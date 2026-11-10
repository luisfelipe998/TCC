"""Reusable individual spectra and Morlet kernels for equal-domain comparisons."""
import hashlib
import json
import copy
from datetime import datetime, timezone
from functools import lru_cache
from numbers import Integral
from pathlib import Path

import numpy as np
from scipy.fft import fft, ifft, fftfreq, next_fast_len

from common import CFG, INTER, SCRIPTS, context_ok, image_path, load_image, sha, write_json
from wavelet_core import OMEGA0, FOURIER_FACTOR, correlation, prepare, safe_coherence

CACHE = INTER / "precalculo_wavelets"
BACKEND_SHA = sha(Path(__file__))
PREPARATION = "recorte, 255-I, detrend linear, remocao de media, desvio padrao unitario"


def validate_method():
    """Fail explicitly instead of silently ignoring unsupported scientific options."""
    supported = dict(estimador="biased", K=1, omega0=OMEGA0,
                     comprimento_onda_maximo="N/2 por padrao")
    for name, expected in supported.items():
        if CFG["wavelet"].get(name) != expected:
            raise ValueError(f"Parametro nao suportado: {name}; este metodo exige {expected!r}.")
    scales = CFG["wavelet"]["escalas"]
    low = CFG["wavelet"]["comprimento_onda_minimo"]
    if not isinstance(scales, Integral) or isinstance(scales, bool) or scales < 2:
        raise ValueError("Numero de escalas invalido.")
    if not np.isfinite(low) or low <= 0:
        raise ValueError("Comprimento de onda minimo invalido.")


def validate_length(length):
    validate_method()
    if not isinstance(length, Integral) or isinstance(length, bool):
        raise ValueError("Comprimento deve ser inteiro, sem arredondamento implicito.")
    length = int(length)
    step = CFG["contexto"]["passo"]
    if length < CFG["contexto"]["comprimento_minimo"] or length % 2 or length % step:
        raise ValueError(f"Comprimento deve ser par e multiplo de {step}, respeitando o minimo.")
    if length < 2*CFG["wavelet"]["comprimento_onda_minimo"]:
        raise ValueError("Janela insuficiente para a faixa de comprimentos de onda.")
    return length


def implementation():
    return dict(nucleo_sha256=sha(SCRIPTS / "wavelet_core.py"),
                entrada_sha256=sha(SCRIPTS / "common.py"), backend_sha256=BACKEND_SHA)


def signature(spec):
    return hashlib.sha256(json.dumps(spec, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _fingerprint(folder, meta):
    files = [folder / "metadata.json"] + [folder / f"{name}.npy" for name in sorted(meta["arrays_sha256"])]
    return tuple((p.name, p.stat().st_size, p.stat().st_mtime_ns, p.stat().st_ctime_ns) for p in files)


@lru_cache(maxsize=128)
def _load_checked(folder_name, fingerprint):
    folder = Path(folder_name)
    meta = json.loads((folder / "metadata.json").read_text())
    if signature(meta["spec"]) != folder.name:
        raise ValueError("Assinatura invalida do precalculo.")
    for field, expected in implementation().items():
        if meta["spec"].get(field) != expected:
            raise ValueError("Precalculo de outra implementacao; regenere os sinais individuais.")
    result = dict(key=folder.name, meta=meta)
    for name, expected in meta["arrays_sha256"].items():
        path = folder / f"{name}.npy"
        if sha(path) != expected:
            raise ValueError(f"Precalculo corrompido: {path}")
        result[name] = np.load(path, mmap_mode="r", allow_pickle=False)
    if _fingerprint(folder, meta) != fingerprint:
        raise ValueError("Precalculo mudou durante a leitura.")
    return result


def load_entry(folder_name):
    validate_method()
    folder = Path(folder_name)
    meta = json.loads((folder / "metadata.json").read_text())
    result = _load_checked(str(folder), _fingerprint(folder, meta))
    # Returned metadata cannot modify the object reused by the process.
    return dict(result, meta=copy.deepcopy(result["meta"]))


def write_entry(folder, spec, arrays):
    folder.mkdir(parents=True, exist_ok=True)
    hashes = {}
    for name, values in arrays.items():
        path = folder / f"{name}.npy"
        temporary = folder / f"{name}.npy.tmp"
        with temporary.open("wb") as file:
            np.save(file, values, allow_pickle=False)
        temporary.replace(path)
        hashes[name] = sha(path)
    metadata = dict(spec=spec, arrays_sha256=hashes, data_utc=datetime.now(timezone.utc).isoformat())
    temporary = folder / "metadata.json.tmp"
    write_json(temporary, metadata)
    temporary.replace(folder / "metadata.json")
    return load_entry(str(folder))


def grid_spec(length):
    length = validate_length(length)
    return dict(tipo="grade", schema=2, comprimento=length, dx_px=1,
                lambda_min=CFG["wavelet"]["comprimento_onda_minimo"], lambda_max=length/2,
                escalas=CFG["wavelet"]["escalas"], omega0=OMEGA0, padding_sigma=8,
                estimador="biased", K=1, preparo=PREPARATION, **implementation())


def get_grid(length):
    spec = grid_spec(length)
    folder = CACHE / "grades" / signature(spec)
    if (folder / "metadata.json").exists():
        return load_entry(str(folder))
    waves = np.geomspace(spec["lambda_min"], spec["lambda_max"], spec["escalas"])
    scales = waves / FOURIER_FACTOR
    lags = np.arange(-length+1, length)
    support = length-1-np.abs(lags)[None, :] >= np.sqrt(2)*scales[:, None]
    nfft = next_fast_len(len(lags)+2*int(np.ceil(spec["padding_sigma"]*scales.max())))
    omega = 2*np.pi*fftfreq(nfft)
    kernels = np.empty((len(scales), nfft), dtype=float)
    for index, scale in enumerate(scales):
        kernels[index] = np.sqrt(2*np.pi*scale)*np.pi**(-.25)*np.exp(-.5*(scale*omega-OMEGA0)**2)
    return write_entry(folder, spec, dict(wavelengths=waves, scales=scales, lags=lags,
                                        support=support, kernels=kernels))


def transform_with_grid(values, grid):
    if len(values) != len(grid["lags"]):
        raise ValueError("Correlacao incompativel com a grade.")
    spectrum = fft(values, grid["kernels"].shape[1])
    transformed = np.empty((len(grid["scales"]), len(values)), dtype=complex)
    for index, kernel in enumerate(grid["kernels"]):
        transformed[index] = ifft(spectrum*kernel)[:len(values)]
    return transformed


def context_record(record, length):
    from anchor_diagnostics import context_diagnostic_fields
    length = validate_length(length)
    center = record["ancora"]
    if length < record["comprimento_minimo"] or not context_ok(center, record["tamanho_panorama"], length):
        raise ValueError(f"Contexto de {length} px inviavel para {record['anchor_id']}; ancora preservada.")
    result = dict(record, comprimento=length, segmento_inicio=center-length//2, segmento_fim=center+length//2)
    if record.get("diagnosticos_comprimento") != length:
        for field in context_diagnostic_fields():
            result.pop(field, None)
    return result


def individual_spec(record, length=None):
    length = record["comprimento"] if length is None else length
    record = context_record(record, length)
    source = image_path(record["aquisicao"])
    return dict(tipo="individual", schema=2, aquisicao=record["aquisicao"], ancora=record["ancora"],
                comprimento=length, inicio=record["segmento_inicio"], fim=record["segmento_fim"],
                linhas=list(CFG["linhas"]), panorama_sha256=sha(source), grade=signature(grid_spec(length)),
                preparo=PREPARATION, **implementation())


def get_auto(record, length=None):
    spec = individual_spec(record, length)
    folder = CACHE / "individuais" / signature(spec)
    if (folder / "metadata.json").exists():
        return load_entry(str(folder))
    grid = get_grid(spec["comprimento"])
    _, raw = load_image(record["aquisicao"])
    raw = raw[spec["inicio"]:spec["fim"]]
    x = prepare(raw, linear=True)
    rxx = correlation(x, x, estimator="biased")
    wxx = transform_with_grid(rxx, grid)
    mean = np.sum(np.where(grid["support"], np.abs(wxx), 0), axis=1) / grid["support"].sum(axis=1)
    return write_entry(folder, spec, dict(raw=raw, x=x, rxx=rxx, wxx=wxx, power_mean=mean))


def compatible_grid(left, right):
    validate_method()
    a, b = left["meta"]["spec"], right["meta"]["spec"]
    for item in (a,b):
        if any(item.get(key) != expected for key, expected in implementation().items()):
            raise ValueError("Implementacao incompativel com os precalculos.")
    for field in ("comprimento", "grade", "linhas", "preparo", "nucleo_sha256", "entrada_sha256", "backend_sha256"):
        if a[field] != b[field]:
            raise ValueError(f"Comparacao incompativel ({field}); use o mesmo perfil de contexto e escala.")
    grid = get_grid(a["comprimento"])
    if grid["key"] != a["grade"]:
        raise ValueError("Parametros da grade mudaram; regenere os precalculos.")
    return grid


def pair_from_autos(left, right):
    """Pattern-independent engine. No partner selection or label-based changes."""
    grid = compatible_grid(left, right)
    length = left["meta"]["spec"]["comprimento"]
    cross = correlation(left["x"], right["x"], estimator="biased")
    wxy = transform_with_grid(cross, grid)
    wco, valid = safe_coherence(wxy, left["wxx"], right["wxx"], length-1)
    return dict(wxx=left["wxx"], wyy=right["wxx"], wxy=wxy, wco=wco,
                scales=grid["scales"], lags=grid["lags"], denominator_valid=valid,
                wavelengths=grid["wavelengths"], support=grid["support"],
                x=left["x"], y=right["x"], raw_x=left["raw"], raw_y=right["raw"],
                power_mean_x=left["power_mean"], power_mean_y=right["power_mean"])
