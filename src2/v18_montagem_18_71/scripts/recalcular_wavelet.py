"""Calcula sinais, perfis individuais e campos de pares da metodologia atual.

Entradas: 01_panoramas e config/ensaios.json. Saídas: 02_wavelet.
A comparação de Pearson e a seleção dos exemplos pertencem à etapa 04.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import hashlib
from itertools import combinations
import json
import os
from pathlib import Path
import sys
import time

sys.dont_write_bytecode = True
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ.setdefault("MPLCONFIGDIR", "/tmp/tcc-v18-wavelet-mpl")
import numpy as np
from PIL import Image

from runtime import OUTPUT_ROOT, SCRIPTS_ROOT, CONFIG_ROOT, artifact_path, reference
import wavelet_core as core

BASE = OUTPUT_ROOT
ROOT = BASE / "02_wavelet"
INTER = ROOT / "intermediario"
TABLES = ROOT / "04_estatistica/dados"
PATTERNS = ["churn", "wavy annular", "plug", "slug", "dispersed bubble"]
WAVES = np.geomspace(50, 3000, 200)
ROW_START, ROW_STOP, ROW_END, ROW_COUNT = 18, 72, 71, 54
INPUT_MANIFEST = BASE / "01_panoramas/manifesto.json"


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def source_record(path):
    return {"arquivo": reference(path), "sha256": sha(path)}


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


CORE_PATH = SCRIPTS_ROOT / "wavelet_core.py"
CORE_HASH = sha(CORE_PATH)
SCRIPT_HASH = sha(Path(__file__).resolve())


def load_raw(run):
    with np.load(INTER / "entradas" / f"{run}.npz") as archive:
        return archive["raw"].copy()


def load_image(row):
    with Image.open(artifact_path(row["wavelet_path"])) as image:
        return np.asarray(image.convert("L"))


def prepare_inputs():
    manifest = json.loads(INPUT_MANIFEST.read_text())
    assert manifest["complete"] and manifest["alignment_slice"] == [18, 72]
    supplied = manifest if isinstance(manifest, list) else manifest["records"]
    records = []
    for entry in supplied:
        row = dict(entry)
        row["wavelet_path"] = reference(artifact_path(row["wavelet_path"]))
        row["optica_path"] = reference(artifact_path(row["optica_path"]))
        row["padrao"] = row["padrao"].replace("_", " ")
        if row["padrao"] == "dispersed bubbles":
            row["padrao"] = "dispersed bubble"
        image = load_image(row)
        assert sha(artifact_path(row["wavelet_path"])) == row["hashes"]["wavelet"]
        raw = image[18:72].sum(axis=0, dtype=np.uint64) / 54
        np.testing.assert_array_equal(raw, image[18:72].mean(axis=0))
        assert image.shape[0] == 80 and len(raw) == int(row["comprimento"])
        path = INTER / "entradas" / f"{row['ensaio']}.npz"
        np.savez_compressed(path, raw=raw)
        row.update(panorama_sha256=sha(artifact_path(row["wavelet_path"])),
                   raw_sha256=hashlib.sha256(raw.astype("<f8").tobytes()).hexdigest())
        records.append(row)
    records.sort(key=lambda r: (PATTERNS.index(r["padrao"]), r["ensaio"]))
    assert len(records) == 26 and len({r["ensaio"] for r in records}) == 26
    selected = json.loads((CONFIG_ROOT / "ensaios.json").read_text())["runs"]
    assert [(r["ensaio"], r["padrao"]) for r in records] == [(r["run"], r["pattern"]) for r in selected]
    for number, row in enumerate(records, 1):
        row["numero"] = number
    write_json(INTER / "aquisicoes.json", records)
    write_csv(TABLES / "indice_ensaios.csv", [dict(numero=r["numero"], ensaio=r["ensaio"],
              padrao=r["padrao"], comprimento_px=r["comprimento"]) for r in records])
    write_json(INTER / "configuracao.json", dict(
        ensaios=26, pares=325, intra=87, inter=238, linhas_inclusivas=[18, 71],
        wavelengths_px=[50, 3000], escalas=200, omega0=6, amostragem_px=1,
        autocorrelacao="enviesada, divisor N", dominio="todos os atrasos, sem COI ou máscara geométrica",
        selecao_par="prefixo min(Na,Nb), preparado após o corte",
        wps="abs(CWT(Rxx))", wcs="CWT(Rxy)", wco="abs(WCS)^2/(WPSx(0)*WPSy(0))",
        potencia_media="média da WPS em todos os atrasos",
        descritor="potência no atraso zero do sinal completo; Pearson direto na etapa 04",
        core=source_record(CORE_PATH),
        selecao_ensaios=source_record(CONFIG_ROOT / "ensaios.json"),
        manifesto_montagem=source_record(INPUT_MANIFEST)))
    return records


def zero_direct(rxx, wavelengths=WAVES):
    lags = np.arange(-(len(rxx) // 2), len(rxx) // 2 + 1)
    values = []
    for scale in wavelengths / core.FOURIER_FACTOR:
        u = lags / scale
        kernel = np.pi ** (-.25) / np.sqrt(scale) * np.exp(-u * u / 2) * np.cos(6 * u)
        values.append(abs(np.sum(rxx * kernel)))
    return np.asarray(values)


def individual(row):
    run = row["ensaio"]
    marker = INTER / "individuais" / f"{run}.json"
    target = marker.with_suffix(".npz")
    signatures = {"input": row["raw_sha256"], "core": CORE_HASH, "script": SCRIPT_HASH}
    if marker.exists() and target.exists():
        previous = json.loads(marker.read_text())
        if previous["assinaturas"] == signatures and previous["npz_sha256"] == sha(target):
            return previous
    raw = load_raw(run)
    prepared = core.prepare(raw, linear=True)
    rxx = core.correlation(prepared, prepared, estimator="biased")
    wxx, scales = core.transform_correlation(rxx, WAVES, engine="morlet")
    zero = np.abs(wxx[:, len(raw) - 1])
    mean = np.abs(wxx).mean(axis=1)
    independent = zero_direct(rxx)
    np.testing.assert_allclose(zero, independent, atol=1e-11, rtol=1e-11)
    assert abs(prepared.mean()) < 1e-12 and abs(prepared.std() - 1) < 1e-12
    np.savez_compressed(target, wavelengths=WAVES, prepared=prepared, wps_zero=zero, wps_mean_all=mean)
    result = dict(ensaio=run, padrao=row["padrao"], comprimento=len(raw),
                  pico_wps_zero_lambda_px=float(WAVES[np.argmax(zero)]),
                  pico_potencia_media_lambda_px=float(WAVES[np.argmax(mean)]),
                  erro_wps_independente=float(abs(zero - independent).max()),
                  npz_sha256=sha(target), assinaturas=signatures)
    write_json(marker, result)
    return result


def specifications(records):
    specs = []
    for a, b in combinations(sorted(records, key=lambda r: r["ensaio"]), 2):
        same = a["padrao"] == b["padrao"]
        prefix = a["padrao"].replace(" ", "_")
        if not same:
            prefix += "_" + b["padrao"].replace(" ", "_")
        stem = prefix + "_" + a["ensaio"] + "_" + b["ensaio"]
        specs.append(dict(run_x=a["ensaio"], run_y=b["ensaio"], padrao_x=a["padrao"],
                          padrao_y=b["padrao"], tipo="intra" if same else "inter",
                          padrao=a["padrao"] if same else a["padrao"] + " × " + b["padrao"],
                          arquivo="04_analise_simplificada/figuras/" + stem + ".png",
                          input_x_sha256=a["raw_sha256"], input_y_sha256=b["raw_sha256"]))
    assert len(specs) == 325 and sum(s["tipo"] == "intra" for s in specs) == 87
    return specs


def calculate_pair(spec, save_maps=False):
    a, b = spec["run_x"], spec["run_y"]
    stem = INTER / "pares" / f"{a}__{b}"
    marker = stem.with_suffix(".json")
    signatures = dict(core=CORE_HASH, script=SCRIPT_HASH,
                      input_x=spec["input_x_sha256"], input_y=spec["input_y_sha256"])
    if not save_maps and marker.exists():
        previous = json.loads(marker.read_text())
        if previous["assinaturas"] == signatures:
            return previous
    raw_a, raw_b = load_raw(a), load_raw(b)
    n = min(len(raw_a), len(raw_b))
    x, y = core.prepare(raw_a[:n], linear=True), core.prepare(raw_b[:n], linear=True)
    maps = core.pair(x, y, WAVES, estimator="biased", engine="morlet")
    co, j0 = maps["wco"], n - 1
    np.testing.assert_array_equal(maps["lags"], np.arange(-n + 1, n))
    assert maps["denominator_valid"].all() and np.isfinite(co).all()
    assert float(co.min()) >= 0 and float(co.max()) <= 1 + 1e-8
    zx, zy = np.abs(maps["wxx"][:, j0]), np.abs(maps["wyy"][:, j0])
    independent_errors = []
    indices = [0, 31, 79, 127, 165, 199]
    for prepared, zero in ((x, zx), (y, zy)):
        independent = zero_direct(core.correlation(prepared, prepared, estimator="biased"), WAVES[indices])
        independent_errors.append(float(abs(independent - zero[indices]).max()))
        np.testing.assert_allclose(independent, zero[indices], atol=1e-11, rtol=1e-11)
    rxy = core.correlation(x, y, estimator="biased")
    direct_errors = []
    for lag in (-n + 1, -101, -1, 0, 1, 101, n - 1):
        direct = np.dot(x[lag:], y[:n - lag]) / n if lag >= 0 else np.dot(x[:n + lag], y[-lag:]) / n
        direct_errors.append(float(abs(rxy[j0 + lag] - direct)))
    assert max(direct_errors) < 1e-12
    peak_scale, peak_lag = np.unravel_index(np.argmax(co), co.shape)
    metrics = dict(comprimento=n, comprimento_x=len(raw_a), comprimento_y=len(raw_b),
                   descarte_x=len(raw_a) - n, descarte_y=len(raw_b) - n,
                   escalas_validas=int(maps["denominator_valid"].sum()),
                   wco_maxima=float(co.max()), wco_minima=float(co.min()),
                   wco_media_integral=float(co.mean()),
                   wco_maxima_zero=float(co[:, j0].max()), wco_mediana_zero=float(np.median(co[:, j0])),
                   pico_lag_px=int(maps["lags"][peak_lag]), pico_lambda_px=float(WAVES[peak_scale]))
    result = dict(**spec, **metrics, erro_soma_direta=max(direct_errors),
                  erro_wps_independente=max(independent_errors), assinaturas=signatures)
    if save_maps:
        previous = json.loads(marker.read_text())
        assert all(previous[key] == value for key, value in metrics.items())
        step = max(1, int(np.ceil(len(maps["lags"]) / 1800)))
        idx = np.unique(np.r_[np.arange(0, len(maps["lags"]), step), j0, len(maps["lags"]) - 1])
        compact = {key: maps[key][:, idx] for key in ("wxx", "wyy", "wxy", "wco")}
        compact.update(lags=maps["lags"][idx], wavelengths=WAVES, scales=maps["scales"], x=x, y=y,
                       mean_x_all=np.abs(maps["wxx"]).mean(axis=1),
                       mean_y_all=np.abs(maps["wyy"]).mean(axis=1),
                       max_wxx=np.abs(maps["wxx"]).max(), max_wyy=np.abs(maps["wyy"]).max(),
                       max_wxy=np.abs(maps["wxy"]).max(), wps_zero_x=zx, wps_zero_y=zy,
                       scale_peak_wco=co.max(axis=1), full_lag_count=len(maps["lags"]))
        np.savez_compressed(stem.with_suffix(".npz"), **compact)
        result["npz_sha256"] = sha(stem.with_suffix(".npz"))
    write_json(marker, result)
    return result


def finalize(records, pairs):
    write_json(INTER / "pares.json", pairs)
    summary = dict(ensaios=records, execucao_pares_wavelet_completa=True,
                   extensao_minima=min(r["comprimento"] for r in records),
                   extensao_maxima=max(r["comprimento"] for r in records),
                   counts=dict(individuals=len(records), pairs=len(pairs), intra=87, inter=238, wavelengths=200),
                   estatistica="04_analise_simplificada: Pearson direto dos perfis completos")
    write_json(ROOT / "resumo_tcc.json", summary)
    write_json(INTER / "auditoria/wavelet.json", dict(status="PASS", pares=len(pairs),
        maior_erro_correlacao_soma_direta=max(p["erro_soma_direta"] for p in pairs),
        maior_erro_wps_soma_direta=max(p["erro_wps_independente"] for p in pairs),
        menor_coerencia=min(p["wco_minima"] for p in pairs), maior_coerencia=max(p["wco_maxima"] for p in pairs),
        todas_200_escalas_validas=all(p["escalas_validas"] == 200 for p in pairs),
        mapas_completos_usados_em_metricas=True, manifesto_montagem=source_record(INPUT_MANIFEST)))


def main(workers):
    started = time.perf_counter()
    for relative in ("entradas", "individuais", "pares", "auditoria"):
        (INTER / relative).mkdir(parents=True, exist_ok=True)
    for relative in ("01_panoramas", "02_wavelets", "03_wavelets_inter", "04_estatistica/dados"):
        (ROOT / relative).mkdir(parents=True, exist_ok=True)
    records = prepare_inputs()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(individual, row) for row in records]
        individuals = []
        for f in as_completed(futures):
            value = f.result()
            individuals.append(value)
            print(f"Individual {len(individuals)}/26 {value['ensaio']}", flush=True)
    individuals.sort(key=lambda row: row["ensaio"])
    write_json(INTER / "individuais.json", individuals)
    specs = specifications(records)
    pairs = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(calculate_pair, spec) for spec in specs]
        for future in as_completed(futures):
            pairs.append(future.result())
            if len(pairs) % 10 == 0 or len(pairs) == 325:
                print(f"Pares {len(pairs)}/325, {time.perf_counter()-started:.1f}s", flush=True)
    pairs.sort(key=lambda p: (p["run_x"], p["run_y"]))
    finalize(records, pairs)
    print(f"WAVELET_COMPLETO {time.perf_counter()-started:.1f}s", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=2)
    main(parser.parse_args().workers)
