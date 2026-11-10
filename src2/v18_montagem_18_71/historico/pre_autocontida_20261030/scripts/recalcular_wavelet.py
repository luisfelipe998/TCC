"""Recalcula a análise v16 sobre panoramas remontados com a faixa 18:72.

python3 src2/v18_montagem_18_71/scripts/recalcular_wavelet.py --workers 2

O núcleo Morlet é carregado da cópia histórica sem executar outros programas.
Os 325 pares usam todos os atrasos e escalas. Apenas os cinco pares escolhidos
pelos critérios do manuscrito conservam mapas compactos para apresentação.
"""
import argparse
import ast
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import hashlib
import importlib.util
from itertools import combinations
import json
import os
from pathlib import Path
import sys
import time

sys.dont_write_bytecode = True
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MPLCONFIGDIR"] = "/tmp/tcc-v18-wavelet-mpl"
import numpy as np
from PIL import Image
from scipy.stats import pearsonr

BASE = Path(__file__).resolve().parents[1]
REPO = BASE.parents[1]
OLD = REPO / "src2/v16_spatial2_18_71"
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
    return {"arquivo": str(path.relative_to(REPO)), "sha256": sha(path)}


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


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


core = load_module(OLD / "scripts/wavelet_core.py", "v18_morlet_historico")
CORE_HASH = sha(OLD / "scripts/wavelet_core.py")
SCRIPT_HASH = sha(Path(__file__).resolve())


def load_functions(path, names, scope):
    nodes = [node for node in ast.parse(path.read_text()).body
             if isinstance(node, ast.FunctionDef) and node.name in names]
    assert {node.name for node in nodes} == set(names)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), "exec"), scope)
    return scope


STAT_FUNCTIONS = load_functions(OLD / "scripts/estatistica.py",
                                {"residuals", "pearson_matrix", "describe"}, {"np": np})


def load_raw(run):
    with np.load(INTER / "entradas" / f"{run}.npz") as archive:
        return archive["raw"].copy()


def load_image(row):
    with Image.open(REPO / row["wavelet_path"]) as image:
        return np.asarray(image.convert("L"))


def prepare_inputs():
    manifest = json.loads(INPUT_MANIFEST.read_text())
    assert manifest["complete"] and manifest["alignment_slice"] == [18, 72]
    supplied = manifest if isinstance(manifest, list) else manifest["records"]
    records = []
    for entry in supplied:
        row = dict(entry)
        row["padrao"] = row["padrao"].replace("_", " ")
        if row["padrao"] == "dispersed bubbles":
            row["padrao"] = "dispersed bubble"
        image = load_image(row)
        assert sha(REPO / row["wavelet_path"]) == row["hashes"]["wavelet"]
        raw = image[18:72].sum(axis=0, dtype=np.uint64) / 54
        np.testing.assert_array_equal(raw, image[18:72].mean(axis=0))
        assert image.shape[0] == 80 and len(raw) == int(row["comprimento"])
        path = INTER / "entradas" / f"{row['ensaio']}.npz"
        np.savez_compressed(path, raw=raw)
        row.update(panorama_sha256=sha(REPO / row["wavelet_path"]),
                   raw_sha256=hashlib.sha256(raw.astype("<f8").tobytes()).hexdigest())
        records.append(row)
    records.sort(key=lambda r: (PATTERNS.index(r["padrao"]), r["ensaio"]))
    assert len(records) == 26 and len({r["ensaio"] for r in records}) == 26
    old_records = json.loads((OLD / "intermediario/aquisicoes.json").read_text())
    assert [(r["ensaio"], r["padrao"]) for r in records] == [(r["ensaio"], r["padrao"]) for r in old_records]
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
        pearson="resíduos de log(WPS(0)+1e-12), reta contra 1 e log(lambda), ensaios completos",
        escores_ponderados="não calculados, retirados do manuscrito",
        core=source_record(OLD / "scripts/wavelet_core.py"),
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
    prepared = core.prepare(raw)
    rxx = core.correlation(prepared, prepared)
    wxx, scales = core.transform_correlation(rxx, WAVES)
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
                          arquivo=("02_wavelets/" if same else "03_wavelets_inter/") + stem + ".png",
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
    x, y = core.prepare(raw_a[:n]), core.prepare(raw_b[:n])
    maps = core.pair(x, y, WAVES)
    co, j0 = maps["wco"], n - 1
    np.testing.assert_array_equal(maps["lags"], np.arange(-n + 1, n))
    assert maps["denominator_valid"].all() and np.isfinite(co).all()
    assert float(co.min()) >= 0 and float(co.max()) <= 1 + 1e-8
    zx, zy = np.abs(maps["wxx"][:, j0]), np.abs(maps["wyy"][:, j0])
    independent_errors = []
    indices = [0, 31, 79, 127, 165, 199]
    for prepared, zero in ((x, zx), (y, zy)):
        independent = zero_direct(core.correlation(prepared, prepared), WAVES[indices])
        independent_errors.append(float(abs(independent - zero[indices]).max()))
        np.testing.assert_allclose(independent, zero[indices], atol=1e-11, rtol=1e-11)
    rxy = core.correlation(x, y)
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


def calculate_statistics(records, specs):
    runs = [r["ensaio"] for r in records]
    powers, means = [], []
    for run in runs:
        with np.load(INTER / "individuais" / f"{run}.npz") as saved:
            np.testing.assert_array_equal(saved["wavelengths"], WAVES)
            powers.append(saved["wps_zero"])
            means.append(saved["wps_mean_all"])
    powers, means = np.asarray(powers), np.asarray(means)
    z, coefficients = STAT_FUNCTIONS["residuals"](powers, WAVES)
    r, unit = STAT_FUNCTIONS["pearson_matrix"](z)
    np.testing.assert_allclose(r, np.corrcoef(z), atol=1e-13, rtol=1e-13)
    errors = [abs(r[i, j] - pearsonr(z[i], z[j]).statistic) for i, j in combinations(range(26), 2)]
    assert max(errors) < 1e-12 and np.linalg.eigvalsh(r).min() > -1e-10
    np.savez_compressed(INTER / "descritores.npz", runs=np.asarray(runs), wavelengths=WAVES,
                        power=powers, power_mean_all=means, residuals=z, unit=unit,
                        coefficients=coefficients, r=r)
    write_csv(TABLES / "perfis_espectrais.csv", [dict(ensaio=run, padrao=records[i]["padrao"],
              lambda_px=float(w), wps_zero=float(powers[i, k]), potencia_media_todos_lags=float(means[i, k]),
              residuo=float(z[i, k]), descritor_unitario=float(unit[i, k]))
              for i, run in enumerate(runs) for k, w in enumerate(WAVES)])
    write_csv(TABLES / "matriz_pearson_todos.csv", [dict(ensaio=run,
              **{other: float(r[i, j]) for j, other in enumerate(runs)}) for i, run in enumerate(runs)])
    rows = []
    for pair in specs:
        i, j = runs.index(pair["run_x"]), runs.index(pair["run_y"])
        rows.append(dict(numero_x=i + 1, numero_y=j + 1, ensaio_x=runs[i], ensaio_y=runs[j],
                         padrao_x=pair["padrao_x"], padrao_y=pair["padrao_y"], tipo=pair["tipo"],
                         comprimento_comum=min(records[i]["comprimento"], records[j]["comprimento"]),
                         pearson_r=float(r[i, j]), arquivo=pair["arquivo"]))
    write_csv(TABLES / "pares_pearson_wco.csv", rows)
    groups = {}
    for row in rows:
        key = tuple(sorted((row["padrao_x"], row["padrao_y"]), key=PATTERNS.index))
        groups.setdefault(key, []).append(row)
    group_rows = [dict(padrao_x=a, padrao_y=b, tipo="intra" if a == b else "inter", pares=len(group),
                      pearson_media=float(np.mean([v["pearson_r"] for v in group])),
                      pearson_mediana=float(np.median([v["pearson_r"] for v in group])))
                  for (a, b), group in sorted(groups.items(), key=lambda v: (PATTERNS.index(v[0][0]), PATTERNS.index(v[0][1])))]
    write_csv(TABLES / "resumo_por_padrao.csv", group_rows)
    aggregate = [dict(tipo=kind, metrica="pearson_r", **STAT_FUNCTIONS["describe"](
                     [v["pearson_r"] for v in rows if v["tipo"] == kind])) for kind in ("intra", "inter")]
    write_csv(TABLES / "resumo_intra_inter.csv", aggregate)
    averages = []
    for i, acquisition in enumerate(records):
        others = np.arange(26) != i
        same = np.array([v["padrao"] == acquisition["padrao"] for v in records]) & others
        different = others & ~same
        averages.append(dict(numero=i + 1, ensaio=runs[i], padrao=acquisition["padrao"],
            pearson_media_todos=float(r[i, others].mean()),
            pearson_media_intra=float(r[i, same].mean()) if same.any() else None,
            pearson_media_inter=float(r[i, different].mean()),
            pares_intra=int(same.sum()), pares_inter=int(different.sum())))
    write_csv(TABLES / "medias_por_ensaio.csv", averages)
    selected = [max([row for row in rows if row["tipo"] == "intra" and row["padrao_x"] == pattern],
                    key=lambda row: row["pearson_r"]) for pattern in PATTERNS[:-1]]
    write_json(TABLES / "selecao_intrapadrao.json", selected)
    summary = dict(ensaios=records, grupos=group_rows, resumos=aggregate,
                   diferenca_global=aggregate[0]["media"] - aggregate[1]["media"],
                   medias_por_ensaio=averages, pares_intrapadrao_selecionados=selected,
                   extensao_minima=min(v["comprimento"] for v in records),
                   extensao_maxima=max(v["comprimento"] for v in records),
                   menor_anular=min((v for v in averages if v["padrao"] == "wavy annular"),
                                    key=lambda v: v["pearson_media_intra"]),
                   m03q40=next(v for v in averages if v["ensaio"] == "m03q40_001"),
                   execucao_pares_wavelet_completa=False)
    write_json(ROOT / "resumo_tcc.json", summary)
    write_json(INTER / "auditoria/estatistica.json", dict(status="PASS", pares_scipy=325,
               maior_erro_scipy=max(errors), menor_autovalor=float(np.linalg.eigvalsh(r).min()),
               funcoes_historicas=source_record(OLD / "scripts/estatistica.py")))
    return summary


def finalize(records, specs, summary, pairs):
    by_pair = {(p["run_x"], p["run_y"]): p for p in pairs}
    rows = read_csv(TABLES / "pares_pearson_wco.csv")
    for row in rows:
        p = by_pair[row["ensaio_x"], row["ensaio_y"]]
        for key in ("wco_maxima", "wco_media_integral", "wco_maxima_zero", "wco_mediana_zero", "pico_lag_px", "pico_lambda_px"):
            row[key] = p[key]
    write_csv(TABLES / "pares_pearson_wco.csv", rows)
    inter = max((p for p in pairs if p["tipo"] == "inter"), key=lambda p: p["wco_maxima"])
    lookup = {(p["run_x"], p["run_y"]): p for p in specs}
    selected_specs = [lookup[p["ensaio_x"], p["ensaio_y"]] for p in summary["pares_intrapadrao_selecionados"]]
    selected_specs.append(lookup[inter["run_x"], inter["run_y"]])
    chosen = []
    for spec in selected_specs:
        value = calculate_pair(spec, save_maps=True)
        p = next(p for p in rows if (p["ensaio_x"], p["ensaio_y"]) == (spec["run_x"], spec["run_y"]))
        chosen.append(dict(**value, pearson_r=float(p["pearson_r"])))
    write_json(INTER / "pares.json", pairs)
    write_json(TABLES / "pares_ilustrados.json", chosen)
    summary.update(execucao_pares_wavelet_completa=True, pares_ilustrados=chosen,
                   maximo_wco_interpadrões=chosen[-1],
                   comprimento_minimo_sobre_lambda_maxima=summary["extensao_minima"] / 3000)
    old_aggregate = read_csv(OLD / "04_estatistica/dados/resumo_intra_inter.csv")
    old_means = {v["tipo"]: float(v["media"]) for v in old_aggregate if v["metrica"] == "pearson_r"}
    summary["comparacao_antigo_novo"] = {
        "media_intra": {"antigo": old_means["intra"], "novo": summary["resumos"][0]["media"]},
        "media_inter": {"antigo": old_means["inter"], "novo": summary["resumos"][1]["media"]},
        "diferenca_global": {"antigo": old_means["intra"] - old_means["inter"], "novo": summary["diferenca_global"]}}
    write_json(ROOT / "resumo_tcc.json", summary)
    write_json(INTER / "auditoria/wavelet.json", dict(status="PASS", pares=len(pairs),
        maior_erro_correlacao_soma_direta=max(p["erro_soma_direta"] for p in pairs),
        maior_erro_wps_soma_direta=max(p["erro_wps_independente"] for p in pairs),
        menor_coerencia=min(p["wco_minima"] for p in pairs), maior_coerencia=max(p["wco_maxima"] for p in pairs),
        todas_200_escalas_validas=all(p["escalas_validas"] == 200 for p in pairs),
        criterio_intra="maior Pearson de cada padrão com pelo menos dois ensaios",
        criterio_inter="maior WCO pontual entre os 238 pares interpadrões",
        mapas_completos_usados_em_metricas=True, compactacao_somente_para_apresentacao=True,
        mapas_compactos_salvos=len(chosen), manifesto_montagem=source_record(INPUT_MANIFEST)))


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
    write_json(INTER / "individuais.json", individuals)
    specs = specifications(records)
    summary = calculate_statistics(records, specs)
    print("PEARSON_PRONTO " + json.dumps({"diferenca": summary["diferenca_global"],
          "pares": summary["pares_intrapadrao_selecionados"]}), flush=True)
    pairs = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(calculate_pair, spec) for spec in specs]
        for future in as_completed(futures):
            pairs.append(future.result())
            if len(pairs) % 10 == 0 or len(pairs) == 325:
                print(f"Pares {len(pairs)}/325, {time.perf_counter()-started:.1f}s", flush=True)
    pairs.sort(key=lambda p: (p["run_x"], p["run_y"]))
    finalize(records, specs, summary, pairs)
    print(f"WAVELET_COMPLETO {time.perf_counter()-started:.1f}s", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=2)
    main(parser.parse_args().workers)
