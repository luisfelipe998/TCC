#!/usr/bin/env python3
"""Complementa as estatísticas dos pares sem alterar as etapas anteriores.

Recalcula a WCO completa para obter sua mediana exata. Pearson é lido da
etapa 04 e conferido com sua matriz. Não executa inferência ou sensibilidade.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import csv
import os
from pathlib import Path
import platform
import sys
import time

sys.dont_write_bytecode = True
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ.setdefault("MPLCONFIGDIR", "/tmp/tcc-v18-resumos-mpl")

import numpy as np
import scipy

from runtime import OUTPUT_ROOT, SCRIPTS_ROOT, CONFIG_ROOT, reference
import wavelet_core as core

OUT = OUTPUT_ROOT / "06_resumos"
INTER = OUTPUT_ROOT / "02_wavelet/intermediario"
PREVIOUS = OUTPUT_ROOT / "04_analise_simplificada"
PATTERNS = ["churn", "wavy annular", "plug", "slug", "dispersed bubble"]
METRICS = ["pearson_r", "wco_media_integral"]
WAVES = np.geomspace(50., 3000., 200)
COMPARE = ["wco_minima", "wco_media_integral", "wco_maxima",
           "wco_maxima_zero", "wco_mediana_zero", "pico_lambda_px", "pico_lag_px"]


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def tree_hashes():
    """Prove that this supplementary stage did not change existing products."""
    found = {}
    for directory in (OUTPUT_ROOT / "02_wavelet", OUTPUT_ROOT / "03_optica",
                      PREVIOUS, CONFIG_ROOT):
        for path in sorted(directory.rglob("*")):
            if path.is_file() and not {"historico", "__pycache__"}.intersection(path.parts):
                found[reference(path)] = sha(path)
    return found


def calculate_pair(task):
    old, pearson, signatures = task
    a, b = old["run_x"], old["run_y"]
    marker = OUT / "intermediario/pares" / f"{a}__{b}.json"
    if marker.is_file():
        saved = json.loads(marker.read_text())
        if saved.get("assinaturas") == signatures:
            assert saved["valores"]["pearson_r"] == pearson
            for key in COMPARE:
                np.testing.assert_allclose(saved["conferencia"][key]["obtido"], old[key],
                                           atol=1e-12, rtol=1e-12, err_msg=key)
            return saved
    raw = []
    for run, hash_key in ((a, "input_x_sha256"), (b, "input_y_sha256")):
        with np.load(INTER / "entradas" / f"{run}.npz") as archive:
            values = archive["raw"].copy()
        assert hashlib.sha256(values.astype("<f8").tobytes()).hexdigest() == old[hash_key]
        raw.append(values)
    n = min(map(len, raw))
    assert n == old["comprimento"]
    x, y = (core.prepare(values[:n], linear=True) for values in raw)
    maps = core.pair(x, y, WAVES, estimator="biased", engine="morlet")
    co, center = maps["wco"], n - 1
    np.testing.assert_array_equal(maps["lags"], np.arange(-n + 1, n))
    assert co.shape == (200, 2*n - 1)
    assert maps["denominator_valid"].all() and np.isfinite(co).all()
    assert co.min() >= 0 and co.max() <= 1 + 1e-8
    peak = np.unravel_index(np.argmax(co), co.shape)
    computed = dict(wco_minima=float(co.min()), wco_media_integral=float(co.mean()),
                    wco_maxima=float(co.max()), wco_maxima_zero=float(co[:, center].max()),
                    wco_mediana_zero=float(np.median(co[:, center])),
                    pico_lambda_px=float(WAVES[peak[0]]), pico_lag_px=int(maps["lags"][peak[1]]))
    checked = {}
    for key in COMPARE:
        np.testing.assert_allclose(computed[key], old[key], atol=1e-12, rtol=1e-12, err_msg=key)
        checked[key] = dict(obtido=computed[key], referencia=old[key],
                            erro_absoluto=abs(computed[key] - old[key]))
    # Retain only WCO before the median partition, releasing the three complex maps.
    del maps, x, y, raw
    median = float(np.median(co))
    assert computed["wco_minima"] <= median <= computed["wco_maxima"]
    values = dict(ensaio_x=a, ensaio_y=b, padrao_x=old["padrao_x"], padrao_y=old["padrao_y"],
                  tipo=old["tipo"], comprimento_x=old["comprimento_x"], comprimento_y=old["comprimento_y"],
                  comprimento_comum=n, atraso_minimo=-n+1, atraso_maximo=n-1,
                  numero_atrasos=2*n-1, numero_escalas=200, escalas_validas=200,
                  numero_valores=co.size, pearson_r=pearson,
                  wco_minima=computed["wco_minima"], wco_media_integral=computed["wco_media_integral"],
                  wco_mediana_integral=median, wco_maxima=computed["wco_maxima"],
                  pico_lambda_px=computed["pico_lambda_px"], pico_lag_px=computed["pico_lag_px"])
    del co
    result = dict(status="PASS", valores=values, conferencia=checked, assinaturas=signatures,
                  mediana="numpy.median de todas as células do campo completo, sem redução de colunas")
    write_json(marker, result)
    return result


def describe(rows, metric):
    values = np.asarray([row[metric] for row in rows])
    acquisitions = {row[key] for row in rows for key in ("ensaio_x", "ensaio_y")}
    q25, q75 = np.quantile(values, [.25, .75], method="linear")
    return dict(metrica=metric, pares=len(rows), ensaios=len(acquisitions),
                minimo=float(values.min()), media=float(values.mean()), mediana=float(np.median(values)),
                maximo=float(values.max()), q25=float(q25), q75=float(q75))


def summarize(rows, records, examples):
    groups = {}
    for row in rows:
        key = tuple(sorted((row["padrao_x"], row["padrao_y"]), key=PATTERNS.index))
        groups.setdefault(key, []).append(row)
    assert len(groups) == 14
    combinations = []
    counts = {pattern: sum(r["padrao"] == pattern for r in records) for pattern in PATTERNS}
    for (a, b), members in sorted(groups.items(), key=lambda item: tuple(map(PATTERNS.index, item[0]))):
        for metric in METRICS:
            combinations.append(dict(padrao_x=a, padrao_y=b, tipo="intra" if a == b else "inter",
                ensaios_padrao_x=counts[a], ensaios_padrao_y=counts[b], **describe(members, metric)))
    global_rows = [dict(tipo=kind, **describe([row for row in rows if row["tipo"] == kind], metric))
                   for metric in METRICS for kind in ("intra", "inter")]
    by_pair = {(row["ensaio_x"], row["ensaio_y"]): row for row in rows}
    selected = [by_pair[item["run_x"], item["run_y"]] for item in examples]
    assert len(selected) == 5
    differences = {}
    for metric in METRICS:
        inside, outside = [row for row in global_rows if row["metrica"] == metric]
        differences[metric] = inside["media"] - outside["media"]
    write_csv(OUT / "dados/pares_pearson_coerencia.csv", rows)
    write_csv(OUT / "dados/resumo_por_combinacao.csv", combinations)
    write_csv(OUT / "dados/resumo_intra_inter.csv", global_rows)
    write_csv(OUT / "dados/pares_ilustrados.csv", selected)
    summary = dict(ensaios=26, pares=325, intra=87, inter=238, combinacoes=14,
                   resumo_intra_inter=global_rows, resumo_por_combinacao=combinations,
                   pares_ilustrados=selected, diferenca_media_intra_menos_inter=differences)
    write_json(OUT / "resumo_tcc.json", summary)
    return summary


def main(workers):
    if workers < 1:
        raise ValueError("O número de processos deve ser positivo.")
    started = time.monotonic()
    for directory in ("dados", "intermediario/pares", "figuras"):
        (OUT / directory).mkdir(parents=True, exist_ok=True)
    before = tree_hashes()
    inputs = {}

    def track(path):
        inputs[reference(path)] = sha(path)
        return path

    old_pairs = json.loads(track(INTER / "pares.json").read_text())
    records = json.loads(track(INTER / "aquisicoes.json").read_text())
    examples = json.loads(track(PREVIOUS / "dados/pares_ilustrados.json").read_text())
    configuration = json.loads(track(INTER / "configuracao.json").read_text())
    assert configuration["escalas"] == 200 and configuration["wavelengths_px"] == [50, 3000]
    assert len(old_pairs) == 325 and len(records) == 26
    with track(PREVIOUS / "dados/pares_pearson_wco.csv").open(newline="") as stream:
        pearson_rows = list(csv.DictReader(stream))
    pearson = {(row["ensaio_x"], row["ensaio_y"]): float(row["pearson_r"]) for row in pearson_rows}
    assert len(pearson_rows) == len(pearson) == 325
    with np.load(track(PREVIOUS / "intermediario/descritores.npz")) as archive:
        np.testing.assert_array_equal(archive["wavelengths"], WAVES)
        matrix, runs = archive["r"].copy(), archive["runs"].tolist()
    index = {run: i for i, run in enumerate(runs)}
    for (a, b), value in pearson.items():
        assert value == float(matrix[index[a], index[b]])
    scripts = {reference(path): sha(path) for path in
               (Path(__file__), SCRIPTS_ROOT / "wavelet_core.py", SCRIPTS_ROOT / "figuras_resumos.py",
                SCRIPTS_ROOT / "runtime.py")}
    task_signatures = dict(core=sha(SCRIPTS_ROOT / "wavelet_core.py"), script=sha(Path(__file__)))
    tasks = []
    raw_hashes = {row["ensaio"]: sha(track(INTER / "entradas" / f"{row['ensaio']}.npz")) for row in records}
    for old in old_pairs:
        a, b = old["run_x"], old["run_y"]
        assert old["assinaturas"]["core"] == task_signatures["core"]
        marker = INTER / "pares" / f"{a}__{b}.json"
        assert json.loads(track(marker).read_text()) == old
        signatures = dict(**task_signatures, input_x=raw_hashes[a], input_y=raw_hashes[b],
                          referencia=sha(marker), pearson=pearson[a, b])
        tasks.append((old, pearson[a, b], signatures))
    outputs = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(calculate_pair, task) for task in tasks]
        for future in as_completed(futures):
            outputs.append(future.result())
            if len(outputs) % 25 == 0 or len(outputs) == 325:
                print(f"WCO completa: {len(outputs)}/325 pares, {time.monotonic()-started:.1f}s", flush=True)
    outputs.sort(key=lambda row: (row["valores"]["ensaio_x"], row["valores"]["ensaio_y"]))
    rows = [row["valores"] for row in outputs]
    assert sum(row["tipo"] == "intra" for row in rows) == 87
    summary = summarize(rows, records, examples)
    print("RESUMOS_PRONTOS " + json.dumps(summary["resumo_intra_inter"], ensure_ascii=False), flush=True)
    from figuras_resumos import generate
    generate()
    after = tree_hashes()
    assert before == after, "As saídas anteriores ou a configuração foram alteradas durante a etapa."
    assert all(sha(Path(__file__) if label == reference(Path(__file__)) else path) == expected
               for path in (Path(__file__), SCRIPTS_ROOT / "wavelet_core.py", SCRIPTS_ROOT / "figuras_resumos.py",
                            SCRIPTS_ROOT / "runtime.py")
               for label, expected in [(reference(path), scripts[reference(path)])])
    errors = {key: max(row["conferencia"][key]["erro_absoluto"] for row in outputs) for key in COMPARE}
    validation = dict(status="PASS", pares=325, todas_escalas_validas=True,
                      erros_maximos_contra_etapa02=errors, coeficientes_pearson_preservados=True,
                      produtos_anteriores_preservados=True, hashes_anteriores=before,
                      limites_wco=[min(row["wco_minima"] for row in rows), max(row["wco_maxima"] for row in rows)],
                      dominio="200 escalas por todos os 2N-1 atrasos, sem redução de colunas",
                      inferencia=False, sensibilidade=False)
    write_json(OUT / "verificacao.json", validation)
    metadata = dict(status="PASS", complete=True, stage="06_resumos",
        generated_utc=datetime.now(timezone.utc).isoformat(),
        inputs=inputs, scripts=scripts,
        counts=dict(ensaios=26, pares=325, intra=87, inter=238, combinacoes=14),
        definitions=dict(
            pearson="Um coeficiente por par, copiado da etapa04 e conferido com sua matriz, sobre sinais completos.",
            wco="Campos de trechos iniciais comuns, cortados antes do preparo, com o núcleo local inalterado.",
            estatisticas_campo="Mínimo, média aritmética, mediana e máximo de todas as células da WCO.",
            dominio="200 comprimentos de onda geométricos de50a3000px, atrasos inteiros de-(N-1)aN-1.",
            pesos_campo="Peso igual por célula na grade logarítmica de comprimentos de onda, sem peso por potência.",
            mascara="Normalização exige ambas as potências no atraso zero acima de1e-10dos respectivos máximos. Todas as200escalas válidas.",
            resumos_grupos="Cada par recebe o mesmo peso. Não se concatenam células de mapas de comprimentos diferentes.",
            quartis="numpy.quantile, method=linear, sobre os valores de cada grupo de pares.",
            leitura="Descritiva, sem inferência ou sensibilidade. Pares compartilham aquisições."),
        runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                     workers=workers, seconds=time.monotonic()-started),
        memory=dict(maior_campo_celulas=max(row["numero_valores"] for row in rows),
                    maior_wco_bytes=max(row["numero_valores"] for row in rows)*8,
                    maior_tres_complexos_e_wco_bytes=max(row["numero_valores"] for row in rows)*56,
                    observacao="Valores estimados dos arrays, não pico medido. Temporários e processos acrescentam memória. Campos completos não são salvos."),
        files={reference(path): sha(path) for path in sorted(OUT.rglob("*"))
               if path.is_file() and path.name != "manifesto.json"})
    write_json(OUT / "manifesto.json", metadata)
    print(f"RESUMOS_COMPLETOS {time.monotonic()-started:.1f}s", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=2)
    main(parser.parse_args().workers)
