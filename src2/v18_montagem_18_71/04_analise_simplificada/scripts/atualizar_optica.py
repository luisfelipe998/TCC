#!/usr/bin/env python3
"""Atualiza a associação óptica local sem recalcular Pearson, mapas ou figuras."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from runtime import PACKAGE_ROOT, OUTPUT_ROOT, artifact_path, reference

OUT = OUTPUT_ROOT / "04_analise_simplificada"
DATA = OUT / "dados"
OPTICAL = OUTPUT_ROOT / "03_optica"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source(path):
    return {"arquivo": reference(path), "sha256": sha(path)}


def read_json(path):
    return json.loads(path.read_text())


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def write_csv(path, rows, fields):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def input_path(label):
    if label.startswith("package/"):
        path = (PACKAGE_ROOT / label.removeprefix("package/")).resolve()
        if not path.is_relative_to(PACKAGE_ROOT):
            raise ValueError(label)
        return path
    return artifact_path(label)


def main():
    manifest_path = OUT / "manifesto.json"
    manifest = read_json(manifest_path)
    assert manifest["status"] == "complete"
    summary_path = OUT / "resumo_tcc.json"
    summary = read_json(summary_path)
    spectral_summary = {k: v for k, v in summary.items() if k != "discussion_examples"}
    optical_path = OPTICAL / "results/void_fraction/void_fraction_flip_plot_data.csv"
    optical_name = reference(optical_path)
    spectral_inputs = {p: h for p, h in manifest["input_sha256"].items() if p != optical_name}
    assert all(sha(input_path(p)) == h for p, h in spectral_inputs.items())
    mutable = {DATA / "optica_por_ensaio.csv", DATA / "comparacao_optica_pares_ilustrados.csv",
               OUT / "intermediario/atualizacao_optica.json"}
    protected = [p for folder in (DATA, OUT / "intermediario", OUT / "figuras")
                 for p in folder.rglob("*") if p.is_file() and p not in mutable]
    protected_hashes = {reference(p): sha(p) for p in protected}
    with optical_path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == len({r["run"] for r in rows})
    assert {r["run"] for r in rows} == {r["ensaio"] for r in summary["ensaios"]}
    for row in rows:
        assert row["eligible_pure_unsaturated"] in ("True", "False")
        assert int(row["lumen_start_row"]) == 18
        assert int(row["lumen_stop_row_exclusive"]) == 72
        assert row["candidate_mask_rows"] == "18:72"
    eligible = {r["run"]: r for r in rows if r["eligible_pure_unsaturated"] == "True"}
    optical_rows = [{"ensaio": run, "padrao": r["source_void_label"],
                     "alpha": float(r["conditional_circle_alpha_spatial_mean"]),
                     "F1000": float(r["candidate_pixel_flips_per1000px"])}
                    for run, r in eligible.items()]
    examples, omitted = [], []
    for pair in summary["pares_ilustrados"]:
        if pair["tipo"] != "intra" or pair["padrao"] not in ("churn", "wavy annular"):
            continue
        a, b = pair["run_x"], pair["run_y"]
        if a not in eligible or b not in eligible:
            omitted.append({"ensaio_x": a, "ensaio_y": b, "motivo": "ensaio óptico inelegível"})
            continue
        examples.append({"padrao": pair["padrao"], "ensaio_x": a, "ensaio_y": b,
                         "pearson_r": pair["pearson_r"],
                         "alpha_x": float(eligible[a]["conditional_circle_alpha_spatial_mean"]),
                         "alpha_y": float(eligible[b]["conditional_circle_alpha_spatial_mean"]),
                         "F1000_x": float(eligible[a]["candidate_pixel_flips_per1000px"]),
                         "F1000_y": float(eligible[b]["candidate_pixel_flips_per1000px"])})
    write_csv(DATA / "optica_por_ensaio.csv", optical_rows, ["ensaio", "padrao", "alpha", "F1000"])
    write_csv(DATA / "comparacao_optica_pares_ilustrados.csv", examples,
              ["padrao", "ensaio_x", "ensaio_y", "pearson_r", "alpha_x", "alpha_y", "F1000_x", "F1000_y"])
    summary["discussion_examples"] = examples
    write_json(summary_path, summary)
    assert {k: v for k, v in read_json(summary_path).items() if k != "discussion_examples"} == spectral_summary
    assert all(sha(artifact_path(p)) == h for p, h in protected_hashes.items())
    assert all(sha(input_path(p)) == h for p, h in spectral_inputs.items())
    report = dict(status="PASS", script=source(Path(__file__).resolve()),
        scope="Atualização exclusiva da associação óptica na faixa 18:72.", optical_input=source(optical_path),
        eligible_runs=list(eligible), eligible_count=len(eligible),
        excluded_runs=[r["run"] for r in rows if r["run"] not in eligible],
        discussion_examples=examples, omitted_discussion_examples=omitted,
        pearson_recomputed=False, wavelet_recomputed=False, figures_redrawn=False,
        spectral_summary_unchanged=True, protected_files_sha256=protected_hashes,
        spectral_inputs_sha256=spectral_inputs)
    report_path = OUT / "intermediario/atualizacao_optica.json"
    write_json(report_path, report)
    manifest["optical_update_script"] = source(Path(__file__).resolve())
    manifest["optical_update_report"] = source(report_path)
    manifest["input_sha256"][optical_name] = sha(optical_path)
    manifest["spectral_analysis_unchanged_in_optical_update"] = True
    manifest["output_sha256"] = {reference(p): sha(p) for p in sorted(OUT.rglob("*"))
        if p.is_file() and p != manifest_path
        and not {"__pycache__", "historico", "scripts"}.intersection(p.relative_to(OUT).parts)}
    write_json(manifest_path, manifest)
    print(json.dumps({"status": "PASS", "eligible": len(eligible), "examples": examples,
                      "protected_files": len(protected_hashes)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
