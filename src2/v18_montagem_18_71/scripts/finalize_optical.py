#!/usr/bin/env python3
"""Join the v18 spectral examples and optical measures with local provenance."""
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path

from runtime import OUTPUT_ROOT, SCRIPTS_ROOT, artifact_path, reference

OPTICAL = OUTPUT_ROOT / "03_optica"
SPECTRAL = OUTPUT_ROOT / "04_analise_simplificada/dados"
FIGURES = OPTICAL / "figuras"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--verification", type=Path)
    args = parser.parse_args()
    selection_path = SPECTRAL / "selecao_intrapadrao.json"
    matrix_path = SPECTRAL / "matriz_pearson_todos.csv"
    plot_path = OPTICAL / "results/void_fraction/void_fraction_flip_plot_data.csv"
    inputs = [selection_path, matrix_path, plot_path,
              OPTICAL / "results/summary/optical_group_summary.csv",
              FIGURES / "origens.json"]
    source_hashes = {reference(path): sha(path) for path in inputs}
    selection = json.loads(selection_path.read_text())
    matrix = {row["ensaio"]: row for row in read_csv(matrix_path)}
    optical = {row["run"]: row for row in read_csv(plot_path)}
    assert len(matrix) == len(optical) == 26 and set(matrix) == set(optical)
    selected = [row for row in selection if row["padrao_x"] in ("churn", "wavy annular")]
    assert {row["padrao_x"] for row in selected} == {"churn", "wavy annular"} and len(selected) == 2
    examples = []
    for row in selected:
        pattern = row["padrao_x"]
        assert row["padrao_y"] == pattern and row["tipo"] == "intra"
        runs = sorted(run for run in optical if optical[run]["supplied_pattern"] == pattern)
        candidates = [float(matrix[a][b]) for a, b in itertools.combinations(runs, 2)]
        x, y = row["ensaio_x"], row["ensaio_y"]
        value = float(row["pearson_r"])
        assert value == float(matrix[x][y]) == max(candidates)
        examples.append(dict(
            pattern=pattern, run_x=x, run_y=y, pearson_r=value,
            alpha_x=float(optical[x]["conditional_circle_alpha_spatial_mean"]),
            alpha_y=float(optical[y]["conditional_circle_alpha_spatial_mean"]),
            F1000_x=float(optical[x]["candidate_pixel_flips_per1000px"]),
            F1000_y=float(optical[y]["candidate_pixel_flips_per1000px"]),
            eligible_x=optical[x]["plotted"] == "True", eligible_y=optical[y]["plotted"] == "True",
            selection="maximum within-pattern Pearson across all retained spectral acquisitions"))
    examples_path = OPTICAL / "results/summary/discussion_example_pairs.csv"
    with examples_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(examples[0]))
        writer.writeheader()
        writer.writerows(examples)
    runs = {row[field] for row in examples for field in ("run_x", "run_y")}
    individual_path = OPTICAL / "results/summary/discussion_example_runs.csv"
    with individual_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(next(iter(optical.values()))))
        writer.writeheader()
        writer.writerows(optical[run] for run in sorted(runs))
    summary_path = OPTICAL / "resumo_tcc.json"
    summary = json.loads(summary_path.read_text())
    summary["discussion_examples"] = examples
    summary["spectral_example_source"] = dict(selection=reference(selection_path),
                                              pearson_matrix=reference(matrix_path),
                                              maxima_independently_verified=True)
    figures = json.loads((FIGURES / "origens.json").read_text())
    for item in figures["figures"]:
        assert sha(artifact_path(item["file"])) == item["sha256"]
    summary["figure_manifest"] = dict(file=reference(FIGURES / "origens.json"),
                                      sha256=sha(FIGURES / "origens.json"),
                                      verification=figures["verification"])
    if args.verification:
        verification = artifact_path(args.verification)
        report = json.loads(verification.read_text())
        assert report["status"] == "PASS"
        summary["independent_verification"] = dict(file=reference(verification),
                                                    sha256=sha(verification), status="PASS")
        source_hashes[reference(verification)] = sha(verification)
    summary["status"] = "complete"
    overlap = summary["hulls"]["overlap_pairs"]
    summary["interpretation"] = [
        "Maior mediana de fração de vazio: " + summary["descriptive_extrema"]["highest_alpha_median"] + ".",
        "Menor mediana de densidade de transições: " + summary["descriptive_extrema"]["lowest_F1000_median"] + ".",
        "Interseções entre envoltórias observadas: " + "; ".join(" / ".join(pair) for pair in overlap) + ".",
        "A leitura conjunta dos dois pares permaneceu descritiva. Não foi calculada correlação entre grandezas ópticas e espectrais."]
    for path in inputs:
        assert sha(path) == source_hashes[reference(path)]
    write_json(summary_path, summary)
    write_json(OPTICAL / "results/validation/joint_summary_input_hashes.json", source_hashes)
    manifest_path = OPTICAL / "results/validation/artifact_manifest.json"
    files = {str(path.relative_to(OPTICAL)): sha(path) for path in sorted(OPTICAL.rglob("*"))
             if path.is_file() and path != manifest_path and "__pycache__" not in path.parts}
    scripts = [SCRIPTS_ROOT / name for name in ("process_optical.py", "optical_model.py",
               "optical_plot.py", "gerar_figuras_opticas.py", "finalize_optical.py", "runtime.py")]
    write_json(manifest_path, dict(hash_algorithm="sha256", files=files,
                                   scripts={reference(path): sha(path) for path in scripts}))
    print(json.dumps(dict(status="complete", examples=examples, groups=summary["groups"],
                           interpretation=summary["interpretation"]), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
