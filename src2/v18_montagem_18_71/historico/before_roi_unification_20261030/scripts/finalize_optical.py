#!/usr/bin/env python3
"""Link the new spectral examples to optical measures and close their provenance."""
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parents[1]
OPTICAL = ROOT / "03_optica"
SPECTRAL = ROOT / "02_wavelet/04_estatistica/dados"
FIGURES = PROJECT / "TCC/figuras/v18_optica"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--verification", type=Path)
    args = parser.parse_args()
    selection_path = SPECTRAL / "selecao_intrapadrao.json"
    matrix_path = SPECTRAL / "matriz_pearson_todos.csv"
    plot_path = OPTICAL / "results/void_fraction/void_fraction_flip_plot_data.csv"
    previous_pairs_path = PROJECT / "src2/v16_spatial2_18_71/04_estatistica/dados/pares_pearson_wco.csv"
    inputs = [selection_path, matrix_path, plot_path, previous_pairs_path,
              OPTICAL / "results/summary/optical_group_summary.csv",
              OPTICAL / "results/validation/change_from_v14_2.csv", FIGURES / "origens.json"]
    source_hashes = {str(path.relative_to(PROJECT)): sha(path) for path in inputs}
    selection = json.loads(selection_path.read_text())
    matrix = {row["ensaio"]: row for row in read_csv(matrix_path)}
    optical = {row["run"]: row for row in read_csv(plot_path)}
    previous_pairs = read_csv(previous_pairs_path)
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
        old = max((pair for pair in previous_pairs if pair["tipo"] == "intra"
                   and pair["padrao_x"] == pair["padrao_y"] == pattern),
                  key=lambda pair: float(pair["pearson_r"]))
        examples.append(dict(
            pattern=pattern, run_x=x, run_y=y, pearson_r=value,
            alpha_x=float(optical[x]["conditional_circle_alpha_spatial_mean"]),
            alpha_y=float(optical[y]["conditional_circle_alpha_spatial_mean"]),
            F1000_x=float(optical[x]["candidate_pixel_flips_per1000px"]),
            F1000_y=float(optical[y]["candidate_pixel_flips_per1000px"]),
            eligible_x=optical[x]["plotted"] == "True", eligible_y=optical[y]["plotted"] == "True",
            previous_run_x=old["ensaio_x"], previous_run_y=old["ensaio_y"],
            previous_pearson_r=float(old["pearson_r"]),
            example_pair_changed={x, y} != {old["ensaio_x"], old["ensaio_y"]},
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
    summary["spectral_example_source"] = dict(selection=str(selection_path.relative_to(PROJECT)),
                                              pearson_matrix=str(matrix_path.relative_to(PROJECT)),
                                              maxima_independently_verified=True)
    for group in summary["previous_groups"]:
        for key in group:
            if key != "pattern":
                group[key] = int(group[key]) if key == "n" else float(group[key])
    previous = {group["pattern"]: group for group in summary["previous_groups"]}
    summary["group_changes"] = [dict(pattern=group["pattern"],
                                     **{key: group[key]-previous[group["pattern"]][key]
                                        for key in ("n", "alpha_median", "F1000_median",
                                                    "alpha_arithmetic_mean", "F1000_arithmetic_mean")})
                                for group in summary["groups"]]
    figures = json.loads((FIGURES / "origens.json").read_text())
    for item in figures["figures"]:
        assert sha(PROJECT / item["file"]) == item["sha256"]
    summary["figure_manifest"] = dict(file=str((FIGURES / "origens.json").relative_to(PROJECT)),
                                      sha256=sha(FIGURES / "origens.json"),
                                      verification=figures["verification"])
    if args.verification:
        verification = args.verification.resolve()
        report = json.loads(verification.read_text())
        assert report["status"] == "PASS"
        summary["independent_verification"] = dict(file=str(verification.relative_to(PROJECT)),
                                                    sha256=sha(verification), status="PASS")
        source_hashes[str(verification.relative_to(PROJECT))] = sha(verification)
    summary["status"] = "complete"
    summary["changes"]["discussion_example_pairs_changed"] = [row["pattern"] for row in examples
                                                                 if row["example_pair_changed"]]
    assert summary["changes"]["highest_alpha_median"] == "wavy annular"
    assert summary["changes"]["lowest_F1000_median"] == "wavy annular"
    assert summary["hulls"]["overlap_pairs"] == [["plug", "slug"]]
    summary["interpretation"] = [
        "O grupo anular apresentou a maior mediana de fração de vazio e a menor mediana de densidade de transições.",
        "As envoltórias de plug e slug se intersectaram. As de churn e anular permaneceram separadas nesta amostra.",
        "A leitura conjunta dos dois pares permaneceu descritiva. Não foi calculada correlação entre grandezas ópticas e espectrais."]
    for path in inputs:
        assert sha(path) == source_hashes[str(path.relative_to(PROJECT))]
    write_json(summary_path, summary)
    write_json(OPTICAL / "results/validation/joint_summary_input_hashes.json", source_hashes)
    manifest_path = OPTICAL / "results/validation/artifact_manifest.json"
    files = {str(path.relative_to(OPTICAL)): sha(path) for path in sorted(OPTICAL.rglob("*"))
             if path.is_file() and path != manifest_path and "__pycache__" not in path.parts}
    scripts = [ROOT / "scripts/process_optical.py", ROOT / "scripts/optical_model.py", Path(__file__).resolve()]
    write_json(manifest_path, dict(hash_algorithm="sha256", files=files,
                                   scripts={str(path.relative_to(PROJECT)): sha(path) for path in scripts}))
    print(json.dumps(dict(status="complete", examples=examples, groups=summary["groups"],
                           interpretation=summary["interpretation"]), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
