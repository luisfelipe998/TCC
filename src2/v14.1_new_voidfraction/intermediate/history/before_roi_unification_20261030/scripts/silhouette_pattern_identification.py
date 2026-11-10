#!/usr/bin/env python3
"""Fold-specific liquid-template silhouette classification comparison.

Rebuild the liquid-like reference from training acquisitions only. Segmentation
never receives the test label; each fold uses identical train/test segmentation.
The silhouette is an assumption-conditional optical reconstruction, not gas
ground truth. Models/features are fixed without parameter search.
"""
from __future__ import annotations

from v14_paths import ROOT as HERE, DATA, CACHE, artifact_path, write_markdown
from v14_geometry import LUMEN_SLICE, LUMEN_ROWS, FIXED_PROBE_ROW, BOUNDARY_PROVENANCE
from v14_selection import selected_folders, selected_runs, PATTERN_ORDER, scope_metadata, is_eligible
import hashlib
import json
import numpy as np
from panorama_void_alternatives import make_reference, REFERENCE_RUNS, silhouette, geometry, read, DIAMETER_PX
from raw_flip_analysis import DATA, HERE, ORDER, labels, state_statistics, write_csv
from pattern_identification import predict, metrics
import matplotlib.pyplot as plt

FEATURES = ["conditional_alpha_visible_mean", "center_gas_duty",
            "mask_pixel_flip_rate53", "center_flip_rate",
            "visible_width_std_over_D", "gas_interval_mean_length_over_D",
            "center_longest_gas_interval_fraction"]


def interval_lengths(state):
    state = np.asarray(state, dtype=bool)
    edges = np.flatnonzero(np.diff(np.r_[False, state, False].astype(np.int8)))
    starts, stops = edges[::2], edges[1::2]
    lengths = stops - starts
    censored = (starts == 0) | (stops == len(state))
    return lengths, censored


def extract(mask):
    core = mask[LUMEN_SLICE]
    center = mask[FIXED_PROBE_ROW]
    lengths, censored = interval_lengths(center)
    widths = mask.sum(axis=0)
    alpha = geometry(mask)
    return {"conditional_alpha_visible_mean": float(alpha.mean()),
            "center_gas_duty": float(center.mean()),
            "mask_pixel_flip_rate53": state_statistics(core, axis=1)["rate"],
            "center_flip_rate": state_statistics(center[:, None])["rate"],
            "visible_width_std_over_D": float(widths.std() / DIAMETER_PX),
            "gas_interval_mean_length_over_D": float(lengths.mean() / DIAMETER_PX) if len(lengths) else 0.,
            "center_longest_gas_interval_fraction": float(lengths.max() / len(center)) if len(lengths) else 0.,
            "n_center_gas_intervals": len(lengths),
            "n_center_gas_intervals_edge_censored": int(censored.sum()),
            "radial_crop_touch_column_fraction": float((mask[LUMEN_ROWS[0]] | mask[LUMEN_ROWS[1] - 1]).mean())}


def checks():
    lengths, censored = interval_lengths([1, 1, 0, 1, 1, 1, 0, 1])
    assert lengths.tolist() == [2, 3, 1]
    assert censored.tolist() == [True, False, True]
    lengths, censored = interval_lengths([1] * 5)
    assert lengths.tolist() == [5] and censored.tolist() == [True]
    lengths, censored = interval_lengths([0] * 5)
    assert len(lengths) == 0 and len(censored) == 0
    full = np.zeros((80, 10), bool); full[LUMEN_SLICE] = True
    result = extract(full)
    assert np.isclose(result["conditional_alpha_visible_mean"], 1.0)
    assert result["center_gas_duty"] == 1
    assert result["center_flip_rate"] == result["mask_pixel_flip_rate53"] == 0
    assert result["n_center_gas_intervals_edge_censored"] == 1


def main():
    checks()
    assert tuple(ORDER) == tuple(PATTERN_ORDER)
    assert set(REFERENCE_RUNS) <= set(selected_runs()), "out-of-scope reference candidate"
    images, runs, patterns, masses = {}, [], [], []
    for folder in selected_folders():
        im = read(folder.name)
        pattern = labels(folder)
        if not is_eligible(folder.name, pattern, float((im[20:60] == 255).mean())):
            continue
        images[folder.name] = im
        runs.append(folder.name); patterns.append(pattern); masses.append(int(folder.name[1:3]))
    y = np.array([ORDER.index(p) for p in patterns])
    masses = np.array(masses)
    cache, references, cached_rows = {}, {}, []

    def matrix(ref_names):
        ref_names = tuple(ref_names)
        if ref_names not in references:
            references[ref_names] = make_reference(ref_names)[0]
        reference = references[ref_names]
        out = []
        for run in runs:
            key = (ref_names, run)
            if key not in cache:
                mask, accepted, _ = silhouette(images[run], reference)
                cache[key] = extract(mask)
                cached_rows.append({"run": run, "pattern": patterns[runs.index(run)],
                                    "reference_training_runs": ";".join(ref_names),
                                    "paired_shadow_acceptance_fraction": float(accepted.mean()), **cache[key]})
            out.append([cache[key][name] for name in FEATURES])
        return np.asarray(out)

    results, prediction_rows, fold_rows = [], [], []
    for validation, mass_out in [("leave_one_run_out", False), ("leave_one_mass_flux_out", True)]:
        model_predictions = {method: np.empty(len(runs), dtype=int) for method in ["centroid", "knn1", "knn3"]}
        for i, run in enumerate(runs):
            train = masses != masses[i] if mass_out else np.arange(len(runs)) != i
            train_names = {runs[j] for j in np.flatnonzero(train)}
            ref_names = tuple(name for name in REFERENCE_RUNS if name in train_names)
            assert ref_names and set(ref_names) <= train_names
            assert run not in ref_names, "held-out reference leakage"
            assert train.any(), "no training acquisitions"
            present = set(y[train])
            present_patterns = ";".join(ORDER[k] for k in range(len(ORDER)) if k in present)
            absent_patterns = ";".join(ORDER[k] for k in range(len(ORDER)) if k not in present)
            x = matrix(ref_names)
            fold_rows.append({"test_run": run, "validation": validation,
                              "n_training_runs": int(train.sum()),
                              "training_present_patterns": present_patterns,
                              "training_absent_patterns": absent_patterns,
                              "test_pattern_has_training_examples": int(y[i]) in present,
                              "reference_training_runs": ";".join(ref_names),
                              "heldout_run_excluded_from_reference": run not in ref_names})
            for method in model_predictions:
                guess = int(predict(x[train], y[train], x[i:i+1], method)[0])
                model_predictions[method][i] = guess
                prediction_rows.append({"run": run, "truth": ORDER[y[i]],
                                        "prediction": ORDER[guess], "correct": bool(guess == y[i]),
                                        "model": method, "validation": validation,
                                        "n_training_runs": int(train.sum()),
                                        "training_present_patterns": present_patterns,
                                        "training_absent_patterns": absent_patterns,
                                        "test_pattern_has_training_examples": int(y[i]) in present,
                                        "reference_training_runs": ";".join(ref_names),
                                        "test_conditional_alpha_visible": float(x[i, 0])})
        for method, pred in model_predictions.items():
            results.append({"model": method, "validation": validation, **metrics(y, pred)})
    write_csv(artifact_path("silhouette_pattern_features.csv"), cached_rows)
    write_csv(artifact_path("silhouette_pattern_predictions.csv"), prediction_rows)
    write_csv(artifact_path("silhouette_pattern_fold_references.csv"), fold_rows)
    write_csv(artifact_path("silhouette_pattern_scores.csv"),
              [{k: v for k, v in r.items() if not isinstance(v, (dict, list))} for r in results])
    source = artifact_path("panorama_void_alternatives.py")
    data = {"n_runs": len(runs), "included_runs": runs, "scope": scope_metadata(),
            "class_counts": dict(zip(ORDER, np.bincount(y, minlength=len(ORDER)).tolist())),
            "features": FEATURES, "lumen_diameter_px": DIAMETER_PX, "boundary_source": BOUNDARY_PROVENANCE, "fixed_reference_candidates": list(REFERENCE_RUNS),
            "unique_fold_templates": len(references), "results": results,
            "segmentation_source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "all_heldout_reference_exclusions_passed": all(r["heldout_run_excluded_from_reference"] for r in fold_rows),
            "folds_with_absent_training_classes": {
                validation: sum(bool(row["training_absent_patterns"]) for row in fold_rows
                                if row["validation"] == validation)
                for validation in ("leave_one_run_out", "leave_one_mass_flux_out")},
            "absent_training_class_policy": "Fit only observed training classes; retain all test acquisitions. An absent class cannot be predicted.",
            "geometry_and_censoring_checks_passed": True}
    (artifact_path("silhouette_pattern_results.json")).write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    review(data, prediction_rows)
    plot(data)
    print(json.dumps(data, indent=2))


def review(data, predictions):
    baseline_data = json.loads((artifact_path("pattern_identification.json")).read_text())
    previous = {(r["model"], r["validation"]): r for r in baseline_data["results"] if r["variant"] == "baseline"}
    better = sum(r["balanced_accuracy"] > previous[(r["model"], r["validation"])]["balanced_accuracy"] for r in data["results"])
    tied = sum(r["balanced_accuracy"] == previous[(r["model"], r["validation"])]["balanced_accuracy"] for r in data["results"])
    worse = len(data["results"]) - better - tied
    centroid_run = next(r for r in data["results"] if r["model"] == "centroid" and r["validation"] == "leave_one_run_out")
    centroid_mass = next(r for r in data["results"] if r["model"] == "centroid" and r["validation"] == "leave_one_mass_flux_out")
    counts = "; ".join(f"{pattern}: {count}" for pattern, count in data["class_counts"].items())
    lines = ["# Silhouette features: internal pattern-identification comparison", "",
             f"This comparison tests whether the panorama silhouette improves identification within the current user-selected scope. It does not validate conditional void values as true gas volume. All {data['n_runs']} selected, acquisitions passing the saturation criterion are evaluated: {counts}. Excluded acquisitions do not enter features, references, training, tests or figures.", "",
             "Seven features are fixed before modeling: conditional visible circular area mean, center candidate-gas duty, candidate-mask pixel flip-rate over rows 18:71, center-row flip-rate, visible radial-width standard deviation /53, mean center gas-interval length /53, and longest connected center gas interval divided by record length. A candidate interval touching either axial record edge is explicitly counted as censored. Interval summaries include those visible censored lengths and are not complete bubble sizes.", "",
             "## Leakage prevention", "",
             f"For every held-out run or mass-flux group, rebuild the liquid-like template from only the fixed reference candidates remaining among that fold's training acquisitions. Feature scaling, class centroids and nearest neighbors likewise use training acquisitions only. The same fold-specific template segments both training and test images. Every fold verifies that the held-out acquisition is absent from its reference. {data['unique_fold_templates']} different templates are used across the folds.", "",
             "Reference candidates were previously chosen from visually reviewed liquid-containing panoramas. The detector and its settings were developed using this dataset. Removing a held-out acquisition from template construction and classifier fitting prevents direct fitting leakage but cannot undo that prior inspection. These scores are internal cross-validation conditional on the current detector, not untouched external validation or an independent liquid-only calibration.", "",
             "Models are the same fixed nearest-centroid, k=1 and k=3 nearest-neighbor methods as the earlier morphology comparison. No hyperparameters or feature subsets are optimized, and every model is reported.", "",
             "| Model | Validation | Earlier morphology accuracy / balanced | Silhouette accuracy / balanced |",
             "|---|---|---:|---:|"]
    for r in data["results"]:
        old = previous[(r["model"], r["validation"])]
        lines.append(f"| {r['model']} | {r['validation']} | {old['accuracy']:.3f} / {old['balanced_accuracy']:.3f} | {r['accuracy']:.3f} / {r['balanced_accuracy']:.3f} |")
    lines += ["", "## Recall and all mistakes", ""]
    for r in data["results"]:
        lines.append(f"**{r['model']}, {r['validation']}**: " + "; ".join(f"{p}: {v:.3f}" for p, v in r["recall_by_pattern"].items()) + ".")
        failures = [p for p in predictions if p["model"] == r["model"] and p["validation"] == r["validation"] and not p["correct"]]
        lines.append("Errors: " + "; ".join(f"{p['run']} ({p['truth']}→{p['prediction']})" for p in failures) + ".")
        lines.append("")
    lines += ["## Interpretation", "",
              f"For the {len(data['results'])} fixed model/validation pairs, silhouette balanced accuracy is higher in {better}, equal in {tied}, and lower in {worse} than the morphology comparison on the same selected acquisitions. Nearest-centroid leave-one-run-out balanced accuracy is {previous[('centroid', 'leave_one_run_out')]['balanced_accuracy']:.3f} for morphology and {centroid_run['balanced_accuracy']:.3f} for silhouettes; leave-one-mass-flux-out is {previous[('centroid', 'leave_one_mass_flux_out')]['balanced_accuracy']:.3f} and {centroid_mass['balanced_accuracy']:.3f}. These observed results do not establish performance for future classifiers or experiments.", "",
              "A score improvement would show that silhouette information helps match these supplied labels under this validation scheme. It would not establish correct absolute void fraction. The user confirmed all 53 rows in 18:71 as the complete tube interior; D=53 and a fully filled lumen gives α=1. Near-unity estimates may still include unresolved liquid films or optical shadows. No extra fluid area beyond these confirmed boundaries is inferred.", "",
              "The candidate-mask flip denominator is 53*(W−1), whereas the legacy dark-pixel baseline uses 40*(W−1). Fixed-row y=40 flips use W−1; legacy center_* field names denote this row, whereas the confirmed lumen midpoint is y=44. Raw switch numerators therefore have different denominators. All features derive from a spatial stitched record, without a conversion to seconds.", "",
              f"Only one dispersed-bubble acquisition remains. Its held-out fold has no bubble training example; absent classes have no fitted centroid and cannot be predicted. Every test remains in the score, including unsupported-class tests. Folds missing at least one training class: leave-one-run-out {data['folds_with_absent_training_classes']['leave_one_run_out']}, leave-one-mass-flux-out {data['folds_with_absent_training_classes']['leave_one_mass_flux_out']}. Per-fold availability is saved. Plug and slug have three acquisitions each, so one error changes their recall by one third. No labels are changed to improve a score.", "",
              "Checks cover exact gas-interval length/censoring examples, the fully filled 53-row lumen with α=1, zero flip-rate for a constant mask, and every fold's test-reference exclusion. Segmentation parameters are inherited unchanged from `panorama_void_alternatives.py`; its source hash is saved.", "",
              "Reproduce with `python3 src2/v14.1_new_voidfraction/scripts/silhouette_pattern_identification.py`. Outputs: `silhouette_pattern_features.csv`, `silhouette_pattern_predictions.csv`, `silhouette_pattern_fold_references.csv`, `silhouette_pattern_scores.csv`, `silhouette_pattern_results.json`, and `silhouette_pattern_confusion_matrices.png`."]
    write_markdown(artifact_path("silhouette_pattern_review.md"), "\n".join(lines) + "\n")


def plot(data):
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    for ax, result in zip(axes.flat, data["results"]):
        conf = np.asarray(result["confusion_true_rows_pred_columns"])
        ax.imshow(conf, cmap="Blues", vmin=0, vmax=10)
        ax.set_xticks(range(len(ORDER)), ["bubbles", "plug", "slug", "churn", "wavy annular"], rotation=35)
        ax.set_yticks(range(len(ORDER)), ["bubbles", "plug", "slug", "churn", "wavy annular"])
        for i in range(len(ORDER)):
            for j in range(len(ORDER)):
                ax.text(j, i, str(conf[i, j]), ha="center", va="center", color="white" if conf[i, j] >= 6 else "black")
        ax.set_title(f"{result['model']} / {'leave mass out' if 'mass' in result['validation'] else 'leave run out'}\nbalanced accuracy={result['balanced_accuracy']:.3f}")
        ax.set_xlabel("predicted"); ax.set_ylabel("provided label")
    fig.text(.5, .015,
             f"Internal cross-validation: n={data['n_runs']} eligible acquisitions; no external validation.\n"
             "The single dispersed-bubble acquisition has no bubble training example in its held-out test.",
             ha="center", va="bottom", fontsize=10)
    fig.tight_layout(rect=[0, .075, 1, 1])
    fig.savefig(artifact_path("silhouette_pattern_confusion_matrices.png"), dpi=180); plt.close(fig)


if __name__ == "__main__":
    main()
