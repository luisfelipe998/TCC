#!/usr/bin/env python3
"""Image-only regime classification feasibility with run-wise validation.

Feature extraction never receives a regime label. Models are intentionally
simple, use no parameter search, and fit scaling on training runs only. Every
selected pure label is retained. The single remaining bubble acquisition cannot
supply a bubble training example when held out. Image crops are sensitivity
checks, never independent replicates.
"""
from __future__ import annotations

from v14_paths import ROOT as HERE, DATA, CACHE, artifact_path, write_markdown
from v14_selection import selected_folders, selected_runs, PATTERN_ORDER, scope_metadata, is_eligible
import json
import re
from pathlib import Path
import numpy as np
from PIL import Image
from scipy.ndimage import binary_closing, binary_fill_holes, label as components
from raw_flip_analysis import DATA, HERE, ORDER, COLORS, otsu, labels, state_statistics, write_csv
import matplotlib.pyplot as plt

FEATURES = ["spatial_flip_rate", "dark_projected_fraction",
            "radial_dark_fraction_std", "center_bright_width_mean",
            "center_bright_width_std", "center_bright_absent_fraction",
            "radial_contrast_median", "radial_contrast_std",
            "radially_flat_column_fraction", "enclosed_hole_fraction",
            "largest_enclosed_hole_fraction", "enclosed_holes_per1000px"]


def features(panorama, start=10, stop=70, offset=0):
    """Predeclared morphology; no label, mass flux, heat flux or raw fps."""
    roi = panorama[start:stop].astype(float)
    height, width = roi.shape
    center = 40 - start
    threshold = np.clip(otsu(roi.ravel()) + offset, 0, 255)
    dark = roi <= threshold
    bright = ~dark
    up = np.cumprod(bright[:center + 1][::-1], axis=0).sum(axis=0)
    down = np.cumprod(bright[center:], axis=0).sum(axis=0)
    widths = (up + down - bright[center]) / height
    p10, p90 = np.percentile(roi, [10, 90], axis=0)
    contrast = (p90 - p10) / np.maximum(p90, 1)
    closed = binary_closing(dark, structure=np.ones((3, 3)))
    holes = binary_fill_holes(closed) & ~closed
    comp, ncomp = components(holes)
    areas = np.bincount(comp.ravel())[1:]
    # Exclude one/two-pixel enclosed gaps from object-density summary.
    substantial = areas[areas >= 9]
    return {"spatial_flip_rate": state_statistics(dark, axis=1)["rate"],
            "dark_projected_fraction": float(dark.mean()),
            "radial_dark_fraction_std": float(dark.mean(axis=0).std()),
            "center_bright_width_mean": float(widths.mean()),
            "center_bright_width_std": float(widths.std()),
            "center_bright_absent_fraction": float((widths == 0).mean()),
            "radial_contrast_median": float(np.median(contrast)),
            "radial_contrast_std": float(contrast.std()),
            "radially_flat_column_fraction": float((contrast < .25).mean()),
            "enclosed_hole_fraction": float(holes.mean()),
            "largest_enclosed_hole_fraction": float(areas.max() / holes.size) if areas.size else 0.,
            "enclosed_holes_per1000px": float(len(substantial) / width * 1000)}


def predict(train_x, train_y, test_x, method):
    # Every scale and centroid uses training data only.
    assert len(train_x) and len(train_x) == len(train_y)
    assert np.isfinite(train_x).all() and np.isfinite(test_x).all()
    present_classes = np.unique(train_y)
    mean = train_x.mean(axis=0)
    sd = train_x.std(axis=0)
    sd = np.where(sd > 1e-9, sd, 1.)
    x = (train_x - mean) / sd
    test = (test_x - mean) / sd
    if method == "centroid":
        # An absent class has no centroid. Fit only observed training classes
        # and map the nearest available centroid back to its actual label.
        centroids = np.stack([x[train_y == k].mean(axis=0) for k in present_classes])
        distances = np.sum((test[:, None, :] - centroids[None, :, :]) ** 2, axis=2)
        return present_classes[np.argmin(distances, axis=1)]
    distances = np.sum((test[:, None, :] - x[None, :, :]) ** 2, axis=2)
    k = int(method[-1])
    nearest = np.argsort(distances, axis=1)[:, :k]
    out = []
    for neighbors in nearest:
        votes = np.bincount(train_y[neighbors], minlength=len(ORDER))
        tied = np.flatnonzero(votes == votes.max())
        # In a tied vote the closest neighbor among tied classes wins.
        out.append(next(int(train_y[j]) for j in neighbors if train_y[j] in tied))
    return np.array(out)


def evaluate(x, y, groups, method="centroid", target_x=None, leave_mass=False):
    target_x = x if target_x is None else target_x
    pred = np.empty(len(y), dtype=int)
    for i in range(len(y)):
        train = groups != groups[i] if leave_mass else np.arange(len(y)) != i
        assert train.any(), "no training acquisitions"
        pred[i] = predict(x[train], y[train], target_x[i:i+1], method)[0]
    return pred


def training_fold_metadata(y, groups, runs, leave_mass=False):
    """Disclose class availability for every held-out test, including singletons."""
    rows = []
    for i, run in enumerate(runs):
        train = groups != groups[i] if leave_mass else np.arange(len(y)) != i
        present = set(y[train])
        rows.append({
            "test_run": run, "n_training_runs": int(train.sum()),
            "training_present_patterns": ";".join(ORDER[k] for k in range(len(ORDER)) if k in present),
            "training_absent_patterns": ";".join(ORDER[k] for k in range(len(ORDER)) if k not in present),
            "test_pattern_has_training_examples": int(y[i]) in present,
        })
    return rows


def metrics(y, pred):
    confusion = np.zeros((len(ORDER), len(ORDER)), dtype=int)
    for true, guess in zip(y, pred):
        confusion[true, guess] += 1
    assert np.all(confusion.sum(axis=1) > 0), "missing test class"
    recall = confusion.diagonal() / confusion.sum(axis=1)
    return {"accuracy": float(np.mean(y == pred)), "balanced_accuracy": float(recall.mean()),
            "recall_by_pattern": dict(zip(ORDER, recall.tolist())),
            "confusion_true_rows_pred_columns": confusion.tolist()}


def main():
    assert tuple(ORDER) == tuple(PATTERN_ORDER)
    images, feature_rows, eligible = {}, [], []
    for folder in selected_folders():
        im = np.asarray(Image.open(folder / "__panorama_v11.png").convert("L"))
        images[folder.name] = im
        pattern = labels(folder)
        sat = float((im[20:60] == 255).mean())
        row = {"run": folder.name, "pattern": pattern, "baseline_saturation_fraction": sat,
               "usable": is_eligible(folder.name, pattern, sat), **features(im)}
        feature_rows.append(row)
        if row["usable"]:
            eligible.append(row)
    write_csv(artifact_path("pattern_morphology_features.csv"), feature_rows)
    runs = [r["run"] for r in eligible]
    x = np.array([[r[k] for k in FEATURES] for r in eligible])
    y = np.array([ORDER.index(r["pattern"]) for r in eligible])
    masses = np.array([int(re.match(r"m(\d+)", run)[1]) for run in runs])
    assert set(images) == set(selected_runs())
    coverage = {
        "leave_one_run_out": training_fold_metadata(y, masses, runs),
        "leave_one_mass_flux_out": training_fold_metadata(y, masses, runs, leave_mass=True),
    }
    coverage_by_run = {validation: {r['test_run']: r for r in entries}
                       for validation, entries in coverage.items()}
    results, prediction_rows = [], []
    base_predictions = {}
    for method in ["centroid", "knn1", "knn3"]:
        for cv, leave_mass in [("leave_one_run_out", False), ("leave_one_mass_flux_out", True)]:
            pred = evaluate(x, y, masses, method, leave_mass=leave_mass)
            results.append({"model": method, "validation": cv, "variant": "baseline", **metrics(y, pred)})
            for run, true, guess in zip(runs, y, pred):
                prediction_rows.append({"run": run, "truth": ORDER[true], "prediction": ORDER[guess],
                                        "correct": bool(true == guess), "model": method,
                                        "validation": cv, "variant": "baseline",
                                        **{k: v for k, v in coverage_by_run[cv][run].items() if k != 'test_run'}})
            if not leave_mass:
                base_predictions[method] = pred
    # Test each unseen run under image perturbations; its training runs retain
    # baseline extraction, preventing a sensitivity variant becoming a refit.
    variants = {"otsu_minus10": dict(offset=-10), "otsu_plus10": dict(offset=10),
                "radial_crop15_65": dict(start=15, stop=65),
                "radial_crop20_60": dict(start=20, stop=60),
                "first_axial_half": {}, "second_axial_half": {}}
    for variant, options in variants.items():
        target_features = []
        for run in runs:
            im = images[run]
            if variant == "first_axial_half":
                im = im[:, :im.shape[1] // 2]
            elif variant == "second_axial_half":
                im = im[:, im.shape[1] // 2:]
            feat = features(im, **options)
            target_features.append([feat[k] for k in FEATURES])
        target_x = np.asarray(target_features)
        for method in ["centroid", "knn1", "knn3"]:
            pred = evaluate(x, y, masses, method, target_x=target_x)
            results.append({"model": method, "validation": "leave_one_run_out", "variant": variant, **metrics(y, pred)})
            for run, true, guess in zip(runs, y, pred):
                prediction_rows.append({"run": run, "truth": ORDER[true], "prediction": ORDER[guess],
                                        "correct": bool(true == guess), "model": method,
                                        "validation": "leave_one_run_out", "variant": variant,
                                        **{k: v for k, v in coverage_by_run['leave_one_run_out'][run].items() if k != 'test_run'}})
            consistent = evaluate(target_x, y, masses, method)
            results.append({"model": method, "validation": "leave_one_run_out_consistent_variant",
                            "variant": variant, **metrics(y, consistent)})
            for run, true, guess in zip(runs, y, consistent):
                prediction_rows.append({"run": run, "truth": ORDER[true], "prediction": ORDER[guess],
                                        "correct": bool(true == guess), "model": method,
                                        "validation": "leave_one_run_out_consistent_variant", "variant": variant,
                                        **{k: v for k, v in coverage_by_run['leave_one_run_out'][run].items() if k != 'test_run'}})
    write_csv(artifact_path("pattern_identification_predictions.csv"), prediction_rows)
    write_csv(artifact_path("pattern_identification_scores.csv"),
              [{k: v for k, v in r.items() if not isinstance(v, (dict, list))} for r in results])
    report_data = {"n_selected_runs": len(feature_rows), "n_eligible_runs": len(y),
                   "included_runs": runs, "scope": scope_metadata(),
                   "classes": ORDER, "features": FEATURES,
                   "class_counts": dict(zip(ORDER, np.bincount(y, minlength=len(ORDER)).tolist())),
                   "validation_training_coverage": coverage,
                   "absent_training_class_policy": "Fit only observed training classes; retain all test acquisitions. An absent class cannot be predicted.",
                   "results": results}
    # Fixed model, fixed features, run-label permutation checks association;
    # this is not proof of calibration or external generalization.
    rng = np.random.default_rng(92377)
    null = []
    for _ in range(499):
        shuffled = rng.permutation(y)
        guessed = evaluate(x, shuffled, masses, "centroid")
        null.append(metrics(shuffled, guessed)["balanced_accuracy"])
    observed = metrics(y, base_predictions["centroid"])["balanced_accuracy"]
    report_data["centroid_run_label_permutation"] = {"n_permutations": 499,
                "observed_balanced_accuracy": observed,
                "p": float((1 + np.count_nonzero(np.asarray(null) >= observed)) / 500),
                "null_balanced_accuracy_mean": float(np.mean(null))}
    (artifact_path("pattern_identification.json")).write_text(json.dumps(report_data, indent=2, allow_nan=False) + "\n")
    write_report(report_data, prediction_rows)
    plot(images, feature_rows, report_data)
    print(json.dumps({"n_runs": len(y), "baseline_results": results[:6],
                      "centroid_permutation": report_data["centroid_run_label_permutation"]}, indent=2))


def write_report(data, predictions):
    counts = "; ".join(f"{pattern}: {count}" for pattern, count in data["class_counts"].items())
    majority = max(data["class_counts"].values()) / data["n_eligible_runs"]
    absent_counts = {validation: sum(bool(r["training_absent_patterns"]) for r in rows)
                     for validation, rows in data["validation_training_coverage"].items()}
    text = ["# Image-only pattern identification: feasibility audit", "",
            "Morphology can help identify regimes without calibrating absolute void fraction. This audit extracts 12 optical/morphology features from each saved panorama. The extraction function receives no labels, mass flux, heat flux or camera timing. Folder labels are used only to train and validate simple classifiers.", "",
            f"The current user-selected scope contains {data['n_selected_runs']} acquisitions. The existing saturation rule leaves {data['n_eligible_runs']} five-class tests: {counts}. Selection is applied before feature extraction and plotting; only selected runs appear in the feature CSV and contact sheet. Saturated selected images are excluded from classifier fitting and scores. Each acquisition is one replicate; axial crops do not increase sample size.", "",
            "Features describe pixel flip density, dark projected area, radial dark-fraction variability, central bright-width mean/variability/absence, radial contrast median/variability/flat-column fraction, and enclosed dark-contour hole area/size/density. They are not validated gas areas. All masks use fixed per-panorama Otsu; contrast-flatness threshold 0.25 and 3×3 closing are fixed across runs.", "",
            "## Validation results", "",
            f"Feature standardization and class centroids/neighbors are fit on training runs only. Models are a nearest class centroid and k-nearest neighbors at fixed k=1 and k=3; there is no parameter or feature search. All tested models are reported. Balanced accuracy is mean per-class recall; five-class uniform-guess chance is 0.20. The majority-class ordinary-accuracy baseline is {majority:.3f}. Missing training classes have no fitted centroid and cannot be predicted; every held-out test is retained in the score.", "",
            "| Model | Validation | Variant | Accuracy | Balanced accuracy |",
            "|---|---|---|---:|---:|"]
    for r in data["results"]:
        text.append(f"| {r['model']} | {r['validation']} | {r['variant']} | {r['accuracy']:.3f} | {r['balanced_accuracy']:.3f} |")
    text += ["", "## Per-pattern recall and mistakes", ""]
    for r in data["results"][:6]:
        text.append(f"**{r['model']}, {r['validation']}**: " + "; ".join(f"{p} {v:.3f}" for p, v in r["recall_by_pattern"].items()) + ".")
        mistakes = [p for p in predictions if p["model"] == r["model"] and p["validation"] == r["validation"] and p["variant"] == "baseline" and not p["correct"]]
        text.append("Errors: " + "; ".join(f"{p['run']} ({p['truth']}→{p['prediction']})" for p in mistakes) + ".")
        text.append("")
    perm = data["centroid_run_label_permutation"]
    text += [f"Fixed-centroid model run-label permutation check (499 permutations): observed balanced accuracy {perm['observed_balanced_accuracy']:.3f}, exploratory p={perm['p']:.3f}; null mean {perm['null_balanced_accuracy_mean']:.3f}. This tests label association within this dataset, not gas calibration or performance on new experiments.", "",
             "## Limits of identification", "",
             f"The dispersed-bubble class now has one selected acquisition. Holding it out leaves no bubble training example, so these data cannot validate a transferable bubble classifier. The absent class is never assigned a fabricated centroid and the unsupported test remains in the denominator. Folds with at least one missing training class: leave-one-run-out {absent_counts['leave_one_run_out']}, leave-one-mass-flux-out {absent_counts['leave_one_mass_flux_out']}. Exact present/absent classes are saved for every fold and prediction. Plug and slug each have three examples, so one error moves their recall by one third.", "",
             "Leave-one-mass-flux-out checks whether signatures transfer to an unseen acquisition condition instead of merely recognizing correlated brightness/advection conditions. Radial/threshold tests perturb the unseen test image without refitting its training images. Half-panorama crops test whether a local segment changes the assigned class. These are adversarial diagnostic checks and are not additional independent validation sets.", "",
             "The consistent-variant rows additionally fit and validate each ROI/threshold/crop definition on all training and held-out acquisitions using the same definition. They separate changed-feature robustness from the harsher test-only perturbation distribution shift. None of these model/variant scores is a selected optimum.", "",
             "A classifier can assign a label to every run even when its error rate is high. Report the full confusion matrices and per-run errors, and do not use its predicted class to manufacture the expected void ordering. Correct image-only absolute void fraction and reliable universal five-pattern identification are separate questions.", "",
             "Reproduce: `python3 src2/v14_new_voidfraction/scripts/pattern_identification.py`. Outputs include all features, per-run predictions, confusion matrices, sensitivity scores, the panorama contact sheet and confusion plots."]
    write_markdown(artifact_path("pattern_identification.md"), "\n".join(text) + "\n")


def plot(images, rows, data):
    sheet_rows = (len(rows) + 3) // 4
    fig, axes = plt.subplots(sheet_rows, 4, figsize=(22, 2.5 * sheet_rows))
    for ax, row in zip(axes.flat, rows):
        im = images[row["run"]]
        center = im.shape[1] // 2
        crop = im[:, max(0, center - 1200):min(im.shape[1], center + 1200)]
        ax.imshow(crop, cmap="gray", vmin=0, vmax=255, aspect="auto")
        ax.set_title(f"{row['run']} / {row['pattern']}" + (" [saturated; outside model]" if not row["usable"] else ""), fontsize=10)
        ax.set_xticks([]); ax.set_yticks([])
    for ax in axes.flat[len(rows):]:
        ax.axis("off")
    fig.suptitle(f"All {len(rows)} selected panoramas: central 2400 axial pixels, original radial rows")
    fig.tight_layout(rect=[0, 0, 1, .97]); fig.savefig(artifact_path("pattern_panorama_contact_sheet.png"), dpi=160); plt.close(fig)
    chosen = data["results"][:6]
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    for ax, result in zip(axes.flat, chosen):
        conf = np.asarray(result["confusion_true_rows_pred_columns"])
        ax.imshow(conf, cmap="Blues", vmin=0, vmax=10)
        ax.set_xticks(range(len(ORDER)), ["bubbles", "plug", "slug", "churn", "wavy annular"], rotation=35)
        ax.set_yticks(range(len(ORDER)), ["bubbles", "plug", "slug", "churn", "wavy annular"])
        for i in range(len(ORDER)):
            for j in range(len(ORDER)):
                ax.text(j, i, str(conf[i, j]), ha="center", va="center", color="white" if conf[i, j] >= 6 else "black")
        ax.set_title(f"{result['model']} / {'leave mass out' if 'mass' in result['validation'] else 'leave run out'}\nbalanced accuracy={result['balanced_accuracy']:.3f}")
        ax.set_xlabel("predicted"); ax.set_ylabel("provided label")
    fig.tight_layout(); fig.savefig(artifact_path("pattern_confusion_matrices.png"), dpi=180); plt.close(fig)


if __name__ == "__main__":
    main()
