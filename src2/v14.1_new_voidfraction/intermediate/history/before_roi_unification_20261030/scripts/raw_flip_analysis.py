#!/usr/bin/env python3
"""Audit spatial and native-time optical flip rates without assigning phases.

Reads only src2/data_processed; writes through the organized v14 artifact layout. No panorama is
treated as a uniformly sampled time record. Pixel intensity is an optical state,
not a validated gas/liquid segmentation. All 501 native frames are retained.
"""
from __future__ import annotations

from v14_paths import panorama_path, original_panorama_path, BASELINE
from v14_paths import ROOT as HERE, DATA, CACHE, artifact_path, write_markdown
from v14_selection import iter_runs, is_eligible, PATTERNS
import csv
import importlib.util
import json
import os
import re
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.ndimage import median_filter
from scipy.stats import mannwhitneyu, spearmanr

os.environ.setdefault("MPLCONFIGDIR", str(CACHE / "matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ORDER = PATTERNS
COLORS = ["#e69f00", "#0072b2", "#009e73", "#d55e00", "#cc79a7"]


def otsu(values):
    """Same inclusive dark class as v11.3; retain exact arithmetic for audit."""
    hist = np.histogram(values, bins=256, range=(0, 256))[0].astype(float)
    total = hist.sum()
    sum_all = np.dot(np.arange(256), hist)
    wb = sb = 0.0
    best = -1.0
    threshold = 0
    for t in range(256):
        wb += hist[t]
        if wb == 0:
            continue
        wf = total - wb
        if wf == 0:
            break
        sb += t * hist[t]
        between = wb * wf * (sb / wb - (sum_all - sb) / wf) ** 2
        if between > best:
            best, threshold = between, t
    return threshold


def state_statistics(states, axis=0):
    """Pooled optical switches per comparison; finite-sample lag-one null."""
    states = np.moveaxis(np.asarray(states, dtype=bool), axis, 0)
    n = states.shape[0]
    if n < 2:
        raise ValueError("Need at least two ordered samples")
    changes = states[1:] != states[:-1]
    count = int(changes.sum())
    rate = count / changes.size
    occupancy = states.mean(axis=0)
    # Each pixel retains its own occupancy: static radial structure is not
    # confused with temporal switching. n/(n-1) is the exact shuffle correction.
    shuffle_null = float(np.mean(2 * occupancy * (1 - occupancy) * n / (n - 1)))
    return {"n_flips": count, "n_comparisons": int(changes.size),
            "rate": float(rate), "dark_occupancy": float(states.mean()),
            "shuffle_expected_rate": shuffle_null,
            "rate_over_shuffle": float(rate / shuffle_null) if shuffle_null else None}


def labels(folder):
    markers = sorted(folder.glob("__*.txt"))
    return markers[0].stem.lstrip("_").replace("_", " ") if markers else "unknown"


def config(folder):
    cfg = {}
    for line in (folder / "_config.xsv").read_text(encoding="latin-1").splitlines():
        if "=" in line:
            key, val = line.split("=", 1)
            cfg[key.strip()] = val.strip()
    return cfg


def write_csv(path, rows):
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def self_checks():
    constant = np.ones((10, 4), dtype=bool)
    assert state_statistics(constant)["rate"] == 0
    alternating = np.tile(np.arange(10)[:, None] % 2, (1, 4)).astype(bool)
    assert state_statistics(alternating)["rate"] == 1
    rng = np.random.default_rng(1944)
    sample = rng.integers(0, 2, (101, 7)).astype(bool)
    assert state_statistics(sample)["rate"] == state_statistics(~sample)["rate"]
    assert state_statistics(sample)["rate"] == state_statistics(sample[::-1])["rate"]
    assert state_statistics(np.array([[0], [1], [1]]))["rate"] == 0.5
    # Pixel-specific null = zero when every pixel is individually static,
    # even though half the spatial pixels are bright and half are dark.
    static_radial = np.tile([False, True], (10, 1))
    assert state_statistics(static_radial)["shuffle_expected_rate"] == 0


def process(folder, reference):
    files = sorted(folder.glob("ImgA*.jpg"))
    cfg = config(folder)
    indices = np.array([int(re.search(r"(\d+)\.jpg$", p.name)[1]) for p in files])
    assert np.all(np.diff(indices) == int(cfg.get("Step", 1))), folder.name
    assert len(files) == int(cfg["Frames"]), folder.name
    native_rate = float(cfg["Rate"])
    index_step = int(cfg.get("Step", 1))
    fs = native_rate / index_step
    frames = np.stack([np.asarray(Image.open(p).convert("L")) for p in files])
    n, h, width = frames.shape
    cx = width // 2
    probes = np.array([width // 4, 3 * width // 8, cx, 5 * width // 8, 3 * width // 4])
    core = frames[:, 20:60, :]
    threshold = otsu(core.ravel())
    sample = core[:, :, probes]
    result = state_statistics(sample <= threshold)
    centered = state_statistics(core[:, :, cx] <= threshold)
    delta = np.abs(np.diff(sample.astype(np.int16), axis=0))
    changes = (sample[1:] <= threshold) != (sample[:-1] <= threshold)
    small_share = float(np.mean(delta[changes] <= 5)) if changes.any() else None
    pano = np.asarray(Image.open(panorama_path(folder)).convert("L"))
    spatial_core = pano[20:60]
    spatial_threshold = otsu(spatial_core.ravel())
    spatial = state_statistics(spatial_core <= spatial_threshold, axis=1)
    old = reference.temporal_flip_rate(pano)
    assert spatial["rate"] == old["flip_rate"], (folder.name, spatial, old)
    sat = float((spatial_core == 255).mean())
    raw_sat = float((core == 255).mean())
    pattern = labels(folder)
    row = {"run": folder.name, "pattern": pattern, "n_native_frames": n,
           "native_camera_fps": native_rate, "frame_index_step": index_step,
           "effective_native_sample_hz": fs,
           "elapsed_s": (n - 1) / fs, "sample_exposure_window_s": n / fs,
           "raw_otsu_threshold": threshold, "spatial_otsu_threshold": spatial_threshold,
           "spatial_pixels": pano.shape[1],
           "spatial_flip_count": spatial["n_flips"],
           "spatial_n_comparisons": spatial["n_comparisons"],
           "spatial_flip_rate": spatial["rate"],
           "spatial_shuffle_expected_rate": spatial["shuffle_expected_rate"],
           "spatial_rate_over_shuffle": spatial["rate_over_shuffle"],
           "raw_x508_flip_count": centered["n_flips"],
           "raw_x508_n_comparisons": centered["n_comparisons"],
           "raw_x508_flip_rate": centered["rate"],
           "raw_x508_flip_hz_per_pixel": centered["rate"] * fs,
           "raw_5probe_flip_count": result["n_flips"],
           "raw_5probe_n_comparisons": result["n_comparisons"],
           "raw_5probe_flip_rate": result["rate"],
           "raw_5probe_flip_hz_per_pixel": result["rate"] * fs,
           "raw_dark_occupancy": result["dark_occupancy"],
           "raw_shuffle_expected_rate": result["shuffle_expected_rate"],
           "raw_rate_over_shuffle": result["rate_over_shuffle"],
           "raw_flip_share_intensity_change_le5": small_share,
           "raw_core_saturation_fraction": raw_sat,
           "spatial_core_saturation_fraction": sat,
           "pure_pattern": pattern in ORDER,
           "usable_existing_20pct_saturation_rule": is_eligible(folder.name, pattern, sat)}
    station_rates = [state_statistics(sample[:, :, j] <= threshold)["rate"] for j in range(5)]
    row["raw_station_rate_min"] = min(station_rates)
    row["raw_station_rate_max"] = max(station_rates)
    variants = []
    for start, stop in [(10, 70), (20, 60), (25, 55)]:
        roi = frames[:, start:stop, :]
        values = roi[:, :, probes]
        base = otsu(roi.ravel())
        med = float(np.median(values))
        for method, thresh in [("otsu", base), ("otsu_minus10", max(0, base - 10)),
                               ("otsu_plus10", min(255, base + 10)),
                               ("fixed70", 70), ("fixed90", 90),
                               ("fixed110", 110), ("median", med)]:
            stats = state_statistics(values <= thresh)
            variants.append({"run": folder.name, "pattern": pattern,
                             "usable_existing_20pct_saturation_rule": row["usable_existing_20pct_saturation_rule"],
                             "rows": f"{start}:{stop}", "method": method,
                             "threshold": thresh, "flip_rate": stats["rate"],
                             "flip_hz_per_pixel": stats["rate"] * fs,
                             "dark_occupancy": stats["dark_occupancy"]})
    filtered = median_filter(core, size=(1, 3, 3), mode="nearest")
    row["raw_spatial_median3_flip_rate"] = state_statistics(filtered[:, :, probes] <= threshold)["rate"]
    # A real chronological permutation preserves marginal optical occupancy
    # while destroying lag-one temporal coherence, including static wall bands.
    perm = np.random.default_rng(4721).permutation(n)
    row["raw_shuffled_flip_rate"] = state_statistics(sample[perm] <= threshold)["rate"]
    # Retain a single-station native kymograph for review, without phase labels.
    diag = {"core": core[:, :, cx].T, "fps": fs, "threshold": threshold,
            "switches_per_frame": ((core[1:, :, cx] <= threshold) != (core[:-1, :, cx] <= threshold)).mean(axis=1)}
    return row, variants, diag


def two_sided_mann_whitney(first, second, method="auto"):
    """Resolve SciPy's small-sample rule and save the method actually used.

    Exact ranks require no pooled ties. The explicit asymptotic option is
    retained for the legacy optical comparisons so their method is unchanged.
    """
    first = np.asarray(first, dtype=float)
    second = np.asarray(second, dtype=float)
    if not first.size or not second.size or not np.isfinite(first).all() or not np.isfinite(second).all():
        raise ValueError("Mann–Whitney samples must be nonempty and finite")
    combined = np.concatenate([first, second])
    ties = np.unique(combined).size != combined.size
    requested = method
    if method == "auto":
        method = "exact" if min(first.size, second.size) <= 8 and not ties else "asymptotic"
    if method not in {"exact", "asymptotic"} or (method == "exact" and ties):
        raise ValueError("Exact Mann–Whitney is used only without pooled ties")
    continuity = method == "asymptotic"
    result = mannwhitneyu(first, second, alternative="two-sided",
                          method=method, use_continuity=continuity)
    metadata = {"mwu_requested_method": requested, "mwu_method": method,
                "mwu_has_pooled_ties": bool(ties),
                "mwu_use_continuity": continuity,
                "mwu_tie_correction_applied": bool(ties and continuity)}
    return result, metadata


def compare_groups(rows, key):
    c = np.array([r[key] for r in rows if r["pattern"] == "churn" and r["usable_existing_20pct_saturation_rule"]])
    a = np.array([r[key] for r in rows if r["pattern"] == "wavy annular" and r["usable_existing_20pct_saturation_rule"]])
    auc = float(np.mean((c[:, None] > a[None, :]) + .5 * (c[:, None] == a[None, :])))
    u, metadata = two_sided_mann_whitney(c, a, method="asymptotic")
    return {"metric": key, "n_churn": len(c), "n_annular": len(a),
            "churn_median": float(np.median(c)), "annular_median": float(np.median(a)),
            "auc_churn_greater": auc, "mwu_two_sided_p": float(u.pvalue),
            **metadata}


def main():
    self_checks()
    spec = importlib.util.spec_from_file_location("old_void_reference", HERE.parent / "v11.3" / "void_common_v11_3.py")
    reference = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reference)
    rows, variants, diags = [], [], {}
    for folder in iter_runs():
        row, sweep, diagnostic = process(folder, reference)
        rows.append(row)
        variants.extend(sweep)
        diags[folder.name] = diagnostic
        print(f"{folder.name} {row['pattern']}: spatial={row['spatial_flip_rate']:.6f}, raw={row['raw_5probe_flip_rate']:.6f}, Hz={row['raw_5probe_flip_hz_per_pixel']:.3f}", flush=True)
    write_csv(artifact_path("flip_rate_per_run.csv"), rows)
    write_csv(artifact_path("flip_rate_sensitivity.csv"), variants)
    keys = ["spatial_flip_rate", "raw_x508_flip_rate", "raw_x508_flip_hz_per_pixel",
            "raw_5probe_flip_rate", "raw_5probe_flip_hz_per_pixel",
            "raw_spatial_median3_flip_rate", "raw_rate_over_shuffle"]
    summary = []
    for pattern in ORDER:
        subset = [r for r in rows if r["pattern"] == pattern and r["usable_existing_20pct_saturation_rule"]]
        for key in keys:
            vals = [r[key] for r in subset if r[key] is not None]
            summary.append({"pattern": pattern, "metric": key, "n": len(vals),
                            "median": float(np.median(vals)), "min": min(vals), "max": max(vals),
                            "mean": float(np.mean(vals))})
    write_csv(artifact_path("flip_rate_by_pattern.csv"), summary)
    comparisons = [compare_groups(rows, key) for key in keys]
    sweep_results = []
    for roi in ["10:70", "20:60", "25:55"]:
        for method in ["otsu", "otsu_minus10", "otsu_plus10", "fixed70", "fixed90", "fixed110", "median"]:
            subset = [v for v in variants if v["rows"] == roi and v["method"] == method]
            comp = compare_groups(subset, "flip_hz_per_pixel")
            sweep_results.append({"roi": roi, "threshold_method": method, **comp})
    write_csv(artifact_path("flip_rate_group_sensitivity.csv"), sweep_results)
    usable = [r for r in rows if r["usable_existing_20pct_saturation_rule"]]
    rho, p = spearmanr([r["spatial_flip_rate"] for r in usable], [r["raw_5probe_flip_hz_per_pixel"] for r in usable])
    audit = {"n_acquisition_runs": len(rows),
             "n_pure_existing_saturation_rule": len(usable),
             "exact_v11_3_flip_reproductions": len(rows),
             "fps_counts": {str(int(fps)): sum(r["native_camera_fps"] == fps for r in rows) for fps in sorted({r["native_camera_fps"] for r in rows})},
             "elapsed_range_s": [min(r["elapsed_s"] for r in rows), max(r["elapsed_s"] for r in rows)],
             "probe_columns": [254, 381, 508, 635, 762],
             "group_comparisons": comparisons,
             "sensitivity_hz_direction_churn_gt_annular": sum(v["churn_median"] > v["annular_median"] for v in sweep_results),
             "sensitivity_variants": len(sweep_results),
             "spatial_vs_raw_hz_spearman": {"rho": float(rho), "p": float(p)},
             "invariants_passed": ["constant", "alternating", "complement", "time_reversal", "N_minus_1_denominator", "pixel_specific_shuffle_null", "native_filename_step", "camera_config_frame_count", "all_run_exact_v11_3_reproduction"]}
    (artifact_path("flip_rate_audit.json")).write_text(json.dumps(audit, indent=2) + "\n")
    plot_results(rows, diags)
    write_review(rows, summary, audit, sweep_results)
    print(json.dumps(audit, indent=2))


def plot_results(rows, diags):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for g, (pattern, color) in enumerate(zip(ORDER, COLORS)):
        subset = [r for r in rows if r["pattern"] == pattern and r["usable_existing_20pct_saturation_rule"]]
        jit = np.linspace(-.12, .12, len(subset)) if len(subset) > 1 else np.zeros(len(subset))
        for ax, key in zip(axes, ["spatial_flip_rate", "raw_5probe_flip_rate", "raw_5probe_flip_hz_per_pixel"]):
            vals = [r[key] for r in subset]
            ax.scatter(g + jit, vals, color=color, s=42, edgecolor="k", linewidth=.4)
            ax.plot([g - .2, g + .2], [np.median(vals)] * 2, color="k", lw=2)
            ax.set_xticks(range(len(ORDER)), ["bubbles", "plug", "slug", "churn", "wavy annular"], rotation=30)
            ax.grid(axis="y", alpha=.2)
    axes[0].set(title="Panorama: spatial pixel changes", ylabel="switches / adjacent spatial comparison")
    axes[1].set(title="Native frames: five fixed stations", ylabel="switches / adjacent frame comparison")
    axes[2].set(title="Native frames: camera timing applied", ylabel="optical threshold switches / pixel / second")
    fig.suptitle("Optical flip rates; phase segmentation is uncalibrated")
    fig.tight_layout()
    fig.savefig(artifact_path("flip_rate_spatial_vs_native.png"), dpi=180)
    plt.close(fig)
    selected = ["m07q01_001", "m06q01_001", "m05q05_001", "m06q20_001", "m06q80_001"]
    fig, axes = plt.subplots(len(selected), 2, figsize=(13, 2 * len(selected) + 3))
    for i, run in enumerate(selected):
        d = diags[run]
        stop_ms = (d["core"].shape[1] - 1) / d["fps"] * 1000
        axes[i, 0].imshow(d["core"], cmap="gray", vmin=0, vmax=255, aspect="auto", extent=[0, stop_ms, 60, 20])
        axes[i, 0].set(title=f"{run}: raw x=508, rows 20:60", ylabel="camera row")
        t_ms = np.arange(1, d["core"].shape[1]) / d["fps"] * 1000
        axes[i, 1].plot(t_ms, d["switches_per_frame"], lw=.8)
        axes[i, 1].set(title=f"Optical switching at fixed x; threshold={d['threshold']}", ylim=(0, 1), ylabel="fraction of rows that switch")
        axes[i, 1].grid(alpha=.2)
    for ax in axes[-1]:
        ax.set_xlabel("elapsed native camera time [ms]")
    fig.tight_layout()
    fig.savefig(artifact_path("flip_rate_native_kymographs.png"), dpi=170)
    plt.close(fig)


def write_review(rows, summary, audit, sweep):
    lookup = {(s["pattern"], s["metric"]): s for s in summary}
    group_counts = ", ".join(f"{p} {lookup[(p, 'spatial_flip_rate')]['n']}" for p in ORDER)
    lines = ["# Flip-rate audit: spatial record versus native camera time", "",
             "All calculations use `src2/data_processed` only. A flip is a change in an optical threshold state. These states have not been independently validated as gas and liquid.", "",
             "## Exact definitions", "",
             f"The v11.3 quantity is `sum(abs(diff(binary_mask, axis=1))) / [40*(panorama_width-1)]`. It is a spatial, per-pixel adjacent-column change fraction. It does not use a Schmitt trigger, and its values near 0.02–0.04 are not counts per 1000 px. All {len(rows)} retained panoramas reproduce the existing v11.3 function exactly.", "",
             "The native-time alternative uses all ordered acquisition frames, camera column x=508 and, separately, x={254,381,508,635,762}, rows 20:60. A fixed Otsu threshold is calculated from all raw core pixels within a run. `r_t = flips/[40*C*(N-1)]`; `r_t*fps` is optical threshold switches per pixel per second, not bubble-passage frequency. No image registration or panorama stitching enters this estimate.", "",
             "Camera `_config.xsv` files provide `Rate`, `Frames` and `Step`; filenames are checked against the recorded step and frame count. All runs contain 501 consecutive frames. Rates are 4000, 5000 or 7000 fps, and elapsed first-to-last-frame duration is 0.125, 0.100 or 0.0714286 s. The zero-based time interval has N−1 comparisons.", "",
             "## Run-level group medians", "",
             f"The current user scope retains {len(rows)} acquisitions; the existing >20% panorama-core saturation rule leaves {audit['n_pure_existing_saturation_rule']} statistical acquisitions: {group_counts}. Mixed labels, smooth annular and the excluded replicate do not enter any per-run computation or grouped statistic. Supplied retained labels are unchanged; a single bubble run does not establish within-pattern variability.", "",
             "| Pattern | n | spatial rate | native x508 rate | native 5-station rate | native optical switches/pixel/s |",
             "|---|---:|---:|---:|---:|---:|"]
    for p in ORDER:
        values = [lookup[(p, k)]["median"] for k in ["spatial_flip_rate", "raw_x508_flip_rate", "raw_5probe_flip_rate", "raw_5probe_flip_hz_per_pixel"]]
        lines.append(f"| {p} | {lookup[(p, 'spatial_flip_rate')]['n']} | {values[0]:.6f} | {values[1]:.6f} | {values[2]:.6f} | {values[3]:.3f} |")
    lines += ["", "## Interpretation and sensitivity", "",
              "Mann–Whitney p-values below retain the two-sided asymptotic method with continuity correction and tie correction where ties occur. Each comparison saves its method, sample sizes and tie flag. This rank test compares distributions; it is not generally a test of medians alone.", ""]
    for comp in audit["group_comparisons"]:
        lines.append(f"- `{comp['metric']}`: churn median {comp['churn_median']:.6g}, annular {comp['annular_median']:.6g}; AUC(churn>annular)={comp['auc_churn_greater']:.3f}, exploratory two-sided run-level Mann–Whitney p={comp['mwu_two_sided_p']:.4g} (n={comp['n_churn']} versus {comp['n_annular']}).")
    lines += ["", f"The churn>annular median direction in native optical switches/pixel/s survives {audit['sensitivity_hz_direction_churn_gt_annular']}/{audit['sensitivity_variants']} fixed ROI/threshold variants (ROIs 10:70, 20:60, 25:55; Otsu, Otsu±10, fixed 70/90/110, median). This sweep evaluates optical activity, not true phase intermittency. Group p-values are exploratory and not corrected for the many alternatives.",
              "", f"Spatial versus native-time-Hz run-level rank correlation: Spearman rho={audit['spatial_vs_raw_hz_spearman']['rho']:.3f}, p={audit['spatial_vs_raw_hz_spearman']['p']:.4g}. Different flow advection speeds and different frame rates make their numerical values non-interchangeable.",
              "", "`raw_flip_share_intensity_change_le5` reports the share of switches caused by changes of at most five gray levels. Such switches may include noise, JPEG effects and weak optical fluctuations; the value diagnoses threshold fragility but does not prove that any individual switch is noise. A 3×3 spatial median filter is an explicit nuisance-control variant.",
              "", "Time shuffling does not make flip-rate disappear. It usually increases changes by destroying temporal correlation. The correct pixel-specific shuffle null is the mean of `2*p_pixel*(1-p_pixel)*N/(N-1)`; the old global-alpha Bernoulli null confounds a static radial profile with temporal occupancy. A complement of the binary mask preserves flip-rate exactly, so flipping gas/bright polarity cannot fix this descriptor.",
              "", "The observation window is short. Native pixel activity is directly measurable, but long-term stationarity and low-frequency passage statistics are not established. Camera exposure, built-in filtering, intensity clipping and JPEG compression remain measurement limits; no phase-velocity estimate or assumed Taylor conversion is used here.",
              "", "## Reproduce", "",
              "`python3 src2/v14.1_new_voidfraction/scripts/raw_flip_analysis.py` regenerates the per-run values, ROI/threshold sensitivity, group statistics and diagnostic plots. Invariants cover constant and alternating records, complement and time-reversal invariance, N−1 normalization, pixel-specific shuffle null, actual filename/config consistency, and exact all-run reproduction of v11.3.",
              "", "Outputs: `flip_rate_per_run.csv`, `flip_rate_by_pattern.csv`, `flip_rate_sensitivity.csv`, `flip_rate_group_sensitivity.csv`, `flip_rate_audit.json`, `flip_rate_spatial_vs_native.png`, `flip_rate_native_kymographs.png`."]
    write_markdown(artifact_path("flip_rate_review.md"), "\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
