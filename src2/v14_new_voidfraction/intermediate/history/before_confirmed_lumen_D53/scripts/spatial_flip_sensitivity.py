#!/usr/bin/env python3
"""Primary panorama-only flip-rate sensitivity, without assigning phases."""
from v14_paths import ROOT as HERE, DATA, CACHE, artifact_path, write_markdown
from v14_selection import iter_runs, is_eligible
from pathlib import Path
import csv
import json
import numpy as np
from PIL import Image
from scipy.ndimage import median_filter
from raw_flip_analysis import DATA, HERE, ORDER, otsu, state_statistics, labels, write_csv, compare_groups


def main():
    rows = []
    for folder in iter_runs():
        pano = np.asarray(Image.open(folder / "__panorama_v11.png").convert("L"))
        core = pano[20:60]
        sat = float((core == 255).mean())
        pattern = labels(folder)
        usable = is_eligible(folder.name, pattern, sat)
        for start, stop in [(10, 70), (20, 60), (25, 55)]:
            roi = pano[start:stop]
            base = otsu(roi.ravel())
            cuts = [("otsu", base), ("otsu_minus10", max(0, base - 10)),
                    ("otsu_plus10", min(255, base + 10)), ("fixed70", 70),
                    ("fixed90", 90), ("fixed110", 110),
                    ("median", float(np.median(roi)))]
            for method, threshold in cuts:
                stats = state_statistics(roi <= threshold, axis=1)
                changes = (roi[:, 1:] <= threshold) != (roi[:, :-1] <= threshold)
                delta = np.abs(np.diff(roi.astype(np.int16), axis=1))
                rows.append({"run": folder.name, "pattern": pattern,
                             "usable_existing_20pct_saturation_rule": usable,
                             "rows": f"{start}:{stop}", "method": method,
                             "threshold": threshold, "flip_count": stats["n_flips"],
                             "n_comparisons": stats["n_comparisons"],
                             "flip_rate": stats["rate"],
                             "flips_per1000px_per_row": stats["rate"] * 1000,
                             "dark_occupancy": stats["dark_occupancy"],
                             "flip_share_intensity_change_le5": float(np.mean(delta[changes] <= 5)) if changes.any() else None})
            filtered = median_filter(roi, size=(3, 3), mode="nearest")
            stats = state_statistics(filtered <= base, axis=1)
            rows.append({"run": folder.name, "pattern": pattern,
                         "usable_existing_20pct_saturation_rule": usable,
                         "rows": f"{start}:{stop}", "method": "otsu_spatial_median3",
                         "threshold": base, "flip_count": stats["n_flips"],
                         "n_comparisons": stats["n_comparisons"],
                         "flip_rate": stats["rate"],
                         "flips_per1000px_per_row": stats["rate"] * 1000,
                         "dark_occupancy": stats["dark_occupancy"]})
    write_csv(artifact_path("flip_rate_spatial_sensitivity.csv"), rows)
    summary = []
    for roi in ["10:70", "20:60", "25:55"]:
        for method in ["otsu", "otsu_minus10", "otsu_plus10", "fixed70", "fixed90", "fixed110", "median", "otsu_spatial_median3"]:
            subset = [r for r in rows if r["rows"] == roi and r["method"] == method]
            comp = compare_groups(subset, "flip_rate")
            summary.append({"roi": roi, "method": method, **comp})
    write_csv(artifact_path("flip_rate_spatial_group_sensitivity.csv"), summary)
    n_positive = sum(r["churn_median"] > r["annular_median"] for r in summary)
    base = next(r for r in summary if r["roi"] == "20:60" and r["method"] == "otsu")
    base_runs = [r for r in rows if r["rows"] == "20:60" and r["method"] == "otsu"]
    raw_audit = json.loads(artifact_path("flip_rate_audit.json").read_text())
    native_comparison = next(c for c in raw_audit["group_comparisons"] if c["metric"] == "raw_5probe_flip_rate")
    bubbles = [r for r in base_runs if r["pattern"] == "dispersed bubbles" and r["usable_existing_20pct_saturation_rule"]]
    report = ["# Optical flip-rate sensitivity: saved panorama pixel boundaries", "",
              f"The current scope contains {len(base_runs)} retained acquisitions and {sum(r["usable_existing_20pct_saturation_rule"] for r in base_runs)} statistical acquisitions passing the saturation criterion. Mixed labels, smooth annular and the user-excluded replicate are omitted before processing.", "",
              "This analysis uses only the saved `src2/data_processed/*/__panorama_v11.png` images. `flip_rate` counts changes in a fixed binary optical threshold mask between neighboring axial pixels, averaged across rows. It is a spatial boundary density and is invariant when the dark/bright labels are swapped. Its value does not certify that those pixels are gas or liquid.", "",
              "`F = number_of_pixel_switches / [H*(W-1)]`. The quantity `1000*F` is mean pixel-state changes per 1000 adjacent axial comparisons per row; the numerator counts all rows, so it is not a count of bubbles or liquid slugs. The v11.3 baseline uses rows 20:60 and a fixed, per-panorama Otsu threshold. There is no Schmitt trigger and no native-time sampling rate in this calculation.", "",
              "| Pattern | n | median F | changes/1000 axial comparisons per row |",
              "|---|---:|---:|---:|"]
    for pattern in ORDER:
        vals = [r["flip_rate"] for r in base_runs if r["pattern"] == pattern and r["usable_existing_20pct_saturation_rule"]]
        median = float(np.median(vals))
        report.append(f"| {pattern} | {len(vals)} | {median:.6f} | {median*1000:.3f} |")
    report += ["", "Mann–Whitney p-values below retain the two-sided asymptotic method with continuity correction and tie correction where ties occur. Saved comparisons record the method, sample sizes and tie flag. These exploratory rank tests compare distributions, not medians alone; no multiple-comparison correction is applied.",
               "", f"Churn versus wavy annular: AUC={base['auc_churn_greater']:.4f}; exploratory two-sided run-level Mann–Whitney p={base['mwu_two_sided_p']:.5f}, n={base['n_churn']} and {base['n_annular']}. The churn>annular median direction holds in {n_positive}/{len(summary)} ROI/threshold/filter variants.",
               "", f"The variants use three radial ROIs (10:70, 20:60, 25:55), Otsu, Otsu±10 gray levels, fixed thresholds 70/90/110, the image median, and Otsu after a 3×3 median filter. None was selected to enforce a regime ordering. Sensitivity values and actual pixel switch counts are recorded for all {len(base_runs)} retained acquisitions. Saturation exclusion is the existing >20% threshold on baseline core pixels and is kept constant across variants.",
               "", f"The saved panorama is stitched from native-frame strips with unequal measured widths and feathering at seams. Spatial density can therefore depend on advection speed, stitching and image noise as well as interface geometry. The supplemental raw-frame audit shows that churn has a higher median native optical flip activity than annular, but substantial run-level overlap; the native five-station comparison has two-sided p={native_comparison['mwu_two_sided_p']:.4f}. A claim that churn is universally more temporally intermittent does not follow from the stronger spatial result.",
               "", f"The retained dispersed-bubble group contains {len(bubbles)} acquisition. Its baseline spatial optical flip rate is {float(np.median([r['flip_rate'] for r in bubbles])):.6f}; one acquisition does not characterize within-pattern variability. The group total is the current unsaturated retained scope, not the earlier unfiltered count.",
               "", "No absolute gas volume fraction, phase velocity or liquid-reference calibration is inferred by this script. The computed spatial flip values are reproducible optical observables. `flip_share_intensity_change_le5` diagnoses small-amplitude boundary crossings but does not identify which boundaries are noise.",
               "", "Run `python3 src2/v14_new_voidfraction/scripts/spatial_flip_sensitivity.py` to regenerate `flip_rate_spatial_sensitivity.csv`, `flip_rate_spatial_group_sensitivity.csv` and this review. Native time checks and exact all-run v11.3 reproduction are documented in `flip_rate_review.md`."]
    write_markdown(artifact_path("spatial_flip_review.md"), "\n".join(report) + "\n")
    print(json.dumps({"n_run_variants": len(rows), "churn_gt_annular_median_variants": n_positive, "total_group_variants": len(summary), "baseline": base}, indent=2))


if __name__ == "__main__":
    main()
