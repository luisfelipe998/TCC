#!/usr/bin/env python3
"""Plot saved v14 void and spatial flip results without rerunning segmentation.

The plotted void is a conditional visible circular reconstruction. Flip rates
count neighboring spatial mask changes, rather than temporal event frequencies.
"""
from __future__ import annotations

import ast
from collections import Counter
import csv
import hashlib
import json
import math
import os
from pathlib import Path

from v14_paths import ROOT, CACHE, artifact_path
from v14_selection import PATTERN_ORDER, selected_runs, scope_metadata, is_eligible

os.environ.setdefault("MPLCONFIGDIR", str(CACHE / "matplotlib"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


PATTERNS = PATTERN_ORDER
PURE_PATTERNS = set(PATTERNS)
COLORS = {
    "dispersed bubbles": "#e69f00", "plug": "#0072b2",
    "slug": "#009e73", "churn": "#d55e00", "wavy annular": "#cc79a7",
}
MARKERS = {
    "dispersed bubbles": "o", "plug": "o", "slug": "o", "churn": "o",
    "wavy annular": "o",
}
SATURATION_THRESHOLD = 0.20
DEFAULT_ROWS = 53
APPARENT_DIAMETER = 60
CROP_CEILING = (DEFAULT_ROWS / APPARENT_DIAMETER) ** 2
INPUT_NAMES = (
    "panorama_void_per_run.csv", "phase_mask_flip_per_run.csv",
    "data_inventory.csv", "adversarial_checks.json",
)
OUTPUT_NAMES = (
    "void_fraction_vs_flip_rate.png", "void_fraction_vs_flip_rate.svg",
    "void_fraction_vs_center_flips.png", "void_fraction_vs_center_flips.svg",
    "void_crop_sensitivity.png", "void_crop_sensitivity.svg",
    "void_fraction_flip_plot_data.csv", "summary_chart_checks.json",
)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def csv_rows(name):
    with artifact_path(name).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    runs = [row["run"] for row in rows]
    assert len(runs) == len(set(runs)), f"Duplicate acquisition in {name}"
    return {row["run"]: row for row in rows}


def boolean(value):
    assert value in {"True", "False"}, value
    return value == "True"


def checked_plot_rows():
    """Join exact saved acquisitions and retain their source arithmetic."""
    void = csv_rows("panorama_void_per_run.csv")
    phase = csv_rows("phase_mask_flip_per_run.csv")
    inventory = csv_rows("data_inventory.csv")
    assert void.keys() == phase.keys() == inventory.keys()
    assert set(void) == set(selected_runs()), "Chart tables must match the selected dataset"
    rows, pixel_errors, center_errors = [], [], []
    for run in sorted(void):
        v, p, item = void[run], phase[run], inventory[run]
        pattern = p["pattern"]
        assert pattern == v["label"].replace("_", " ") == item["pattern"], run
        assert pattern in PATTERNS, pattern
        width = int(p["panorama_width"])
        assert width == int(v["panorama_width_px"]) == int(item["panorama_width"])
        n_rows = int(p["candidate_mask_n_rows"])
        assert n_rows == DEFAULT_ROWS and p["candidate_mask_rows"] == "18:71"
        count = int(p["candidate_pixel_flip_count"])
        comparisons = int(p["candidate_pixel_flip_comparisons"])
        assert comparisons == n_rows * (width - 1)
        assert 0 <= count <= comparisons
        pixel_rate = float(p["candidate_pixel_flip_rate"])
        pixel_error = abs(pixel_rate - count / comparisons)
        assert pixel_error < 1e-15, run
        pixel_errors.append(pixel_error)
        center_count = int(p["center_candidate_flip_count"])
        center_comparisons = width - 1
        center_rate = float(p["center_candidate_flip_rate"])
        center_scaled = float(p["center_candidate_flips_per1000px"])
        center_error = abs(center_rate - center_count / center_comparisons)
        assert 0 <= center_count <= center_comparisons
        assert center_error < 1e-15, run
        assert abs(center_scaled - 1000 * center_count / center_comparisons) < 1e-12
        assert center_count == (int(p["center_candidate_gas_entries"])
                                + int(p["center_candidate_gas_exits"]))
        center_errors.append(center_error)
        saturation = float(v["core_saturated_fraction"])
        assert abs(saturation - float(item["core_saturation_fraction"])) < 1e-15
        saturated = saturation > SATURATION_THRESHOLD
        pure = pattern in PURE_PATTERNS
        eligible = boolean(p["eligible"])
        assert eligible == boolean(item["pure_unsaturated"]) == is_eligible(run, pattern, saturation)
        assert pure == boolean(item["pure_pattern"])
        touch = float(v["radial_crop_touch_column_fraction"])
        alpha = float(v["conditional_circle_alpha_spatial_mean"])
        assert 0 <= touch <= 1 and 0 <= alpha <= CROP_CEILING + 1e-12
        assert saturated == (v["validity"] == "unsupported_saturation")
        rows.append({
            "run": run, "supplied_pattern": pattern, "source_void_label": v["label"],
            "eligible_pure_unsaturated": eligible, "pure_pattern": pure,
            "saturated": saturated, "saturation_threshold": SATURATION_THRESHOLD,
            "core_saturated_fraction": v["core_saturated_fraction"],
            "void_validity": v["validity"],
            "radial_crop_touch_column_fraction": v["radial_crop_touch_column_fraction"],
            "radial_crop_touched": touch > 0,
            "all_columns_touch_radial_crop": math.isclose(touch, 1.0, rel_tol=0, abs_tol=1e-12),
            "default_crop_ceiling": CROP_CEILING,
            "conditional_circle_alpha_spatial_mean": v["conditional_circle_alpha_spatial_mean"],
            "conditional_missing_width_alpha_lower": v["conditional_missing_width_alpha_lower"],
            "conditional_missing_width_alpha_upper": v["conditional_missing_width_alpha_upper"],
            "panorama_width_px": width, "candidate_mask_rows": p["candidate_mask_rows"],
            "candidate_mask_n_rows": n_rows,
            "candidate_pixel_flip_count": count,
            "candidate_pixel_flip_comparisons": comparisons,
            "candidate_pixel_flip_rate": p["candidate_pixel_flip_rate"],
            "candidate_pixel_flips_per1000px": 1000 * pixel_rate,
            "center_candidate_flip_count": center_count,
            "center_candidate_flip_comparisons": center_comparisons,
            "center_candidate_flip_rate": p["center_candidate_flip_rate"],
            "center_candidate_flips_per1000px": p["center_candidate_flips_per1000px"],
            "center_candidate_gas_entries": int(p["center_candidate_gas_entries"]),
            "center_candidate_gas_exits": int(p["center_candidate_gas_exits"]),
            "phase_status": p["phase_status"],
        })
    assert set(row["supplied_pattern"] for row in rows) == set(PATTERNS)
    return rows, {
        "joined_acquisitions": len(rows), "unique_supplied_patterns": len(PATTERNS),
        "pattern_counts_all_shown": dict(Counter(r["supplied_pattern"] for r in rows)),
        "eligible_pure_unsaturated": sum(r["eligible_pure_unsaturated"] for r in rows),
        "saturated_acquisitions": [r["run"] for r in rows if r["saturated"]],
        "unsaturated_transition_acquisitions": [r["run"] for r in rows
                                                if not r["pure_pattern"] and not r["saturated"]],
        "any_radial_crop_touch_acquisitions": sum(r["radial_crop_touched"] for r in rows),
        "all_columns_touch_crop_acquisitions": [r["run"] for r in rows
                                                if r["all_columns_touch_radial_crop"]],
        "max_pixel_rate_count_denominator_error": max(pixel_errors),
        "max_center_rate_count_denominator_error": max(center_errors),
        "labels_match_void_phase_and_inventory": True,
        "saturation_eligibility_and_widths_match_inventory": True,
        "entry_plus_exit_equals_center_flip_count": True,
    }


def hull_indices(points):
    """Return the observed convex polygon; reject fewer than three/collinear points."""
    unique = {}
    for index, (x, y) in enumerate(points):
        unique.setdefault((x, y), index)
    ordered = [(x, y, index) for (x, y), index in sorted(unique.items())]
    if len(ordered) < 3:
        return []

    def cross(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    lower, upper = [], []
    for point in ordered:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    for point in reversed(ordered):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    hull = lower[:-1] + upper[:-1]
    return [point[2] for point in hull] if len(hull) >= 3 else []


def save_figure(fig, basename):
    fig.savefig(artifact_path(basename + ".png"), dpi=180, facecolor="white")
    fig.savefig(artifact_path(basename + ".svg"), facecolor="white")
    plt.close(fig)


def alternation_chart(rows, metric, basename, center=False):
    fig, ax = plt.subplots(figsize=(13.7, 8.0))
    fig.subplots_adjust(left=0.09, right=0.71, top=0.83, bottom=0.18)
    all_y = [float(r[metric]) for r in rows]
    ymax = max(all_y) * 1.12
    ax.set_xlim(-0.025, 0.835)
    ax.set_ylim(-0.025 * ymax, ymax)
    ax.set_xticks(np.arange(0, 0.81, 0.1))
    ax.grid(color="#d5d5d5", linestyle="--", linewidth=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    ax.axvspan(CROP_CEILING, 0.835, color="#efefef", zorder=0)
    ax.axvline(CROP_CEILING, color="#666666", linestyle="--", linewidth=1.25)
    ax.text(CROP_CEILING + 0.011, 0.98 * ymax, "Default crop ceiling\n0.7803",
            fontsize=10, color="#555555", ha="left", va="top", rotation=90)
    hull_checks = {}
    for pattern in PATTERNS:
        selected = [r for r in rows if r["supplied_pattern"] == pattern]
        eligible = [r for r in selected if r["eligible_pure_unsaturated"]]
        points = [(float(r["conditional_circle_alpha_spatial_mean"]), float(r[metric]))
                  for r in eligible]
        hull = hull_indices(points) if pattern in PURE_PATTERNS else []
        if hull:
            polygon = np.asarray([points[i] for i in hull])
            ax.fill(polygon[:, 0], polygon[:, 1], facecolor=COLORS[pattern],
                    edgecolor=COLORS[pattern], linewidth=1.4, alpha=0.16, zorder=1)
        hull_checks[pattern] = {
            "eligible_runs": [r["run"] for r in eligible],
            "polygon_vertex_runs": [eligible[i]["run"] for i in hull],
            "shaded": bool(hull),
        }
        for row in selected:
            x = float(row["conditional_circle_alpha_spatial_mean"])
            y = float(row[metric])
            ax.scatter(x, y, marker=MARKERS[pattern], s=83,
                       facecolors="white" if row["saturated"] else COLORS[pattern],
                       edgecolors=COLORS[pattern] if row["saturated"] else "#222222",
                       linewidths=1.8 if row["saturated"] else 0.75, zorder=4)
    by_run = {r["run"]: r for r in rows}
    for row in [r for r in rows if r["supplied_pattern"] == "dispersed bubbles"]:
        run, offset = row["run"], (13, 10)
        ax.annotate(run, (float(row["conditional_circle_alpha_spatial_mean"]),
                          float(row[metric])), xytext=offset, textcoords="offset points",
                    color=COLORS["dispersed bubbles"], fontsize=10, fontweight="bold")
    saturated_note = ["Saturated acquisitions:"] + [
        f"{row['run']} — {row['supplied_pattern']}"
        for row in rows if row['saturated']
    ]
    fig.text(0.735, 0.36, "\n".join(saturated_note), fontsize=9.5,
             color="#555555", linespacing=1.7)
    handles = []
    for pattern in PATTERNS:
        n = sum(r["supplied_pattern"] == pattern for r in rows)
        handles.append(Line2D([], [], color="none", marker=MARKERS[pattern],
                              markerfacecolor=COLORS[pattern], markeredgecolor="#222222",
                              markersize=8, label=f"{pattern} (n={n})"))
    handles.append(Line2D([], [], color="none", marker="o", markerfacecolor="white",
                          markeredgecolor="#777777", markersize=8,
                          label=f"Hollow: saturated image (n={sum(r['saturated'] for r in rows)})"))
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.025, 1.02),
              title="Supplied pattern — all acquisitions", fontsize=10.2,
              title_fontsize=10.2, frameon=False, labelspacing=0.85)
    ax.set_xlabel("Conditional visible void fraction, α [–]", fontsize=13)
    ax.set_ylabel("Center-line changes / 1,000 adjacent positions" if center else
                  "Gas-mask changes / 1,000 adjacent comparisons\n(53-row average)",
                  fontsize=12)
    ax.tick_params(labelsize=11)
    fig.suptitle("Void fraction × center-row phase alternation" if center else
                 "Void fraction × spatial phase alternation", fontsize=18, y=0.96)
    fig.text(0.09, 0.895, f"{len(rows)} retained acquisitions • {len(PATTERNS)} single patterns • spatial image estimates",
             fontsize=11.5, color="#444444")
    footnote = (
        "Shading = observed spread of eligible pure-pattern runs; not classification boundaries.\n"
        "Void estimates depend on gas shape and optical mapping. The dashed line marks the 53-row crop limit.\n"
        + ("Center-row flips can miss off-center bubbles. They describe space, not events per second."
           if center else "Spatial mask changes are not events per second; saturation and radial clipping are separate limits.")
    )
    fig.text(0.09, 0.045, footnote, fontsize=10.2, color="#444444", linespacing=1.55)
    save_figure(fig, basename)
    return hull_checks


def crop_chart(rows, adversarial):
    by_run = {r["run"]: r for r in rows}
    records = adversarial["crop_sensitivity"]
    assert len(records) == len({r["run"] for r in records}) and all(r["run"] in by_run for r in records)
    fig, ax = plt.subplots(figsize=(13.4, 7.3))
    fig.subplots_adjust(left=0.085, right=0.66, top=0.83, bottom=0.2)
    crop_checks = []
    seen_patterns = Counter()
    for index, record in enumerate(records):
        run = record["run"]
        pattern = by_run[run]["supplied_pattern"]
        points = []
        for key, values in record["crop_variants"].items():
            start, stop = ast.literal_eval(key)
            assert isinstance(start, int) and isinstance(stop, int) and 0 <= start < stop
            width = stop - start
            ceiling = (width / APPARENT_DIAMETER) ** 2
            alpha = float(values["conditional_alpha"])
            assert abs(ceiling - values["analytic_crop_ceiling"]) < 1e-12
            assert 0 <= alpha <= ceiling + 1e-12
            if (start, stop) == (18, 71):
                assert abs(alpha - float(by_run[run]["conditional_circle_alpha_spatial_mean"])) < 1e-15
            points.append((width, alpha, (start, stop)))
        points.sort()
        assert [p[0] for p in points] == [50, 53, 57, 60]
        marker = "D" if seen_patterns[pattern] else "o"
        seen_patterns[pattern] += 1
        ax.plot([p[0] for p in points], [p[1] for p in points],
                color=COLORS[pattern], marker=marker, markersize=6, linewidth=1.8,
                label=f"{run} — {pattern}")
        crop_checks.append({"run": run, "supplied_pattern": pattern,
                            "crop_widths_px": [p[0] for p in points],
                            "conditional_alpha_values": [p[1] for p in points],
                            "crop_row_intervals": [list(p[2]) for p in points]})
    widths = np.asarray([50, 53, 57, 60])
    ax.plot(widths, (widths / APPARENT_DIAMETER) ** 2, linestyle="--",
            color="#666666", linewidth=1.4, label="Crop-imposed maximum")
    ax.axvline(DEFAULT_ROWS, color="#666666", linewidth=1.1, linestyle=":")
    ax.text(DEFAULT_ROWS + 0.13, 0.08, "Default crop\n53 rows", color="#555555", fontsize=10)
    ax.set_xticks(widths)
    ax.set_xlim(49.4, 60.6)
    ax.set_ylim(-0.03, 1.04)
    ax.set_xlabel("Radial crop width [image rows / pixels]", fontsize=12)
    ax.set_ylabel("Conditional visible void fraction, α [–]", fontsize=12)
    ax.grid(color="#d5d5d5", linestyle="--", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=11)
    ax.legend(loc="upper left", bbox_to_anchor=(1.025, 1.015), frameon=False,
              title="Saved crop-sensitivity examples", fontsize=9.7, title_fontsize=10.2,
              labelspacing=0.85)
    fig.suptitle("Changing the crop changes the apparent void estimate", fontsize=17, y=0.965)
    fig.text(0.085, 0.895, "Same saved panoramas, reference and threshold; only included rows change.",
             fontsize=11, color="#444444")
    fig.text(0.085, 0.055,
             f"{len(records)} saved examples; sensitivity, not a confidence interval.\n"
             "A wider crop can include wall/reflection shadows; none of these widths is calibrated physical truth.",
             fontsize=10.3, color="#444444", linespacing=1.55)
    save_figure(fig, "void_crop_sensitivity")
    return crop_checks


def main():
    # Fail before reading/writing results if integration has not registered all destinations.
    for name in OUTPUT_NAMES:
        artifact_path(name)
    inputs = {name: sha256(artifact_path(name)) for name in INPUT_NAMES}
    rows, checks = checked_plot_rows()
    adversarial = json.loads(artifact_path("adversarial_checks.json").read_text())
    with artifact_path("void_fraction_flip_plot_data.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with plt.rc_context({"font.family": "DejaVu Sans", "svg.fonttype": "none"}):
        pixel_hulls = alternation_chart(rows, "candidate_pixel_flips_per1000px",
                                       "void_fraction_vs_flip_rate")
        center_hulls = alternation_chart(rows, "center_candidate_flips_per1000px",
                                        "void_fraction_vs_center_flips", center=True)
        crop_checks = crop_chart(rows, adversarial)
    assert inputs == {name: sha256(artifact_path(name)) for name in INPUT_NAMES}
    metadata = {
        "passed": True,
        "source_inputs": [{"artifact": artifact_path(name).relative_to(ROOT).as_posix(),
                           "sha256": digest} for name, digest in inputs.items()],
        "script_sha256": sha256(Path(__file__)),
        "checks": checks,
        "dataset_scope": scope_metadata(),
        "axis_definitions": {
            "x": "Saved conditional_circle_alpha_spatial_mean; spatial mean of separate circular gas-width areas divided by 60-pixel apparent diameter squared.",
            "primary_y": "1000 × saved candidate_pixel_flip_rate = 1000 × candidate_pixel_flip_count / [53 × (panorama_width−1)].",
            "secondary_y": "Saved center_candidate_flips_per1000px = 1000 × center_candidate_flip_count / (panorama_width−1), at image row 40.",
            "sampling": "Neighboring panorama columns; spatial phase alternation, not temporal event frequency.",
            "default_crop": "Image rows 18:71; 53 rows; apparent diameter 60 pixels; conditional ceiling (53/60)^2.",
        },
        "default_crop_ceiling": CROP_CEILING,
        "saturation_threshold": SATURATION_THRESHOLD,
        "shading": "Observed convex hulls of eligible pure-pattern points, at least three non-collinear observations. No singleton-bubble hull or classifier boundaries; user exclusions apply before plotting.",
        "primary_chart_hulls": pixel_hulls,
        "center_chart_hulls": center_hulls,
        "crop_sensitivity": crop_checks,
        "limits": [
            "Visible optical/geometric estimates are not independently calibrated absolute gas-volume fractions.",
            "Only the selected five single patterns are analyzed; saturation-unsupported retained points are omitted from hulls.",
            "Radial crop clipping differs from saturated intensity and can censor annular estimates.",
            "Both axes derive from the same mask and are not independent validation.",
            "Filled gas envelopes can hide liquid fragments; center-row flips can miss off-center bubbles.",
            "Saved crop perturbations are heuristic alternatives, not physical confidence intervals.",
        ],
        "outputs": [artifact_path(name).relative_to(ROOT).as_posix() for name in OUTPUT_NAMES],
        "saved_csv_json_inputs_unchanged_during_chart_generation": True,
    }
    artifact_path("summary_chart_checks.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({"passed": True, "acquisitions": len(rows),
                      "patterns": len(PATTERNS), "crop_examples": len(crop_checks)}, indent=2))


if __name__ == "__main__":
    main()
