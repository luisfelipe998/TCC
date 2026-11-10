#!/usr/bin/env python3
"""Recompute optical references, masks and measures after panorama reconstruction on 54 rows."""
from __future__ import annotations

import csv
import hashlib
import itertools
import json
from pathlib import Path
import platform
import sys

import numpy as np
import PIL
from PIL import Image
import scipy
from scipy.optimize import linprog
from scipy.spatial import ConvexHull

from optical_model import (DIAMETER_PX, FUNCTION_NAMES, LUMEN_HEIGHT_PX, LUMEN_ROWS,
                           REFERENCE_RUNS, SOURCE_SHA256, THRESHOLD, load_model)

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parents[1]
SOURCE = ROOT.parent / "v14.1_new_voidfraction"
BASELINE = ROOT.parent / "v14.2_54_linhas"
OUTROOT = ROOT / "03_optica"
PANORAMA_MANIFEST = ROOT / "01_panoramas/manifesto.json"
DATA = ROOT.parent / "data_processed"
PATTERNS = ("dispersed bubbles", "plug", "slug", "churn", "wavy annular")
BOUNDARY_SOURCE = "User correction on 2026-10-29: rows 18 through 71 inclusive, Python slice 18:72, 54 pixel cells, D=54 px."


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path):
    return path.relative_to(PROJECT).as_posix()


def output(name):
    path = OUTROOT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def write_json(name, value):
    output(name).write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def write_csv(name, rows):
    with output(name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def independent_geometry(mask):
    """Detect each contiguous vertical interval independently, column by column."""
    values = []
    for column in mask[18:72].T:
        transitions = np.flatnonzero(np.diff(np.r_[False, column, False]))
        widths = transitions[1::2] - transitions[::2]
        values.append(sum(int(width) ** 2 for width in widths) / (54 ** 2))
    return np.asarray(values)


def hull_checks(plot_rows):
    polygons, halfspaces, details = {}, {}, {}
    for pattern in PATTERNS:
        rows = [row for row in plot_rows if row["plotted"] and row["supplied_pattern"] == pattern]
        points = np.array([[r["conditional_circle_alpha_spatial_mean"],
                            r["candidate_pixel_flips_per1000px"]] for r in rows])
        vertices = ConvexHull(points).vertices if len(points) >= 3 else []
        details[pattern] = dict(eligible_runs=[r["run"] for r in rows],
                                polygon_vertex_runs=[rows[i]["run"] for i in vertices],
                                shaded=len(vertices) >= 3)
        if len(vertices) >= 3:
            polygons[pattern] = points[vertices]
            halfspaces[pattern] = ConvexHull(points).equations

    def cross(a, b, c):
        return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])

    def inside(point, polygon):
        signs = [cross(a, b, point) for a, b in zip(polygon, np.roll(polygon, -1, axis=0))]
        return min(signs) >= -1e-12 or max(signs) <= 1e-12

    def crosses(a, b, c, d):
        if any(max(min(a[k], b[k]), min(c[k], d[k])) >
               min(max(a[k], b[k]), max(c[k], d[k])) + 1e-12 for k in (0, 1)):
            return False
        return cross(a, b, c)*cross(a, b, d) <= 1e-12 and cross(c, d, a)*cross(c, d, b) <= 1e-12

    intersections = []
    for left, right in itertools.combinations(polygons, 2):
        a, b = polygons[left], polygons[right]
        geometric = any(inside(p, b) for p in a) or any(inside(p, a) for p in b)
        geometric |= any(crosses(p, q, r, s) for p, q in zip(a, np.roll(a, -1, axis=0))
                         for r, s in zip(b, np.roll(b, -1, axis=0)))
        constraints = np.vstack([halfspaces[left], halfspaces[right]])
        feasible = linprog([0., 0.], A_ub=constraints[:, :2], b_ub=-constraints[:, 2],
                           bounds=[(None, None)]*2, method="highs")
        assert feasible.success or feasible.status == 2, feasible.message
        assert bool(geometric) == feasible.success, (left, right)
        if geometric:
            intersections.append([left, right])
    return dict(primary_chart_hulls=details, overlap_pairs=intersections,
                churn_annular_intersect=["churn", "wavy annular"] in intersections,
                independent_polygon_and_halfspace_tests_agree=True)


def main():
    reference_path = BASELINE / "intermediate/features/panorama_liquid_reference.json"
    old_void_path = BASELINE / "results/void_fraction/panorama_void_per_run.csv"
    old_phase_path = BASELINE / "results/flip_rate/phase_mask_flip_per_run.csv"
    reconstruction = json.loads(PANORAMA_MANIFEST.read_text())
    assert reconstruction["complete"] is True
    records = {record["ensaio"]: record for record in reconstruction["records"]}
    runs = sorted(records)
    panorama_paths = {run: PROJECT / records[run]["optica_path"] for run in runs}
    assert len(reconstruction["records"]) == 26
    for run, path in panorama_paths.items():
        assert sha256(path) == records[run]["hashes"]["optica"], run
    assert len(runs) == len(set(runs)) == 26
    old_void = {row["run"]: row for row in read_csv(old_void_path)}
    old_phase = {row["run"]: row for row in read_csv(old_phase_path)}
    assert set(runs) == set(old_void) == set(old_phase)
    inputs = [PANORAMA_MANIFEST, reference_path, old_void_path, old_phase_path,
              BASELINE / "results/summary/optical_group_summary.csv",
              BASELINE / "results/validation/hull_intersections.json",
              SOURCE / "scripts/panorama_void_alternatives.py",
              SOURCE / "scripts/v14_geometry.py", SOURCE / "scripts/v14_selection.py"]
    markers = {}
    for run in runs:
        candidates = list((DATA / run).glob("__*.txt"))
        assert len(candidates) == 1, run
        markers[run] = candidates[0]
        inputs.extend([candidates[0], panorama_paths[run]])
    before_hashes = {relative(path): sha256(path) for path in inputs}

    model = load_model(SOURCE, panorama_paths)
    reference, selections = model.make_reference()
    historical_reference = json.loads(reference_path.read_text())
    reference_comparison = dict(
        recomputed_from_new_panoramas=True,
        algorithm_and_three_plug_runs_unchanged=True,
        exactly_equal_to_previous=bool(np.array_equal(reference, historical_reference["radial_reference"])),
        max_absolute_radial_reference_difference=float(np.max(abs(reference-np.asarray(historical_reference["radial_reference"])))),
        old_selected_columns={run: historical_reference["liquid_candidates"][run]["n_selected"] for run in REFERENCE_RUNS},
        new_selected_columns={run: selections[run]["n_selected"] for run in REFERENCE_RUNS})
    write_json("intermediate/features/panorama_liquid_reference.json", dict(
        reference_runs=REFERENCE_RUNS, liquid_candidates=selections, radial_reference=reference.tolist(),
        recomputed_from_new_panoramas=True, comparison=reference_comparison, lumen_rows=LUMEN_ROWS, lumen_diameter_px=DIAMETER_PX,
        boundary_source=BOUNDARY_SOURCE))

    void_rows, phase_rows, plot_rows, comparisons, checks = [], [], [], [], []
    for run in runs:
        im = model.read(run)
        label = markers[run].stem[2:]
        pattern = label.replace("_", " ")
        assert pattern in PATTERNS and pattern == old_phase[run]["pattern"]
        assert label == old_void[run]["label"]
        record_pattern = records[run]["padrao"].replace("_", " ")
        if record_pattern == "dispersed bubble":
            record_pattern = "dispersed bubbles"
        assert pattern == record_pattern, run
        assert im.shape[1] == records[run]["comprimento"], run
        assert im.shape[0] == 80 and im.shape[1] > 1
        width = im.shape[1]
        mask, accepted, _ = model.silhouette(im, reference)
        assert mask.dtype == bool and mask.shape == im.shape
        assert not mask[:18].any() and not mask[72:].any()
        alpha = model.geometry(mask)
        independent = independent_geometry(mask)
        assert np.array_equal(alpha, independent), run
        assert np.all((0 <= alpha) & (alpha <= 1))
        mask_path = output("intermediate/phase_masks/" + run + "_mask.png")
        Image.fromarray(mask.astype(np.uint8)*255).save(mask_path)
        with Image.open(mask_path) as image:
            saved_mask = np.asarray(image) == 255
        assert np.array_equal(saved_mask, mask)
        core = saved_mask[18:72]
        count = int(np.count_nonzero(core[:, 1:] != core[:, :-1]))
        independent_count = sum(int(np.count_nonzero(np.diff(row.astype(np.int8)))) for row in core)
        assert count == independent_count
        assert count == np.count_nonzero(core[:, ::-1][:, 1:] != core[:, ::-1][:, :-1])
        comparisons_count = 54*(width-1)
        rate = count / comparisons_count
        saturation = float((im[20:60] >= 255).mean())
        eligible = saturation <= 0.20
        center = mask[40]
        center_count = int(np.count_nonzero(center[1:] != center[:-1]))
        entries = int(np.count_nonzero(~center[:-1] & center[1:]))
        exits = int(np.count_nonzero(center[:-1] & ~center[1:]))
        assert entries + exits == center_count
        touch = float((mask[18] | mask[71]).mean())
        alpha_mean = float(alpha.mean())
        void_rows.append(dict(
            run=run, label=label, panorama_width_px=width, lumen_start_row=18,
            lumen_stop_row_exclusive=72, lumen_diameter_px=54.0, boundary_source=BOUNDARY_SOURCE,
            gas_silhouette_projected_fraction=float(mask.sum()/(54*width)),
            conditional_circle_alpha_spatial_mean=alpha_mean,
            conditional_circle_alpha_spatial_median=float(np.median(alpha)),
            detected_center_gas_column_fraction=float(accepted.mean()),
            center_gas_silhouette_duty=float(center.mean()),
            radial_crop_touch_column_fraction=touch, core_saturated_fraction=saturation,
            validity="conditional_optical_geometry" if eligible else "unsupported_saturation"))
        phase_rows.append(dict(
            run=run, pattern=pattern, panorama_width=width, candidate_mask_rows="18:72",
            candidate_mask_n_rows=54, candidate_pixel_flip_count=count,
            candidate_pixel_flip_comparisons=comparisons_count, candidate_pixel_flip_rate=rate,
            center_candidate_flip_count=center_count, center_candidate_flip_rate=center_count/(width-1),
            center_candidate_flips_per1000px=1000*center_count/(width-1),
            center_candidate_gas_entries=entries, center_candidate_gas_exits=exits,
            center_candidate_gas_duty=float(center.mean()), paired_shadow_detection_fraction=float(accepted.mean()),
            eligible=eligible, phase_status="inferred gas silhouette; optical calibration conditional"))
        plot_rows.append(dict(
            run=run, supplied_pattern=pattern, source_void_label=label,
            eligible_pure_unsaturated=eligible, pure_pattern=True, plotted=eligible,
            plot_exclusion_reason="" if eligible else "unsupported_saturation",
            saturated=not eligible, saturation_threshold=0.20, core_saturated_fraction=saturation,
            void_validity=void_rows[-1]["validity"], radial_crop_touch_column_fraction=touch,
            radial_crop_touched=touch > 0, all_columns_touch_radial_crop=touch == 1.0,
            full_lumen_maximum=1.0, lumen_start_row=18, lumen_stop_row_exclusive=72,
            lumen_diameter_px=54.0, boundary_source=BOUNDARY_SOURCE,
            conditional_circle_alpha_spatial_mean=alpha_mean, panorama_width_px=width,
            candidate_mask_rows="18:72", candidate_mask_n_rows=54,
            candidate_pixel_flip_count=count, candidate_pixel_flip_comparisons=comparisons_count,
            candidate_pixel_flip_rate=rate, candidate_pixel_flips_per1000px=1000*rate,
            center_candidate_flip_count=center_count, center_candidate_flip_comparisons=width-1,
            center_candidate_flip_rate=center_count/(width-1),
            center_candidate_flips_per1000px=1000*center_count/(width-1),
            center_candidate_gas_entries=entries, center_candidate_gas_exits=exits,
            phase_status=phase_rows[-1]["phase_status"]))
        comparisons.append(dict(
            run=run, pattern=pattern, old_eligible=old_phase[run]["eligible"] == "True", new_eligible=eligible,
            old_width=int(old_void[run]["panorama_width_px"]), new_width=width,
            old_lumen_height_px=54, new_lumen_height_px=54,
            old_saturated_fraction=float(old_void[run]["core_saturated_fraction"]), new_saturated_fraction=saturation,
            old_alpha=float(old_void[run]["conditional_circle_alpha_spatial_mean"]), new_alpha=alpha_mean,
            old_F1000=1000*float(old_phase[run]["candidate_pixel_flip_rate"]), new_F1000=1000*rate,
            old_flip_count=int(old_phase[run]["candidate_pixel_flip_count"]), new_flip_count=count))
        checks.append(dict(run=run, width=width, independent_geometry_max_error=float(np.max(abs(alpha-independent))),
                           independent_flip_count=count, denominator=comparisons_count,
                           saved_mask_identical=True, no_foreground_outside_roi=True,
                           saturation_recomputed_from_new_panorama=True, eligible=eligible))
        print(f"{run}: alpha={alpha_mean:.8f}, 1000F={1000*rate:.6f}, eligible={eligible}", flush=True)

    void_groups, phase_groups, group_summary = [], [], []
    for pattern in PATTERNS:
        subset = [row for row in plot_rows if row["plotted"] and row["supplied_pattern"] == pattern]
        alphas = [row["conditional_circle_alpha_spatial_mean"] for row in subset]
        rates = [row["candidate_pixel_flip_rate"] for row in subset]
        void_groups.append(dict(pattern=pattern, n=len(subset), conditional_alpha_median=float(np.median(alphas)),
                                conditional_alpha_min=min(alphas), conditional_alpha_max=max(alphas),
                                conditional_alpha_mean=float(np.mean(alphas))))
        phase_groups.append(dict(pattern=pattern, n=len(subset),
                                 candidate_pixel_flip_rate_median=float(np.median(rates)),
                                 candidate_pixel_flip_rate_mean=float(np.mean(rates))))
        group_summary.append(dict(pattern=pattern, n=len(subset), alpha_median=float(np.median(alphas)),
                                  F1000_median=float(np.median(rates))*1000,
                                  alpha_arithmetic_mean=float(np.mean(alphas)),
                                  F1000_arithmetic_mean=float(np.mean(rates))*1000))
    eligible_runs = [row["run"] for row in plot_rows if row["plotted"]]
    excluded_runs = [row["run"] for row in plot_rows if not row["plotted"]]
    assert len(eligible_runs) + len(excluded_runs) == 26
    write_csv("results/void_fraction/panorama_void_per_run.csv", void_rows)
    write_csv("results/void_fraction/panorama_void_by_pattern.csv", void_groups)
    write_csv("results/flip_rate/phase_mask_flip_per_run.csv", phase_rows)
    write_csv("results/flip_rate/phase_mask_flip_by_pattern.csv", phase_groups)
    write_csv("results/void_fraction/void_fraction_flip_plot_data.csv", plot_rows)
    write_csv("results/summary/optical_group_summary.csv", group_summary)
    write_csv("results/validation/change_from_v14_2.csv", comparisons)
    hulls = hull_checks(plot_rows)
    write_json("results/validation/hull_intersections.json", hulls)
    old_groups = read_csv(BASELINE / "results/summary/optical_group_summary.csv")
    old_hulls = json.loads((BASELINE / "results/validation/hull_intersections.json").read_text())
    interpretation_changes = dict(
        eligibility_changed=[row["run"] for row in comparisons if row["old_eligible"] != row["new_eligible"]],
        old_hull_intersections=old_hulls["overlap_pairs"], new_hull_intersections=hulls["overlap_pairs"],
        highest_alpha_median=max(group_summary, key=lambda row: row["alpha_median"])["pattern"],
        lowest_F1000_median=min(group_summary, key=lambda row: row["F1000_median"])["pattern"])
    write_json("resumo_tcc.json", dict(
        status="optical_complete_pending_spectral_examples", processed_acquisitions=26,
        eligible_acquisitions=len(eligible_runs), eligible_runs=eligible_runs,
        excluded=[dict(run=row["run"], reason="core saturation above 20 percent in rows 20:60",
                       saturated_fraction=row["core_saturated_fraction"]) for row in plot_rows if not row["plotted"]],
        geometry=dict(rows="18:72", included_rows="18 through 71", height=54, diameter_px=54),
        saturation_rows="20:60", groups=group_summary, previous_groups=old_groups,
        reference=reference_comparison, hulls=hulls, changes=interpretation_changes,
        discussion_examples=[],
        figure_paths=["TCC/figuras/v18_optica/panorama_phase_masks.png",
                      "TCC/figuras/v18_optica/void_fraction_vs_flip_rate.png"]))
    write_json("results/validation/independent_arithmetic_checks.json", checks)
    after_hashes = {relative(path): sha256(path) for path in inputs}
    assert before_hashes == after_hashes
    write_json("results/validation/input_hashes.json", before_hashes)
    write_json("results/validation/metadata.json", dict(
        version="v18_montagem_18_71", purpose="Optical reanalysis after panorama reconstruction using rows 18:72",
        command="python3 -B src2/v18_montagem_18_71/scripts/process_optical.py",
        source_version=relative(ROOT / "01_panoramas"), previous_optical_version=relative(BASELINE),
        source_algorithm=relative(SOURCE), source_algorithm_sha256=SOURCE_SHA256,
        isolated_unchanged_functions=FUNCTION_NAMES, historical_module_imported=False,
        geometry=dict(start=18, stop_exclusive=72, included_rows="18 through 71", height=54, diameter_px=54),
        geometry_default_arguments_verified=True, boundary_source=BOUNDARY_SOURCE,
        segmentation=dict(shadow_threshold=0.70, radial_gaussian_sigma=1, morphological_closing_shape=[3, 3],
                          minimum_open_interval=20, minimum_isolated_component_pixels=4,
                          upper_search_slice=[19, 38], lower_search_slice=[42, 71],
                          search_slices_unchanged=True, radial_extension_can_include_row_71=True),
        reference=dict(runs=REFERENCE_RUNS, column_normalization_slice=[30, 50], percentile=75,
                       selection_std_slice=[20, 68], selection_quantile=0.05, comparison=reference_comparison),
        saturation=dict(slice=[20, 60], saturated_value=255, maximum_fraction=0.20,
                        rule_unchanged_from_v14_2=True, values_recomputed=True, excluded_runs=excluded_runs),
        selected_runs=runs, processed_acquisitions=26, eligible_acquisitions=len(eligible_runs),
        group_counts={row["pattern"]: row["n"] for row in group_summary},
        alpha_formula="mean_over_columns(sum_over_contiguous_radial_intervals((width/54)**2))",
        flip_formula="longitudinal_neighbor_changes / (54*(panorama_width-1))",
        mask_encoding="255=inferred gas silhouette, 0=unmarked; 80-row image, foreground confined to rows 18:72",
        table_summary="Median across eligible acquisitions of individual spatial means and flip densities",
        chart_stars="Arithmetic coordinate means with equal weight for each eligible acquisition",
        no_parameter_sweep=True, no_new_statistical_inference=True, no_wavelet_processing=True,
        independent_geometry_and_flip_arithmetic_verified=True, source_inputs_hashes_unchanged=True,
        runtime=dict(python=platform.python_version(), executable=sys.executable, numpy=np.__version__,
                     scipy=scipy.__version__, pillow=PIL.__version__)))
    manifest_path = OUTROOT / "results/validation/artifact_manifest.json"
    outputs = {path.relative_to(OUTROOT).as_posix(): sha256(path) for path in sorted(OUTROOT.rglob("*"))
               if path.is_file() and path != manifest_path and "__pycache__" not in path.parts}
    script_hashes = {relative(path): sha256(path) for path in [Path(__file__).resolve(), ROOT / "scripts/optical_model.py"]}
    write_json("results/validation/artifact_manifest.json", dict(hash_algorithm="sha256", files=outputs,
                                                                 scripts=script_hashes))
    print(json.dumps(group_summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
