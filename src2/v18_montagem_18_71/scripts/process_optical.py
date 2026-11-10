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
                           REFERENCE_RUNS, THRESHOLD, load_model)

from runtime import DATA_ROOT, OUTPUT_ROOT, SCRIPTS_ROOT, artifact_path, reference

OUTROOT = OUTPUT_ROOT / "03_optica"
PANORAMA_MANIFEST = OUTPUT_ROOT / "01_panoramas/manifesto.json"
DATA = DATA_ROOT
PATTERNS = ("dispersed bubbles", "plug", "slug", "churn", "wavy annular")
BOUNDARY_SOURCE = "User correction2026-10-30: all current optical operations use only rows18 through71 inclusive, Python slice18:72; H=D=54px."


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path):
    return reference(path)


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
    reconstruction = json.loads(PANORAMA_MANIFEST.read_text())
    assert reconstruction["complete"] is True
    records = {record["ensaio"]: record for record in reconstruction["records"]}
    runs = sorted(records)
    panorama_paths = {run: artifact_path(records[run]["optica_path"]) for run in runs}
    assert len(reconstruction["records"]) == 26
    for run, path in panorama_paths.items():
        assert sha256(path) == records[run]["hashes"]["optica"], run
    assert len(runs) == len(set(runs)) == 26
    inputs = [PANORAMA_MANIFEST, SCRIPTS_ROOT / "optical_model.py"]
    markers = {}
    for run in runs:
        candidates = list((DATA / run).glob("__*.txt"))
        assert len(candidates) == 1, run
        markers[run] = candidates[0]
        inputs.extend([candidates[0], panorama_paths[run]])
    before_hashes = {relative(path): sha256(path) for path in inputs}

    model = load_model(panorama_paths)
    full = np.zeros((80, 25), dtype=bool)
    full[18:72] = True
    assert np.array_equal(model.geometry(full), np.ones(25))
    probe = np.full((80, 25), 100.0)
    probe[71] = 255
    assert model.roi_saturation(probe) == 1/54
    assert np.all(model.normalize(probe)[-1] == 2.55)
    probe_before = model.normalize(probe)
    probe[71] = 100
    assert not np.array_equal(probe_before, model.normalize(probe))
    reference, selections = model.make_reference()
    reference_metadata = dict(reference_runs=REFERENCE_RUNS,
                              selected_columns={run: selections[run]["n_selected"] for run in REFERENCE_RUNS},
                              height=len(reference), roi=list(LUMEN_ROWS))
    write_json("intermediate/features/panorama_liquid_reference.json", dict(
        reference_runs=REFERENCE_RUNS, liquid_candidates=selections, radial_reference=reference.tolist(),
        lumen_rows=LUMEN_ROWS, lumen_diameter_px=DIAMETER_PX,
        boundary_source=BOUNDARY_SOURCE))

    void_rows, phase_rows, plot_rows, checks = [], [], [], []
    for run in runs:
        im = model.read(run)
        label = markers[run].stem[2:]
        pattern = label.replace("_", " ")
        assert pattern in PATTERNS
        record_pattern = records[run]["padrao"].replace("_", " ")
        if record_pattern == "dispersed bubble":
            record_pattern = "dispersed bubbles"
        assert pattern == record_pattern, run
        assert im.shape[1] == records[run]["comprimento"], run
        assert im.shape[0] == 80 and im.shape[1] > 1
        width = im.shape[1]
        mask, accepted, ratio = model.silhouette(im, reference)
        perturbed = im.copy()
        perturbed[:18] = 255 - perturbed[:18]
        perturbed[72:] = 255 - perturbed[72:]
        alt_mask, alt_accepted, alt_ratio = model.silhouette(perturbed, reference)
        assert np.array_equal(model.normalize(im), model.normalize(perturbed)), run
        assert model.roi_saturation(im) == model.roi_saturation(perturbed), run
        assert np.array_equal(mask, alt_mask) and np.array_equal(accepted, alt_accepted), run
        assert np.array_equal(ratio, alt_ratio), run
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
        assert np.array_equal(im[18:72] >= 255, im[18:72] == 255), run
        saturation = model.roi_saturation(im)
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
        checks.append(dict(run=run, width=width, independent_geometry_max_error=float(np.max(abs(alpha-independent))),
                           independent_flip_count=count, denominator=comparisons_count,
                           saved_mask_identical=True, no_foreground_outside_roi=True,
                           saturation_recomputed_from_roi=True, outside_roi_intensity_invariance=True,
                           reference_candidates_outside_roi_invariant=True, eligible=eligible))
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
    hulls = hull_checks(plot_rows)
    write_json("results/validation/hull_intersections.json", hulls)
    extrema = dict(
        highest_alpha_median=max(group_summary, key=lambda row: row["alpha_median"])["pattern"],
        lowest_F1000_median=min(group_summary, key=lambda row: row["F1000_median"])["pattern"])
    write_json("resumo_tcc.json", dict(
        status="optical_complete_pending_spectral_examples", processed_acquisitions=26,
        eligible_acquisitions=len(eligible_runs), eligible_runs=eligible_runs,
        excluded=[dict(run=row["run"], reason="saturation above20percent in the sole ROI18:72",
                       saturated_fraction=row["core_saturated_fraction"]) for row in plot_rows if not row["plotted"]],
        geometry=dict(rows="18:72", included_rows="18 through 71", height=54, diameter_px=54),
        saturation_rows="18:72", groups=group_summary,
        reference=reference_metadata, hulls=hulls, descriptive_extrema=extrema,
        discussion_examples=[],
        figure_paths=["03_optica/figuras/panorama_phase_masks.png",
                      "03_optica/figuras/void_fraction_vs_flip_rate.png"]))
    write_json("results/validation/independent_arithmetic_checks.json", checks)
    write_json("results/validation/single_roi_checks.json", dict(
        status="PASS", roi=[18, 72], full_lumen_geometry_equals_one=True,
        row71_included_in_saturation_and_normalization=True,
        row71_only_saturation_fraction=1/54, checked_acquisitions=len(checks),
        saturation_equality_255_matches_greater_equal_for_uint8_inputs=True,
        saturation_equality_checked_runs=runs,
        outside_roi_invariance_of_normalization_reference_candidates_segmentation_saturation=True,
        independent_interval_geometry_and_neighbor_counts=True))
    after_hashes = {relative(path): sha256(path) for path in inputs}
    assert before_hashes == after_hashes
    write_json("results/validation/input_hashes.json", before_hashes)
    write_json("results/validation/metadata.json", dict(
        version="v18_montagem_18_71", purpose="All optical intensity operations and measures confined to the single ROI18:72",
        command="python3 -B scripts/process_optical.py",
        source_panoramas=relative(OUTPUT_ROOT / "01_panoramas"),
        source_algorithm=relative(SCRIPTS_ROOT / "optical_model.py"),
        source_algorithm_sha256=sha256(SCRIPTS_ROOT / "optical_model.py"),
        optical_functions=FUNCTION_NAMES,
        geometry=dict(start=18, stop_exclusive=72, included_rows="18 through 71", height=54, diameter_px=54),
        geometry_default_arguments_verified=True, boundary_source=BOUNDARY_SOURCE,
        segmentation=dict(shadow_threshold=0.70, radial_gaussian_sigma=1, morphological_closing_shape=[3, 3],
                          minimum_open_interval=20, minimum_isolated_component_pixels=4,
                          upper_search_slice=[18, 45], lower_search_slice=[45, 72],
                          search_partition="equal halves derived from ROI height",
                          smoothing_input_slice=[18, 72], gaussian_mode="reflect",
                          radial_extension_can_include_row_71=True),
        reference=dict(runs=REFERENCE_RUNS, column_normalization_slice=[18, 72], percentile=75,
                       selection_std_slice=[18, 72], selection_quantile=0.05, reference_height=54),
        saturation=dict(slice=[18, 72], saturated_value=255, maximum_fraction=0.20,
                        comparison="equality to255 in8bit inputs",
                        excluded_runs=excluded_runs),
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
    script_hashes = {relative(path): sha256(path) for path in [Path(__file__).resolve(), SCRIPTS_ROOT / "optical_model.py"]}
    write_json("results/validation/artifact_manifest.json", dict(hash_algorithm="sha256", files=outputs,
                                                                 scripts=script_hashes))
    print(json.dumps(group_summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
