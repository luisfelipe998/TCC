"""Audit image-only void alternatives; no regime labels enter any estimator.

Outputs are optical descriptors/conditional geometric reconstructions, not
calibrated volumetric void measurements. Uses only src2/data_processed images.
"""
from v14_paths import panorama_path, original_panorama_path, BASELINE
from v14_paths import ROOT as HERE, DATA, CACHE, artifact_path, write_markdown
from v14_selection import iter_runs, is_eligible
from pathlib import Path
import csv
import json
import os
import sys
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter1d, binary_closing, binary_fill_holes

sys.path.insert(0, str(HERE.parent / "v11.3"))
from void_common_v11_3 import manual_otsu
os.environ.setdefault("MPLCONFIGDIR", str(CACHE / "matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

N_FRAMES = 120
AXIAL_STEP = 8
# These are apparent pixel diameters, not an optical calibration.
D_CANDIDATES = (56.0, 60.0, 64.0)


def files_for(run):
    files = sorted(run.glob("ImgA*.jpg"))
    idx = np.linspace(0, len(files) - 1, min(N_FRAMES, len(files))).round().astype(int)
    return [files[i] for i in idx]


def center_width(mask, center):
    # Shape (radial rows, axial columns), includes the centre exactly once.
    upper = np.cumprod(mask[:center + 1][::-1], axis=0).sum(axis=0)
    lower = np.cumprod(mask[center:], axis=0).sum(axis=0)
    return upper + lower - mask[center]


def apparent_interfaces(profile):
    top = 19 + np.argmin(profile[19:38], axis=0)
    bottom = 42 + np.argmin(profile[42:70], axis=0)
    cols = np.arange(profile.shape[1])
    center = np.median(profile[35:46], axis=0)
    contrast = np.minimum(center - profile[top, cols], center - profile[bottom, cols])
    interior = (top > 19) & (top < 37) & (bottom > 42) & (bottom < 69)
    return top, bottom, contrast, interior


def expanded_interfaces(profile):
    # Explicit window-sensitivity check: allow minima near the lower reflection
    # band, still reject the boundary pixel at row 72. This is not calibration.
    top = 19 + np.argmin(profile[19:38], axis=0)
    bottom = 42 + np.argmin(profile[42:73], axis=0)
    cols = np.arange(profile.shape[1])
    center = np.median(profile[35:46], axis=0)
    contrast = np.minimum(center - profile[top, cols], center - profile[bottom, cols])
    gate = (top > 19) & (top < 37) & (bottom > 42) & (bottom < 72) & (contrast >= 15)
    return top, bottom, gate


def baseline_profile():
    # Liquid-like negative-control run has some real dispersed bubbles. Robust
    # profile is a pseudo-reference, never asserted to be an all-liquid photo.
    ref = DATA / "m07q01_001"
    prof = []
    for f in files_for(ref)[::4]:
        im = np.asarray(Image.open(f).convert("L"), dtype=float)
        scaled = im / np.maximum(np.median(im[:8], axis=0), 1)[None, :]
        prof.append(np.median(scaled, axis=1))
    return gaussian_filter1d(np.median(prof, axis=0), 1)


def summarize(run, reference):
    bc_otsu, bc_mid, contour, dark, geo, tmpl = [], [], [], [], [], []
    contrasts, gates, template_gates = [], [], []
    expanded_widths, expanded_gates = [], []
    frame_means, saturation, profiles = [], [], []
    all_files = sorted(run.glob("ImgA*.jpg"))
    for f in files_for(run):
        im = np.asarray(Image.open(f).convert("L"), dtype=float)
        win = im[10:70, ::AXIAL_STEP]
        # Match v11.3 Otsu histogram over all axial pixels.
        otsu = manual_otsu(im[10:70].ravel())
        w_otsu = center_width(win > otsu, 30)
        midpoint = (win.min(axis=0) + win.max(axis=0)) / 2
        w_mid = center_width(win > midpoint[None, :], 30)
        bc_otsu.append(w_otsu)
        bc_mid.append(w_mid)
        frame_means.append(w_mid.mean() / 60)
        saturation.append(float((im[20:60] >= 255).mean()))
        profile = gaussian_filter1d(im, 1, axis=0)[:, ::AXIAL_STEP]
        profiles.append(np.median(profile, axis=1))
        top, bottom, contrast, interior = apparent_interfaces(profile)
        geo.append(bottom - top)
        contrasts.append(contrast)
        gates.append(interior)
        etop, ebottom, egate = expanded_interfaces(profile)
        expanded_widths.append(ebottom - etop)
        expanded_gates.append(egate)
        scale = np.maximum(np.median(im[:8, ::AXIAL_STEP], axis=0), 1)
        ratio = profile / scale[None, :] / np.maximum(reference[:, None], .01)
        u = 19 + np.argmin(ratio[19:38], axis=0)
        b = 42 + np.argmin(ratio[42:70], axis=0)
        cols = np.arange(ratio.shape[1])
        gate = (ratio[u, cols] < .65) & (ratio[b, cols] < .65)
        gate &= (u > 19) & (u < 37) & (b > 42) & (b < 69)
        tmpl.append(b - u)
        template_gates.append(gate)
        # Filling dark closed contours can include gas interiors, but open
        # annular interfaces cannot be filled by this topology operation.
        crop = im[18:71]
        mask = crop <= manual_otsu(crop.ravel())
        closed = binary_closing(mask, structure=np.ones((3, 3)))
        dark.append(mask.mean())
        contour.append(binary_fill_holes(closed).mean())

    widths = np.concatenate(bc_mid)
    otsu_widths = np.concatenate(bc_otsu)
    gwidth = np.concatenate(geo).astype(float)
    contrast = np.concatenate(contrasts)
    gate_interior = np.concatenate(gates)
    twidth = np.concatenate(tmpl).astype(float)
    tgate = np.concatenate(template_gates)
    ewidth = np.concatenate(expanded_widths).astype(float)
    egate = np.concatenate(expanded_gates)
    marker = [x.stem.lstrip("_") for x in run.glob("__*.txt")]
    label = marker[0] if marker else "unknown"
    row = dict(run=run.name, label=label, n_frames_total=len(all_files),
               n_frames_sampled=len(frame_means), n_axial_profiles=len(widths),
               saturated_core_fraction_mean=np.mean(saturation),
               bright_mid_chord_mean=np.mean(widths / 60),
               bright_mid_chord_median=np.median(widths / 60),
               bright_otsu_chord_mean=np.mean(otsu_widths / 60),
               bright_otsu_chord_median=np.median(otsu_widths / 60),
               dark_projected_fraction_mean=np.mean(dark),
               closed_contour_projected_fraction_mean=np.mean(contour),
               template_interface_detection_fraction=np.mean(tgate),
               expanded_interface_detection_fraction=np.mean(egate),
               expanded_interface_zero_fill_circle_D60=np.mean(egate * (ewidth / 60) ** 2),
               expanded_interface_detected_circle_D60=np.mean((ewidth[egate] / 60) ** 2) if egate.any() else float("nan"))
    for diameter in D_CANDIDATES:
        key = int(diameter)
        row[f"bright_width_conditional_circle_D{key}_mean"] = np.mean((widths / diameter) ** 2)
        row[f"template_envelope_zero_fill_circle_D{key}"] = np.mean(tgate * (twidth / diameter) ** 2)
    for threshold in (10, 15, 20):
        gate = gate_interior & (contrast >= threshold)
        row[f"interface_detection_contrast{threshold}"] = gate.mean()
        row[f"interface_envelope_zero_fill_circle_D60_contrast{threshold}"] = np.mean(gate * (gwidth / 60) ** 2)
        row[f"interface_envelope_detected_circle_D60_contrast{threshold}"] = np.mean((gwidth[gate] / 60) ** 2) if gate.any() else float("nan")
    blocks = np.array_split(np.asarray(frame_means), 5)
    row["bright_mid_five_block_mean_min"] = min(x.mean() for x in blocks)
    row["bright_mid_five_block_mean_max"] = max(x.mean() for x in blocks)
    print(run.name, label, "bright mean/median", round(row["bright_mid_chord_mean"], 3),
          round(row["bright_mid_chord_median"], 3), "interface detection", round(row["interface_detection_contrast15"], 3), flush=True)
    return row, np.median(profiles, axis=0)


def write_report(rows):
    lines = ["# Image-only void alternatives: numerical audit", "",
             "**Historical descriptor comparison:** D60 and D56/D64 fields below reproduce the earlier raw-frame alternatives. They are not current void estimates. The primary panorama analysis now uses the user-confirmed complete interior 18:71 and D=53 pixels; see the [normalization correction](lumen_normalization_review.md). Raw-frame mappings are not silently recalibrated by this panorama boundary identification.", "",
             f"The current scope contains {len(rows)} retained acquisitions. Mixed labels, smooth annular and the user-excluded replicate are omitted before every raw-frame estimator and graphic.", "",
             "All numerical estimators use raw frames from `src2/data_processed`. Regime labels are attached after computation. No estimator is tuned to a desired regime ranking. 120 uniformly spaced frames and every eighth axial column are sampled per sequence (501 original frames). These correlated observations do not constitute independent experimental replicates.", "",
             "**Result:** replacing the median with a mean is necessary for a time/space average, but neither statistic calibrates brightness into gas volume. Dark-interface envelopes and closed-contour filling recover additional gas-outline information. They remain optical proxies because detected shadow positions are refracted interface images; a geometric diameter cannot yet be assigned to those shadows.", "",
             "## Per-run results", "",
             "BC values are center bright-width / 60 pixels. Envelope values below are zero-filled *conditional circular reconstructions*, E[detected × (apparent interface separation/60)²]. Zero means failed detection or no detected center gas, not a proven all-liquid section. They are not calibrated void fractions.", "",
             "| Run | Label | BC mean | BC median | Interface detected | Envelope surrogate | Filled dark projected area | Template detected |",
             "|---|---|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        lines.append(f"| {r['run']} | {r['label']} | {r['bright_mid_chord_mean']:.3f} | {r['bright_mid_chord_median']:.3f} | {r['interface_detection_contrast15']:.3f} | {r['interface_envelope_zero_fill_circle_D60_contrast15']:.3f} | {r['closed_contour_projected_fraction_mean']:.3f} | {r['template_interface_detection_fraction']:.3f} |")
    lines += ["", "## What each alternative resolves", "",
              "1. **Mean instead of median.** Physical averaged void is an arithmetic expectation of instantaneous gas area. v11.3 returns the median bright width across frames and columns, which can ignore the less common phase. This audit records both statistics and five consecutive frame-block means. The bright-width numerator still does not equal gas diameter.",
              "2. **Dark-interface envelope.** Smooth radial intensity with σ=1 pixel; locate upper/lower dark extrema within rows 19:38 and 42:70. Accept interior extrema when both are ≥15 gray levels below the center. Repeat at 10 and 20 levels. All gates are identical across runs. This uses the boundary, including a bright interior, rather than declaring every shadow pixel to be gas. An absence of accepted interfaces is an unresolved observation.",
              "The CSV also records an expanded lower search, rows 42:73, rejecting row 72. It exposes substantial search-window dependence when gas shadows overlap the lower wall/reflection band. Better detection coverage does not establish that the optical shadow is the physical gas boundary.",
              "3. **Circular reconstruction.** If an undistorted side view identifies the physical connected gas-core diameter d_g, and the core is concentric and axisymmetric, α(t,z)=(d_g/D)². Averaging the square matters: E[(d_g/D)²] is generally different from E[d_g/D]². Absolute px→mm scale cancels from this ratio. A second camera is not mathematically mandatory under that assumption. The present bright band and dark shadows have no validated mapping to d_g, however; in churn, off-axis fragments further violate this reconstruction. The former diameter candidates 56, 60 and 64 pixels remain historical descriptor sensitivity tests only, superseded for the primary panorama estimate by the confirmed 53-pixel interior.",
              "4. **Filled dark contours.** Otsu-segment rows 18:71, close 3×3 gaps, then fill enclosed holes. This captures interiors of some dark-outlined bubbles. Annular interfaces may connect to image edges and fail to enclose any region, while merging shadows can falsely enclose liquid. The projected image area is not the circular cross-sectional area.",
              "5. **Liquid-like template.** Build a robust normalized radial profile from `m07q01_001`, normalize each column using its upper eight background rows, and detect paired deficits below 65% of the template. This explicitly exposes the reference dependence. `m07q01_001` contains bubbles and is not a validated liquid-only calibration photo; cross-run exposure, gamma and reflections can differ. The resulting template detector is a feasibility probe.", "",
              "## Negative control and identifiability", "",
              "The liquid-like `m07q01_001` sequence has visible dispersed bubbles but a nearly featureless bright center. Its bright-core midpoint median can be high because that operation measures brightness shape, not phase occupancy. A per-run adaptive threshold necessarily finds image classes even without a separable gas/liquid signal. Likewise an own-run temporal percentile is not automatically a liquid baseline: a persistent annular core can occupy the same pixel throughout the record.", "",
              "The 80-row acquisition ROI includes the wall/reflection bands. The windows 20:60 or 10:70 are analysis crops, not demonstrated physical tube diameters. Even an apparent wall spacing is affected by cylindrical refraction. A proper inverse model requires the tube's optical geometry (including glass outer diameter and refractive index), camera/illumination geometry, fluid refractive index, and validation against liquid-only and known-volume states or another independent void measurement.", "",
              "## Source check", "",
              "Schepperle et al. (2024), Sensors 24, 3363, DOI 10.3390/s24113363, uses a **rectangular** channel (depth 0.5 mm, width 1.5 mm; §2.1), contour extraction, and histogram comparison against a **separate liquid-only image of the same subsection acquired before boiling** (§2.3, pp. 5–7). Its Eq. (3) counts pixels *after phase labeling*. It does not validate replacing contour/phase classification with uncalibrated bright-pixel counting or an own-run 20th percentile in a circular tube.", "",
              "**Current next step:** use the retained panorama liquid intervals and gas outlines described in the main reconstruction report. These earlier descriptors remain optical diagnostics. Validate masks and circular/refraction assumptions; correcting a rank order does not establish physical calibration."]
    write_markdown(artifact_path("image_void_alternatives.md"), "\n".join(lines) + "\n")


def main():
    reference = baseline_profile()
    rows, profiles = [], {}
    for run in iter_runs():
        row, profile = summarize(run, reference)
        rows.append(row)
        profiles[run.name] = profile
    with (artifact_path("image_void_alternatives.csv")).open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    # Match the primary analysis selection using panorama-core saturation.
    # The raw-frame per-run values above are retained as diagnostic records.
    plot_rows = []
    for row in rows:
        panorama = np.asarray(Image.open(panorama_path(row["run"])).convert("L"))
        saturation = float((panorama[20:60] == 255).mean())
        if is_eligible(row["run"], row["label"].replace("_", " "), saturation):
            plot_rows.append(row)
    (artifact_path("image_void_alternatives_metadata.json")).write_text(json.dumps(dict(
        n_acquisition_runs=len(rows), included_acquisitions=[r["run"] for r in rows],
        n_frames_sampled=N_FRAMES, axial_step=AXIAL_STEP,
        diameter_candidates_px=D_CANDIDATES,
        bright_width_plot_runs=[r["run"] for r in plot_rows],
        bright_width_plot_excluded_runs=[r["run"] for r in rows if r not in plot_rows],
        bright_width_plot_selection="Single-pattern acquisitions with at most 20% saturated pixels in panorama rows 20:60; individual categorical points, no connecting trajectory",
        interpretation="Optical proxy/conditional geometric reconstruction; not calibrated physical void"), indent=2) + "\n")
    write_report(rows)
    chosen = ["m07q01_001", "m06q01_001", "m06q05_001", "m06q20_001", "m06q80_001"]
    fig, axes = plt.subplots(3, 2, figsize=(12, 9), sharex=True, sharey=True)
    for ax, run in zip(axes.flat, chosen):
        ax.plot(np.arange(80), profiles[run], color="black")
        ax.axvspan(0, 18, alpha=.12, color="red")
        ax.axvspan(71, 80, alpha=.12, color="red")
        ax.set_title(run); ax.set_xlabel("Radial image row (px)"); ax.set_ylabel("Median gray level")
        ax.grid(alpha=.2)
    for ax in axes.flat[len(chosen):]:
        ax.set_visible(False)
    fig.suptitle("Radial optical profiles: bright interior and dark boundaries are different observables")
    fig.legend(handles=[Patch(facecolor="red", alpha=.12,
                              label="Outside the nominal 18:71 row window; not calibrated wall locations")],
               loc="lower center", frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, .045, 1, .965)); fig.savefig(artifact_path("radial_interface_profiles.png"), dpi=180); plt.close(fig)
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(plot_rows))
    ax.scatter(x - .12, [r["bright_mid_chord_mean"] for r in plot_rows], marker="o", label="Bright-width arithmetic mean")
    ax.scatter(x + .12, [r["bright_mid_chord_median"] for r in plot_rows], marker="s", label="Bright-width median (v11.3 statistic)")
    ax.set_xticks(x); ax.set_xticklabels([r["run"].replace("_001", "") for r in plot_rows], rotation=90, fontsize=8)
    ax.set_title(f"Bright-width comparison: {len(plot_rows)} eligible acquisitions")
    ax.set_ylabel("Bright width / 60 px: optical descriptor"); ax.grid(alpha=.2); ax.legend()
    fig.tight_layout(); fig.savefig(artifact_path("bright_width_mean_vs_median.png"), dpi=180); plt.close(fig)
    fig, axes = plt.subplots(len(chosen), 1, figsize=(13, 2 * len(chosen)), sharex=True)
    for ax, run in zip(axes, chosen):
        im = np.asarray(Image.open(DATA / run / "ImgA000250.jpg").convert("L"), dtype=float)
        p = gaussian_filter1d(im, 1, axis=0)
        t, b, gate = expanded_interfaces(p)
        x = np.arange(im.shape[1])
        ax.imshow(im, cmap="gray", vmin=0, vmax=255, aspect="auto")
        ax.plot(x, np.where(gate, t, np.nan), color="cyan", lw=.6)
        ax.plot(x, np.where(gate, b, np.nan), color="cyan", lw=.6)
        ax.set_ylabel(run, fontsize=8); ax.set_yticks([0, 40, 79])
    axes[-1].set_xlabel("Axial image column (px)")
    fig.suptitle("Dark-interface detector: accepted apparent boundaries in cyan (uncalibrated)")
    fig.tight_layout(); fig.savefig(artifact_path("interface_detector_overlay.png"), dpi=180); plt.close(fig)

    fig, axes = plt.subplots(len(chosen), 1, figsize=(13, 1.7 * len(chosen)))
    for ax, run in zip(axes, chosen):
        im = np.asarray(Image.open(DATA / run / 'ImgA000250.jpg').convert('L'))
        ax.imshow(im, cmap='gray', vmin=0, vmax=255, aspect='auto')
        ax.set_ylabel(run, fontsize=9); ax.set_xticks([]); ax.set_yticks([0, 40, 79])
    fig.suptitle('Retained representative raw frames: ImgA000250.jpg')
    fig.tight_layout(); fig.savefig(artifact_path('representative_raw_images.png'), dpi=180); plt.close(fig)
    fig, axes = plt.subplots(len(chosen), 2, figsize=(14, 1.8 * len(chosen)))
    for pair, run in zip(axes, chosen):
        raw = np.asarray(Image.open(DATA / run / 'ImgA000250.jpg').convert('L'))
        panorama = np.asarray(Image.open(panorama_path(run)).convert('L'))
        pair[0].imshow(raw, cmap='gray', vmin=0, vmax=255, aspect='auto')
        pair[1].imshow(panorama[:, :2500], cmap='gray', vmin=0, vmax=255, aspect='auto')
        pair[0].set_ylabel(run, fontsize=9)
        for ax in pair:
            ax.set_xticks([]); ax.set_yticks([0, 40, 79])
    axes[0, 0].set_title('Raw frame ImgA000250.jpg')
    axes[0, 1].set_title('Saved panorama, first up to 2500 columns')
    fig.suptitle('Retained acquisitions: two different spatial records')
    fig.tight_layout(); fig.savefig(artifact_path('raw_comparison.png'), dpi=180); plt.close(fig)


if __name__ == "__main__":
    main()
