"""Panorama gas silhouettes using inferred liquid intervals and shadow envelopes.

The same segmentation rule is used for every selected acquisition. Liquid
reference candidates are the least radially structured 5% of columns in three
visually reviewed plug panoramas. The dispersed-bubble panorama is held out from
reference construction; it was inspected during method development and is not
an independent validation dataset. Report apparent geometric reconstructions
conditionally, not calibrated physical void.
"""
from v14_paths import ROOT as HERE, DATA, CACHE, artifact_path, write_markdown
from v14_selection import iter_runs, is_eligible, PATTERNS
from pathlib import Path
import csv
import json
import os
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter1d, binary_closing, binary_fill_holes, label

os.environ.setdefault("MPLCONFIGDIR", str(CACHE / "matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REFERENCE_RUNS = ("m05q01_001", "m06q01_001", "m07q05_001")
THRESHOLD = .70
DIAMETER_PX = 60.0


def read(run):
    return np.asarray(Image.open(DATA / run / "__panorama_v11.png").convert("L"), dtype=float)


def normalize(im):
    # Gas and liquid centers can have similar brightness. Normalize only to
    # compare radial shape; this is not a Beer-Lambert volume calibration.
    return im / np.maximum(np.percentile(im[30:50], 75, axis=0), 1)[None, :]


def make_reference(names=REFERENCE_RUNS, percentile=.05):
    profiles, selections = [], {}
    for run in names:
        im = read(run)
        norm = normalize(im)
        score = norm[20:68].std(axis=0)
        selected = score <= np.quantile(score, percentile)
        profiles.append(np.median(norm[:, selected], axis=1))
        selections[run] = dict(n_selected=int(selected.sum()),
                               score_cutoff=float(np.quantile(score, percentile)),
                               selected_position_px=np.flatnonzero(selected).tolist())
    return gaussian_filter1d(np.median(profiles, axis=0), 1), selections


def close_shadow_crop(crop):
    """Close on the zero-extended plane, then restore the observed crop.

    Closing a finite array with SciPy's default border value removes its outer
    row/column during erosion. Padding by the 3x3 element's radius preserves
    real boundary contact, which the radial-component rejection below needs.
    Axial record-edge contact remains allowed.
    """
    padded = np.pad(np.asarray(crop, dtype=bool), 1, mode='constant')
    return binary_closing(padded, structure=np.ones((3, 3)))[1:-1, 1:-1]


def silhouette(im, reference, threshold=THRESHOLD, rows=(18, 71), min_open_interval=20):
    ratio = gaussian_filter1d(normalize(im), 1, axis=0) / np.maximum(reference[:, None], .01)
    cols = np.arange(im.shape[1])
    upper = 19 + np.argmin(ratio[19:38], axis=0)
    lower = 42 + np.argmin(ratio[42:71], axis=0)
    accepted = (ratio[upper, cols] < threshold) & (ratio[lower, cols] < threshold)
    shadow = ratio < threshold
    crop = shadow[rows[0]:rows[1]]
    closed = close_shadow_crop(crop)
    components, count = label(closed)
    top_component = components[np.clip(upper - rows[0], 0, crop.shape[0] - 1), cols]
    bottom_component = components[np.clip(lower - rows[0], 0, crop.shape[0] - 1), cols]
    same_outline = (top_component > 0) & (top_component == bottom_component)
    axial_components, _ = label(accepted)
    axial_sizes = np.bincount(axial_components)
    # A short paired deficit from different objects can bridge liquid between
    # unrelated bubbles. Keep closed outlines, and keep long open interfaces
    # such as an annular core. This rule is global and is never label-tuned.
    accepted &= same_outline | (axial_sizes[axial_components] >= min_open_interval)
    outer_upper, outer_lower = upper.copy(), lower.copy()
    for _ in range(32):
        advance = (outer_upper > rows[0]) & shadow[np.maximum(outer_upper - 1, 0), cols]
        outer_upper[advance] -= 1
        advance = (outer_lower < rows[1] - 1) & shadow[np.minimum(outer_lower + 1, 79), cols]
        outer_lower[advance] += 1
    out = np.zeros_like(im, dtype=bool)
    for y in range(rows[0], rows[1]):
        out[y] = accepted & (y >= outer_upper) & (y <= outer_lower)
    # Small/off-center shadow candidates: retain components away from radial
    # crop edges and fill any enclosed holes. A >=4-pixel component need not
    # enclose a hole, so unresolved dark spots are candidates too. This rule
    # can miss true bubbles connected to a crop edge and retain liquid shadows.
    # Lateral record-edge contact is permitted.
    rejected = np.unique(np.concatenate([components[0], components[-1]]))
    sizes = np.bincount(components.ravel(), minlength=count + 1)
    retained = (~np.isin(components, rejected)) & (sizes[components] >= 4)
    filled = binary_fill_holes(retained)
    out[rows[0]:rows[1]] |= filled
    return out, accepted, ratio


def geometry(mask, diameter=DIAMETER_PX):
    # Every separate vertical gas-silhouette interval is treated as a disk
    # cross-section. This is conditional on axisymmetry/no overlapping depth.
    counts = np.zeros(mask.shape[1], dtype=float)
    square_sum = np.zeros(mask.shape[1], dtype=float)
    for y in range(mask.shape[0]):
        active = mask[y]
        square_sum += np.where(~active, counts ** 2, 0)
        counts = np.where(active, counts + 1, 0)
    square_sum += counts ** 2
    return square_sum / diameter ** 2


def emit_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main():
    reference, selections = make_reference()
    top_wall = int(8 + np.argmin(reference[8:22]))
    bottom_wall = int(65 + np.argmin(reference[65:80]))
    (artifact_path("panorama_liquid_reference.json")).write_text(json.dumps(dict(
        liquid_candidates=selections, radial_reference=reference.tolist(),
        apparent_wall_minima=[top_wall, bottom_wall],
        apparent_wall_separation_px=bottom_wall - top_wall,
        independent_negative_control="m07q01_001",
        control_status="Legacy field name: held out from reference construction, but inspected during detector development; not independent validation",
        status="Inferred liquid-like reference; not a separately acquired liquid-only calibration"), indent=2) + "\n")
    alternatives = {"baseline": reference}
    for removed in REFERENCE_RUNS:
        alternatives["without_" + removed] = make_reference(tuple(r for r in REFERENCE_RUNS if r != removed))[0]
    alternatives["reference_1percent"] = make_reference(percentile=.01)[0]
    alternatives["reference_10percent"] = make_reference(percentile=.10)[0]
    control_rows = []
    for run in REFERENCE_RUNS:
        heldout_reference = alternatives["without_" + run]
        control_mask, _, _ = silhouette(read(run), heldout_reference)
        selected = selections[run]["selected_position_px"]
        control_rows.append(dict(run=run, n_liquid_like_candidate_columns=len(selected),
                                 heldout_reference_projected_gas_fraction=float(control_mask[:, selected].sum() / (DIAMETER_PX * len(selected))),
                                 heldout_reference_any_gas_column_fraction=float(control_mask[:, selected].any(axis=0).mean()),
                                 heldout_reference_conditional_alpha=float(geometry(control_mask)[selected].mean())))
    emit_csv(artifact_path("panorama_liquid_holdout_checks.csv"), control_rows)
    rows, sensitivity = [], []
    for folder in iter_runs():
        run = folder.name; im = read(run)
        marker = [p.stem.lstrip("_") for p in folder.glob("__*.txt")]
        pattern = marker[0] if marker else "unknown"
        mask, accepted, _ = silhouette(im, reference)
        alpha = geometry(mask)
        radial_clipped = mask[18] | mask[70]
        # Edge-contact extension scenario: retain unflagged columns and allow
        # flagged columns to extend up to a full tube. This omits gas wholly
        # unseen outside the crop and rejected/missed objects; it is a bound
        # only with the additional completeness assumption, correct phase
        # assignment and undistorted circular geometry.
        missing_width_upper = alpha + radial_clipped * (1 - alpha)
        sat = float((im[20:60] >= 255).mean())
        row = dict(run=run, label=pattern, panorama_width_px=im.shape[1],
                   gas_silhouette_projected_fraction=float(mask.sum() / (DIAMETER_PX * im.shape[1])),
                   conditional_circle_alpha_spatial_mean=float(alpha.mean()),
                   conditional_circle_alpha_spatial_median=float(np.median(alpha)),
                   detected_center_gas_column_fraction=float(accepted.mean()),
                   center_gas_silhouette_duty=float(mask[40].mean()),
                   radial_crop_touch_column_fraction=float(radial_clipped.mean()),
                   conditional_missing_width_alpha_lower=float(alpha.mean()),
                   conditional_missing_width_alpha_upper=float(missing_width_upper.mean()),
                   core_saturated_fraction=sat,
                   validity="unsupported_saturation" if sat > .20 else "conditional_optical_geometry",
                   conditional_alpha_D56=float(geometry(mask, 56).mean()),
                   conditional_alpha_D64=float(geometry(mask, 64).mean()))
        values = []
        for threshold in (.60, .65, .70, .75, .80):
            smask, saccepted, _ = silhouette(im, reference, threshold=threshold)
            value = float(geometry(smask).mean()); values.append(value)
            sensitivity.append(dict(run=run, variant=f"shadow_ratio_{threshold:.2f}",
                                    alpha_conditional=value, center_gas_fraction=float(saccepted.mean())))
        row["threshold_060_080_min"] = min(values)
        row["threshold_060_080_max"] = max(values)
        refs = []
        for name, ref in alternatives.items():
            smask, saccepted, _ = silhouette(im, ref)
            value = float(geometry(smask).mean()); refs.append(value)
            sensitivity.append(dict(run=run, variant=name, alpha_conditional=value,
                                    center_gas_fraction=float(saccepted.mean())))
        row["reference_variants_min"] = min(refs)
        row["reference_variants_max"] = max(refs)
        for span in (10, 30):
            smask, saccepted, _ = silhouette(im, reference, min_open_interval=span)
            sensitivity.append(dict(run=run, variant=f"open_interface_min_span_{span}",
                                    alpha_conditional=float(geometry(smask).mean()),
                                    center_gas_fraction=float(saccepted.mean())))
        for bounds in ((19, 70), (20, 70)):
            smask, saccepted, _ = silhouette(im, reference, rows=bounds)
            sensitivity.append(dict(run=run, variant=f"crop_{bounds[0]}_{bounds[1]}",
                                    alpha_conditional=float(geometry(smask).mean()),
                                    center_gas_fraction=float(saccepted.mean())))
        blocks = np.array_split(alpha, 5)
        row["five_spatial_block_means_min"] = float(min(x.mean() for x in blocks))
        row["five_spatial_block_means_max"] = float(max(x.mean() for x in blocks))
        rows.append(row)
        print(run, pattern, "conditional alpha", round(alpha.mean(), 4),
              "threshold range", round(min(values), 4), round(max(values), 4), flush=True)
    emit_csv(artifact_path("panorama_void_per_run.csv"), rows)
    emit_csv(artifact_path("panorama_void_sensitivity.csv"), sensitivity)
    pure = PATTERNS
    group = []
    for pattern in pure:
        values = [r["conditional_circle_alpha_spatial_mean"] for r in rows
                  if r["label"].replace("_", " ") == pattern
                  and is_eligible(r["run"], pattern, r["core_saturated_fraction"])]
        group.append(dict(pattern=pattern, n=len(values),
                          conditional_alpha_median=float(np.median(values)),
                          conditional_alpha_min=min(values), conditional_alpha_max=max(values)))
    emit_csv(artifact_path("panorama_void_by_pattern.csv"), group)
    report = ["# Panorama shadow-envelope reconstruction", "",
              f"The current scope contains {len(rows)} acquisitions across {len(PATTERNS)} retained patterns, with {sum(g["n"] for g in group)} acquisitions passing the saturation criterion in grouped statistics. Mixed labels, smooth annular and the user-excluded replicate are omitted before processing. Supplied retained labels are unchanged.", "",
              "**A usable image solution is possible as a conditional geometric estimate.** This method identifies dark gas outlines against a liquid-like radial template, fills the gas interior, and includes separately enclosed off-center bubbles. It does not equate dark pixels alone, or bright pixels alone, with gas. The result is a physically meaningful reconstruction formula applied to an optical silhouette; its absolute calibration remains conditional.", "",
              "## Reference inferred from available panoramas", "",
              "The three visually reviewed plug panoramas `m05q01_001`, `m06q01_001`, `m07q05_001` provide liquid-looking candidate intervals. In each, normalize radial profiles by the 75th percentile of rows 30:50; select columns at or below the 5th percentile of radial standard deviation in rows 20:68. Median selected profiles give the liquid-like template. Candidate columns are selected from profile structure; the three reference panoramas were chosen by inspection. `m07q01_001` is excluded from reference construction but was inspected during detector development, so it is an internal check rather than independent validation. The JSON saves exact candidate column coordinates.", "",
              f"The inferred reference has wall/reflection minima at image rows **{top_wall} and {bottom_wall}**, an apparent separation of **{bottom_wall-top_wall} px**. We use D_app=60 px, and repeat with 56 and 64 px. These are apparent optical dimensions: the minima have not been validated as physical lumen edges.", "",
              "## Segmentation and physical reconstruction", "",
              "Smooth normalized radial profiles by σ=1 px and divide by the inferred template. Locate paired upper/lower shadow deficits below ratio 0.70, expand each contiguous deficit to its outer shadow edge, and fill between them. Reject paired fills from different connected shadow components when their accepted axial interval is shorter than 20 px; this reduces false bridges but does not prove phase identity. Long open interfaces remain eligible. Separately, apply 3×3 morphological closing on a zero-padded crop, restore the crop, reject components touching a radial edge, retain components of at least 4 pixels, and fill enclosed holes. Axial record-edge contact is allowed. Components need not contain a hole: unresolved dark spots remain gas candidates, and genuine bubbles connected to a radial edge can be missed. Identical parameters apply to every panorama.", "",
              "The projected gas-silhouette fraction is A_g,projected/(D_app L_record). The **conditional circular void estimate** is E_z[Σ_i(d_i(z)/D_app)²], where d_i is each separate silhouette's radial width. This follows from A_g(z)=Σ_i π d_i²/4 and A_tube=πD²/4. It requires circular/axisymmetric gas cross-sections and negligible overlap in viewing depth. The arithmetic mean of squares is used, not the median width or square of its mean. The px→mm factor cancels in the ratio.", "",
              "## Per-run estimates and sensitivity", "",
              "These numbers are **conditional optical reconstructions**, not independently calibrated measurements. Saturated runs are retained and explicitly marked unsupported. The threshold interval is a perturbation range, not a statistical confidence interval.", "",
              "| Run | Pattern | Conditional α | Edge-contact extension scenario | Threshold 0.60–0.80 range | Reference variant range | Fixed-row y=40 gas duty | Status |",
              "|---|---|---:|---:|---:|---:|---:|---|"]
    for r in rows:
        report.append(f"| {r['run']} | {r['label']} | {r['conditional_circle_alpha_spatial_mean']:.3f} | {r['conditional_missing_width_alpha_lower']:.3f}–{r['conditional_missing_width_alpha_upper']:.3f} | {r['threshold_060_080_min']:.3f}–{r['threshold_060_080_max']:.3f} | {r['reference_variants_min']:.3f}–{r['reference_variants_max']:.3f} | {r['center_gas_silhouette_duty']:.3f} | {r['validity']} |")
    report += ["", "## Adversarial checks", "",
               "The sensitivity CSV repeats thresholds 0.60/0.65/0.70/0.75/0.80, leaves out each of the three reference panoramas in turn, changes reference candidate selection to 1% and 10%, narrows the radial crop, and changes the open-interface minimum interval from 20 to 10/30 px. The dispersed panorama is never a reference. Five spatial-block means expose within-record variability. These checks measure internal robustness; they cannot substitute for independent void calibration.", "",
               "Liquid-like candidate checks are saved separately, with reference construction excluding the tested panorama. Projected foreground fractions: " + ", ".join(f"{r['run']} {100*r['heldout_reference_projected_gas_fraction']:.4f}%" for r in control_rows) + ". Conditional area contributions: " + ", ".join(f"{r['run']} {r['heldout_reference_conditional_alpha']:.7f}" for r in control_rows) + ". The candidates were inferred from the same dataset and are not certified all-liquid ground truth.", "",
               "The overlay figure makes possible false phase assignments reviewable. The dispersed internal check should show small bubbles with little foreground between them. Plug/slug images should show both shadow boundaries and bright interiors, with liquid slugs empty. Annular should keep the connected core, while churn retains uncertain fragmented outlines. Fixed-row duty uses image row y=40, not a validated tube center; the midpoint between apparent wall minima y=15 and y=75 is y=45.", "",
               "**A censoring limit matters for annular:** filling all 53 rows in the interior crop imposes a computational maximum (53/60)²=0.7803. This follows from the chosen crop and diameter and is not a physical maximum of void fraction. Several annular runs approach it; apparent agreement or low threshold scatter must not be mistaken for precision. The legacy CSV columns `conditional_missing_width_alpha_lower/upper` now explicitly denote an edge-contact extension scenario: retain unflagged columns and let each flagged column extend up to a full tube, E[α_visible+1_edge(1−α_visible)]. This is only an upper bound if no gas is wholly unseen outside the crop and no gas objects are rejected or missed, in addition to correct phase labels, circular geometry and an undistorted diameter. A bubble entirely outside the crop can increase real void while leaving both scenario values unchanged. The interval is not a confidence interval or a general bound on total gas.", "",
               "## Remaining physical assumptions", "",
               "A single corrected side view can reconstruct an axisymmetric gas body; a second view and absolute mm scale are not invariably necessary. The present cylindrical tube's refraction can move shadow edges away from physical interface positions. The glass/illumination inverse map remains unknown, and irregular churn or overlapping dispersed bubbles can violate circular reconstruction. Those are the principal limits on calling the numerical result the true void fraction.", "",
               "A stitched panorama supplies a spatial record. Its mean is spatial occupancy under the panorama reconstruction's frozen-advection assumption; it is not automatically a time average at a fixed point. A physical time average requires native frame timestamps and fixed-column phase masks. No conversion of panorama length to independent sample count or seconds is made here.", "",
               "Source: Schepperle et al. (2024), Sensors 24, 3363, DOI 10.3390/s24113363, supports contour extraction and comparison against liquid appearance (§2.3); its experiment is rectangular and uses a separate liquid-only photograph. Our inferred reference and circular reconstruction are extensions and must be reported as such."]
    write_markdown(artifact_path("panorama_void_alternatives.md"), "\n".join(report) + "\n")
    selected = ["m07q01_001", "m06q01_001", "m06q05_001", "m06q20_001", "m06q80_001"]
    fig, axes = plt.subplots(len(selected), 2, figsize=(15, 2 * len(selected)))
    for pair, run in zip(axes, selected):
        im = read(run); mask, _, _ = silhouette(im, reference)
        stop = min(im.shape[1], 2500)
        pair[0].imshow(im[:, :stop], cmap="gray", vmin=0, vmax=255, aspect="auto")
        pair[1].imshow(im[:, :stop], cmap="gray", vmin=0, vmax=255, aspect="auto")
        color = np.zeros((80, stop, 4)); color[..., 0]=.1; color[..., 1]=.8; color[..., 2]=1
        color[..., 3]=mask[:, :stop]*.42
        pair[1].imshow(color, aspect="auto")
        pair[0].set_ylabel(run, fontsize=8)
        for ax in pair: ax.set_yticks([0, 40, 79]); ax.set_xticks([])
    axes[0, 0].set_title("Panorama original (primeiros ≤2.500 px)")
    axes[0, 1].set_title("Gás estimado em azul")
    fig.suptitle("Segmentação dos panoramas com referência líquida inferida")
    fig.tight_layout(); fig.savefig(artifact_path("panorama_phase_masks.png"), dpi=180); plt.close(fig)
    fig, axes = plt.subplots(len(REFERENCE_RUNS), 1, figsize=(12, 5))
    for ax, run in zip(axes, REFERENCE_RUNS):
        im = read(run); selection = selections[run]["selected_position_px"]
        ax.imshow(im, cmap="gray", vmin=0, vmax=255, aspect="auto")
        ax.scatter(selection, np.full(len(selection), 40), s=2, color="cyan")
        ax.set_ylabel(run, fontsize=8); ax.set_yticks([0,40,79])
    fig.suptitle("Colunas de referência líquida inferida, em ciano")
    fig.tight_layout(); fig.savefig(artifact_path("panorama_liquid_reference_candidates.png"), dpi=180); plt.close(fig)


if __name__ == "__main__":
    main()
