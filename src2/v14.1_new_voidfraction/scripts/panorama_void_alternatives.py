"""Current optical model: every intensity operation is confined to rows 18:72.

The earlier D53 model and its supplemental analyses are preserved in
intermediate/history/before_roi_unification_20261030. This module implements
conditional optical geometry; it does not calibrate the physical gas volume.
"""
from pathlib import Path
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter1d, binary_closing, binary_fill_holes, label
from v14_geometry import LUMEN_ROWS, LUMEN_HEIGHT_PX, DIAMETER_PX

REFERENCE_RUNS = ("m05q01_001", "m06q01_001", "m07q05_001")
THRESHOLD = .70


def read(run):
    project = Path(__file__).resolve().parents[3]
    path = project / "src2/v18_montagem_18_71/01_panoramas/optica" / run / "panorama.png"
    return np.asarray(Image.open(path).convert("L"), dtype=float)


def crop_roi(im):
    """Discard all out-of-lumen intensities before any numerical operation."""
    array = np.asarray(im)
    if array.ndim != 2 or array.shape[0] < LUMEN_ROWS[1] or array.shape[1] < 1:
        raise ValueError("Expected a panorama containing rows 18 through 71 inclusive")
    return array[LUMEN_ROWS[0]:LUMEN_ROWS[1]].astype(float, copy=True)


def roi_saturation(im):
    return float(np.mean(crop_roi(im) == 255))


def normalize(im):
    roi = crop_roi(im)
    return roi / np.maximum(np.percentile(roi, 75, axis=0), 1)[None, :]


def make_reference(names=REFERENCE_RUNS, percentile=.05):
    profiles, selections = [], {}
    for run in names:
        norm = normalize(read(run))
        score = norm.std(axis=0)
        cutoff = np.quantile(score, percentile)
        selected = score <= cutoff
        profiles.append(np.median(norm[:, selected], axis=1))
        selections[run] = dict(n_selected=int(selected.sum()), score_cutoff=float(cutoff),
                               selected_position_px=np.flatnonzero(selected).tolist())
    return gaussian_filter1d(np.median(profiles, axis=0), 1), selections


def close_shadow_crop(crop):
    # Zero padding preserves actual contact with the observed radial boundary.
    padded = np.pad(np.asarray(crop, dtype=bool), 1, mode="constant")
    return binary_closing(padded, structure=np.ones((3, 3)))[1:-1, 1:-1]


def silhouette(im, reference, threshold=THRESHOLD, rows=LUMEN_ROWS, min_open_interval=20):
    if tuple(rows) != LUMEN_ROWS:
        raise ValueError("The current optical analysis uses the single full ROI18:72")
    reference = np.asarray(reference)
    if reference.shape != (LUMEN_HEIGHT_PX,):
        raise ValueError("The liquid reference must contain exactly the 54 ROI rows")
    ratio = gaussian_filter1d(normalize(im), 1, axis=0) / np.maximum(reference[:, None], .01)
    height, width = ratio.shape
    midpoint = height // 2
    cols = np.arange(width)
    upper = np.argmin(ratio[:midpoint], axis=0)
    lower = midpoint + np.argmin(ratio[midpoint:], axis=0)
    accepted = (ratio[upper, cols] < threshold) & (ratio[lower, cols] < threshold)
    shadow = ratio < threshold
    closed = close_shadow_crop(shadow)
    components, count = label(closed)
    top_component, bottom_component = components[upper, cols], components[lower, cols]
    same_outline = (top_component > 0) & (top_component == bottom_component)
    axial_components, _ = label(accepted)
    axial_sizes = np.bincount(axial_components)
    accepted &= same_outline | (axial_sizes[axial_components] >= min_open_interval)
    outer_upper, outer_lower = upper.copy(), lower.copy()
    for _ in range(height):
        advance = (outer_upper > 0) & shadow[np.maximum(outer_upper - 1, 0), cols]
        outer_upper[advance] -= 1
        advance = (outer_lower < height - 1) & shadow[np.minimum(outer_lower + 1, height - 1), cols]
        outer_lower[advance] += 1
    core = np.zeros_like(shadow)
    for y in range(height):
        core[y] = accepted & (y >= outer_upper) & (y <= outer_lower)
    rejected = np.unique(np.concatenate([components[0], components[-1]]))
    sizes = np.bincount(components.ravel(), minlength=count + 1)
    retained = (~np.isin(components, rejected)) & (sizes[components] >= 4)
    core |= binary_fill_holes(retained)
    out = np.zeros_like(im, dtype=bool)
    out[LUMEN_ROWS[0]:LUMEN_ROWS[1]] = core
    return out, accepted, ratio


def geometry(mask, diameter=DIAMETER_PX):
    # Each contiguous radial interval is treated as a circular gas section.
    # Circular gas sections/no overlapping depth are additional assumptions.
    core = crop_roi(mask).astype(bool)
    counts = np.zeros(core.shape[1], dtype=float)
    square_sum = np.zeros(core.shape[1], dtype=float)
    for active in core:
        square_sum += np.where(~active, counts ** 2, 0)
        counts = np.where(active, counts + 1, 0)
    square_sum += counts ** 2
    return square_sum / diameter ** 2


if __name__ == "__main__":
    from run_all import main
    raise SystemExit(main())
