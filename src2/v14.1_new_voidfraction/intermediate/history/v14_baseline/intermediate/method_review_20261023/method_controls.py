"""Read-only production audit: analytic controls and boundary/probe sensitivity."""
import csv
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'scripts'))
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / 'intermediate/cache/matplotlib'))
import numpy as np
from PIL import Image
from scipy.ndimage import binary_closing, binary_fill_holes, gaussian_filter1d, label
from panorama_void_alternatives import geometry, make_reference, normalize, read, silhouette
from v14_selection import iter_runs


def boundary_preserving_silhouette(im, reference, threshold=.70, rows=(18, 71), min_open_interval=20):
    """Identical estimator except zero-padded closing before component tests."""
    ratio = gaussian_filter1d(normalize(im), 1, axis=0) / np.maximum(reference[:, None], .01)
    cols = np.arange(im.shape[1])
    upper = 19 + np.argmin(ratio[19:38], axis=0)
    lower = 42 + np.argmin(ratio[42:71], axis=0)
    accepted = (ratio[upper, cols] < threshold) & (ratio[lower, cols] < threshold)
    shadow = ratio < threshold
    crop = shadow[rows[0]:rows[1]]
    # Closing on the zero-extended image, then restrict to the observed crop.
    # Unlike closing on the cropped finite array, this preserves its boundary.
    closed = binary_closing(np.pad(crop, 1), structure=np.ones((3, 3)))[1:-1, 1:-1]
    components, count = label(closed)
    top_component = components[np.clip(upper - rows[0], 0, crop.shape[0] - 1), cols]
    bottom_component = components[np.clip(lower - rows[0], 0, crop.shape[0] - 1), cols]
    same_outline = (top_component > 0) & (top_component == bottom_component)
    axial_components, _ = label(accepted)
    axial_sizes = np.bincount(axial_components)
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
    rejected = np.unique(np.concatenate([components[0], components[-1]]))
    sizes = np.bincount(components.ravel(), minlength=count + 1)
    retained = (~np.isin(components, rejected)) & (sizes[components] >= 4)
    out[rows[0]:rows[1]] |= binary_fill_holes(retained)
    return out


def flips(arr):
    return int(np.count_nonzero(arr[..., 1:] != arr[..., :-1]))


def main():
    reference, _ = make_reference()
    measurements = []
    originals = {}
    alternatives = {}
    for folder in iter_runs():
        im = read(folder.name)
        mask, _, ratio = silhouette(im, reference)
        alt = boundary_preserving_silhouette(im, reference)
        originals[folder.name], alternatives[folder.name] = mask, alt
        crop = ratio[18:71] < .70
        closed = binary_closing(crop, structure=np.ones((3, 3)))
        comp, _ = label(closed)
        row = dict(run=folder.name, pattern=next(folder.glob('__*.txt')).stem.lstrip('_'),
                   raw_radial_edge_foreground=int(crop[[0,-1]].sum()),
                   closed_radial_edge_foreground=int(closed[[0,-1]].sum()),
                   n_rejected_original_foreground_components=int(np.count_nonzero(np.unique(np.r_[comp[0],comp[-1]]))),
                   alpha_original=float(geometry(mask).mean()),
                   alpha_padded_closing=float(geometry(alt).mean()),
                   alpha_change=float(geometry(alt).mean()-geometry(mask).mean()),
                   changed_mask_pixels=int((mask != alt).sum()),
                   flip_original=flips(mask[18:71])/(53*(im.shape[1]-1)),
                   flip_padded_closing=flips(alt[18:71])/(53*(im.shape[1]-1)),
                   flip_row40_per1000=1000*flips(mask[40])/(im.shape[1]-1),
                   flip_row45_per1000=1000*flips(mask[45])/(im.shape[1]-1),
                   duty_row40=float(mask[40].mean()), duty_row45=float(mask[45].mean()))
        measurements.append(row)
    with (OUT/'boundary_and_probe_sensitivity.csv').open('w',newline='') as fp:
        writer=csv.DictWriter(fp, fieldnames=list(measurements[0])); writer.writeheader(); writer.writerows(measurements)
    # A mathematically valid gas disk lying wholly outside the observed crop.
    full=np.zeros((80,10),bool); full[15:18]=True
    visible=np.zeros_like(full); visible[18:71]=full[18:71]
    seen=geometry(visible); clipped=visible[18] | visible[70]
    upper=seen+clipped*(1-seen)
    controls=dict(
        hidden_disk=dict(diameter_px=3, full_alpha=float(geometry(full).mean()),
                         observed_alpha=float(seen.mean()), reported_upper=float(upper.mean()),
                         claimed_interval_covers_true=float(geometry(full).mean())<=float(upper.mean())),
        all_original_boundary_checks_reject_no_foreground=all(r['n_rejected_original_foreground_components']==0 for r in measurements),
        runs_changed_by_boundary_fix=sum(r['changed_mask_pixels']>0 for r in measurements),
        largest_abs_alpha_change=max(measurements,key=lambda r:abs(r['alpha_change'])),
        probe_row45=dict(inferred_wall_rows=[15,75], apparent_midpoint=45, implemented_probe=40),
    )
    (OUT/'analytic_controls.json').write_text(json.dumps(controls,indent=2)+'\n')
    # Keep reviewable mask-change previews, clearly marked diagnostic alternatives.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    selected=['m07q01_001','m06q01_001','m06q05_001','m06q20_001','m06q80_001']
    fig, axes=plt.subplots(len(selected),3,figsize=(18,10))
    for axrow, run in zip(axes,selected):
        im=read(run); stop=min(im.shape[1],2000)
        for ax in axrow:
            ax.imshow(im[:,:stop],cmap='gray',vmin=0,vmax=255,aspect='auto'); ax.set_xticks([]); ax.set_yticks([0,40,79])
        axrow[0].set_ylabel(run)
        for ax, mask in zip(axrow[1:],[originals[run],alternatives[run]]):
            rgba=np.zeros((80,stop,4)); rgba[:,:,1:3]=1; rgba[:,:,3]=mask[:,:stop]*.40
            ax.imshow(rgba,aspect='auto')
    for ax,title in zip(axes[0],['Original input','Production mask','Diagnostic padded-closing alternative']): ax.set_title(title)
    fig.tight_layout(); fig.savefig(OUT/'boundary_fix_diagnostic.png',dpi=140); plt.close(fig)
    print(json.dumps(controls,indent=2))


if __name__=='__main__': main()
