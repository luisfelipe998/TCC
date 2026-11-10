#!/usr/bin/env python3
"""Verify and report the targeted contour restriction and object-track probe."""
from pathlib import Path
from collections import defaultdict
import csv
import json
import os
import numpy as np
from PIL import Image

BASE=Path(__file__).resolve().parents[1]
OUT=BASE/'results'
os.environ.setdefault('MPLCONFIGDIR',str(BASE/'intermediate'/'closed_contour_mplcache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


def read(name):return list(csv.DictReader((OUT/name).open()))


def run():
    summary=read('closed_contour_summary.csv');frames=read('closed_contour_frames.csv')
    links=read('closed_contour_bubble_links.csv');bsummary=read('closed_contour_bubble_summary.csv')
    numeric=0
    for row in summary:
        selected=[r for r in frames if r['run']==row['run']]
        suffix=f"{round(float(row['threshold'])*100):03d}"
        for key,source in [('baseline_projected_gas_fraction',f'masks_{suffix}_baseline_fraction'),
                           ('closed_projected_gas_fraction',f'masks_{suffix}_closed_fraction')]:
            assert np.isclose(float(row[key]),np.mean([float(r[source]) for r in selected]),atol=1e-12,rtol=0)
            numeric+=1
        assert int(row['added_pixels'])==0
    for row in bsummary:
        current=[r for r in links if r['threshold']==row['threshold']]
        good=[r for r in current if r['stable_chain']=='True']
        assert len(current)==int(row['matched_intervals']) and len(good)==int(row['stable_intervals'])
        values=np.array([float(r['dx_px']) for r in good])
        for key,value in [('mean_px_frame',values.mean()),('median_px_frame',np.median(values)),
                          ('median_px_s',np.median(values)*7000),('q25_px_frame',np.quantile(values,.25)),
                          ('q75_px_frame',np.quantile(values,.75))]:
            assert np.isclose(float(row[key]),value,atol=1e-12,rtol=0);numeric+=1
        chains=defaultdict(list)
        for r in current:
            assert np.isclose(float(r['dx_px']),float(r['next_x'])-float(r['x']))
            assert np.isclose(float(r['dy_px']),float(r['next_y'])-float(r['y']))
            assert np.isclose(float(r['velocity_px_s']),float(r['dx_px'])*7000)
            numeric+=3;chains[int(r['chain_id'])].append(r)
        for rows in chains.values():
            x=np.array([float(r['dx_px']) for r in rows]);y=np.array([float(r['dy_px']) for r in rows])
            passed=bool(len(rows)>=5 and np.max(np.abs(np.diff(x)))<=2 and np.max(np.abs(np.diff(y)))<=2)
            assert all((r['stable_chain']=='True')==passed for r in rows)
    groups=defaultdict(list)
    for r in links:
        if r['threshold']=='070' and r['stable_chain']=='True':groups[int(r['chain_id'])].append(r)
    chosen=[k for k,rows in groups.items() if min(int(r['frame_id']) for r in rows)<=200 and max(int(r['frame_id']) for r in rows)>=210]
    selected_frames=(200,205,210)
    fig,axes=plt.subplots(3,1,figsize=(14,4))
    colors=('#00e5ff','#ff5b5b','#ffdf00')
    data=BASE.parent.parent/'data_processed'/'m07q01_001'
    files=sorted(data.glob('ImgA*.jpg'))
    for frame,ax in zip(selected_frames,axes):
        raw=np.asarray(Image.open(files[frame]).convert('L'))
        ax.imshow(raw,cmap='gray',vmin=0,vmax=255)
        for c,chain in enumerate(chosen):
            r=next(r for r in groups[chain] if int(r['frame_id'])==frame)
            x,y=float(r['x']),float(r['y'])
            ax.add_patch(Circle((x,y),10,fill=False,color=colors[c%3],lw=1))
            ax.text(x+12,y-10,f'candidate {chain}',color=colors[c%3],fontsize=8)
        ax.set_title(f'Frame {frame}; {(frame-200)/7000*1000:.3f} ms after frame 200',fontsize=9)
        ax.set_xlim(0,1015);ax.set_ylim(79,0);ax.axis('off')
    fig.suptitle('Selected small-object tracks after contour restriction (threshold 0.70)\nCentroid correspondence illustrates feasibility; no liquid motion or phase ground truth',fontsize=11)
    fig.tight_layout(rect=(0,0,1,.88))
    fig.savefig(BASE/'figures'/'closed_contour_bubble_tracks.png',dpi=180)
    plt.close(fig)
    m07=[r for r in summary if r['run']=='m07q01_001']
    lines=['# Targeted closed-contour and dispersed-object probe','',
        '## Executive finding','',
        '**Disabling the long open-envelope fill partially repairs a conspicuous native-frame dispersed-bubble artifact.** It removes the large apparent right-edge gas block in the inspected frames, while preserving several recognizable small bubble candidates. Residual reflection/shadow candidates and missed bubbles remain. The same rule removes large amounts of legitimate open/cropped gas structure in other patterns, so it is not a universal replacement for v14.',
        '',
        'An exploratory small-object tracker then produces persistent bubble-candidate motion near −6.1 to −6.3 px/frame in the dispersed recording. This supports a pattern-specific object-tracking path. It does not provide liquid parcel motion or a gas/liquid material velocity ratio.',
        '',
        '## Controlled topology change','',
        'The current v14 silhouette function and frozen liquid-like reference were retained. Only min_open_interval was increased from 20 to 1,000,000 pixels, disabling acceptance solely because an open paired-shadow interval is long. Same-outline acceptance and the separate closed-component branch remain. Thresholds 0.60, 0.70 and 0.80 are tested. All 501 frames of m07q01_001 and 20 fixed equally spaced frames in each of the other 25 retained recordings were processed: 1,001 frames and 3,003 detector evaluations.',
        '',
        'These are identical-frame baseline/alternative comparisons, not comparisons against phase truth. Every alternative mask is asserted to be a subset of its baseline mask. Lower area alone is not evidence of correctness.',
        '',
        '| Threshold | Baseline projected gas (%) | Open-fill-disabled projected gas (%) | Baseline gas pixels removed (%) |',
        '|---|---:|---:|---:|']
    for r in m07:
        lines.append(f"| {float(r['threshold']):.2f} | {100*float(r['baseline_projected_gas_fraction']):.3f} | {100*float(r['closed_projected_gas_fraction']):.3f} | {100*float(r['removed_fraction_of_baseline']):.2f} |")
    lines += ['',
        'Across the 501 dispersed frames, threshold disagreement decreases from **24.369% to 4.360% of lumen pixels**. This is improved internal consistency, not segmentation accuracy.',
        '',
        '![Dispersed-frame topology comparison](../figures/closed_contour_m07q01_frame200.png)',
        '',
        'In inspected frames 200 and 210, the large right-edge block disappears. Thin lower reflection/shadow structures remain labeled at higher thresholds. Several obvious isolated bubble outlines remain, while some neighboring visible bubbles are missed or only partly outlined. No manual pixel-level ground truth was constructed.',
        '',
        '## Why this requires pattern-specific handling','',
        '| Example | Baseline projected gas (%) | Open-fill-disabled projected gas (%) |',
        '|---|---:|---:|']
    for run in ('m05q01_001','m06q01_001','m05q20_001','m04q80_001'):
        r=next(r for r in summary if r['run']==run and float(r['threshold'])==.7)
        lines.append(f"| {run} | {100*float(r['baseline_projected_gas_fraction']):.3f} | {100*float(r['closed_projected_gas_fraction']):.3f} |")
    lines += ['',
        'These examples use 20 sampled frames, not complete-record averages. Long gas objects, open annular cores and interfaces intersecting the crop can fail the closed-outline requirement. The large decreases therefore expose lost coverage as well as artifact suppression. A dispersed-object method may use this restriction, but it must not turn unmarked portions of plug/churn/annular images into certified liquid.',
        '',
        '## Exploratory dispersed bubble-candidate trajectories','',
        'On the full 501-frame dispersed sequence, connected candidates must have area 15–1,000 pixels, width 3–60, height 3–35 and width/height ≤4, inside rows 20:69 and columns 10:1006. Adjacent matches use reciprocal nearest centroid distance, displacement ≤30 pixels, vertical displacement ≤5 pixels, area ratio 0.5–2.0, and a ≥2-pixel margin over competitors. A retained chain has at least five matched intervals and consecutive axial/vertical step changes no larger than 2 pixels. These transparent exploratory restrictions deliberately select persistent small objects; they are not validated detector-error bounds.',
        '',
        '| Threshold | Candidate observations | Stable chains | Retained matched intervals | Median px/frame | Median px/s |',
        '|---|---:|---:|---:|---:|---:|']
    for r in bsummary:
        lines.append(f"| {float(r['threshold'])/100:.2f} | {r['candidate_objects']} | {r['stable_chains']} | {r['stable_intervals']} | {float(r['median_px_frame']):.3f} | {float(r['median_px_s']):,.0f} |")
    lines += ['',
        'Multiple objects can be matched in one frame pair, so the interval counts can exceed 500. Chains can be fragments of the same physical bubble and are not proven independent experimental observations. Reported medians weight retained matched intervals, not whole bubbles or phase mass. The three thresholds select different candidates and chains.',
        '',
        '![Selected candidate trajectories](../figures/closed_contour_bubble_tracks.png)',
        '',
        'The figure shows all threshold-0.70 stable chains that cover frames 200–210, with circles on their current centroid positions. These selected candidates visibly correspond to small moving features in the raw images. The consistency with previous broad-profile motion is not independent validation because both methods use the same images.',
        '',
        'This special handling can recover a bubble-translation descriptor where a whole-cross-section gas-dominance gate misses dispersed events. **No liquid velocity is inferred from the background, and no material slip ratio is established.** Dark specks may remain ambiguous; touching, coalescing, cropped, tiny and deformed bubbles are omitted or fragmented. A vetted annotation set would be needed to quantify identification errors.',
        '',
        '## Verification and reproduction','',
        'Mask-subset and packed-array roundtrip assertions pass. Saved per-frame gas fractions reproduce the per-run summaries. A known two-object translation verifies centroid direction and reciprocal matching. Saved link displacements, pixel rates, chain-stability gates and summary statistics are reaggregated by this report script. These checks establish computation consistency, not physical phase accuracy.',
        '',
        'Run scripts/closed_contour_probe.py, scripts/closed_contour_bubble_tracks.py, then scripts/closed_contour_report.py using python3 -B.',
        '',
        'Data: [topology summary](../results/closed_contour_summary.csv), [per-frame comparison](../results/closed_contour_frames.csv), [bubble-candidate summary](../results/closed_contour_bubble_summary.csv), [candidate links](../results/closed_contour_bubble_links.csv), [verification](../results/closed_contour_verification.json). Source images, reference and baseline masks were preserved.']
    (BASE/'findings'/'closed_contour_review.md').write_text('\n'.join(lines)+'\n')
    verification=dict(passed=True,n_numeric_reaggregations=numeric,mask_summary_rows=len(summary),
        per_frame_rows=len(frames),bubble_link_rows=len(links),bubble_summary_rows=len(bsummary),
        displayed_tracks=chosen,scope='Saved computation and selection consistency; no physical phase ground truth.')
    (OUT/'closed_contour_verification.json').write_text(json.dumps(verification,indent=2)+'\n')
    print(json.dumps(verification))


if __name__=='__main__':run()
