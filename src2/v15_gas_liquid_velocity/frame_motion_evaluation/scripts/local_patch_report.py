#!/usr/bin/env python3
"""Reaggregate and report the local 2-D image-motion experiment."""
from pathlib import Path
import json
import os
import numpy as np
import pandas as pd
BASE = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR', '/tmp/local_patch_mplconfig')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def run():
    out = BASE/'results'
    t = pd.read_csv(out/'local_patch_tracks.csv')
    p = pd.read_csv(out/'local_patch_pairs.csv')
    s = pd.read_csv(out/'local_patch_summary.csv')
    v = json.loads((out/'local_patch_validation.json').read_text())
    assert len(t)==15600 and len(p)==2600 and len(s)==52
    err = float(np.nanmax(abs(t.axial_px_s-t.dx_px*t.fps)))
    assert err<1e-7
    expected = t.forward_pass & t.reverse_pass & (t.roundtrip_px<=1.) & (t.local_saturation<=.2)
    assert np.array_equal(expected, t.accepted)
    for key, q in t.groupby(['run','method','frame_index']):
        saved = p[(p.run==key[0])&(p.method==key[1])&(p.frame_index==key[2])].iloc[0]
        good=q[q.accepted]
        assert len(good)==saved.accepted_patches
        if len(good):
            med=float(good.dx_px.median())
            assert np.isclose(med,saved.median_axial_px_frame,atol=1e-10)
            distributed=len(good)>=3 and good.window_x.nunique()==2 and good.window_y.nunique()>=2
            coherent=distributed and max(abs(good.dx_px-med))<=max(2.,.2*abs(med))
            assert bool(coherent)==saved.coherent
    for key,q in p.groupby(['run','method']):
        saved=s[(s.run==key[0])&(s.method==key[1])].iloc[0]
        assert q.coherent.sum()==saved.coherent_pairs
        assert q.accepted_patches.sum()==saved.accepted_patch_tracks
        d=q[q.coherent].median_axial_px_frame
        if len(d): assert np.isclose(d.median(),saved.median_coherent_axial_px_frame)
        assert (q.coherent.mean()>=.5 and saved.core_saturation<=.2)==saved.supported
    totals=[]
    for method in ['raw','background_subtracted']:
        a=t[t.method==method]; b=p[p.method==method]; c=s[s.method==method]
        both=b[b.coherent & b.global_1d_pass]
        diff=abs(both.median_axial_px_frame-both.global_1d_dx_px)
        totals.append(dict(method=method,attempted_tracks=len(a),accepted_tracks=int(a.accepted.sum()),
                           sampled_pairs=len(b),distributed_pairs=int(b.distributed.sum()),
                           coherent_pairs=int(b.coherent.sum()),supported_runs=int(c.supported.sum()),
                           annular_accepted_tracks=int(a[a.pattern=='wavy annular'].accepted.sum()),
                           annular_attempted_tracks=int(sum(a.pattern=='wavy annular')),
                           annular_coherent_pairs=int(b[b.pattern=='wavy annular'].coherent.sum()),
                           annular_supported_runs=int(c[c.pattern=='wavy annular'].supported.sum()),
                           joint_global_local_pairs=len(both),median_joint_difference_px=float(diff.median()),
                           p90_joint_difference_px=float(diff.quantile(.9))))
    checks=dict(px_s_max_arithmetic_residual=err,all_saved_track_gates_reproduced=True,
                all_pair_gates_reproduced=True,all_run_summaries_reproduced=True,aggregates=totals)
    (out/'local_patch_saved_verification.json').write_text(json.dumps(checks,indent=2))
    rows=s[s.method=='background_subtracted'].reset_index(drop=True)
    fig,axes=plt.subplots(1,2,figsize=(12,10),sharey=True)
    y=np.arange(len(rows))
    for k,method in enumerate(['raw','background_subtracted']):
        a=s[s.method==method].set_index('run').loc[rows.run]
        offset=-.18 if k==0 else .18
        label='Raw images' if k==0 else 'Background subtracted'
        color='#a5a5a5' if k==0 else '#2672ae'
        axes[0].barh(y+offset,100*a.accepted_patch_tracks/a.attempted_patch_tracks,height=.32,label=label,color=color)
        axes[1].barh(y+offset,100*a.coherent_pair_fraction,height=.32,label=label,color=color)
    axes[0].set_yticks(y,[r.run+(' *' if r.core_saturation>.2 else '') for r in rows.itertuples()])
    axes[0].invert_yaxis()
    axes[0].set_xlabel('Accepted local tracks (%)')
    axes[1].set_xlabel('Pairs with distributed agreement (%)')
    for ax in axes:
        ax.set_xlim(0,100);ax.grid(axis='x',alpha=.2);ax.set_axisbelow(True)
    axes[1].axvline(50,color='#a33',linestyle='--',linewidth=1,label='50% reporting gate')
    axes[0].legend(loc='lower right',fontsize=9);axes[1].legend(loc='lower right',fontsize=9)
    fig.suptitle('Local 2-D frame correspondence: 50 pairs x 6 windows per run',fontsize=13)
    fig.text(.02,.012,'* Native core saturation >20%. Agreement gates assess image correspondence, not phase-velocity accuracy.',fontsize=9)
    fig.tight_layout(rect=(0,.035,1,.965))
    (BASE/'figures').mkdir(exist_ok=True)
    for ext in ['png','svg']:fig.savefig(BASE/'figures'/('local_patch_support.'+ext),dpi=180)
    plt.close(fig)
    lines=['# Local two-dimensional frame-motion evaluation','',
           '**Executive finding:** local 2-D tracking recovers individual annular interface/texture displacements, but does not establish a representative annular speed. It adds spatial-consistency evidence; it does not identify separate liquid and gas material velocities.','',
           '## Method and scope','',
           'All 501 native frames in each of the 26 retained recordings supply the temporal mean background. The test uses 50 evenly spaced adjacent pairs (indices 0–499) and six fixed 160×16-pixel windows: x origins 220 and 636, y origins 12, 32 and 52. Coordinates are zero-based. Raw images and temporal-background-subtracted images are compared at identical locations and times. There are 15,600 local tracking attempts across the two methods.','',
           'Full 2-D normalized correlation searches ±96 px axially and ±4 px vertically. Positive displacement means increasing image coordinate. Gates fixed before inspecting results are peak correlation ≥0.60, peak margin ≥0.03 over competitors outside ±5 axial/±2 vertical pixels, texture standard deviation ≥1 gray level, and no search-boundary hit. Reverse tracking starts from the matched integer location; both directions must pass and return within 1 px. Saturated pixels (gray level ≥250) may occupy at most 20% of either matched patch. The vertical search is intentionally narrow for the axial geometry, and can reject larger transverse motions.','',
           '**Distributed agreement:** at least three accepted patches cover both axial positions and at least two radial positions; every accepted axial shift is within max(2 px, 20% of the absolute median) of the median. The illustrative run-reporting gate requires ≥50% of sampled pairs with distributed agreement and whole-run core saturation ≤20%. These are heuristic correspondence checks, not confidence levels or physical validation.','',
           '## Results','',
           '| Method | Accepted tracks | Agreement pairs | Runs passing reporting gate | Annular accepted tracks | Annular agreement pairs |',
           '|---|---:|---:|---:|---:|---:|']
    for a in totals:
        lines.append(f"| {a['method']} | {a['accepted_tracks']}/{a['attempted_tracks']} | {a['coherent_pairs']}/{a['sampled_pairs']} | {a['supported_runs']}/26 | {a['annular_accepted_tracks']}/{a['annular_attempted_tracks']} | {a['annular_coherent_pairs']}/450 |")
    lines+=['',
            'Only m07q05_001 (plug) meets this distributed-agreement reporting gate: 28/50 raw pairs and 29/50 background-subtracted pairs. Its medians over those selected pairs are −101,377 and −101,407 px/s, respectively. The negative sign means movement to smaller axial coordinates.','',
            '**No annular recording meets this gate.** Background subtraction provides 434 accepted local annular tracks out of 2,700 attempts, with ten agreement pairs: m03q10 (3), m04q40 (1), m06q40 (5), m06q80 (1). These sparse local events do not establish a representative annular speed. Saturated m03q40 remains flagged.','',
            'This spatial-agreement criterion is substantially stricter than the previous global-profile criterion. The 1/26 result does not mean only one recording has trackable motion, and does not directly invalidate the previous 17/26 profile-supported recordings. Distinct bubbles, fronts and waves may move differently, while a global profile returns an aggregate shift. Small windows also contain less information and lose features more readily.','',
            '![Local patch support](../figures/local_patch_support.png)','',
            '| Run | Pattern | Background-subtracted accepted tracks / 300 | Agreement pairs / 50 |',
            '|---|---|---:|---:|']
    for r in rows.itertuples():lines.append(f'| {r.run} | {r.pattern} | {r.accepted_patch_tracks} | {r.coherent_pairs} |')
    lines+=['','## Comparison with the previous global profile','',
            'The following comparison uses only the identical sampled pairs passing both local and global gates. It selects easier examples; the global estimate is not independent ground truth.','',
            '| Local method | Jointly passing pairs | Median absolute difference (px/frame) | 90th percentile (px/frame) |',
            '|---|---:|---:|---:|']
    for a in totals:lines.append(f"| {a['method']} | {a['joint_global_local_pairs']} | {a['median_joint_difference_px']:.3f} | {a['p90_joint_difference_px']:.3f} |")
    maxerr=max(r['error_px'] for r in v['translation_tests'])
    lines+=['','## Verification and limitations','',
            f'Five known synthetic 2-D translations, including fractional shifts, pass with maximum Euclidean displacement error {maxerr:.6f} px. Uniform images are rejected. One FFT peak agrees with direct Pearson correlation to {v["peak_direct_pearson_difference"]:.3g}. These limited translation tests verify sign, normalization and arithmetic; they do not recreate multiphase deformation or reflections.','',
            'The saved-result verifier reproduces every local-track acceptance gate, pair gate, run summary and px/s conversion from CSVs. Consecutive frame IDs and agreement between both camera-config FPS values are checked. Camera clock and magnification are not independently calibrated.','',
            'Background subtraction can alter persistent structures. A patch can follow refractive texture, an interface or a wave; these need not follow gas or liquid material. The same moving boundary can occur in both phase masks. Spatial agreement provides no phase mass/volume weighting. Reported spreads describe selected motion events, not measurement uncertainty.','',
            'Recommendation: keep global-profile results as candidate structure-motion summaries subject to temporal/ROI checks, use local patches to reveal spatial inconsistency and illustrate particular events, and leave representative annular speeds unsupported. Do not relabel image regions as measured liquid and gas velocities.','',
            'Reproduce by running scripts/local_patch_motion.py, then scripts/local_patch_report.py, with python3 -B from the repository environment.','',
            'Outputs: [tracks](../results/local_patch_tracks.csv), [pairs](../results/local_patch_pairs.csv), [run summary](../results/local_patch_summary.csv), [metadata and synthetic checks](../results/local_patch_validation.json), [saved arithmetic verification](../results/local_patch_saved_verification.json).','']
    (BASE/'findings/local_patch_review.md').write_text('\n'.join(lines))
    print(json.dumps(checks,indent=2))


if __name__=='__main__':run()
