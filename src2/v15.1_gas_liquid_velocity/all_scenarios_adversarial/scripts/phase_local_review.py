#!/usr/bin/env python3
"""Reconstruct phase-local decisions and selected NCC peaks without FFT code."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw
from scipy.ndimage import minimum_filter

BASE=Path(__file__).resolve().parents[1]
SRC2=BASE.parent.parent
MASKS=SRC2/'v15_gas_liquid_velocity/phase_conditioned_motion/intermediate/masks'

def native(name):
    files=sorted((SRC2/'data_processed'/name).glob('ImgA*.jpg'))
    images=np.stack([np.asarray(Image.open(p).convert('L'),dtype=float) for p in files])
    with np.load(MASKS/(name+'.npz')) as z:
        masks=[np.unpackbits(z[k],axis=-1,count=1016,bitorder='big').astype(bool) for k in ('masks_060','masks_070','masks_080')]
    g=masks[0]&masks[1]&masks[2]
    l=~(masks[0]|masks[1]|masks[2]);l[:,:18]=False;l[:,71:]=False
    return images,{'gas_candidate':g,'liquid_candidate':l}

def verify():
    d=pd.read_csv(BASE/'results/phase_local_summary.csv')
    b=pd.read_csv(BASE/'results/phase_local_blocks.csv')
    q=pd.read_csv(BASE/'results/phase_local_lags.csv')
    z=pd.read_csv(BASE/'results/phase_local_nulls.csv')
    pairs=pd.read_csv(BASE/'results/phase_local_pairs.csv')
    prev=pd.read_csv(BASE.parent/'results/scenario_inventory.csv')
    assert set(d.run)==set(prev[prev.included].run) and len(d)==104 and len(pairs)==52
    checks=0
    for r in d.itertuples():
        bb=b[(b.run==r.run)&(b.erosion_radius==r.erosion_radius)&(b.phase==r.phase)]
        qq=q[(q.run==r.run)&(q.erosion_radius==r.erosion_radius)&(q.phase==r.phase)].set_index('lag')
        zz=z[(z.run==r.run)&(z.erosion_radius==r.erosion_radius)&(z.phase==r.phase)]
        v=r.signed_px_frame;tol=max(2,.2*abs(v)) if np.isfinite(v) else np.nan
        assert np.isclose(r.signed_px_s,v*r.fps,equal_nan=True)
        count=sum(np.isfinite(t.peak) and t.peak>=.2 and t.std>=1 and not t.boundary and abs(t.velocity_px_frame-v)<=tol for t in bb.itertuples() if t.lag==1)
        assert count==r.agreeing_blocks
        lagok=all(np.isfinite(qq.loc[k,'peak']) and qq.loc[k,'peak']>=.2 and qq.loc[k,'std']>=1 and not qq.loc[k,'boundary'] and abs(qq.loc[k,'shift']/k-v)<=tol for k in (1,2,4))
        assert lagok==r.lag_consistent
        nullmax=zz.peak.max()
        assert np.isclose(nullmax,r.null_max_peak,equal_nan=True)
        reasons=[]
        if r.frames_with_100_phase_pixels<20:reasons.append('insufficient_phase_area')
        if r.raw_phase_std<1:reasons.append('no_raw_phase_texture')
        if not np.isfinite(v):reasons.append('undefined_correlation')
        if not lagok:reasons.append('lag_or_texture_inconsistency')
        if count<4:reasons.append('block_inconsistency')
        if not np.isfinite(r.margin) or r.margin<.02:reasons.append('ambiguous_peak')
        if not np.isfinite(nullmax) or r.peak<nullmax+.05:reasons.append('null_not_separated')
        if r.saturation_fraction>.2:reasons.append('phase_saturation')
        actual=[] if pd.isna(r.reasons) else [s for s in r.reasons.split(';') if s]
        assert actual==reasons and r.optical_supported==(not reasons),(r.run,r.phase,reasons,actual)
        checks+=1
    direct=[]
    for name in ('m05q01_001','m04q80_001','m07q01_001'):
        raw,mm=native(name);dyn=raw-raw.mean(axis=0,keepdims=True)
        for r in d[d.run==name].itertuples():
            mask=minimum_filter(mm[r.phase].astype(np.uint8),size=(1,2*r.erosion_radius+1,2*r.erosion_radius+1),mode='constant',cval=0).astype(bool)
            m=mask[:,18:71]
            assert int((m.sum(axis=(1,2))>=100).sum())==r.frames_with_100_phase_pixels
            assert abs(m.mean()-r.raw_phase_area_fraction)<1e-14
            if not np.isfinite(r.signed_px_frame):continue
            displacement=int(np.rint(r.signed_px_frame))
            xs=slice(max(0,-displacement),min(1016,1016-displacement))
            xt=slice(max(0,displacement),min(1016,1016+displacement))
            n=sa=sb=saa=sbb=sab=0.
            for ids in np.array_split(np.arange(501),5):
                ma=mask[ids[:-1],18:71,xs];mb=mask[ids[1:],18:71,xt]
                chosen=ma&mb
                a=dyn[ids[:-1],18:71,xs][chosen];c=dyn[ids[1:],18:71,xt][chosen]
                n+=len(a);sa+=a.sum();sb+=c.sum();saa+=a@a;sbb+=c@c;sab+=a@c
            corr=(sab-sa*sb/n)/np.sqrt((saa-sa*sa/n)*(sbb-sb*sb/n))
            error=abs(corr-r.peak)
            assert error<2e-5,(name,r.phase,error)
            # FFT products are complex64; aggregated pixel counts therefore
            # have float32 relative rounding, not exact integer arithmetic.
            count_error=abs(n-r.matched_pixel_pairs)
            assert count_error<max(2.,4*np.finfo(np.float32).eps*n),(name,r.phase,n,r.matched_pixel_pairs)
            direct.append(dict(run=name,phase=r.phase,erosion_radius=r.erosion_radius,integer_shift=displacement,
                direct_peak=corr,saved_peak=r.peak,error=error,matched_pixels=int(n),count_error=count_error))
        print('Direct reconstruction',name,flush=True)
    out=dict(status='passed',scope_runs=26,support_decisions_reconstructed=checks,direct_native_peaks=direct,
        max_direct_peak_error=max(r['error'] for r in direct),
        limitation='Arithmetic checks using a separate direct implementation, not experimental validation of phase identity.')
    (BASE/'results/phase_local_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    return d,pairs,out

def report(d,pairs,verification):
    meta=json.loads((BASE/'results/phase_local_metadata.json').read_text())
    challenge=json.loads((BASE/'results/phase_local_challenge.json').read_text())
    ab=[]
    for r in d.itertuples():
        reasons=[] if pd.isna(r.reasons) else [x for x in r.reasons.split(';') if x]
        ab.append(dict(run=r.run,erosion_radius=r.erosion_radius,phase=r.phase,
            original_supported=r.optical_supported,supported_if_fixed_margin_removed=not [x for x in reasons if x!='ambiguous_peak'],
            finite_candidate=np.isfinite(r.signed_px_s)))
    ab=pd.DataFrame(ab);ab.to_csv(BASE/'results/phase_local_gate_ablation.csv',index=False)
    abpair=[]
    for (run,radius),rr in ab.groupby(['run','erosion_radius']):
        abpair.append(dict(run=run,erosion_radius=radius,pair_after_margin_removal=bool(rr.supported_if_fixed_margin_removed.all())))
    abpair=pd.DataFrame(abpair);abpair.to_csv(BASE/'results/phase_local_pair_ablation.csv',index=False)
    sensitivity=[]
    for run in d.run.unique():
        rr=d[d.run==run];vals={}
        for phase in ('gas_candidate','liquid_candidate'):
            r3=rr[(rr.phase==phase)&(rr.erosion_radius==3)].iloc[0]
            r1=rr[(rr.phase==phase)&(rr.erosion_radius==1)].iloc[0]
            difference=abs(r3.signed_px_frame-r1.signed_px_frame)
            vals[phase+'_erosion_consistent']=np.isfinite(difference) and difference<=max(2,.2*abs(r3.signed_px_frame))
        sensitivity.append(dict(run=run,**vals,both_consistent=all(vals.values())))
    sensitivity=pd.DataFrame(sensitivity)
    sensitivity.to_csv(BASE/'results/phase_local_erosion_sensitivity.csv',index=False)
    stable=sensitivity.set_index('run').both_consistent
    both_support=pairs.groupby('run').pair_optical_supported.all()
    consensus=pairs[pairs.erosion_radius==3].copy()
    consensus['both_erosions_optical_supported']=consensus.run.map(both_support)
    consensus['both_erosions_speed_consistent']=consensus.run.map(stable)
    consensus['erosion_consensus_supported']=consensus.both_erosions_optical_supported & consensus.both_erosions_speed_consistent
    consensus.to_csv(BASE/'results/phase_local_erosion_consensus_pairs.csv',index=False)
    lines=['# Distributed phase-interior motion experiment','',
        'This experiment replaces contiguous pure-phase boxes with the union of all available inferred phase-interior pixels. A small bubble or thin irregular region can contribute even when no 32×12-pixel rectangle fits. It used every one of the 501 native frames in all 26 included recordings.','',
        '## Result','',
        '| Spatial erosion | Finite two-coordinate candidates | Two-coordinate optical support | Support after removing only fixed peak-margin gate |',
        '|---|---:|---:|---:|']
    for radius in (3,1):
        pp=pairs[pairs.erosion_radius==radius];aa=abpair[abpair.erosion_radius==radius]
        lines.append(f'| Radius {radius} px ({2*radius+1}×{2*radius+1}) | {int(pp.both_finite.sum())}/26 | {int(pp.pair_optical_supported.sum())}/26 | {int(aa.pair_after_margin_removal.sum())}/26 |')
    lines+=['',f'Requiring support at both erosion radii and compatible velocities leaves {int(consensus.erosion_consensus_supported.sum())}/26 pairs. The consensus table keeps the radius-3 coordinates for every run; it never switches radius to fill a missing coordinate.','',
        'None is independently identified as a gas/liquid material-velocity pair. The new optical matches show that the previous rectangular-patch failure was not a universal limit on image correspondence; they do not establish phase-material identity.','',
        '## Observable and controls','',
        'For every axial trial displacement, only pixels belonging to the same eroded consensus phase candidate at both source and destination enter the correlation. Gas uses the intersection of the 0.60/0.70/0.80 masks. The liquid candidate is the lumen complement of their union. Rows 18:71 are retained. Radius 3 is the initial reference and radius 1 is a declared interface-distance sensitivity. These masks remain optical candidates, not independently verified material labels.','',
        'The estimator subtracts the full-record temporal mean at each pixel. Six summed products (count, source/target sums, squares and cross-product) define one normalized correlation over all matched pixels. Larger visible phase areas receive greater weight; this is neither equal frame weighting nor a volume-weighted phase velocity. No contiguous-patch requirement or tracking from another recording is imposed.','',
        'Five disjoint chronological blocks use lags of 1, 2 and 4 frames. Pairs crossing block boundaries are omitted: 496, 491 and 481 comparisons respectively, while all 501 images contribute. Acceptance requires at least 100 phase pixels in 20 frames, original within-mask gray standard deviation ≥1, dynamic matched-pixel standard deviation ≥1, correlation ≥0.20 at all lags, compatible lag speeds within max(2 px/frame,20%), four of five compatible lag-1 blocks, a peak margin ≥0.02 outside ±5 pixels, and a lag-1 peak at least 0.05 above three blockwise shuffled-order controls. Trial shifts with fewer than max(1000,25% of the maximum overlap count) corresponding pixels are undefined; ±240-pixel search-edge maxima are rejected. More than 20% saturated phase-candidate pixels also rejects the estimate.','',
        'These numerical gates are operational feasibility criteria, not experimentally calibrated error rates. The full-record temporal background is shared; block consistency is not an independent out-of-sample experiment. Three null replicates are a diagnostic and do not establish a p-value. Pure axial matching omits transverse/out-of-plane motion.','',
        '## Adversarial checks','',
        f'Two independent synthetic textures were translated at (−17.3,−7.4) and (+13.25,−11.6) px/frame. The four estimates had maximum absolute error {meta["synthetic"]["maximum_translation_error"]:.6f} px/frame. The opposing-direction case confirms that the method does not force common phase direction or a prescribed velocity ratio. A stationary sequence and untextured liquid interior were rejected.','',
        'A pure moving binary interface with uniform phase interiors was also rejected for both phases. Temporal-background subtraction alone produced approximately −0.006 px/frame in both interiors, illustrating how background structure can produce numbers without material tracers. The raw-phase-variance guard prevents accepting those synthetic numbers.','',
        'A harder synthetic test moved the phase boundary at −12 px/frame while independent interior grayscale textures moved at −21.3 and −5.7 px/frame. Both recovered texture speeds passed, with errors below 0.015 px/frame. Adding static radial shading to otherwise untextured moving phases defeated the raw-variance guard alone, but the shuffled-order control still rejected both candidates. Independently randomizing grayscale texture each frame while retaining the moving phase mask also rejected both. These controls demonstrate that moving masks alone need not determine the reported grayscale speed; they do not certify native liquid texture as material tracers.','',
        f'A separate direct spatial calculation reconstructed {len(verification["direct_native_peaks"])} native correlation peaks from the original images and masks without importing the FFT estimator. Maximum peak difference was {verification["max_direct_peak_error"]:.3g}; all 104 saved support decisions were reconstructed. Synthetic accuracy does not measure native-frame physical accuracy.','',
        'The fixed ±5-pixel competitor rule can reject a broad, repeatable peak. The gate-ablation table therefore removes only that rule while preserving every other criterion. This diagnoses a potentially overconservative rejection; it is not evidence that the added candidates are physically correct. Even this relaxation does not complete all 26 pairs.','',
        '## Per-recording result','',
        '| Recording | Pattern | Radius 3 pair | Radius 1 pair | Radius 1 gas candidate (px/s) | Radius 1 liquid candidate (px/s) |',
        '|---|---|---|---|---:|---:|']
    for run in pairs.run.unique():
        pp=pairs[pairs.run==run].set_index('erosion_radius');r=pp.loc[1]
        fmt=lambda v: f'{v:.1f}' if np.isfinite(v) else 'undefined'
        lines.append(f'| {run} | {r["pattern"]} | {"supported" if pp.loc[3,"pair_optical_supported"] else "withheld"} | {"supported" if r.pair_optical_supported else "withheld"} | {fmt(r.gas_candidate_px_s)} | {fmt(r.liquid_candidate_px_s)} |')
    lines+=['','All finite values are retained in this diagnostic table, including rejected candidates. They must not all be plotted as accepted phase velocities.','',
        '## Native visual challenge','',
        '![Original frames and candidate-mask overlays](../figures/phase_local_masks_examples.png)','',
        'The fixed frame 250 examples show why optical support still falls short of independent phase identification. In m04q10 the liquid-candidate support can lie in narrow strips next to the moving silhouette; refracted interface structure can remain after finite erosion. In m05q05, the liquid-looking area contains small visible objects that may carry texture without tracing the continuous liquid itself. In m04q40 the liquid-candidate support is almost absent. Blue is the gas candidate; orange is the liquid candidate; erosion radius is 1 pixel. These are illustrative frame inspections, not ground-truth annotations.','',
        '## Reproduction','',
        'Run phase_local_experiment.py, phase_local_challenge.py, then phase_local_review.py from this investigation’s scripts directory. Inputs and hashes are in results/phase_local_inputs.csv. Complete estimates and reasons are in phase_local_summary.csv and phase_local_pairs.csv. Block, lag, null and verification outputs preserve all tested alternatives. No earlier v15/v15.1 output was overwritten.']
    (BASE/'findings/phase_local_review.md').write_text('\n'.join(lines)+'\n')

def figure():
    names=['m04q10_001','m05q05_001','m04q40_001']
    canvas=Image.new('RGB',(1016,180*len(names)),(255,255,255));draw=ImageDraw.Draw(canvas)
    for i,name in enumerate(names):
        frame=np.array(Image.open(SRC2/'data_processed'/name/'ImgA000250.jpg').convert('L'))
        with np.load(MASKS/(name+'.npz')) as z:
            mm=[np.unpackbits(z[k][250],axis=-1,count=1016,bitorder='big').astype(bool) for k in ('masks_060','masks_070','masks_080')]
        g=np.logical_and.reduce(mm);l=~np.logical_or.reduce(mm);l[:18]=False;l[71:]=False
        g=minimum_filter(g.astype(np.uint8),size=3,mode='constant').astype(bool)
        l=minimum_filter(l.astype(np.uint8),size=3,mode='constant').astype(bool)
        rgb=np.repeat(frame[:,:,None],3,axis=2).astype(float)
        rgb[g]=.45*rgb[g]+.55*np.array([0,130,255]);rgb[l]=.45*rgb[l]+.55*np.array([255,120,0])
        draw.text((8,i*180+2),name+' frame250, erosion radius1; blue: gas candidate, orange: liquid candidate',fill=(0,0,0))
        canvas.paste(Image.fromarray(frame).convert('RGB'),(0,i*180+18));canvas.paste(Image.fromarray(rgb.astype(np.uint8)),(0,i*180+98))
    canvas.save(BASE/'figures/phase_local_masks_examples.png')

if __name__=='__main__':
    d,pairs,v=verify();report(d,pairs,v);figure();print(json.dumps(v,indent=2))
