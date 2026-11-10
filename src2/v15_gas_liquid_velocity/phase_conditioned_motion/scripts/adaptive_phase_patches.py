#!/usr/bin/env python3
"""Adaptive-location smaller patches within eroded consensus phase candidates.

Source selection uses texture only before matching; selected liquid regions are
not verified liquid material and correlation is not phase-velocity truth.
"""
from pathlib import Path
from collections import defaultdict
import csv
import hashlib
import importlib.util
import json
import math
import os
import time
import numpy as np
from PIL import Image
from scipy.ndimage import binary_erosion, gaussian_filter, shift

BASE=Path(__file__).resolve().parents[1]
V15=BASE.parent
DATA=V15.parent/'data_processed'
OUT=BASE/'results'
FIG=BASE/'figures'
os.environ.setdefault('MPLCONFIGDIR',str(BASE/'intermediate'/'adaptive_mplcache'))
SOURCE=V15/'frame_motion_evaluation'/'scripts'/'local_patch_motion.py'
spec=importlib.util.spec_from_file_location('local_ncc',SOURCE)
L=importlib.util.module_from_spec(spec);spec.loader.exec_module(L)
L.DX,L.DY=64,4
SIZES=((64,12),(32,12))
SAMPLED=np.linspace(0,499,50,dtype=int)
GRID_X=tuple(range(96,921,32))
GRID_Y=(24,36,48,56)
PHASES=('gas_candidate','liquid_candidate')


def read(name):
    with (OUT/name).open() as f:return list(csv.DictReader(f))


def save(name,rows):
    with (OUT/name).open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)


def mean(values):return float(np.mean(values)) if len(values) else float('nan')


def in_search(x,y,w,h):
    return x-64>=0 and x+w+64<=1016 and y-4>=0 and y+h+4<=80


def match(a,b,x,y,w,h):
    L.WIDTH,L.HEIGHT=w,h
    f=L.match(a,b,x,y);r=None;rt=float('nan')
    bx=x+f['ix'];by=y+f['iy']
    rb=in_search(bx,by,w,h)
    if f['passed'] and rb:
        r=L.match(b,a,bx,by)
        rt=float(np.hypot(f['ix']+r['dx'],f['iy']+r['dy']))
    ok=bool(f['passed'] and r and r['passed'] and rt<=1.)
    return f,r,rt,rb,ok


def synthetic():
    rng=np.random.default_rng(20261024)
    a=gaussian_filter(rng.normal(size=(80,1016)),.8)
    a=a/a.std()*30
    tests=[]
    for w,h in SIZES:
        for dx,dy in [(-13,0),(17,2),(-6.25,-1.2),(0,0)]:
            b=shift(a,(dy,dx),order=3,mode='constant')
            f,r,rt,rb,ok=match(a,b,320,36,w,h)
            error=float(np.hypot(f['dx']-dx,f['dy']-dy))
            assert ok and error<.16,(w,dx,dy,f,error)
            tests.append(dict(width=w,height=h,true_dx=dx,true_dy=dy,
                measured_dx=f['dx'],measured_dy=f['dy'],error_px=error,roundtrip_px=rt))
        f,*_=match(np.ones_like(a),np.ones_like(a),320,36,w,h)
        assert not f['passed']
    return dict(translations=tests,constant_rejected=True,
                maximum_error_px=max(r['error_px'] for r in tests))


def run():
    start=time.monotonic();validation=synthetic()
    tracks=[];opportunities=[];summary=[];ratios=[]
    for inv in read('native_phase_inventory.csv'):
        run=inv['run'];pattern=inv['pattern'];fps=float(inv['fps'])
        high_sat=float(inv['core_saturation_fraction'])>.20
        tids=sorted(set(map(int,SAMPLED))|set(map(int,SAMPLED+1)))
        lookup={t:i for i,t in enumerate(tids)}
        with np.load(BASE/'intermediate'/'masks'/(run+'.npz'),allow_pickle=False) as z:
            assert tuple(z['mask_shape'])==(501,80,1016)
            masks=[np.unpackbits(z[k][tids],axis=-1,count=1016,bitorder='big').astype(bool)
                   for k in ('masks_060','masks_070','masks_080')]
        gas=np.logical_and.reduce(masks)
        lumen=np.zeros((80,1016),bool);lumen[18:71]=True
        liquid=(~np.logical_or.reduce(masks))&lumen
        interiors={p:np.stack([binary_erosion(m,structure=np.ones((7,7),bool)) for m in arr])
                   for p,arr in zip(PHASES,(gas,liquid))}
        files=sorted((DATA/run).glob('ImgA*.jpg'))
        assert len(files)==501
        frames=np.stack([np.asarray(Image.open(f).convert('L')) for f in files])
        background=frames.mean(axis=0)
        runtracks=[]
        for t in SAMPLED:
            t=int(t);a=frames[t].astype(float)-background;b=frames[t+1].astype(float)-background
            for w,h in SIZES:
                grid=[(x,y) for x in GRID_X for y in GRID_Y if in_search(x,y,w,h)]
                for phase in PHASES:
                    source_mask=interiors[phase][lookup[t]]
                    target_mask=interiors[phase][lookup[t+1]]
                    candidates=[];area_count=0;low_texture=0;saturated_source=0
                    for x,y in grid:
                        if not source_mask[y:y+h,x:x+w].all():continue
                        area_count+=1
                        std=float(a[y:y+h,x:x+w].std())
                        if std<1.:
                            low_texture+=1;continue
                        sat=float(np.mean(frames[t,y:y+h,x:x+w]>=250))
                        if sat>.20:
                            saturated_source+=1;continue
                        candidates.append((std,x,y,sat))
                    # Stable deterministic order; no target-frame match or speed
                    # enters source selection. Avoid overlapping source boxes.
                    candidates.sort(key=lambda v:(-v[0],v[1],v[2]))
                    selected=[]
                    for c in candidates:
                        _,x,y,_=c
                        if any(x<sx+w and x+w>sx and y<sy+h and y+h>sy for _,sx,sy,_ in selected):continue
                        selected.append(c)
                        if len(selected)==2:break
                    counts=defaultdict(int)
                    for rank,(std,x,y,source_sat) in enumerate(selected,1):
                        f,r,rt,rb,motion_ok=match(a,b,x,y,w,h)
                        bx=x+f['ix'];by=y+f['iy']
                        inside=bx>=0 and bx+w<=1016 and by>=18 and by+h<=71
                        phase_ok=bool(inside and target_mask[by:by+h,bx:bx+w].all())
                        target_sat=float(np.mean(frames[t+1,by:by+h,bx:bx+w]>=250)) if inside else float('nan')
                        accepted=bool(motion_ok and phase_ok and target_sat<=.20)
                        reason='accepted' if accepted else 'motion_failure' if not motion_ok else 'destination_phase_failure' if not phase_ok else 'target_saturation'
                        counts[reason]+=1
                        row=dict(run=run,pattern=pattern,width=w,height=h,phase=phase,frame_index=t,next_frame_index=t+1,
                            selected_rank=rank,source_x=x,source_y=y,matched_x=bx,matched_y=by,
                            source_std=std,source_saturation=source_sat,target_saturation=target_sat,
                            dx_px=f['dx'],dy_px=f['dy'],velocity_px_s=f['dx']*fps,
                            peak=f['peak'],margin=f['margin'],texture=f['texture'],boundary=f['boundary'],
                            forward_pass=f['passed'],reverse_search_in_bounds=rb,
                            reverse_dx=r['dx'] if r else float('nan'),reverse_dy=r['dy'] if r else float('nan'),
                            reverse_peak=r['peak'] if r else float('nan'),reverse_margin=r['margin'] if r else float('nan'),
                            reverse_pass=r['passed'] if r else False,roundtrip_px=rt,
                            motion_pass=motion_ok,destination_phase_pass=phase_ok,accepted=accepted,reason=reason)
                        tracks.append(row);runtracks.append(row)
                    opportunities.append(dict(run=run,pattern=pattern,width=w,height=h,phase=phase,frame_index=t,
                        grid_positions=len(grid),pure_interior_positions=area_count,
                        low_texture_positions=low_texture,saturated_source_positions=saturated_source,
                        texture_eligible_positions=len(candidates),selected_patches=len(selected),
                        no_interior_area=area_count==0,no_texture_eligible=len(candidates)==0,
                        motion_failures=counts['motion_failure'],destination_phase_failures=counts['destination_phase_failure'],
                        target_saturation_failures=counts['target_saturation'],accepted_patches=counts['accepted']))
        for w,h in SIZES:
            stats={}
            for phase in PHASES:
                attempted=[r for r in runtracks if r['width']==w and r['phase']==phase]
                good=[r for r in attempted if r['accepted']]
                groups=defaultdict(list)
                for r in good:groups[r['frame_index']].append(r['dx_px'])
                v=[mean(x) for x in groups.values()]
                r=dict(run=run,pattern=pattern,width=w,height=h,phase=phase,selected_tracks=len(attempted),
                    accepted_tracks=len(good),accepted_frames=len(groups),mean_px_frame=mean(v),
                    median_px_frame=float(np.median(v)) if v else float('nan'),
                    mean_px_s=mean(v)*fps,patch_weighted_mean_px_frame=mean([r['dx_px'] for r in good]),
                    positive_fraction=float(np.mean(np.array(v)>0)) if v else float('nan'),
                    near_zero_fraction=float(np.mean(np.abs(v)<1)) if v else float('nan'),
                    high_saturation=high_sat)
                stats[phase]=r;summary.append(r)
            g,l=stats['gas_candidate'],stats['liquid_candidate'];gm,lm=g['mean_px_frame'],l['mean_px_frame']
            reasons=[]
            for label,r in [('gas',g),('liquid',l)]:
                if r['accepted_tracks']<20:reasons.append(label+'_too_few_tracks')
                if r['accepted_frames']<10:reasons.append(label+'_too_few_frames')
            finite=math.isfinite(gm) and math.isfinite(lm)
            if not finite:reasons.append('missing_mean')
            elif abs(lm)<1:reasons.append('near_zero_liquid_mean')
            if finite and gm*lm<=0:reasons.append('opposed_or_zero_means')
            if high_sat:reasons.append('saturation')
            ratios.append(dict(run=run,pattern=pattern,width=w,height=h,gas_mean_px_frame=gm,liquid_mean_px_frame=lm,
                gas_tracks=g['accepted_tracks'],liquid_tracks=l['accepted_tracks'],gas_frames=g['accepted_frames'],liquid_frames=l['accepted_frames'],
                unfiltered_ratio=gm/lm if finite and lm!=0 else float('nan'),
                descriptive_ratio=gm/lm if not reasons else float('nan'),eligible_descriptive_ratio=not reasons,
                withheld_reasons=';'.join(reasons)))
        print(run,[(r['width'],r['phase'],r['accepted_tracks']) for r in summary[-4:]],flush=True)
    assert len(opportunities)==26*50*2*2
    save('adaptive_tracks.csv',tracks);save('adaptive_opportunities.csv',opportunities)
    save('adaptive_phase_summary.csv',summary);save('adaptive_ratios.csv',ratios)
    metadata=dict(sizes=SIZES,frame_indices=SAMPLED.tolist(),grid_x=GRID_X,grid_y=GRID_Y,
        source_selection='Up to two nonoverlapping highest-dynamic-texture wholly interior grid patches; no match outcome enters selection.',
        source_masks='AND of three gas masks; liquid candidate is NOT OR within18:71; both eroded with square7x7.',
        thresholds=[.6,.7,.8],erosion_width=7,search=[64,4],min_peak=.6,min_margin=.03,min_texture=1.,
        maximum_roundtrip_px=1.,maximum_saturation=.2,method='temporal-background-subtracted 2D NCC',
        ratio_definition='Ratio of phase means, where each phase mean averages accepted-window means per sampled frame.',
        ratio_gates=dict(minimum_tracks=20,minimum_frames=10,minimum_abs_liquid_mean_px_frame=1.,same_nonzero_direction=True),
        original_ncc_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),synthetic=validation,
        rows=dict(opportunities=len(opportunities),tracks=len(tracks),phase_summary=len(summary),ratios=len(ratios)),
        elapsed_seconds=time.monotonic()-start,
        limitations=['Phase-candidate purity is not independently established phase identity.',
            'Texture-maximizing selection may favor residual bubbles, reflections or refractive features.',
            'Smaller windows can be more ambiguous and less representative.',
            'Image descriptors do not establish material gas/liquid velocity or slip.'])
    (OUT/'adaptive_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    report(summary,ratios,opportunities,tracks,metadata)


def report(summary,ratios,opportunities,tracks,metadata):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    lines=['# Adaptive interior-patch feasibility','',
        '**This experiment tests whether the earlier fixed patch size/placement unnecessarily removed liquid-candidate observations.** It uses smaller source patches selected wholly inside eroded consensus masks before matching. These remain optical phase candidates, not material tracers.','',
        'Each recording supplies50 fixed adjacent pairs. Source grid x=96..920 in steps of32, y=24/36/48/56, filtered so the search window fits. Per phase and size, up to two nonoverlapping highest-texture patches are selected before looking at the target match. Gas requires agreement of all three masks; liquid candidates are outside all three within rows18:71. Both masks are eroded by a7×7 square. The entire source and matched destination patch must stay inside the same eroded phase candidate.','',
        'NCC uses ±64 axial/±4 radial pixels, peak≥0.60, competitor margin≥0.03, texture standard deviation≥1, no boundary hit, reciprocal matching within1px, and≤20% saturated pixels. Reverse searches that do not fit the frame are rejected. Source texture uses original intensity minus the full501-frame temporal background. No velocity enters source selection.','',
        '| Patch size | Gas-candidate accepted tracks | Liquid-candidate accepted tracks | Recordings with eligible two-phase descriptor |',
        '|---|---:|---:|---:|']
    for w,h in SIZES:
        totals=[sum(r['accepted_tracks'] for r in summary if r['width']==w and r['phase']==p) for p in PHASES]
        eligible=[r for r in ratios if r['width']==w and r['eligible_descriptive_ratio']]
        lines.append(f'| {w}×{h} | {totals[0]} | {totals[1]} | {len(eligible)}/26 |')
    lines+=['','A descriptor requires≥20 accepted tracks and≥10 distinct sampled frames for each phase candidate, same nonzero mean direction, |liquid mean|≥1px/frame and no high whole-run saturation. The mean for each phase is the mean of its per-frame accepted-window means. These are information gates, not independent sample sizes or accuracy bounds.','',
        '| Recording | Size | Gas tracks / frames | Liquid tracks / frames | Ratio, if eligible |',
        '|---|---|---:|---:|---:|']
    for r in ratios:
        ratio=f'{r["descriptive_ratio"]:.4f}' if r['eligible_descriptive_ratio'] else '—'
        lines.append(f'| {r["run"]} | {r["width"]}×{r["height"]} | {r["gas_tracks"]}/{r["gas_frames"]} | {r["liquid_tracks"]}/{r["liquid_frames"]} | {ratio} |')
    lines+=['','## Interpretation','',
        'Smaller adaptive patches test geometric availability separately from motion correspondence. Opportunity CSVs retain source-area absence, low texture, saturated source, matching failure, destination-phase failure and acceptance. Comparisons with the earlier160×16 fixed patches change both size and placement; they do not isolate one causal factor.','',
        'A surviving liquid-candidate match can follow an unresolved bubble or optical structure missed by segmentation. A high-purity inferred mask does not certify liquid parcel identity. Texture-based source selection is explicitly selective and can favor such residual objects. No material phase-velocity ratio is established here.','',
        f'Synthetic translated textures at both patch sizes pass, with maximum displacement error {metadata["synthetic"]["maximum_error_px"]:.6f}px. Constant inputs are rejected. These tests verify numerical behavior for constructed translations, not real-data phase accuracy.','',
        'Reproduce: python3 -B src2/v15_gas_liquid_velocity/phase_conditioned_motion/scripts/adaptive_phase_patches.py. Full outputs use the adaptive_ prefix.']
    lines += ['', '## Observed support and visual limitations', '',
        'Adaptive placement recovers more liquid-candidate matches than the earlier fixed160×16 boxes, so the earlier four matches were not a universal limit on possible image correspondence. Nevertheless, no recording satisfies the two-phase reporting criteria at either tested size. Most recovered liquid-candidate matches belong to the dispersed recording, for which no gas-candidate patch qualifies at these sizes.', '',
        '| Patch size | Liquid-candidate phase/frame opportunities | No wholly interior source patch | Accepted liquid tracks in dispersed run | Accepted liquid tracks in all other runs |',
        '|---|---:|---:|---:|---:|']
    for w,h in SIZES:
        oo=[r for r in opportunities if r['width']==w and r['phase']=='liquid_candidate']
        ss=[r for r in summary if r['width']==w and r['phase']=='liquid_candidate']
        dispersed=sum(r['accepted_tracks'] for r in ss if r['run']=='m07q01_001')
        other=sum(r['accepted_tracks'] for r in ss if r['run']!='m07q01_001')
        lines.append(f'| {w}×{h} | {len(oo)} | {sum(r["no_interior_area"] for r in oo)} | {dispersed} | {other} |')
    lines += ['',
        'The masks leave no eligible liquid interior at most tested phase/frame opportunities. Where candidate interiors exist, correspondence rejection also matters. This is a limitation of the tested segmentation/patch definition and observations; it does not prove that liquid is physically absent or stationary.', '',
        'The surviving-example figure was visually inspected. Several nominal liquid-interior boxes visibly include part of a bubble-like contour or a nearby blurred interface: m04q01 frame152, m05q01 frame61, m07q01 frame0 and m07q05 frame203 are examples. Other boxes contain weak unidentified spots or brightness structure. This shows that even consensus plus erosion of these masks does not certify actual phase purity. These are targeted visual observations, not an annotated error-rate measurement.', '',
        '![Adaptive liquid-candidate examples](../figures/adaptive_liquid_examples.png)', '',
        '![Adaptive support](../figures/adaptive_support.png)', '',
        'Saved arithmetic and independently reconstructed mask-interior checks are available in [adaptive_verification.json](../results/adaptive_verification.json), reproduced by scripts/adaptive_verify.py.']
    (BASE/'findings'/'adaptive_patch_review.md').write_text('\n'.join(lines)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(10,4.5))
    for ax,phase in zip(axes,PHASES):
        counts=[sum(r['accepted_tracks'] for r in summary if r['width']==w and r['phase']==phase) for w,h in SIZES]
        bars=ax.bar([0,1],counts,color='#326bc2' if phase=='gas_candidate' else '#298c63')
        ax.bar_label(bars,padding=3);ax.set_xticks([0,1],['64×12','32×12']);ax.set_title(phase.replace('_',' '))
        ax.set_ylabel('Accepted adaptive interior tracks');ax.set_ylim(0,max(counts+[1])*1.2);ax.grid(axis='y',alpha=.2)
    fig.suptitle('Smaller patches selected entirely inside eroded consensus candidates')
    fig.tight_layout(rect=(0,0,1,.92));fig.savefig(FIG/'adaptive_support.png',dpi=180);plt.close(fig)
    examples=[]
    for run in sorted({r['run'] for r in tracks}):
        good=[r for r in tracks if r['run']==run and r['phase']=='liquid_candidate' and r['accepted'] and r['width']==32]
        if good:examples.append(sorted(good,key=lambda r:(r['frame_index'],r['selected_rank']))[0])
    examples=examples[:8]
    if examples:
        fig,axes=plt.subplots(len(examples),2,figsize=(11,1.7*len(examples)),squeeze=False)
        for r,axs in zip(examples,axes):
            files=sorted((DATA/r['run']).glob('ImgA*.jpg'))
            for ax,t,x,y in [(axs[0],r['frame_index'],r['source_x'],r['source_y']),
                              (axs[1],r['next_frame_index'],r['matched_x'],r['matched_y'])]:
                a=np.asarray(Image.open(files[t]).convert('L'))
                lo=max(0,x-32);hi=min(1016,x+r['width']+32);top=max(0,y-12);bot=min(80,y+r['height']+12)
                ax.imshow(a[top:bot,lo:hi],cmap='gray',vmin=0,vmax=255,interpolation='nearest',extent=[lo,hi,bot,top],aspect='auto')
                ax.add_patch(Rectangle((x,y),r['width'],r['height'],fill=False,ec='#29b580',lw=1.2))
                ax.set_title(f'{r["run"]}, frame{t}, dx={r["dx_px"]:.2f}px',fontsize=8)
        fig.suptitle('Surviving liquid-candidate matches: material identity is unverified',fontsize=11)
        fig.tight_layout(rect=(0,0,1,.95));fig.savefig(FIG/'adaptive_liquid_examples.png',dpi=180);plt.close(fig)
        save('adaptive_visual_examples.csv',examples)
    print('Finished',metadata['rows'],'elapsed',metadata['elapsed_seconds'],flush=True)


if __name__=='__main__':run()
