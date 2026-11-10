#!/usr/bin/env python3
"""Test passage-conditioned image motion and stricter within-mask patch motion.

The outputs are conditional image descriptors, not measured material slip.
Every classification, rejection and denominator is retained for review.
"""
from __future__ import annotations
import csv
from collections import defaultdict
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np
from PIL import Image
from scipy.ndimage import binary_erosion, gaussian_filter1d

BASE=Path(__file__).resolve().parents[1]
V15=BASE.parent
DATA=V15.parent/'data_processed'
FRAME=V15/'frame_motion_evaluation'
OUT=BASE/'results'
OUT.mkdir(parents=True,exist_ok=True)
LUMEN=(18,71)
PROBE=(496,520)
KEYS=('masks_060','masks_070','masks_080')
LABELS={-1:'unresolved',0:'liquid_candidate',1:'gas_candidate'}
MIN_EPISODE_FRAMES=5
MIN_PHASE_PAIRS=20
MIN_PHASE_EPISODES=3
MIN_LOCAL_FRAMES=10
MIN_LOCAL_TRACKS=20
MIN_DENOM_PX_FRAME=1.
CONFIGS=('threshold_060','threshold_070','threshold_080','consensus','consensus_reference_guard')
GATES=('adjacent','lag2_confirmed')
PATCH_VARIANTS=('majority_070','consensus_95','interior_consensus_100')


def read_csv(path):
    with path.open() as f: return list(csv.DictReader(f))


def save_csv(name,rows):
    assert rows,name
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def mean(x): return float(np.mean(x)) if len(x) else float('nan')
def median(x): return float(np.median(x)) if len(x) else float('nan')
def frac(x): return float(np.mean(x)) if len(x) else float('nan')
def truth(v): return v=='True'


def episodes(states):
    cleaned=states.copy()
    ids=np.full(len(states),-1,dtype=int)
    rows=[]
    start=0
    for stop in range(1,len(states)+1):
        if stop<len(states) and states[stop]==states[start]: continue
        state=int(states[start]); n=stop-start
        accepted=state>=0 and n>=MIN_EPISODE_FRAMES
        if accepted:
            idx=len(rows);ids[start:stop]=idx
            rows.append(dict(episode_id=idx,phase=LABELS[state],first_frame=start,
                             last_frame=stop-1,n_frames=n))
        elif state>=0: cleaned[start:stop]=-1
        start=stop
    return cleaned,ids,rows


def phase_stats(values,byepisode):
    em=[mean(v) for v in byepisode.values() if v]
    return dict(n_accepted_pairs=len(values),n_accepted_episodes=len(em),
        mean_px_frame=mean(values),median_px_frame=median(values),
        q25_px_frame=float(np.quantile(values,.25)) if values else float('nan'),
        q75_px_frame=float(np.quantile(values,.75)) if values else float('nan'),
        episode_equal_mean_px_frame=mean(em),positive_fraction=frac(np.array(values)>0),
        near_zero_fraction=frac(np.abs(values)<MIN_DENOM_PX_FRAME))


def ratio_row(gas,liquid,high_sat,local=False):
    gm,lm=gas['mean_px_frame'],liquid['mean_px_frame']
    finite=math.isfinite(gm) and math.isfinite(lm)
    denominator_safe=finite and abs(lm)>=MIN_DENOM_PX_FRAME
    same_direction=finite and gm*lm>0
    reasons=[]
    if local:
        for label,r in [('gas',gas),('liquid',liquid)]:
            if r['n_accepted_frames']<MIN_LOCAL_FRAMES:reasons.append(label+'_too_few_frames')
            if r['n_accepted_tracks']<MIN_LOCAL_TRACKS:reasons.append(label+'_too_few_tracks')
    else:
        for label,r in [('gas',gas),('liquid',liquid)]:
            if r['n_accepted_pairs']<MIN_PHASE_PAIRS:reasons.append(label+'_too_few_pairs')
            if r['n_accepted_episodes']<MIN_PHASE_EPISODES:reasons.append(label+'_too_few_episodes')
    if high_sat:reasons.append('saturation')
    if not finite:reasons.append('missing_phase_mean')
    elif not denominator_safe:reasons.append('near_zero_liquid_mean')
    if finite and not same_direction:reasons.append('opposed_or_zero_means')
    return dict(gas_mean_px_frame=gm,liquid_mean_px_frame=lm,
        unfiltered_ratio_of_means=gm/lm if finite and lm!=0 else float('nan'),
        descriptive_ratio=gm/lm if not reasons else float('nan'),
        eligible_descriptive_ratio=not reasons,withheld_reasons=';'.join(reasons))


def complement_control():
    # Complementary masks contain the same interfaces and hence the same
    # normalized-correlation displacement even without any liquid tracers.
    spec=importlib.util.spec_from_file_location('ncc',V15/'scripts'/'optical_investigation.py')
    ncc=importlib.util.module_from_spec(spec);spec.loader.exec_module(ncc)
    gas=np.zeros(1016);gas[290:405]=1;gas[440:610]=1;gas[670:731]=1
    next_gas=np.zeros_like(gas);next_gas[:-12]=gas[12:]
    rg=ncc.ncc(gas,next_gas);rl=ncc.ncc(1-gas,1-next_gas)
    assert abs(rg['shift']-rl['shift'])<1e-12 and abs(rg['shift']+12)<.01
    return dict(true_interface_shift=-12.,gas_mask_shift=rg['shift'],
        liquid_complement_shift=rl['shift'],ratio=rg['shift']/rl['shift'],
        interpretation='Ratio 1 arises from shared moving boundaries; no liquid parcel motion was measured.')


def run():
    complement=complement_control()
    baseline={}
    for r in read_csv(V15/'results'/'optical_preferred_candidate_pairs.csv'):
        baseline[(r['run'],int(r['frame_index']))]=r
    temporal={}
    for r in read_csv(FRAME/'results'/'temporal_pairs.csv'):
        if int(r['lag_frames'])==2:temporal[(r['run'],int(r['frame_index']))]=r
    local=defaultdict(list)
    for r in read_csv(FRAME/'results'/'local_patch_tracks.csv'):local[r['run']].append(r)
    inventory=read_csv(BASE/'results'/'native_phase_inventory.csv')
    reference=np.asarray(json.loads((OUT/'native_phase_reference.json').read_text())['radial_reference'])
    frame_rows,pair_rows,episode_rows,phase_rows,ratio_rows=[],[],[],[],[]
    patch_rows,patch_phase_rows,patch_ratio_rows=[],[],[]
    for inv in inventory:
        run=inv['run'];pattern=inv['pattern'];fps=float(inv['fps'])
        high_sat=float(inv['core_saturation_fraction'])>.20
        archive=BASE/'intermediate'/'masks'/(run+'.npz')
        with np.load(archive,allow_pickle=False) as data:
            assert tuple(data['mask_shape'])==(501,80,1016)
            assert np.array_equal(data['frame_ids'],np.arange(501))
            packed={k:data[k] for k in KEYS}
        # Probe boundaries align exactly to byte boundaries: 496/8=62,520/8=65.
        occupancy=np.array([np.unpackbits(packed[k][:,18:71,62:65],axis=-1,
            count=24,bitorder='big').mean(axis=(1,2)) for k in KEYS])
        ref_mae=[];ref_p95=[]
        files=sorted((DATA/run).glob('ImgA*.jpg'))
        for f in files:
            im=np.asarray(Image.open(f).convert('L'),dtype=float)[:,496:520]
            normalized=im/np.maximum(np.percentile(im[30:50],75,axis=0),1)[None,:]
            ratio=gaussian_filter1d(normalized,1,axis=0)/np.maximum(reference[:,None],.01)
            err=np.abs(ratio[18:71]-1)
            ref_mae.append(float(err.mean()));ref_p95.append(float(np.quantile(err,.95)))
        liquid_reference_like=(np.array(ref_mae)<=.12)&(np.array(ref_p95)<=.30)
        states_by_config={}
        for c,config in enumerate(CONFIGS):
            lower=occupancy[c] if c<3 else occupancy.min(axis=0)
            upper=occupancy[c] if c<3 else occupancy.max(axis=0)
            states=np.full(501,-1,dtype=np.int8)
            states[lower>=.60]=1
            states[upper<=.05]=0
            if config=='consensus_reference_guard':states[(states==0)&~liquid_reference_like]=-1
            raw_states=states.copy()
            states,epids,eps=episodes(states)
            states_by_config[config]=(states,epids)
            for i in range(501):
                frame_rows.append(dict(run=run,pattern=pattern,config=config,frame_id=i,fps=fps,
                    occupancy_060=occupancy[0,i],occupancy_070=occupancy[1,i],occupancy_080=occupancy[2,i],
                    raw_state=LABELS[int(raw_states[i])],state=LABELS[int(states[i])],episode_id=int(epids[i]),
                    reference_mean_abs_error=ref_mae[i],reference_p95_abs_error=ref_p95[i],
                    liquid_reference_like=bool(liquid_reference_like[i])))
            for gate in GATES:
                grouped={p:[] for p in (0,1)}
                byepisode={p:defaultdict(list) for p in (0,1)}
                opportunities={p:0 for p in (0,1)}
                for i in range(500):
                    phase=int(states[i]) if states[i]>=0 and states[i]==states[i+1] else -1
                    if phase>=0:opportunities[phase]+=1
                    b=baseline[(run,i)];motion_pass=truth(b['heuristic_pass'])
                    lg=temporal.get((run,i))
                    temporal_pass=bool(lg and truth(lg['temporal_within_tolerance']))
                    phase_through_lag2=bool(i<499 and phase>=0 and states[i+2]==phase and epids[i+2]==epids[i])
                    accepted=bool(phase>=0 and motion_pass and
                        (gate=='adjacent' or (temporal_pass and phase_through_lag2)))
                    v=float(b['displacement_px_per_frame'])
                    if accepted:
                        grouped[phase].append(v);byepisode[phase][int(epids[i])].append(v)
                    pair_rows.append(dict(run=run,pattern=pattern,config=config,motion_gate=gate,
                        frame_index=i,phase=LABELS[phase],episode_id=int(epids[i]) if phase>=0 else -1,
                        motion_pass=motion_pass,lag2_temporal_pass=temporal_pass,
                        phase_through_lag2=phase_through_lag2,accepted=accepted,
                        displacement_px_frame=v,velocity_px_s=v*fps))
                stats={}
                for phase in (0,1):
                    s=phase_stats(grouped[phase],byepisode[phase]);stats[phase]=s
                    all_eps=[e for e in eps if e['phase']==LABELS[phase]]
                    phase_rows.append(dict(run=run,pattern=pattern,config=config,motion_gate=gate,
                        phase=LABELS[phase],fps=fps,n_tagged_frames=int(np.sum(states==phase)),
                        n_tagged_episodes=len(all_eps),n_phase_pair_opportunities=opportunities[phase],
                        **s,mean_px_s=s['mean_px_frame']*fps,
                        acceptance_fraction=s['n_accepted_pairs']/opportunities[phase] if opportunities[phase] else float('nan'),
                        high_saturation=high_sat))
                    for e in all_eps:
                        vs=byepisode[phase].get(e['episode_id'],[])
                        episode_rows.append(dict(run=run,pattern=pattern,config=config,motion_gate=gate,
                            **e,n_accepted_pairs=len(vs),mean_px_frame=mean(vs),mean_px_s=mean(vs)*fps))
                rr=ratio_row(stats[1],stats[0],high_sat)
                eg,el=stats[1]['episode_equal_mean_px_frame'],stats[0]['episode_equal_mean_px_frame']
                ratio_rows.append(dict(run=run,pattern=pattern,config=config,motion_gate=gate,
                    **rr,gas_pairs=stats[1]['n_accepted_pairs'],liquid_pairs=stats[0]['n_accepted_pairs'],
                    gas_episodes=stats[1]['n_accepted_episodes'],liquid_episodes=stats[0]['n_accepted_episodes'],
                    episode_equal_ratio=eg/el if math.isfinite(eg) and math.isfinite(el) and abs(el)>=1 else float('nan'),
                    high_saturation=high_sat))
        # Stronger test: classify the actual source AND matched target patches.
        # Only sampled local-track frames need full masks and erosion.
        tids=sorted({int(r['frame_index']) for r in local[run]}|
                    {int(r['next_frame_index']) for r in local[run]})
        tidmap={t:i for i,t in enumerate(tids)}
        masks=[np.unpackbits(packed[k][tids],axis=-1,count=1016,bitorder='big').astype(bool) for k in KEYS]
        intersection=np.logical_and.reduce(masks)
        union=np.logical_or.reduce(masks)
        lumen=np.zeros((80,1016),dtype=bool);lumen[18:71]=True
        candidate_liquid=(~union)&lumen
        gas_interior=np.stack([binary_erosion(m,structure=np.ones((7,7),bool)) for m in intersection])
        liquid_interior=np.stack([binary_erosion(m,structure=np.ones((7,7),bool)) for m in candidate_liquid])
        variant_masks={'majority_070':(masks[1],(~masks[1])&lumen,.80),
            'consensus_95':(intersection,candidate_liquid,.95),
            'interior_consensus_100':(gas_interior,liquid_interior,1.)}
        run_patch_rows=[]
        for source_index,r in enumerate(local[run]):
            i=int(r['frame_index']);j=int(r['next_frame_index']);x=int(r['window_x']);y=int(r['window_y'])
            w=int(r['width']);h=int(r['height']);dx=float(r['dx_px']);dy=float(r['dy_px'])
            finite=math.isfinite(dx) and math.isfinite(dy)
            bx=x+round(dx) if finite else -10000;by=y+round(dy) if finite else -10000
            src_inside=y>=18 and y+h<=71 and x>=0 and x+w<=1016
            dst_inside=finite and by>=18 and by+h<=71 and bx>=0 and bx+w<=1016
            for variant,(gm,lm,purity) in variant_masks.items():
                def fractions(t,xx,yy,inside):
                    return (float(gm[tidmap[t],yy:yy+h,xx:xx+w].mean()),
                            float(lm[tidmap[t],yy:yy+h,xx:xx+w].mean())) if inside else (float('nan'),float('nan'))
                sg,sl=fractions(i,x,y,src_inside);dg,dl=fractions(j,bx,by,dst_inside)
                source_phase='gas_candidate' if sg>=purity else 'liquid_candidate' if sl>=purity else 'unresolved'
                target_phase='gas_candidate' if dg>=purity else 'liquid_candidate' if dl>=purity else 'unresolved'
                assigned=source_phase if source_phase==target_phase and source_phase!='unresolved' else 'unresolved'
                accepted=truth(r['accepted']) and assigned!='unresolved'
                row=dict(run=run,pattern=pattern,variant=variant,method=r['method'],frame_index=i,
                    next_frame_index=j,window_x=x,window_y=y,width=w,height=h,matched_x=bx,matched_y=by,
                    dx_px=dx,dy_px=dy,velocity_px_s=dx*fps,source_in_lumen=src_inside,target_in_lumen=dst_inside,
                    source_gas_fraction=sg,source_liquid_fraction=sl,target_gas_fraction=dg,target_liquid_fraction=dl,
                    source_phase=source_phase,target_phase=target_phase,assigned_phase=assigned,
                    original_motion_accepted=truth(r['accepted']),accepted=accepted)
                patch_rows.append(row);run_patch_rows.append(row)
        for variant in PATCH_VARIANTS:
            for method in ('raw','background_subtracted'):
                stats={}
                for phase in ('gas_candidate','liquid_candidate'):
                    allrows=[r for r in run_patch_rows if r['variant']==variant and r['method']==method]
                    source=[r for r in allrows if r['source_phase']==phase]
                    good=[r for r in allrows if r['assigned_phase']==phase and r['accepted']]
                    perframe=defaultdict(list)
                    for r in good:perframe[r['frame_index']].append(r['dx_px'])
                    frame_means=[mean(v) for v in perframe.values()]
                    s=dict(n_source_candidate_tracks=len(source),n_phase_consistent_tracks=sum(r['assigned_phase']==phase for r in allrows),
                        n_accepted_tracks=len(good),n_accepted_frames=len(perframe),mean_px_frame=mean(frame_means),
                        median_px_frame=median(frame_means),patch_weighted_mean_px_frame=mean([r['dx_px'] for r in good]),
                        axial_locations=','.join(str(v) for v in sorted({r['window_x'] for r in good})),
                        radial_locations=','.join(str(v) for v in sorted({r['window_y'] for r in good})))
                    stats[phase]=s
                    patch_phase_rows.append(dict(run=run,pattern=pattern,variant=variant,method=method,phase=phase,
                        **s,mean_px_s=s['mean_px_frame']*fps,high_saturation=high_sat))
                rr=ratio_row(stats['gas_candidate'],stats['liquid_candidate'],high_sat,local=True)
                patch_ratio_rows.append(dict(run=run,pattern=pattern,variant=variant,method=method,**rr,
                    gas_tracks=stats['gas_candidate']['n_accepted_tracks'],liquid_tracks=stats['liquid_candidate']['n_accepted_tracks'],
                    gas_frames=stats['gas_candidate']['n_accepted_frames'],liquid_frames=stats['liquid_candidate']['n_accepted_frames']))
        print(run,'passage and within-mask checks complete',flush=True)
    assert len(inventory)==26 and len(pair_rows)==26*5*2*500
    for name,rows in [('passage_frames.csv',frame_rows),('passage_pairs.csv',pair_rows),
        ('passage_episodes.csv',episode_rows),('passage_phase_summary.csv',phase_rows),('passage_ratios.csv',ratio_rows),
        ('phase_patch_tracks.csv',patch_rows),('phase_patch_summary.csv',patch_phase_rows),('phase_patch_ratios.csv',patch_ratio_rows)]:
        save_csv(name,rows)
    metadata=dict(probe_columns=list(PROBE),lumen_rows=list(LUMEN),gas_min_occupancy=.60,
        liquid_max_occupancy=.05,minimum_episode_frames=MIN_EPISODE_FRAMES,configs=list(CONFIGS),motion_gates=list(GATES),
        liquid_reference_guard=dict(mean_abs_ratio_error_max=.12,p95_abs_ratio_error_max=.30),
        passage_ratio_gates=dict(minimum_accepted_pairs_per_phase=MIN_PHASE_PAIRS,
            minimum_accepted_episodes_per_phase=MIN_PHASE_EPISODES,min_abs_liquid_mean_px_frame=MIN_DENOM_PX_FRAME,
            same_nonzero_direction=True,maximum_core_saturation=.20),
        local_variants=dict(majority_070=dict(purity=.80,threshold=.70,erosion=0),
            consensus_95=dict(purity=.95,thresholds=[.6,.7,.8],erosion=0),
            interior_consensus_100=dict(purity=1.,thresholds=[.6,.7,.8],erosion_square_width=7)),
        local_ratio_gates=dict(minimum_tracks_per_phase=MIN_LOCAL_TRACKS,minimum_distinct_frames_per_phase=MIN_LOCAL_FRAMES,
            min_abs_liquid_mean_px_frame=MIN_DENOM_PX_FRAME,same_nonzero_direction=True,maximum_core_saturation=.20),
        complement_control=complement,rows=dict(frames=len(frame_rows),pairs=len(pair_rows),episodes=len(episode_rows),
            passage_phase=len(phase_rows),passage_ratios=len(ratio_rows),patches=len(patch_rows),
            patch_phase=len(patch_phase_rows),patch_ratios=len(patch_ratio_rows)),
        definition='Image-motion means conditional on inferred phase appearance; no material phase-velocity identification.',
        limitations=['Unmarked v14 pixels are not verified liquid.',
            'A 24-pixel station labels a 536-pixel global correlation template only by time; it does not determine template phase purity.',
            'Interior consensus is segmentation robustness, not gas/liquid ground truth or tracer identification.',
            'Means of local frame means weight sampled frames equally; accepted windows can vary by phase and time.',
            'Count thresholds are reporting diagnostics, not effective independent sample sizes or confidence levels.'])
    (OUT/'conditioned_metadata.json').write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n')
    print(json.dumps(metadata['rows']),flush=True)


if __name__=='__main__':run()
