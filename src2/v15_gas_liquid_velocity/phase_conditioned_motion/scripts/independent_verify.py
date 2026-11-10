#!/usr/bin/env python3
"""Independent reaggregation of phase-conditioned image-motion evidence.

No estimator implementation is imported. Passing checks establish arithmetic
and provenance consistency, not true material phase identity or velocity.
"""
from pathlib import Path
from collections import defaultdict
from itertools import groupby
import csv
import hashlib
import json
import math
import re
import numpy as np
from PIL import Image
from scipy.ndimage import convolve1d, minimum_filter

BASE=Path(__file__).resolve().parents[1]
V15=BASE.parent
SRC2=V15.parent
OUT=BASE/'results'


def read_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def true(value):
    assert value in ('True','False'),value
    return value=='True'


class Checks:
    def __init__(self):
        self.numeric_cells=0
        self.max_difference=0.

    def same(self,actual,expected):
        a,b=float(actual),float(expected)
        self.numeric_cells+=1
        if math.isnan(a) or math.isnan(b):
            assert math.isnan(a) and math.isnan(b),(a,b)
            return
        d=abs(a-b)
        self.max_difference=max(self.max_difference,d)
        assert math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-9),(a,b)


def verify_native():
    checks=Checks()
    metadata=json.loads((OUT/'native_phase_metadata.json').read_text())
    inventory=read_csv(OUT/'native_phase_inventory.csv')
    frames=defaultdict(dict)
    for row in read_csv(OUT/'native_phase_frame_fractions.csv'):
        run,index=row['run'],int(row['frame_id'])
        assert index not in frames[run]
        frames[run][index]=row
    declared={r['run']:r['pattern'] for r in metadata['included_runs']}
    expected={}
    excluded=[]
    valid={'dispersed bubbles','plug','slug','churn','wavy annular'}
    for folder in sorted((SRC2/'data_processed').glob('m??q??_???')):
        markers=list(folder.glob('__*.txt'))
        assert len(markers)==1
        pattern=markers[0].stem[2:].replace('_',' ')
        if folder.name=='m05q05_002' or pattern not in valid:
            excluded.append(folder.name)
        else:
            expected[folder.name]=pattern
    assert declared==expected and len(expected)==26 and len(excluded)==5
    assert sorted(excluded)==sorted(r['run'] for r in metadata['excluded_runs'])
    for r in metadata['source_files']:
        assert hashlib.sha256((SRC2/r['path']).read_bytes()).hexdigest()==r['sha256']
    all_pixels=0
    for row in inventory:
        run=row['run']
        assert row['pattern']==expected[run]
        folder=SRC2/'data_processed'/run
        files=sorted(folder.glob('ImgA*.jpg'))
        ids=[int(re.search(r'([0-9]+)[.]jpg$',p.name).group(1)) for p in files]
        assert ids==list(range(501))
        mii=(folder/'_config.mii').read_text()
        xsv=(folder/'_config.xsv').read_text()
        fps=float(re.search(r'Film speed 1\s*:([0-9.]+)',mii).group(1))
        assert fps==float(re.search(r'^Rate=([0-9.]+)',xsv,re.M).group(1))
        checks.same(row['fps'],fps)
        checks.same(row['duration_s'],500/fps)
        assert set(frames[run])==set(range(501))
        file=BASE/row['packed_file']
        assert hashlib.sha256(file.read_bytes()).hexdigest()==row['packed_sha256']
        assert int(row['packed_bytes'])==file.stat().st_size
        with np.load(file,allow_pickle=False) as npz:
            assert tuple(npz['mask_shape'])==(501,80,1016)
            assert np.array_equal(npz['frame_ids'],np.arange(501))
            assert tuple(npz['lumen_rows'])==(18,71)
            assert int(npz['original_width'])==1016 and str(npz['bitorder'])=='big'
            checks.same(npz['fps'],fps)
            masks=[]
            for key in ('masks_060','masks_070','masks_080'):
                packed=npz[key]
                assert packed.dtype==np.uint8 and packed.shape==(501,80,127)
                mask=np.unpackbits(packed,axis=-1,count=1016,bitorder='big').astype(bool)
                assert not mask[:,:18].any() and not mask[:,71:].any()
                core=mask[:,18:71]
                total=int(core.sum())
                checks.same(row[key+'_gas_pixels'],total)
                checks.same(row[key+'_lumen_fraction'],total/(501*53*1016))
                fractions=core.mean(axis=(1,2))
                for i,v in enumerate(fractions):
                    checks.same(frames[run][i][key+'_lumen_fraction'],v)
                masks.append(mask)
                all_pixels+=mask.size
            union=masks[0]|masks[1]|masks[2]
            intersection=masks[0]&masks[1]&masks[2]
            disagreement=(union&~intersection)[:,18:71]
            checks.same(row['threshold_disagreement_lumen_fraction'],disagreement.mean())
            for i,v in enumerate(disagreement.mean(axis=(1,2))):
                checks.same(frames[run][i]['threshold_disagreement_lumen_fraction'],v)
            nonmonotonic=((masks[0]&~masks[1])|(masks[1]&~masks[2]))[:,18:71]
            checks.same(row['nonmonotonic_lumen_fraction'],nonmonotonic.mean())
    assert len(inventory)==26
    return dict(runs=26,frames=26*501,threshold_masks=26*501*3,
        unpacked_pixels_checked=all_pixels,numeric_cells=checks.numeric_cells,
        max_difference=checks.max_difference,all_mask_inventory_checks_passed=True)


def mean(values):
    return float(np.mean(values)) if len(values) else float('nan')


def expected_ratio(gas,liquid,high_sat,local=False):
    gm,lm=gas['mean_px_frame'],liquid['mean_px_frame']
    finite=math.isfinite(gm) and math.isfinite(lm)
    reasons=[]
    for phase,r in [('gas',gas),('liquid',liquid)]:
        if local:
            if r['n_accepted_frames']<10:reasons.append(phase+'_too_few_frames')
            if r['n_accepted_tracks']<20:reasons.append(phase+'_too_few_tracks')
        else:
            if r['n_accepted_pairs']<20:reasons.append(phase+'_too_few_pairs')
            if r['n_accepted_episodes']<3:reasons.append(phase+'_too_few_episodes')
    if high_sat:reasons.append('saturation')
    if not finite:reasons.append('missing_phase_mean')
    elif abs(lm)<1:reasons.append('near_zero_liquid_mean')
    if finite and gm*lm<=0:reasons.append('opposed_or_zero_means')
    raw=gm/lm if finite and lm!=0 else float('nan')
    return raw,raw if not reasons else float('nan'),reasons


def verify_passage():
    checks=Checks()
    inventory=read_csv(OUT/'native_phase_inventory.csv')
    inv={r['run']:r for r in inventory}
    frame_groups=defaultdict(list)
    frame_rows=read_csv(OUT/'passage_frames.csv')
    for r in frame_rows:frame_groups[(r['run'],r['config'])].append(r)
    baseline={(r['run'],int(r['frame_index'])):r for r in read_csv(V15/'results'/'optical_preferred_candidate_pairs.csv')}
    temporal={(r['run'],int(r['frame_index'])):r for r in read_csv(V15/'frame_motion_evaluation'/'results'/'temporal_pairs.csv') if int(r['lag_frames'])==2}
    ref=np.asarray(json.loads((OUT/'native_phase_reference.json').read_text())['radial_reference'])
    kernel=np.exp(-np.arange(-4,5,dtype=float)**2/2)
    kernel/=kernel.sum()
    labels={-1:'unresolved',0:'liquid_candidate',1:'gas_candidate'}
    configs=('threshold_060','threshold_070','threshold_080','consensus','consensus_reference_guard')
    states_all={}
    episodes_all={}
    reference_images=0
    for run,row in inv.items():
        with np.load(BASE/row['packed_file'],allow_pickle=False) as z:
            # Reconstruct the full-width mask before selecting probe columns;
            # this does not depend on the analysis byte-alignment shortcut.
            occupancy=np.array([np.unpackbits(z[key],axis=-1,count=1016,bitorder='big')[:,18:71,496:520].mean(axis=(1,2)) for key in ('masks_060','masks_070','masks_080')])
        ref_mean=[];ref_p95=[]
        for path in sorted((SRC2/'data_processed'/run).glob('ImgA*.jpg')):
            image=np.asarray(Image.open(path).convert('L'),float)[:,496:520]
            den=np.maximum(np.quantile(image[30:50],.75,axis=0),1)
            norm=image/den
            smoothed=convolve1d(norm,kernel,axis=0,mode='reflect')
            errors=np.abs(smoothed[18:71]/np.maximum(ref[18:71,None],.01)-1)
            ref_mean.append(float(errors.mean()));ref_p95.append(float(np.quantile(errors,.95)))
            reference_images+=1
        guard=(np.asarray(ref_mean)<=.12)&(np.asarray(ref_p95)<=.30)
        for c,config in enumerate(configs):
            saved=sorted(frame_groups[(run,config)],key=lambda r:int(r['frame_id']))
            assert [int(r['frame_id']) for r in saved]==list(range(501))
            lo=occupancy[c] if c<3 else occupancy.min(axis=0)
            hi=occupancy[c] if c<3 else occupancy.max(axis=0)
            raw=np.where(lo>=.60,1,np.where(hi<=.05,0,-1))
            if config=='consensus_reference_guard':raw[(raw==0)&~guard]=-1
            cleaned=raw.copy();epids=np.full(501,-1,int);eps=[]
            for state,indices in groupby(range(501),key=lambda i:int(raw[i])):
                ids=list(indices)
                if state>=0 and len(ids)>=5:
                    eid=len(eps);epids[ids]=eid
                    eps.append(dict(episode_id=eid,phase=labels[state],first_frame=ids[0],last_frame=ids[-1],n_frames=len(ids)))
                elif state>=0:cleaned[ids]=-1
            states_all[(run,config)]=(cleaned,epids)
            episodes_all[(run,config)]=eps
            for i,r in enumerate(saved):
                for j,key in enumerate(('occupancy_060','occupancy_070','occupancy_080')):checks.same(r[key],occupancy[j,i])
                checks.same(r['reference_mean_abs_error'],ref_mean[i]);checks.same(r['reference_p95_abs_error'],ref_p95[i])
                assert true(r['liquid_reference_like'])==bool(guard[i])
                assert r['raw_state']==labels[int(raw[i])]
                assert r['state']==labels[int(cleaned[i])]
                assert int(r['episode_id'])==epids[i]
                checks.same(r['fps'],row['fps'])
    pairs=read_csv(OUT/'passage_pairs.csv')
    accepted=defaultdict(list);accepted_episode=defaultdict(lambda:defaultdict(list))
    opportunities=defaultdict(int);pair_counts=defaultdict(int)
    for r in pairs:
        run,config,gate=r['run'],r['config'],r['motion_gate']
        i=int(r['frame_index']);states,epids=states_all[(run,config)]
        phase=int(states[i]) if states[i]>=0 and states[i]==states[i+1] else -1
        b=baseline[(run,i)];lag=temporal.get((run,i))
        motion=true(b['heuristic_pass'])
        lag_pass=bool(lag and true(lag['temporal_within_tolerance']))
        phase_through=bool(i<499 and phase>=0 and states[i+2]==phase and epids[i+2]==epids[i])
        ok=bool(phase>=0 and motion and (gate=='adjacent' or (lag_pass and phase_through)))
        assert r['phase']==labels[phase]
        assert int(r['episode_id'])==(epids[i] if phase>=0 else -1)
        assert true(r['motion_pass'])==motion
        assert true(r['lag2_temporal_pass'])==lag_pass
        assert true(r['phase_through_lag2'])==phase_through
        assert true(r['accepted'])==ok
        v=float(b['displacement_px_per_frame'])
        checks.same(r['displacement_px_frame'],v)
        checks.same(r['velocity_px_s'],v*float(inv[run]['fps']))
        key=(run,config,gate,labels[phase])
        pair_counts[(run,config,gate)]+=1
        if phase>=0:opportunities[key]+=1
        if ok:
            accepted[key].append(v)
            accepted_episode[key][int(epids[i])].append(v)
    assert len(pairs)==130000 and all(n==500 for n in pair_counts.values())
    phase_rows=read_csv(OUT/'passage_phase_summary.csv')
    summary={}
    for r in phase_rows:
        run,config,gate,phase=r['run'],r['config'],r['motion_gate'],r['phase']
        key=(run,config,gate,phase);values=accepted[key]
        episodes=accepted_episode[key]
        episode_means=[mean(v) for v in episodes.values()]
        lower,median,upper=np.quantile(values,[.25,.5,.75]) if values else [float('nan')]*3
        data=dict(n_accepted_pairs=len(values),n_accepted_episodes=len(episodes),mean_px_frame=mean(values),
            median_px_frame=median,q25_px_frame=lower,q75_px_frame=upper,episode_equal_mean_px_frame=mean(episode_means),
            positive_fraction=mean(np.asarray(values)>0),near_zero_fraction=mean(np.abs(values)<1))
        for k,v in data.items():checks.same(r[k],v)
        phase_id=0 if phase=='liquid_candidate' else 1
        checks.same(r['n_tagged_frames'],int(np.sum(states_all[(run,config)][0]==phase_id)))
        checks.same(r['n_tagged_episodes'],sum(e['phase']==phase for e in episodes_all[(run,config)]))
        checks.same(r['n_phase_pair_opportunities'],opportunities[key])
        checks.same(r['mean_px_s'],data['mean_px_frame']*float(inv[run]['fps']))
        checks.same(r['acceptance_fraction'],len(values)/opportunities[key] if opportunities[key] else float('nan'))
        assert true(r['high_saturation'])==(float(inv[run]['core_saturation_fraction'])>.20)
        summary[key]=data
    episode_rows=read_csv(OUT/'passage_episodes.csv')
    expected_episode_count=sum(len(e) for e in episodes_all.values())*2
    assert len(episode_rows)==expected_episode_count
    for r in episode_rows:
        run,config,gate,phase=r['run'],r['config'],r['motion_gate'],r['phase']
        e=episodes_all[(run,config)][int(r['episode_id'])]
        for k in ('first_frame','last_frame','n_frames'):checks.same(r[k],e[k])
        assert r['phase']==e['phase']
        values=accepted_episode[(run,config,gate,phase)].get(int(r['episode_id']),[])
        checks.same(r['n_accepted_pairs'],len(values));checks.same(r['mean_px_frame'],mean(values))
        checks.same(r['mean_px_s'],mean(values)*float(inv[run]['fps']))
    ratios=read_csv(OUT/'passage_ratios.csv');eligible=defaultdict(list)
    for r in ratios:
        run,config,gate=r['run'],r['config'],r['motion_gate']
        gas=summary[(run,config,gate,'gas_candidate')];liquid=summary[(run,config,gate,'liquid_candidate')]
        high_sat=float(inv[run]['core_saturation_fraction'])>.20
        raw,ratio,reasons=expected_ratio(gas,liquid,high_sat)
        checks.same(r['gas_mean_px_frame'],gas['mean_px_frame']);checks.same(r['liquid_mean_px_frame'],liquid['mean_px_frame'])
        checks.same(r['unfiltered_ratio_of_means'],raw);checks.same(r['descriptive_ratio'],ratio)
        assert true(r['eligible_descriptive_ratio'])==(not reasons)
        assert r['withheld_reasons']==';'.join(reasons)
        for label,data in [('gas',gas),('liquid',liquid)]:
            checks.same(r[label+'_pairs'],data['n_accepted_pairs']);checks.same(r[label+'_episodes'],data['n_accepted_episodes'])
        eg,el=gas['episode_equal_mean_px_frame'],liquid['episode_equal_mean_px_frame']
        er=eg/el if math.isfinite(eg) and math.isfinite(el) and abs(el)>=1 else float('nan')
        checks.same(r['episode_equal_ratio'],er)
        if not reasons:eligible[config+' / '+gate].append(dict(run=run,ratio=ratio))
    return dict(frame_rows=len(frame_rows),pairs=len(pairs),episodes=len(episode_rows),phase_summaries=len(phase_rows),ratios=len(ratios),
        source_images_used_for_independent_reference_guard=reference_images,numeric_cells=checks.numeric_cells,
        max_difference=checks.max_difference,eligible_ratios=dict(eligible),passed=True)


def verify_local():
    checks=Checks()
    inv={r['run']:r for r in read_csv(OUT/'native_phase_inventory.csv')}
    original={}
    original_by_run=defaultdict(list)
    for r in read_csv(V15/'frame_motion_evaluation'/'results'/'local_patch_tracks.csv'):
        key=(r['run'],r['method'],int(r['frame_index']),int(r['window_x']),int(r['window_y']))
        assert key not in original
        original[key]=r;original_by_run[r['run']].append(r)
    tracks=read_csv(OUT/'phase_patch_tracks.csv')
    byrun=defaultdict(list)
    for r in tracks:byrun[r['run']].append(r)
    groups=defaultdict(list)
    for run,rows in byrun.items():
        tids=sorted({int(r['frame_index']) for r in original_by_run[run]}|{int(r['next_frame_index']) for r in original_by_run[run]})
        tidmap={t:i for i,t in enumerate(tids)}
        with np.load(BASE/inv[run]['packed_file'],allow_pickle=False) as z:
            masks=[np.unpackbits(z[key][tids],axis=-1,count=1016,bitorder='big').astype(bool) for key in ('masks_060','masks_070','masks_080')]
        both=masks[0]&masks[1]&masks[2]
        anygas=masks[0]|masks[1]|masks[2]
        valid=np.zeros((80,1016),bool);valid[18:71]=True
        liquid=(~anygas)&valid
        # Minimum filtering over a 7x7 spatial neighborhood is independent of
        # the analysis binary_erosion call and has the same stated semantics.
        inside_gas=minimum_filter(both.astype(np.uint8),size=(1,7,7),mode='constant',cval=0).astype(bool)
        inside_liquid=minimum_filter(liquid.astype(np.uint8),size=(1,7,7),mode='constant',cval=0).astype(bool)
        variants={'majority_070':(masks[1],(~masks[1])&valid,.80),
            'consensus_95':(both,liquid,.95),'interior_consensus_100':(inside_gas,inside_liquid,1.)}
        fps=float(inv[run]['fps'])
        for r in rows:
            i,j=int(r['frame_index']),int(r['next_frame_index'])
            x,y,w,h=[int(r[k]) for k in ('window_x','window_y','width','height')]
            orig=original[(run,r['method'],i,x,y)]
            dx,dy=float(orig['dx_px']),float(orig['dy_px'])
            finite=math.isfinite(dx) and math.isfinite(dy)
            bx=x+round(dx) if finite else -10000
            by=y+round(dy) if finite else -10000
            checks.same(r['matched_x'],bx);checks.same(r['matched_y'],by)
            checks.same(r['dx_px'],dx);checks.same(r['dy_px'],dy)
            checks.same(r['velocity_px_s'],dx*fps)
            assert j==i+1
            src=x>=0 and x+w<=1016 and y>=18 and y+h<=71
            dst=finite and bx>=0 and bx+w<=1016 and by>=18 and by+h<=71
            assert true(r['source_in_lumen'])==src and true(r['target_in_lumen'])==dst
            gm,lm,purity=variants[r['variant']]
            fractions=[]
            for index,xx,yy,inside in [(i,x,y,src),(j,bx,by,dst)]:
                if inside:
                    fractions.extend([float(gm[tidmap[index],yy:yy+h,xx:xx+w].sum())/(w*h),
                                      float(lm[tidmap[index],yy:yy+h,xx:xx+w].sum())/(w*h)])
                else:fractions.extend([float('nan'),float('nan')])
            for k,v in zip(('source_gas_fraction','source_liquid_fraction','target_gas_fraction','target_liquid_fraction'),fractions):checks.same(r[k],v)
            sg,sl,dg,dl=fractions
            sp='gas_candidate' if sg>=purity else 'liquid_candidate' if sl>=purity else 'unresolved'
            dp='gas_candidate' if dg>=purity else 'liquid_candidate' if dl>=purity else 'unresolved'
            assigned=sp if sp==dp and sp!='unresolved' else 'unresolved'
            assert r['source_phase']==sp and r['target_phase']==dp and r['assigned_phase']==assigned
            assert true(r['original_motion_accepted'])==true(orig['accepted'])
            assert true(r['accepted'])==(true(orig['accepted']) and assigned!='unresolved')
            groups[(run,r['variant'],r['method'])].append(r)
    assert len(tracks)==46800
    phases=read_csv(OUT/'phase_patch_summary.csv')
    summaries={}
    totals=defaultdict(lambda:dict(gas_tracks=0,liquid_tracks=0))
    for r in phases:
        run,variant,method,phase=r['run'],r['variant'],r['method'],r['phase']
        rows=groups[(run,variant,method)]
        source=[q for q in rows if q['source_phase']==phase]
        consistent=[q for q in rows if q['assigned_phase']==phase]
        good=[q for q in consistent if true(q['accepted'])]
        perframe=defaultdict(list)
        for q in good:perframe[int(q['frame_index'])].append(float(q['dx_px']))
        means=[mean(v) for v in perframe.values()]
        data=dict(n_source_candidate_tracks=len(source),n_phase_consistent_tracks=len(consistent),
            n_accepted_tracks=len(good),n_accepted_frames=len(perframe),mean_px_frame=mean(means),
            median_px_frame=float(np.median(means)) if means else float('nan'),
            patch_weighted_mean_px_frame=mean([float(q['dx_px']) for q in good]))
        for k,v in data.items():checks.same(r[k],v)
        checks.same(r['mean_px_s'],data['mean_px_frame']*float(inv[run]['fps']))
        for key,column in [('axial_locations','window_x'),('radial_locations','window_y')]:
            expected=','.join(str(v) for v in sorted({int(q[column]) for q in good}))
            assert r[key]==expected
        assert true(r['high_saturation'])==(float(inv[run]['core_saturation_fraction'])>.20)
        summaries[(run,variant,method,phase)]=data
        totals[variant+' / '+method]['gas_tracks' if phase=='gas_candidate' else 'liquid_tracks']+=len(good)
    ratios=read_csv(OUT/'phase_patch_ratios.csv')
    eligible=defaultdict(list)
    for r in ratios:
        run,variant,method=r['run'],r['variant'],r['method']
        gas=summaries[(run,variant,method,'gas_candidate')]
        liquid=summaries[(run,variant,method,'liquid_candidate')]
        raw,ratio,reasons=expected_ratio(gas,liquid,float(inv[run]['core_saturation_fraction'])>.20,local=True)
        checks.same(r['gas_mean_px_frame'],gas['mean_px_frame']);checks.same(r['liquid_mean_px_frame'],liquid['mean_px_frame'])
        checks.same(r['unfiltered_ratio_of_means'],raw);checks.same(r['descriptive_ratio'],ratio)
        assert true(r['eligible_descriptive_ratio'])==(not reasons)
        assert r['withheld_reasons']==';'.join(reasons)
        for label,data in [('gas',gas),('liquid',liquid)]:
            checks.same(r[label+'_tracks'],data['n_accepted_tracks']);checks.same(r[label+'_frames'],data['n_accepted_frames'])
        if not reasons:eligible[variant+' / '+method].append(dict(run=run,ratio=ratio))
    return dict(tracks=len(tracks),phase_summaries=len(phases),ratios=len(ratios),
        numeric_cells=checks.numeric_cells,max_difference=checks.max_difference,
        eligible_ratios=dict(eligible),accepted_track_totals=dict(totals),passed=True)


def verify_sensitivity():
    checks=Checks()
    inventory=read_csv(OUT/'native_phase_inventory.csv')
    baseline={(r['run'],int(r['frame_index'])):r for r in read_csv(V15/'results'/'optical_preferred_candidate_pairs.csv')}
    temporal={(r['run'],int(r['frame_index'])):r for r in read_csv(V15/'frame_motion_evaluation'/'results'/'temporal_pairs.csv') if int(r['lag_frames'])==2}
    ratios=read_csv(OUT/'passage_sensitivity_ratios.csv')
    saved_ratios={(r['run'],r['variant'],r['motion_gate']):r for r in ratios}
    saved_episodes=read_csv(OUT/'passage_sensitivity_episodes.csv')
    episode_lookup={(r['run'],r['variant'],r['motion_gate'],int(r['episode_id'])):r for r in saved_episodes}
    main_ratios={(r['run'],r['motion_gate']):r for r in read_csv(OUT/'passage_ratios.csv') if r['config']=='consensus'}
    variants=[('center_strict',496,520,.05),('left_strict',208,232,.05),('right_strict',624,648,.05),
              ('center_liquidrich_020',496,520,.20),('center_liquidrich_040',496,520,.40)]
    eligible={v+' / '+g:[] for v,*_ in variants for g in ('adjacent','lag2_confirmed')}
    episode_count=0
    labels={0:'liquid_candidate',1:'gas_candidate'}
    for inv in inventory:
        run=inv['run']
        with np.load(BASE/inv['packed_file'],allow_pickle=False) as z:
            masks=[np.unpackbits(z[k],axis=-1,count=1016,bitorder='big') for k in ('masks_060','masks_070','masks_080')]
        for variant,lo,hi,liqmax in variants:
            occ=np.asarray([m[:,18:71,lo:hi].mean(axis=(1,2)) for m in masks])
            raw=np.where(occ.min(axis=0)>=.60,1,np.where(occ.max(axis=0)<=liqmax,0,-1))
            states=raw.copy();epids=np.full(501,-1,int);eps=[]
            for state,indices in groupby(range(501),key=lambda i:int(raw[i])):
                ids=list(indices)
                if state>=0 and len(ids)>=5:
                    eid=len(eps);epids[ids]=eid
                    eps.append(dict(episode_id=eid,phase=labels[state],first_frame=ids[0],last_frame=ids[-1],n_frames=len(ids)))
                elif state>=0:states[ids]=-1
            for gate in ('adjacent','lag2_confirmed'):
                values={0:[],1:[]};byepisode={0:defaultdict(list),1:defaultdict(list)}
                for i in range(500):
                    phase=int(states[i]) if states[i]>=0 and states[i]==states[i+1] else -1
                    motion=true(baseline[(run,i)]['heuristic_pass'])
                    accepted=phase>=0 and motion
                    if gate=='lag2_confirmed':
                        t=temporal.get((run,i))
                        accepted=bool(accepted and t and true(t['temporal_within_tolerance']) and i<499 and states[i+2]==phase and epids[i+2]==epids[i])
                    if accepted:
                        value=float(baseline[(run,i)]['displacement_px_per_frame'])
                        values[phase].append(value);byepisode[phase][int(epids[i])].append(value)
                stats={p:dict(mean_px_frame=mean(values[p]),n_accepted_pairs=len(values[p]),n_accepted_episodes=len(byepisode[p])) for p in (0,1)}
                unfiltered,ratio,reasons=expected_ratio(stats[1],stats[0],float(inv['core_saturation_fraction'])>.20)
                row=saved_ratios[(run,variant,gate)]
                for k,v in [('probe_start',lo),('probe_stop',hi),('liquid_max_occupancy',liqmax),('gas_min_occupancy',.60),
                            ('gas_mean_px_frame',stats[1]['mean_px_frame']),('liquid_mean_px_frame',stats[0]['mean_px_frame']),
                            ('unfiltered_ratio_of_means',unfiltered),('descriptive_ratio',ratio),
                            ('gas_pairs',stats[1]['n_accepted_pairs']),('liquid_pairs',stats[0]['n_accepted_pairs']),
                            ('gas_episodes',stats[1]['n_accepted_episodes']),('liquid_episodes',stats[0]['n_accepted_episodes'])]:checks.same(row[k],v)
                assert true(row['eligible_descriptive_ratio'])==(not reasons)
                assert row['withheld_reasons']==';'.join(reasons)
                if not reasons:eligible[variant+' / '+gate].append(dict(run=run,ratio=ratio))
                if variant=='center_strict':
                    prior=main_ratios[(run,gate)]
                    assert prior['eligible_descriptive_ratio']==row['eligible_descriptive_ratio']
                    for key in ('gas_mean_px_frame','liquid_mean_px_frame','unfiltered_ratio_of_means','gas_pairs','liquid_pairs','gas_episodes','liquid_episodes'):checks.same(row[key],prior[key])
                for e in eps:
                    erow=episode_lookup[(run,variant,gate,e['episode_id'])]
                    assert erow['phase']==e['phase']
                    for key in ('first_frame','last_frame','n_frames'):checks.same(erow[key],e[key])
                    phase=1 if e['phase']=='gas_candidate' else 0
                    vals=byepisode[phase].get(e['episode_id'],[])
                    checks.same(erow['accepted_pairs'],len(vals));checks.same(erow['mean_px_frame'],mean(vals))
                    episode_count+=1
    assert len(ratios)==260 and len(saved_episodes)==episode_count
    return dict(ratio_rows=len(ratios),episode_rows=len(saved_episodes),numeric_cells=checks.numeric_cells,
        max_difference=checks.max_difference,center_baseline_reproduced=True,eligible_ratios=eligible,passed=True)


def main():
    result=dict(native=verify_native(),passage=verify_passage(),local=verify_local(),sensitivity=verify_sensitivity(),
        conditional_analysis_status='all_saved_decisions_and_aggregates_verified',
        physical_interpretation='Arithmetic consistency does not establish true phase identity or material velocities.')
    (OUT/'independent_verification.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
