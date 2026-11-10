"""Audit robust endpoint slopes, resolved motion and common-window ratios."""
from pathlib import Path
import gzip,hashlib,json,re
import numpy as np
import pandas as pd
import independent_object_verify as audit
BASE=audit.BASE;SRC2=audit.SRC2

def pair_slope(t,x):
    i,j=np.triu_indices(len(t),1)
    velocity=float(np.median((x[j]-x[i])/(t[j]-t[i])))
    intercept=np.median(x-velocity*t);error=abs(x-intercept-velocity*t)
    return velocity,float(np.median(error)),float(np.percentile(error,90))

def assess(nodes):
    t=np.array([n['frame'] for n in nodes],float);left=np.array([n['xmin'] for n in nodes],float);right=np.array([n['xmax'] for n in nodes],float)
    vl,el,ql=pair_slope(t,left);vr,er,qr=pair_slope(t,right);v=.5*(vl+vr);tol=max(1.5,.2*abs(v));halves=[]
    for ids in np.array_split(np.arange(len(t)),2):
        halves.append(.5*(pair_slope(t[ids],left[ids])[0]+pair_slope(t[ids],right[ids])[0]) if len(ids)>=3 else np.nan)
    reasons=[]
    if t[-1]-t[0]<20:reasons.append('fewer_than_20_intervals')
    if abs(v)<=2 or abs(v)*(t[-1]-t[0])<50:reasons.append('insufficient_displacement')
    if abs(vl-vr)>tol:reasons.append('front_rear_inconsistent')
    if max(el,er)>3 or max(ql,qr)>12:reasons.append('endpoint_fit_residual')
    if any(not np.isfinite(h) or abs(h-v)>tol for h in halves):reasons.append('half_velocity_inconsistent')
    return dict(robust_px_frame=v,robust_left_px_frame=vl,robust_right_px_frame=vr,
        endpoint_median_error_px=max(el,er),endpoint_q90_error_px=max(ql,qr),half1_px_frame=halves[0],half2_px_frame=halves[1],
        robust_accepted=not reasons,robust_rejection=';'.join(reasons))

def motion_resolved(t):return t['accepted'] and abs(t['signed_px_frame'])>2 and abs(t['signed_px_frame'])*(t['last_frame']-t['first_frame'])>=20

def run():
    inputs=['robust_meniscus_tracks.csv','robust_meniscus_summary.csv','motion_resolved_object_summary.csv','motion_resolved_small_tracks.csv',
        'paired_object_variants.csv','paired_object_intervals.csv','paired_object_consensus.csv','inter_recording_object_ratios.csv','zero_size_extrapolation.csv']
    tables={name:pd.read_csv(BASE/'results'/name) for name in inputs}
    inv=pd.read_csv(SRC2/'v15.1_gas_liquid_velocity/results/scenario_inventory.csv');inv=inv[inv.included]
    derived=[];resolved_rows=[];paired_rows=[];interval_rows=[];robust_missing_correspondence=0;manifest={}
    for r in inv.to_dict('records'):
        name=r['run'];fps=r['fps'];exposure=float(re.search(r'^Shutter time 1\s*:([\d.]+)$',(SRC2/'data_processed'/name/'_config.mii').read_text(),re.M).group(1))
        for th in audit.THRESHOLDS:
            pp=BASE/'intermediate'/f'meniscus_{name}_{th:.1f}.json.gz';manifest[str(pp.relative_to(BASE))]=audit.digest(pp)
            with gzip.open(pp,'rt') as f:menisci=json.load(f)
            robust=tables['robust_meniscus_tracks.csv'];subset=robust[(robust.run==name)&(robust.threshold==th)].set_index('track_id')
            assert len(subset)==len(menisci)
            accepted=[]
            for track in menisci:
                value=assess(track['nodes']);row=subset.loc[track['track_id']]
                for key,v in value.items():audit.check(row[key] if not isinstance(v,str) else (row[key] if pd.notna(row[key]) else ''),v,'robust_'+key,tol=1e-7)
                audit.check(row.robust_px_s,value['robust_px_frame']*fps,'robust_native_fps',tol=1e-6)
                if value['robust_accepted']:
                    accepted.append(dict(track,**value))
                    robust_missing_correspondence+=int(track['lag2_correspondence_fraction']<.8 or not track['reverse_all'])
                audit.bump('independent_robust_tracks')
            intervals={f for t in accepted for f in range(t['first_frame'],t['last_frame'])};blocks={f//100 for f in intervals}
            values=np.array([t['robust_px_frame'] for t in accepted]);rs=tables['robust_meniscus_summary.csv'];rs=rs[(rs.run==name)&(rs.threshold==th)].iloc[0]
            for field,value in dict(accepted_tracks=len(accepted),unique_intervals=len(intervals),observed_blocks=len(blocks),
                median_px_frame=np.median(values) if len(values) else np.nan,median_px_s=np.median(values)*fps if len(values) else np.nan,
                q25_px_s=np.quantile(values,.25)*fps if len(values) else np.nan,q75_px_s=np.quantile(values,.75)*fps if len(values) else np.nan,
                supported=len(accepted)>=3 and len(intervals)>=50 and len(blocks)>=3).items():audit.check(rs[field],value,'robust_summary_'+field,tol=1e-6)
            pp=BASE/'intermediate'/f'tracks_{name}_{th:.1f}.json.gz';manifest[str(pp.relative_to(BASE))]=audit.digest(pp)
            with gzip.open(pp,'rt') as f:objects=json.load(f)
            resolved=[dict(t,accepted=bool(motion_resolved(t))) for t in objects]
            small=[t for t in resolved if t['accepted'] and t['object_class']=='small']
            for group in ('small','isolated_small'):
                values=audit.summary(resolved,group,fps);table=tables['motion_resolved_object_summary.csv'];row=table[(table.run==name)&(table.threshold==th)&(table.group==group)].iloc[0]
                for key,value in values.items():audit.check(row[key],value,'resolved_summary_'+key,tol=1e-6)
                derived.append(dict(run=name,threshold=th,group=group,**values))
            if th==.7:
                table=tables['motion_resolved_small_tracks.csv'];subset=table[table.run==name].set_index('track_id');assert len(subset)==len(small)
                for t in small:
                    row=subset.loc[t['track_id']];nodes=t['nodes'];height=np.median([n['height'] for n in nodes])
                    fields=dict(signed_px_s=t['signed_px_frame']*fps,height_px=height,width_px=np.median([n['width'] for n in nodes]),
                        y_px=np.median([n['y'] for n in nodes]),block=min(int((t['first_frame']+t['last_frame'])/2)//100,4),
                        blur_travel_over_height=abs(t['signed_px_frame'])*fps*exposure/height,isolated_small=t['isolated_small'])
                    for key,value in fields.items():audit.check(row[key],value,'resolved_track_'+key,tol=1e-6)
                    resolved_rows.append({'run':name,**t,**fields})
            ss={};gg={};small_ids=set();large_ids=set()
            for t in small:
                for f in range(t['first_frame'],t['last_frame']):ss.setdefault(f,[]).append(t)
            for t in accepted:
                for f in range(t['first_frame'],t['last_frame']):gg.setdefault(f,[]).append(t)
            common=sorted(set(ss)&set(gg));sv=[];gv=[]
            for f in common:
                smean=float(np.mean([t['signed_px_frame'] for t in ss[f]])*fps);gmean=float(np.mean([t['robust_px_frame'] for t in gg[f]])*fps)
                sv.append(smean);gv.append(gmean);small_ids.update(t['track_id'] for t in ss[f]);large_ids.update(t['track_id'] for t in gg[f])
                interval_rows.append(dict(run=name,threshold=th,interval=f,small_px_s=smean,meniscus_px_s=gmean,n_small=len(ss[f]),n_menisci=len(gg[f])))
            blocks={f//100 for f in common};supported=len(common)>=50 and len(blocks)>=3 and len(small_ids)>=3 and len(large_ids)>=3
            sm=np.mean(sv) if sv else np.nan;gm=np.mean(gv) if gv else np.nan
            fields=dict(common_intervals=len(common),common_blocks=len(blocks),small_tracks=len(small_ids),meniscus_tracks=len(large_ids),
                supported=supported,small_common_px_s=sm,meniscus_common_px_s=gm,ratio=sm/gm if supported and abs(gm)>2*fps else np.nan)
            table=tables['paired_object_variants.csv'];row=table[(table.run==name)&(table.threshold==th)].iloc[0]
            for key,value in fields.items():audit.check(row[key],value,'paired_'+key,tol=1e-6)
            paired_rows.append(dict(run=name,pattern=r['pattern'],fps=fps,threshold=th,**fields))
        print(name,'additional comparisons verified',flush=True)
    intervals=pd.DataFrame(interval_rows).set_index(['run','threshold','interval']).sort_index();saved=tables['paired_object_intervals.csv'].set_index(['run','threshold','interval']).sort_index()
    pd.testing.assert_index_equal(intervals.index,saved.index)
    assert np.allclose(intervals.to_numpy(),saved[intervals.columns].to_numpy(),atol=1e-7,rtol=1e-12)
    pairs=pd.DataFrame(paired_rows);supported_names=[]
    for name,d in pairs.groupby('run'):
        a=d[d.threshold==.7].iloc[0];stable=bool(d.supported.all())
        for col in ('small_common_px_s','meniscus_common_px_s'):
            stable=stable and bool(np.isfinite(d[col]).all() and (np.sign(d[col])==np.sign(a[col])).all()
                and (abs(d[col]-a[col])<=max(1.5*a.fps,.2*abs(a[col]))).all())
        row=tables['paired_object_consensus.csv'];row=row[row.run==name].iloc[0]
        fields=dict(baseline_supported=bool(a.supported),threshold_stable=stable,common_intervals=a.common_intervals,
            small_px_s=a.small_common_px_s if a.supported else np.nan,meniscus_px_s=a.meniscus_common_px_s if a.supported else np.nan,
            ratio=a.ratio if stable else np.nan,ratio_sensitivity_min=d.ratio.min() if stable else np.nan,ratio_sensitivity_max=d.ratio.max() if stable else np.nan,
            material_phase_pair=False)
        for key,value in fields.items():audit.check(row[key],value,'consensus_'+key,tol=1e-6)
        if stable:supported_names.append(name)
    sums=pd.DataFrame(derived);inter=tables['inter_recording_object_ratios.csv'];names=set()
    for (name,group),d in sums.groupby(['run','group']):
        if group!='small':continue
        v=d.loc[d.threshold==.7,'median_px_frame'].iloc[0]
        if d.supported.all() and (np.sign(d.median_px_frame)==np.sign(v)).all() and (abs(d.median_px_frame-v)<=max(1.5,.2*abs(v))).all():names.add(name)
    assert len(inter)==len(names)*(len(names)-1)//2
    for row in inter.to_dict('records'):
        a=sums[(sums.run==row['numerator_run'])&(sums.group=='small')];b=sums[(sums.run==row['denominator_run'])&(sums.group=='small')]
        va=abs(a.loc[a.threshold==.7,'median_px_s'].iloc[0]);vb=abs(b.loc[b.threshold==.7,'median_px_s'].iloc[0])
        audit.check(row['image_speed_ratio'],va/vb,'inter_recording_ratio')
        audit.check(row['segmentation_envelope_min'],abs(a.median_px_s).min()/abs(b.median_px_s).max(),'inter_recording_min')
        audit.check(row['segmentation_envelope_max'],abs(a.median_px_s).max()/abs(b.median_px_s).min(),'inter_recording_max')
    points=tables['motion_resolved_small_tracks.csv'];overlap_rows=[];fits=tables['zero_size_extrapolation.csv']
    for name,d in points.groupby('run'):
        if len(d)<12 or d.block.nunique()<3:continue
        y=d.signed_px_s.to_numpy();h=d.height_px.to_numpy()/20;rad=(d.y_px.to_numpy()-44.5)/26.5
        train_by_block={}
        for block in sorted(d.block.unique()):
            test=d[d.block==block];train=d[d.block!=block]
            test_intervals={i for t in test.to_dict('records') for i in range(t['first_frame'],t['last_frame'])}
            train_intervals={i for t in train.to_dict('records') for i in range(t['first_frame'],t['last_frame'])}
            retain=np.array([row['block']!=block and not test_intervals.intersection(range(row['first_frame'],row['last_frame'])) for row in d.to_dict('records')],bool)
            purged_intervals={i for t in d[retain].to_dict('records') for i in range(t['first_frame'],t['last_frame'])}
            assert not test_intervals.intersection(purged_intervals)
            train_by_block[block]=retain
            overlap_rows.append(dict(run=name,heldout_midpoint_block=int(block),test_intervals=len(test_intervals),
                pre_purge_overlapping_train_intervals=len(test_intervals&train_intervals),
                removed_training_tracks=len(train)-int(retain.sum()),remaining_training_tracks=int(retain.sum()),
                overlapping_train_intervals=len(test_intervals&purged_intervals),
                overlapping_test_fraction=len(test_intervals&purged_intervals)/len(test_intervals) if test_intervals else np.nan))
        for model in ('constant','height','height_squared','height_squared_radial'):
            x=np.ones((len(d),1))
            if model=='height':x=np.column_stack([x,h])
            if 'squared' in model:x=np.column_stack([x,h**2])
            if 'radial' in model:x=np.column_stack([x,rad**2])
            coef=np.linalg.pinv(x)@y;pred=np.full(len(d),np.nan);intercepts=[]
            for block in sorted(d.block.unique()):
                is_test=d.block.to_numpy()==block;is_train=train_by_block[block]
                if is_train.sum()<x.shape[1] or np.linalg.matrix_rank(x[is_train])<x.shape[1]:continue
                fold=np.linalg.pinv(x[is_train])@y[is_train];pred[is_test]=x[is_test]@fold;intercepts.append(fold[0])
            row=fits[(fits.run==name)&(fits.model==model)].iloc[0]
            for key,value in dict(zero_size_intercept_px_s=coef[0],heldout_mae_px_s=np.nanmean(abs(y-pred)),
                leave_block_intercept_min_px_s=min(intercepts) if intercepts else np.nan,
                leave_block_intercept_max_px_s=max(intercepts) if intercepts else np.nan,
                design_condition_number=np.linalg.cond(x),heldout_predicted_tracks=int(np.isfinite(pred).sum())).items():audit.check(row[key],value,'extrapolation_'+key,tol=1e-6)
            audit.bump('extrapolation_models_recalculated')
    pd.DataFrame(overlap_rows).to_csv(BASE/'results/extrapolation_holdout_overlap_audit.csv',index=False)
    for file in inputs:manifest['results/'+file]=audit.digest(BASE/'results'/file)
    for path in [Path(__file__),BASE/'scripts/robust_meniscus.py',BASE/'scripts/object_comparisons.py',BASE/'scripts/paired_object_motion.py']:
        manifest[str(path.relative_to(BASE))]=audit.digest(path)
    result=dict(status='passed_numerical_audit',counts=audit.COUNTS,max_errors=audit.ERRORS,
        accepted_robust_tracks_missing_old_correspondence_gates=robust_missing_correspondence,
        motion_resolved_small_baseline_tracks=len(resolved_rows),common_intervals_checked=len(interval_rows),
        threshold_stable_common_pairs=supported_names,threshold_stable_small_object_runs=sorted(names),inter_recording_pairs=len(inter),
        extrapolation_midpoint_folds_with_temporal_overlap=sum(r['overlapping_train_intervals']>0 for r in overlap_rows),
        extrapolation_midpoint_folds_with_overlap_before_purge=sum(r['pre_purge_overlapping_train_intervals']>0 for r in overlap_rows),
        extrapolation_midpoint_folds=len(overlap_rows),
        physical_phase_pair_verified=False,manifest=manifest,
        limitations=['Paired interval coordinates average fitted whole-trajectory slopes; fits also use nodes outside the common interval set.',
            'The size-extrapolation holdouts assign tracks by midpoint block and purge overlapping training intervals; they remain dependent samples from one short recording.',
            'Robust endpoint method uses new gates rather than requiring original OLS acceptance; current accepted chains nevertheless pass the original lag/reversal correspondence diagnostics.'])
    (BASE/'results/comparison_independent_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('manifest','max_errors')},indent=2))

if __name__=='__main__':run()
