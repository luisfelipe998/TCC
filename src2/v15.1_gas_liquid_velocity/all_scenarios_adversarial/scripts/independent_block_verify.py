#!/usr/bin/env python3
"""Independent arithmetic, source and sampled-Pearson audit; imports no generator."""
from pathlib import Path
import hashlib
import json
import re
import time
import numpy as np
import pandas as pd
from PIL import Image
from scipy.signal import correlate
from scipy.ndimage import gaussian_filter1d

BASE=Path(__file__).resolve().parents[1]
OLD=BASE.parent
SEED=20261024
REGIONS=('core','upper_outer','lower_outer')
ROIS={'core':(35,54),'upper_outer':(18,29),'lower_outer':(60,71)}

def audit():
    start=time.monotonic()
    detail=pd.read_csv(BASE/'results/block_ncc_details.csv')
    perrun=pd.read_csv(BASE/'results/block_ncc_per_run.csv')
    selected=pd.read_csv(BASE/'results/selected_wall_per_run.csv')
    split_meta=json.loads((BASE/'results/split_background_metadata.json').read_text())
    split_script_hash=hashlib.sha256((BASE/'scripts/split_background_ncc.py').read_bytes()).hexdigest()
    assert split_meta['script_sha256']==split_script_hash,'Training-background generator run has not finished for current source'
    split_detail=pd.read_csv(BASE/'results/split_background_details.csv')
    split_pairs=pd.read_csv(BASE/'results/split_background_pairs.csv')
    assert len(split_detail)==4680 and len(split_pairs)==156
    assert not split_detail.duplicated(['run','variant','selection_direction','region','block']).any()
    scope=pd.read_csv(OLD/'results/scenario_inventory.csv')
    scope=scope[scope.included].reset_index(drop=True)
    assert len(scope)==26 and set(detail.run)==set(scope.run)
    assert not detail.duplicated(['run','variant','region','block']).any()
    assert len(detail)==5070 and len(perrun)==156 and len(selected)==208
    numeric_fields=0;max_summary_error=0.
    def close(actual,expected,label,tol=1e-7):
        nonlocal numeric_fields,max_summary_error
        numeric_fields+=1
        if pd.isna(expected):
            assert pd.isna(actual),(label,actual,expected);return
        assert np.isfinite(actual) and np.isfinite(expected),(label,actual,expected)
        difference=float(abs(actual-expected));max_summary_error=max(max_summary_error,difference)
        assert difference<=tol+1e-12*abs(expected),(label,actual,expected,difference)
    sources=[];archive=np.load(OLD/'intermediate/annular_profiles.npz')
    for z in scope.itertuples():
        folder=OLD.parent/'data_processed'/z.run
        xsv=(folder/'_config.xsv').read_text();mii=(folder/'_config.mii').read_text()
        fps=float(re.search(r'^Rate=([0-9.]+)',xsv,re.M)[1])
        alternate=float(re.search(r'Film speed 1[ ]*:([0-9.]+)',mii)[1])
        assert fps==alternate==z.fps
        assert np.all(detail.loc[detail.run==z.run,'fps']==fps)
        assert np.all(split_detail.loc[split_detail.run==z.run,'fps']==fps)
        files=sorted(folder.glob('ImgA*.jpg'))
        assert len(files)==501
        ids=[int(re.search(r'([0-9]+)[.]jpg$',f.name)[1]) for f in files]
        assert np.array_equal(np.diff(ids),np.ones(500))
        maximum=0.
        for frame_id in (0,250,500):
            native=np.asarray(Image.open(files[frame_id]).convert('L'))
            assert native.shape==(80,1016)
            for region,(lo,hi) in ROIS.items():
                difference=np.max(abs(native[lo:hi].mean(axis=0)-archive[f'{z.run}__{region}'][frame_id]))
                maximum=max(maximum,float(difference))
                assert difference<2e-5
        sources.append(dict(run=z.run,fps=fps,frames=501,profile_sample_max_error=maximum))
    # Reconstruct every local decision from recorded measurements, with no
    # imported gate function or aggregate result.
    for z in pd.concat([detail,split_detail],ignore_index=True).itertuples():
        tol=max(2.,.2*abs(z.velocity_px_frame))
        close(z.tolerance_px_frame,tol,'tolerance')
        reasons=[]
        if z.texture_coverage<.6:reasons.append('texture')
        if not np.isfinite(z.peak) or z.peak<.3:reasons.append('weak_peak')
        if z.boundary:reasons.append('search_boundary')
        if not np.isfinite(z.null_max_peak) or z.peak<z.null_max_peak+.05:reasons.append('shuffle')
        for name in ('lag2','reverse','x160','x320'):
            value=getattr(z,name+'_velocity_px_frame')
            expected=-z.velocity_px_frame if name=='reverse' else z.velocity_px_frame
            ok=(np.isfinite(value) and abs(value-expected)<=tol and
                not getattr(z,name+'_boundary') and getattr(z,name+'_peak')>=.15 and
                getattr(z,name+'_texture_coverage')>=.6)
            assert bool(ok)==getattr(z,name+'_ok'),(z.run,z.region,z.variant,z.block,name)
            if not ok:reasons.append(name)
        assert bool(not reasons)==z.supported
        assert ';'.join(reasons)==('' if pd.isna(z.reasons) else z.reasons)
        assert z.last_frame-z.first_frame==z.intervals
    reconstructed=[]
    for (run,variant),g in detail.groupby(['run','variant'],sort=False):
        groups={r:g[g.region==r].sort_values('block') for r in REGIONS}
        core=groups['core'];nb=len(core)
        expected_ids=np.array_split(np.arange(501),nb)
        for region,q in groups.items():
            assert np.array_equal(q.block,np.arange(nb))
            assert np.array_equal(q.first_frame,[ids[0] for ids in expected_ids])
            assert np.array_equal(q.last_frame,[ids[-1] for ids in expected_ids])
            assert np.array_equal(q.intervals,core.intervals)
        valid=np.array([all(bool(groups[r].iloc[k].supported) for r in REGIONS) for k in range(nb)])
        w=core.intervals.to_numpy(float)
        mids=(core.first_frame.to_numpy()+core.last_frame.to_numpy())*.5
        fifths=len(set(np.minimum(np.floor(mids[valid]/501*5).astype(int),4)))
        fps=float(core.fps.iloc[0])
        values=np.stack([groups[r].velocity_px_frame.to_numpy() for r in REGIONS])
        whole=np.sum(values*w[None,:],axis=1)/w.sum()*fps
        mean=np.sum(values[:,valid]*w[valid][None,:],axis=1)/w[valid].sum()*fps if valid.any() else np.full(3,np.nan)
        wanted=dict(blocks=nb,joint_accepted_blocks=int(valid.sum()),joint_block_fraction=valid.mean(),
            temporal_fifths=fifths,joint_accepted_intervals=int(w[valid].sum()),total_internal_intervals=int(w.sum()),
            joint_interval_fraction_of_500=w[valid].sum()/500,core_signed_px_s=mean[0],
            wall_signed_px_s=(mean[1]+mean[2])/2,upper_signed_px_s=mean[1],lower_signed_px_s=mean[2],
            all_blocks_core_candidate_px_s=whole[0],all_blocks_wall_candidate_px_s=(whole[1]+whole[2])/2)
        reported=perrun[(perrun.run==run)&(perrun.variant==variant)]
        assert len(reported)==1
        reported=reported.iloc[0]
        for key,expected in wanted.items():close(reported[key],expected,'aggregate:'+key)
        support=valid.sum()/nb>=.6 and fifths>=4 and np.isfinite(mean).all()
        assert bool(support)==reported.supported and not reported.material_phase_pair
        reconstructed.append(dict(run=run,variant=variant,supported=bool(support),blocks=int(valid.sum())))
    # Reconstruct training ranking and all assessment means using the full
    # detailed input table; no generator helper is imported.
    selected_reconstruction=[]
    for z in selected.itertuples():
        q=detail[(detail.run==z.run)&(detail.variant==z.variant)]
        groups={r:q[q.region==r].set_index('block').sort_index() for r in REGIONS}
        core=groups['core']
        train=list(range(5)) if z.selection_direction=='forward' else list(range(5,10))
        test=list(range(5,10)) if z.selection_direction=='forward' else list(range(5))
        def rank(region):
            values=groups[region].loc[train]
            available=values.peak.to_numpy()[np.isfinite(values.peak)]
            peak_mean=float(available.mean()) if len(available) else -np.inf
            return int(values.supported.sum()),peak_mean
        upper_rank=rank('upper_outer');lower_rank=rank('lower_outer')
        choice='upper_outer' if upper_rank>=lower_rank else 'lower_outer'
        assert choice==z.selected_wall
        wall=groups[choice]
        train_ids=[k for k in train if core.loc[k,'supported'] and wall.loc[k,'supported']]
        test_ids=[k for k in test if core.loc[k,'supported'] and wall.loc[k,'supported']]
        assert ','.join(map(str,test_ids))==('' if pd.isna(z.accepted_block_ids) else str(z.accepted_block_ids))
        assert len(train_ids)==z.training_joint_blocks and len(test_ids)==z.test_joint_blocks
        assert bool(len(train_ids)>=3 and len(test_ids)>=3)==z.supported
        assert z.test_blocks==5 and not z.material_phase_pair
        for ids in (train,test,test_ids):
            assert np.array_equal(core.loc[ids,'first_frame'],wall.loc[ids,'first_frame'])
            assert np.array_equal(core.loc[ids,'last_frame'],wall.loc[ids,'last_frame'])
        def weighted(region,indices):
            if not indices:return np.nan
            weights=core.loc[indices,'intervals'].to_numpy(float)
            speeds=groups[region].loc[indices,'velocity_px_frame'].to_numpy(float)
            return float(np.dot(speeds,weights)/weights.sum()*z.fps)
        other='lower_outer' if choice=='upper_outer' else 'upper_outer'
        expect=dict(test_intervals_accepted=int(core.loc[test_ids,'intervals'].sum()),
            test_intervals_total=int(core.loc[test,'intervals'].sum()),
            core_signed_px_s=weighted('core',test_ids),wall_signed_px_s=weighted(choice,test_ids),
            all_test_core_candidate_px_s=weighted('core',test),all_test_wall_candidate_px_s=weighted(choice,test),
            other_wall_same_epochs_candidate_px_s=weighted(other,test_ids))
        for key,value in expect.items():close(getattr(z,key),value,'selection:'+key)
        selected_reconstruction.append(dict(run=z.run,variant=z.variant,direction=z.selection_direction,
            selected_wall=choice,training_joint=len(train_ids),assessment_joint=len(test_ids),supported=bool(z.supported)))
    # Final corrected study: background estimated only from the selection half.
    # Reconstruct every final selection from the newly computed detail table.
    for z in split_pairs.itertuples():
        q=split_detail[(split_detail.run==z.run)&(split_detail.variant==z.variant)&
                       (split_detail.selection_direction==z.selection_direction)]
        groups={r:q[q.region==r].set_index('block').sort_index() for r in REGIONS}
        core=groups['core']
        train=list(range(5)) if z.selection_direction=='forward' else list(range(5,10))
        test=list(range(5,10)) if z.selection_direction=='forward' else list(range(5))
        expected=np.array_split(np.arange(501),10)
        for region,g in groups.items():
            assert np.array_equal(g.index,np.arange(10))
            assert np.array_equal(g.first_frame,[x[0] for x in expected])
            assert np.array_equal(g.last_frame,[x[-1] for x in expected])
        assert z.background_first_frame==int(core.loc[train,'first_frame'].min())
        assert z.background_last_frame==int(core.loc[train,'last_frame'].max())
        background_ids=set(range(z.background_first_frame,z.background_last_frame+1))
        test_frames=set(np.concatenate([expected[k] for k in test]).tolist())
        assert not background_ids.intersection(test_frames)
        def rank(region):
            arr=groups[region].loc[train,'peak'].to_numpy()
            finite=arr[np.isfinite(arr)]
            return int(groups[region].loc[train,'supported'].sum()),float(finite.mean()) if len(finite) else -np.inf
        choice='upper_outer' if rank('upper_outer')>=rank('lower_outer') else 'lower_outer'
        assert choice==z.selected_wall
        wall=groups[choice]
        training=[k for k in train if core.loc[k,'supported'] and wall.loc[k,'supported']]
        assessment=[k for k in test if core.loc[k,'supported'] and wall.loc[k,'supported']]
        assert len(training)==z.training_joint_blocks and len(assessment)==z.test_joint_blocks
        assert bool(len(training)>=3 and len(assessment)>=3)==z.supported
        assert ','.join(map(str,assessment))==('' if pd.isna(z.accepted_block_ids) else str(z.accepted_block_ids))
        assert not z.material_phase_pair
        def mean(region,ids):
            if not ids:return np.nan
            w=core.loc[ids,'intervals'].to_numpy(float)
            return float(np.dot(groups[region].loc[ids,'velocity_px_frame'],w)/w.sum()*z.fps)
        check=dict(test_intervals_accepted=int(core.loc[assessment,'intervals'].sum()),
            test_intervals_total=int(core.loc[test,'intervals'].sum()),
            core_signed_px_s=mean('core',assessment),wall_signed_px_s=mean(choice,assessment),
            all_test_core_candidate_px_s=mean('core',test),all_test_wall_candidate_px_s=mean(choice,test))
        for key,value in check.items():close(getattr(z,key),value,'training_only:'+key)
    # Underlying correlations: independent scipy.signal linear cross-correlation
    # and explicit Pearson products at selected shifts/pairs.
    max_pearson_error=0.;pearson_values=0;curves_checked=0;max_peak_error=0.;nulls_checked=0
    def independent_curve(sequence,lag=1,source=240,direct=False):
        nonlocal max_pearson_error,pearson_values,curves_checked,max_peak_error
        a=sequence[:-lag];b=sequence[lag:];width=536
        template=a[:,source:source+width]
        centred=template-template.mean(axis=1)[:,None]
        ss=np.cumsum(b,axis=1);ss=np.column_stack([np.zeros(len(b)),ss])
        ss2=np.cumsum(b*b,axis=1);ss2=np.column_stack([np.zeros(len(b)),ss2])
        target_sums=ss[:,width:]-ss[:,:-width]
        target_energy=np.maximum(ss2[:,width:]-ss2[:,:-width]-target_sums**2/width,0)
        denom=np.sqrt(np.sum(centred**2,axis=1)[:,None]*target_energy)
        covariance=np.stack([correlate(bb,tt,mode='valid',method='fft') for bb,tt in zip(b,centred)])
        curves=np.divide(covariance,denom,out=np.full_like(covariance,np.nan),where=denom>1e-9)
        valid=np.minimum(a.std(axis=1),b.std(axis=1))>=1
        mean=np.nanmean(curves[valid],axis=0) if valid.any() else np.full(481,np.nan)
        if not np.isfinite(mean).any():return dict(speed=np.nan,peak=np.nan,boundary=True,coverage=float(valid.mean()))
        k=int(np.nanargmax(mean));offset=0.
        if 0<k<480:
            curvature=mean[k-1]-2*mean[k]+mean[k+1]
            if np.isfinite(curvature) and curvature<-1e-12:
                offset=float(np.clip((mean[k-1]-mean[k+1])/(2*curvature),-.5,.5))
        if direct:
            for i in (0,len(a)//2,len(a)-1):
                for j in sorted(set([0,160,240,320,480,max(k-1,0),k,min(k+1,480)])):
                    x=template[i]-template[i].mean()
                    y=b[i,j:j+width]-b[i,j:j+width].mean()
                    denominator=np.linalg.norm(x)*np.linalg.norm(y)
                    if denominator>1e-9:
                        actual=float(x@y/denominator)
                        error=abs(actual-curves[i,j]);max_pearson_error=max(max_pearson_error,error)
                        assert error<1e-9,(source,lag,i,j,error)
                        pearson_values+=1
        curves_checked+=1
        return dict(speed=(k-source+offset)/lag,peak=float(mean[k]),boundary=k<=2 or k>=478,coverage=float(valid.mean()))
    # Representative original and sensitivity cases: all3 regions, diverse
    # patterns, saturated annular, problematic churn/slug, both time halves.
    sample=[]
    target_runs=('m03q05_001','m03q40_001','m04q01_001','m05q01_001','m07q01_001')
    for run in target_runs:
        for region in REGIONS:
            for block in (0,5,9):sample.append((run,region,'blocks10',block))
    for variant,block in [('raw_blocks10',5),('highpass8_blocks10',5),('highpass32_blocks10',5),('blocks5',2),('blocks20',10)]:
        for run in ('m03q05_001','m03q40_001','m04q01_001'):
            sample.append((run,'core',variant,block))
    cached={}
    for run,region,variant,block in sample:
        key=(run,region,variant)
        if key not in cached:
            raw=archive[f'{run}__{region}'].astype(float)
            sequence=raw if variant=='raw_blocks10' else raw-raw.mean(axis=0)
            if variant.startswith('highpass'):
                sigma=8 if variant.startswith('highpass8_') else 32
                sequence-=gaussian_filter1d(sequence,sigma,axis=1)
            cached[key]=sequence
        row=detail[(detail.run==run)&(detail.region==region)&(detail.variant==variant)&(detail.block==block)].iloc[0]
        sequence=cached[key][int(row.first_frame):int(row.last_frame)+1]
        jobs=[('',sequence,1,240),('lag2',sequence,2,240),('reverse',sequence[::-1],1,240),
              ('x160',sequence,1,160),('x320',sequence,1,320)]
        for name,seq,lag,source in jobs:
            recovered=independent_curve(seq,lag,source,direct=True)
            field='velocity_px_frame' if not name else name+'_velocity_px_frame'
            error=abs(recovered['speed']-row[field]);max_peak_error=max(max_peak_error,error)
            assert error<2e-8,(run,region,variant,block,name,error)
            prefix=name+'_' if name else ''
            close(row[prefix+'peak'],recovered['peak'],'direct_peak',tol=1e-9)
            assert row[prefix+'boundary']==recovered['boundary']
            close(row[prefix+'texture_coverage'],recovered['coverage'],'direct_texture',tol=1e-12)
        # Regenerate exactly the9 saved permutations for each sampled block.
        index=int(scope.index[scope.run==run][0])
        rng=np.random.default_rng(SEED+index*100+block)
        null=[]
        for _ in range(9):null.append(independent_curve(sequence[rng.permutation(len(sequence))])['peak'])
        close(row.null_max_peak,float(np.max(null)),'recomputed_shuffle',tol=1e-9)
        nulls_checked+=9
    split_sample=[]
    for run in target_runs:
        for region in REGIONS:
            for direction in ('forward','reverse'):
                for block in (0,5):split_sample.append((run,region,'training_background',direction,block))
    for run in ('m03q05_001','m03q40_001','m04q01_001'):
        for variant in ('training_background_hp8','training_background_hp32'):
            for direction in ('forward','reverse'):
                split_sample.append((run,'core',variant,direction,5 if direction=='forward' else 0))
    split_cache={}
    for run,region,variant,direction,block in split_sample:
        key=(run,region,variant,direction)
        if key not in split_cache:
            raw=archive[f'{run}__{region}'].astype(float)
            training=raw[:251] if direction=='forward' else raw[251:]
            sequence=raw-training.mean(axis=0)
            if variant.endswith('_hp8') or variant.endswith('_hp32'):
                sigma=8 if variant.endswith('_hp8') else 32
                sequence-=gaussian_filter1d(sequence,sigma,axis=1)
            split_cache[key]=sequence
        row=split_detail[(split_detail.run==run)&(split_detail.region==region)&(split_detail.variant==variant)&
                         (split_detail.selection_direction==direction)&(split_detail.block==block)].iloc[0]
        sequence=split_cache[key][int(row.first_frame):int(row.last_frame)+1]
        for name,seq,lag,source in [('',sequence,1,240),('lag2',sequence,2,240),
                ('reverse',sequence[::-1],1,240),('x160',sequence,1,160),('x320',sequence,1,320)]:
            recovered=independent_curve(seq,lag,source,direct=True)
            field='velocity_px_frame' if not name else name+'_velocity_px_frame'
            error=abs(recovered['speed']-row[field]);max_peak_error=max(max_peak_error,error)
            assert error<2e-8,(run,region,variant,direction,block,name,error)
            prefix=name+'_' if name else ''
            close(row[prefix+'peak'],recovered['peak'],'training_only_peak',tol=1e-9)
            assert row[prefix+'boundary']==recovered['boundary']
            close(row[prefix+'texture_coverage'],recovered['coverage'],'training_only_texture',tol=1e-12)
        index=int(scope.index[scope.run==run][0])
        rng=np.random.default_rng(SEED+index*100+block)
        null=[independent_curve(sequence[rng.permutation(len(sequence))])['peak'] for _ in range(9)]
        close(row.null_max_peak,float(np.max(null)),'training_only_shuffle',tol=1e-9)
        nulls_checked+=9
    # Interpretation-sensitive summaries, not additional acceptance rules.
    a=selected[(selected.variant=='blocks10')&(selected.selection_direction=='forward')].set_index('run')
    z=selected[(selected.variant=='blocks10')&(selected.selection_direction=='reverse')].set_index('run')
    both=a.supported&z.supported
    differences={}
    for field in ('core_signed_px_s','wall_signed_px_s'):
        changed=both&(abs(a[field]-z[field])>np.maximum(2*a.fps,.2*abs(a[field])))
        differences[field]=a.index[changed].tolist()
    block_counts=perrun.groupby('variant').supported.sum().astype(int).to_dict()
    selected_counts={v+'__'+d:int(k) for (v,d),k in selected.groupby(['variant','selection_direction']).supported.sum().items()}
    split_counts={v+'__'+d:int(k) for (v,d),k in split_pairs.groupby(['variant','selection_direction']).supported.sum().items()}
    final=split_pairs[(split_pairs.variant=='training_background')&(split_pairs.selection_direction=='forward')].set_index('run')
    reverse=split_pairs[(split_pairs.variant=='training_background')&(split_pairs.selection_direction=='reverse')].set_index('run')
    alternatives={'reverse':reverse,
        'raw':selected[(selected.variant=='raw_blocks10')&(selected.selection_direction=='forward')].set_index('run'),
        'hp32':split_pairs[(split_pairs.variant=='training_background_hp32')&(split_pairs.selection_direction=='forward')].set_index('run'),
        'hp8':split_pairs[(split_pairs.variant=='training_background_hp8')&(split_pairs.selection_direction=='forward')].set_index('run')}
    stable=final.supported.copy();sensitivity={}
    for label,other in alternatives.items():
        both=final.supported&other.supported
        agree=both.copy()
        for field in ('core_signed_px_s','wall_signed_px_s'):
            agree &= abs(final[field]-other[field])<=np.maximum(2*final.fps,.2*abs(final[field]))
        stable &= agree
        sensitivity[label]=dict(both_supported=int(both.sum()),both_coordinates_agree=int(agree.sum()),
            changed_or_unsupported_runs=final.index[final.supported&~agree].tolist())
    hashes={str(path.relative_to(BASE)):hashlib.sha256(path.read_bytes()).hexdigest() for path in
        [BASE/'scripts/block_ncc.py',BASE/'scripts/selected_wall.py',BASE/'PROTOCOL.md',
         BASE/'results/block_ncc_details.csv',BASE/'results/block_ncc_per_run.csv',BASE/'results/selected_wall_per_run.csv']}
    for path in [BASE/'scripts/split_background_ncc.py',BASE/'results/split_background_details.csv',
                 BASE/'results/split_background_pairs.csv',BASE/'results/split_background_metadata.json']:
        hashes[str(path.relative_to(BASE))]=hashlib.sha256(path.read_bytes()).hexdigest()
    assert hashlib.sha256((BASE/'scripts/split_background_ncc.py').read_bytes()).hexdigest()==split_script_hash
    result=dict(status='passed_with_interpretation_limits',runtime_s=time.monotonic()-start,
        numerical_errors_found=0,detail_rows_checked=len(detail),summary_rows_checked=len(perrun),
        selection_rows_checked=len(selected),numeric_fields_checked=numeric_fields,max_numeric_difference=max_summary_error,
        final_split_detail_rows_checked=len(split_detail),final_split_pair_rows_checked=len(split_pairs),
        final_background_train_assessment_disjoint=True,
        representative_blocks_recomputed=len(sample)+len(split_sample),averaged_curves_recomputed=curves_checked,
        direct_centered_Pearson_values=pearson_values,max_direct_Pearson_difference=max_pearson_error,
        max_recomputed_velocity_difference_px_frame=max_peak_error,shuffled_controls_recomputed=nulls_checked,
        input_sources=sources,block_counts=block_counts,selected_counts=selected_counts,
        final_split_counts=split_counts,final_sensitivity=sensitivity,
        final_all_alternatives_stable_count=int(stable.sum()),final_all_alternatives_stable_runs=final.index[stable].tolist(),
        final_forward_reverse_wall_changes=final.index[final.selected_wall!=reverse.selected_wall].tolist(),
        final_forward_interval_fraction_of_full_record=[float(final.loc[final.supported,'test_intervals_accepted'].min()/500),float(final.loc[final.supported,'test_intervals_accepted'].max()/500)],
        forward_reverse_selected_wall_changes=a.index[a.selected_wall!=z.selected_wall].tolist(),
        forward_reverse_differences_over_max_2px_or20percent=differences,
        forward_reported_interval_fraction_of_full_record=[float(a.loc[a.supported,'test_intervals_accepted'].min()/500),float(a.loc[a.supported,'test_intervals_accepted'].max()/500)],
        warnings=[
            'The ARCHIVED selected_wall experiment shared the501-frame mean; FINAL split_background uses training-only mean, frozen for assessment. That processing leakage is fixed.',
            'Selection was developed after inspecting the recordings; no external experimental validation.',
            'Selected-wall means describe accepted assessment-half blocks and one chosen wall, not whole-record time averages or mean liquid-film speed.',
            'm03q05 has finite forward coordinates from one accepted block but fails training2/5 and assessment1/5 support; never mark this point accepted.',
            'm03q40 passes selected-wall reference but fails raw and both high-pass alternatives; retaining it must state preprocessing dependence and original saturation.',
            'High-pass failures expose frequency/content sensitivity and do not alone establish that reference velocities are physically false.',
            'Passing numerical and transport tests does not identify gas/liquid material or superficial velocities.',
            'The other-wall candidate on assessment epochs need not pass its own support gate and must not be presented as a validated second-wall measurement.'
        ],input_sha256=hashes)
    (BASE/'results/independent_block_verification.json').write_text(json.dumps(result,indent=2,allow_nan=False))
    report=f'''# Revisão independente — correlação por blocos e parede selecionada

**A aritmética, os critérios salvos, a seleção e as médias foram reproduzidos sem erro detectado.** A versão final corrigida, com fundo estimado somente na metade de seleção, contém **25/26 pares ópticos** com uma parede selecionada, em cada direção. **22/25** concordam ao inverter as metades e **7/26** passam e concordam em todas as alternativas declaradas. O resultado de referência anterior, usando ambas as paredes, contém20/26 pares. Nenhum desses números estabelece26 pares aceitos ou velocidades materiais de gás/líquido.

## Correção de compartilhamento entre as metades

A auditoria inicial identificou que `selected_wall_per_run.csv` utilizava perfis com fundo estimado pela média das501 imagens. Esse experimento foi preservado como **arquivado, com fundo compartilhado**. A análise final está em `split_background_details.csv` e `split_background_pairs.csv`: na direção inicial→final, o fundo usa somente os quadros0–250 e é congelado para avaliar251–500; na direção oposta, usa somente251–500 e avalia0–250. A seleção da parede também usa exclusivamente sua metade de seleção.

Foi conferido que não existe interseção entre os quadros usados para o fundo e os de avaliação. Os perfis, picos e controles foram recalculados por implementação independente com o fundo correto. **O compartilhamento do fundo foi corrigido na versão final.** Isso ainda não constitui validação experimental externa: o desenho do método explorou essas mesmas gravações. A filtragem8/32 foi repetida com o fundo restrito à seleção, sem transportar automaticamente os resultados anteriores.

## Escopo da conferência

Foram reconstruídas as decisões dos **{len(detail)} registros regionais/blocos**, as **{len(perrun)} linhas por gravação/alternativa** e as **{len(selected)} linhas de seleção de parede** diretamente do CSV detalhado, sem importar os geradores. As regras conferidas incluem textura, picos, controle embaralhado, atraso2, reversão, posições160/320, tolerância, índices de quadros, blocos simultâneos, cobertura temporal, pesos por intervalos, escolha no conjunto de seleção e médias apenas nos blocos aceitos de avaliação.

Adicionalmente, foram reconstruídas **todas as {len(split_detail)} decisões regionais/blocos e {len(split_pairs)} seleções da versão final corrigida**, incluindo as duas direções e três alternativas de processamento. A correção não alterou as contagens25/25 da referência, mas suas coordenadas e verificações foram recalculadas e conferidas; não se presumiu equivalência.

FPS foi comparado com dois campos nativos em cada uma das26 gravações; foram conferidos501 nomes consecutivos e perfis de três imagens nativas por gravação. A média das duas paredes utiliza blocos simultâneos e pesos idênticos: a média das duas médias é matematicamente correta para esse descritor regional.

Recalcularam-se **{len(sample)+len(split_sample)} blocos representativos** ({len(sample)} arquivados e {len(split_sample)} da versão final), incluindo m03q05, m03q40, m04q01, m05q01 e m07q01, com outra implementação de correlação linear e normalização de Pearson. Foram verificadas **{pearson_values} avaliações de Pearson por produtos centrados explícitos**, **{curves_checked} curvas médias** e **{nulls_checked} controles embaralhados regenerados**. A maior diferença de Pearson foi {max_pearson_error:.3g}; a maior diferença de velocidade foi {max_peak_error:.3g} px/quadro. Esses números verificam implementação, não erro experimental.

## Resultados reproduzidos

| Alternativa de uma parede | Seleção inicial→final | Seleção final→inicial |
|---|---:|---:|
| **Final: fundo somente da seleção**,10 blocos | {split_counts['training_background__forward']} | {split_counts['training_background__reverse']} |
| Sem subtração temporal | {selected_counts['raw_blocks10__forward']} | {selected_counts['raw_blocks10__reverse']} |
| Fundo somente da seleção + passa-altas32 | {split_counts['training_background_hp32__forward']} | {split_counts['training_background_hp32__reverse']} |
| Fundo somente da seleção + passa-altas8 | {split_counts['training_background_hp8__forward']} | {split_counts['training_background_hp8__reverse']} |
| Arquivado: fundo compartilhado501 quadros | {selected_counts['blocks10__forward']} | {selected_counts['blocks10__reverse']} |

As alternativas de cinco e vinte blocos forneceram {block_counts['blocks5']} e {block_counts['blocks20']} pares com as duas paredes. A cobertura maior da parede selecionada acompanha uma mudança explícita da grandeza: uma região escolhida por observabilidade substitui a média de ambas.

## Pontos que precisam aparecer no resumo

1. **m03q05 continua rejeitado.** Na seleção inicial→final, a parede inferior foi escolhida, mas somente2/5 blocos conjuntos na seleção e1/5 na avaliação passaram. O ponto finito usa apenas esse bloco de49 intervalos. Ele pode aparecer como candidato rejeitado; não completa26 pontos aceitos.
2. **A versão final separa o fundo e a escolha da parede entre seleção/avaliação, mas não é validação externa.** O vazamento de fundo existia no experimento arquivado e foi corrigido. A escolha do desenho analítico após explorar o conjunto permanece como limitação; não se deve afirmar que os25 resultados foram confirmados em um conjunto experimental independente.
3. **As coordenadas reportadas pertencem à metade de avaliação e a seus blocos aceitos.** Na direção inicial→final, elas usam29,4%–49,0% dos500 intervalos do vídeo. Isso não é a média integral do ensaio. Núcleo e parede usam os mesmos blocos, o que evita uma mistura temporal entre as duas coordenadas.
4. **Inversão das metades e processamento:** na versão final, **{sensitivity['reverse']['both_coordinates_agree']}/25** pares concordaram em ambas as coordenadas ao inverter as metades, dentro de max(2px/quadro,20% do módulo da referência). Os casos discordantes foram {', '.join(sensitivity['reverse']['changed_or_unsupported_runs'])}. Exigindo aprovação e concordância também com a variante sem fundo e com as duas filtragens, restaram **{int(stable.sum())}** gravações: {', '.join(final.index[stable])}. Isso combina variabilidade temporal, mudança de região, escala óptica e mudança dos blocos selecionados; não é um erro experimental calibrado.
5. **m03q40 depende do processamento.** Ele passa na referência de parede selecionada e falha na variante sem subtração e em ambas as variantes passa-altas. A saturação original e essa dependência precisam acompanhar seu ponto. Os controles sintéticos de recorte mostram possibilidade de recuperar movimentos sob saturação, mas não validam automaticamente esta gravação.
6. **As falhas passa-altas não provam velocidades fisicamente falsas.** Alterar a escala óptica pode mudar a estrutura dominante, sua textura e a seleção dos blocos. A evidência limita a representatividade de um único valor e precisa ser mantida como sensibilidade.
7. **Não houve identificação material.** Uma correlação coerente pode acompanhar ondas, interfaces ou brilho compartilhado entre as regiões. A localização núcleo/parede e a escolha da região mais observável não fornecem, por si, gás/líquido material.
8. **A outra parede foi salva como candidata no experimento arquivado.** `other_wall_same_epochs_candidate_px_s` não exige suporte dessa outra parede nos mesmos blocos. Não deve ser apresentada como uma segunda medida de parede validada nem confundida com a parede da versão final.

## Disposição

Os cálculos finais, com fundo restrito à metade de seleção, podem ser utilizados como **descritores de propagação óptica sob o protocolo declarado**, com25 pares na referência exploratória e um candidato explicitamente rejeitado. A sensibilidade ampliada retém apenas7 gravações aprovadas e concordantes em todas as alternativas; os25 pontos não devem ser apresentados como invariantes ao método. O gráfico completo precisa preservar o status do26º ponto, o caráter regional dos eixos, a escolha de uma parede e o período efetivamente amostrado. Nenhuma das conferências converte esses resultados em26 pares de velocidades materiais.

Reprodução da auditoria: `python3 -B src2/v15.1_gas_liquid_velocity/all_scenarios_adversarial/scripts/independent_block_verify.py`. Evidência estruturada: [independent_block_verification.json](../results/independent_block_verification.json).
'''
    (BASE/'findings/block_independent_review.md').write_text(report)
    print(json.dumps({k:v for k,v in result.items() if k not in ('input_sources','input_sha256','warnings')},indent=2))

if __name__=='__main__':audit()
