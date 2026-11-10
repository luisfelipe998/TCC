#!/usr/bin/env python3
"""Additional negative controls and independent aggregation of new methods."""
from pathlib import Path
import hashlib
import importlib.util
import json
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

BASE=Path(__file__).resolve().parents[1]
OUT=BASE/'results'
SEED=20261024

def module(name):
    path=BASE/'scripts'/name
    spec=importlib.util.spec_from_file_location(name.replace('.py',''),path)
    loaded=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded

def clean(v):
    if isinstance(v,dict):return {str(k):clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple,np.ndarray)):return [clean(x) for x in v]
    if isinstance(v,(bool,np.bool_)):return bool(v)
    if isinstance(v,(int,np.integer)):return int(v)
    if isinstance(v,(float,np.floating)):return float(v) if np.isfinite(v) else None
    return v

def nuisance_controls():
    phase=module('phase_local_experiment.py')
    block=module('block_ncc.py')
    rng=np.random.default_rng(SEED)
    n,h,w=101,20,1016
    x=np.arange(w);times=np.arange(n)
    gas=np.zeros((n,h,w),bool);gas[:,:10]=True;liquid=~gas
    texture=gaussian_filter1d(rng.normal(size=12000),1.5)*25
    illumination=np.stack([np.interp(x-12*t+3000,np.arange(12000),texture) for t in times])
    images=(np.where(gas,160.,70.)+illumination[:,None,:]).astype(np.float32)
    records=[]
    for label,mask in [('gas_candidate',gas),('liquid_candidate',liquid)]:
        result,_,_,_=phase.analyze(images,mask,nulls=True)
        assert result['optical_supported'] and abs(result['signed_px_frame']-12)<.01
        records.append(dict(control='common_moving_illumination_without_material_texture',phase=label,
            illumination_speed_px_frame=12.,estimated_px_frame=result['signed_px_frame'],
            optical_supported=result['optical_supported'],peak=result['peak'],
            null_max_peak=result['null_max_peak'],
            interpretation='Correct optical speed despite no material tracer. This is not a phase-velocity validation or an algorithm failure.'))
    pd.DataFrame(records).to_csv(OUT/'adversarial_shared_illumination.csv',index=False)
    periodic=[]
    t=np.arange(51)[:,None]
    for period in (24,32,48,64,80,100):
        for truth in (-12.,-35.,22.75):
            sequence=50*np.sin(2*np.pi*(x[None,:]-truth*t)/period)+15*np.cos(4*np.pi*(x[None,:]-truth*t)/period+.3)
            sequence-=sequence.mean(axis=0)
            result=block.local(sequence,SEED)
            periodic.append(dict(control='ambiguous_periodic_translation',period_px=period,
                truth_px_frame=truth,estimated_px_frame=result['velocity_px_frame'],
                supported=result['supported'],reasons=result['reasons']))
    p=pd.DataFrame(periodic)
    p.to_csv(OUT/'adversarial_periodic_block_ncc.csv',index=False)
    return dict(illumination=records,periodic_cases=len(p),periodic_rejected=int((~p.supported).sum()),
        interpretation='These tested periodic aliases are rejected by the joint temporal/reversal/position gates. This finite test is not a general anti-aliasing proof.')

def independent_aggregation():
    results={}
    details=pd.read_csv(OUT/'block_ncc_details.csv')
    saved=pd.read_csv(OUT/'block_ncc_per_run.csv').set_index(['run','variant'])
    errors=[];support=[]
    for (name,variant),group in details.groupby(['run','variant']):
        velocity=group.pivot(index='block',columns='region',values='velocity_px_frame')
        flags=group.pivot(index='block',columns='region',values='supported').all(axis=1)
        temporal=group[group.region=='core'].set_index('block').sort_index()
        weights=temporal.intervals.to_numpy()
        good=flags.to_numpy()
        ids=temporal.index.to_numpy()[good]
        fifths=len(set(np.minimum(((temporal.first_frame.to_numpy()[good]+temporal.last_frame.to_numpy()[good])/2/501*5).astype(int),4)))
        eligible=good.mean()>=.6 and fifths>=4
        row=saved.loc[(name,variant)]
        assert bool(row.supported)==bool(eligible)
        if good.any():
            core=np.average(velocity.core.to_numpy()[good],weights=weights[good])*row.fps
            walls=np.average(velocity[['upper_outer','lower_outer']].mean(axis=1).to_numpy()[good],weights=weights[good])*row.fps
            errors.extend([abs(core-row.core_signed_px_s),abs(walls-row.wall_signed_px_s)])
        support.append(dict(run=name,variant=variant,supported=bool(eligible),selected_block_ids=';'.join(map(str,ids))))
    pd.DataFrame(support).to_csv(OUT/'adversarial_block_independent_status.csv',index=False)
    results['block_ncc']=dict(rows=len(saved),max_coordinate_difference_px_s=max(errors),
        support_by_variant=pd.DataFrame(support).groupby('variant').supported.sum().astype(int).to_dict())
    assert max(errors)<1e-7
    blocks=pd.read_csv(OUT/'spectral_blocks.csv')
    pairs=pd.read_csv(OUT/'spectral_pairs.csv')
    errors=[];supported=0;allfive=0
    assert 'joint_supported' in blocks and 'accepted_core_signed_px_s' in pairs
    for row in pairs.itertuples():
        group=blocks[blocks.run==row.run]
        speeds=group.pivot(index='block',columns='region',values='signed_px_s')
        joint=group.pivot(index='block',columns='region',values='joint_supported').all(axis=1)
        count=int(joint.sum())
        assert bool(row.robust_supported)==(count>=4)
        supported+=count>=4;allfive+=count==5
        errors.extend([abs(speeds.core.mean()-row.core_signed_px_s),
                       abs(speeds[['upper_outer','lower_outer']].to_numpy().mean()-row.wall_signed_px_s)])
        if count:
            errors.extend([abs(speeds.loc[joint,'core'].mean()-row.accepted_core_signed_px_s),
                           abs(speeds.loc[joint,['upper_outer','lower_outer']].to_numpy().mean()-row.accepted_wall_signed_px_s)])
    assert max(errors)<1e-7
    results['spectral']=dict(rows=len(pairs),max_coordinate_difference_px_s=max(errors),
        simultaneous_supported_pairs=supported,all_five_joint_blocks_pairs=allfive)
    summaries=pd.read_csv(OUT/'phase_local_summary.csv')
    lags=pd.read_csv(OUT/'phase_local_lags.csv')
    blocks=pd.read_csv(OUT/'phase_local_blocks.csv')
    nulls=pd.read_csv(OUT/'phase_local_nulls.csv')
    pairtable=pd.read_csv(OUT/'phase_local_pairs.csv')
    errors=[]
    for row in summaries.itertuples():
        def select(table):
            return table[(table.run==row.run)&(table.erosion_radius==row.erosion_radius)&(table.phase==row.phase)]
        lg=select(lags);bl=select(blocks);nl=select(nulls)
        v=float(lg.loc[lg.lag==1,'shift'].iloc[0])
        tol=max(2.,.2*abs(v)) if np.isfinite(v) else np.nan
        lagok=bool((np.isfinite(lg.peak)&(lg.peak>=.2)&(lg['std']>=1)&~lg.boundary&
                    (abs(lg['shift']/lg.lag-v)<=tol)).all())
        one=bl[bl.lag==1]
        blockcount=int((np.isfinite(one.peak)&(one.peak>=.2)&(one['std']>=1)&~one.boundary&
                        (abs(one.velocity_px_frame-v)<=tol)).sum())
        nullmax=nl.peak.max()
        expected=bool(row.frames_with_100_phase_pixels>=20 and row.raw_phase_std>=1 and
            np.isfinite(v) and lagok and blockcount>=4 and np.isfinite(row.margin) and row.margin>=.02 and
            np.isfinite(nullmax) and row.peak>=nullmax+.05 and row.saturation_fraction<=.20)
        assert bool(row.optical_supported)==expected,(row.run,row.erosion_radius,row.phase)
        assert bool(row.lag_consistent)==lagok and int(row.agreeing_blocks)==blockcount
        if np.isfinite(v):errors.append(abs(row.signed_px_s-v*row.fps))
    for row in pairtable.itertuples():
        phases=summaries[(summaries.run==row.run)&(summaries.erosion_radius==row.erosion_radius)]
        assert len(phases)==2
        assert bool(row.pair_optical_supported)==bool(phases.optical_supported.all())
        assert bool(row.both_finite)==bool(np.isfinite(phases.signed_px_s).all())
        assert not row.material_pair_identified
    assert max(errors)<1e-7
    results['phase_local']=dict(region_rows=len(summaries),pair_rows=len(pairtable),
        finite_by_erosion_radius=pairtable.groupby('erosion_radius').both_finite.sum().astype(int).to_dict(),
        optical_pairs_by_erosion_radius=pairtable.groupby('erosion_radius').pair_optical_supported.sum().astype(int).to_dict(),
        max_fps_conversion_difference_px_s=max(errors),all_acceptance_rules_recomputed=True,
        independently_confirmed_material_pairs=0)
    return results

def main():
    result=dict(status='additional_controls_and_aggregations_checked',
        controls=nuisance_controls(),independent_aggregation=independent_aggregation())
    files=['scripts/phase_local_experiment.py','scripts/block_ncc.py','scripts/spectral_transport.py',
           'results/block_ncc_details.csv','results/block_ncc_per_run.csv',
           'results/spectral_blocks.csv','results/spectral_pairs.csv',
           'results/phase_local_summary.csv','results/phase_local_lags.csv',
           'results/phase_local_blocks.csv','results/phase_local_nulls.csv','results/phase_local_pairs.csv']
    result['sources_sha256']={name:hashlib.sha256((BASE/name).read_bytes()).hexdigest() for name in files}
    result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (OUT/'adversarial_new_methods.json').write_text(json.dumps(clean(result),indent=2,allow_nan=False)+'\n')
    print(json.dumps(clean(result),indent=2))

if __name__=='__main__':main()
