#!/usr/bin/env python3
"""Time-local Fourier transport, without gas/liquid material assignment.

Every choice below is fixed for every recording. This experiment changes the
observable to the equal-block mean of five local dominant optical speeds.
It does not require one constant speed over the entire recording.
"""
from pathlib import Path
from functools import lru_cache
import hashlib
import json
import time
import numpy as np
import pandas as pd
from scipy.fft import rfft2
from scipy.ndimage import gaussian_filter1d

BASE=Path(__file__).resolve().parents[1]
OLD=BASE.parent
VEL=np.arange(-80.,80.0001,.25)
BANDS={'main':(1/256,1/32),'coarse':(1/512,1/64),'fine':(1/128,1/16)}
SEED=20261024
NULLS=19
REGIONS=('core','upper_outer','lower_outer')

def serial(v):
    if isinstance(v,dict):return {k:serial(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [serial(x) for x in v]
    if isinstance(v,np.ndarray):return serial(v.tolist())
    if isinstance(v,np.bool_):return bool(v)
    if isinstance(v,(np.integer,)):return int(v)
    if isinstance(v,(float,np.floating)):return float(v) if np.isfinite(v) else None
    return v

@lru_cache(None)
def coordinates(nt,nx,band):
    fx=np.fft.rfftfreq(nx)
    indices=np.flatnonzero((fx>=BANDS[band][0])&(fx<=BANDS[band][1]))
    p=(-VEL[:,None]*fx[None,indices])%1*nt
    ip=np.floor(p).astype(int)%nt
    frac=p-np.floor(p)
    return indices,ip,frac

def spectrum(a):
    a=np.asarray(a,float)
    dynamic=a-a.mean(axis=0,keepdims=True)
    dynamic-=dynamic.mean(axis=1,keepdims=True)
    texture=float(dynamic.std())
    window=np.hanning(len(a))[:,None]*np.hanning(a.shape[1])[None,:]
    power=abs(rfft2(dynamic*window,workers=1))**2
    # Remove DC and its immediate temporal leakage neighbours.
    power[[0,1,-1]]=0.
    return power,texture

def ridge(power,nx,band='main'):
    nt=len(power);idx,ip,frac=coordinates(nt,nx,band)
    p=power[:,idx]
    total=p.sum(axis=0)
    good=total>max(total.max()*1e-10,1e-15)
    if not good.any():return dict(speed=np.nan,score=0.,margin=0.,boundary=True),np.zeros(len(VEL))
    # Sum power in a three-bin temporal tube around f_t=-v f_x, including aliases.
    norm=np.divide(p,total[None,:],out=np.zeros_like(p),where=total[None,:]>0)
    j=np.arange(len(idx))[None,:]
    value=np.zeros_like(ip,dtype=float)
    for off in (-1,0,1):
        value+=(1-frac)*norm[(ip+off)%nt,j]+frac*norm[(ip+off+1)%nt,j]
    weights=np.sqrt(total)*good
    curve=value@weights/weights.sum()
    k=int(np.argmax(curve));speed=VEL[k]
    # Piecewise-linear interpolation makes a discrete grid maximum appropriate.
    outside=abs(VEL-speed)>max(3.,.15*abs(speed))
    margin=float(curve[k]-curve[outside].max())
    return dict(speed=float(speed),score=float(curve[k]),margin=margin,boundary=k in (0,len(VEL)-1)),curve

def estimate(a,band='main'):
    p,tex=spectrum(a);r,c=ridge(p,a.shape[1],band);r['texture']=tex
    return r,c,p

def base_gate(r,null_max):
    reasons=[]
    if r['texture']<1.:reasons.append('low_dynamic_texture')
    if not np.isfinite(r['speed']) or abs(r['speed'])<1.:reasons.append('missing_or_near_zero')
    if r['boundary']:reasons.append('velocity_search_boundary')
    if r['score']<.10:reasons.append('weak_directional_energy')
    if r['score']<null_max+.03:reasons.append('weak_shuffle_separation')
    if r['margin']<.03:reasons.append('ambiguous_directional_peak')
    return reasons

def synthetic():
    rng=np.random.default_rng(SEED)
    nx=1016;nt=101;x=np.arange(nx);tt=np.arange(nt)[:,None]
    rows=[]
    def moving(speed,seed):
        rr=np.random.default_rng(seed);freq=np.arange(4,32)/nx
        phase=rr.uniform(-np.pi,np.pi,len(freq));amp=rr.uniform(.5,1.5,len(freq))
        return sum(a*np.cos(2*np.pi*f*(x[None,:]-speed*tt)+ph) for f,ph,a in zip(freq,phase,amp))*10
    for speed in (-51.,-31.8,-17.3,-6.7,0.,23.75):
        seq=moving(speed,int(1000+speed*10))
        r,_,_=estimate(seq)
        null=[]
        for _ in range(NULLS):null.append(estimate(seq[rng.permutation(nt)])[0]['score'])
        rev=estimate(seq[::-1])[0]
        rows.append(dict(test='known_transport',truth=speed,**r,error=abs(r['speed']-speed),
                         reverse_error=abs(rev['speed']+r['speed']),null_max=max(null),
                         supported=not base_gate(r,max(null))))
    # Two simultaneous speeds challenge the dominant-ridge interpretation.
    for name,a,b in [('equal_multiplier_competitors',1.,1.),('unequal_multiplier_competitors',1.,.5)]:
        seq=a*moving(-12.,234)+b*moving(-38.,456)
        r,c,_=estimate(seq)
        rows.append(dict(test=name,truth=-12.,other_truth=-38.,**r,
                         interpretation='dominant wave or ambiguity, not mean phase/material speed'))
    for name,seq in [('constant',np.ones((nt,nx))*40),
                     ('stationary_texture',np.tile(rng.normal(size=nx)*20,(nt,1))),
                     ('random_frames',rng.normal(size=(nt,nx))*20)]:
        r,_,_=estimate(seq);null=[]
        for _ in range(NULLS):null.append(estimate(seq[rng.permutation(nt)])[0]['score'])
        rows.append(dict(test=name,**r,null_max=max(null),supported=not base_gate(r,max(null))))
    # Independently moving optical layers: regional separation can recover two
    # optical speeds when separate texture really exists.
    layers=[moving(-35.,120),moving(-14.,121)]
    for n,seq in zip(('separate_core_layer','separate_wall_layer'),layers):
        truth=-35 if 'core' in n else -14
        r,_,_=estimate(seq);rows.append(dict(test=n,truth=truth,**r,error=abs(r['speed']-truth)))
    # Explicit nonstationarity: five blocks with different known velocities.
    truth=np.array([-12.,-20.,-30.,-18.,-25.])
    fitted=np.array([estimate(moving(v,100+i))[0]['speed'] for i,v in enumerate(truth)])
    rows.append(dict(test='piecewise_speed_time_mean',truth=float(truth.mean()),speed=float(fitted.mean()),
                     error=abs(fitted.mean()-truth.mean()),block_truth=truth.tolist(),block_estimates=fitted.tolist()))
    return rows

def main():
    start=time.monotonic()
    controls=synthetic()
    pd.DataFrame(controls).to_csv(BASE/'results/spectral_synthetic.csv',index=False)
    check=[r for r in controls if r['test']=='known_transport' and r['truth']!=0]
    assert all(r['error']<=.5 and r['reverse_error']<=.25 for r in check),check
    assert all(not r['supported'] for r in controls if r['test'] in ('constant','stationary_texture','random_frames'))
    scope=pd.read_csv(OLD/'results/all_run_status.csv')
    regions=pd.read_csv(OLD/'results/annular_region_summary.csv').set_index(['run','region'])
    archive=np.load(OLD/'intermediate/annular_profiles.npz')
    blocks=[];null_rows=[];summaries=[];curves={}
    rng=np.random.default_rng(SEED)
    for run in scope.run:
        for region in REGIONS:
            seq=archive[f'{run}__{region}'].astype(float)
            old=regions.loc[(run,region)];fps=float(old.fps)
            br=[]
            for b,ids in enumerate(np.array_split(np.arange(501),5)):
                a=seq[ids];r,curve,p=estimate(a)
                null=[]
                for trial in range(NULLS):
                    nr=estimate(a[rng.permutation(len(a))])[0]
                    null.append(nr['score'])
                    null_rows.append(dict(run=run,region=region,block=b,trial=trial,score=nr['score'],speed=nr['speed']))
                rev=estimate(a[::-1])[0]
                coarse=ridge(p,a.shape[1],'coarse')[0];fine=ridge(p,a.shape[1],'fine')[0]
                left=estimate(a[:,:508])[0];right=estimate(a[:,508:])[0]
                tol=max(3.,.25*abs(r['speed']))
                reasons=base_gate(r,max(null))
                reverse_ok=abs(r['speed']+rev['speed'])<=.25
                band_ok=abs(coarse['speed']-r['speed'])<=tol and abs(fine['speed']-r['speed'])<=tol
                axial_ok=abs(left['speed']-r['speed'])<=tol and abs(right['speed']-r['speed'])<=tol
                if not reverse_ok:reasons.append('reversal_inconsistent')
                row=dict(run=run,region=region,block=b,first_frame=int(ids[0]),last_frame=int(ids[-1]),
                         fps=fps,signed_px_frame=r['speed'],signed_px_s=r['speed']*fps,
                         score=r['score'],margin=r['margin'],texture=r['texture'],boundary=r['boundary'],
                         shuffle_max=max(null),reverse_signed_px_frame=rev['speed'],reverse_ok=reverse_ok,
                         coarse_px_frame=coarse['speed'],fine_px_frame=fine['speed'],band_ok=band_ok,
                         left_px_frame=left['speed'],right_px_frame=right['speed'],axial_ok=axial_ok,
                         tolerance_px_frame=tol,transport_supported=not reasons,reasons=';'.join(reasons))
                row['saturation_fraction']=float(old.saturation_fraction)
                row['joint_supported']=bool(not reasons and band_ok and axial_ok and old.saturation_fraction<=.20)
                blocks.append(row);br.append(row)
                curves[f'{run}__{region}__block{b}']=curve.astype(np.float32)
            velocities=np.array([r['signed_px_frame'] for r in br])
            supported=sum(r['transport_supported'] for r in br)
            band=sum(r['band_ok'] for r in br);axial=sum(r['axial_ok'] for r in br)
            direction=max(sum(velocities<0),sum(velocities>0))
            sat=float(old.saturation_fraction)
            reasons=[]
            if supported<4:reasons.append('fewer_than_four_supported_blocks')
            if band<4:reasons.append('frequency_band_sensitivity')
            if axial<4:reasons.append('axial_half_sensitivity')
            if direction<4:reasons.append('mixed_propagation_directions')
            if sat>.20:reasons.append('region_saturation')
            joint=[r for r in br if r['joint_supported']]
            summaries.append(dict(run=run,pattern=old['pattern'],region=region,fps=fps,
                signed_px_frame=float(velocities.mean()),signed_px_s=float(velocities.mean()*fps),
                block_min_px_frame=float(velocities.min()),block_max_px_frame=float(velocities.max()),
                supported_blocks=supported,band_consistent_blocks=band,axial_consistent_blocks=axial,
                same_direction_blocks=direction,saturation_fraction=sat,
                temporal_supported=supported>=4 and direction>=4,
                marginal_robust_supported=not reasons,joint_supported_blocks=len(joint),
                joint_accepted_mean_px_s=float(np.mean([r['signed_px_s'] for r in joint])) if joint else np.nan,
                joint_reasons='' if len(joint)>=4 and direction>=4 else 'fewer_than_four_joint_blocks_or_mixed_sign',
                robust_supported=len(joint)>=4 and direction>=4,reasons=';'.join(reasons)))
        print(run,flush=True)
    blockdf=pd.DataFrame(blocks);summary=pd.DataFrame(summaries)
    blockdf.to_csv(BASE/'results/spectral_blocks.csv',index=False)
    pd.DataFrame(null_rows).to_csv(BASE/'results/spectral_nulls.csv',index=False)
    summary.to_csv(BASE/'results/spectral_regions.csv',index=False)
    np.savez_compressed(BASE/'results/spectral_score_curves.npz',velocities=VEL,**curves)
    pairs=[]
    for run,g in summary.groupby('run',sort=False):
        d=g.set_index('region');c=d.loc['core'];u=d.loc['upper_outer'];l=d.loc['lower_outer']
        wall=(u.signed_px_s+l.signed_px_s)/2
        joint=blockdf[blockdf.run==run].groupby('block').joint_supported.all()
        temporal_joint=blockdf[blockdf.run==run].groupby('block').transport_supported.all()
        good_blocks=joint.index[joint].tolist()
        selected=blockdf[(blockdf.run==run)&blockdf.block.isin(good_blocks)]
        core_mean=selected[selected.region=='core'].signed_px_s.mean()
        wall_mean=selected[selected.region.isin(['upper_outer','lower_outer'])].signed_px_s.mean()
        pairs.append(dict(run=run,pattern=c['pattern'],fps=c.fps,core_signed_px_s=c.signed_px_s,
            wall_signed_px_s=wall,upper_signed_px_s=u.signed_px_s,lower_signed_px_s=l.signed_px_s,
            upper_lower_difference_px_s=abs(u.signed_px_s-l.signed_px_s),
            core_magnitude_px_s=abs(c.signed_px_s),wall_magnitude_px_s=abs(wall),
            marginal_temporal_supported=bool(g.temporal_supported.all()),
            temporal_joint_blocks=int(temporal_joint.sum()),temporal_supported=int(temporal_joint.sum())>=4,
            marginal_robust_supported=bool(g.marginal_robust_supported.all()),
            joint_supported_blocks=len(good_blocks),joint_block_ids=';'.join(map(str,good_blocks)),
            accepted_core_signed_px_s=core_mean,accepted_wall_signed_px_s=wall_mean,
            robust_supported=len(good_blocks)>=4,full_candidate_mean_supported=len(good_blocks)==5,
            joint_reasons='' if len(good_blocks)>=4 else 'fewer_than_four_simultaneous_joint_blocks',
            reasons='|'.join(r.region+':'+r.reasons for r in g.itertuples() if r.reasons),
            observable='Equal-block mean of five local dominant optical speeds; wall averages two separately defined strips; no phase material assignment'))
    pd.DataFrame(pairs).to_csv(BASE/'results/spectral_pairs.csv',index=False)
    metadata=dict(runtime_s=time.monotonic()-start,seed=SEED,runs=len(scope),regions=len(summary),blocks=len(blocks),
        null_trials=len(null_rows),temporal_pair_support=sum(r['temporal_supported'] for r in pairs),
        robust_pair_support=sum(r['robust_supported'] for r in pairs),
        profile_archive_sha256=hashlib.sha256((OLD/'intermediate/annular_profiles.npz').read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        parameters=dict(velocity_min=-80,velocity_max=80,velocity_step=.25,frequency_bands=BANDS,
            five_blocks='consecutive, disjoint 101/100/100/100/100 frames; no frame pair crossing needed',
            preprocessing='Subtract per-block temporal mean and each frame spatial mean; Hann time/space; remove temporal DC and neighbour bins',
            spectrum='sqrt-energy weighted mean of each spatial-frequency normalized power in a 3-bin tube at f_t=-v*f_x modulo 1',
            block_gates='dynamic std>=1; abs(v)>=1; not boundary; score>=0.10; score>=max(19 shuffles)+0.03; competitor margin>=0.03 outside max(3,15% speed); reversal<=0.25',
            run_gates='At least 4/5 SIMULTANEOUS blocks passing transport, both alternate-frequency bands, both axial halves, and saturation<=20% in ALL THREE regions. Original separate marginal 4/5 gates retained as sensitivity. Full all-five-block mean supported only if all 5 simultaneous blocks pass.',
            aggregation='Candidate coordinates average all five blocks, including failures. Separate accepted_block coordinates use only simultaneously supported blocks for all three regions; >=4 required for pair support. Full candidate mean requires all5. Upper/lower wall speeds averaged without equality assumption; disagreement is reported.',
            selection='No tuning, band selection, run-specific method switch, class labels or inter-run prediction used'),
        limitations=['All acceptance criteria are empirical transport diagnostics, not calibrated error or formal significance.',
            'The temporal blocks are not independent experimental repeats.',
            'Frequency-band/axial sensitivity measures different optical scales/content.',
            'A regional image wave is not an identified phase-material parcel.',
            'Nonstationary optical wave speed is represented by five block maxima; sub-block variability and weaker coexisting waves remain unresolved.',
            '19 temporal shuffles are negative controls, not a claimed inferential p-value.',
            'Native profile archive is reused with hash; no new experimental data or spatial calibration.'])
    (BASE/'results/spectral_metadata.json').write_text(json.dumps(serial(metadata),indent=2))
    print(json.dumps(serial(metadata),indent=2),flush=True)

if __name__=='__main__':main()
