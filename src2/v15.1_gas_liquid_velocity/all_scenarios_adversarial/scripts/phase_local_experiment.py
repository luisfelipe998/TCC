#!/usr/bin/env python3
"""Distributed phase-interior NCC; no contiguous rectangular patch required.

Correlations use grayscale at source AND target consensus-phase pixels.
The estimator measures phase-candidate optical motion, not material velocity.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import time
import numpy as np
from PIL import Image
from scipy.fft import rfft, irfft
from scipy.ndimage import binary_erosion, gaussian_filter1d

BASE = Path(__file__).resolve().parents[1]
V151 = BASE.parent
SRC2 = V151.parent
MASKS = SRC2/'v15_gas_liquid_velocity/phase_conditioned_motion/intermediate/masks'
LAGS = (1, 2, 4)
SHIFTS = np.arange(-240, 241)
NFFT = 2048
SEED = 20261024
PHASES = ('gas_candidate', 'liquid_candidate')

def clean(value):
    if isinstance(value, dict): return {k:clean(v) for k,v in value.items()}
    if isinstance(value, (list, tuple)): return [clean(v) for v in value]
    if isinstance(value, np.ndarray): return clean(value.tolist())
    if isinstance(value, (np.bool_,)): return bool(value)
    if isinstance(value, (np.integer,)): return int(value)
    if isinstance(value, (float, np.floating)): return float(value) if np.isfinite(value) else None
    return value

def write_csv(name, rows):
    with (BASE/'results'/name).open('w', newline='') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def channels(images, masks):
    m=masks.astype(np.float32)
    return [rfft(c, n=NFFT, axis=-1, workers=1) for c in (m,m*images,m*images*images)]

def sufficient(ch, ia, ib):
    # r_ba[d] = sum_x a[x] b[x+d], independent of FFT sign conventions.
    fields=((0,0),(1,0),(0,1),(2,0),(0,2),(1,1))
    outputs=[]
    for ca,cb in fields:
        spectrum=np.sum(ch[cb][ib]*np.conj(ch[ca][ia]),axis=(0,1),dtype=np.complex128)
        outputs.append(irfft(spectrum,n=NFFT)[SHIFTS%NFFT])
    return np.array(outputs)

def curve_from_sums(s):
    n,sa,sb,saa,sbb,sab=s
    n=np.maximum(n,0)
    safe=np.maximum(n,1)
    va=np.maximum(saa-sa*sa/safe,0); vb=np.maximum(sbb-sb*sb/safe,0)
    den=np.sqrt(va*vb)
    corr=np.full(len(n),np.nan)
    valid=(n>=max(1000,.25*float(n.max())))&(den>1e-6)
    corr[valid]=(sab[valid]-sa[valid]*sb[valid]/safe[valid])/den[valid]
    return corr,n,np.sqrt(va/safe),np.sqrt(vb/safe)

def peak(s):
    corr,n,sta,stb=curve_from_sums(s)
    if not np.isfinite(corr).any():
        return dict(shift=np.nan,peak=np.nan,margin=np.nan,overlap=0.,std=0.,boundary=True)
    k=int(np.nanargmax(corr)); offset=0.
    if 0<k<len(corr)-1 and np.isfinite(corr[k-1:k+2]).all():
        curvature=corr[k-1]-2*corr[k]+corr[k+1]
        if curvature< -1e-12: offset=float(np.clip(.5*(corr[k-1]-corr[k+1])/curvature,-.5,.5))
    competitors=corr[np.abs(SHIFTS-SHIFTS[k])>5]
    margin=float(corr[k]-np.nanmax(competitors)) if np.isfinite(competitors).any() else np.nan
    return dict(shift=float(SHIFTS[k]+offset),peak=float(corr[k]),margin=margin,
                overlap=float(n[k]),std=float(min(sta[k],stb[k])),boundary=bool(k<=2 or k>=len(corr)-3))

def analyze(images, mask, nulls=True):
    blocks=[]; null_sums=[]; sums_by_lag={k:[] for k in LAGS}
    dynamic=images-images.mean(axis=0,keepdims=True)
    raw_std=float(images[mask].std()) if mask.any() else 0.
    area_frames=int((mask.sum(axis=(1,2))>=100).sum())
    rng=np.random.default_rng(SEED)
    for bi, ids in enumerate(np.array_split(np.arange(len(images)),5)):
        ch=channels(dynamic[ids],mask[ids]); n=len(ids)
        for lag in LAGS:
            sums=sufficient(ch,np.arange(n-lag),np.arange(lag,n))
            sums_by_lag[lag].append(sums)
            p=peak(sums)
            blocks.append(dict(block=bi,lag=lag,comparisons=n-lag,
                first_frame=int(ids[0]),last_frame=int(ids[-1]),
                velocity_px_frame=p['shift']/lag,**p))
        if nulls:
            for repeat in range(3):
                perm=rng.permutation(n)
                null_sums.append((repeat,sufficient(ch,perm[:-1],perm[1:])))
        del ch
    stats={k:peak(np.sum(sums_by_lag[k],axis=0)) for k in LAGS}
    control=[peak(np.sum([s for r,s in null_sums if r==repeat],axis=0)) for repeat in range(3)] if nulls else []
    finite_null=[p['peak'] for p in control if np.isfinite(p['peak'])]
    null_max=max(finite_null) if finite_null else np.nan
    v=stats[1]['shift']; tol=max(2.,.2*abs(v)) if np.isfinite(v) else np.nan
    block_ok=sum(np.isfinite(b['peak']) and b['peak']>=.2 and b['std']>=1 and not b['boundary'] and
                 abs(b['velocity_px_frame']-v)<=tol for b in blocks if b['lag']==1)
    lag_ok=all(np.isfinite(stats[k]['peak']) and stats[k]['peak']>=.2 and stats[k]['std']>=1 and
               not stats[k]['boundary'] and abs(stats[k]['shift']/k-v)<=tol for k in LAGS)
    reasons=[]
    if area_frames<20: reasons.append('insufficient_phase_area')
    if raw_std<1: reasons.append('no_raw_phase_texture')
    if not np.isfinite(v): reasons.append('undefined_correlation')
    if not lag_ok: reasons.append('lag_or_texture_inconsistency')
    if block_ok<4: reasons.append('block_inconsistency')
    if not np.isfinite(stats[1]['margin']) or stats[1]['margin']<.02: reasons.append('ambiguous_peak')
    if nulls and (not np.isfinite(null_max) or stats[1]['peak']<null_max+.05): reasons.append('null_not_separated')
    summary=dict(signed_px_frame=v,lag2_px_frame=stats[2]['shift']/2,lag4_px_frame=stats[4]['shift']/4,
        peak=stats[1]['peak'],margin=stats[1]['margin'],dynamic_std=stats[1]['std'],raw_phase_std=raw_std,
        raw_phase_area_fraction=float(mask.mean()),frames_with_100_phase_pixels=area_frames,
        matched_pixel_pairs=stats[1]['overlap'],null_max_peak=null_max,lag_consistent=lag_ok,
        agreeing_blocks=block_ok,optical_supported=not reasons,reasons=';'.join(reasons))
    return summary,blocks,stats,control

def synth():
    rng=np.random.default_rng(SEED)
    n=61;h=20;w=1016;xx=np.arange(w);t=np.arange(n)
    gas=np.zeros((n,h,w),bool);gas[:,:10]=True;liq=~gas
    # Distinct grayscale textures translate independently in disjoint fixed rows.
    textures=[gaussian_filter1d(rng.normal(size=(h,8000)),1.1,axis=-1)*40 for _ in range(2)]
    results=[]
    for speeds in [(-17.3,-7.4),(13.25,-11.6)]:
        images=np.zeros((n,h,w),np.float32)
        for j,v in enumerate(speeds):
            for row in range(j*10,(j+1)*10):
                images[:,row]=np.stack([np.interp(xx-v*i+3000,np.arange(8000),textures[j][row]) for i in t])
        for j,mask in enumerate([gas,liq]):
            out,_,_,_=analyze(images,mask,nulls=False)
            err=abs(out['signed_px_frame']-speeds[j])
            assert err<.15 and out['optical_supported'],(speeds,j,out)
            results.append(dict(case='two_independent_textures',phase=PHASES[j],truth=speeds[j],estimate=out['signed_px_frame'],error=err))
    blank=images.copy();blank[:,10:]=37.
    out,_,_,_=analyze(blank,liq,nulls=False)
    assert not out['optical_supported'] and out['raw_phase_std']==0
    results.append(dict(case='untextured_liquid',phase=PHASES[1],truth=None,estimate=out['signed_px_frame'],error=None))
    static=np.repeat(images[:1],n,axis=0)
    out,_,_,_=analyze(static,gas,nulls=False)
    assert not out['optical_supported']
    # Pure moving binary interface has two constant interiors; masks alone
    # must not manufacture two identifiable material tracers.
    shift_mask=np.stack([((xx-12*i)%340)<170 for i in t])[:,None,:]
    shift_mask=np.repeat(shift_mask,h,axis=1)
    shared=np.where(shift_mask,170.,45.).astype(np.float32)
    for label,m in zip(PHASES,(shift_mask,~shift_mask)):
        out,_,_,_=analyze(shared,m,nulls=False)
        assert not out['optical_supported'] and out['raw_phase_std']==0
        results.append(dict(case='shared_binary_interface',phase=label,truth=None,estimate=out['signed_px_frame'],error=None))
    # Direct masked sufficient-statistic check at three shifts.
    a=rng.normal(size=(3,8,1016)).astype(np.float32)
    m=rng.random(a.shape)>.3;ch=channels(a,m);s=sufficient(ch,[0],[1]);max_error=0.
    for d in [-17,0,23]:
        xs=slice(max(0,-d),min(1016,1016-d));xt=slice(max(0,d),min(1016,1016+d))
        use=m[0,:,xs]&m[1,:,xt];av=a[0,:,xs][use];bv=a[1,:,xt][use]
        direct=np.array([len(av),av.sum(dtype=float),bv.sum(dtype=float),np.sum(av.astype(float)**2),np.sum(bv.astype(float)**2),np.sum(av.astype(float)*bv)])
        err=float(np.max(abs(s[:,d+240]-direct)));max_error=max(max_error,err)
        assert err<.002,(d,err)
    return dict(cases=results,static_rejected=True,blank_liquid_rejected=True,shared_interface_rejected=True,
        direct_sufficient_max_error=max_error,maximum_translation_error=max(r['error'] for r in results if r['error'] is not None))

def run(pilot=False):
    started=time.monotonic();validation=synth();print('synthetic passed',flush=True)
    with (V151/'results/scenario_inventory.csv').open() as f:scope=[r for r in csv.DictReader(f) if r['included']=='True']
    if pilot:scope=[r for r in scope if r['run'] in ('m05q01_001','m04q80_001','m07q01_001')]
    rows=[];blocks=[];lagrows=[];controls=[];provenance=[]
    for inv in scope:
        rs=time.monotonic();name=inv['run'];fps=float(inv['fps'])
        files=sorted((SRC2/'data_processed'/name).glob('ImgA*.jpg'))
        assert len(files)==501
        frames=np.stack([np.array(Image.open(f).convert('L'),dtype=np.float32) for f in files])
        assert frames.shape==(501,80,1016)
        with np.load(MASKS/(name+'.npz')) as z:
            assert float(z['fps'])==fps
            masks=[np.unpackbits(z[k],axis=-1,count=1016,bitorder='big').astype(bool) for k in ('masks_060','masks_070','masks_080')]
        g=np.logical_and.reduce(masks);l=~np.logical_or.reduce(masks);l[:,:18]=False;l[:,71:]=False
        del masks
        digest=hashlib.sha256()
        for f in files:digest.update(hashlib.sha256(f.read_bytes()).digest())
        provenance.append(dict(run=name,frames=len(files),ordered_image_sha256=digest.hexdigest(),mask_sha256=hashlib.sha256((MASKS/(name+'.npz')).read_bytes()).hexdigest()))
        for radius in (3,1):
            for phase,base_mask in zip(PHASES,(g,l)):
                mask=binary_erosion(base_mask,structure=np.ones((1,2*radius+1,2*radius+1),bool))[:,18:71]
                out,bs,st,co=analyze(frames[:,18:71],mask)
                sat=float(np.mean(frames[:,18:71][mask]>=250)) if mask.any() else 0.
                if sat>.20:
                    out['optical_supported']=False;out['reasons']+=';phase_saturation'
                row=dict(run=name,pattern=inv['pattern'],fps=fps,erosion_radius=radius,phase=phase,
                    signed_px_s=out['signed_px_frame']*fps,saturation_fraction=sat,**out)
                rows.append(row)
                for b in bs:blocks.append(dict(run=name,erosion_radius=radius,phase=phase,**b))
                for lag,p in st.items():lagrows.append(dict(run=name,erosion_radius=radius,phase=phase,lag=lag,**p))
                for repeat,p in enumerate(co):controls.append(dict(run=name,erosion_radius=radius,phase=phase,repeat=repeat,**p))
        print(name,'seconds',round(time.monotonic()-rs,1),[(r['erosion_radius'],r['phase'],round(r['signed_px_frame'],2),r['optical_supported']) for r in rows[-4:]],flush=True)
        del frames,g,l
    pairs=[]
    for inv in scope:
        rr=[r for r in rows if r['run']==inv['run']]
        for radius in (3,1):
            g=next(r for r in rr if r['phase']==PHASES[0] and r['erosion_radius']==radius)
            l=next(r for r in rr if r['phase']==PHASES[1] and r['erosion_radius']==radius)
            finite=np.isfinite(g['signed_px_s']) and np.isfinite(l['signed_px_s'])
            pairs.append(dict(run=inv['run'],pattern=inv['pattern'],erosion_radius=radius,
                gas_candidate_px_s=g['signed_px_s'],liquid_candidate_px_s=l['signed_px_s'],both_finite=finite,
                gas_optical_supported=g['optical_supported'],liquid_optical_supported=l['optical_supported'],
                pair_optical_supported=g['optical_supported'] and l['optical_supported'],
                gas_reasons=g['reasons'],liquid_reasons=l['reasons'],material_pair_identified=False))
    prefix='phase_local_pilot_' if pilot else 'phase_local_'
    for name,data in [('summary.csv',rows),('blocks.csv',blocks),('lags.csv',lagrows),('nulls.csv',controls),('pairs.csv',pairs),('inputs.csv',provenance)]:write_csv(prefix+name,data)
    metadata=dict(status='complete',runtime_s=time.monotonic()-started,runs=len(scope),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),synthetic=validation,
        method='Whole-lumen distributed phase-interior grayscale NCC using source/target mask intersections at every tested shift. Aggregate raw sufficient statistics across matched pixels, then normalize.',
        preprocessing='Subtract the full-record per-pixel temporal mean; require also nonzero raw within-phase intensity variance.',
        masks='Gas AND(060,070,080); liquid candidate lumen AND NOT OR(060,070,080); spatial erosion radius3 primary/radius1 sensitivity.',
        source_and_destination_masks=True,lags=LAGS,search_px=[-240,240],blocks=5,
        pairs_by_lag={k:sum(len(ids)-k for ids in np.array_split(np.arange(501),5)) for k in LAGS},
        null='Three independently shuffled orders in each temporal block, concatenated via sufficient statistics; feasibility control, no p-value.',
        gates=dict(min_phase_pixels=100,min_frames=20,min_raw_std=1.,min_dynamic_std=1.,min_corr=.2,min_margin=.02,
            min_overlap=max(1000,0),min_relative_overlap=.25,min_null_gap=.05,min_blocks=4,
            tolerance='max(2 px/frame,20% absolute lag1 speed)',max_saturation=.20),
        limitations=['Consensus silhouette masks do not establish material phase purity.',
            'Mask-conditioned NCC pools scattered pixels and can track shared optical structure.',
            'Full-record temporal background and five blocks are consistency diagnostics, not out-of-sample experimental validation.',
            'Raw pixel counts are not independent observations.',
            'Pure axial matching does not account for transverse motion or out-of-plane changes.',
            'No calibration or independent phase-velocity truth is supplied.'])
    (BASE/'results'/(prefix+'metadata.json')).write_text(json.dumps(clean(metadata),indent=2)+'\n')
    print('DONE',prefix,'runtime',metadata['runtime_s'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--pilot',action='store_true');a=p.parse_args();run(a.pilot)
