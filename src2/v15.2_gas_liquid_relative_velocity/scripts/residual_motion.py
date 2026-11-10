#!/usr/bin/env python3
"""Held-out interface-distance subtraction and phase-candidate optical motion.

This script never identifies material phase velocities. Masks are optical
candidates; residual texture may remain refracted/shared illumination.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import importlib.util
import json
import time
import numpy as np
from PIL import Image
from scipy.ndimage import binary_erosion, distance_transform_edt, gaussian_filter1d

BASE = Path(__file__).resolve().parents[1]
SRC2 = BASE.parent
MASKS = SRC2 / 'v15_gas_liquid_velocity/phase_conditioned_motion/intermediate/masks'
OLD_SCRIPT = SRC2 / 'v15.1_gas_liquid_velocity/all_scenarios_adversarial/scripts/phase_local_experiment.py'
spec = importlib.util.spec_from_file_location('phase_local_math', OLD_SCRIPT)
math = importlib.util.module_from_spec(spec)
spec.loader.exec_module(math)
PHASES = ('gas_candidate', 'liquid_candidate')
METHODS = ('raw', 'fixed_background', 'interface_distance')
SEED = 20261024
LAGS = (1, 2, 4)

def write_csv(name, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with (BASE / 'results' / name).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)

def dump(name, obj):
    (BASE / 'results' / name).write_text(json.dumps(math.clean(obj), indent=2) + '\n')

def halves(n):
    return (np.arange(n // 2), np.arange(n // 2, n))

def block_indices(n):
    return [b for half in halves(n) for b in np.array_split(half, 3)]

def distance_bins(mask):
    """Signed distance to the segmented boundary, bins width 2, clipped +/-24.

    The tube-wall boundary can also contribute to this segmentation distance.
    Radial position is a separate predictor. No optical-flow estimate is used.
    """
    result = np.empty(mask.shape, np.int16)
    for i, m in enumerate(mask):
        signed = distance_transform_edt(m) - distance_transform_edt(~m)
        result[i] = np.rint(np.clip(signed, -24, 24) / 2).astype(np.int16) + 12
    return result

def crossfit_predictions(images, bins, min_count=100):
    """Each half is predicted with means learned only from the other half.

    Interface model: E[I | row, binned signed distance], smoothed in distance
    (sigma one bin). Sparse bins interpolate within their row, never from the
    evaluation images. The fixed model is the training-half per-pixel mean.
    """
    n, h, w = images.shape
    fixed = np.empty_like(images)
    interface = np.empty_like(images)
    diagnostics = []
    row_offsets = (np.arange(h) * 25)[None, :, None]
    cat = bins + row_offsets
    split = halves(n)
    for side in range(2):
        train, test = split[1-side], split[side]
        fixed[test] = images[train].mean(axis=0, dtype=np.float64).astype(np.float32)
        counts = np.bincount(cat[train].ravel(), minlength=h*25).reshape(h, 25)
        sums = np.bincount(cat[train].ravel(), weights=images[train].ravel(), minlength=h*25).reshape(h, 25)
        means = sums / np.maximum(counts, 1)
        populated = counts >= min_count
        for row in range(h):
            good = np.flatnonzero(populated[row])
            if not len(good):
                means[row] = images[train, row].mean()
            else:
                means[row] = np.interp(np.arange(25), good, means[row, good])
        means = gaussian_filter1d(means, 1., axis=-1, mode='nearest')
        interface[test] = means.ravel()[cat[test]]
        unsupported = ~populated.ravel()[cat[test]]
        diagnostics.append(dict(evaluation_half=side, training_first=int(train[0]),
            training_last=int(train[-1]), evaluation_first=int(test[0]), evaluation_last=int(test[-1]),
            supported_categories=int(populated.sum()), total_categories=h*25,
            extrapolated_pixel_fraction=float(unsupported.mean()),
            raw_std=float(images[test].std()), fixed_residual_std=float((images[test]-fixed[test]).std()),
            interface_residual_std=float((images[test]-interface[test]).std())))
    return {'raw':images, 'fixed_background':images-fixed, 'interface_distance':images-interface}, diagnostics

def analyze(images, mask, original, nulls=True):
    """Six independent temporal blocks; no requirement of equal block speeds.

    Optical support is a per-block gate. A run needs >=5/6 blocks, at least two
    in each half. Final pairs average only common accepted blocks. This is a
    selected-subsequence mean, not a full-video phase mean.
    """
    blocks=[]; lagrows=[]; nullrows=[]
    rng=np.random.default_rng(SEED)
    for bi, ids in enumerate(block_indices(len(images))):
        a=images[ids]; m=mask[ids]; n=len(ids)
        ch=math.channels(a, m)
        pp={}
        for lag in LAGS:
            sums=math.sufficient(ch,np.arange(n-lag),np.arange(lag,n))
            p=math.peak(sums);pp[lag]=p
            lagrows.append(dict(block=bi,lag=lag,comparisons=n-lag,**p))
        controls=[]
        if nulls:
            for repeat in range(3):
                perm=rng.permutation(n)
                p=math.peak(math.sufficient(ch,perm[:-1],perm[1:]))
                nullrows.append(dict(block=bi,repeat=repeat,**p))
                if np.isfinite(p['peak']):controls.append(p['peak'])
        null_peak=max(controls) if controls else np.nan
        del ch
        p=pp[1];v=p['shift'];tol=max(2.,.2*abs(v)) if np.isfinite(v) else np.nan
        phase_area=(m.sum(axis=(1,2))>=100)
        rawstd=float(original[ids][m].std()) if m.any() else 0.
        sat=float(np.mean(original[ids][m]>=250)) if m.any() else 0.
        reasons=[]
        if int(phase_area.sum()) < 20:reasons.append('insufficient_phase_area')
        if rawstd < 1:reasons.append('no_raw_phase_texture')
        if not np.isfinite(v):reasons.append('undefined_correlation')
        if not all(np.isfinite(pp[k]['peak']) and pp[k]['peak']>=.2 and pp[k]['std']>=1
            and not pp[k]['boundary'] and abs(pp[k]['shift']/k-v)<=tol for k in LAGS):
            reasons.append('lag_or_texture_inconsistency')
        if not np.isfinite(p['margin']) or p['margin']<.02:reasons.append('ambiguous_peak')
        if nulls and (not np.isfinite(null_peak) or p['peak'] < null_peak+.05):
            reasons.append('null_not_separated')
        if sat>.2:reasons.append('phase_saturation')
        blocks.append(dict(block=bi,half=bi//3,first_frame=int(ids[0]),last_frame=int(ids[-1]),
            frame_intervals=n-1,signed_px_frame=v,lag2_px_frame=pp[2]['shift']/2,lag4_px_frame=pp[4]['shift']/4,
            peak=p['peak'],margin=p['margin'],texture_std=p['std'],raw_phase_std=rawstd,
            phase_area_fraction=float(m.mean()),frames_with_100_pixels=int(phase_area.sum()),
            overlap_pixel_pairs=p['overlap'],null_max_peak=null_peak,saturation_fraction=sat,
            optical_supported=not reasons,reasons=';'.join(reasons)))
    accepted=[b for b in blocks if b['optical_supported']]
    support=len(accepted)>=5 and all(sum(b['half']==h for b in accepted)>=2 for h in range(2))
    velocity=float(np.average([b['signed_px_frame'] for b in accepted],weights=[b['frame_intervals'] for b in accepted])) if accepted else np.nan
    summary=dict(signed_px_frame=velocity,optical_supported=support,accepted_blocks=len(accepted),
        accepted_intervals=sum(b['frame_intervals'] for b in accepted),available_intervals=sum(b['frame_intervals'] for b in blocks),
        material_phase_identified=False,reason='' if support else 'fewer_than_5_of_6_blocks')
    return summary,blocks,lagrows,nullrows

def erode(mask, radius):
    return binary_erosion(mask, structure=np.ones((1,2*radius+1,2*radius+1),bool))

def synthetic():
    rng=np.random.default_rng(SEED)
    n,h,w=181,24,1016
    xx=np.arange(w); tt=np.arange(n)
    g=np.zeros((n,h,w),bool);g[:,:h//2]=True
    cases={}
    tex=gaussian_filter1d(rng.normal(size=(h,14000)),1.2,axis=-1)*45
    speeds=(-17.3,-7.4)
    independent=np.empty((n,h,w),np.float32)
    for row in range(h):
        v=speeds[int(row>=h//2)]
        independent[:,row]=np.stack([np.interp(xx-v*i+3000,np.arange(14000),tex[row])+120 for i in tt])
    cases['independent_material_textures']=(independent,g,speeds)
    shared=np.empty_like(independent)
    illumination=gaussian_filter1d(rng.normal(size=14000),1.2)*45
    for i in tt:shared[i]=120+np.interp(xx-12*i+4000,np.arange(14000),illumination)[None,:]
    cases['shared_moving_illumination']=(shared,g,(12.,12.))
    moving=np.stack([((xx-12*i)%500)<250 for i in tt])[:,None,:]
    moving=np.repeat(moving,h,axis=1)
    bins=distance_bins(moving)
    # Exactly one brightness-vs-distance law, no material texture.
    halo=(90+3*(bins-12)+30*np.exp(-((bins-12)/3.)**2)).astype(np.float32)
    cases['interface_halo_without_material_texture']=(halo,moving.copy(),(None,None))
    for i in tt:
        for row in range(h):
            moving[i,row]=((xx-12*i-12*np.sin(i*.14+row*.25))%500)<250
    bins=distance_bins(moving)
    deformed=(90+3*(bins-12)+30*np.exp(-((bins-12)/3.)**2)).astype(np.float32)
    cases['deforming_interface_without_material_texture']=(deformed,moving,(None,None))
    cases['blank_interiors']=(np.where(g,160,45).astype(np.float32),g,(None,None))
    rows=[]
    for name,(im,mask,truths) in cases.items():
        prepared,_=crossfit_predictions(im,distance_bins(mask))
        for method in METHODS:
            for pi,phase in enumerate(PHASES):
                m=erode(mask if pi==0 else ~mask,3)
                out,blocks,_,_=analyze(prepared[method],m,im,nulls=True)
                truth=truths[pi]
                error=abs(out['signed_px_frame']-truth) if truth is not None and np.isfinite(out['signed_px_frame']) else np.nan
                rows.append(dict(case=name,method=method,phase=phase,known_optical_speed=truth,
                    recovered_px_frame=out['signed_px_frame'],absolute_error=error,
                    supported=out['optical_supported'],accepted_blocks=out['accepted_blocks'],
                    physical_material_tracer=(name=='independent_material_textures')))
                if name=='independent_material_textures':
                    assert out['optical_supported'] and error<.15,(name,method,phase,out,error)
                if name=='blank_interiors':assert not out['optical_supported']
                if name=='shared_moving_illumination':
                    assert out['optical_supported'] and error<.15,(name,method,phase,out,error)
    return rows

def pair_rows(scope, summaries, blocks):
    output=[]
    for inv in scope:
        run=inv['run'];fps=float(inv['fps'])
        for radius in (3,1):
            for method in METHODS:
                b=[x for x in blocks if x['run']==run and x['erosion_radius']==radius and x['method']==method]
                pp={p:{x['block']:x for x in b if x['phase']==p} for p in PHASES}
                common=[i for i in range(6) if all(pp[p][i]['optical_supported'] for p in PHASES)]
                support=len(common)>=5 and all(sum(i//3==h for i in common)>=2 for h in range(2))
                weights=[pp[PHASES[0]][i]['frame_intervals'] for i in common]
                vv=[float(np.average([pp[p][i]['signed_px_frame'] for i in common],weights=weights)) if common else np.nan for p in PHASES]
                ratio=vv[0]/vv[1] if support and np.isfinite(vv[1]) and abs(vv[1])>2 else np.nan
                output.append(dict(run=run,pattern=inv['pattern'],fps=fps,method=method,erosion_radius=radius,
                    common_blocks=len(common),common_block_ids=';'.join(map(str,common)),
                    matched_intervals=sum(weights),pair_optical_supported=support,
                    gas_candidate_px_frame=vv[0],liquid_candidate_px_frame=vv[1],
                    gas_candidate_px_s=vv[0]*fps,liquid_candidate_px_s=vv[1]*fps,
                    signed_ratio_gas_candidate_over_liquid_candidate=ratio,material_pair_identified=False))
    return output

def run(pilot=False, synth_only=False):
    start=time.monotonic()
    synth=synthetic();write_csv('residual_synthetic.csv',synth)
    print('synthetic controls complete',flush=True)
    if synth_only:return
    with (SRC2/'v15.1_gas_liquid_velocity/results/scenario_inventory.csv').open() as f:
        scope=[r for r in csv.DictReader(f) if r['included']=='True']
    if pilot:scope=[r for r in scope if r['run'] in ('m05q01_001','m04q80_001','m07q01_001')]
    summaries=[];blocks=[];lags=[];nulls=[];provenance=[];predictions=[]
    for inv in scope:
        stamp=time.monotonic();run=inv['run'];fps=float(inv['fps'])
        files=sorted((SRC2/'data_processed'/run).glob('ImgA*.jpg'))
        assert len(files)==501 and all(f.stem==f'ImgA{i:06d}' for i,f in enumerate(files))
        images=np.stack([np.array(Image.open(f).convert('L'),dtype=np.float32) for f in files])
        assert images.shape==(501,80,1016)
        with np.load(MASKS/(run+'.npz'),allow_pickle=False) as z:
            assert float(z['fps'])==fps
            mm=[np.unpackbits(z[k],axis=-1,count=1016,bitorder='big').astype(bool) for k in ('masks_060','masks_070','masks_080')]
        gas=np.logical_and.reduce(mm); liquid=~np.logical_or.reduce(mm)
        liquid[:,:18]=False;liquid[:,71:]=False
        bins=distance_bins(mm[1])[:,18:71]
        del mm
        original=images[:,18:71]
        prepared,diag=crossfit_predictions(original,bins)
        for d in diag:predictions.append(dict(run=run,**d))
        for radius in (3,1):
            for phase,mask in zip(PHASES,(gas,liquid)):
                mask=erode(mask,radius)[:,18:71]
                for method in METHODS:
                    out,bs,ls,ns=analyze(prepared[method],mask,original)
                    keys=dict(run=run,pattern=inv['pattern'],fps=fps,erosion_radius=radius,phase=phase,method=method)
                    summaries.append(dict(**keys,signed_px_s=out['signed_px_frame']*fps,**out))
                    blocks.extend(dict(**keys,**x) for x in bs)
                    lags.extend(dict(**keys,**x) for x in ls)
                    nulls.extend(dict(**keys,**x) for x in ns)
        digest=hashlib.sha256()
        for f in files:digest.update(f.name.encode());digest.update(hashlib.sha256(f.read_bytes()).digest())
        provenance.append(dict(run=run,frames=len(files),ordered_named_image_sha256=digest.hexdigest(),
            mask_sha256=hashlib.sha256((MASKS/(run+'.npz')).read_bytes()).hexdigest()))
        print(run,round(time.monotonic()-stamp,1),'s',[(s['erosion_radius'],s['method'],s['phase'],s['accepted_blocks']) for s in summaries[-12:]],flush=True)
        for name,data in [('summary',summaries),('blocks',blocks),('lags',lags),('nulls',nulls),('inputs',provenance),('predictor',predictions)]:
            write_csv('residual_'+('pilot_' if pilot else '')+name+'.csv',data)
        del images,original,prepared,gas,liquid,bins
    pairs=pair_rows(scope,summaries,blocks)
    write_csv('residual_'+('pilot_' if pilot else '')+'pairs.csv',pairs)
    dump('residual_'+('pilot_' if pilot else '')+'metadata.json',dict(status='complete',runs=len(scope),
        runtime_s=time.monotonic()-start,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        reused_math_sha256=hashlib.sha256(OLD_SCRIPT.read_bytes()).hexdigest(),
        inventory_sha256=hashlib.sha256((SRC2/'v15.1_gas_liquid_velocity/results/scenario_inventory.csv').read_bytes()).hexdigest(),
        methods=METHODS,erosion_radii=[3,1],lags=LAGS,blocks=6,
        block_policy='Three blocks per half; >=5/6 supported, >=2 per half. Pair uses identical common blocks.',
        predictor='Other-half training only; radial row x signed segmentation-distance bin, 25 bins of width2 clipped +/-24px; sigma1-bin smooth; sparse bins (<100 pixels) interpolated using training only.',
        velocity_average='Mean of accepted per-block lag1 NCC velocities weighted by compared frame intervals; not peak of pooled correlations.',
        gates=dict(min_area_pixels_per_frame=100,min_area_frames_per_block=20,min_raw_std=1,min_residual_std=1,min_corr=.2,
            min_margin=.02,min_null_gap=.05,max_saturation=.2,min_overlap_pairs=1000,min_overlap_relative=.25,
            lag_tolerance='max(2 px/frame, 20% lag1 magnitude)',min_common_blocks=5),
        material_pairs_identified=0,limitations=[
            'The predictor removes one interface-conditioned brightness law, not all interface optics.',
            'No new material phase labels or tracer identity are created by residualization.',
            'Shared moving illumination passes the optical validation controls.',
            'Intervals across block boundaries are omitted; selected subsets can bias means.',
            'Other-half temporal training is separated, but halves are not independent experimental replicates.',
            'Static wall texture can dominate raw correlations; fixed subtraction can imprint training structures.',
            'Physical pixel calibration is absent; no certified m/s or cross-run magnification is assumed.']))
    print('DONE',time.monotonic()-start,flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--pilot',action='store_true');parser.add_argument('--synth-only',action='store_true')
    args=parser.parse_args();run(args.pilot,args.synth_only)
