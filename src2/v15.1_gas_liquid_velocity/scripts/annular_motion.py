#!/usr/bin/env python3
"""Spatial image-motion alternatives, with special annular diagnostics.

No phase-material velocity assignment or spatial calibration is made here.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import re
import time
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter, gaussian_filter1d
from scipy.fft import rfft, irfft

BASE = Path(__file__).resolve().parents[1]
DATA = BASE.parent / 'data_processed'
V15 = BASE.parent / 'v15_gas_liquid_velocity'
os.environ.setdefault('MPLCONFIGDIR', str(BASE/'intermediate/annular_mplcache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LAGS = (1, 2, 4)
MAX_SHIFT = 240
WIDTH = 1016 - 2*MAX_SHIFT
SHIFTS = np.arange(-MAX_SHIFT, MAX_SHIFT+1)
ROIS = {'core': (35,54), 'upper_inner': (24,35), 'lower_inner': (54,65),
        'upper_outer': (18,29), 'lower_outer': (60,71)}
NULL_REPEATS = 19
SEED = 20261023


def save_csv(name, rows):
    with (BASE/'results'/name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def to_json(value):
    if isinstance(value, dict): return {k:to_json(v) for k,v in value.items()}
    if isinstance(value, (list, tuple)): return [to_json(v) for v in value]
    if isinstance(value, np.ndarray): return to_json(value.tolist())
    if isinstance(value, (np.bool_,)): return bool(value)
    if isinstance(value, (np.integer,)): return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    return value


def correlations(a, b, source_x=MAX_SHIFT):
    """One fixed-template normalized correlation curve per pair."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    t = a[:, source_x:source_x+WIDTH]
    t = t - t.mean(axis=1, keepdims=True)
    c = irfft(rfft(b, n=2048, axis=1)*rfft(t, n=2048, axis=1).conj(), n=2048, axis=1)[:, :481]
    sums = np.pad(np.cumsum(b, axis=1), ((0,0),(1,0)))
    sums2 = np.pad(np.cumsum(b*b, axis=1), ((0,0),(1,0)))
    s = sums[:, WIDTH:] - sums[:, :-WIDTH]
    var = np.maximum(sums2[:, WIDTH:] - sums2[:, :-WIDTH] - s*s/WIDTH, 0)
    energy = np.sum(t*t, axis=1)
    den = np.sqrt(energy[:, None]*var)
    corr = np.full_like(c, np.nan)
    np.divide(c, den, out=corr, where=den>1e-9)
    texture = np.minimum(np.std(a, axis=1), np.std(b, axis=1))
    return corr, texture


def peak(curve, shifts=SHIFTS):
    if not np.isfinite(curve).any():
        return dict(shift=np.nan, peak=np.nan, margin=np.nan, boundary=True)
    k = int(np.nanargmax(curve))
    off = 0.
    if 0<k<len(curve)-1:
        den = curve[k-1]-2*curve[k]+curve[k+1]
        if np.isfinite(den) and den < -1e-12:
            off = float(np.clip(.5*(curve[k-1]-curve[k+1])/den,-.5,.5))
    return dict(shift=float(shifts[k]+off), peak=float(curve[k]),
                margin=float(curve[k]-np.nanmax(curve[np.abs(shifts-shifts[k])>5])),
                boundary=bool(k<=2 or k>=len(curve)-3))


def mean_curve(corr, valid):
    if not np.any(valid): return np.full(481, np.nan)
    finite = np.isfinite(corr[valid])
    return np.divide(np.nansum(corr[valid],axis=0), finite.sum(axis=0),
                     out=np.full(481,np.nan), where=finite.sum(axis=0)>0)


def temporal_analysis(profiles, texture_min, saturation, seed=SEED):
    p = profiles - profiles.mean(axis=0, keepdims=True)
    stats, curves, pair_arrays = {}, {}, {}
    for lag in LAGS:
        corr, texture = correlations(p[:-lag], p[lag:])
        valid = texture>=texture_min
        detections = [peak(curve) for curve in corr]
        shifts = np.array([d['shift'] for d in detections])
        accepted = valid & np.array([d['peak']>=.6 and d['margin']>=.03 and not d['boundary'] for d in detections])
        curve = mean_curve(corr, valid)
        ensemble = peak(curve)
        stats[lag] = dict(lag=lag, comparisons=len(corr), valid_texture_pairs=int(valid.sum()),
            individual_accepted_pairs=int(accepted.sum()), individual_coverage=float(accepted.mean()),
            individual_median_px_frame=float(np.median(shifts[accepted])/lag) if accepted.any() else np.nan,
            ensemble_shift_px=ensemble['shift'], ensemble_velocity_px_frame=ensemble['shift']/lag,
            ensemble_peak=ensemble['peak'], ensemble_margin=ensemble['margin'], boundary=ensemble['boundary'])
        curves[lag] = curve
        pair_arrays[lag] = np.column_stack([shifts, np.array([d['peak'] for d in detections]),
                                           np.array([d['margin'] for d in detections]), texture, accepted])
    # Full temporal order is deliberately destroyed; all 501 frames remain used.
    rng = np.random.default_rng(seed)
    null_peaks_by_lag={lag:[] for lag in LAGS}
    null_speeds_by_lag={lag:[] for lag in LAGS}
    for _ in range(NULL_REPEATS):
        indices = rng.permutation(len(p))
        for lag in LAGS:
            corr, texture = correlations(p[indices[:-lag]],p[indices[lag:]])
            result=peak(mean_curve(corr,texture>=texture_min))
            null_peaks_by_lag[lag].append(result['peak'])
            null_speeds_by_lag[lag].append(result['shift']/lag)
    null_peaks=null_peaks_by_lag[1];null_speeds=null_speeds_by_lag[1]
    finite_null = np.asarray(null_peaks)[np.isfinite(null_peaks)]
    null_max=float(finite_null.max()) if len(finite_null) else np.nan
    for lag in LAGS:
        finite=np.asarray(null_peaks_by_lag[lag])[np.isfinite(null_peaks_by_lag[lag])]
        maximum=float(finite.max()) if len(finite) else np.nan
        stats[lag]['shuffle_max_peak']=maximum
        stats[lag]['peak_and_shuffle_supported']=bool(np.isfinite(maximum) and
            stats[lag]['ensemble_peak']>=.10 and stats[lag]['ensemble_peak']>=maximum+.03 and
            stats[lag]['valid_texture_pairs']>=250 and not stats[lag]['boundary'])
    block_speeds, block_peaks = [], []
    # Five disjoint temporal windows; their comparisons do not cross boundaries.
    for indices in np.array_split(np.arange(len(p)),5):
        block=p[indices]
        corr, texture=correlations(block[:-1],block[1:])
        result=peak(mean_curve(corr,texture>=texture_min))
        block_speeds.append(result['shift']);block_peaks.append(result['peak'])
    speed=stats[1]['ensemble_velocity_px_frame']
    tolerance=max(2.,.2*abs(speed)) if np.isfinite(speed) else np.nan
    lag_consistent=all(np.isfinite(stats[k]['ensemble_velocity_px_frame']) and
                       abs(stats[k]['ensemble_velocity_px_frame']-speed)<=tolerance and
                       not stats[k]['boundary'] for k in (2,4))
    block_agree=sum(np.isfinite(v) and abs(v-speed)<=tolerance for v in block_speeds)
    reasons=[]
    if stats[1]['valid_texture_pairs']<250: reasons.append('low_texture_coverage')
    if not np.isfinite(speed) or abs(speed)<1: reasons.append('missing_or_near_zero_motion')
    if stats[1]['boundary']: reasons.append('search_boundary')
    if not np.isfinite(stats[1]['ensemble_peak']) or stats[1]['ensemble_peak']<.10: reasons.append('weak_ensemble_peak')
    if not np.isfinite(null_max) or stats[1]['ensemble_peak']<null_max+.03: reasons.append('weak_shuffle_separation')
    if not lag_consistent: reasons.append('lag_inconsistency')
    if block_agree<4: reasons.append('block_inconsistency')
    if saturation>.20: reasons.append('region_saturation')
    summary=dict(ensemble_supported=not reasons, withheld_reasons=';'.join(reasons),
        signed_velocity_px_frame=speed, ensemble_peak=stats[1]['ensemble_peak'],
        individual_accepted_pairs=stats[1]['individual_accepted_pairs'],
        individual_coverage=stats[1]['individual_coverage'],
        individual_median_px_frame=stats[1]['individual_median_px_frame'],
        valid_texture_pairs=stats[1]['valid_texture_pairs'], saturation_fraction=saturation,
        shuffle_max_peak=null_max, shuffle_median_peak=float(np.median(finite_null)) if len(finite_null) else np.nan,
        lag2_velocity_px_frame=stats[2]['ensemble_velocity_px_frame'],
        lag4_velocity_px_frame=stats[4]['ensemble_velocity_px_frame'],
        lag_consistent=lag_consistent, block_agree_count=block_agree,
        min_block_px_frame=float(np.nanmin(block_speeds)) if np.isfinite(block_speeds).any() else np.nan,
        max_block_px_frame=float(np.nanmax(block_speeds)) if np.isfinite(block_speeds).any() else np.nan,
        temporal_tolerance_px_frame=tolerance)
    # Fixed-window reversal and alternate axial templates are sensitivity checks.
    # They do not follow an identical material parcel or turn optical labels into phases.
    sensitivity=[]
    alternatives=[('time_reversed',p[::-1],240),('source_x160',p,160),('source_x320',p,320),('raw_profiles',profiles,240),
        ('highpass_sigma8',p-gaussian_filter1d(p,8.,axis=1),240),
        ('highpass_sigma32',p-gaussian_filter1d(p,32.,axis=1),240)]
    for label,sequence,start in alternatives:
        corr,texture=correlations(sequence[:-1],sequence[1:],source_x=start)
        shifts=np.arange(481)-start
        result=peak(mean_curve(corr,texture>=texture_min),shifts)
        difference=abs(result['shift']+speed) if label=='time_reversed' else abs(result['shift']-speed)
        sensitivity.append(dict(check=label,source_x=start,source_width=WIDTH,
            min_shift_px=int(shifts[0]),max_shift_px=int(shifts[-1]),signed_px_frame=result['shift'],
            peak=result['peak'],boundary=result['boundary'],difference_px_frame=difference,
            within_tolerance=bool(np.isfinite(difference) and difference<=tolerance and not result['boundary'])))
    summary['time_reversal_consistent']=sensitivity[0]['within_tolerance']
    summary['axial_template_consistent']=all(r['within_tolerance'] for r in sensitivity[1:3])
    summary['all_lag_peaks_supported']=all(stats[lag]['peak_and_shuffle_supported'] for lag in LAGS)
    summary['robust_ensemble_supported']=bool(summary['ensemble_supported'] and summary['time_reversal_consistent'] and summary['axial_template_consistent'] and summary['all_lag_peaks_supported'])
    return summary, stats, curves, pair_arrays, dict(null_peaks=null_peaks,null_speeds=null_speeds,
        null_peaks_by_lag=null_peaks_by_lag,null_speeds_by_lag=null_speeds_by_lag,
        block_speeds=block_speeds,block_peaks=block_peaks,sensitivity=sensitivity)


def edge_profiles(frames, lo, hi, polarity):
    smooth=gaussian_filter(frames.astype(np.float32),sigma=(0,.8,1.2))
    gradient=np.diff(smooth,axis=1)*polarity
    local=gradient[:,lo:hi]
    at=np.argmax(local,axis=1)
    height=lo+at.astype(float)+.5
    strength=np.take_along_axis(local,at[:,None,:],axis=1)[:,0]
    boundary=(at==0)|(at==hi-lo-1)
    return height, strength, boundary


def synthetic_controls():
    rng=np.random.default_rng(SEED)
    tests=[]
    for truth in (-35.,-17.3,23.75):
        texture=gaussian_filter1d(rng.normal(size=12000),2.)*80
        x=np.arange(1016)
        seq=np.stack([np.interp(x-truth*t+4000,np.arange(len(texture)),texture) for t in range(61)])
        seq+=10*np.sin(x/75)[None,:]
        dynamic=seq-seq.mean(axis=0)
        estimates=[]
        for lag in LAGS:
            c,t=correlations(dynamic[:-lag],dynamic[lag:])
            result=peak(mean_curve(c,t>=1.))
            error=abs(result['shift']/lag-truth)
            assert error<.15,(truth,lag,error)
            tests.append(dict(test='known_translation',truth_px_frame=truth,lag=lag,
                              estimate_px_frame=result['shift']/lag,error_px_frame=error))
            estimates.append(result['shift']/lag)
    constant=np.ones((51,1016))
    c,t=correlations(constant,constant)
    assert np.isnan(c).all() and not (t>=1).any()
    fixed=np.tile(rng.normal(size=1016),(51,1))
    fixed-=fixed.mean(axis=0)
    c,t=correlations(fixed[:-1],fixed[1:])
    assert not(t>=1).any()
    # Correlation is checked against direct Pearson at a chosen displacement.
    a=rng.normal(size=(2,1016));b=rng.normal(size=(2,1016))
    c,_=correlations(a,b)
    direct=np.corrcoef(a[0,240:776],b[0,253:789])[0,1]
    assert abs(c[0,253]-direct)<1e-12
    return dict(known_translations=tests,constant_rejected=True,static_background_rejected=True,
                direct_pearson_difference=abs(c[0,253]-direct))


def run(pilot=False):
    start=time.monotonic();controls=synthetic_controls()
    scope=list(csv.DictReader((V15/'results/scenario_inventory.csv').open()))
    scope=[r for r in scope if r['included']=='True']
    if pilot: scope=[r for r in scope if r['run'] in ('m03q10_001','m06q80_001')]
    summaries,lag_rows,inputs,blocks,null_trials,additional=[],[],[],[],[],[]
    arrays={};profiles_archive={}
    for entry in scope:
        folder=DATA/entry['run'];files=sorted(folder.glob('ImgA*.jpg'))
        ids=[int(re.search(r'([0-9]+)[.]jpg$',f.name)[1]) for f in files]
        assert len(files)==501 and np.all(np.diff(ids)==1)
        xsv=(folder/'_config.xsv').read_text();mii=(folder/'_config.mii').read_text()
        fps=float(re.search(r'^Rate=([0-9.]+)',xsv,re.M)[1])
        assert fps==float(re.search(r'Film speed 1[ ]*:([0-9.]+)',mii)[1])==float(entry['native_fps'])
        frames=np.stack([np.asarray(Image.open(f).convert('L')) for f in files])
        assert frames.shape==(501,80,1016)
        digest=hashlib.sha256()
        for f in files:digest.update(hashlib.sha256(f.read_bytes()).digest())
        inputs.append(dict(run=entry['run'],pattern=entry['pattern'],fps=fps,frames=501,
            duration_s=500/fps,ordered_image_sha256=digest.hexdigest(),
            config_sha256=hashlib.sha256((folder/'_config.xsv').read_bytes()).hexdigest(),
            mii_sha256=hashlib.sha256((folder/'_config.mii').read_bytes()).hexdigest()))
        variants={name:(frames[:,lo:hi].mean(axis=1),1.,float(np.mean(frames[:,lo:hi]>=250)),np.nan,np.nan)
                  for name,(lo,hi) in ROIS.items()}
        if entry['pattern']=='wavy annular':
            for name,lo,hi,polarity in [('upper_edge',22,42,1),('lower_edge',45,67,-1)]:
                profile,strength,boundary=edge_profiles(frames,lo,hi,polarity)
                variants[name]=(profile,.25,float(np.mean(frames[:,lo:hi+1]>=250)),
                    float(np.mean(strength<2.)),float(np.mean(boundary)))
        for name,(profile,texture,saturation,weak,boundary) in variants.items():
            summary,bylag,curves,pairs,diagnostics=temporal_analysis(profile,texture,saturation)
            row=dict(run=entry['run'],pattern=entry['pattern'],region=name,fps=fps,
                profile_units='image_row_px' if 'edge' in name else 'gray_level',
                weak_edge_fraction=weak,edge_search_boundary_fraction=boundary,**summary)
            row['signed_velocity_px_s']=row['signed_velocity_px_frame']*fps
            if 'edge' in name and (weak>.2 or boundary>.2):
                row['ensemble_supported']=False
                row['robust_ensemble_supported']=False
                row['withheld_reasons']+=';ambiguous_edge_location'
            summaries.append(row)
            for lag,values in bylag.items():
                lag_rows.append(dict(run=entry['run'],region=name,**values))
                arrays[f"{entry['run']}__{name}__lag{lag}"]=pairs[lag].astype(np.float32)
                arrays[f"{entry['run']}__{name}__curve{lag}"]=curves[lag].astype(np.float32)
            for k,(v,p) in enumerate(zip(diagnostics['block_speeds'],diagnostics['block_peaks'])):
                blocks.append(dict(run=entry['run'],region=name,block_index=k,signed_px_frame=v,peak=p))
            for lag in LAGS:
                for k,(v,p) in enumerate(zip(diagnostics['null_speeds_by_lag'][lag],diagnostics['null_peaks_by_lag'][lag])):
                    null_trials.append(dict(run=entry['run'],region=name,lag=lag,shuffle_index=k,signed_px_frame=v,peak=p))
            for check in diagnostics['sensitivity']:
                additional.append(dict(run=entry['run'],region=name,**check))
            profiles_archive[f"{entry['run']}__{name}"]=profile.astype(np.float32)
        print(entry['run'],[(r['region'],round(r['signed_velocity_px_frame'],2),r['ensemble_supported'])
                           for r in summaries if r['run']==entry['run']],flush=True)
    prefix='annular_pilot_' if pilot else 'annular_'
    save_csv(prefix+'region_summary.csv',summaries);save_csv(prefix+'lag_summary.csv',lag_rows)
    save_csv(prefix+'block_summary.csv',blocks)
    save_csv(prefix+'null_trials.csv',null_trials)
    save_csv(prefix+'additional_checks.csv',additional)
    np.savez_compressed(BASE/'results'/(prefix+'pair_results.npz'),**arrays)
    np.savez_compressed(BASE/'intermediate'/(prefix+'profiles.npz'),**profiles_archive)
    paired=[]
    for entry in scope:
        regs={r['region']:r for r in summaries if r['run']==entry['run']}
        core,upper,lower=(regs[k] for k in ('core','upper_inner','lower_inner'))
        nearwall=(upper['signed_velocity_px_frame']+lower['signed_velocity_px_frame'])/2
        same_sign=np.sign(upper['signed_velocity_px_frame'])==np.sign(lower['signed_velocity_px_frame'])
        agree=bool(same_sign and abs(upper['signed_velocity_px_frame']-lower['signed_velocity_px_frame'])<=max(2.,.25*abs(nearwall)))
        reasons=[k+':'+regs[k]['withheld_reasons'] for k in ('core','upper_inner','lower_inner') if not regs[k]['ensemble_supported']]
        if not agree:reasons.append('upper_lower_disagreement')
        if not np.isfinite(nearwall) or abs(nearwall)<1:reasons.append('nearwall_missing_or_near_zero')
        pair=dict(run=entry['run'],pattern=entry['pattern'],fps=float(entry['native_fps']),
            eligible_spatial_descriptor=not reasons,withheld_reasons='|'.join(reasons),
            core_signed_px_frame=core['signed_velocity_px_frame'],nearwall_signed_px_frame=nearwall,
            upper_signed_px_frame=upper['signed_velocity_px_frame'],lower_signed_px_frame=lower['signed_velocity_px_frame'],
            upper_lower_agreement=agree,core_supported=core['ensemble_supported'],
            upper_supported=upper['ensemble_supported'],lower_supported=lower['ensemble_supported'],
            robust_spatial_descriptor=bool(not reasons and all(r['robust_ensemble_supported'] for r in (core,upper,lower))),
            core_signed_px_s=core['signed_velocity_px_s'],nearwall_signed_px_s=nearwall*float(entry['native_fps']),
            interpretation='Spatial image-structure descriptor; no gas/liquid material assignment')
        paired.append(pair)
    save_csv(prefix+'spatial_pairs.csv',paired)
    metadata=dict(status='pilot' if pilot else 'complete',script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        runtime_s=time.monotonic()-start,inputs=inputs,synthetic_controls=controls,roi_half_open=ROIS,
        annular_edge_search=dict(upper=[22,42],lower=[45,67],smoothing_sigma_y=.8,smoothing_sigma_x=1.2),
        correlation=dict(source_x=[240,776],search_px=[-240,240],lags=LAGS,
            preprocessing='Subtract 501-frame mean at each axial profile coordinate',
            pair_gates=dict(peak=.6,competitor_margin=.03,competitor_exclusion_px=5,texture_gray=1.,texture_edge_px=.25)),
        ensemble_gates=dict(min_texture_pairs=250,min_peak=.10,shuffle_repeats=NULL_REPEATS,
            min_peak_above_max_shuffle=.03,lag_tolerance='max(2 px/frame,20% of lag1 speed)',
            min_agreeing_blocks=4,blocks=5,max_region_saturation=.20,min_absolute_speed_px_frame=1.,
            max_edge_boundary_or_weak_fraction=.20,upper_lower_tolerance='max(2 px/frame,25% of mean wall speed)'),
        additional_checks=dict(time_reversal='Compare fixed-template reversed sequence with negative reference shift',
            axial_source_x=[160,240,320],source_width=536,
            axial_search='Full valid target windows; shift limits are [-source_x,480-source_x]',
            tolerance='max(2 px/frame,20% of reference speed)',
            robust_flag='Original ensemble gates plus both alternate templates and reversal within tolerance, and peak>=0.10 and peak>=shuffle_max+0.03 at each lag1/2/4',
            raw_profiles='Diagnostic of temporal-background subtraction; no raw-based eligibility gate',
            highpass='Additional lag1 sensitivity at Gaussian sigma8 and sigma32; changes the tracked optical scales, not a reporting gate'),
        limitations=['The gates are feasibility diagnostics, not accuracy levels or significance tests.',
            'Core and near-wall labels are spatial, not material phase labels.',
            'Upper/lower edge coordinates track selected brightness gradients, not validated film thickness.',
            'Ensemble peaks can favor the dominant coherent optical structure and suppress minority motion.',
            'ROI, temporal averaging, search window and rapid deformation limit representativeness.',
            'The same interfaces can influence core and wall profiles; agreement is not phase independence.'],
        counts=dict(runs=len(scope),annular=sum(r['pattern']=='wavy annular' for r in scope),
            regions=len(summaries),supported_regions=sum(r['ensemble_supported'] for r in summaries),
            eligible_spatial_pairs=sum(r['eligible_spatial_descriptor'] for r in paired)))
    (BASE/'results'/(prefix+'metadata.json')).write_text(json.dumps(to_json(metadata),indent=2)+chr(10))
    print(json.dumps(metadata['counts']), 'seconds',metadata['runtime_s'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--pilot',action='store_true')
    run(parser.parse_args().pilot)
