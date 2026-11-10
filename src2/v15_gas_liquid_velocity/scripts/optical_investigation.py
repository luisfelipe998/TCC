#!/usr/bin/env python3
"""Native-frame displacement investigation. No conversion to phase velocity.

Run from anywhere with python3 -B optical_investigation.py. Outputs are all
under v15_gas_liquid_velocity and do not alter source acquisitions.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import re
from collections import Counter
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter1d, shift as image_shift
from scipy.signal import correlate

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT.parent / 'data_processed'
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / 'intermediate' / 'optical_mplcache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROIS = ((10, 70), (20, 60), (25, 55))
LAGS = (1, 2, 4)
METHODS = ('ncc_raw', 'ncc_dynamic', 'ncc_dynamic_highpass', 'phase_dynamic_highpass')
MAX_SHIFT = 240
SAMPLE_SIZE = 100
PATTERNS = ('dispersed bubbles', 'plug', 'slug', 'churn', 'wavy annular')
COLORS = dict(zip(PATTERNS, ('#b58900', '#4269d0', '#d55e00', '#009e73', '#9c6ade')))


def save_csv(name, rows):
    path = ROOT / 'results' / name
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def peak_details(values, shifts):
    """Quadratic subpixel interpolation; competitor outside +/-5 pixels."""
    if not np.isfinite(values).any():
        return dict(shift=np.nan, peak=np.nan, margin=np.nan, boundary=True)
    k = int(np.nanargmax(values))
    offset = 0.0
    if 0 < k < len(values) - 1:
        den = values[k-1] - 2*values[k] + values[k+1]
        if np.isfinite(den) and den < -1e-12:
            offset = float(np.clip(0.5*(values[k-1]-values[k+1])/den, -0.5, 0.5))
    other = values[np.abs(shifts-shifts[k]) > 5]
    return dict(shift=float(shifts[k]+offset), peak=float(values[k]),
                margin=float(values[k]-np.nanmax(other)),
                boundary=bool(abs(shifts[k]) >= MAX_SHIFT-2))


def ncc(a, b):
    """Fixed-template, signed NCC: positive shift means A advects right in B.

    The central len(a)-2*MAX_SHIFT samples of A are compared against equally
    sized B windows. Each candidate has its own mean/variance normalization.
    Equal template length prevents correlation from rising at tiny overlaps.
    """
    t = a[MAX_SHIFT:-MAX_SHIFT].astype(float)
    t = t-t.mean()
    n = len(t)
    b = np.asarray(b, float)
    c = correlate(b, t, mode='valid', method='fft')
    sums = np.r_[0., np.cumsum(b)]
    sums2 = np.r_[0., np.cumsum(b*b)]
    s = sums[n:]-sums[:-n]
    var = np.maximum(sums2[n:]-sums2[:-n]-s*s/n, 0.)
    den = np.sqrt(np.sum(t*t)*var)
    values = np.full_like(c, np.nan)
    np.divide(c, den, out=values, where=den > 1e-9)
    return peak_details(values, np.arange(-MAX_SHIFT, MAX_SHIFT+1))


def phase(a, b):
    """Zero-padded phase-only correlation with a Hann window."""
    n = len(a)
    fft_n = 2**int(np.ceil(np.log2(2*n-1)))
    win = np.hanning(n)
    fa = np.fft.rfft((a-np.mean(a))*win, n=fft_n)
    fb = np.fft.rfft((b-np.mean(b))*win, n=fft_n)
    cross = fb*np.conj(fa)
    absolute = np.abs(cross)
    normalized = np.zeros_like(cross)
    np.divide(cross, absolute, out=normalized, where=absolute > 1e-8*absolute.max())
    c = np.fft.irfft(normalized, n=fft_n)
    shifts = np.arange(-MAX_SHIFT, MAX_SHIFT+1)
    return peak_details(c[shifts % fft_n], shifts)


def preprocess(profiles):
    dynamic = profiles - np.mean(profiles, axis=0, keepdims=True)
    return dynamic-gaussian_filter1d(dynamic, 8., axis=-1, mode='reflect')


def quality(result, a, b, method):
    texture = min(float(np.std(a)), float(np.std(b)))
    # Heuristic evidence thresholds fixed before inspecting the acquisitions.
    threshold = 0.60 if method.startswith('ncc') else 0.08
    margin = 0.03 if method.startswith('ncc') else 0.02
    passes = bool(np.isfinite(result['shift']) and result['peak'] >= threshold
                  and result['margin'] >= margin and not result['boundary']
                  and texture >= 1.)
    return passes, texture


def quantify(profiles, method, lag, indices, fps):
    fn = phase if method.startswith('phase') else ncc
    rows = []
    for i in indices:
        a, b = profiles[i], profiles[i+lag]
        r = fn(a, b)
        passed, texture = quality(r, a, b, method)
        rows.append(dict(frame_index=int(i), lag_frames=lag,
                         shift_px=r['shift'], displacement_px_per_frame=r['shift']/lag,
                         structure_celerity_px_per_s=r['shift']*fps/lag,
                         peak=r['peak'], competitor_margin=r['margin'],
                         boundary=r['boundary'], min_texture_std=texture,
                         heuristic_pass=passed))
    return rows


def summarize(rows):
    shifts = np.array([r['displacement_px_per_frame'] for r in rows])
    passed = np.array([r['heuristic_pass'] for r in rows])
    good = shifts[passed]
    return dict(pairs=len(rows), finite_pairs=int(np.isfinite(shifts).sum()),
                heuristic_pass_pairs=int(passed.sum()), heuristic_pass_fraction=float(passed.mean()),
                median_all_px_per_frame=float(np.nanmedian(shifts)),
                median_pass_px_per_frame=float(np.median(good)) if len(good) else np.nan,
                q25_pass_px_per_frame=float(np.quantile(good, .25)) if len(good) else np.nan,
                q75_pass_px_per_frame=float(np.quantile(good, .75)) if len(good) else np.nan,
                negative_fraction_all=float(np.mean(shifts < -0.5)),
                near_zero_fraction_all=float(np.mean(np.abs(shifts) < .5)),
                boundary_fraction=float(np.mean([r['boundary'] for r in rows])),
                low_texture_fraction=float(np.mean([r['min_texture_std'] < 1 for r in rows])))


def synthetic_tests():
    rng = np.random.default_rng(20261023)
    rows = []
    for seed in range(5):
        base = gaussian_filter1d(rng.normal(size=1016), 3.)
        base = (base-base.mean())/base.std()*25
        for truth in (-170., -31.4, -7., 0., 1., 7.35, 31., 93.6, 170.):
            for noise in (0., 2., 10.):
                a = base+rng.normal(0, noise, len(base))
                b = image_shift(base, truth, order=3, mode='constant')+rng.normal(0, noise, len(base))
                for name, fn in [('ncc', ncc), ('phase', phase)]:
                    r = fn(a, b)
                    rows.append(dict(test='translated_random_texture', seed=seed, method=name,
                                     true_shift_px=truth, noise_std=noise,
                                     estimated_shift_px=r['shift'], abs_error_px=abs(r['shift']-truth)))
    # Both estimators must fail to infer motion in a completely featureless pair.
    nulls = {}
    for name, fn in [('ncc', ncc), ('phase', phase)]:
        r = fn(np.ones(1016), np.ones(1016))
        ok, texture = quality(r, np.ones(1016), np.ones(1016), name)
        assert not ok
        nulls[name] = dict(heuristic_pass=ok, texture_std=texture)
    save_csv('optical_synthetic_trials.csv', rows)
    summary = []
    for method in ('ncc', 'phase'):
        for noise in (0., 2., 10.):
            e = np.array([r['abs_error_px'] for r in rows if r['method']==method and r['noise_std']==noise])
            summary.append(dict(method=method, noise_std=noise, n=len(e), mae_px=float(e.mean()),
                                median_error_px=float(np.median(e)), p95_error_px=float(np.quantile(e,.95)),
                                max_error_px=float(e.max()), fraction_within_1px=float(np.mean(e<=1))))
    # Sign and subpixel arithmetic validation, independent of data results.
    assert next(s for s in summary if s['method']=='ncc' and s['noise_std']==0)['max_error_px'] < .1
    return summary, nulls


def run():
    for sub in ('results', 'findings', 'figures', 'intermediate'):
        (ROOT/sub).mkdir(parents=True, exist_ok=True)
    synthetic, nulls = synthetic_tests()
    per_run, all_pairs, preferred_pairs, sensitivity, excluded, full_methods = [], [], [], [], [], []
    sample_agreement = []
    for folder in sorted(DATA.glob('m??q??_???')):
        markers = list(folder.glob('__*.txt'))
        if len(markers) != 1:
            raise ValueError(f'Expected exactly one marker: {folder}')
        pattern = markers[0].stem[2:].replace('_', ' ')
        if folder.name=='m05q05_002' or pattern not in PATTERNS:
            excluded.append(dict(run=folder.name, pattern=pattern))
            continue
        config_path = folder/'_config.mii'
        config = config_path.read_text()
        fps = float(re.search(r'Film speed 1\s*:([0-9.]+)', config).group(1))
        expected = int(re.search(r'Number of images 1\s*:([0-9]+)', config).group(1))
        files = sorted(folder.glob('ImgA*.jpg'))
        ids = [int(re.search(r'(\d+)\.jpg$', p.name).group(1)) for p in files]
        assert len(files)==expected and np.all(np.diff(ids)==1), folder.name
        digest = hashlib.sha256(config_path.read_bytes())
        profiles = [[] for _ in ROIS]
        saturated = []
        for file in files:
            raw = file.read_bytes()
            digest.update(file.name.encode())
            digest.update(hashlib.sha256(raw).digest())
            with Image.open(io.BytesIO(raw)) as im:
                a = np.asarray(im.convert('L'))
            assert a.shape == (80, 1016), (folder.name, a.shape)
            for k, (lo, hi) in enumerate(ROIS):
                profiles[k].append(a[lo:hi].mean(axis=0))
            saturated.append(float(np.mean(a[20:60]>=250)))
        profiles = np.asarray(profiles)
        dynamic = profiles - np.mean(profiles, axis=1, keepdims=True)
        processed = np.asarray([preprocess(p) for p in profiles])
        sample = np.linspace(0, len(files)-1-max(LAGS), SAMPLE_SIZE, dtype=int)
        full = quantify(processed[1], 'ncc_dynamic_highpass', 1, range(len(files)-1), fps)
        full_summary = summarize(full)
        for method in METHODS:
            p = profiles[1] if method=='ncc_raw' else dynamic[1] if method=='ncc_dynamic' else processed[1]
            frows = full if method=='ncc_dynamic_highpass' else quantify(p, method, 1, range(len(files)-1), fps)
            full_methods.append(dict(run=folder.name, pattern=pattern, method=method, **summarize(frows)))
            if method=='ncc_dynamic':
                preferred_pairs.extend([dict(run=folder.name, pattern=pattern, **r) for r in frows])
        all_pairs.extend([dict(run=folder.name, pattern=pattern, **r) for r in full])
        native_saturation = float(np.mean(saturated))
        summary = dict(run=folder.name, pattern=pattern, fps=fps, frames=len(files),
                       duration_between_first_last_s=(len(files)-1)/fps,
                       mean_native_core_saturation=native_saturation,
                       high_saturation=native_saturation>.20, input_sha256=digest.hexdigest(),
                       **full_summary)
        summary['median_pass_structure_celerity_px_per_s'] = summary['median_pass_px_per_frame']*fps
        estimates = {}
        for k, (lo,hi) in enumerate(ROIS):
            for lag in LAGS:
                for method in METHODS:
                    p = profiles[k] if method=='ncc_raw' else dynamic[k] if method=='ncc_dynamic' else processed[k]
                    measured = quantify(p, method, lag, sample, fps)
                    sensitivity.append(dict(run=folder.name, pattern=pattern, method=method,
                                            roi=f'{lo}:{hi}', lag_frames=lag, **summarize(measured)))
                    estimates[k,lag,method] = measured
        for method in ('ncc_raw', 'ncc_dynamic', 'phase_dynamic_highpass'):
            reference = estimates[1,1,'ncc_dynamic_highpass']
            alternative = estimates[1,1,method]
            both = [(a,b) for a,b in zip(reference, alternative) if a['heuristic_pass'] and b['heuristic_pass']]
            err = np.array([abs(a['displacement_px_per_frame']-b['displacement_px_per_frame']) for a,b in both])
            sample_agreement.append(dict(run=folder.name, pattern=pattern, alternative=method,
                                         sampled_pairs=SAMPLE_SIZE, jointly_passing_pairs=len(both),
                                         median_abs_difference_px_per_frame=float(np.median(err)) if len(err) else np.nan,
                                         fraction_within_2px=float(np.mean(err<=2)) if len(err) else np.nan))
        variants = [s for s in sensitivity if s['run']==folder.name and s['method']=='ncc_dynamic_highpass'
                    and s['heuristic_pass_pairs']>=20]
        medians = [s['median_pass_px_per_frame'] for s in variants]
        summary['sensitivity_supported_variants'] = len(variants)
        summary['sensitivity_median_min_px_per_frame'] = min(medians) if medians else np.nan
        summary['sensitivity_median_max_px_per_frame'] = max(medians) if medians else np.nan
        summary['sensitivity_median_range_px_per_frame'] = max(medians)-min(medians) if medians else np.nan
        summary['low_coverage'] = full_summary['heuristic_pass_fraction'] < .5
        per_run.append(summary)
        print(folder.name, pattern, f"pass={full_summary['heuristic_pass_fraction']:.3f}",
              f"median={full_summary['median_pass_px_per_frame']:.3f} px/frame", flush=True)
    assert len(per_run)==26 and len(excluded)==5
    save_csv('optical_per_run.csv', per_run)
    save_csv('optical_primary_pairs.csv', all_pairs)
    save_csv('optical_preferred_candidate_pairs.csv', preferred_pairs)
    save_csv('optical_sensitivity.csv', sensitivity)
    save_csv('optical_method_agreement.csv', sample_agreement)
    save_csv('optical_full_method_summary.csv', full_methods)
    preferred = []
    for r in per_run:
        candidate = next(s for s in full_methods if s['run']==r['run'] and s['method']=='ncc_dynamic')
        preferred.append(dict(run=r['run'], pattern=r['pattern'], fps=r['fps'],
                              **{k:v for k,v in candidate.items() if k not in ('run','pattern')},
                              median_pass_structure_celerity_px_per_s=candidate['median_pass_px_per_frame']*r['fps'],
                              high_saturation=r['high_saturation'],
                              low_coverage=candidate['heuristic_pass_fraction']<.5,
                              note='Selected for numerical coverage; physical accuracy unvalidated.'))
    save_csv('optical_preferred_candidate_per_run.csv', preferred)
    metrics = dict(selected_runs=len(per_run), patterns=dict(Counter(r['pattern'] for r in per_run)),
                   excluded=excluded, primary_pairs=len(all_pairs), preferred_candidate_pairs=len(preferred_pairs),
                   primary_heuristic_pass_pairs=sum(r['heuristic_pass_pairs'] for r in per_run),
                   full_method_quality={m:dict(passing_pairs=sum(r['heuristic_pass_pairs'] for r in full_methods if r['method']==m),
                                               runs_with_at_least_half_passing=sum(r['heuristic_pass_fraction']>=.5 for r in full_methods if r['method']==m)) for m in METHODS},
                   low_coverage_runs=[r['run'] for r in per_run if r['low_coverage']],
                   high_saturation_runs=[r['run'] for r in per_run if r['high_saturation']],
                   synthetic=synthetic, featureless_pair_checks=nulls,
                   fixed_method_settings=dict(max_displacement_px=MAX_SHIFT, rois=ROIS, lags=LAGS,
                                              sensitivity_pairs_per_variant=SAMPLE_SIZE, methods=METHODS,
                                              primary='ncc_dynamic_highpass; ROI 20:60; lag 1',
                                              highpass_gaussian_sigma_px=8,
                                              ncc_peak_threshold=.60, ncc_competitor_margin=.03,
                                              phase_peak_threshold=.08, phase_competitor_margin=.02,
                                              min_profile_std_graylevels=1., saturation_intensity=250,
                                              max_native_saturation_fraction=.20),
                   caveat='Optical structure celerity only. Heuristic quality is not calibrated measurement confidence. No physical velocity ground truth or verified axial scale.')
    (ROOT/'results'/'optical_metadata.json').write_text(json.dumps(metrics, indent=2, allow_nan=False)+'\n')
    figures(per_run, sensitivity)
    visual_correspondence()
    report(per_run, sensitivity, metrics, sample_agreement)
    verify_saved_outputs()


def figures(rows, sensitivity):
    fig, axes = plt.subplots(1, 2, figsize=(12, 8), gridspec_kw={'width_ratios':[3,1]}, sharey=True)
    for k,r in enumerate(rows):
        color = COLORS[r['pattern']]
        x = r['median_pass_px_per_frame']
        if np.isfinite(x):
            axes[0].plot([r['q25_pass_px_per_frame'],r['q75_pass_px_per_frame']],[k,k],color=color, lw=3)
            axes[0].plot(x,k,'o',color=color,markerfacecolor='white' if r['low_coverage'] or r['high_saturation'] else color)
        axes[1].barh(k,r['heuristic_pass_fraction'],color=color,alpha=.8)
    axes[0].set_yticks(range(len(rows)),[r['run'] for r in rows],fontsize=8)
    axes[0].set_xlabel('Signed optical displacement (px/frame); median and IQR of passing pairs')
    axes[1].set_xlabel('Passing fraction')
    axes[1].axvline(.5,color='black',ls=':',lw=1)
    axes[0].axvline(0,color='black',ls=':',lw=1)
    axes[0].invert_yaxis()
    fig.suptitle('High-pass NCC diagnostic: no calibrated phase velocity\nHollow: <50% passing pairs or >20% native-core saturation', fontsize=12)
    fig.tight_layout()
    fig.savefig(ROOT/'figures'/'optical_displacement_quality.png',dpi=180)
    plt.close(fig)

    fig, ax=plt.subplots(figsize=(8,5))
    for method,marker in [('ncc_raw','s'),('ncc_dynamic','x'),('phase_dynamic_highpass','^')]:
        for r in rows:
            base=next(s for s in sensitivity if s['run']==r['run'] and s['roi']=='20:60' and s['lag_frames']==1 and s['method']=='ncc_dynamic_highpass')
            alt=next(s for s in sensitivity if s['run']==r['run'] and s['roi']=='20:60' and s['lag_frames']==1 and s['method']==method)
            if base['heuristic_pass_pairs']>=20 and alt['heuristic_pass_pairs']>=20:
                ax.plot(base['median_pass_px_per_frame'],alt['median_pass_px_per_frame'],marker,color=COLORS[r['pattern']],alpha=.8)
    low,high=ax.get_xlim()
    ax.plot([low,high],[low,high],':',color='black')
    ax.set(xlabel='Dynamic high-pass NCC median (px/frame)',ylabel='Alternative median (px/frame)',
           title='Method agreement is not ground-truth accuracy\nSquares: raw NCC; crosses: dynamic NCC; triangles: phase; ≥20 passing pairs each')
    fig.tight_layout()
    fig.savefig(ROOT/'figures'/'optical_method_comparison.png',dpi=180)
    plt.close(fig)


def visual_correspondence():
    """Show native images and the actual fixed-template correspondence."""
    folder = DATA / 'm05q01_001'
    files = sorted(folder.glob('ImgA*.jpg'))
    frames = []
    for k in (200, 210):
        with Image.open(files[k]) as im:
            frames.append(np.asarray(im.convert('L')))
    a,b = [f[20:60].mean(axis=0) for f in frames]
    result = ncc(a,b)
    shift = result['shift']
    fig, axes = plt.subplots(3,1,figsize=(12,6), gridspec_kw={'height_ratios':[1,1,2]})
    for ax,frame,k in zip(axes[:2],frames,(200,210)):
        ax.imshow(frame,cmap='gray',vmin=0,vmax=255,aspect='auto')
        ax.axhline(20,color='#4269d0',lw=.8)
        ax.axhline(59,color='#4269d0',lw=.8)
        ax.set(ylim=(79,0),ylabel=f'Frame {k}',xlim=(0,1015))
    axes[0].axvspan(MAX_SHIFT,1016-MAX_SHIFT,color='#4269d0',alpha=.08)
    axes[1].axvspan(MAX_SHIFT+shift,1016-MAX_SHIFT+shift,color='#d55e00',alpha=.08)
    x=np.arange(1016)
    axes[2].plot(x,a,label='Frame 200 profile',color='#4269d0',lw=1)
    axes[2].plot(x,b,label='Frame 210, native coordinates',color='#999999',alpha=.7,lw=1)
    axes[2].plot(x-shift,b,label='Frame 210, shifted back for comparison',color='#d55e00',lw=1)
    axes[2].axvspan(MAX_SHIFT,1016-MAX_SHIFT,color='#4269d0',alpha=.08)
    axes[2].set(xlabel='Native axial pixel coordinate',ylabel='Mean intensity, rows 20:60',xlim=(0,1015))
    axes[2].legend(fontsize=8,loc='lower right')
    fig.suptitle(f'm05q01_001: native correspondence over 10 frames; signed raw NCC shift {shift:.2f} px\nNegative shift = leftward movement. Correspondence sanity check, not independent ground truth.',fontsize=12)
    fig.tight_layout()
    fig.savefig(ROOT/'figures'/'optical_native_correspondence.png',dpi=180)
    plt.close(fig)
    (ROOT/'results'/'optical_visual_correspondence.json').write_text(json.dumps(dict(
        run=folder.name, first_file=files[200].name, second_file=files[210].name,lag_frames=10,
        roi='20:60',method='ncc_raw',shift_px=shift,peak=result['peak'],margin=result['margin'],
        passes_adjacent_pair_gate=quality(result,a,b,'ncc_raw')[0],
        caveat='Image/profile correspondence check only, not independent displacement ground truth.'),indent=2)+'\n')


def verify_saved_outputs():
    """Independently recompute saved arithmetic and one brute-force NCC peak."""
    def read(name):
        with (ROOT/'results'/name).open() as file:
            return list(csv.DictReader(file))
    summaries=read('optical_per_run.csv')
    pairs=read('optical_primary_pairs.csv')
    sensitivity=read('optical_sensitivity.csv')
    preferred=read('optical_preferred_candidate_per_run.csv')
    preferred_pairs=read('optical_preferred_candidate_pairs.csv')
    assert len(pairs)==13000 and len(preferred_pairs)==13000 and len(sensitivity)==26*36 and len(preferred)==26
    for r in summaries:
        run_pairs=[s for s in pairs if s['run']==r['run']]
        passed=[s for s in run_pairs if s['heuristic_pass']=='True']
        assert len(passed)==int(r['heuristic_pass_pairs'])
        if passed:
            median=float(np.median([float(s['shift_px']) for s in passed]))
            assert abs(median-float(r['median_pass_px_per_frame']))<1e-12
        for s in run_pairs:
            assert abs(float(s['structure_celerity_px_per_s'])-float(s['shift_px'])*float(r['fps']))<1e-7
    preferred_pass_total=0
    for r in preferred:
        run_pairs=[s for s in preferred_pairs if s['run']==r['run']]
        assert len(run_pairs)==int(r['pairs'])==500
        assert sorted(int(s['frame_index']) for s in run_pairs)==list(range(500))
        passed=[]
        for s in run_pairs:
            lag=int(s['lag_frames'])
            displacement=float(s['shift_px'])
            assert lag==1 and s['pattern']==r['pattern']
            assert abs(float(s['displacement_px_per_frame'])-displacement/lag)<1e-12
            assert abs(float(s['structure_celerity_px_per_s'])-displacement*float(r['fps'])/lag)<1e-7
            numerical_gate=(np.isfinite(displacement) and float(s['peak'])>=.60
                            and float(s['competitor_margin'])>=.03
                            and s['boundary']=='False' and float(s['min_texture_std'])>=1.)
            assert numerical_gate==(s['heuristic_pass']=='True')
            if numerical_gate:
                passed.append(float(s['displacement_px_per_frame']))
        preferred_pass_total+=len(passed)
        assert len(passed)==int(r['heuristic_pass_pairs'])
        assert abs(len(passed)/len(run_pairs)-float(r['heuristic_pass_fraction']))<1e-12
        assert (len(passed)/len(run_pairs)<.5)==(r['low_coverage']=='True')
        if passed:
            for q,key in ((.25,'q25_pass_px_per_frame'),(.5,'median_pass_px_per_frame'),(.75,'q75_pass_px_per_frame')):
                assert abs(float(np.quantile(passed,q))-float(r[key]))<1e-12
            assert abs(float(np.median(passed))*float(r['fps'])-float(r['median_pass_structure_celerity_px_per_s']))<1e-7
    metadata=json.loads((ROOT/'results'/'optical_metadata.json').read_text())
    assert preferred_pass_total==metadata['full_method_quality']['ncc_dynamic']['passing_pairs']
    accepted=[r for r in preferred if r['high_saturation']=='False' and r['low_coverage']=='False']
    assert len(accepted)==17 and all(r['pattern']!='wavy annular' for r in accepted)
    rng=np.random.default_rng(42)
    a=rng.normal(size=1016)
    b=np.roll(a,37)+.5*rng.normal(size=1016)
    result=ncc(a,b)
    template=a[240:-240]
    coefficients=np.array([np.corrcoef(template,b[i:i+len(template)])[0,1] for i in range(481)])
    assert np.argmax(coefficients)-240==37 and abs(result['peak']-coefficients.max())<1e-12
    out=dict(passed=True,primary_pair_arithmetic_verified=len(pairs),
             preferred_candidate_pair_arithmetic_verified=len(preferred_pairs),
             preferred_candidate_passing_pairs_recomputed=preferred_pass_total,
             preferred_candidate_run_medians_and_iqr_recomputed=len(preferred),
             preferred_candidate_quality_flags_and_coverage_recomputed=True,
             sensitivity_rows_verified=len(sensitivity),preferred_runs_passing_coverage_and_saturation=len(accepted),
             independent_bruteforce_ncc_peak_matches=True,independent_bruteforce_test_signed_integer_shift_px=37)
    (ROOT/'results'/'optical_validation.json').write_text(json.dumps(out,indent=2)+'\n')


def report(rows, sensitivity, metrics, agreement):
    lines=['# Native-frame optical displacement investigation','',
           '**These are structure/interface displacements, not independently measured liquid or gas phase velocities.** No verified axial pixel scale was found. Multiplying by the per-acquisition camera rate gives px/s, not m/s.','',
           f"The analysis retains {len(rows)} acquisitions, excludes the five requested runs/categories before processing, and uses all {metrics['primary_pairs']:,} adjacent native-frame pairs for the primary estimate. No stitched panorama is used. The supplied labels are unchanged.",'',
           '## Methods and fixed checks','',
           '- Primary: signed, normalized correlation of a fixed central template with equal-length candidate windows. Positive is motion to the right from frame A to B. Search ±240 px; no forced positive shift. Parabolic peak interpolation supplies a subpixel estimate.','- Profile preprocessing: average rows 20:60, subtract the per-pixel temporal mean over the acquisition, then subtract a Gaussian blur (σ=8 px). This suppresses stationary illumination and broad intensity gradients, but may also remove persistent structures.','- Alternatives: raw-profile NCC, temporal-mean-subtracted NCC without high-pass filtering, and Hann-windowed, zero-padded phase-only correlation on the dynamic high-pass profile. Sensitivity uses rows 10:70, 20:60 and 25:55; temporal lags 1, 2 and 4; 100 evenly spaced starting frames per acquisition and variant. At lag 4 the same ±240 px search permits at most ±60 px/frame: boundary hits flag this restricted range.','- A passing NCC pair has peak ≥0.60, a ≥0.03 margin over the strongest competitor outside ±5 px, profile standard deviation ≥1 gray level, and a peak more than 2 px from the search boundary. Phase-only thresholds are peak ≥0.08 and margin ≥0.02 with the same texture/boundary rules. These fixed heuristic rules are **not** probability/confidence estimates and have different acceptance rates across methods.','- Native-core saturation means gray level ≥250 in rows 20:60, averaged over all frames; >20% is flagged. A correlation can pass despite saturated images. Runs with <50% passing primary pairs are flagged, retained, and not silently given representative velocity status.','',
           '## Reference high-pass results','',
           f"The reference high-pass NCC diagnostic passes on {metrics['primary_heuristic_pass_pairs']:,}/{metrics['primary_pairs']:,} pairs ({metrics['primary_heuristic_pass_pairs']/metrics['primary_pairs']:.1%}). {len(metrics['low_coverage_runs'])}/{len(rows)} acquisitions have <50% passing pairs; {len(metrics['high_saturation_runs'])} have >20% native-core saturation. These rates assess numerical support, not physical accuracy.",'',
           '![Displacement and numerical support](../figures/optical_displacement_quality.png)','',
           '| Acquisition | Pattern | Passing pairs | Median passing px/frame | IQR px/frame | Median passing px/s | Flags |','|---|---|---:|---:|---:|---:|---|']
    for r in rows:
        flags=[]
        if r['low_coverage']: flags.append('low coverage')
        if r['high_saturation']: flags.append('saturation')
        lines.append(f"| {r['run']} | {r['pattern']} | {r['heuristic_pass_pairs']}/{r['pairs']} | {r['median_pass_px_per_frame']:.3f} | {r['q25_pass_px_per_frame']:.3f}–{r['q75_pass_px_per_frame']:.3f} | {r['median_pass_structure_celerity_px_per_s']:.1f} | {', '.join(flags) or '—'} |")
    lines += ['', '| Method (all adjacent pairs) | Passing pairs / 13,000 | Runs with ≥50% passing pairs |', '|---|---:|---:|']
    for method, entry in metrics['full_method_quality'].items():
        lines.append(f"| {method} | {entry['passing_pairs']} | {entry['runs_with_at_least_half_passing']} / 26 |")
    lines += ['', 'Different passing rates do not establish relative accuracy: high-pass filtering discards broad profiles, while raw profiles can follow illumination and large interfaces. Phase and NCC thresholds are not on the same scale. All signed medians are given in the saved full-method table.']
    lines += ['', 'IQR is observed temporal spread among accepted pairs, not an uncertainty interval for the true motion. Temporal-background subtraction, the texture threshold and rejection can preferentially select short interface passages, especially in bubbly or intermittent flows. A median of accepted displacements is not a volume- or mass-weighted phase velocity.', '',
              '## Synthetic estimator verification','',
              'Five independent Gaussian-smoothed random profiles, nine signed translations (including fractions of a pixel and ±170 px), and additive independent noise levels of 0, 2 and 10 gray levels were tested. Profile signal standard deviation is 25 gray levels. Ground truth is known **only in these constructed translations**; deformation, occlusion and multiphase topology are absent. Featureless pairs fail the quality gate for both methods.', '',
              '| Estimator | Noise σ | n | MAE px | 95th percentile absolute error px | Within 1 px |','|---|---:|---:|---:|---:|---:|']
    for s in metrics['synthetic']:
        lines.append(f"| {s['method']} | {s['noise_std']:.0f} | {s['n']} | {s['mae_px']:.4f} | {s['p95_error_px']:.4f} | {s['fraction_within_1px']:.1%} |")
    lines += ['', 'These tests verify implementation/sign and illustrate behavior under translation. They cannot rank physical velocity accuracy on the acquisitions.', '',
              '## Method and crop/lag sensitivity','',
              '![Cross-method medians](../figures/optical_method_comparison.png)','',
              '| Alternative to dynamic high-pass NCC | Jointly passing sampled pairs | Median of per-run median absolute differences (px/frame) |','|---|---:|---:|']
    for method in ('ncc_raw','ncc_dynamic','phase_dynamic_highpass'):
        subset=[r for r in agreement if r['alternative']==method]
        values=[r['median_abs_difference_px_per_frame'] for r in subset if r['jointly_passing_pairs']>0]
        lines.append(f"| {method} | {sum(r['jointly_passing_pairs'] for r in subset)} | {np.median(values):.3f} |")
    lines += ['', 'Agreement uses the intersection of each method’s passing pairs. This selects easier pairs and must not be read as an unbiased accuracy estimate. Each row in [optical_method_agreement.csv](../results/optical_method_agreement.csv) preserves run-specific support and differences. [optical_sensitivity.csv](../results/optical_sensitivity.csv) preserves all 36 method/crop/lag variants per run. [optical_full_method_summary.csv](../results/optical_full_method_summary.csv) also compares all four methods on every adjacent frame pair.', '',
              '## Preferred optical candidate and limits','',
              'Prefer temporal-mean-subtracted NCC without high-pass filtering as the candidate for further optical validation: it gives the broadest numerical support under the common NCC gate (7,528/13,000 pairs; 18/26 runs with at least half passing, including saturated m03q40_001). Only 17/26 runs meet both coverage and the saturation rule; none is wavy annular. This selection is based on coverage, not proven physical accuracy. High-pass filtering reduces support to 3,928 pairs and 6/26 runs; retain that baseline as a sensitivity diagnostic rather than declaring its more selective estimates more accurate. Raw NCC and phase-only correlation expose preprocessing dependence. The [preferred-candidate per-run CSV](../results/optical_preferred_candidate_per_run.csv) contains its signed medians, frame rates, IQRs and quality flags. All 13,000 underlying adjacent-pair records, including rejected pairs, are retained in [optical_preferred_candidate_pairs.csv](../results/optical_preferred_candidate_pairs.csv). No method is declared the most accurate real-data estimator because there is no displacement or phase-velocity ground truth.', '',
              '![Native frames and aligned intensity profiles](../figures/optical_native_correspondence.png)', '',
              'The native-frame comparison uses m05q01_001 frames 200 and 210. Their broad intensity structure moves to smaller axial pixel coordinates; the translated profiles overlay closely. The exact signed shift is recorded in [optical_visual_correspondence.json](../results/optical_visual_correspondence.json). This confirms correspondence and sign only; the same images generate the NCC estimate, so the overlay is not independent ground truth. This 10-frame illustration has NCC=0.959 but competitor margin=0.0228, below the adjacent-pair margin threshold 0.03; it is not presented as passing that gate.', '',
              'A tracked shadow, bubble nose, interface wave, or shared brightness pattern travels at a structure celerity. It need not equal either phase-average velocity, and a single axial brightness profile cannot separate simultaneous liquid and vapor velocity fields. Liquid-phase velocity would require tracers/phase-resolved velocimetry or independently measured superficial flux and cross-sectional void fraction. A literature camera pixel scale is not proof that this dataset has the same magnification; the native metadata must be calibrated independently. The use of an apparent tube diameter would also require distortion and axial/radial scale validation.', '',
              'Existing pipeline_v11.py uses the reverse correlation orientation relative to its positive-right docstring and clamps shifts to ≥1 pixel. Therefore saved panorama strip widths are unsuitable as an unbiased signed velocity record; this analysis recomputes displacements directly from native images.', '',
              'Reproduce: `python3 -B src2/v15_gas_liquid_velocity/scripts/optical_investigation.py`. Input fingerprints, actual FPS values, sequence counts, fixed thresholds and synthetic results are saved in [metadata](../results/optical_metadata.json) and [per-run results](../results/optical_per_run.csv). Hashes cover config bytes, sorted filenames and native-frame content hashes.']
    lines += ['', '[Saved-result verification](../results/optical_validation.json) independently recomputes all 13,000 reference and 13,000 preferred-candidate pixel-rate conversions, reaggregates preferred-candidate quality gates, coverage, medians and IQRs from the saved pairs, checks the 936 sensitivity rows and 17 unsaturated preferred candidates, and compares an FFT NCC result with brute-force NumPy correlation on a known signed shift.']
    (ROOT/'findings'/'optical_investigation.md').write_text('\n'.join(lines)+'\n')


if __name__ == '__main__':
    run()
