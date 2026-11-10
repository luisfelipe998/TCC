#!/usr/bin/env python3
"""Fixed 2-D patch NCC feasibility experiment; image motion, not phase speed.

Reads native frames only. All outputs stay within this evaluation. Thresholds and
windows are fixed globally, with no pattern-specific fitting.
"""
from __future__ import annotations
import csv
import hashlib
import json
import re
import time
from pathlib import Path
import numpy as np
from PIL import Image
from scipy.signal import correlate

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / 'results'
V15 = BASE.parent
DATA = V15.parent / 'data_processed'
PATTERNS = {'dispersed bubbles', 'plug', 'slug', 'churn', 'wavy annular'}
DX, DY = 96, 4
WIDTH, HEIGHT = 160, 16
WINDOWS = [(x, y) for x in (220, 636) for y in (12, 32, 52)]
METHODS = ('raw', 'background_subtracted')
N_PAIRS = 50
MIN_CORR, MIN_MARGIN, MIN_TEXTURE = .60, .03, 1.
ROUNDTRIP_LIMIT = 1.
SAT_LIMIT = .20


def box_sums(a, height, width):
    c = np.pad(np.cumsum(np.cumsum(a, axis=0), axis=1), ((1, 0), (1, 0)))
    return c[height:, width:] - c[:-height, width:] - c[height:, :-width] + c[:-height, :-width]


def match(a, b, x, y):
    """Template at (x,y) in A, search around same position in B."""
    t = np.asarray(a[y:y+HEIGHT, x:x+WIDTH], dtype=float)
    s = np.asarray(b[y-DY:y+HEIGHT+DY, x-DX:x+WIDTH+DX], dtype=float)
    assert t.shape == (HEIGHT, WIDTH) and s.shape == (HEIGHT+2*DY, WIDTH+2*DX)
    t = t - t.mean()
    std = float(t.std())
    numerator = correlate(s, t, mode='valid', method='fft')
    n = HEIGHT*WIDTH
    sums, sums2 = box_sums(s, HEIGHT, WIDTH), box_sums(s*s, HEIGHT, WIDTH)
    var = np.maximum(sums2-sums*sums/n, 0.)
    den = np.sqrt(np.sum(t*t)*var)
    corr = np.full_like(numerator, np.nan)
    np.divide(numerator, den, out=corr, where=den>1e-9)
    if not np.isfinite(corr).any():
        return dict(dx=np.nan, dy=np.nan, ix=0, iy=0, peak=np.nan, margin=np.nan,
                    texture=std, boundary=True, passed=False)
    iy, ix = np.unravel_index(np.nanargmax(corr), corr.shape)
    peak = float(corr[iy, ix])
    yy, xx = np.indices(corr.shape)
    competitor = corr[(np.abs(xx-ix)>5) | (np.abs(yy-iy)>2)]
    margin = float(peak-np.nanmax(competitor))
    offsets = []
    for axis, k in ((1, ix), (0, iy)):
        values = corr[iy, :] if axis == 1 else corr[:, ix]
        d = 0.
        if 0<k<len(values)-1:
            curvature = values[k-1]-2*values[k]+values[k+1]
            if np.isfinite(curvature) and curvature < -1e-12:
                d = float(np.clip(.5*(values[k-1]-values[k+1])/curvature, -.5, .5))
        offsets.append(d)
    boundary = bool(abs(ix-DX)>=DX-1 or abs(iy-DY)>=DY)
    texture = min(std, float(np.sqrt(var[iy, ix]/n)))
    return dict(dx=float(ix-DX+offsets[0]), dy=float(iy-DY+offsets[1]),
                ix=int(ix-DX), iy=int(iy-DY), peak=peak, margin=margin,
                texture=texture, boundary=boundary,
                passed=bool(peak>=MIN_CORR and margin>=MIN_MARGIN and
                            texture>=MIN_TEXTURE and not boundary))


def track(a, b, x, y):
    f = match(a, b, x, y)
    r = None
    rt = np.nan
    # The reverse template is sampled at the matched INTEGER location in B.
    # Its correct return displacement is therefore minus that integer shift.
    if f['passed']:
        r = match(b, a, x+f['ix'], y+f['iy'])
        rt = float(np.hypot(f['ix']+r['dx'], f['iy']+r['dy']))
    ok = bool(f['passed'] and r and r['passed'] and rt<=ROUNDTRIP_LIMIT)
    return f, r, rt, ok


def write_csv(name, rows):
    with (OUT/('local_'+name)).open('w', newline='') as h:
        w = csv.DictWriter(h, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)


def sanity():
    from scipy.ndimage import gaussian_filter, shift
    rng = np.random.default_rng(20261023)
    a = gaussian_filter(rng.normal(size=(80, 1016)), .8)
    a = a/a.std()*30
    tests = []
    for dx, dy in [(-45, 0), (13, 2), (0, 0), (31.25, -1.2), (-13.4, 1.3)]:
        b = shift(a, (dy, dx), order=3, mode='constant')
        f, r, rt, passed = track(a, b, 220, 32)
        err = float(np.hypot(f['dx']-dx, f['dy']-dy))
        tests.append(dict(dx_true=dx, dy_true=dy, dx=f['dx'], dy=f['dy'],
                          error_px=err, roundtrip_px=rt, passed=passed))
        assert err < .15 and passed, tests[-1]
    f, _, _, ok = track(np.ones_like(a), np.ones_like(a), 220, 32)
    assert not ok and not f['passed']
    # Independent direct Pearson correlation for one candidate window.
    b = shift(a, (2, 13), order=0, mode='constant')
    f = match(a, b, 220, 32)
    t = a[32:48, 220:380].ravel()
    s = b[34:50, 233:393].ravel()
    direct = float(np.corrcoef(t, s)[0, 1])
    assert abs(f['peak']-direct)<1e-10
    return dict(translation_tests=tests, constant_rejected=True,
                peak_direct_pearson_difference=abs(f['peak']-direct))


def run():
    start = time.monotonic()
    validation = sanity()
    old = {}
    with (V15/'results/optical_preferred_candidate_pairs.csv').open() as h:
        for r in csv.DictReader(h):
            old[(r['run'], int(r['frame_index']))] = r
    tracks, pairs, summaries, inputs, excluded = [], [], [], [], []
    for folder in sorted(DATA.glob('m??q??_???')):
        markers = list(folder.glob('__*.txt'))
        assert len(markers)==1, folder
        pattern = markers[0].stem[2:].replace('_', ' ')
        if folder.name=='m05q05_002' or pattern not in PATTERNS:
            excluded.append(dict(run=folder.name, pattern=pattern)); continue
        files = sorted(folder.glob('ImgA*.jpg'))
        ids = [int(re.search(r'(\d+)\.jpg$', p.name).group(1)) for p in files]
        assert len(files)==501 and np.all(np.diff(ids)==1), folder
        config = (folder/'_config.mii').read_text()
        fps = float(re.search(r'Film speed 1\s*:([0-9.]+)', config).group(1))
        xsv = (folder/'_config.xsv').read_text()
        assert fps == float(re.search(r'^Rate=([0-9.]+)', xsv, re.M).group(1))
        frames = np.stack([np.asarray(Image.open(f).convert('L')) for f in files])
        assert frames.shape == (501, 80, 1016)
        background = frames.mean(axis=0)
        saturation = float(np.mean(frames[:, 20:60]>=250))
        sampled = np.linspace(0, 499, N_PAIRS, dtype=int)
        inp = dict(run=folder.name, pattern=pattern, fps=fps, nframes=len(files),
                   duration_s=500/fps, core_saturation=saturation,
                   config_sha256=hashlib.sha256((folder/'_config.mii').read_bytes()).hexdigest(),
                   frame_content_sha256=hashlib.sha256(frames.tobytes()).hexdigest())
        inputs.append(inp)
        run_pairs = []
        for i in sampled:
            for method in METHODS:
                a, b = frames[i].astype(float), frames[i+1].astype(float)
                if method=='background_subtracted': a-=background; b-=background
                local = []
                for x, y in WINDOWS:
                    f, r, rt, ok = track(a, b, x, y)
                    bx, by = x+f['ix'], y+f['iy']
                    local_sat = max(float(np.mean(frames[i,y:y+HEIGHT,x:x+WIDTH]>=250)),
                                    float(np.mean(frames[i+1,by:by+HEIGHT,bx:bx+WIDTH]>=250)))
                    accepted = bool(ok and local_sat<=SAT_LIMIT)
                    row = dict(run=folder.name, pattern=pattern, method=method,
                               frame_index=int(i), next_frame_index=int(i+1),
                               window_x=x, window_y=y, width=WIDTH, height=HEIGHT,
                               fps=fps, dx_px=f['dx'], dy_px=f['dy'],
                               axial_px_s=f['dx']*fps, peak=f['peak'],
                               margin=f['margin'], texture_std=f['texture'],
                               search_boundary=f['boundary'], forward_pass=f['passed'],
                               reverse_dx_px=r['dx'] if r else np.nan,
                               reverse_dy_px=r['dy'] if r else np.nan,
                               reverse_peak=r['peak'] if r else np.nan,
                               reverse_margin=r['margin'] if r else np.nan,
                               reverse_pass=r['passed'] if r else False,
                               roundtrip_px=rt, local_saturation=local_sat, accepted=accepted)
                    tracks.append(row); local.append(row)
                good = [r for r in local if r['accepted']]
                d = np.array([r['dx_px'] for r in good])
                med = float(np.median(d)) if len(d) else np.nan
                dev = float(np.max(np.abs(d-med))) if len(d) else np.nan
                nx = len({r['window_x'] for r in good}); ny = len({r['window_y'] for r in good})
                tol = max(2., .20*abs(med)) if len(d) else np.nan
                distributed = len(d)>=3 and nx==2 and ny>=2
                coherent = bool(distributed and dev<=tol)
                baseline = old[(folder.name, int(i))]
                prow = dict(run=folder.name, pattern=pattern, method=method, frame_index=int(i),
                            accepted_patches=len(d), available_patches=len(WINDOWS),
                            distinct_axial_windows=nx, distinct_radial_windows=ny,
                            median_axial_px_frame=med, median_axial_px_s=med*fps,
                            max_patch_deviation_px=dev, agreement_tolerance_px=tol,
                            distributed=distributed, coherent=coherent,
                            global_1d_pass=baseline['heuristic_pass']=='True',
                            global_1d_dx_px=float(baseline['displacement_px_per_frame']))
                pairs.append(prow); run_pairs.append(prow)
        for method in METHODS:
            p = [r for r in run_pairs if r['method']==method]
            good = [r for r in p if r['coherent']]
            d = np.array([r['median_axial_px_frame'] for r in good])
            summary = dict(run=folder.name, pattern=pattern, method=method, fps=fps,
                           sampled_pairs=len(p), attempted_patch_tracks=len(p)*len(WINDOWS),
                           accepted_patch_tracks=sum(r['accepted_patches'] for r in p),
                           pairs_with_any_patch=sum(r['accepted_patches']>0 for r in p),
                           distributed_pairs=sum(r['distributed'] for r in p), coherent_pairs=len(good),
                           coherent_pair_fraction=len(good)/len(p),
                           median_coherent_axial_px_frame=float(np.median(d)) if len(d) else np.nan,
                           q25_coherent_axial_px_frame=float(np.quantile(d,.25)) if len(d) else np.nan,
                           q75_coherent_axial_px_frame=float(np.quantile(d,.75)) if len(d) else np.nan,
                           core_saturation=saturation, supported=bool(len(good)>=len(p)/2 and saturation<=SAT_LIMIT),
                           global_1d_passing_sampled_pairs=sum(r['global_1d_pass'] for r in p))
            summary['median_coherent_axial_px_s']=summary['median_coherent_axial_px_frame']*fps
            summaries.append(summary)
        print(folder.name, [(r['method'], r['coherent_pairs'], r['accepted_patch_tracks'])
                            for r in summaries[-2:]], flush=True)
    assert len(inputs)==26 and len(tracks)==26*N_PAIRS*len(WINDOWS)*len(METHODS)
    write_csv('patch_tracks.csv', tracks); write_csv('patch_pairs.csv', pairs)
    write_csv('patch_summary.csv', summaries)
    validation.update(inputs=inputs, excluded=excluded,
                      methods=list(METHODS), windows=WINDOWS, samples=N_PAIRS,
                      search_dx=DX, search_dy=DY, min_correlation=MIN_CORR,
                      min_margin=MIN_MARGIN, min_texture=MIN_TEXTURE,
                      roundtrip_limit_px=ROUNDTRIP_LIMIT, saturation_limit=SAT_LIMIT,
                      coherent_definition='At least 3 accepted patches spanning both axial windows and >=2 radial windows; all axial shifts within max(2 px, 20% of absolute median).',
                      support_definition='>=50% sampled pairs coherent and whole-run native core saturation <=20%. Heuristic consistency, not physical validation.',
                      rows=dict(tracks=len(tracks), pairs=len(pairs), summaries=len(summaries)),
                      elapsed_s=time.monotonic()-start)
    (OUT/'local_patch_validation.json').write_text(json.dumps(validation, indent=2))
    print('complete',validation['rows'],'seconds',validation['elapsed_s'],flush=True)


if __name__=='__main__': run()
