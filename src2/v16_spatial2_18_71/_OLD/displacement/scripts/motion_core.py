"""Signed image motion with explicit correspondence checks; no forced motion."""
from __future__ import annotations

import numpy as np
from scipy.signal import correlate

# Fixed before the experimental comparison. These are correspondence heuristics,
# not calibrated uncertainties. Settings follow the earlier native-frame study.
SEARCH = 200
MIN_CORRELATION = 0.60
MIN_MARGIN = 0.03
MIN_TEXTURE = 1.0
COMPETITOR_RADIUS = 5
ABS_TOLERANCE = 1.0
REL_TOLERANCE = 0.15
ROI = slice(18, 72)
FIT_ROWS = np.arange(18, 72, 2)
CHECK_ROWS = np.arange(19, 72, 2)


def ncc(a, b, search=SEARCH):
    """A fixed central template in A searches B; positive dx means rightward.

    Every displacement uses the same template length. Quadratic interpolation
    refines the integer peak, retaining the sign and allowing genuine zero.
    """
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.ndim != 1 or a.shape != b.shape or len(a) <= 2 * search + 8:
        raise ValueError('Matching finite profiles of sufficient width required.')
    if not (np.isfinite(a).all() and np.isfinite(b).all()):
        raise ValueError('Non-finite profile.')
    template = a[search:-search].copy()
    template -= template.mean()
    n = len(template)
    numerator = correlate(b, template, mode='valid', method='fft')
    c1 = np.r_[0., np.cumsum(b)]
    c2 = np.r_[0., np.cumsum(b * b)]
    sums = c1[n:] - c1[:-n]
    variances = np.maximum(c2[n:] - c2[:-n] - sums * sums / n, 0.)
    template_energy = np.dot(template, template)
    denominator = np.sqrt(template_energy * variances)
    values = np.full(numerator.shape, np.nan)
    np.divide(numerator, denominator, out=values, where=denominator > 1e-9)
    empty = dict(dx=np.nan, peak=np.nan, margin=np.nan, texture=0.,
                 boundary=True, passed=False)
    if not np.isfinite(values).any():
        return empty
    k = int(np.nanargmax(values))
    offsets = np.arange(-search, search + 1)
    subpixel = 0.
    if 0 < k < len(values) - 1:
        curvature = values[k - 1] - 2 * values[k] + values[k + 1]
        if np.isfinite(curvature) and curvature < -1e-12:
            subpixel = float(np.clip(.5 * (values[k - 1] - values[k + 1]) / curvature, -.5, .5))
    peak = float(values[k])
    competitor = values[np.abs(offsets - offsets[k]) > COMPETITOR_RADIUS]
    margin = float(peak - np.nanmax(competitor))
    texture = float(min(np.sqrt(template_energy / n), np.sqrt(variances[k] / n)))
    boundary = bool(abs(offsets[k]) >= search - 2)
    passed = bool(peak >= MIN_CORRELATION and margin >= MIN_MARGIN
                  and texture >= MIN_TEXTURE and not boundary)
    return dict(dx=float(offsets[k] + subpixel), peak=peak, margin=margin,
                texture=texture, boundary=boundary, passed=passed)


def dynamic_profiles(frames, rows):
    profiles = frames[:, rows, :].mean(axis=1, dtype=np.float64)
    # Offline background estimate uses all observed frames, explicitly. No
    # smoothing or expected positive velocity is imposed on the trajectory.
    return profiles - profiles.mean(axis=0, keepdims=True)


def estimate(frames):
    """Fit on alternate rows; check other rows, reversed order, and frame gap 2.

    Reversing central templates is a spatial-consistency test, not a same-patch
    round trip. Neighboring rows/frames are correlated, not independent trials.
    """
    fit = dynamic_profiles(frames, FIT_ROWS)
    check = dynamic_profiles(frames, CHECK_ROWS)
    forward = [ncc(a, b) for a, b in zip(fit[:-1], fit[1:])]
    reverse = [ncc(b, a) for a, b in zip(fit[:-1], fit[1:])]
    other_rows = [ncc(a, b) for a, b in zip(check[:-1], check[1:])]
    long = [ncc(a, b) for a, b in zip(fit[:-2], fit[2:])]
    count = len(forward)
    dx = np.array([r['dx'] for r in forward])
    tolerances = np.maximum(ABS_TOLERANCE, REL_TOLERANCE * abs(dx))
    row_error = abs(dx - np.array([r['dx'] for r in other_rows]))
    reverse_error = abs(dx + np.array([r['dx'] for r in reverse]))
    base = np.array([a['passed'] and b['passed'] and c['passed']
                     for a, b, c in zip(forward, reverse, other_rows)])
    base &= (row_error <= tolerances) & (reverse_error <= tolerances)
    checks = np.zeros(count, dtype=int)
    failures = np.zeros(count, dtype=int)
    temporal = []
    for i, result in enumerate(long):
        total = dx[i] + dx[i + 1]
        tol = max(2 * ABS_TOLERANCE, REL_TOLERANCE * abs(total))
        error = abs(result['dx'] - total)
        comparable = bool(base[i] and base[i + 1] and result['passed'])
        agrees = bool(comparable and error <= tol)
        if comparable:
            checks[i:i + 2] += 1
            if not agrees:
                failures[i:i + 2] += 1
        temporal.append(dict(frame_i=i, frame_j=i + 2, direct_dx=result['dx'],
                             adjacent_sum_dx=float(total), absolute_difference=error,
                             tolerance=tol, comparable=comparable, agrees=agrees,
                             direct_peak=result['peak'], direct_pass=result['passed']))
    supported = base & (checks > 0) & (failures == 0)
    records = []
    for i, (f, r, c) in enumerate(zip(forward, reverse, other_rows)):
        reasons = []
        if not f['passed']: reasons.append('fit_quality')
        if not r['passed'] or reverse_error[i] > tolerances[i]: reasons.append('reverse_consistency')
        if not c['passed'] or row_error[i] > tolerances[i]: reasons.append('row_consistency')
        if checks[i] == 0: reasons.append('no_temporal_support')
        if failures[i]: reasons.append('temporal_disagreement')
        records.append(dict(frame_i=i, frame_j=i + 1, candidate_dx=f['dx'],
                            new_dx=f['dx'] if supported[i] else np.nan,
                            supported=bool(supported[i]), fit_pass=f['passed'],
                            fit_peak=f['peak'], peak_margin=f['margin'],
                            texture_std=f['texture'], search_boundary=f['boundary'],
                            check_rows_dx=c['dx'], check_rows_peak=c['peak'],
                            check_rows_pass=c['passed'], row_difference=float(row_error[i]),
                            reverse_dx=r['dx'], reverse_peak=r['peak'],
                            reverse_pass=r['passed'], reverse_sum=float(reverse_error[i]),
                            agreement_tolerance=float(tolerances[i]),
                            base_supported=bool(base[i]), temporal_checks=int(checks[i]),
                            temporal_failures=int(failures[i]), reason=';'.join(reasons)))
    return records, temporal


def sample_columns(image, x):
    x = np.asarray(x, float)
    if np.min(x) < 0 or np.max(x) > image.shape[1] - 1:
        raise ValueError('Sample outside the native image.')
    left = np.floor(x).astype(int)
    right = np.minimum(left + 1, image.shape[1] - 1)
    fraction = x - left
    return image[:, left] * (1 - fraction) + image[:, right] * fraction


def normalized_correlation(a, b):
    a = np.asarray(a, float).ravel(); b = np.asarray(b, float).ravel()
    a = a - a.mean(); b = b - b.mean()
    norm = np.linalg.norm(a) * np.linalg.norm(b)
    return float(np.dot(a, b) / norm) if norm > 1e-9 else np.nan


def image_alignment(frames, background, i, dx, rows=CHECK_ROWS):
    """2-D check on rows unused for fitting, over identical source coordinates."""
    if not np.isfinite(dx):
        return dict(dynamic_correlation=np.nan, raw_correlation=np.nan)
    x = np.arange(SEARCH, frames.shape[2] - SEARCH, dtype=float)
    a = frames[i, rows].astype(float)
    b = frames[i + 1, rows].astype(float)
    source = a[:, SEARCH:-SEARCH]
    shifted = sample_columns(b, x + dx)
    da = source - background[rows, SEARCH:-SEARCH]
    db = shifted - sample_columns(background[rows], x + dx)
    return dict(dynamic_correlation=normalized_correlation(da, db),
                raw_correlation=normalized_correlation(source, shifted))


def longest_segment(records):
    """Longest supported monotone segment; first one wins ties. No gap filling."""
    best = None; start = None; direction = None
    for i in range(len(records) + 1):
        r = records[i] if i < len(records) else None
        ok = r is not None and r['supported'] and abs(r['new_dx']) >= .5
        sign = int(np.sign(r['new_dx'])) if ok else None
        if start is not None and (not ok or sign != direction):
            candidate = (start, i)  # inclusive frame endpoints, exclusive transition end
            if best is None or candidate[1] - candidate[0] > best[1] - best[0]:
                best = candidate
            start = None
        if ok and start is None:
            start = i; direction = sign
    return best


def reconstruct_segment(frames, displacements):
    """Cumulative subpixel positions and midpoint cells, sampled every 1 px.

    Support stops at the first/last observed frame center. No extrapolated final
    frame shift and no artificial strip width, feathering or bridge across gaps.
    """
    dx = np.asarray(displacements, float)
    if len(frames) != len(dx) + 1 or not np.isfinite(dx).all():
        raise ValueError('A finite displacement for every included interval is required.')
    positions = np.r_[0., np.cumsum(-dx)]
    if not (np.all(np.diff(positions) > 0) or np.all(np.diff(positions) < 0)):
        raise ValueError('Reconstruction requires a monotone supported segment.')
    order = np.argsort(positions)
    positions = positions[order]
    positions -= positions[0]
    grid = np.arange(0., positions[-1], 1.)
    if len(grid) < 8:
        raise ValueError('Less than eight reconstructed pixels.')
    midpoints = (positions[:-1] + positions[1:]) / 2
    owners_sorted = np.searchsorted(midpoints, grid, side='right')
    owners = order[owners_sorted]
    columns = (frames.shape[2] - 1) / 2 + grid - positions[owners_sorted]
    if columns.min() < 0 or columns.max() > frames.shape[2] - 1:
        raise ValueError('Estimated motion exceeds available image support.')
    lo = np.floor(columns).astype(int); hi = np.minimum(lo + 1, frames.shape[2] - 1)
    frac = columns - lo
    canvas = (frames[owners, :, lo].T * (1 - frac)
              + frames[owners, :, hi].T * frac)
    rendered = np.rint(np.clip(canvas, 0, 255)).astype(np.uint8)
    return rendered, dict(position_px=grid, source_frame=owners,
                          source_column=columns, frame_positions=positions,
                          total_travel_px=float(positions[-1]))
