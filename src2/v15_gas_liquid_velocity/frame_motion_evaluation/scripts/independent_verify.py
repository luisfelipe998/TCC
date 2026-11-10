#!/usr/bin/env python3
"""Independently reaggregate saved frame-motion evidence; no source mutation.

This verifier does not import the analysis implementations. It checks arithmetic,
scope and declared decision rules, not experimentally accurate material velocity.
"""
from __future__ import annotations
import csv
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np

BASE = Path(__file__).resolve().parents[1]
V15 = BASE.parent
DATA = V15.parent / 'data_processed'
OUT = BASE / 'results'
PATTERNS = {'dispersed bubbles', 'plug', 'slug', 'churn', 'wavy annular'}


def read_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def true(value):
    assert value in ('True', 'False'), value
    return value == 'True'


class Checks:
    def __init__(self):
        self.numeric_cells = 0
        self.maximum_difference = 0.

    def same(self, actual, expected):
        a, b = float(actual), float(expected)
        self.numeric_cells += 1
        if math.isnan(a) or math.isnan(b):
            assert math.isnan(a) and math.isnan(b), (a, b)
            return
        difference = abs(a-b)
        self.maximum_difference = max(self.maximum_difference, difference)
        assert math.isclose(a, b, rel_tol=1e-11, abs_tol=1e-9), (a, b)


def quantiles(values):
    return np.quantile(values, [.25, .5, .75]) if len(values) else [np.nan]*3


def inventory():
    selected, excluded = [], []
    for folder in sorted(DATA.glob('m??q??_???')):
        markers = list(folder.glob('__*.txt'))
        assert len(markers) == 1
        pattern = markers[0].stem[2:].replace('_', ' ')
        if folder.name == 'm05q05_002' or pattern not in PATTERNS:
            excluded.append(dict(run=folder.name, pattern=pattern))
            continue
        files = sorted(folder.glob('ImgA*.jpg'))
        ids = [int(re.search(r'([0-9]+)[.]jpg$', p.name).group(1)) for p in files]
        mii = (folder / '_config.mii').read_text()
        fps = float(re.search(r'Film speed 1\s*:([0-9.]+)', mii).group(1))
        nf = int(re.search(r'Number of images 1\s*:([0-9]+)', mii).group(1))
        acquisition = (folder / '_config.xsv').read_text().split('[', 2)[1]
        values = dict(line.split('=', 1) for line in acquisition.splitlines() if '=' in line)
        assert len(files) == nf == int(values['Images']) == 501
        assert all(b-a == 1 for a, b in zip(ids, ids[1:]))
        assert fps == float(values['Rate'])
        selected.append(dict(run=folder.name, pattern=pattern, frames=nf, fps=fps,
                             first_id=ids[0], last_id=ids[-1], duration_s=500/fps))
    assert len(selected) == 26 and len(excluded) == 5
    return selected, excluded


def verify_baseline(rates):
    checks = Checks()
    pairs = read_csv(V15 / 'results/optical_preferred_candidate_pairs.csv')
    groups = defaultdict(list)
    for row in pairs:
        run = row['run']
        shift, lag = float(row['shift_px']), int(row['lag_frames'])
        checks.same(row['displacement_px_per_frame'], shift/lag)
        checks.same(row['structure_celerity_px_per_s'], shift*rates[run]/lag)
        groups[run].append(row)
    summary = read_csv(V15 / 'results/optical_preferred_candidate_per_run.csv')
    for row in summary:
        current = groups[row['run']]
        assert len(current) == int(row['pairs']) == 500
        good = [float(r['displacement_px_per_frame']) for r in current if true(r['heuristic_pass'])]
        assert len(good) == int(row['heuristic_pass_pairs'])
        lo, med, hi = quantiles(good)
        checks.same(row['median_pass_px_per_frame'], med)
        checks.same(row['q25_pass_px_per_frame'], lo)
        checks.same(row['q75_pass_px_per_frame'], hi)
        checks.same(row['median_pass_structure_celerity_px_per_s'], med*rates[row['run']])
    assert len(pairs) == 13000
    return dict(pair_rows=len(pairs), numeric_cells=checks.numeric_cells,
                maximum_difference=checks.maximum_difference,
                passing_pairs=sum(true(r['heuristic_pass']) for r in pairs),
                supported_runs=sum(not true(r['low_coverage']) and not true(r['high_saturation']) for r in summary))


def verify_local(rates):
    checks = Checks()
    tracks = read_csv(OUT / 'local_patch_tracks.csv')
    pairs = read_csv(OUT / 'local_patch_pairs.csv')
    summary = read_csv(OUT / 'local_patch_summary.csv')
    groups, pairgroups = defaultdict(list), defaultdict(list)
    for row in tracks:
        assert float(row['fps']) == rates[row['run']]
        assert int(row['next_frame_index'])-int(row['frame_index']) == 1
        dx, dy = float(row['dx_px']), float(row['dy_px'])
        checks.same(row['axial_px_s'], dx*rates[row['run']])
        forward = float(row['peak']) >= .60 and float(row['margin']) >= .03 and float(row['texture_std']) >= 1 and not true(row['search_boundary'])
        assert true(row['forward_pass']) == forward
        if forward:
            residual = math.hypot(round(dx)+float(row['reverse_dx_px']), round(dy)+float(row['reverse_dy_px']))
            checks.same(row['roundtrip_px'], residual)
        accepted = forward and true(row['reverse_pass']) and float(row['roundtrip_px']) <= 1 and float(row['local_saturation']) <= .20
        assert true(row['accepted']) == accepted
        groups[(row['run'], row['method'], row['frame_index'])].append(row)
    for row in pairs:
        local = groups[(row['run'], row['method'], row['frame_index'])]
        assert len(local) == int(row['available_patches']) == 6
        good = [r for r in local if true(r['accepted'])]
        assert len(good) == int(row['accepted_patches'])
        nx, ny = len({r['window_x'] for r in good}), len({r['window_y'] for r in good})
        assert nx == int(row['distinct_axial_windows']) and ny == int(row['distinct_radial_windows'])
        shifts = np.array([float(r['dx_px']) for r in good])
        med = float(np.median(shifts)) if len(shifts) else np.nan
        deviation = float(np.max(np.abs(shifts-med))) if len(shifts) else np.nan
        tolerance = max(2., .20*abs(med)) if len(shifts) else np.nan
        checks.same(row['median_axial_px_frame'], med)
        checks.same(row['median_axial_px_s'], med*rates[row['run']])
        checks.same(row['max_patch_deviation_px'], deviation)
        checks.same(row['agreement_tolerance_px'], tolerance)
        distributed = len(good) >= 3 and nx == 2 and ny >= 2
        coherent = distributed and deviation <= tolerance
        assert true(row['distributed']) == distributed
        assert true(row['coherent']) == coherent
        pairgroups[(row['run'], row['method'])].append(row)
    for row in summary:
        current = pairgroups[(row['run'], row['method'])]
        good = [r for r in current if true(r['coherent'])]
        assert len(current) == int(row['sampled_pairs']) == 50
        assert int(row['attempted_patch_tracks']) == 300
        assert int(row['accepted_patch_tracks']) == sum(int(r['accepted_patches']) for r in current)
        assert int(row['pairs_with_any_patch']) == sum(int(r['accepted_patches']) > 0 for r in current)
        assert int(row['distributed_pairs']) == sum(true(r['distributed']) for r in current)
        assert int(row['coherent_pairs']) == len(good)
        checks.same(row['coherent_pair_fraction'], len(good)/len(current))
        lo, med, hi = quantiles([float(r['median_axial_px_frame']) for r in good])
        checks.same(row['median_coherent_axial_px_frame'], med)
        checks.same(row['q25_coherent_axial_px_frame'], lo)
        checks.same(row['q75_coherent_axial_px_frame'], hi)
        checks.same(row['median_coherent_axial_px_s'], med*rates[row['run']])
        assert true(row['supported']) == (len(good) >= 25 and float(row['core_saturation']) <= .20)
    assert len(tracks) == 15600 and len(pairs) == 2600 and len(summary) == 52
    methods = {}
    for method in ['raw', 'background_subtracted']:
        selected = [r for r in summary if r['method'] == method]
        methods[method] = dict(accepted_tracks=sum(int(r['accepted_patch_tracks']) for r in selected),
                               coherent_pairs=sum(int(r['coherent_pairs']) for r in selected),
                               supported_runs=[r['run'] for r in selected if true(r['supported'])])
    return dict(track_rows=len(tracks), pair_rows=len(pairs), summary_rows=len(summary),
                numeric_cells=checks.numeric_cells, maximum_difference=checks.maximum_difference, methods=methods)


def verify_temporal(rates):
    checks = Checks()
    pairs = read_csv(OUT / 'temporal_pairs.csv')
    summary = read_csv(OUT / 'temporal_per_run.csv')
    metadata = json.loads((OUT / 'temporal_metadata.json').read_text())
    groups = defaultdict(list)
    for row in pairs:
        groups[(row['run'], int(row['lag_frames']))].append(row)
    stats = {}
    for run in rates:
        adjacent = groups[(run, 1)]
        assert [int(r['frame_index']) for r in adjacent] == list(range(500))
        short_shifts = np.array([float(r['shift_px']) for r in adjacent])
        short_passes = np.array([true(r['heuristic_pass']) for r in adjacent])
        reversal_errors = []
        for lag in (1, 2, 4):
            current = groups[(run, lag)]
            assert [int(r['frame_index']) for r in current] == list(range(501-lag))
            differences, relative_differences, agreements, speeds = [], [], [], []
            for i, row in enumerate(current):
                assert float(row['fps']) == rates[run]
                assert int(row['second_frame_id'])-int(row['first_frame_id']) == lag
                assert int(row['first_frame_id']) == i
                d = float(row['shift_px'])
                speed = d/lag
                checks.same(row['displacement_px_per_frame'], speed)
                checks.same(row['structure_celerity_px_per_s'], speed*rates[run])
                short_mean = float(np.mean(short_shifts[i:i+lag]))
                path_pass = bool(np.all(short_passes[i:i+lag]))
                joint = true(row['heuristic_pass']) and path_pass
                difference = abs(speed-short_mean)
                relative = difference/max(abs(short_mean), 1.)
                tolerance = max(1., .15*abs(short_mean))
                checks.same(row['adjacent_path_mean_px_frame'], short_mean)
                checks.same(row['temporal_abs_difference_px_frame'], difference)
                checks.same(row['temporal_relative_difference'], relative)
                checks.same(row['temporal_tolerance_px_frame'], tolerance)
                assert true(row['adjacent_path_pass']) == path_pass
                assert true(row['temporal_joint_pass']) == joint
                assert true(row['temporal_within_tolerance']) == (joint and difference <= tolerance)
                if joint:
                    differences.append(difference)
                    relative_differences.append(relative)
                    agreements.append(difference <= tolerance)
                    speeds.append(speed)
                if lag == 1:
                    reverse_error = abs(d+float(row['reversal_shift_px']))
                    reverse_joint = joint and true(row['reversal_pass'])
                    checks.same(row['reversal_abs_sum_px'], reverse_error)
                    assert true(row['reversal_joint_pass']) == reverse_joint
                    if reverse_joint:
                        reversal_errors.append(reverse_error)
                else:
                    assert row['reversal_shift_px'] == row['reversal_pass'] == row['reversal_joint_pass'] == row['reversal_abs_sum_px'] == ''
            stats[(run, lag)] = dict(pairs=len(current),
                passing=sum(true(r['heuristic_pass']) for r in current),
                joint=len(differences), joint_fraction=len(differences)/len(current),
                med_abs=np.median(differences) if differences else np.nan,
                med_relative=np.median(relative_differences) if relative_differences else np.nan,
                agreement=np.mean(agreements) if agreements else np.nan,
                med_speed=np.median(speeds) if speeds else np.nan)
        stats[(run, 'reversal')] = (len(reversal_errors), np.median(reversal_errors) if reversal_errors else np.nan)
    stable_runs = []
    for row in summary:
        run = row['run']
        assert int(row['frames']) == 501 and float(row['fps']) == rates[run]
        checks.same(row['first_to_last_s'], 500/rates[run])
        current = groups[(run, 1)]
        good = [float(r['shift_px']) for r in current if true(r['heuristic_pass'])]
        lo, med, hi = quantiles(good)
        assert int(row['adjacent_pairs']) == 500
        assert int(row['adjacent_pass_pairs']) == len(good)
        checks.same(row['adjacent_pass_fraction'], len(good)/500)
        for prefix, expected in [('q25', lo), ('median', med), ('q75', hi)]:
            checks.same(row[prefix+'_signed_px_frame'], expected)
            checks.same(row[prefix+'_signed_px_s'], expected*rates[run])
        checks.same(row['estimated_motion_during_exposure_px'], abs(med)*rates[run]*float(row['exposure_us'])*1e-6)
        saturated = float(row['core_saturation_fraction']) > .20
        assert true(row['high_saturation']) == saturated
        base_supported = len(good) >= 250 and not saturated
        assert true(row['base_coverage_and_saturation_pass']) == base_supported
        lag_flags = []
        for lag in (2, 4):
            st = stats[(run, lag)]
            prefix = 'lag'+str(lag)+'_'
            for key, skey in [('pairs', 'pairs'), ('pass_pairs', 'passing'), ('joint_pairs', 'joint'), ('joint_fraction', 'joint_fraction'),
                              ('median_abs_difference_px_frame', 'med_abs'), ('median_relative_difference', 'med_relative'),
                              ('within_tolerance_fraction', 'agreement'), ('joint_median_px_frame', 'med_speed')]:
                checks.same(row[prefix+key] if row[prefix+key] else np.nan, st[skey])
            lag_flags.append(st['joint_fraction'] >= .25 and st['agreement'] >= .80)
        stable = base_supported and all(lag_flags)
        assert true(row['stable_under_lag_diagnostics']) == stable
        if stable:
            stable_runs.append(run)
        count, revmed = stats[(run, 'reversal')]
        assert int(row['reversal_joint_pairs']) == count
        checks.same(row['reversal_median_abs_sum_px'] if row['reversal_median_abs_sum_px'] else np.nan, revmed)
    totals = dict(runs=26, adjacent_pairs=13000,
                  adjacent_pass=sum(stats[(r, 1)]['passing'] for r in rates),
                  base_supported=sum(true(r['base_coverage_and_saturation_pass']) for r in summary),
                  lag_stable=len(stable_runs), temporal_pairs=len(pairs))
    assert metadata['totals'] == totals
    assert len(pairs) == 38896 and len(summary) == 26
    return dict(pair_rows=len(pairs), summary_rows=len(summary), totals=totals,
                stable_runs=stable_runs, numeric_cells=checks.numeric_cells,
                maximum_difference=checks.maximum_difference)


def summarize_existing_sensitivity():
    existing = read_csv(V15 / 'results/optical_sensitivity.csv')
    baseline = {r['run']: r for r in read_csv(V15 / 'results/optical_preferred_candidate_per_run.csv')}
    groups = defaultdict(list)
    for row in existing:
        if row['method'] == 'ncc_dynamic':
            groups[row['run']].append(row)
    result = []
    for run, rows in sorted(groups.items()):
        assert len(rows) == 9
        reference = float(baseline[run]['median_pass_px_per_frame'])
        values = [float(r['median_pass_px_per_frame']) for r in rows if float(r['heuristic_pass_fraction']) >= .5]
        regions = [float(r['median_pass_px_per_frame']) for r in rows if r['lag_frames'] == '1' and float(r['heuristic_pass_fraction']) >= .5]
        lags = [float(r['median_pass_px_per_frame']) for r in rows if r['roi'] == '20:60' and float(r['heuristic_pass_fraction']) >= .5]
        def span(values):
            return 100*(max(values)-min(values))/abs(reference) if len(values) > 1 else None
        result.append(dict(run=run, pattern=rows[0]['pattern'],
            baseline_support=not true(baseline[run]['low_coverage']) and not true(baseline[run]['high_saturation']),
            baseline_px_frame=reference, supported_variants=len(values),
            min_support=min(float(r['heuristic_pass_fraction']) for r in rows),
            max_support=max(float(r['heuristic_pass_fraction']) for r in rows),
            max_boundary_fraction=max(float(r['boundary_fraction']) for r in rows),
            supported_range_pct_baseline=span(values), lag1_supported_roi_count=len(regions),
            lag1_roi_range_pct=span(regions), center_supported_lag_count=len(lags),
            center_lag_range_pct=span(lags)))
    (OUT / 'independent_existing_sensitivity.json').write_text(json.dumps(result, indent=2))
    return dict(variants=26*9, all_nine_supported=sum(r['baseline_support'] and r['supported_variants'] == 9 for r in result),
                supported_span_over_10pct=[r['run'] for r in result if r['baseline_support'] and r['supported_range_pct_baseline'] is not None and r['supported_range_pct_baseline'] > 10])


def main():
    selected, excluded = inventory()
    rates = {r['run']: r['fps'] for r in selected}
    result = dict(scope=selected, excluded=excluded,
                  rate_counts=dict(Counter(r['fps'] for r in selected)),
                  baseline=verify_baseline(rates), local_patch=verify_local(rates),
                  temporal=verify_temporal(rates),
                  existing_sensitivity=summarize_existing_sensitivity(),
                  meaning='Scope, arithmetic and declared gates verified; these checks do not establish physical accuracy or phase identity.')
    (OUT / 'independent_verification.json').write_text(json.dumps(result, indent=2))
    print(json.dumps({k: v for k, v in result.items() if k not in ('scope', 'excluded')}, indent=2))


if __name__ == '__main__':
    main()
