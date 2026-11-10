#!/usr/bin/env python3
"""Independent saved-evidence checks; no analysis/plot generator is imported."""
from collections import defaultdict
from pathlib import Path
import hashlib
import itertools
import json
import math
import re
import numpy as np
import pandas as pd
from CoolProp.CoolProp import PropsSI

BASE = Path(__file__).resolve().parents[1]
SRC2 = BASE.parent
OLD = SRC2 / 'v15_gas_liquid_velocity'
PHASE = OLD / 'phase_conditioned_motion'
RESULT = BASE / 'results'
CONFIG = 'consensus_reference_guard'
SCALE = 1e-3 / 53
COUNTS = defaultdict(int)
ERRORS = defaultdict(float)
INPUTS = {}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def remember(path):
    INPUTS[str(path.relative_to(SRC2))] = sha(path)
    return path


def table(path):
    return pd.read_csv(remember(path))


def same(actual, expected, group, atol=1e-9, rtol=2e-12):
    a, b = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
    assert a.shape == b.shape, (group, a.shape, b.shape)
    assert np.array_equal(np.isnan(a), np.isnan(b)), group
    assert np.allclose(a, b, atol=atol, rtol=rtol, equal_nan=True), (group, a, b)
    finite = np.isfinite(a) & np.isfinite(b)
    if finite.any():
        ERRORS[group] = max(ERRORS[group], float(np.max(np.abs(a[finite]-b[finite]))))
    COUNTS[group] += a.size


def verify_scope():
    inv = table(RESULT / 'scenario_inventory.csv').set_index('run')
    selected = inv[inv.included]
    assert len(inv) == 31 and len(selected) == 26
    assert set(inv[~inv.included].index) == {
        'm03q01_001', 'm03q80_001', 'm04q05_001', 'm05q05_002', 'm07q10_001'}
    calibration = table(RESULT / 'calibration_inventory.csv').set_index('run')
    for run, row in inv.iterrows():
        folder = SRC2 / 'data_processed' / run
        xsv = remember(folder / '_config.xsv')
        mii = remember(folder / '_config.mii')
        assert sha(xsv) == row.config_sha256 and sha(mii) == row.mii_sha256
        fps = float(re.search(r'^Rate=([0-9.]+)', xsv.read_text(), re.M)[1])
        other = float(re.search(r'Film speed 1[ ]*:([0-9.]+)', mii.read_text())[1])
        same([row.fps, row.frames, row.first_to_last_s, row.conditional_axial_scale_m_px],
             [fps, 501, 500/fps, SCALE], 'scope')
        assert other == fps
        ids = sorted(int(p.stem.split('ImgA')[-1]) for p in folder.glob('ImgA*.jpg'))
        assert ids == list(range(501))
        if run in calibration.index:
            c = calibration.loc[run]
            same([c.native_fps, c.frames, c.interval_count, c.duration_s],
                 [fps, 501, 500, 500/fps], 'calibration_crosscheck')
            assert c.config_sha256 == row.config_sha256
            assert c.mii_sha256 == row.mii_sha256
            assert not c.documented_camera_m_per_px
    return selected


def verify_passage(scope):
    frames = table(PHASE / 'results/passage_frames.csv')
    frames = frames[frames.config == CONFIG].sort_values(['run', 'frame_id'])
    pairs = table(RESULT / 'passage_pairs_reference.csv')
    summary = table(RESULT / 'passage_per_run.csv').set_index('run')
    temporal = table(OLD / 'frame_motion_evaluation/results/temporal_pairs.csv')
    temporal = temporal.set_index(['run', 'lag_frames', 'frame_index'])
    optical = table(OLD / 'frame_motion_evaluation/results/temporal_per_run.csv').set_index('run')
    assert len(pairs) == 13000 and set(pairs.run) == set(scope.index)
    assert set(summary.index) == set(scope.index) and len(summary) == 26
    assert set(pairs.config) == {CONFIG} and set(pairs.motion_gate) == {'lag2_confirmed'}
    frame_groups = {}
    for run, inv in scope.iterrows():
        f = frames[frames.run == run].copy()
        assert list(f.frame_id) == list(range(501))
        frame_groups[run] = f
        raw = []
        for row in f.itertuples():
            reference = row.reference_mean_abs_error <= .12 and row.reference_p95_abs_error <= .30
            assert bool(row.liquid_reference_like) == reference
            occupancy = [row.occupancy_060, row.occupancy_070, row.occupancy_080]
            value = 'gas_candidate' if min(occupancy) >= .6 else 'unresolved'
            if max(occupancy) <= .05 and reference:
                value = 'liquid_candidate'
            raw.append(value)
        assert list(f.raw_state) == raw
        states = np.array(['unresolved']*501, dtype=object)
        epids = np.full(501, -1, dtype=int)
        next_id = 0
        for label, group in itertools.groupby(enumerate(raw), key=lambda item: item[1]):
            indices = [item[0] for item in group]
            if label != 'unresolved' and len(indices) >= 5:
                states[indices] = label
                epids[indices] = next_id
                next_id += 1
        assert np.array_equal(states, f.state)
        same(f.episode_id.to_numpy(), epids, 'episode_identity')
        p = pairs[pairs.run == run].sort_values('frame_index')
        assert list(p.frame_index) == list(range(500))
        fps = float(inv.fps)
        adj = temporal.loc[(run, 1)].sort_index()
        lag2 = temporal.loc[(run, 2)].sort_index()
        for row in p.itertuples():
            i = row.frame_index
            label = states[i] if states[i] != 'unresolved' and states[i] == states[i+1] else 'unresolved'
            episode = epids[i] if label != 'unresolved' else -1
            through = i < 499 and label != 'unresolved' and states[i+2] == label and epids[i+2] == epids[i]
            temporal_ok = False
            if i < 499:
                a = adj.loc[[i, i+1]]
                long = lag2.loc[i]
                joint = bool(long.heuristic_pass) and bool(a.heuristic_pass.all())
                delta = abs(float(long.displacement_px_per_frame) - float(a.displacement_px_per_frame.mean()))
                tolerance = max(1., .15*abs(float(a.displacement_px_per_frame.mean())))
                temporal_ok = joint and delta <= tolerance
                assert bool(long.temporal_within_tolerance) == temporal_ok
            motion = bool(adj.loc[i].heuristic_pass)
            accepted = label != 'unresolved' and motion and temporal_ok and through
            assert (row.phase, row.episode_id, bool(row.motion_pass), bool(row.lag2_temporal_pass),
                    bool(row.phase_through_lag2), bool(row.accepted)) == (
                        label, episode, motion, temporal_ok, through, accepted)
            same([row.displacement_px_frame, row.velocity_px_s],
                 [adj.loc[i].displacement_px_per_frame, adj.loc[i].displacement_px_per_frame*fps], 'pair_displacement')
            COUNTS['pair_gate_checks'] += 6
        s = summary.loc[run]
        reasons, values = [], {}
        for phase, key in [('gas_candidate', 'gas'), ('liquid_candidate', 'liquid')]:
            good = p[p.accepted & (p.phase == phase)]
            speed = list(good.velocity_px_s)
            per_episode = defaultdict(list)
            for row in good.itertuples():
                per_episode[row.episode_id].append(row.velocity_px_s)
            mean = math.fsum(speed)/len(speed) if speed else math.nan
            episode_means = [math.fsum(v)/len(v) for v in per_episode.values()]
            eq = math.fsum(episode_means)/len(episode_means) if episode_means else math.nan
            q25, q75 = np.quantile(speed, [.25, .75]) if speed else (math.nan, math.nan)
            fraction = sum(v > 0 for v in speed)/len(speed) if speed else math.nan
            same([s[key+'_pairs'], s[key+'_episodes'], s[key+'_mean_signed_px_s'],
                  s[key+'_episode_equal_mean_signed_px_s'], s[key+'_positive_fraction'],
                  s[key+'_q25_signed_px_s'], s[key+'_q75_signed_px_s']],
                 [len(speed), len(episode_means), mean, eq, fraction, q25, q75], 'passage_summary')
            values[key], values[key+'_episode'] = mean, eq
            if len(speed) < 20: reasons.append(key+'_too_few_pairs')
            if len(episode_means) < 3: reasons.append(key+'_too_few_episodes')
        g, l = values['gas'], values['liquid']
        if not (math.isfinite(g) and math.isfinite(l)):
            reasons.append('missing_phase_mean')
        elif abs(l)/fps < 1:
            reasons.append('near_zero_liquid_mean')
        elif g*l <= 0:
            reasons.append('opposed_or_zero_mean_directions')
        high_sat = bool(optical.loc[run].core_saturation_fraction > .2)
        assert bool(s.high_saturation) == high_sat
        if high_sat: reasons.append('high_saturation')
        assert bool(s.eligible_passage_comparison) == (not reasons)
        assert (s.withheld_reasons if pd.notna(s.withheld_reasons) else '') == ';'.join(reasons)
        expected = [g/l, values['gas_episode']/values['liquid_episode'], abs(g), abs(l), g*SCALE, l*SCALE] if not reasons else [math.nan]*6
        same([s.ratio_of_signed_means, s.episode_equal_ratio, s.gas_plot_magnitude_px_s,
              s.liquid_plot_magnitude_px_s, s.gas_conditional_m_s, s.liquid_conditional_m_s],
             expected, 'passage_plot_conversion')
        if not reasons:
            assert s.gas_positive_fraction == s.liquid_positive_fraction == 0
    assert int(summary.eligible_passage_comparison.sum()) == 4
    return frame_groups


def radial_run_squares(mask):
    """Independent edge-pair algorithm, unlike the generator's row counter."""
    radial = mask[:, 18:71, :].transpose(0, 2, 1).reshape(-1, 53)
    padded = np.pad(radial.astype(np.int8), ((0, 0), (1, 1)))
    edges = np.diff(padded, axis=1)
    starts_row, starts_y = np.where(edges == 1)
    ends_row, ends_y = np.where(edges == -1)
    assert np.array_equal(starts_row, ends_row)
    lengths = ends_y-starts_y
    assert np.all((lengths > 0) & (lengths <= 53))
    square_sum = np.bincount(starts_row, weights=lengths.astype(np.int64)**2,
                             minlength=radial.shape[0])
    alpha = square_sum.reshape(mask.shape[0], mask.shape[2]) / 2809
    return alpha


def verify_geometry_stress(scope, frame_groups):
    geometry = table(RESULT / 'native_geometry_stress.csv').set_index(['run', 'threshold'])
    trials = table(RESULT / 'mass_balance_stress.csv')
    summaries = table(RESULT / 'mass_balance_stress_per_run.csv').set_index('run')
    flow = table(OLD / 'results/massflow_per_run.csv').set_index('run')
    optical = table(OLD / 'frame_motion_evaluation/results/temporal_per_run.csv').set_index('run')
    inv = table(PHASE / 'results/native_phase_inventory.csv').set_index('run')
    model = json.loads(remember(OLD / 'results/validation/model_metadata.json').read_text())
    rho_l = PropsSI('D', 'T', 300.15, 'Q', 0, 'IsoButane')
    rho_g = PropsSI('D', 'T', 300.15, 'Q', 1, 'IsoButane')
    same([model['properties']['temperature_c'], model['properties']['rho_l'], model['properties']['rho_g']],
         [27, rho_l, rho_g], 'property_assumptions')
    assert len(geometry) == 78 and len(trials) == 702 and len(summaries) == 26
    assert geometry.index.is_unique and set(geometry.index.get_level_values('run')) == set(scope.index)
    assert not trials.accepted_as_phase_velocity.any() and not summaries.accepted_as_phase_velocity.any()
    for run in scope.index:
        path = remember(PHASE / inv.loc[run].packed_file)
        assert sha(path) == inv.loc[run].packed_sha256
        with np.load(path, allow_pickle=False) as archive:
            assert tuple(archive['mask_shape']) == (501, 80, 1016)
            assert tuple(archive['lumen_rows']) == (18, 71)
            assert np.array_equal(archive['frame_ids'], np.arange(501))
            for threshold, suffix in [(.6, '060'), (.7, '070'), (.8, '080')]:
                mask = np.unpackbits(archive['masks_'+suffix], axis=-1, bitorder='big', count=1016).astype(bool)
                assert not mask[:, :18].any() and not mask[:, 71:].any()
                values = radial_run_squares(mask)
                alpha = float(values.mean())
                g = geometry.loc[(run, threshold)]
                same([g.native_circle_alpha, g.native_projected_fraction, g.full_height_column_fraction,
                      g.alpha_frame_min, g.alpha_frame_max],
                     [alpha, mask[:,18:71].mean(), (mask[:,18:71].all(axis=1)).mean(),
                      values.mean(axis=1).min(), values.mean(axis=1).max()], 'native_geometry', atol=2e-14, rtol=0)
                probe = mask[:,18:71,496:520].mean(axis=(1,2))
                same(frame_groups[run]['occupancy_'+suffix].to_numpy(), probe, 'native_probe_occupancy', atol=2e-14, rtol=0)
                COUNTS['native_masks_reconstructed'] += mask.shape[0]
                COUNTS['radial_columns_reconstructed'] += values.size
        print(run, 'independent radial runs checked', flush=True)
    for row in trials.itertuples():
        g = geometry.loc[(row.run, row.threshold)]
        assert row.pattern == g.pattern
        assert row.massflow_statistic in ('mean', 'p05', 'p95')
        assert row.axial_scale_multiplier in (.9, 1., 1.1)
        alpha = float(g.native_circle_alpha)
        massflux = float(flow.loc[row.run, 'massflow_'+row.massflow_statistic+'_g_s'])*.001/(math.pi*1e-6/4)
        optical_value = float(optical.loc[row.run, 'median_signed_px_s'])
        scale = SCALE*row.axial_scale_multiplier
        ug = abs(optical_value)*scale
        # Reconstruct through phase mass fluxes, not the generator's direct expression.
        gas_massflux = rho_g*alpha*ug
        liquid_massflux = massflux-gas_massflux
        jl = liquid_massflux/rho_l
        ul = jl/(1-alpha) if alpha < 1 else math.nan
        jg = gas_massflux/rho_g
        same([row.G_kg_m2_s, row.optical_proxy_signed_px_s, row.assumed_axial_scale_m_px,
              row.substituted_ug_m_s, row.inferred_ul_m_s, row.implied_Jg_m_s, row.implied_Jl_m_s,
              row.rho_l_kg_m3, row.rho_g_kg_m3],
             [massflux, optical_value, scale, ug, ul, jg, jl, rho_l, rho_g], 'stress_arithmetic')
        same(row.conservation_residual_kg_m2_s, rho_g*row.implied_Jg_m_s+rho_l*row.implied_Jl_m_s-massflux,
             'stress_conservation', atol=1e-9, rtol=0)
        reasons = ['unverified_material_identity', 'conditional_native_geometry', 'condition_only_massflow',
                   'product_of_nonmatched_means', 'assumed_axial_scale_and_properties']
        if not optical.loc[row.run, 'base_coverage_and_saturation_pass']:
            reasons.append('low_optical_support_or_saturation')
        if alpha >= .98: reasons.append('small_liquid_area_denominator')
        if alpha == 1: reasons.append('singular_alpha')
        if math.isfinite(ul) and ul < 0: reasons.append('negative_liquid_inference')
        assert row.limitation_flags == ';'.join(reasons)
        COUNTS['stress_rejection_checks'] += 1
    assert not trials.duplicated(['run', 'threshold', 'massflow_statistic', 'axial_scale_multiplier']).any()
    assert (trials.groupby('run').size() == 27).all()
    threshold_only = trials[(trials.massflow_statistic == 'mean') & (trials.axial_scale_multiplier == 1)]
    for run, s in summaries.iterrows():
        group = threshold_only[threshold_only.run == run]
        values = group.inferred_ul_m_s.dropna()
        factor = values.max()/values.min() if len(values) and values.min() > 0 else math.nan
        same([s.alpha_min, s.alpha_max, s.candidate_liquid_min_m_s, s.candidate_liquid_max_m_s, s.threshold_max_min_factor],
             [group.native_circle_alpha.min(), group.native_circle_alpha.max(), values.min(), values.max(), factor],
             'stress_run_summary')
        assert bool(s.any_singular_alpha) == bool((group.native_circle_alpha == 1).any())
    meta = json.loads(remember(RESULT / 'mass_balance_stress_checks.json').read_text())
    baseline = threshold_only[threshold_only.threshold == .7]
    same([meta['baseline_candidate_liquid_min_m_s'], meta['baseline_candidate_liquid_max_m_s'],
          meta['baseline_alpha_min'], meta['baseline_alpha_max'], meta['max_threshold_liquid_factor'],
          meta['singular_rows'], meta['negative_liquid_rows']],
         [baseline.inferred_ul_m_s.min(), baseline.inferred_ul_m_s.max(), baseline.native_circle_alpha.min(),
          baseline.native_circle_alpha.max(), summaries.threshold_max_min_factor.max(),
          trials.inferred_ul_m_s.isna().sum(), (trials.inferred_ul_m_s < 0).sum()], 'stress_metadata')
    return dict(geometry_rows=78, trial_rows=702, accepted_material_pairs=0,
                baseline_max_inferred_liquid_m_s=float(baseline.inferred_ul_m_s.max()),
                interpretation='Rejected assumption stress test; arithmetic identity is not physical validation.')


def main():
    scope = verify_scope()
    frames = verify_passage(scope)
    stress = verify_geometry_stress(scope, frames)
    for relative, digest in INPUTS.items():
        assert sha(SRC2/relative) == digest, ('Input changed during verification', relative)
    result = dict(status='passed', scope_recordings=len(scope), passage_comparisons=4,
                  checks=dict(COUNTS), max_absolute_differences=dict(ERRORS), stress=stress,
                  input_sha256=INPUTS, verifier_sha256=sha(Path(__file__)),
                  independence='No generating module imported; independent temporal gates, math.fsum means, edge-paired radial runs, phase-flux mass reconstruction and CoolProp property call.',
                  limitations=['The camera clock, axial scale, phase identities and physical velocity accuracy are not independently validated.',
                               'Raw-frame fingerprints and the liquid-reference appearance computation are supported by separately recorded audits; the present run reconstructs occupancy from every saved mask.',
                               'Regional/annular outputs and final figures will be audited separately after their final support criteria are saved.'])
    (RESULT/'independent_verification.json').write_text(json.dumps(result, indent=2)+chr(10))
    print(json.dumps({k:v for k,v in result.items() if k != 'input_sha256'}, indent=2))


if __name__ == '__main__':
    main()
