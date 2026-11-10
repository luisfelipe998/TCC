#!/usr/bin/env python3
"""Reaggregate saved native-frame passage evidence without new phase claims."""
from pathlib import Path
import csv
import hashlib
import json
import math
import re
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
SRC2 = BASE.parent
OLD = SRC2 / 'v15_gas_liquid_velocity'
PHASE = OLD / 'phase_conditioned_motion'
DATA = SRC2 / 'data_processed'
CONFIG, GATE = 'consensus_reference_guard', 'lag2_confirmed'
SCALE = .001 / 53


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(name, rows):
    pd.DataFrame(rows).to_csv(BASE / 'results' / name, index=False)


def run():
    prior_inventory = pd.read_csv(OLD/'results/scenario_inventory.csv').set_index('run')
    inventory = []
    for folder in sorted(DATA.glob('m??q??_???')):
        markers = list(folder.glob('__*.txt'))
        assert len(markers) == 1
        pattern = markers[0].stem[2:].replace('_', ' ')
        files = sorted(folder.glob('ImgA*.jpg'))
        ids = [int(re.search(r'(\d+)\.jpg$', p.name)[1]) for p in files]
        assert ids == list(range(501))
        xsv, mii = (folder/'_config.xsv').read_text(), (folder/'_config.mii').read_text()
        fps = float(re.search(r'^Rate=([0-9.]+)', xsv, re.M)[1])
        assert fps == float(re.search(r'Film speed 1\s*:([0-9.]+)', mii)[1])
        assert int(re.search(r'^Images=([0-9]+)', xsv, re.M)[1]) == len(files)
        old = prior_inventory.loc[folder.name]
        assert pattern == old['pattern'] and fps == old.native_fps
        assert sha(folder/'_config.xsv') == old.config_sha256
        assert sha(folder/'_config.mii') == old.mii_sha256
        inventory.append(dict(run=folder.name, pattern=pattern, included=bool(old.included),
                              exclusion_reason=old.exclusion_reason if pd.notna(old.exclusion_reason) else '',
                              fps=fps, frames=len(files), first_to_last_s=500/fps,
                              config_sha256=old.config_sha256, mii_sha256=old.mii_sha256,
                              conditional_axial_scale_m_px=SCALE,
                              axial_scale_status='assumed D=1mm over 53 radial pixel cells; optical mapping uncalibrated'))
    assert len(inventory) == 31 and sum(r['included'] for r in inventory) == 26
    write('scenario_inventory.csv', inventory)
    pairs = pd.read_csv(PHASE/'results/passage_pairs.csv')
    pairs = pairs[(pairs.config == CONFIG) & (pairs.motion_gate == GATE)].copy()
    assert len(pairs) == 13000
    prior = pd.read_csv(PHASE/'results/passage_ratios.csv')
    prior = prior[(prior.config == CONFIG) & (prior.motion_gate == GATE)].set_index('run')
    temporal = pd.read_csv(OLD/'frame_motion_evaluation/results/temporal_per_run.csv').set_index('run')
    rows = []
    for inv in inventory:
        if not inv['included']:
            continue
        run, fps = inv['run'], inv['fps']
        selected = pairs[(pairs.run == run) & pairs.accepted]
        np.testing.assert_allclose(pairs.loc[pairs.run == run, 'velocity_px_s'],
                                   pairs.loc[pairs.run == run, 'displacement_px_frame']*fps,
                                   rtol=1e-12)
        row = dict(run=run, pattern=inv['pattern'], fps=fps, config=CONFIG, motion_gate=GATE)
        reasons = []
        for phase, short in [('gas_candidate', 'gas'), ('liquid_candidate', 'liquid')]:
            values = selected[selected.phase == phase]
            episode_means = values.groupby('episode_id').velocity_px_s.mean()
            n, ne = len(values), len(episode_means)
            mean = float(values.velocity_px_s.mean())
            row.update({f'{short}_pairs': n, f'{short}_episodes': ne,
                        f'{short}_mean_signed_px_s': mean,
                        f'{short}_episode_equal_mean_signed_px_s': float(episode_means.mean()),
                        f'{short}_positive_fraction': float((values.velocity_px_s > 0).mean()) if n else math.nan,
                        f'{short}_q25_signed_px_s': float(values.velocity_px_s.quantile(.25)),
                        f'{short}_q75_signed_px_s': float(values.velocity_px_s.quantile(.75))})
            if n < 20:
                reasons.append(short+'_too_few_pairs')
            if ne < 3:
                reasons.append(short+'_too_few_episodes')
        g, l = row['gas_mean_signed_px_s'], row['liquid_mean_signed_px_s']
        saturated = bool(temporal.loc[run, 'high_saturation'])
        if not (math.isfinite(g) and math.isfinite(l)):
            reasons.append('missing_phase_mean')
        elif abs(l)/fps < 1:
            reasons.append('near_zero_liquid_mean')
        elif g*l <= 0:
            reasons.append('opposed_or_zero_mean_directions')
        if saturated:
            reasons.append('high_saturation')
        eligible = not reasons
        row.update(eligible_passage_comparison=eligible, withheld_reasons=';'.join(reasons),
                   high_saturation=saturated,
                   ratio_of_signed_means=g/l if eligible else math.nan,
                   episode_equal_ratio=row['gas_episode_equal_mean_signed_px_s']/row['liquid_episode_equal_mean_signed_px_s'] if eligible else math.nan,
                   gas_plot_magnitude_px_s=abs(g) if eligible else math.nan,
                   liquid_plot_magnitude_px_s=abs(l) if eligible else math.nan,
                   gas_conditional_m_s=g*SCALE if eligible else math.nan,
                   liquid_conditional_m_s=l*SCALE if eligible else math.nan,
                   scale_m_px=SCALE,
                   observable='global image motion conditioned on local phase appearance; not material or superficial velocity')
        assert eligible == bool(prior.loc[run, 'eligible_descriptive_ratio']), run
        if eligible:
            for field, ref in [('ratio_of_signed_means','descriptive_ratio'), ('episode_equal_ratio','episode_equal_ratio')]:
                assert math.isclose(row[field], prior.loc[run,ref], rel_tol=1e-12), (run,field)
            assert row['gas_positive_fraction'] == row['liquid_positive_fraction'] == 0
        rows.append(row)
    assert len(rows) == 26 and sum(r['eligible_passage_comparison'] for r in rows) == 4
    write('passage_per_run.csv', rows)
    pairs.to_csv(BASE/'results/passage_pairs_reference.csv', index=False)
    bubble = pd.read_csv(PHASE/'results/closed_contour_bubble_summary.csv')
    bubble.insert(0, 'run', 'm07q01_001')
    bubble['mean_signed_px_s'] = bubble.mean_px_frame*7000
    bubble['conditional_mean_m_s'] = bubble.mean_signed_px_s*SCALE
    bubble['observable'] = 'selected candidate bubble-object translation; no independent liquid speed'
    bubble.to_csv(BASE/'results/bubble_object_reference.csv', index=False)
    sources = [OLD/'results/scenario_inventory.csv', PHASE/'results/passage_pairs.csv',
               PHASE/'results/passage_ratios.csv', PHASE/'results/closed_contour_bubble_summary.csv',
               OLD/'frame_motion_evaluation/results/temporal_per_run.csv']
    metadata = dict(status='passed', selected_recordings=26, excluded_recordings=5,
                    pair_rows_reaggregated=13000, primary_comparisons=4,
                    config=CONFIG, motion_gate=GATE, velocity_units='px/s',
                    conditional_scale_m_px=SCALE,
                    averaging='signed arithmetic mean of accepted adjacent comparisons, separately per appearance; absolute value only for plotting',
                    source_mode='Reaggregated previous native-frame image evidence; segmentation/tracking was not refitted.',
                    inputs={str(p.relative_to(SRC2)):sha(p) for p in sources},
                    script_sha256=sha(Path(__file__)),
                    limitations=['Phase appearance is inferred, not independently verified.',
                                 'Motion template is wider than appearance probe.',
                                 'Four selected pairs do not describe all five flow patterns.',
                                 'Record duration is 71–125 ms; overlapping comparisons are not independent repeats.'])
    (BASE/'results/passage_reaggregation_checks.json').write_text(json.dumps(metadata, indent=2)+chr(10))
    print(json.dumps({k:metadata[k] for k in ['status','selected_recordings','pair_rows_reaggregated','primary_comparisons']}, indent=2))


if __name__ == '__main__':
    run()
