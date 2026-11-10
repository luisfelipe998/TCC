#!/usr/bin/env python3
"""Plot the current reporting reference directly from saved passage means."""
from pathlib import Path
import csv
import math
import os

BASE = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR', str(BASE / 'intermediate/figure_mplcache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter


def read(name):
    with (BASE / 'results' / name).open() as stream:
        return list(csv.DictReader(stream))


def build():
    config, gate = 'consensus_reference_guard', 'lag2_confirmed'
    phase = {(r['run'], r['phase']): r for r in read('passage_phase_summary.csv')
             if r['config'] == config and r['motion_gate'] == gate}
    selected = [r for r in read('passage_ratios.csv')
                if r['config'] == config and r['motion_gate'] == gate]
    assert len(selected) == 26
    rows = []
    for r in selected:
        if r['eligible_descriptive_ratio'] != 'True':
            continue
        gas = phase[r['run'], 'gas_candidate']
        liquid = phase[r['run'], 'liquid_candidate']
        fps = float(gas['fps'])
        assert fps == float(liquid['fps'])
        vg, vl = float(gas['mean_px_s']), float(liquid['mean_px_s'])
        assert math.isclose(vg, float(r['gas_mean_px_frame']) * fps, rel_tol=1e-12)
        assert math.isclose(vl, float(r['liquid_mean_px_frame']) * fps, rel_tol=1e-12)
        ratio = float(r['descriptive_ratio'])
        assert math.isclose(ratio, vg / vl, rel_tol=1e-12)
        rows.append(dict(run=r['run'], pattern=r['pattern'], config=config, motion_gate=gate,
                         fps=fps, gas_passage_mean_px_s=vg, liquid_gap_mean_px_s=vl,
                         ratio_of_means=ratio, episode_equal_ratio=float(r['episode_equal_ratio']),
                         gas_pairs=int(r['gas_pairs']), liquid_pairs=int(r['liquid_pairs']),
                         gas_episodes=int(r['gas_episodes']), liquid_episodes=int(r['liquid_episodes'])))
    assert len(rows) == 4
    with (BASE / 'results/conditional_passage_velocities.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, (vel, ratio_ax) = plt.subplots(1, 2, figsize=(13, 6.2), sharey=True,
                                      gridspec_kw={'width_ratios': [1.55, 1]})
    for i, r in enumerate(rows):
        for key, dy, color, marker, label in [
            ('gas_passage_mean_px_s', -.13, '#2769a7', 'o', 'During gas-like passages'),
            ('liquid_gap_mean_px_s', .13, '#218563', 's', 'During liquid-looking gaps')]:
            value = r[key]
            vel.scatter(value, i + dy, color=color, marker=marker, s=55,
                        label=label if i == 0 else None, zorder=3)
            vel.annotate(f'{value:,.0f}', (value, i + dy), xytext=(7, 0),
                         textcoords='offset points', va='center', fontsize=9, color=color)
        for key, dy, color, marker, label in [
            ('ratio_of_means', -.13, '#4b4b9b', 'o', 'Ratio of pair-weighted means'),
            ('episode_equal_ratio', .13, '#b16a24', 'D', 'Equal weight per episode')]:
            value = r[key]
            ratio_ax.scatter(value, i + dy, color=color, marker=marker, s=45,
                             label=label if i == 0 else None, zorder=3)
            ratio_ax.annotate(f'{value:.3f}', (value, i + dy), xytext=(7, 0),
                              textcoords='offset points', va='center', fontsize=9, color=color)
    vel.set_yticks(range(len(rows)), [r['run'] + chr(10) + r['pattern'] for r in rows])
    vel.set_ylim(len(rows) - .5, -.6)
    vel.set_xlim(-118000, -65000)
    vel.set_xticks([-115000, -105000, -95000, -85000, -75000])
    vel.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    vel.set_xlabel('Signed mean image motion (px/s)')
    vel.set_title('Separate passage averages')
    ratio_ax.axvline(1, color='.6', linestyle='--', linewidth=1)
    ratio_ax.set_xlim(.86, 1.05)
    ratio_ax.set_xticks([.88, .92, .96, 1, 1.04])
    ratio_ax.set_xlabel('Gas-passage mean / liquid-gap mean (dimensionless)')
    ratio_ax.set_title('Ratio and averaging sensitivity')
    for ax in (vel, ratio_ax):
        ax.grid(axis='x', alpha=.2)
        ax.legend(loc='upper center', bbox_to_anchor=(.5, -.18), frameon=False, fontsize=9)
    fig.suptitle('Current conditional image-motion results: 4 of 26 recordings' + chr(10)
                 + 'Threshold consensus + liquid appearance guard + lag-2 consistency screen', fontsize=13)
    fig.text(.5, .025, 'Negative = leftward. These are image-structure descriptors; separate material phase velocities remain unverified.',
             ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .06, 1, .91))
    for extension in ('png', 'svg'):
        fig.savefig(BASE / f'figures/conditional_passage_velocities.{extension}', dpi=180)
    plt.close(fig)
    print('Saved four current passage rows and the px/s / dimensionless PNG and SVG.')


if __name__ == '__main__':
    build()
