#!/usr/bin/env python3
"""Conditional thermodynamic velocity comparisons; no ground-truth accuracy claim."""
from __future__ import annotations

import configparser
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import sys

import CoolProp
from CoolProp.CoolProp import PropsSI
import numpy as np
import pandas as pd
from scipy.optimize import brentq

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[1]
DATA = ROOT / 'src2/data_processed'
os.environ.setdefault('MPLCONFIGDIR', str(BASE / 'intermediate/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PATTERNS = ('dispersed bubbles', 'plug', 'slug', 'churn', 'wavy annular')
COLORS = dict(zip(PATTERNS, ('#b58900', '#4269d0', '#d55e00', '#009e73', '#9c6ade')))
MARKERS = dict(zip(PATTERNS, ('P', 'o', 's', '^', 'D')))
G_MAP = {3: 380., 4: 510., 5: 640., 6: 760., 7: 890.}
D, L, GRAVITY = 0.001, 0.265, 9.81
MODELS = ('homogeneous', 'rouhani_axelsson', 'xu_fang_literal_2021', 'xu_fang_fluids')


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def properties(temperature_c=27.):
    t = temperature_c + 273.15
    return dict(temperature_c=float(temperature_c),
                rho_l=PropsSI('D', 'T', t, 'Q', 0, 'R600a'),
                rho_g=PropsSI('D', 'T', t, 'Q', 1, 'R600a'),
                h_l=PropsSI('H', 'T', t, 'Q', 0, 'R600a'),
                h_fg=PropsSI('H', 'T', t, 'Q', 1, 'R600a') - PropsSI('H', 'T', t, 'Q', 0, 'R600a'),
                pressure_pa=PropsSI('P', 'T', t, 'Q', 0, 'R600a'),
                sigma=PropsSI('I', 'T', t, 'Q', 0, 'R600a'))


def inventory():
    rows = []
    for folder in sorted(DATA.glob('m??q??_???')):
        m = re.fullmatch(r'm(\d{2})q(\d{2})_(\d{3})', folder.name)
        markers = sorted(folder.glob('__*.txt'))
        if len(markers) != 1:
            raise ValueError(f'{folder}: expected exactly one pattern marker')
        pattern = markers[0].stem[2:].replace('_', ' ')
        reason = ''
        if folder.name == 'm05q05_002':
            reason = 'explicit_user_exclusion'
        elif pattern == 'smooth annular':
            reason = 'smooth_annular'
        elif '-' in pattern:
            reason = 'double_classification'
        elif pattern not in PATTERNS:
            raise ValueError(f'Unrecognized pattern: {pattern}')
        config = configparser.ConfigParser(strict=False)
        config.read(folder / '_config.xsv')
        acq = config['Acquisition']
        fps, n = float(acq['Rate']), int(acq['Images'])
        mii = (folder / '_config.mii').read_text()
        assert fps == float(re.search(r'Film speed 1\s*:(\S+)', mii)[1])
        assert n == int(re.search(r'Number of images 1\s*:(\S+)', mii)[1])
        assert len(list(folder.glob('ImgA*.jpg'))) == n
        rows.append(dict(run=folder.name, pattern=pattern, included=not reason,
                         exclusion_reason=reason, G_kg_m2_s=G_MAP[int(m[1])],
                         q_kw_m2=float(m[2]), native_fps=fps, frames=n,
                         duration_s=(n - 1) / fps,
                         config_sha256=sha256(folder / '_config.xsv'),
                         mii_sha256=sha256(folder / '_config.mii'),
                         marker_file=str(markers[0].relative_to(ROOT)),
                         marker_sha256=sha256(markers[0])))
    result = pd.DataFrame(rows)
    selected = result[result.included]
    assert len(result) == 31 and len(selected) == 26
    assert selected.pattern.value_counts().to_dict() == {
        'dispersed bubbles': 1, 'plug': 3, 'slug': 3, 'churn': 10, 'wavy annular': 9}
    return result


def alpha_model(x, G, p, model):
    if not 0 <= x <= 1:
        return np.nan
    if x in (0., 1.):
        return x
    vg, vl = x / p['rho_g'], (1 - x) / p['rho_l']
    ah = vg / (vg + vl)
    if model == 'homogeneous':
        return ah
    if model == 'rouhani_axelsson':
        c0 = 1 + 0.12 * (1 - x)
        drift = 1.18 * (1 - x) * (GRAVITY * p['sigma'] * (p['rho_l'] - p['rho_g']))**0.25 / (G * np.sqrt(p['rho_l']))
        return vg / (c0 * (vg + vl) + drift)
    if model in ('xu_fang_literal_2021', 'xu_fang_fluids'):
        # Literal de Oliveira 2021 Eq.5 differs from the fluids implementation.
        fr_lo = G**2 / (p['rho_l']**2 * GRAVITY * D)
        exponent = -2 if model == 'xu_fang_literal_2021' else -.2
        return 1 / (1 + (1 + 2 * fr_lo**exponent * ah**3.5) * vl / vg)
    raise ValueError(model)


def velocities(x, G, p, model):
    if not 0 <= x <= 1:
        return dict(quality_status='outside_two_phase_equilibrium_domain',
                    Jg=np.nan, Jl=np.nan, alpha=np.nan, ug=np.nan, ul=np.nan, slip=np.nan)
    jg, jl = G * x / p['rho_g'], G * (1 - x) / p['rho_l']
    a = alpha_model(x, G, p, model)
    ug = jg / a if a > 0 else np.nan
    ul = jl / (1 - a) if a < 1 else np.nan
    return dict(quality_status='two_phase' if 0 < x < 1 else 'single_phase_endpoint',
                Jg=jg, Jl=jl, alpha=a, ug=ug, ul=ul, slip=ug / ul if 0 < a < 1 else np.nan)


def experiments(selected):
    # One-at-a-time stress scenarios, not measured uncertainty or calibrated bounds.
    scenarios = [dict(scenario='saturated_inlet', temperature_c=27., subcool_k=0., xin=0., eta=1.)]
    scenarios += [dict(scenario=f'subcool_{k}K', temperature_c=27., subcool_k=float(k), xin=0., eta=1.) for k in (1, 2, 5, 10)]
    scenarios += [dict(scenario=f'inlet_quality_{x:.2f}', temperature_c=27., subcool_k=0., xin=x, eta=1.) for x in (0.01, 0.05)]
    scenarios += [dict(scenario=f'Tsat_{t}C', temperature_c=float(t), subcool_k=0., xin=0., eta=1.) for t in (25, 29)]
    scenarios += [dict(scenario='net_heat_90pct', temperature_c=27., subcool_k=0., xin=0., eta=0.9)]
    scenarios += [dict(scenario='recorded_flow_' + stat, temperature_c=27., subcool_k=0., xin=0., eta=1.,
                       flow_column=column) for stat, column in
                  [('mean', 'G_mean_kg_m2_s'), ('median', 'G_median_kg_m2_s'), ('trimmed10', 'G_trimmed10_kg_m2_s')]]
    flows = pd.read_csv(BASE / 'results/massflow_per_run.csv').set_index('run')
    assert set(flows.index) == set(selected.run)
    out = []
    for spec in scenarios:
        p = properties(spec['temperature_c'])
        xin = spec['xin']
        if spec['subcool_k']:
            hin = PropsSI('H', 'T', spec['temperature_c'] + 273.15 - spec['subcool_k'],
                          'P', p['pressure_pa'], 'R600a')
            xin = (hin - p['h_l']) / p['h_fg']
        for r in selected.to_dict('records'):
            G, q = r['G_kg_m2_s'], r['q_kw_m2'] * 1000
            if 'flow_column' in spec:
                G = float(flows.loc[r['run'], spec['flow_column']])
                assert flows.loc[r['run'], 'pattern'] == r['pattern']
            dx = spec['eta'] * 4 * q * L / (G * D * p['h_fg'])
            x = xin + dx
            for model in MODELS:
                out.append(dict(run=r['run'], pattern=r['pattern'], model=model,
                                scenario=spec['scenario'], G_kg_m2_s=G, q_kw_m2=q/1000,
                                G_basis=spec.get('flow_column', 'nominal_literature_lookup'),
                                temperature_c=spec['temperature_c'], net_heat_fraction=spec['eta'],
                                inlet_equilibrium_quality=xin, quality_increment=dx, x_eq=x,
                                rho_l=p['rho_l'], rho_g=p['rho_g'],
                                Fr_lo=G**2 / (p['rho_l']**2 * GRAVITY * D),
                                XF_secondary_Fr_domain=(0.02 <= G**2 / (p['rho_l']**2 * GRAVITY * D) <= 145),
                                **velocities(x, G, p, model)))
    return pd.DataFrame(out), scenarios


def image_diagnostics(baseline):
    src = ROOT / 'src2/v14_new_voidfraction/results/void_fraction/panorama_void_per_run.csv'
    snapshot = BASE / 'intermediate/v14_void_input_snapshot.csv'
    # Freeze this dependency once; --refresh-inputs explicitly updates the snapshot.
    if not snapshot.exists() or '--refresh-inputs' in sys.argv:
        shutil.copyfile(src, snapshot)
    v14 = pd.read_csv(snapshot).set_index('run')
    output = []
    p = properties()
    for r in baseline[baseline.model == 'rouhani_axelsson'].to_dict('records'):
        v = v14.loc[r['run']]
        assert v['label'] == r['pattern']
        for name, key in [('D60', 'conditional_circle_alpha_spatial_mean'),
                          ('D56', 'conditional_alpha_D56'), ('D64', 'conditional_alpha_D64')]:
            a = float(v[key])
            supported = v['validity'] != 'unsupported_saturation' and 0 < a < 1
            item = dict(run=r['run'], pattern=r['pattern'], variant=name, alpha_image=a,
                        image_validity=v['validity'],
                        diagnostic_status='unvalidated_spatial_alpha_hybrid' if supported else 'unsupported_saturation_or_alpha',
                        ug=np.nan, ul=np.nan, slip=np.nan,
                        RA_implied_x=np.nan, RA_implied_x_in=np.nan)
            if supported:
                ug, ul = r['Jg'] / a, r['Jl'] / (1-a)
                xi = brentq(lambda x: alpha_model(x, r['G_kg_m2_s'], p, 'rouhani_axelsson') - a, 1e-12, 1-1e-12)
                item.update(ug=ug, ul=ul, slip=ug/ul, RA_implied_x=xi, RA_implied_x_in=xi-r['x_eq'])
            output.append(item)
    return pd.DataFrame(output), dict(source=str(src.relative_to(ROOT)),
                                      snapshot=str(snapshot.relative_to(BASE)),
                                      snapshot_sha256=sha256(snapshot))


def validation(all_rows, baseline, selected):
    p = properties()
    mass = all_rows.rho_g * all_rows.Jg + all_rows.rho_l * all_rows.Jl
    valid = all_rows.quality_status == 'two_phase'
    max_mass = float(np.max(np.abs(mass[valid] - all_rows.loc[valid, 'G_kg_m2_s'])))
    assert max_mass < 1e-9
    for model in MODELS:
        assert np.isnan(velocities(-.01, 640, p, model)['ug'])
        assert np.isnan(velocities(1.01, 640, p, model)['ul'])
        assert np.isnan(velocities(0., 640, p, model)['ug'])
        assert np.isnan(velocities(1., 640, p, model)['ul'])
        assert alpha_model(0., 640, p, model) == 0
        assert alpha_model(1., 640, p, model) == 1
    homogeneous = all_rows[(all_rows.model == 'homogeneous') & valid]
    assert np.allclose(homogeneous.ug, homogeneous.ul, rtol=1e-12)
    assert np.allclose(all_rows.loc[valid, 'alpha'] * all_rows.loc[valid, 'ug'], all_rows.loc[valid, 'Jg'])
    assert np.allclose((1-all_rows.loc[valid, 'alpha']) * all_rows.loc[valid, 'ul'], all_rows.loc[valid, 'Jl'])
    assert all_rows[valid].groupby(['scenario', 'run']).Jg.nunique().eq(1).all()
    assert all_rows[valid].groupby(['scenario', 'run']).Jl.nunique().eq(1).all()
    assert np.allclose(baseline.groupby('q_kw_m2').Jg.max(), baseline.groupby('q_kw_m2').Jg.min())
    assert set(all_rows.run) == set(selected.run)
    assert len(baseline) == len(selected) * len(MODELS)
    assert not (all_rows.run == 'm05q05_002').any()
    return dict(status='passed', scope_count=len(selected), closure_rows=len(baseline),
                experiment_rows=len(all_rows), max_mass_balance_residual_kg_m2_s=max_mass,
                checks=['exact scope', 'native config agreement', 'frame count',
                        'mass conservation', 'phase flux identities', 'homogeneous equal speeds',
                        'superficial invariance to closure', 'same-q cancellation under zero inlet',
                        'invalid quality not clipped', 'absent phase velocity undefined'])


def plots(baseline, sensitivity):
    fig, ax = plt.subplots(figsize=(8.5, 6.2))
    b = baseline[baseline.model == 'rouhani_axelsson']
    for pat in PATTERNS:
        s = b[b.pattern == pat]
        ax.scatter(s.Jg, s.Jl, c=COLORS[pat], marker=MARKERS[pat], label=f'{pat} (n={len(s)})', s=65)
    ax.set(xscale='log', yscale='log', xlabel='Gas superficial velocity Jg [m/s]',
           ylabel='Liquid superficial velocity Jl [m/s]',
           title='26 selected scenarios: conditional superficial velocities\nSaturated inlet, nominal G and heat flux; no slip closure required')
    ax.grid(alpha=.2, which='both')
    ax.legend(fontsize=9)
    fig.tight_layout()
    for ext in ('png', 'svg'):
        fig.savefig(BASE / f'figures/superficial_velocity.{ext}', dpi=170)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(11, 10), sharex=True, sharey=True)
    labels = ('Homogeneous: imposed equal phase speeds', 'Rouhani–Axelsson: v12.1 closure',
              'Xu–Fang: literal 2021 exponent −2', 'Xu–Fang: fluids exponent −0.2')
    for ax, model, label in zip(axes.flat, MODELS, labels):
        for pat in PATTERNS:
            s = baseline[(baseline.model == model) & (baseline.pattern == pat)]
            ax.scatter(s.ug, s.ul, c=COLORS[pat], marker=MARKERS[pat], s=45, label=pat)
        ax.plot([.5, 50], [.5, 50], '--', color='.5', lw=1, label='ug = ul')
        ax.set(xscale='log', yscale='log', xlabel='Gas mean phase velocity ug [m/s]', title=label,
               ylabel='Liquid mean phase velocity ul [m/s]')
        ax.grid(alpha=.2, which='both')
    axes[0,0].legend(fontsize=8, loc='upper left')
    fig.suptitle('Model predictions at x_in = 0; differences are model spread, not measured error')
    fig.tight_layout()
    for ext in ('png', 'svg'):
        fig.savefig(BASE / f'figures/phase_model_comparison.{ext}', dpi=170)
    plt.close(fig)
    shown = ['subcool_1K', 'subcool_2K', 'subcool_5K', 'subcool_10K', 'inlet_quality_0.01', 'inlet_quality_0.05', 'Tsat_25C', 'Tsat_29C', 'net_heat_90pct']
    pivot = sensitivity.pivot(index='run', columns='scenario', values='Jg_change_pct')[shown]
    fig, ax = plt.subplots(figsize=(12, 9))
    values = pivot.to_numpy()
    im = ax.imshow(np.ma.masked_invalid(values), cmap='RdBu_r', vmin=-100, vmax=100, aspect='auto')
    for i in range(len(pivot)):
        for j in range(len(shown)):
            v = values[i,j]
            ax.text(j, i, 'n/a' if not np.isfinite(v) else f'{v:.0f}%', ha='center', va='center', fontsize=7, color='black' if not np.isfinite(v) or abs(v)<55 else 'white')
    ax.set_xticks(range(len(shown)), [x.replace('_',' ') for x in shown], rotation=35, ha='right')
    ax.set_yticks(range(len(pivot)), pivot.index)
    ax.set_title('Gas superficial velocity sensitivity relative to saturated inlet\nOne-at-a-time assumed scenarios; n/a = equilibrium quality outside [0, 1]')
    fig.colorbar(im, ax=ax, label='Change [%]; color saturated beyond ±100%, labels show actual change')
    fig.tight_layout()
    fig.savefig(BASE / 'figures/inlet_sensitivity.png', dpi=170)
    plt.close(fig)


def flow_comparison_plot(baseline, recorded):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.6))
    for pat in PATTERNS:
        a = baseline[(baseline.model == 'rouhani_axelsson') & (baseline.pattern == pat)].set_index('run')
        b = recorded[(recorded.model == 'rouhani_axelsson') & (recorded.pattern == pat)].set_index('run').loc[a.index]
        for ax, x, y in [(axes[0], 'Jg', 'Jl'), (axes[1], 'ug', 'ul')]:
            for run in a.index:
                ax.plot([a.loc[run,x], b.loc[run,x]], [a.loc[run,y], b.loc[run,y]], color=COLORS[pat], alpha=.5, lw=.9)
            ax.scatter(a[x], a[y], marker=MARKERS[pat], edgecolors=COLORS[pat], facecolors='none', s=50)
            ax.scatter(b[x], b[y], marker=MARKERS[pat], c=COLORS[pat], s=35, label=pat)
    axes[0].set(xlabel='Gas superficial velocity Jg [m/s]', ylabel='Liquid superficial velocity Jl [m/s]', title='Superficial velocities; no void-fraction model')
    axes[1].set(xlabel='Gas mean phase velocity ug [m/s]', ylabel='Liquid mean phase velocity ul [m/s]', title='Rouhani–Axelsson conditional phase means')
    for ax in axes:
        ax.set(xscale='log', yscale='log')
        ax.grid(alpha=.2, which='both')
    axes[0].legend(fontsize=8)
    fig.suptitle('26 scenarios, saturated inlet: hollow = nominal G; filled = recorded mean flow\nCondition-level matching only; camera synchronization and inlet state unverified', fontsize=11)
    fig.tight_layout()
    for ext in ('png', 'svg'):
        fig.savefig(BASE / f'figures/recorded_flow_comparison.{ext}', dpi=170)
    plt.close(fig)


def main():
    for folder in ('results', 'results/validation', 'figures', 'intermediate'):
        (BASE / folder).mkdir(parents=True, exist_ok=True)
    inv = inventory()
    selected = inv[inv.included]
    all_rows, scenarios = experiments(selected)
    baseline = all_rows[all_rows.scenario == 'saturated_inlet'].copy()
    image, image_manifest = image_diagnostics(baseline)
    check = validation(all_rows, baseline, selected)
    inv.to_csv(BASE / 'results/scenario_inventory.csv', index=False)
    all_rows.to_csv(BASE / 'results/thermodynamic_sensitivity.csv', index=False)
    baseline.to_csv(BASE / 'results/conditional_velocities_per_run.csv', index=False)
    summary = baseline.groupby(['model', 'pattern'], sort=False)[['x_eq','Jg','Jl','alpha','ug','ul','slip']].agg(['min', 'max', 'mean', 'median'])
    summary.columns = ['_'.join(c) for c in summary.columns]
    summary.insert(0, 'n', baseline.groupby(['model','pattern'], sort=False).size())
    summary.reset_index().to_csv(BASE / 'results/conditional_velocities_by_pattern.csv', index=False)
    recorded = all_rows[all_rows.scenario == 'recorded_flow_mean']
    recorded.to_csv(BASE / 'results/recorded_flow_conditional_velocities_per_run.csv', index=False)
    rec_summary = recorded.groupby(['model', 'pattern'], sort=False)[['G_kg_m2_s', 'x_eq', 'Jg', 'Jl', 'alpha', 'ug', 'ul', 'slip']].agg(['min', 'max', 'mean', 'median'])
    rec_summary.columns = ['_'.join(c) for c in rec_summary.columns]
    rec_summary.insert(0, 'n', recorded.groupby(['model','pattern'], sort=False).size())
    rec_summary.reset_index().to_csv(BASE / 'results/recorded_flow_conditional_velocities_by_pattern.csv', index=False)
    image.to_csv(BASE / 'results/image_alpha_diagnostics.csv', index=False)
    b = baseline[baseline.model == 'rouhani_axelsson'].set_index('run')
    s = all_rows[all_rows.model == 'rouhani_axelsson'].copy()
    for key in ('Jg', 'Jl', 'ug', 'ul'):
        s[key+'_change_pct'] = (s[key].to_numpy() / s.run.map(b[key]).to_numpy() - 1) * 100
    s.to_csv(BASE / 'results/sensitivity_relative_to_baseline.csv', index=False)
    closure = baseline.groupby('run').agg(ug_min=('ug','min'),ug_max=('ug','max'),ul_min=('ul','min'),ul_max=('ul','max'),alpha_min=('alpha','min'),alpha_max=('alpha','max'))
    closure['ug_max_to_min'] = closure.ug_max / closure.ug_min
    closure['ul_max_to_min'] = closure.ul_max / closure.ul_min
    closure.to_csv(BASE / 'results/closure_spread.csv')
    metadata = dict(properties=properties(), CoolProp_version=CoolProp.__version__,
                    python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
                    diameter_m=D, heated_length_m=L, G_map=G_MAP, gravity_m_s2=GRAVITY,
                    baseline_temperature_and_G_status='saturated_inlet baseline uses nominal literature/setpoint values; recorded_flow branches use massflow_input',
                    scopes=dict(source=len(inv), selected=len(selected), patterns=dict(Counter(selected.pattern))),
                    scenarios=scenarios, sensitivity_status='assumed scenarios, not confidence intervals',
                    image_input=image_manifest,
                    massflow_input=dict(path='results/massflow_per_run.csv', sha256=sha256(BASE / 'results/massflow_per_run.csv'),
                                        status='owner-documented g/s; original blocks verified; same condition only, camera synchronization unproven'),
                    xu_fang_Fr_definition='G^2/(rho_l^2*g*D); conventional closure',
                    xu_fang_exponents={'literal_2021': -2., 'fluids': -.2},
                    primary_velocity_status='conditional model prediction; no measured truth available',
                    script_sha256=sha256(__file__))
    (BASE / 'results/validation/model_checks.json').write_text(json.dumps(check, indent=2) + '\n')
    (BASE / 'results/validation/model_metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    plots(baseline, s)
    flow_comparison_plot(baseline, recorded)
    print(json.dumps(check, indent=2))
    print('RA baseline pattern medians:')
    print(baseline[baseline.model == 'rouhani_axelsson'].groupby('pattern')[['Jg','Jl','ug','ul','slip']].median().to_string())
    print('Outside two-phase domain by scenario:')
    print(s.groupby('scenario').quality_status.apply(lambda x: (x != 'two_phase').sum()).to_string())


if __name__ == '__main__':
    main()
