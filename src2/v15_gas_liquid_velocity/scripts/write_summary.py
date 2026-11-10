#!/usr/bin/env python3
"""Build historical model findings and publish the saved current pixel report."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from publish_current_summary import publish

BASE = Path(__file__).resolve().parents[1]
PATTERNS = ('dispersed bubbles', 'plug', 'slug', 'churn', 'wavy annular')


def main():
    results = BASE / 'results'
    inv = pd.read_csv(results / 'scenario_inventory.csv')
    nominal = pd.read_csv(results / 'conditional_velocities_per_run.csv')
    recorded = pd.read_csv(results / 'recorded_flow_conditional_velocities_per_run.csv')
    sensitivity = pd.read_csv(results / 'sensitivity_relative_to_baseline.csv')
    flow = pd.read_csv(results / 'massflow_per_run.csv')
    optical = pd.read_csv(results / 'optical_preferred_candidate_per_run.csv')
    all_optical = pd.read_csv(results / 'optical_full_method_summary.csv')
    closure = pd.read_csv(results / 'closure_spread.csv')
    hybrid = pd.read_csv(results / 'image_alpha_diagnostics.csv')
    model_checks = json.loads((results / 'validation/model_checks.json').read_text())
    optical_meta = json.loads((results / 'optical_metadata.json').read_text())
    ra = recorded[recorded.model == 'rouhani_axelsson']
    nominal_ra = nominal[nominal.model == 'rouhani_axelsson']
    scoped = inv[inv.included]
    n = len(scoped)
    acceptable = optical[~optical.low_coverage & ~optical.high_saturation]
    metrics = dict(selected_runs=n, pattern_counts=scoped.pattern.value_counts().to_dict(),
                   massflux_change_pct_min=float(flow.G_mean_delta_pct.min()),
                   massflux_change_pct_max=float(flow.G_mean_delta_pct.max()),
                   maximum_nominal_liquid_model_factor=float(closure.ul_max_to_min.max()),
                   preferred_optical_passing_pairs=int(optical.heuristic_pass_pairs.sum()),
                   optical_pairs=int(optical.pairs.sum()),
                   preferred_optical_supported_runs=len(acceptable),
                   preferred_optical_supported_annular_runs=int((acceptable.pattern == 'wavy annular').sum()),
                   thermodynamic_experiment_rows=model_checks['experiment_rows'],
                   nominal_XF_Fr_extrapolated_runs=int((~nominal_ra.XF_secondary_Fr_domain).sum()),
                   recorded_XF_Fr_extrapolated_runs=int((~ra.XF_secondary_Fr_domain).sum()))
    subcool = sensitivity[sensitivity.scenario == 'subcool_1K']
    metrics['subcool_1K_out_of_domain_runs'] = int((subcool.quality_status != 'two_phase').sum())
    lines = [
        '# Historical conditional model summary', '',
        '**This is the earlier thermodynamic/model investigation.** The [current executive summary](EXECUTIVE_SUMMARY.md) reports the later phase-passage image-motion results in px/s and dimensionless ratios. The model predictions below retain their actual m/s units and explicit inlet/property/closure assumptions; they are not measured phase velocities or converted pixel results.', '',
        '**We corrected the comparison and tested the alternatives. The results are conditional estimates; actual gas and liquid phase velocities are not yet validated.** No defensible experimental accuracy percentage or unique best physical closure can be inferred from these records.', '',
        f'The analysis uses **{n} acquisitions**: 1 dispersed-bubble, 3 plug, 3 slug, 10 churn and 9 wavy-annular. It excludes m05q05_002, smooth annular and all double classifications. Saturated m03q40 remains in the thermodynamic calculation and is flagged as unsupported for optical interpretation.', '',
        '## Conditional gas × liquid model predictions', '',
        '![Conditional velocity comparison with nominal and recorded flow](../figures/recorded_flow_comparison.png)', '',
        '**Use the documented recorded mass flow for a condition-level alternative, alongside the nominal-G baseline.** Processed values match the original LabVIEW files. However, condition names alone do not establish synchronization with each video. The suffix-free m05q05 sensor blocks cannot distinguish camera _001 from excluded _002; their mean is only a condition-level sensitivity input. Both baselines assume saturated liquid at the heated inlet, net nominal heating, R600a at 27 °C and a 1 mm tube.', '',
        'The table reports **medians across acquisitions**, in m/s. Jg/Jl are superficial velocities. ug/ul are **Rouhani–Axelsson (RA) predictions** using recorded mean flow; they are not measured phase speeds. RA is a provisional reporting reference, not an experimentally established winner.', '',
        '| Pattern | n | Jg, either G source | Jl, nominal G | Jl, recorded flow | ug, recorded flow + RA | ul, recorded flow + RA |',
        '|---|---:|---:|---:|---:|---:|---:|']
    medians = []
    for pattern in PATTERNS:
        r = ra[ra.pattern == pattern]
        old = nominal_ra[nominal_ra.pattern == pattern]
        item = dict(pattern=pattern, n=len(r), Jg=float(r.Jg.median()), Jl_nominal=float(old.Jl.median()),
                    Jl_recorded=float(r.Jl.median()), ug_recorded_RA=float(r.ug.median()), ul_recorded_RA=float(r.ul.median()))
        medians.append(item)
        lines.append(f"| {pattern} | {len(r)} | {item['Jg']:.3f} | {item['Jl_nominal']:.3f} | {item['Jl_recorded']:.3f} | {item['ug_recorded_RA']:.3f} | {item['ul_recorded_RA']:.3f} |")
    lines += ['',
        'Jg does not change with G under the assumed zero-inlet-quality energy balance: G cancels algebraically. This is model structure, not experimental confirmation. Between-pattern differences also reflect different heat/mass-flux conditions and unequal sample counts. Three displayed decimals facilitate reproducibility; they do not imply measurement precision.', '',
        'Full tables: [26 acquisitions and both baselines](velocity_results.md), [recorded-flow CSV](../results/recorded_flow_conditional_velocities_per_run.csv), [nominal CSV](../results/conditional_velocities_per_run.csv).', '',
        '## What the investigation changed', '',
        '1. **The legacy arithmetic reproduces.** All 280 displayed per-run values and 90 summary values agree with independent recomputation. The defects are scope, interpretation and unsupported accuracy claims. Superficial liquid velocity is not liquid material speed, and a slip model cannot repair missing inlet enthalpy.',
        f"2. **Measured-flow evidence exists.** The owner identifies Voltage_1 as g/s in the original LabVIEW data; the similarly named src2 image-derived columns are zero placeholders. We verified 241,000 samples against 241 original files. Recorded mean G differs from the nominal lookup by **{metrics['massflux_change_pct_min']:.2f}% to +{metrics['massflux_change_pct_max']:.2f}%**. Temporal camera pairing and instrument calibration remain unresolved.",
        f"3. **Phase velocity depends strongly on the closure.** Four closures and 13 input scenarios produce **{metrics['thermodynamic_experiment_rows']:,} calculations**. At nominal G, the largest liquid-velocity ratio between tested closures is **{metrics['maximum_nominal_liquid_model_factor']:.2f}×**. This is model spread, not an error bar.",
        f"4. **Inlet uncertainty is decisive.** Assuming only 1 K of inlet subcooling makes **{metrics['subcool_1K_out_of_domain_runs']}/{n}** nominal-G cases fall outside the two-phase equilibrium domain. Those cases are flagged, not clipped into apparently valid velocities.",
        '5. **The source literature is not a ground-truth table.** Published x/alpha/slip combinations are inconsistent. The closest paper prints a Xu–Fang exponent of −2, while a pinned secondary implementation uses −0.2. Both are tested separately; the original correlation paper remains unavailable.', '',
        '![Phase model sensitivity at the same nominal conditions](../figures/phase_model_comparison.png)', '',
        f"Xu–Fang's secondary-source Froude limit is exceeded by **{metrics['nominal_XF_Fr_extrapolated_runs']}/{n} nominal-G** and **{metrics['recorded_XF_Fr_extrapolated_runs']}/{n} recorded-flow** cases. No published error percentage is transferred to these predictions.", '',
        '## Which optical option worked best in these tests?', '',
        '**Signed normalized cross-correlation after subtracting the temporal background, without high-pass filtering**, is the preferred candidate for further validation. It offers the broadest numerical support among the tested methods; this does not establish physical accuracy.', '',
        '| Optical method | Passing adjacent pairs / 13,000 | Runs with ≥50% passing pairs |',
        '|---|---:|---:|']
    names = {'ncc_raw':'Raw-profile NCC', 'ncc_dynamic':'Background-subtracted NCC',
             'ncc_dynamic_highpass':'Background-subtracted + high-pass NCC', 'phase_dynamic_highpass':'Phase-only correlation, high-pass'}
    for method in names:
        o = all_optical[all_optical.method == method]
        lines.append(f"| {names[method]} | {int(o.heuristic_pass_pairs.sum()):,} | {int((o.heuristic_pass_fraction >= .5).sum())}/{n} |")
    lines += ['',
        f"The preferred candidate passes **{metrics['preferred_optical_passing_pairs']:,}/{metrics['optical_pairs']:,} pairs**. Excluding saturated images from supported interpretations leaves **{len(acceptable)}/{n} supported recordings and zero wavy-annular recordings**. Different estimator thresholds and preprocessing affect acceptance; coverage is not an accuracy score.", '',
        'Synthetic translations favor NCC over phase-only correlation for the constructed textures/noise, but do not simulate all real interface deformation. All four methods use native timing; crop/lag tests preserve rejected and low-support outcomes. Results stay in **px/frame and px/s**. A structure speed cannot independently supply both liquid and gas phase speeds, even after a pixel-to-mm calibration.', '',
        'Optical results: [preferred candidate](../results/optical_preferred_candidate_per_run.csv), [full investigation and synthetic results](optical_investigation.md). Files named optical_primary_pairs/optical_per_run preserve the high-pass reference experiment, not the preferred candidate.', '',
        '## Unresolved inputs and validation limits of this model investigation', '',
        '- **Unknown inlet enthalpy and net heat input:** obtain these, plus the absolute state at the camera, to determine vapor quality and superficial velocities.',
        '- **Flow-to-video matching:** establish which LabVIEW interval belongs to each acquisition and verify the engineering-unit calibration.',
        '- **Independent phase reference:** validate volumetric void fraction or measure phase-resolved velocity with an appropriate averaging convention. Inserting v14 spatial image alpha is not sufficient; its unvalidated geometry and crop sensitivity can produce extreme velocities.',
        '- **Optical scale and tracking:** calibrate axial length and validate feature correspondence. Annular motion remains poorly supported by the tested optical methods.',
        '- **Source conflict:** obtain the original Xu–Fang equation or author clarification before an authoritative attribution.', '',
        'This historical investigation produced the conditional superficial comparison, an explicitly assumed RA phase estimate with closure sensitivity, and a separate optical structure-motion diagnostic. The later [current executive summary](EXECUTIVE_SUMMARY.md) evaluates image-motion descriptors using the available dataset. **A measured-accuracy winner remains undetermined.**', '',
        'Details: [method and equations](method_comparison.md), [legacy audit](legacy_audit.md), [mass-flow audit](massflow_review.md), [literature](literature_review.md), [reproduction](../README.md).']
    (BASE / 'findings/conditional_model_summary.md').write_text('\n'.join(lines) + '\n')

    detailed = ['# Historical conditional model velocity results', '',
        'These are the earlier model tables. See the [current executive summary](EXECUTIVE_SUMMARY.md) for the phase-passage analysis in px/s and dimensionless ratios, or the [historical model summary](conditional_model_summary.md) for this calculation and its assumptions.', '',
        'All velocities below are conditional predictions in m/s. Reference assumptions: x_in = 0, eta = 1, R600a at 27 °C, D = 1 mm and heated length 265 mm. RA is the reporting reference. Full-precision CSVs retain all four closures, sensitivity scenarios and Froude flags.', '',
        '## Recorded-mean-flow alternative, by acquisition', '',
        'Flow input is documented in g/s and verified against the original LabVIEW values. Matching to the video is by condition name only. The averaging and instrument limitations apply to every row.', '',
        '| Acquisition | Pattern | G [kg/m²s] | x | Jg | Jl | alpha RA | ug RA | ul RA | S RA |',
        '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in ra.to_dict('records'):
        detailed.append(f"| {r['run']} | {r['pattern']} | {r['G_kg_m2_s']:.2f} | {r['x_eq']:.5f} | {r['Jg']:.3f} | {r['Jl']:.3f} | {r['alpha']:.4f} | {r['ug']:.3f} | {r['ul']:.3f} | {r['slip']:.3f} |")
    detailed += ['', '## Nominal-G baseline, by acquisition', '',
        'This reproduces the old model after the requested exclusions. It is not fitted or calibrated to the images.', '',
        '| Acquisition | Pattern | G [kg/m²s] | x | Jg | Jl | alpha RA | ug RA | ul RA | S RA |',
        '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in nominal_ra.to_dict('records'):
        detailed.append(f"| {r['run']} | {r['pattern']} | {r['G_kg_m2_s']:.0f} | {r['x_eq']:.5f} | {r['Jg']:.3f} | {r['Jl']:.3f} | {r['alpha']:.4f} | {r['ug']:.3f} | {r['ul']:.3f} | {r['slip']:.3f} |")
    detailed += ['', '## Model comparison at nominal G: median predicted phase speeds', '',
        '| Pattern | Homogeneous ug = ul | RA ug / ul | XF secondary −0.2 ug / ul | Literal 2021 −2 ug / ul |',
        '|---|---:|---:|---:|---:|']
    for p in PATTERNS:
        vals = nominal[nominal.pattern == p].groupby('model')[['ug','ul']].median()
        pairs = [f"{vals.loc[m,'ug']:.3f} / {vals.loc[m,'ul']:.3f}" for m in ['rouhani_axelsson','xu_fang_fluids','xu_fang_literal_2021']]
        detailed.append(f"| {p} | {vals.loc['homogeneous','ug']:.3f} | " + ' | '.join(pairs) + ' |')
    detailed += ['', '## Input sensitivity relative to nominal saturated-inlet baseline', '',
        'Percent ranges span the 26 acquisitions. Invalid two-phase equilibrium cases are counted and omitted from numerical percent ranges, never treated as zero velocity. These are hypothetical one-at-a-time tests, not statistical uncertainty bounds.', '',
        '| Scenario | Out-of-domain cases | Jg change [%] | Jl change [%] | RA ug change [%] | RA ul change [%] |',
        '|---|---:|---:|---:|---:|---:|']
    for scenario, group in sensitivity.groupby('scenario', sort=False):
        ranges = [f'{group[k].min():+.2f} to {group[k].max():+.2f}' for k in ['Jg_change_pct','Jl_change_pct','ug_change_pct','ul_change_pct']]
        detailed.append(f"| {scenario} | {int((group.quality_status != 'two_phase').sum())}/26 | " + ' | '.join(ranges) + ' |')
    h = hybrid[(hybrid.run == 'm07q01_001') & (hybrid.variant == 'D60')].iloc[0]
    detailed += ['', '![Input sensitivity](../figures/inlet_sensitivity.png)', '',
        '## Image-alpha hybrid stress test', '',
        f"For m07q01_001, combining nominal energy-model Jg with v14 image alpha={h.alpha_image:.6f} yields ug={h.ug:.2f} m/s and S={h.slip:.2f}. These are **unvalidated hybrid outputs, not accepted phase-velocity estimates**. Saturated m03q40 has no hybrid velocity result. Other cases can yield S < 1; this alone is not a proof of impossibility. The test cannot determine whether inlet quality, optical geometry, sampling or closure is responsible.", '',
        '[All hybrid cases and diameter assumptions](../results/image_alpha_diagnostics.csv). See [method comparison](method_comparison.md) for averaging and sensitivity limits.', '',
        '## Reproducibility and exclusions', '',
        '[Scenario inventory](../results/scenario_inventory.csv) preserves all 31 folder identities, the five exclusion reasons, native timing and input hashes. Every reported velocity comparison uses only the same 26 selected acquisitions. Synthetic dataset folders are not experimental scenarios. Group CSVs contain min/max/mean/median; summaries are descriptive statistics over the sampled conditions.']
    (BASE / 'findings/velocity_results.md').write_text('\n'.join(detailed) + '\n')
    metrics['recorded_RA_pattern_medians'] = medians
    (results / 'validation/executive_summary_values.json').write_text(json.dumps(metrics, indent=2) + '\n')
    publish()
    print(json.dumps(metrics, indent=2))


if __name__ == '__main__':
    main()
