#!/usr/bin/env python3
"""Integration/provenance checks; not experimental velocity validation."""
from pathlib import Path
import hashlib
import json
import re
from urllib.parse import unquote
import numpy as np
import pandas as pd
from publish_current_summary import current_summary_text

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_current_summary():
    """Check current units, source selection, displayed means and plot data."""
    phase = BASE / 'phase_conditioned_motion'
    reports = [phase / 'findings/EXECUTIVE_SUMMARY.md', BASE / 'findings/EXECUTIVE_SUMMARY.md']
    assert reports[1].read_text() == current_summary_text(), 'Main summary is stale'
    ratios = pd.read_csv(phase / 'results/passage_ratios.csv')
    reference = ratios[(ratios.config == 'consensus_reference_guard')
                       & (ratios.motion_gate == 'lag2_confirmed')]
    assert len(reference) == 26
    accepted = reference[reference.eligible_descriptive_ratio].set_index('run')
    expected_runs = {'m05q01_001', 'm05q05_001', 'm06q01_001', 'm07q05_001'}
    assert set(accepted.index) == expected_runs and len(accepted) == 4
    phase_rows = pd.read_csv(phase / 'results/passage_phase_summary.csv')
    phase_rows = phase_rows[(phase_rows.config == 'consensus_reference_guard')
                            & (phase_rows.motion_gate == 'lag2_confirmed')]
    plot_path = phase / 'results/conditional_passage_velocities.csv'
    plotted = pd.read_csv(plot_path).set_index('run')
    assert plotted.index.is_unique and set(plotted.index) == expected_runs
    assert (plotted.config == 'consensus_reference_guard').all()
    assert (plotted.motion_gate == 'lag2_confirmed').all()
    expected = {}
    numeric_cells = 0
    for run, row in accepted.iterrows():
        candidates = phase_rows[phase_rows.run == run]
        assert len(candidates) == 2 and candidates.fps.nunique() == 1
        assert set(candidates.phase) == {'gas_candidate', 'liquid_candidate'}
        fps = float(candidates.fps.iloc[0])
        values = dict(fps=fps,
                      gas_passage_mean_px_s=float(row.gas_mean_px_frame) * fps,
                      liquid_gap_mean_px_s=float(row.liquid_mean_px_frame) * fps,
                      ratio_of_means=float(row.descriptive_ratio),
                      episode_equal_ratio=float(row.episode_equal_ratio),
                      gas_pairs=int(row.gas_pairs), liquid_pairs=int(row.liquid_pairs),
                      gas_episodes=int(row.gas_episodes), liquid_episodes=int(row.liquid_episodes))
        assert np.isclose(values['gas_passage_mean_px_s'] / values['liquid_gap_mean_px_s'],
                          values['ratio_of_means'], rtol=1e-13)
        for candidate, field in [('gas_candidate', 'gas_passage_mean_px_s'),
                                 ('liquid_candidate', 'liquid_gap_mean_px_s')]:
            source_value = candidates.loc[candidates.phase == candidate, 'mean_px_s'].iloc[0]
            assert np.isclose(source_value, values[field], rtol=1e-13)
        assert plotted.loc[run, 'pattern'] == row['pattern']
        for key, value in values.items():
            assert np.isclose(plotted.loc[run, key], value, rtol=1e-13), (run, key)
            numeric_cells += 1
        expected[run] = [run, row['pattern'],
                         f"{values['gas_passage_mean_px_s']:,.0f}",
                         f"{values['liquid_gap_mean_px_s']:,.0f}",
                         f"{values['ratio_of_means']:.3f}",
                         f"{values['gas_pairs']}/{values['liquid_pairs']}",
                         f"{values['gas_episodes']}/{values['liquid_episodes']}"]
    table_rows = 0
    embedded_images = set()
    for report in reports:
        markdown = report.read_text()
        assert 'Mean during gas-like passages (px/s)' in markdown
        assert 'Mean during liquid-like gaps (px/s)' in markdown
        checked = set()
        for line in markdown.splitlines():
            cells = [cell.strip() for cell in line.strip('|').split('|')]
            if len(cells) == 7 and cells[0] in expected_runs:
                assert cells == expected[cells[0]], (report.name, cells)
                assert cells[0] not in checked
                checked.add(cells[0])
                table_rows += 1
        assert checked == expected_runs
        images = re.findall(r'!\[[^\]\n]*\]\(([^\)\n]+)\)', markdown)
        assert any('conditional_passage_velocities.png' in image for image in images)
        for image in images:
            target = (report.parent / unquote(image.split('#')[0])).resolve()
            assert target.is_relative_to((phase / 'figures').resolve()), target
            assert target.exists()
            embedded_images.add(str(target.relative_to(BASE)))
    chart = phase / 'figures/conditional_passage_velocities.svg'
    svg = chart.read_text()
    assert 'Signed mean image motion (px/s)' in svg
    assert 'Gas-passage mean / liquid-gap mean (dimensionless)' in svg
    assert 'm/s' not in svg, 'Physical units in the current image-motion chart'
    return dict(status='passed', eligible_recordings=4, report_table_rows=table_rows,
                plotted_numeric_cells=numeric_cells,
                configuration='consensus_reference_guard / lag2_confirmed',
                velocity_units='px/s', ratio_units='dimensionless',
                current_plot_csv_sha256=digest(plot_path),
                current_plot_svg_sha256=digest(chart), embedded_images=sorted(embedded_images))


def main():
    results = BASE / 'results'
    inv = pd.read_csv(results / 'scenario_inventory.csv')
    selected = set(inv.loc[inv.included, 'run'])
    assert len(selected) == 26 and len(inv) == 31
    assert set(inv.loc[~inv.included, 'run']) == {
        'm03q01_001', 'm03q80_001', 'm04q05_001', 'm05q05_002', 'm07q10_001'}
    scoped_tables = []
    for path in sorted(results.glob('*.csv')):
        if path.name == 'scenario_inventory.csv':
            continue
        frame = pd.read_csv(path)
        if 'run' in frame:
            assert set(frame.run) <= selected, f'Out-of-scope result in {path}'
            assert set(frame.run) == selected, f'Missing scenario in {path}'
            scoped_tables.append(path.name)
    for row in inv.to_dict('records'):
        folder = ROOT / 'src2/data_processed' / row['run']
        assert digest(folder / '_config.xsv') == row['config_sha256']
        assert digest(folder / '_config.mii') == row['mii_sha256']
        assert digest(ROOT / row['marker_file']) == row['marker_sha256']
    legacy = json.loads((results / 'validation/legacy_audit.json').read_text())
    assert legacy['original_sources_unchanged']
    for name, value in legacy['original_source_hashes'].items():
        assert digest(ROOT / 'src2' / name) == value
    model = json.loads((results / 'validation/model_metadata.json').read_text())
    assert digest(BASE / 'scripts/model_comparison.py') == model['script_sha256']
    assert digest(BASE / model['image_input']['snapshot']) == model['image_input']['snapshot_sha256']
    assert digest(BASE / model['massflow_input']['path']) == model['massflow_input']['sha256']
    independent = json.loads((results / 'validation/independent_model_review.json').read_text())
    assert independent['status'] == 'passed'
    for name, value in independent['reviewed_input_sha256'].items():
        assert digest(BASE / name) == value, name
    for name, value in independent.get('reviewed_narrative_sha256', {}).items():
        assert digest(BASE / name) == value, name
    assert json.loads((results / 'optical_validation.json').read_text())['passed']
    assert json.loads((results / 'validation/model_checks.json').read_text())['status'] == 'passed'
    assert json.loads((results / 'validation/massflow_audit.json').read_text())['raw_processed_all_channels_match']

    headline = json.loads((results / 'validation/executive_summary_values.json').read_text())
    ra = pd.read_csv(results / 'recorded_flow_conditional_velocities_per_run.csv')
    ra = ra[ra.model == 'rouhani_axelsson']
    for row in headline['recorded_RA_pattern_medians']:
        subset = ra[ra.pattern == row['pattern']]
        assert len(subset) == row['n']
        for field, target in [('Jg','Jg'), ('Jl','Jl_recorded'), ('ug','ug_recorded_RA'), ('ul','ul_recorded_RA')]:
            assert np.isclose(subset[field].median(), row[target], rtol=1e-13)
    pair = pd.read_csv(results / 'optical_preferred_candidate_pairs.csv')
    assert len(pair) == headline['optical_pairs'] == 13000
    assert int(pair.heuristic_pass.sum()) == headline['preferred_optical_passing_pairs'] == 7528
    preferred = pd.read_csv(results / 'optical_preferred_candidate_per_run.csv')
    assert int((~preferred.low_coverage & ~preferred.high_saturation).sum()) == headline['preferred_optical_supported_runs'] == 17
    current_report = verify_current_summary()

    json_count = 0
    for path in sorted(results.rglob('*.json')):
        json.loads(path.read_text())
        json_count += 1
    link_count, markdown_count = 0, 0
    for path in sorted(BASE.rglob('*.md')):
        markdown_count += 1
        for target in re.findall(r'\[[^\]\n]*\]\(([^\)\n]+)\)', path.read_text()):
            target = target.strip('<>')
            if target.startswith(('https://', 'http://', 'mailto:', '#')):
                continue
            target = unquote(target.split('#')[0])
            assert (path.parent / target).exists(), f'Broken link: {path.relative_to(BASE)} -> {target}'
            link_count += 1

    output = dict(status='passed', selected_acquisitions=26, scope_checked_tables=scoped_tables,
                  native_config_and_marker_hashes_checked=31,
                  original_legacy_source_hashes_checked=len(legacy['original_source_hashes']),
                  independent_model_numeric_cells=independent['numeric_cells_compared'],
                  independent_summary_cells=independent['summary_cells_compared'],
                  independent_report_cells=independent.get('narrative_numeric_cells_compared'),
                  checked_json_files=json_count, checked_markdown_files=markdown_count,
                  existing_relative_links=link_count,
                  optical_preferred_pair_count=len(pair),
                  current_phase_summary=current_report,
                  physical_accuracy='unvalidated; arithmetic, scope and provenance only')
    manifest = {}
    excluded_names = {'artifact_manifest.json', 'final_verification.json', 'pipeline_run.json'}
    for path in sorted(BASE.rglob('*')):
        relative = path.relative_to(BASE)
        if not path.is_file() or path.name in excluded_names:
            continue
        if '__pycache__' in relative.parts or any('cache' in part or part in ('matplotlib','logs') for part in relative.parts):
            continue
        manifest[str(relative)] = dict(sha256=digest(path), bytes=path.stat().st_size)
    (results / 'validation/artifact_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    output['manifest_artifacts'] = len(manifest)
    (results / 'validation/final_verification.json').write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
