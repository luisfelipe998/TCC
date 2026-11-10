#!/usr/bin/env python3
"""Generate a short illustrated executive summary and folder indexes."""
import csv
import ast
import hashlib
import json

from v14_geometry import DIAMETER_PX, INTERNAL_CROP_VARIANTS, LUMEN_HEIGHT_PX
from v14_navigation import write_navigation
from v14_paths import ROOT, artifact_path, write_markdown
from v14_selection import selected_runs, scope_metadata


def rows(name):
    with artifact_path(name).open() as fp:
        return list(csv.DictReader(fp))


def main():
    chart_checks = json.loads(artifact_path('summary_chart_checks.json').read_text())
    assert chart_checks['passed'], 'Chart verification failed'
    for item in chart_checks['source_inputs']:
        digest = hashlib.sha256((ROOT / item['artifact']).read_bytes()).hexdigest()
        assert digest == item['sha256'], 'Chart inputs changed; rerun summary_figures.py'
    script_hash = hashlib.sha256(artifact_path('summary_figures.py').read_bytes()).hexdigest()
    assert script_hash == chart_checks['script_sha256'], 'Chart code changed; rerun summary_figures.py'
    void = {r['run']: r for r in rows('panorama_void_per_run.csv')}
    phase = {r['run']: r for r in rows('phase_mask_flip_per_run.csv')}
    groups = rows('panorama_void_by_pattern.csv')
    phase_groups = {r['pattern']: r for r in rows('phase_mask_flip_by_pattern.csv')}
    adversarial = json.loads(artifact_path('adversarial_checks.json').read_text())
    crop_example = next(r for r in adversarial['crop_sensitivity'] if r['run'] == 'm06q80_001')
    assert {ast.literal_eval(key) for key in crop_example['crop_variants']} == set(INTERNAL_CROP_VARIANTS)
    crop_values = [r['conditional_alpha'] for r in crop_example['crop_variants'].values()]
    assert set(void) == set(phase) == set(selected_runs())
    n_runs = len(void)
    n_patterns = len({r['label'] for r in void.values()})
    n_eligible = sum(r['eligible'] == 'True' for r in phase.values())
    radial_rows = {int(r['candidate_mask_n_rows']) for r in phase.values()}
    assert len(radial_rows) == 1
    n_rows = next(iter(radial_rows))
    assert n_rows == LUMEN_HEIGHT_PX == DIAMETER_PX == 53
    assert sum(int(r['n']) for r in groups) == n_eligible
    for figure in ('void_fraction_vs_flip_rate.png', 'void_fraction_vs_center_flips.png',
                   'void_crop_sensitivity.png', 'panorama_phase_masks.png'):
        if not artifact_path(figure).is_file():
            raise FileNotFoundError(f'Generate {figure} with summary_figures.py before writing the summary.')

    from summary_text import build_summary, hull_overlap_pairs
    overlaps = hull_overlap_pairs(chart_checks, rows('void_fraction_flip_plot_data.csv'))
    effect = json.loads(artifact_path('fix_effectiveness.json').read_text())
    assert effect['passed']
    for pt, filename in [(False, 'EXECUTIVE_SUMMARY.md'), (True, 'EXECUTIVE_SUMMARY_PT_BR.md')]:
        report = build_summary(pt=pt, groups=groups, phase_groups=phase_groups,
                               n_runs=n_runs, n_rows=n_rows, diameter=DIAMETER_PX,
                               crop_values=crop_values, effect=effect, overlap_pairs=overlaps)
        write_markdown(artifact_path(filename), report)
    refresh_source_scope()
    write_navigation()
    print('Wrote the illustrated executive summary and folder indexes from saved results.')


def refresh_source_scope():
    import re
    from collections import Counter
    path = artifact_path('source_audit.md')
    content = path.read_text()
    start, stop = '<!-- current-dataset-scope:start -->', '<!-- current-dataset-scope:end -->'
    if start in content:
        before, _, rest = content.partition(start)
        _, _, after = rest.partition(stop)
        content = before.rstrip() + after
    note = (start + '\n**Current dataset scope:** ' + str(len(selected_runs()))
            + ' retained acquisitions in five single patterns; annular means wavy annular. '
            + 'Mixed labels, smooth annular and m05q05_002 are excluded from current calculations. '
            + 'Published literature values below remain source evidence, not current dataset measurements.\n'
            + stop + '\n')
    first, rest = content.split('\n', 1)
    content = first + '\n\n' + note + rest.lstrip('\n')
    native = rows('flip_rate_per_run.csv')
    counts = dict(sorted(Counter(int(float(r['native_camera_fps'])) for r in native).items()))
    frame_counts = sorted({int(r['n_native_frames']) for r in native})
    replacement = (f"Native config evidence for the current {len(native)} retained acquisitions: frame counts {frame_counts}; acquisition counts by native fps {counts}. "
                   + "The actual endpoint-to-endpoint record duration is (N−1)/fps. ")
    pattern = r'Native config evidence (?:in|for the current).*?(?=' + chr(96) + '__timeseries_v9.xlsx' + chr(96) + ')'
    content = re.sub(pattern, replacement, content, count=1)
    content = content.replace('for all 31 ' + chr(96) + 'src2' + chr(96) + ' records', 'for the retained ' + chr(96) + 'src2' + chr(96) + ' records')
    write_markdown(path, content)


if __name__ == '__main__':
    main()
