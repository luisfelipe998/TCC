#!/usr/bin/env python3
"""Generate a short illustrated executive summary and folder indexes."""
import csv
import hashlib
import json

from panorama_void_alternatives import DIAMETER_PX, REFERENCE_RUNS
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
    verification = json.loads(artifact_path('final_verification.json').read_text())
    void = {r['run']: r for r in rows('panorama_void_per_run.csv')}
    phase = {r['run']: r for r in rows('phase_mask_flip_per_run.csv')}
    groups = rows('panorama_void_by_pattern.csv')
    phase_groups = {r['pattern']: r for r in rows('phase_mask_flip_by_pattern.csv')}
    adversarial = json.loads(artifact_path('adversarial_checks.json').read_text())
    crop_example = next(r for r in adversarial['crop_sensitivity'] if r['run'] == 'm06q80_001')
    crop_values = [r['conditional_alpha'] for r in crop_example['crop_variants'].values()]
    assert set(void) == set(phase) == set(selected_runs())
    n_runs = len(void)
    n_patterns = len({r['label'] for r in void.values()})
    n_eligible = sum(r['eligible'] == 'True' for r in phase.values())
    radial_rows = {int(r['candidate_mask_n_rows']) for r in phase.values()}
    assert len(radial_rows) == 1
    n_rows = next(iter(radial_rows))
    ceiling = (n_rows / DIAMETER_PX) ** 2
    assert sum(int(r['n']) for r in groups) == n_eligible
    for figure in ('void_fraction_vs_flip_rate.png', 'void_fraction_vs_center_flips.png',
                   'void_crop_sensitivity.png', 'panorama_phase_masks.png'):
        if not artifact_path(figure).is_file():
            raise FileNotFoundError(f'Generate {figure} with summary_figures.py before writing the summary.')

    report = [
        '# Executive summary', '',
        '**We estimated visible gas and spatial gas/liquid switches. Absolute gas-volume fractions remain uncalibrated.** Results depend on detected outlines and assumed circular geometry.', '',
        '## Void fraction × spatial intermittency', '',
        '![Conditional visible void estimate versus spatial gas-mask flip-rate, by supplied pattern](void_fraction_vs_flip_rate.png)', '',
        f'The current analysis retains **{n_runs} acquisitions and {n_patterns} single patterns**. Colored areas show sample spread, not proven pattern boundaries. The group table uses {n_eligible} pure-pattern, unsaturated acquisitions. Hollow points have unsupported estimates because their images are saturated.', '',
        f'Y-axis: gas/liquid changes per 1,000 neighboring pixel comparisons, averaged over {n_rows} rows. Higher means more changes **along the image**, not per second. The [center-line version](void_fraction_vs_center_flips.png) misses off-center bubbles.', '',
        f'**Limits:** plug, slug and churn overlap. Both axes use the same gas map, so agreement does not verify true void fraction. The artificial crop ceiling at **{ceiling:.4f}** prevents clustered annular points from proving precision.', '',
        '| Pattern / acquisition | n | Visible fraction α | Spatial flips / 1,000 comparisons |',
        '|---|---:|---:|---:|',
    ]
    for run in sorted(r['run'] for r in void.values() if r['label'].replace('_', ' ') == 'dispersed bubbles'):
        report.append(f"| Dispersed bubbles: {run} | 1 | {float(void[run]['conditional_circle_alpha_spatial_mean']):.4f} | {1000 * float(phase[run]['candidate_pixel_flip_rate']):.3f} |")
    for group in groups:
        pattern = group['pattern']
        if pattern != 'dispersed bubbles':
            report.append(f"| {pattern} | {group['n']} | {float(group['conditional_alpha_median']):.4f} | {1000 * float(phase_groups[pattern]['candidate_pixel_flip_rate_median']):.3f} |")
    report += [
        '', 'Group values are medians. Only m07q01_001 represents dispersed bubbles; one run cannot establish within-pattern variability or independently validate that class.', '',
        'Mixed classifications, smooth annular and m05q05_002 are excluded from every calculation, reference check, classifier and chart. Annular means wavy annular throughout. Previous results are [archived](intermediate/history/pre_scope_exclusions/README.md).', '',
        '## How we reached this result', '',
        f'1. **Read images and labels.** We retained {n_runs} single-pattern panoramas under the [selection policy](dataset_selection.json). The one saturated retained image is shown as unsupported and omitted from group statistics.', '',
        f'2. **Find liquid-looking gaps.** Gaps in {len(REFERENCE_RUNS)} plug panoramas provide a liquid-appearance reference. **Limit:** these gaps are inferred, not a confirmed liquid-only calibration.', '',
        '3. **Trace gas outlines and fill their insides.** Like coloring a balloon, the bright middle belongs inside its outline. We kept separate bubbles and rejected short false bridges. **Limit:** broken outlines can still cause errors.', '',
        '![Original panoramas alongside the inferred gas silhouettes in blue](panorama_phase_masks.png)', '',
        f'4. **Turn width into area.** We converted separate gas widths to circular areas, added them at each column, then averaged along the image. Assumed apparent tube diameter: {DIAMETER_PX:g} pixels. **Limit:** the tube edges are unvalidated; a side view hides shape in depth and optical distortion remains uncertain.', '',
        f'5. **Count switches.** We compared neighboring gas/liquid mask pixels across {n_rows} rows, keeping counts and denominators. **Limit:** spacing and stitching affect this spatial measure; it is not a time frequency.', '',
        '6. **Try to break the result.** We changed thresholds, liquid references and crops, and inspected clear gas/liquid patches. This found a false bridge that was corrected. **Limit:** clear patches do not validate every boundary.', '',
        '![Conditional void estimates when only the radial crop changes](void_crop_sensitivity.png)', '',
        f'**The biggest demonstrated uncertainty is crop width.** For annular m06q80, crop alone moves α from **{min(crop_values):.4f} to {max(crop_values):.4f}**. This is sensitivity, not a confidence interval. Wider crops can include wall reflections.', '',
        f'7. **Check calculations and sources.** We independently verified all {verification["saved_binary_masks_verified"]} gas maps, void calculations, flip counts and input hashes. **Limit:** the local papers provide methods, not verified true values for these images.', '',
        '## Findings and next steps', '',
        'The outline method avoids counting a bright liquid middle as a large gas core. Annular images usually show broad gas cores with few switches. Pattern identification remains unresolved; held-out training can lack the only bubble example, and those test cases remain in the reported scores.', '',
        '- Annotate gas/liquid outlines across all panoramas and measure mask errors.',
        '- Validate tube boundaries, optical distortion and circular geometry before claiming absolute void fractions.', '',
        'Use the numbers as **conditional image estimates and spatial flip-rates**. All values come from saved results; flip-rates are multiplied by 1,000 for readability.', '',
        'Details: [per-run void](panorama_void_per_run.csv), [per-run flips](phase_mask_flip_per_run.csv), [method](panorama_void_alternatives.md), [flip definitions](phase_mask_flip_review.md), [adversarial checks](adversarial_review.md), [literature](source_audit.md), [post-regeneration adversarial review](post_regeneration_adversarial_review.md).',
    ]
    write_markdown(artifact_path('EXECUTIVE_SUMMARY.md'), '\n'.join(report) + '\n')
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
