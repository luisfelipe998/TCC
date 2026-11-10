#!/usr/bin/env python3
"""Audit the current results AFTER the analysis and executive summary are written.

Selection is independently derived from the original folder markers. This script
does not import the shared selection helper or rerun an analysis stage.
"""
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import importlib
import json
import math
from pathlib import Path
import re
import warnings

import numpy as np
from PIL import Image

from v14_paths import DATA, ROOT, artifact_path, write_markdown

ALLOWED = ('dispersed bubbles', 'plug', 'slug', 'churn', 'wavy annular')
EXPLICITLY_EXCLUDED = {'m05q05_002'}
RUN_ID = re.compile(r'm\d{2}q\d{2}_\d{3}')
ARCHIVE = ROOT / 'intermediate/history/pre_scope_exclusions'
CSV_AREAS = ('results', 'intermediate/features', 'intermediate/sensitivity',
             'intermediate/exploratory')
LABEL_COLUMNS = {'pattern', 'label', 'truth', 'prediction', 'supplied_pattern',
                 'source_void_label', 'true_pattern', 'predicted_pattern'}
ALL_RUN_TABLES = (
    'data_inventory.csv', 'panorama_void_per_run.csv', 'phase_mask_flip_per_run.csv',
    'flip_rate_per_run.csv', 'image_void_alternatives.csv',
    'pattern_morphology_features.csv', 'void_fraction_flip_plot_data.csv',
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_csv(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def table(name):
    return read_csv(artifact_path(name))


def by_run(name):
    rows = table(name)
    keyed = {row['run']: row for row in rows}
    require(len(keyed) == len(rows), f'Duplicate run in {name}')
    return keyed


def close(actual, expected, context, tolerance=1e-12):
    actual, expected = float(actual), float(expected)
    require(math.isfinite(actual) and math.isfinite(expected), f'Nonfinite value: {context}')
    require(abs(actual - expected) <= tolerance, f'{context}: {actual} != {expected}')


def discover_scope():
    """Read raw markers and original pixels without using v14_selection."""
    selected, excluded, eligible, saturation = {}, {}, set(), {}
    folders = [p for p in sorted(DATA.iterdir()) if p.is_dir() and RUN_ID.fullmatch(p.name)]
    for folder in folders:
        markers = [p for p in folder.iterdir() if p.name.startswith('__') and p.suffix == '.txt']
        require(len(markers) == 1, f'Unexpected pattern markers in {folder.name}')
        pattern = markers[0].name[2:-4].replace('_', ' ')
        if pattern not in ALLOWED or folder.name in EXPLICITLY_EXCLUDED:
            excluded[folder.name] = pattern
            continue
        image = np.asarray(Image.open(folder / '__panorama_v11.png').convert('L'))
        require(image.ndim == 2 and image.shape[0] == 80, f'Unexpected panorama shape: {folder.name}')
        selected[folder.name] = pattern
        saturation[folder.name] = float(np.count_nonzero(image[20:60] == 255) / image[20:60].size)
        if saturation[folder.name] <= .20:
            eligible.add(folder.name)
    require(len(selected) == 26, f'Independent selection found {len(selected)}, expected 26')
    require(len(eligible) == 25, f'Independent saturation rule found {len(eligible)}, expected 25')
    require(set(selected.values()) == set(ALLOWED), 'A requested class is absent from selected data')
    return selected, excluded, eligible, saturation


def scan_csv_scope(selected, excluded, eligible):
    checked, rows_checked = [], 0
    for area in CSV_AREAS:
        for path in sorted((ROOT / area).rglob('*.csv')):
            rows = read_csv(path)
            if path.name == 'dataset_exclusions.csv':
                require({r['run']: r['supplied_pattern'] for r in rows} == excluded,
                        'Exclusion bookkeeping does not match raw markers')
                continue
            for number, row in enumerate(rows, 2):
                for column, value in row.items():
                    require(value is not None, f'Malformed CSV row: {path.name}:{number}')
                    mentioned = set(RUN_ID.findall(value))
                    require(mentioned <= set(selected),
                            f'Out-of-scope run in {path.relative_to(ROOT)}:{number}:{column}: {mentioned - set(selected)}')
                    if column in LABEL_COLUMNS and value:
                        label = value.replace('_', ' ')
                        require(label in ALLOWED,
                                f'Excluded/unknown label in {path.name}:{number}:{column}: {value}')
                    if column in {'training_present_patterns', 'training_absent_patterns'} and value:
                        require(set(value.split(';')) <= set(ALLOWED), f'Unknown training class: {path.name}:{number}')
                rows_checked += 1
            checked.append(path.relative_to(ROOT).as_posix())
    for name in ALL_RUN_TABLES:
        require(set(by_run(name)) == set(selected), f'{name} does not cover exactly the selected runs')
    require({r['run'] for r in table('silhouette_pattern_features.csv')} == eligible,
            'Silhouette features do not cover exactly the statistical runs')
    for name, keys in (
        ('panorama_void_sensitivity.csv', ('variant',)),
        ('phase_mask_flip_sensitivity.csv', ('shadow_threshold',)),
        ('flip_rate_sensitivity.csv', ('rows', 'method')),
        ('flip_rate_spatial_sensitivity.csv', ('rows', 'method')),
    ):
        variants = defaultdict(set)
        for row in table(name):
            variants[tuple(row[key] for key in keys)].add(row['run'])
        require(variants and all(runs == set(selected) for runs in variants.values()),
                f'Incomplete/out-of-scope sensitivity variant in {name}')
    return {'active_csv_files_checked': len(checked), 'active_csv_rows_checked': rows_checked,
            'checked_files': checked, 'selected_acquisitions': len(selected),
            'statistical_acquisitions': len(eligible), 'excluded_acquisitions': excluded}


def interval_area(mask):
    """Match start/end transitions independently of either estimator recurrence."""
    transitions = np.diff(np.pad(mask.astype(np.int8), ((1, 1), (0, 0))), axis=0)
    start_y, start_x = np.nonzero(transitions == 1)
    end_y, end_x = np.nonzero(transitions == -1)
    start_order = np.lexsort((start_y, start_x))
    end_order = np.lexsort((end_y, end_x))
    require(np.array_equal(start_x[start_order], end_x[end_order]), 'Unmatched gas intervals')
    widths = end_y[end_order] - start_y[start_order]
    require((widths > 0).all(), 'Nonpositive interval width')
    return np.bincount(start_x[start_order], weights=widths.astype(float) ** 2 / 60 ** 2,
                       minlength=mask.shape[1])


def mask_arithmetic(selected, eligible, saturation):
    inventory, void, phase = (by_run(name) for name in
                              ('data_inventory.csv', 'panorama_void_per_run.csv', 'phase_mask_flip_per_run.csv'))
    old_inventory = {r['run']: r for r in read_csv(ARCHIVE / 'results/validation/data_inventory.csv')}
    mask_dir = artifact_path('phase_masks')
    require({p.name[:-9] for p in mask_dir.glob('*_mask.png')} == set(selected), 'Unexpected/missing active masks')
    errors = []
    for run, pattern in selected.items():
        source = DATA / run / '__panorama_v11.png'
        require(digest(source) == old_inventory[run]['sha256'] == inventory[run]['sha256'],
                f'Original panorama changed or hash stale: {run}')
        source_image = np.asarray(Image.open(source).convert('L'))
        mask_path = mask_dir / (run + '_mask.png')
        encoded = np.asarray(Image.open(mask_path).convert('L'))
        require(encoded.shape == source_image.shape and set(np.unique(encoded)) <= {0, 255},
                f'Invalid mask shape/encoding: {run}')
        mask = encoded == 255
        require(not mask[:18].any() and not mask[71:].any(), f'Mask outside defined crop: {run}')
        alpha = interval_area(mask)
        error = abs(float(alpha.mean()) - float(void[run]['conditional_circle_alpha_spatial_mean']))
        close(alpha.mean(), void[run]['conditional_circle_alpha_spatial_mean'], f'Independent area: {run}')
        require((alpha <= (53 / 60) ** 2 + 1e-12).all(), f'Crop ceiling violated: {run}')
        comparisons = 53 * (mask.shape[1] - 1)
        flips = int(np.count_nonzero(mask[18:71, 1:] ^ mask[18:71, :-1]))
        require(flips == int(phase[run]['candidate_pixel_flip_count']), f'Wrong mask flips: {run}')
        require(comparisons == int(phase[run]['candidate_pixel_flip_comparisons']), f'Wrong denominator: {run}')
        close(flips / comparisons, phase[run]['candidate_pixel_flip_rate'], f'Mask rate: {run}')
        center_flips = int(np.count_nonzero(mask[40, 1:] ^ mask[40, :-1]))
        require(center_flips == int(phase[run]['center_candidate_flip_count']), f'Wrong center flips: {run}')
        require(center_flips == int(phase[run]['center_candidate_gas_entries']) + int(phase[run]['center_candidate_gas_exits']),
                f'Entries/exits disagree with center flips: {run}')
        close(1000 * center_flips / (mask.shape[1] - 1), phase[run]['center_candidate_flips_per1000px'], f'Center scaling: {run}')
        close(mask[40].mean(), phase[run]['center_candidate_gas_duty'], f'Center occupancy: {run}')
        close(saturation[run], inventory[run]['core_saturation_fraction'], f'Saturation: {run}')
        require((phase[run]['eligible'] == 'True') == (run in eligible), f'Wrong eligibility: {run}')
        require(inventory[run]['pattern'] == phase[run]['pattern'] == void[run]['label'].replace('_', ' ') == pattern,
                f'Authoritative label changed: {run}')
        require(digest(mask_path) == digest(ARCHIVE / 'intermediate/phase_masks' / mask_path.name),
                f'Retained default mask changed after exclusions: {run}')
        errors.append(error)
    return {'input_panorama_hashes_verified': len(selected), 'masks_byte_identical_to_archive': len(selected),
            'independent_area_means_verified': len(errors), 'max_area_error': max(errors),
            'pixel_and_center_flip_counts_verified': len(selected),
            'saturation_excluded_selected_runs': sorted(set(selected) - eligible)}


def same_saved_value(current, old):
    if current == old:
        return True
    try:
        current, old = float(current), float(old)
    except (ValueError, TypeError):
        return False
    if math.isnan(current) and math.isnan(old):
        return True
    return math.isfinite(current) and math.isfinite(old) and abs(current - old) <= 1e-12


def retained_measurements(selected):
    checked = []
    for name in ('panorama_void_per_run.csv', 'phase_mask_flip_per_run.csv', 'flip_rate_per_run.csv',
                 'data_inventory.csv', 'image_void_alternatives.csv', 'pattern_morphology_features.csv'):
        path = artifact_path(name)
        archived_path = ARCHIVE / path.relative_to(ROOT)
        old = {r['run']: r for r in read_csv(archived_path)}
        current = by_run(name)
        cells = 0
        for run in selected:
            for column in current[run].keys() & old[run].keys():
                require(same_saved_value(current[run][column], old[run][column]),
                        f'Retained measurement changed: {name}:{run}:{column}')
                cells += 1
        checked.append({'artifact': path.relative_to(ROOT).as_posix(), 'common_cells_compared': cells,
                        'retained_runs': len(selected)})
    return {'tables_compared': checked, 'scope_and_classification_recalculations_are_intentional': True}


def calibration_and_metadata_inputs(selected, eligible):
    """Reject excluded run IDs in computational metadata, not scope bookkeeping."""
    checked = []
    bookkeeping = {'scope', 'dataset_scope', 'policy', 'excluded_runs', 'excluded_acquisitions',
                   'explicitly_excluded_acquisitions'}

    def inspect(value, location):
        if isinstance(value, dict):
            for key, item in value.items():
                if key in bookkeeping:
                    continue
                require(set(RUN_ID.findall(key)) <= set(selected), f'Excluded metadata key: {location}.{key}')
                inspect(item, f'{location}.{key}')
        elif isinstance(value, list):
            for index, item in enumerate(value):
                inspect(item, f'{location}[{index}]')
        elif isinstance(value, str):
            require(set(RUN_ID.findall(value)) <= set(selected), f'Excluded computational run: {location}')

    for name in ('panorama_liquid_reference.json', 'image_void_alternatives_metadata.json',
                 'phase_mask_flip_metadata.json', 'flip_rate_audit.json', 'legacy_independent_audit.json',
                 'pattern_identification.json', 'silhouette_pattern_results.json',
                 'adversarial_checks.json', 'final_verification.json', 'summary_chart_checks.json'):
        data = json.loads(artifact_path(name).read_text())
        inspect(data, name)
        checked.append(name)
    segmentation_hash = digest(artifact_path('panorama_void_alternatives.py'))
    for name, key in (('phase_mask_flip_metadata.json', 'segmentation_source_sha256'),
                      ('silhouette_pattern_results.json', 'segmentation_source_sha256'),
                      ('adversarial_checks.json', 'estimator_sha256'),
                      ('final_verification.json', 'segmentation_source_sha256')):
        data = json.loads(artifact_path(name).read_text())
        require(data[key] == segmentation_hash, f'Stale segmentation hash in {name}')
    phase = json.loads(artifact_path('phase_mask_flip_metadata.json').read_text())
    require(phase['n_saved_binary_masks'] == len(selected), 'Stale saved-mask count')
    chart = json.loads(artifact_path('summary_chart_checks.json').read_text())
    require(chart['checks']['joined_acquisitions'] == len(selected), 'Stale chart acquisition count')
    require(chart['checks']['unique_supplied_patterns'] == len(ALLOWED), 'Stale chart pattern count')
    require(chart['checks']['eligible_pure_unsaturated'] == len(eligible), 'Stale chart statistical count')
    return {'computational_json_inputs_checked': checked,
            'scope_bookkeeping_subtrees_exempt_from_input_scan': sorted(bookkeeping),
            'segmentation_source_sha256': segmentation_hash, 'all_four_source_hashes_current': True}


def grouped_and_chart_values(selected, eligible):
    void, phase = by_run('panorama_void_per_run.csv'), by_run('phase_mask_flip_per_run.csv')
    counts = Counter(selected[run] for run in eligible)
    void_groups, phase_groups = table('panorama_void_by_pattern.csv'), table('phase_mask_flip_by_pattern.csv')
    require({r['pattern'] for r in void_groups} == {r['pattern'] for r in phase_groups} == set(ALLOWED), 'Missing/extra group')
    for row in void_groups:
        pattern = row['pattern']
        values = [float(void[run]['conditional_circle_alpha_spatial_mean']) for run in eligible if selected[run] == pattern]
        require(int(row['n']) == counts[pattern], f'Wrong void group n: {pattern}')
        for key, value in (('conditional_alpha_median', np.median(values)), ('conditional_alpha_min', min(values)),
                           ('conditional_alpha_max', max(values))):
            close(row[key], value, f'Void group {pattern}:{key}')
    for row in phase_groups:
        pattern = row['pattern']
        require(int(row['n']) == counts[pattern], f'Wrong phase group n: {pattern}')
        for output, source in (('candidate_pixel_flip_rate_median', 'candidate_pixel_flip_rate'),
                               ('center_candidate_flip_rate_median', 'center_candidate_flip_rate'),
                               ('center_candidate_gas_duty_median', 'center_candidate_gas_duty')):
            close(row[output], np.median([float(phase[r][source]) for r in eligible if selected[r] == pattern]),
                  f'Phase group {pattern}:{output}')
    legacy = by_run('flip_rate_per_run.csv')
    for row in table('flip_rate_by_pattern.csv'):
        pattern, metric = row['pattern'], row['metric']
        values = [float(legacy[run][metric]) for run in eligible if selected[run] == pattern]
        require(int(row['n']) == counts[pattern], f'Wrong optical group n: {pattern}')
        for key, value in (('median', np.median(values)), ('min', min(values)), ('max', max(values)), ('mean', np.mean(values))):
            close(row[key], value, f'Optical group {pattern}:{metric}:{key}')
    points = by_run('void_fraction_flip_plot_data.csv')
    require(set(points) == set(selected), 'Wrong plotted acquisition set')
    for run, point in points.items():
        require(point['supplied_pattern'] == selected[run], f'Wrong plotted label: {run}')
        require((point['eligible_pure_unsaturated'] == 'True') == (run in eligible), f'Wrong plotted eligibility: {run}')
        close(point['conditional_circle_alpha_spatial_mean'], void[run]['conditional_circle_alpha_spatial_mean'], f'Chart x: {run}')
        close(point['candidate_pixel_flips_per1000px'], 1000 * float(phase[run]['candidate_pixel_flip_rate']), f'Chart y: {run}')
        close(point['center_candidate_flips_per1000px'], phase[run]['center_candidate_flips_per1000px'], f'Center chart y: {run}')
    meta = json.loads(artifact_path('summary_chart_checks.json').read_text())
    require(meta['passed'], 'Chart generator reported failed checks')
    for source in meta['source_inputs']:
        require(digest(ROOT / source['artifact']) == source['sha256'], f'Stale chart source: {source["artifact"]}')
    require(digest(artifact_path('summary_figures.py')) == meta['script_sha256'], 'Stale chart script provenance')
    for output in meta['outputs']:
        require((ROOT / output).is_file(), f'Missing chart export: {output}')
        if output.endswith('.svg'):
            content = (ROOT / output).read_text()
            require(not (set(RUN_ID.findall(content)) - set(selected)), f'Excluded ID in SVG: {output}')
            require('plug-slug' not in content and 'slug-churn' not in content and 'smooth annular' not in content,
                    f'Excluded pattern in SVG: {output}')
    for chart in ('primary_chart_hulls', 'center_chart_hulls'):
        require(set(meta[chart]) == set(ALLOWED), f'Wrong hull categories: {chart}')
        for group in meta[chart].values():
            require(set(group['eligible_runs']) <= eligible and set(group['polygon_vertex_runs']) <= set(group['eligible_runs']),
                    f'Invalid hull participants: {chart}')
    require(not meta['primary_chart_hulls']['dispersed bubbles']['shaded'], 'Invented singleton bubble area hull')
    return {'group_counts_independently_recomputed': dict(counts), 'chart_coordinate_pairs_verified': len(points),
            'flip_scaling_factor': 1000, 'chart_source_hashes_current': True, 'svg_scope_and_hull_scope_verified': True}


def classifier_audit(selected, eligible):
    controls = importlib.import_module('pattern_identification')
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        prediction = controls.predict(np.array([[0., 0.], [1., 1.], [8., 8.], [9., 9.]]),
                                      np.array([1, 1, 3, 3]), np.array([[.2, .2], [8.5, 8.5]]), 'centroid')
    require(set(prediction) <= {1, 3}, 'Centroid control predicted an absent training class')
    bubble_runs = {run for run in eligible if selected[run] == 'dispersed bubbles'}
    require(len(bubble_runs) == 1, 'The retained bubble support is not singleton')
    checked_models, prediction_rows, absent_counts = 0, 0, {}
    for metadata_name, csv_name in (('pattern_identification.json', 'pattern_identification_predictions.csv'),
                                    ('silhouette_pattern_results.json', 'silhouette_pattern_predictions.csv')):
        data = json.loads(artifact_path(metadata_name).read_text())
        grouped = defaultdict(list)
        absent_by_validation = defaultdict(set)
        for row in table(csv_name):
            run, validation = row['run'], row['validation']
            require(run in eligible and row['truth'] == selected[run], f'Wrong held-out truth: {csv_name}:{run}')
            require(validation in {'leave_one_run_out', 'leave_one_run_out_consistent_variant', 'leave_one_mass_flux_out'}, f'Unknown validation: {validation}')
            training = {r for r in eligible if r[:3] != run[:3]} if validation == 'leave_one_mass_flux_out' else eligible - {run}
            present = {selected[r] for r in training}
            absent = set(ALLOWED) - present
            if absent:
                absent_by_validation[validation].add(run)
            require(row['prediction'] in present, f'Prediction uses unavailable training class: {csv_name}:{run}')
            require(set(filter(None, row['training_present_patterns'].split(';'))) == present, f'Wrong present classes: {run}')
            require(set(filter(None, row['training_absent_patterns'].split(';'))) == absent, f'Wrong absent classes: {run}')
            require(int(row['n_training_runs']) == len(training), f'Wrong training n: {run}')
            require((row['test_pattern_has_training_examples'] == 'True') == (row['truth'] in present), f'Wrong support flag: {run}')
            require((row['correct'] == 'True') == (row['truth'] == row['prediction']), f'Wrong correctness flag: {run}')
            if 'reference_training_runs' in row:
                require(set(filter(None, row['reference_training_runs'].split(';'))) <= training, f'Reference leakage: {run}')
            grouped[(row['model'], validation, row.get('variant', 'baseline'))].append(row)
            prediction_rows += 1
        require(len(grouped) == len(data['results']), f'Missing/extra model results: {metadata_name}')
        for result in data['results']:
            key = (result['model'], result['validation'], result.get('variant', 'baseline'))
            records = grouped[key]
            require(len(records) == len(eligible) and {r['run'] for r in records} == eligible, f'Dropped/duplicated test case: {key}')
            matrix = np.zeros((len(ALLOWED), len(ALLOWED)), dtype=int)
            for row in records:
                matrix[ALLOWED.index(row['truth']), ALLOWED.index(row['prediction'])] += 1
            require(np.array_equal(matrix, result['confusion_true_rows_pred_columns']), f'Wrong confusion matrix: {key}')
            require((matrix.sum(axis=1) > 0).all(), f'Missing test class: {key}')
            close(result['accuracy'], np.trace(matrix) / matrix.sum(), f'Accuracy: {key}')
            close(result['balanced_accuracy'], np.mean(np.diag(matrix) / matrix.sum(axis=1)), f'Balanced accuracy: {key}')
            checked_models += 1
        require('absent_training_class_policy' in data or 'training_class_availability' in data or
                'fold_training_class_availability' in data, f'Missing class-support policy: {metadata_name}')
        observed = {validation: len(runs) for validation, runs in absent_by_validation.items()}
        if 'folds_with_absent_training_classes' in data:
            require(all(observed.get(key, 0) == count for key, count in data['folds_with_absent_training_classes'].items()),
                    f'Wrong reported absent-class counts in {metadata_name}')
        absent_counts[metadata_name] = observed
    folds = table('silhouette_pattern_fold_references.csv')
    require(len(folds) == 2 * len(eligible), 'Wrong number of held-out reference folds')
    for validation in ('leave_one_run_out', 'leave_one_mass_flux_out'):
        require({row['test_run'] for row in folds if row['validation'] == validation} == eligible,
                f'Missing held-out fold: {validation}')
    for row in folds:
        run = row['test_run']
        training = eligible - {run} if row['validation'] == 'leave_one_run_out' else {r for r in eligible if r[:3] != run[:3]}
        refs = set(filter(None, row['reference_training_runs'].split(';')))
        require(refs and refs <= training and run not in refs, f'Reference leakage: {run}')
    for name in ('pattern_identification.md', 'silhouette_pattern_review.md'):
        prose = artifact_path(name).read_text().lower()
        require('one' in prose and 'bubble' in prose and 'training' in prose and 'cannot' in prose,
                f'Singleton-class limitation missing from {name}')
    return {'model_variant_matrices_rebuilt': checked_models, 'heldout_predictions_checked': prediction_rows,
            'heldout_reference_folds_checked': len(folds), 'all_statistical_tests_retained_per_model': len(eligible),
            'singleton_bubble_runs': sorted(bubble_runs), 'absent_class_centroid_control_passed': True,
            'independently_detected_absent_training_class_folds': absent_counts}


def executive_summary_audit(selected, eligible):
    path = artifact_path('EXECUTIVE_SUMMARY.md')
    prose = path.read_text()
    void, phase = by_run('panorama_void_per_run.csv'), by_run('phase_mask_flip_per_run.csv')
    summary_rows = []
    for line in prose.splitlines():
        if not line.startswith('|'):
            continue
        cells = [cell.strip() for cell in line.strip('|').split('|')]
        if len(cells) != 4 or not cells[1].isdigit():
            continue
        label = cells[0].lower()
        pattern = 'dispersed bubbles' if label.startswith('dispersed bubbles') else label
        require(pattern in ALLOWED, f'Excluded/unknown executive-summary row: {label}')
        mentioned = set(RUN_ID.findall(label))
        require(mentioned <= set(selected), f'Excluded example in summary table: {mentioned}')
        runs = {run for run in eligible if selected[run] == pattern}
        require(int(cells[1]) == len(runs), f'Wrong summary group n: {pattern}')
        alpha = float(np.median([float(void[r]['conditional_circle_alpha_spatial_mean']) for r in runs]))
        flips = 1000 * float(np.median([float(phase[r]['candidate_pixel_flip_rate']) for r in runs]))
        for display, actual in ((cells[2], alpha), (cells[3], flips)):
            places = len(display.partition('.')[2])
            close(display, actual, f'Summary rounded value: {pattern}', tolerance=.50001 * 10 ** -places)
        summary_rows.append(pattern)
    require(len(summary_rows) == len(ALLOWED) and set(summary_rows) == set(ALLOWED), 'Missing/duplicated summary groups')
    normalized = prose.lower()
    require('both bubble' not in normalized and 'two bubble' not in normalized, 'Old two-bubble finding remains')
    require('bubble' in normalized and ('one' in normalized or 'single' in normalized), 'Singleton bubble limitation absent from summary')
    require('uncalibrated' in normalized or 'not calibrated' in normalized, 'Physical calibration limit absent')
    require('crop' in normalized and '0.7803' in normalized, 'Crop-ceiling limitation absent')
    require('1,000' in prose and ('not per second' in normalized or 'not a time frequency' in normalized), 'Spatial flip definition/scaling absent')
    require(not re.search(r'(?:all|kept)\s+\*{0,2}31\s+(?:acquisitions|panoramas)', normalized), 'Stale full-dataset count in summary')
    embedded = re.findall(r'!\[[^\]]*\]\(([^)]+)\)', prose)
    require(len(embedded) >= 3 and all((path.parent / href).is_file() for href in embedded), 'Missing summary graphics')
    return {'summary_table_groups_checked': summary_rows, 'summary_table_rounding_checked': True,
            'summary_flip_scaling_checked': 1000, 'embedded_figures_checked': len(embedded),
            'executive_summary_sha256': digest(path), 'limitations_and_singleton_support_disclosed': True}


def main():
    started = datetime.now(timezone.utc).isoformat()
    checks, failures = [], []
    try:
        selected, excluded, eligible, saturation = discover_scope()
    except Exception as exc:
        selected, excluded, eligible, saturation = {}, {}, set(), {}
        failures.append({'check': 'independent_raw_selection', 'error': f'{type(exc).__name__}: {exc}'})
    jobs = (
        ('scope_and_all_active_csvs', lambda: scan_csv_scope(selected, excluded, eligible)),
        ('independent_mask_arithmetic_and_input_hashes', lambda: mask_arithmetic(selected, eligible, saturation)),
        ('retained_measurements_against_archive', lambda: retained_measurements(selected)),
        ('calibration_and_computational_metadata_scope', lambda: calibration_and_metadata_inputs(selected, eligible)),
        ('groups_charts_scaling_and_provenance', lambda: grouped_and_chart_values(selected, eligible)),
        ('heldout_classifier_tests_and_missing_classes', lambda: classifier_audit(selected, eligible)),
        ('written_executive_summary', lambda: executive_summary_audit(selected, eligible)),
    )
    if selected:
        for name, job in jobs:
            try:
                details = job()
                checks.append({'check': name, 'passed': True, 'details': details})
            except Exception as exc:
                record = {'check': name, 'passed': False, 'error': f'{type(exc).__name__}: {exc}'}
                checks.append(record)
                failures.append(record)
    result = {
        'passed': not failures, 'started_utc': started,
        'finished_utc': datetime.now(timezone.utc).isoformat(),
        'script_sha256': digest(Path(__file__)),
        'independent_selection_uses_original_markers_not_shared_helper': True,
        'selected_acquisitions': len(selected), 'statistical_acquisitions': len(eligible),
        'supplied_pattern_counts_selected': dict(Counter(selected.values())),
        'selected_runs': sorted(selected), 'excluded_runs': excluded,
        'checks': checks, 'failures': failures,
        'scope': 'Post-regeneration arithmetic, dataset selection, leakage, chart and written-findings audit; not independent physical void calibration.',
    }
    artifact_path('post_regeneration_adversarial_checks.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    lines = [
        '# Independent adversarial review after regeneration', '',
        f"**Result: {'PASS' if result['passed'] else 'FAIL'}.** This check ran after the analysis and executive summary were written.", '',
        f"Independent reading of original folder labels and pixels finds **{len(selected)} selected acquisitions and {len(eligible)} unsaturated statistical acquisitions**, across the five approved single-pattern classes. Excluded data are allowed only in explicit exclusion bookkeeping and historical outputs.", '',
        '| Check | Result |', '|---|---|',
    ]
    lines.extend(f"| {item['check'].replace('_', ' ')} | {'PASS' if item['passed'] else 'FAIL'} |" for item in checks)
    lines += ['', '## Findings', '',
              'The audit tests selection consistency, saved mask arithmetic, retained default measurements, chart coordinates and scaling, held-out predictions and the numbers written in the executive summary. Excluded runs must stay out of current tables, references and classifier tests.', '',
              '**Limit:** a single retained bubble acquisition cannot validate a transferable bubble classifier. Its unsupported held-out test remains in the scores; absent training classes receive no manufactured centroid.', '',
              '**Limit:** passing these checks establishes computational and scope consistency, not correct three-dimensional gas shape or calibrated physical gas volume. The crop and optical-mapping limitations remain.', '',
              'Detailed observations and failures are recorded in [post_regeneration_adversarial_checks.json](post_regeneration_adversarial_checks.json).', '']
    if failures:
        lines += ['## Failures to resolve', '']
        lines.extend(f"- {item['check']}: {item['error']}" for item in failures)
    write_markdown(artifact_path('post_regeneration_adversarial_review.md'), '\n'.join(lines) + '\n')
    print(json.dumps({'passed': result['passed'], 'completed_checks': len(checks),
                      'selected_acquisitions': len(selected), 'statistical_acquisitions': len(eligible),
                      'failures': failures}, indent=2))
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
