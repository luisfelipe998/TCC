#!/usr/bin/env python3
"""Independently verify selected masks, folds, scores and input panorama hashes."""
from collections import Counter, defaultdict
import csv
import hashlib
import json

import numpy as np
from PIL import Image

from panorama_void_alternatives import DIAMETER_PX, REFERENCE_RUNS
from v14_paths import DATA, artifact_path
from v14_selection import selected_runs, PATTERN_ORDER, scope_metadata


def table(filename):
    with artifact_path(filename).open() as fp:
        return list(csv.DictReader(fp))


def independent_area(mask):
    # Each active radial streak contributes odd increments whose sum is d².
    # This is independent of the estimator's flush-on-run-end arithmetic.
    streak = np.zeros(mask.shape[1], dtype=int)
    result = np.zeros(mask.shape[1], dtype=float)
    for active in mask:
        streak = np.where(active, streak + 1, 0)
        result += np.where(active, 2 * streak - 1, 0)
    return result / DIAMETER_PX ** 2


def split_patterns(value):
    return set(value.split(';')) if value else set()


def expected_training(test_run, validation, eligible):
    if validation == 'leave_one_mass_flux_out':
        return {run for run in eligible if run[:3] != test_run[:3]}
    return eligible - {test_run}


def main():
    selected = set(selected_runs())
    inventory = table('data_inventory.csv')
    inventory_by_run = {r['run']: r for r in inventory}
    void = {r['run']: r for r in table('panorama_void_per_run.csv')}
    phase = {r['run']: r for r in table('phase_mask_flip_per_run.csv')}
    assert len(inventory) == len(inventory_by_run)
    assert set(inventory_by_run) == set(void) == set(phase) == selected
    masks = artifact_path('phase_masks')
    assert {path.stem.removesuffix('_mask') for path in masks.glob('*_mask.png')} == selected
    eligible = {run for run, row in inventory_by_run.items() if row['pure_unsaturated'] == 'True'}
    assert eligible == {run for run, row in phase.items() if row['eligible'] == 'True'}
    expected_counts = dict(Counter(inventory_by_run[run]['pattern'] for run in eligible))
    assert set(expected_counts) == set(PATTERN_ORDER)
    errors = []
    for item in inventory:
        run = item['run']
        panorama_path = DATA / run / '__panorama_v11.png'
        assert hashlib.sha256(panorama_path.read_bytes()).hexdigest() == item['sha256'], run
        original = np.asarray(Image.open(panorama_path).convert('L'))
        saved = np.asarray(Image.open(masks / (run + '_mask.png')).convert('L'))
        assert saved.shape == original.shape, run
        assert set(np.unique(saved)) <= {0, 255}
        mask = saved == 255
        assert not mask[:18].any() and not mask[71:].any()
        area = independent_area(mask)
        assert np.all((area >= 0) & (area <= (53 / DIAMETER_PX) ** 2 + 1e-15))
        error = abs(float(area.mean()) - float(void[run]['conditional_circle_alpha_spatial_mean']))
        errors.append(error)
        assert error < 1e-12, (run, error)
        assert item['pattern'] == void[run]['label'].replace('_', ' ') == phase[run]['pattern']
        assert item['pattern'] in PATTERN_ORDER
        count = int(np.count_nonzero(mask[18:71, 1:] != mask[18:71, :-1]))
        assert count == int(phase[run]['candidate_pixel_flip_count'])
        denominator = 53 * (original.shape[1] - 1)
        assert denominator == int(phase[run]['candidate_pixel_flip_comparisons'])
        assert abs(count / denominator - float(phase[run]['candidate_pixel_flip_rate'])) < 1e-15
        center_count = int(np.count_nonzero(mask[40, 1:] != mask[40, :-1]))
        assert center_count == int(phase[run]['center_candidate_flip_count'])
        assert center_count == int(phase[run]['center_candidate_gas_entries']) + int(phase[run]['center_candidate_gas_exits'])
        assert abs(float(mask[40].mean()) - float(void[run]['center_gas_silhouette_duty'])) < 1e-15
        assert abs(float(mask[40].mean()) - float(phase[run]['center_candidate_gas_duty'])) < 1e-15

    source_hash = hashlib.sha256(artifact_path('panorama_void_alternatives.py').read_bytes()).hexdigest()
    phase_meta = json.loads(artifact_path('phase_mask_flip_metadata.json').read_text())
    classifier = json.loads(artifact_path('silhouette_pattern_results.json').read_text())
    adversarial = json.loads(artifact_path('adversarial_checks.json').read_text())
    assert source_hash == phase_meta['segmentation_source_sha256'] == classifier['segmentation_source_sha256'] == adversarial['estimator_sha256']
    assert phase_meta['n_saved_binary_masks'] == len(selected)
    scores = json.loads(artifact_path('pattern_identification.json').read_text())
    assert scores['n_eligible_runs'] == classifier['n_runs'] == len(eligible)
    assert scores['n_selected_runs'] == len(selected)
    assert scores['class_counts'] == classifier['class_counts'] == expected_counts
    assert set(scores['included_runs']) == set(classifier['included_runs']) == eligible
    assert {row['run'] for row in table('pattern_morphology_features.csv')} == selected
    silhouette_features = table('silhouette_pattern_features.csv')
    assert {row['run'] for row in silhouette_features} == eligible
    for row in silhouette_features:
        assert row['pattern'] == inventory_by_run[row['run']]['pattern']
        assert set(row['reference_training_runs'].split(';')) <= set(REFERENCE_RUNS) & eligible

    folds = table('silhouette_pattern_fold_references.csv')
    assert len(folds) == 2 * len(eligible)
    for validation in ('leave_one_run_out', 'leave_one_mass_flux_out'):
        subset = [row for row in folds if row['validation'] == validation]
        assert len(subset) == len(eligible)
        assert {row['test_run'] for row in subset} == eligible
    missing_folds = Counter()
    for fold in folds:
        run = fold['test_run']
        training = expected_training(run, fold['validation'], eligible)
        refs = set(fold['reference_training_runs'].split(';'))
        assert refs and refs <= training
        assert refs == set(REFERENCE_RUNS) & training
        assert run not in refs
        assert fold['heldout_run_excluded_from_reference'] == 'True'
        assert int(fold['n_training_runs']) == len(training)
        present = {inventory_by_run[name]['pattern'] for name in training}
        absent = set(PATTERN_ORDER) - present
        assert split_patterns(fold['training_present_patterns']) == present
        assert split_patterns(fold['training_absent_patterns']) == absent
        assert (fold['test_pattern_has_training_examples'] == 'True') == (inventory_by_run[run]['pattern'] in present)
        missing_folds[fold['validation']] += bool(absent)
    assert dict(missing_folds) == classifier['folds_with_absent_training_classes']

    prediction_checks = 0
    for collection, filename in ((scores['results'], 'pattern_identification_predictions.csv'),
                                 (classifier['results'], 'silhouette_pattern_predictions.csv')):
        predictions = table(filename)
        grouped = defaultdict(list)
        for row in predictions:
            run = row['run']
            assert run in eligible
            assert row['truth'] == inventory_by_run[run]['pattern']
            assert row['prediction'] in PATTERN_ORDER
            training = expected_training(run, row['validation'], eligible)
            present = {inventory_by_run[name]['pattern'] for name in training}
            assert row['prediction'] in present, 'prediction of an absent training class'
            assert split_patterns(row['training_present_patterns']) == present
            assert split_patterns(row['training_absent_patterns']) == set(PATTERN_ORDER) - present
            assert int(row['n_training_runs']) == len(training)
            assert (row['test_pattern_has_training_examples'] == 'True') == (row['truth'] in present)
            assert (row['correct'] == 'True') == (row['truth'] == row['prediction'])
            key = (row['model'], row['validation'], row.get('variant', 'baseline'))
            grouped[key].append(row)
            prediction_checks += 1
        assert len(grouped) == len(collection)
        for result in collection:
            key = (result['model'], result['validation'], result.get('variant', 'baseline'))
            rows = grouped[key]
            assert len(rows) == len(eligible) and {r['run'] for r in rows} == eligible
            independent_confusion = np.zeros((len(PATTERN_ORDER), len(PATTERN_ORDER)), dtype=int)
            for row in rows:
                independent_confusion[PATTERN_ORDER.index(row['truth']), PATTERN_ORDER.index(row['prediction'])] += 1
            confusion = np.asarray(result['confusion_true_rows_pred_columns'])
            assert np.array_equal(confusion, independent_confusion)
            assert confusion.sum() == len(eligible)
            assert abs(np.trace(confusion) / confusion.sum() - result['accuracy']) < 1e-15
            recall = np.diag(confusion) / confusion.sum(axis=1)
            assert np.isfinite(recall).all()
            assert abs(recall.mean() - result['balanced_accuracy']) < 1e-15
            assert all(abs(float(value) - result['recall_by_pattern'][pattern]) < 1e-15
                       for pattern, value in zip(PATTERN_ORDER, recall))

    result = dict(
        selected_acquisitions=len(selected), eligible_classifier_acquisitions=len(eligible),
        scope=scope_metadata(),
        unchanged_input_panorama_hashes=len(selected), saved_binary_masks_verified=len(selected),
        independently_recomputed_visible_area_means=len(selected), max_visible_area_absolute_error=max(errors),
        independently_recomputed_candidate_flip_counts=len(selected),
        authoritative_labels_preserved=True, source_hash_matches_phase_and_classifier=True,
        source_hash_matches_adversarial=True,
        heldout_reference_exclusions_verified=len(folds),
        heldout_training_class_availability_verified_predictions=prediction_checks,
        folds_with_absent_training_classes=dict(missing_folds),
        confusion_matrices_recomputed_from_per_run_predictions=True,
        confusion_matrix_scores_recomputed=True,
        no_out_of_scope_results_or_masks=True,
        selected_feature_scope_verified=True,
        singleton_class_limit='The only dispersed-bubble test has no bubble training example when held out; these data cannot validate a transferable bubble classifier.',
        segmentation_source_sha256=source_hash,
        interpretation='Selected-scope arithmetic and artifact consistency checks; not physical calibration.',
    )
    artifact_path('final_verification.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
