#!/usr/bin/env python3
"""Independent inventory and arithmetic reproduction from saved panoramas."""
from v14_paths import panorama_path, original_panorama_path, BASELINE
from v14_paths import ROOT as HERE, DATA, CACHE, artifact_path, write_markdown
from pathlib import Path
from v14_selection import iter_runs, selected_runs, PATTERN_ORDER, scope_metadata, is_eligible
from collections import Counter
import csv
import hashlib
import json

import numpy as np
from PIL import Image

PURE = set(PATTERN_ORDER)


def independent_otsu(a):
    histogram = np.bincount(a.ravel(), minlength=256).astype(float)
    weights = np.cumsum(histogram)
    moments = np.cumsum(histogram * np.arange(256))
    total = weights[-1]
    valid = (weights > 0) & (weights < total)
    between = np.full(256, -1.0)
    between[valid] = ((moments[-1] * weights[valid] - moments[valid] * total) ** 2 /
                      (weights[valid] * (total - weights[valid])))
    return int(np.argmax(between))


def main():
    with (artifact_path('flip_rate_per_run.csv')).open() as fp:
        saved = {r['run']: r for r in csv.DictReader(fp)}
    inventory, errors = [], []
    for folder in iter_runs():
        if not folder.is_dir():
            continue
        markers = list(folder.glob('__*.txt'))
        assert len(markers) == 1, folder
        pattern = markers[0].stem[2:].replace('_', ' ')
        path = panorama_path(folder)
        image = np.asarray(Image.open(path).convert('L'))
        core = image[20:60]
        threshold = independent_otsu(core)
        binary = core <= threshold
        count = int(np.count_nonzero(binary[:, 1:] != binary[:, :-1]))
        denominator = int(core.shape[0] * (core.shape[1] - 1))
        rate = count / denominator
        sat = float(np.mean(core == 255))
        original = saved[folder.name]
        assert count == int(original['spatial_flip_count']), folder.name
        assert denominator == int(original['spatial_n_comparisons'])
        errors.append(abs(rate - float(original['spatial_flip_rate'])))
        assert count == np.count_nonzero((~binary)[:, 1:] != (~binary)[:, :-1])
        assert count == np.count_nonzero(binary[:, ::-1][:, 1:] != binary[:, ::-1][:, :-1])
        inventory.append(dict(
            run=folder.name, pattern=pattern, panorama_height=int(image.shape[0]),
            panorama_width=int(image.shape[1]), core_saturation_fraction=sat,
            pure_pattern=pattern in PURE, pure_unsaturated=is_eligible(folder.name, pattern, sat),
            threshold=threshold, flip_count_all_40_rows=count,
            adjacent_comparisons_all_40_rows=denominator, flip_rate=rate,
            sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        ))
    assert len(inventory) == len(saved) == len(selected_runs())
    assert max(errors) < 1e-15
    selected = [r for r in inventory if r['pure_unsaturated']]
    result = dict(
        all_runs=len(inventory), pattern_counts=dict(Counter(r['pattern'] for r in inventory)),
        pure_unsaturated_runs=len(selected),
        pure_unsaturated_counts=dict(Counter(r['pattern'] for r in selected)),
        independently_reproduced_runs=len(errors), max_flip_rate_absolute_error=max(errors),
        independent_threshold_implementation='vectorized cumulative histogram variance',
        complement_and_spatial_reversal_invariance=f'passed on all {len(inventory)} selected panoramas',
        dataset_scope=scope_metadata(),
        data_role='original folder markers retained; no relabeling; panorama read-only',
    )
    with (artifact_path('data_inventory.csv')).open('w', newline='') as fp:
        writer = csv.DictWriter(fp, fieldnames=list(inventory[0]))
        writer.writeheader(); writer.writerows(inventory)
    (artifact_path('legacy_independent_audit.json')).write_text(json.dumps(result, indent=2) + '\n')
    scope = scope_metadata()
    artifact_path('dataset_selection.json').write_text(json.dumps(scope, indent=2) + '\n')
    excluded = scope['excluded_acquisitions']
    with artifact_path('dataset_exclusions.csv').open('w', newline='') as fp:
        writer = csv.DictWriter(fp, fieldnames=['run', 'supplied_pattern', 'reason'])
        writer.writeheader(); writer.writerows(excluded)
    report = [
        '# Independent check under the current dataset selection', '',
        f"The current analysis uses {len(inventory)} acquisitions across {len(PATTERN_ORDER)} single-pattern classes. Mixed labels, smooth annular and m05q05_002 are excluded before computation; only wavy annular represents annular flow. The original input files and their labels are retained.", '',
        f"Independent Otsu arithmetic reproduces all {len(errors)} retained optical flip counts and denominators. Maximum absolute rate difference: {max(errors):.6g}.", '',
        f"The existing >20% saturated-core exclusion leaves {len(selected)} statistical acquisitions. The selected per-pattern counts are {dict(Counter(r['pattern'] for r in selected))}. Only m03q40_001 is saturation-unsupported within the retained scope.", '',
        'The [selection manifest](dataset_selection.json) records the policy and retained runs. [Excluded acquisitions](dataset_exclusions.csv) are bookkeeping only and do not contribute to estimates, references, classifiers or charts.', '',
        'Previous wider-scope results are archived under intermediate/history/v14_baseline/intermediate/history/pre_scope_exclusions. Their counts and statistics are historical.', '',
        'The original file conflates two definitions. v11.1 counts Schmitt-trigger ',
        'switches of a one-dimensional dark occupancy. v11.3 counts individual binary ',
        'pixel changes across neighboring panorama columns, averaged over 40 rows. ',
        'The v11.3 values near 0.02 are fractions of adjacent comparisons, ',
        'not automatically switches per 1000 px or temporal frequencies.', '',
        'Complementing the optical binary mask preserves its flip count exactly. ',
        'Thus swapping dark and bright phase names cannot correct or validate flip-rate. ',
        'Both complement and spatial-reversal invariance pass for every retained image.', '',
        'The one-view/scale argument also needs correction: an axisymmetric physical ',
        'gas diameter can yield alpha=(d_g/D)^2 with absolute scale canceling. ',
        'Unvalidated optical shadows and unknown cross-sectional shapes are the actual ',
        'obstacles here. The legacy 60-pixel window is not a tube calibration; the current ',
        'primary estimator uses the user-confirmed full interior 18:71 and D=53 pixels.', '',
        'The data inventory records SHA-256 for every input panorama to make the source ',
        'images used in this check identifiable. The pipeline does not change them.', '',
        'Reproduce: `python3 src2/v14.1_new_voidfraction/scripts/inventory_and_legacy_audit.py`.',
    ]
    write_markdown(artifact_path('prior_findings_audit.md'), '\n'.join(report) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
