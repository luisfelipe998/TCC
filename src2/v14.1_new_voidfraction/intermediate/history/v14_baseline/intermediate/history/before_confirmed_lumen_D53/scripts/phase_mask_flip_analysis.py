#!/usr/bin/env python3
"""Spatial flips of the reconstructed gas-silhouette mask, not raw shadows."""
from v14_paths import ROOT as HERE, DATA, CACHE, artifact_path, write_markdown
from v14_selection import iter_runs, is_eligible, PATTERNS
from pathlib import Path
import csv
import json
import hashlib

import numpy as np
from raw_flip_analysis import two_sided_mann_whitney
from panorama_void_alternatives import (DATA, HERE, make_reference, read, silhouette)

ORDER = PATTERNS


def flips(mask):
    mask = np.asarray(mask, dtype=bool)
    return int(np.count_nonzero(mask[..., 1:] != mask[..., :-1]))


def write_csv(name, rows):
    with (artifact_path(name)).open('w', newline='') as fp:
        writer = csv.DictWriter(fp, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main():
    segmentation_source = artifact_path('panorama_void_alternatives.py')
    segmentation_hash = hashlib.sha256(segmentation_source.read_bytes()).hexdigest()
    mask_folder = artifact_path('phase_masks')
    mask_folder.mkdir(exist_ok=True)
    reference, _ = make_reference()
    rows, variants = [], []
    assert flips([0, 1, 1, 0]) == 2
    assert flips([1, 1, 1, 1]) == 0
    folders = list(iter_runs())
    for folder in folders:
        im = read(folder.name)
        pattern = next(folder.glob('__*.txt')).stem[2:].replace('_', ' ')
        mask, accepted, _ = silhouette(im, reference)
        from PIL import Image
        Image.fromarray(mask.astype(np.uint8) * 255).save(mask_folder / (folder.name + '_mask.png'))
        core = mask[18:71]
        W = mask.shape[1]
        center = mask[40]
        count = flips(core)
        center_count = flips(center)
        assert count == flips(~core) == flips(core[:, ::-1])
        assert center_count == flips(~center) == flips(center[::-1])
        entries = int(np.count_nonzero(~center[:-1] & center[1:]))
        exits = int(np.count_nonzero(center[:-1] & ~center[1:]))
        assert entries + exits == center_count
        rows.append(dict(
            run=folder.name, pattern=pattern, panorama_width=W,
            candidate_mask_rows='18:71', candidate_mask_n_rows=53,
            candidate_pixel_flip_count=count,
            candidate_pixel_flip_comparisons=53 * (W - 1),
            candidate_pixel_flip_rate=count / (53 * (W - 1)),
            center_candidate_flip_count=center_count,
            center_candidate_flip_rate=center_count / (W - 1),
            center_candidate_flips_per1000px=center_count * 1000 / (W - 1),
            center_candidate_gas_entries=entries,
            center_candidate_gas_exits=exits,
            center_candidate_gas_duty=float(center.mean()),
            paired_shadow_detection_fraction=float(accepted.mean()),
            eligible=is_eligible(folder.name, pattern, float((im[20:60] == 255).mean())),
            phase_status='inferred gas silhouette; optical calibration conditional',
        ))
        for threshold in (.60, .65, .70, .75, .80):
            alt, _, _ = silhouette(im, reference, threshold=threshold)
            variants.append(dict(
                run=folder.name, pattern=pattern, shadow_threshold=threshold,
                candidate_pixel_flip_rate=flips(alt[18:71]) / (53 * (W - 1)),
                center_candidate_flip_rate=flips(alt[40]) / (W - 1),
                eligible=rows[-1]['eligible'],
            ))
    assert {r['run'] for r in rows} == {folder.name for folder in folders}
    assert hashlib.sha256(segmentation_source.read_bytes()).hexdigest() == segmentation_hash
    (artifact_path('phase_mask_flip_metadata.json')).write_text(json.dumps(dict(
        segmentation_source_sha256=segmentation_hash, n_saved_binary_masks=len(rows),
        mask_encoding='255=inferred gas silhouette, 0=unmarked; not independent phase truth',
        fixed_probe_row_zero_based=40,
        center_field_semantics='Legacy center_* names refer to fixed image row y=40; physical centering is unvalidated',
        source_data='src2/data_processed/*/__panorama_v11.png',
    ), indent=2) + '\n')
    groups = []
    for pattern in ORDER:
        subset = [r for r in rows if r['pattern'] == pattern and r['eligible']]
        groups.append(dict(
            pattern=pattern, n=len(subset),
            candidate_pixel_flip_rate_median=float(np.median([r['candidate_pixel_flip_rate'] for r in subset])),
            center_candidate_flip_rate_median=float(np.median([r['center_candidate_flip_rate'] for r in subset])),
            center_candidate_gas_duty_median=float(np.median([r['center_candidate_gas_duty'] for r in subset])),
        ))
    write_csv('phase_mask_flip_per_run.csv', rows)
    write_csv('phase_mask_flip_by_pattern.csv', groups)
    write_csv('phase_mask_flip_sensitivity.csv', variants)
    comparisons = []
    for threshold in (.60, .65, .70, .75, .80):
        subset = [r for r in variants if r['eligible'] and r['shadow_threshold'] == threshold]
        for metric in ('candidate_pixel_flip_rate', 'center_candidate_flip_rate'):
            churn = [r[metric] for r in subset if r['pattern'] == 'churn']
            annular = [r[metric] for r in subset if r['pattern'] == 'wavy annular']
            test, method_metadata = two_sided_mann_whitney(churn, annular, method='auto')
            comparisons.append(dict(
                shadow_threshold=threshold, metric=metric,
                n_churn=len(churn), n_annular=len(annular),
                churn_median=float(np.median(churn)), annular_median=float(np.median(annular)),
                exploratory_mwu_two_sided_p=float(test.pvalue),
                auc_churn_greater=float(test.statistic / (len(churn) * len(annular))),
                **method_metadata,
            ))
    write_csv('phase_mask_flip_group_sensitivity.csv', comparisons)
    lines = [
        '# Flips of the gas-silhouette candidate', '',
        f'The current scope contains {len(rows)} retained acquisitions and {sum(r["eligible"] for r in rows)} statistical acquisitions passing the saturation criterion across {len(ORDER)} patterns. Mixed labels, smooth annular and the excluded replicate are omitted before computing masks.', '',
        'These measurements use the same inferred phase mask as the panorama void ',
        'reconstruction. They are spatial changes in a gas-silhouette candidate, ',
        'whose gas labeling and radial crop remain conditional.', '',
        'Pixel flip rate = mask changes over all 53 rows / [53 × (W−1)]. ',
        'Fixed-row flip rate = changes in mask row y=40 / (W−1). Each complete isolated ',
        'gas segment along that row contributes two interfaces; entry count and flip count ',
        'are saved separately. Neither quantity is a temporal event frequency.', '',
        'The legacy center_* column names mean fixed image row y=40; this row has not been validated as the physical tube center.', '',
        '| Pattern | n | Candidate pixel flip rate | Row y=40 flip rate | Row y=40 gas duty |',
        '|---|---:|---:|---:|---:|',
    ]
    for r in groups:
        lines.append(f"| {r['pattern']} | {r['n']} | {r['candidate_pixel_flip_rate_median']:.6f} | "
                     f"{r['center_candidate_flip_rate_median']:.6f} | {r['center_candidate_gas_duty_median']:.4f} |")
    lines.extend([
        '', 'The original v11.3 optical flip rates use raw dark/bright classes and ',
        '40 rows; these candidate-mask values use reconstructed silhouettes and ',
        '53 rows. They measure different masks and must retain separate names.', '',
        'The same reconstruction that makes the annular core continuous also ',
        'suppresses internal brightness fluctuations. This is an intended phase ',
        'definition, but can suppress real liquid fragments inside churn gas ',
        'envelopes. It is not independent validation of a desired ranking.', '',
        'Thresholds 0.60–0.80 are all reported in the sensitivity CSV. Group tests ',
        'are exploratory comparisons with no multiple-comparison correction. ',
        'Mann–Whitney uses the exact method when the smaller group has at most eight acquisitions and there are no pooled ties; otherwise it uses the asymptotic method with tie and continuity corrections. Each saved comparison records the actual method, tie flag and both sample sizes. The rank test concerns distributions, not medians alone, and does not validate absolute phase calibration.', '',
        'Retained supplied labels are unchanged. The single retained dispersed-bubble acquisition does not establish within-pattern variability.', '',
        'Arithmetic checks reproduce two interfaces for a complete isolated ',
        'segment and zero for a constant gas record, and verify complement/reversal ',
        f'invariance for all {len(rows)} reconstructed masks.', '',
        'Reproduce: `python3 src2/v14_new_voidfraction/scripts/phase_mask_flip_analysis.py`.',
    ])
    write_markdown(artifact_path('phase_mask_flip_review.md'), '\n'.join(lines) + '\n')
    print(json.dumps(groups, indent=2))


if __name__ == '__main__':
    main()
