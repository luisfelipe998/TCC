#!/usr/bin/env python3
"""Verify the current selected scope, retained measurements and navigation.

Uses archived measurements to preserve inputs and unrelated observables.
Checks documented mask corrections against the independent method-change ledger.
"""
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit

from v14_paths import ROOT, SCRIPT_DIR, LAYOUT, artifact_path, write_markdown

ROOT_ENTRIES = {"README.md", "findings", "results", "figures",
                "scripts", "intermediate", "literature"}
FROZEN_LUMEN_SNAPSHOT = ROOT / 'intermediate/history/before_confirmed_lumen_D53'


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def markdown_links():
    """Check local inline/image and reference-definition link destinations."""
    checked, failures = [], []
    inline = re.compile(r"!?\[[^\]\n]*\]\((<[^>\n]+>|[^\s)]+)(?:\s+[\"'][^\n]*?[\"'])?\)")
    reference = re.compile(r"^ {0,3}\[[^\]\n]+\]:\s*(<[^>\n]+>|\S+)", re.MULTILINE)
    for document in sorted(ROOT.rglob("*.md")):
        content = document.read_text()
        # This exact-byte snapshot retains links relative to the original v14
        # location. Audit them in that original layout rather than rewriting
        # the frozen evidence. Current documents use their actual locations.
        frozen = document.is_relative_to(FROZEN_LUMEN_SNAPSHOT)
        link_base = ((ROOT / document.relative_to(FROZEN_LUMEN_SNAPSHOT)).parent
                     if frozen else document.parent)
        for pattern in (inline, reference):
            for match in pattern.finditer(content):
                target = match.group(1).strip("<>")
                parsed = urlsplit(target)
                if parsed.scheme or parsed.netloc or not parsed.path:
                    continue
                path = link_base / unquote(parsed.path)
                item = dict(document=document.relative_to(ROOT).as_posix(),
                            line=content[:match.start()].count("\n") + 1,
                            target=target, exists=path.exists(),
                            resolution_basis='original v14 layout for byte-frozen snapshot' if frozen else 'current document location')
                checked.append(item)
                if not item["exists"]:
                    failures.append(item)
    return checked, failures


def table(path):
    import csv
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def main():
    from v14_selection import selected_runs, PATTERN_ORDER, scope_metadata
    expected = set(selected_runs())
    history = ROOT / 'intermediate/history/pre_scope_exclusions'
    failures, measurement_checks = [], []
    snapshot_manifest = json.loads((FROZEN_LUMEN_SNAPSHOT / 'snapshot_manifest.json').read_text())
    for relative, expected_hash in snapshot_manifest['files'].items():
        frozen_path = FROZEN_LUMEN_SNAPSHOT / relative
        if not frozen_path.is_file() or sha256(frozen_path) != expected_hash:
            failures.append({'kind': 'frozen D60 snapshot changed', 'artifact': relative})
    names = ['results/void_fraction/panorama_void_per_run.csv',
             'results/flip_rate/phase_mask_flip_per_run.csv',
             'results/flip_rate/flip_rate_per_run.csv',
             'results/validation/data_inventory.csv',
             'intermediate/exploratory/image_void_alternatives.csv']
    for relative in names:
        current = {r['run']: r for r in table(ROOT / relative)}
        prior = {r['run']: r for r in table(history / relative)}
        if set(current) != expected:
            failures.append({'kind': 'measurement scope mismatch', 'artifact': relative,
                             'actual': sorted(current), 'expected': sorted(expected)})
        differences = []
        for run in sorted(set(current) & expected):
            if run not in prior:
                differences.append({'run': run, 'field': 'missing historical row'})
                continue
            for field in current[run].keys() & prior[run].keys():
                if current[run][field] != prior[run][field]:
                    differences.append({'run': run, 'field': field,
                                        'before': prior[run][field], 'after': current[run][field]})
        derived = relative in {'results/void_fraction/panorama_void_per_run.csv',
                               'results/flip_rate/phase_mask_flip_per_run.csv'}
        identity = {'run', 'label', 'pattern', 'panorama_width_px', 'panorama_width',
                    'core_saturated_fraction', 'eligible', 'candidate_mask_rows',
                    'candidate_mask_n_rows', 'candidate_pixel_flip_comparisons', 'validity'}
        failures.extend({'kind': 'unrelated retained measurement changed', 'artifact': relative, **item}
                        for item in differences if not derived or item['field'] in identity)
        measurement_checks.append({'artifact': relative, 'retained_rows': len(current),
                                   'differences': differences})
    masks = ROOT / 'intermediate/phase_masks'
    actual_masks = {f.name.removesuffix('_mask.png') for f in masks.glob('*_mask.png')}
    if actual_masks != expected:
        failures.append({'kind': 'mask scope mismatch', 'actual': sorted(actual_masks)})
    identical_masks = 0
    changed_masks = 0
    from PIL import Image
    import numpy as np
    method_archive = ROOT / 'intermediate/history/before_deep_method_review'
    corrections = {r['run']: r for r in table(artifact_path('deep_method_changes.csv'))}
    if set(corrections) != expected:
        failures.append({'kind': 'method-correction ledger scope mismatch'})
    for run in sorted(actual_masks & expected):
        current = masks / f'{run}_mask.png'
        prior = method_archive / 'intermediate/phase_masks' / current.name
        delta = corrections.get(run, {})
        before = np.asarray(Image.open(prior)) > 0
        after = np.asarray(Image.open(current)) > 0
        changed = int(np.count_nonzero(before != after))
        if (delta.get('before_mask_sha256') != sha256(prior)
                or delta.get('after_mask_sha256') != sha256(current)
                or int(delta.get('changed_pixels', -1)) != changed):
            failures.append({'kind': 'mask correction does not match independent ledger', 'run': run})
        if changed:
            changed_masks += 1
        else:
            identical_masks += 1
    source_hash = sha256(artifact_path('panorama_void_alternatives.py'))
    method = json.loads(artifact_path('deep_method_checks.json').read_text())
    if not method['passed'] or method['segmentation_source_sha256'] != source_hash:
        failures.append({'kind': 'method verification failed or source hash stale'})
    if method.get('morphology_comparison_diameter_px') != 60:
        failures.append({'kind': 'historical morphology ledger must compare both masks at D60'})
    from lumen_normalization_checks import verify_normalization
    normalization, _ = verify_normalization(expected_runs=expected)
    recorded_normalization = json.loads(artifact_path('lumen_normalization_checks.json').read_text())
    if (not recorded_normalization['passed']
            or recorded_normalization['segmentation_source_sha256'] != source_hash
            or recorded_normalization['geometry_source_sha256'] != sha256(artifact_path('v14_geometry.py'))):
        failures.append({'kind': 'lumen normalization failed or source provenance stale'})
    current_hashes = {}
    for name, field in [('phase_mask_flip_metadata.json', 'segmentation_source_sha256'),
                        ('silhouette_pattern_results.json', 'segmentation_source_sha256'),
                        ('adversarial_checks.json', 'estimator_sha256'),
                        ('final_verification.json', 'segmentation_source_sha256')]:
        data = json.loads(artifact_path(name).read_text())
        current_hashes[name] = data[field]
        if data[field] != source_hash:
            failures.append({'kind': 'stale source hash', 'artifact': name})
    # Group statistics and classifiers change intentionally; defaults for retained
    # acquisitions must match the archived measurements, rather than old group totals.
    csv_checks = []
    csv_files = list((ROOT / 'results').rglob('*.csv'))
    for folder in ['features', 'sensitivity', 'exploratory']:
        csv_files += list((ROOT / 'intermediate' / folder).rglob('*.csv'))
    for path in sorted(csv_files):
        if path.name == 'dataset_exclusions.csv':
            continue
        records = table(path)
        leaks = []
        for row in records:
            for key, value in row.items():
                if key in {'run', 'test_run'} and value and value not in expected:
                    leaks.append({'field': key, 'value': value})
                if key in {'pattern', 'label', 'true_pattern', 'predicted_pattern',
                           'supplied_pattern', 'actual_pattern'} and value and value.replace('_', ' ') not in PATTERN_ORDER:
                    leaks.append({'field': key, 'value': value})
                # Run-reference columns end in _run/_runs. A Boolean such as
                # heldout_run_excluded_from_reference is not a list of IDs.
                if 'reference' in key and key.endswith(('_run', '_runs')) and value:
                    for run in value.split(';'):
                        if run and run not in expected:
                            leaks.append({'field': key, 'value': run})
        csv_checks.append({'artifact': path.relative_to(ROOT).as_posix(), 'rows': len(records), 'leaks': leaks})
        failures.extend({'kind': 'excluded CSV input leaked', 'artifact': path.relative_to(ROOT).as_posix(), **item} for item in leaks)
    allowed_scripts = {Path(v).name for v in LAYOUT.values() if Path(v).parent.as_posix() == 'scripts'}
    allowed_scripts |= {'README.md', 'artifact_layout.json'}
    unexpected = [f.relative_to(ROOT).as_posix() for f in ROOT.iterdir() if f.name not in ROOT_ENTRIES]
    unexpected += [f.relative_to(ROOT).as_posix() for f in SCRIPT_DIR.iterdir() if f.name not in allowed_scripts or not f.is_file()]
    failures.extend({'kind': 'unexpected output location', 'artifact': item} for item in unexpected)
    result = {'passed': not failures, 'scope': scope_metadata(),
              'measurement_checks': measurement_checks, 'active_csv_scope_checks': csv_checks,
              'retained_masks_byte_identical': identical_masks,
              'masks_changed_by_verified_method_correction': changed_masks,
              'masks_identical_after_D53_normalization': normalization['identical_masks'],
              'flip_metrics_unchanged_after_D53_normalization': normalization['unchanged_flip_metrics'],
              'confirmed_lumen_rows': [18,71], 'current_lumen_diameter_px': 53,
              'boundary_source': normalization['boundary_source'],
              'current_segmentation_source_sha256': source_hash,
              'frozen_D60_files_hash_checked': len(snapshot_manifest['files']),
              'frozen_snapshot_link_policy': 'Exact-byte snapshot links are checked in their original v14 layout; current document links are checked at their current paths.',
              'source_hashes': current_hashes, 'failures': failures,
              'interpretation': 'Dataset scope, historical morphology correction at D60, unchanged masks/flips during D53 normalization, input preservation and navigation; physical phase calibration is not established.'}
    artifact_path('scope_verification.json').write_text(json.dumps(result, indent=2) + '\n')
    lines = ['# Current scope and artifact verification', '',
             f"**{'PASS' if result['passed'] else 'FAIL'}:** current scope contains {len(expected)} acquisitions across five single-pattern classes. Only wavy annular represents annular flow.", '',
             f"Compared {len(measurement_checks)} per-run tables against the [earlier outputs](intermediate/history/pre_scope_exclusions/README.md). Input observables and labels are preserved. Against the [pre-review snapshot](intermediate/history/before_deep_method_review/README.md), {changed_masks} masks change under the independently verified border correction and {identical_masks} remain identical. Derived values, groups and classifiers are recomputed.", '',
             f"The user then confirmed the complete lumen as slice 18:71, D=53. Against the [pre-normalization snapshot](intermediate/history/before_confirmed_lumen_D53/README.md), all {normalization['identical_masks']} masks and {normalization['unchanged_flip_metrics']} sets of flip metrics remain identical. The alpha rescaling by (60/53)^2 is independently verified; the historical morphology comparison remains at D60 on both sides.", '',
             f"Checked {len(csv_checks)} active CSV artifacts for excluded runs, labels and reference inputs. Current source hashes match the segmentation script.", '',
             'Exact checks: [scope_verification.json](scope_verification.json). Independent review after the summary: [post-regeneration adversarial review](post_regeneration_adversarial_review.md).', '',
             'These checks verify exclusion consistency, documented corrections and arithmetic; absolute gas-volume calibration remains unresolved.', '',
             'The D60 snapshot preserves original bytes and relative links. Its file hashes are verified, and its links are evaluated in the original v14 layout. Current report links are checked at their current locations.', '',
             '## Failures', '']
    lines += [json.dumps(item, ensure_ascii=False) for item in failures] or ['None.']
    write_markdown(artifact_path('scope_verification.md'), '\n'.join(lines) + '\n')
    write_markdown(artifact_path('layout_verification.md'),
                   '# Current scope and navigation verification\n\n'
                   'The previous reorganization-only check is historical. The current check follows the user dataset exclusions.\n\n'
                   '[Current verification](scope_verification.md) · [Historical outputs](intermediate/history/pre_scope_exclusions/README.md)\n')
    links, link_failures = markdown_links()
    failures.extend({'kind': 'broken Markdown link', **item} for item in link_failures)
    result.update(passed=not failures, local_markdown_links=links,
                  local_link_count=len(links), broken_link_count=len(link_failures), failures=failures)
    for name in ['scope_verification.json', 'layout_verification.json']:
        artifact_path(name).write_text(json.dumps(result, indent=2) + '\n')
    if link_failures:
        doc = artifact_path('scope_verification.md')
        doc.write_text(doc.read_text().replace('**PASS:**', '**FAIL:**'))
        with doc.open('a') as fp:
            fp.write('\nLocal Markdown link failures: ' + json.dumps(link_failures) + '\n')
    print(json.dumps({'passed': result['passed'], 'selected_runs': len(expected),
                      'identical_retained_masks': identical_masks,
                      'masks_changed_by_verified_method_correction': changed_masks,
                      'active_csv_artifacts_checked': len(csv_checks),
                      'local_links_checked': len(links), 'failures': failures}, indent=2))
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
