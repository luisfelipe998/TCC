#!/usr/bin/env python3
"""Verify the current selected scope, retained measurements and navigation.

Uses archived pre-exclusion measurements to check retained default values.
Group statistics and classifier scores are intentionally recomputed for the selected scope.
"""
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit

from v14_paths import ROOT, SCRIPT_DIR, LAYOUT, artifact_path, write_markdown

ROOT_ENTRIES = {"README.md", "findings", "results", "figures",
                "scripts", "intermediate", "literature"}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def markdown_links():
    """Check local inline/image and reference-definition link destinations."""
    checked, failures = [], []
    inline = re.compile(r"!?\[[^\]\n]*\]\((<[^>\n]+>|[^\s)]+)(?:\s+[\"'][^\n]*?[\"'])?\)")
    reference = re.compile(r"^ {0,3}\[[^\]\n]+\]:\s*(<[^>\n]+>|\S+)", re.MULTILINE)
    for document in sorted(ROOT.rglob("*.md")):
        content = document.read_text()
        for pattern in (inline, reference):
            for match in pattern.finditer(content):
                target = match.group(1).strip("<>")
                parsed = urlsplit(target)
                if parsed.scheme or parsed.netloc or not parsed.path:
                    continue
                path = document.parent / unquote(parsed.path)
                item = dict(document=document.relative_to(ROOT).as_posix(),
                            line=content[:match.start()].count("\n") + 1,
                            target=target, exists=path.exists())
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
        failures.extend({'kind': 'retained measurement changed', 'artifact': relative, **item}
                        for item in differences)
        measurement_checks.append({'artifact': relative, 'retained_rows': len(current),
                                   'differences': differences})
    masks = ROOT / 'intermediate/phase_masks'
    actual_masks = {f.name.removesuffix('_mask.png') for f in masks.glob('*_mask.png')}
    if actual_masks != expected:
        failures.append({'kind': 'mask scope mismatch', 'actual': sorted(actual_masks)})
    identical_masks = 0
    for run in sorted(actual_masks & expected):
        current = masks / f'{run}_mask.png'
        prior = history / 'intermediate/phase_masks' / current.name
        if sha256(current) != sha256(prior):
            failures.append({'kind': 'retained mask changed', 'run': run})
        else:
            identical_masks += 1
    source_hash = sha256(artifact_path('panorama_void_alternatives.py'))
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
              'current_segmentation_source_sha256': source_hash,
              'source_hashes': current_hashes, 'failures': failures,
              'interpretation': 'Dataset scope, retained default measurements and navigation consistency; physical calibration is not established.'}
    artifact_path('scope_verification.json').write_text(json.dumps(result, indent=2) + '\n')
    lines = ['# Current scope and artifact verification', '',
             f"**{'PASS' if result['passed'] else 'FAIL'}:** current scope contains {len(expected)} acquisitions across five single-pattern classes. Only wavy annular represents annular flow.", '',
             f"Compared {len(measurement_checks)} retained per-run measurement tables against the [previous outputs](intermediate/history/pre_scope_exclusions/README.md). Retained default measurements and {identical_masks} masks match the archive. Group statistics, classifier fits and reference perturbations are recomputed for the new selection.", '',
             f"Checked {len(csv_checks)} active CSV artifacts for excluded runs, labels and reference inputs. Current source hashes match the segmentation script.", '',
             'Exact checks: [scope_verification.json](scope_verification.json). Independent review after the summary: [post-regeneration adversarial review](post_regeneration_adversarial_review.md).', '',
             'These checks verify exclusion consistency and arithmetic preservation, not absolute gas-volume calibration.', '',
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
                      'active_csv_artifacts_checked': len(csv_checks),
                      'local_links_checked': len(links), 'failures': failures}, indent=2))
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
