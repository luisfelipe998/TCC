"""Verify current v14.1 provenance, scope, required outputs and report links."""
import ast
import csv
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit

from v14_paths import ROOT, SCRIPT_DIR, BASELINE, LAYOUT, artifact_path, write_markdown
from v14_selection import selected_runs, PATTERN_ORDER, scope_metadata


def sha256(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def markdown_links():
    checked,failures=[],[]
    patterns=[re.compile(r'!?\[[^\]\n]*\]\((<[^>\n]+>|[^\s)]+)(?:\s+[\"\x27][^\n]*?[\"\x27])?\)'),
              re.compile(r'^ {0,3}\[[^\]\n]+\]:\s*(<[^>\n]+>|\S+)',re.MULTILINE)]
    for document in sorted(ROOT.rglob('*.md')):
        if document.is_relative_to(BASELINE):continue
        content=document.read_text()
        for pattern in patterns:
            for match in pattern.finditer(content):
                href=match.group(1).strip('<>');parts=urlsplit(href)
                if parts.scheme or parts.netloc or not parts.path:continue
                resolved=document.parent/unquote(parts.path)
                item=dict(document=document.relative_to(ROOT).as_posix(),target=href,exists=resolved.exists())
                checked.append(item)
                if not item['exists']:failures.append(item)
    return checked,failures


def main():
    failures=[];checks={};selected=set(selected_runs())
    def run_check(name,callback):
        try:checks[name]=dict(passed=True,details=callback())
        except Exception as exc:
            checks[name]=dict(passed=False,error=f'{type(exc).__name__}: {exc}')
            failures.append(dict(check=name,error=checks[name]['error']))
    def evidence():
        names=['input_reconstruction.json','deep_method_checks.json','lumen_normalization_checks.json',
               'final_verification.json','summary_chart_checks.json','fix_effectiveness.json',
               'seam_adversarial_checks.json','post_regeneration_adversarial_checks.json']
        for name in names:assert json.loads(artifact_path(name).read_text())['passed'],name
        for name,script in [('seam_adversarial_checks.json','seam_adversarial.py'),
                            ('post_regeneration_adversarial_checks.json','post_regeneration_adversarial.py')]:
            assert json.loads(artifact_path(name).read_text())['script_sha256']==sha256(artifact_path(script)),name
        source=sha256(artifact_path('panorama_void_alternatives.py'))
        for name,key in [('phase_mask_flip_metadata.json','segmentation_source_sha256'),
                         ('silhouette_pattern_results.json','segmentation_source_sha256'),
                         ('adversarial_checks.json','estimator_sha256'),('final_verification.json','segmentation_source_sha256')]:
            assert json.loads(artifact_path(name).read_text())[key]==source,name
        return names
    def frozen():
        manifest=json.loads((ROOT/'intermediate/history/v14_baseline_manifest.json').read_text())
        for relative,digest in manifest['files'].items():
            assert sha256(BASELINE/relative)==digest,relative
            assert sha256(Path(manifest['original_root'])/relative)==digest,relative
        source=json.loads(artifact_path('source_input_hashes.json').read_text())
        for relative,digest in source.items():assert sha256(ROOT.parent/relative)==digest,relative
        return dict(frozen_and_original_v14_files=len(manifest['files']),unchanged_source_files=len(source))
    def scope():
        tables=[]
        folders=[ROOT/'results',ROOT/'intermediate/features',ROOT/'intermediate/sensitivity',ROOT/'intermediate/exploratory']
        for folder in folders:
            for path in sorted(folder.rglob('*.csv')):
                if path.name=='dataset_exclusions.csv':continue
                with path.open() as stream:rows=list(csv.DictReader(stream))
                for row in rows:
                    for column,value in row.items():
                        assert value is not None,(path,column)
                        assert set(re.findall(r'm\d{2}q\d{2}_\d{3}',value))<=selected,(path,column,value)
                        if column in {'pattern','label','truth','prediction','supplied_pattern'} and value:
                            assert value.replace('_',' ') in PATTERN_ORDER,(path,column,value)
                tables.append(dict(path=path.relative_to(ROOT).as_posix(),rows=len(rows)))
        assert {p.stem.removesuffix('_mask') for p in artifact_path('phase_masks').glob('*_mask.png')}==selected
        assert {p.parent.name for p in (ROOT/'intermediate/reconstructed_inputs').glob('*/__panorama_v11.png')}==selected
        return tables
    def scripts_and_required_outputs():
        names=[]
        for path in SCRIPT_DIR.glob('*.py'):
            ast.parse(path.read_text(),filename=str(path));names.append(path.name)
            assert path.name in LAYOUT,path.name
        required=['EXECUTIVE_SUMMARY.md','EXECUTIVE_SUMMARY_PT_BR.md','fix_effectiveness.md',
                  'void_fraction_vs_flip_rate.png','void_fraction_vs_flip_rate.svg','threshold_sensitivity.png',
                  'fix_comparison.png','fix_comparison_per_run.csv','classifier_comparison.csv',
                  'threshold_hull_intersections.csv','source_input_hashes.json','reconstruction_shifts.csv']
        for name in required:assert artifact_path(name).is_file(),name
        expected={'README.md','scripts','results','findings','figures','intermediate','literature'}
        assert {p.name for p in ROOT.iterdir()}==expected
        return dict(parsed_scripts=names,required_artifacts=required)
    for name,job in [('verification_evidence',evidence),('original_and_baseline_preservation',frozen),
                     ('active_dataset_scope',scope),('scripts_and_required_outputs',scripts_and_required_outputs)]:
        run_check(name,job)
    for filename in ('scope_verification.md','layout_verification.md'):
        write_markdown(artifact_path(filename),'# v14.1 scope, provenance and navigation\n\n'
            'Checks cover all active result tables, current corrected inputs, fixed sample selection, analysis source hashes, '
            'original/source-file preservation and current report links. Frozen v14 snapshot documents retain original bytes; '
            'their historical links are not rewritten or interpreted as current results.\n\n'
            '[Machine-readable evidence](layout_verification.json) · [Adversarial correction audit](seam_adversarial_review.md) · '
            '[Independent result audit](post_regeneration_adversarial_review.md)\n')
    # Create the JSON destinations before checking links that point to them.
    for name in ('scope_verification.json','layout_verification.json'):
        artifact_path(name).write_text('{}\n')
    links,broken=markdown_links()
    failures.extend(dict(check='local_markdown_link',**r) for r in broken)
    result=dict(passed=not failures,scope=scope_metadata(),checks=checks,
                local_markdown_links=links,local_link_count=len(links),broken_link_count=len(broken),failures=failures,
                interpretation='Implementation, input preservation, scope and reporting checks; physical gas-volume accuracy is not established.')
    for name in ('scope_verification.json','layout_verification.json'):
        artifact_path(name).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    status='PASS' if result['passed'] else 'FAIL'
    for name in ('scope_verification.md','layout_verification.md'):
        path=artifact_path(name);path.write_text(path.read_text()+'\n**'+status+'** — '+str(len(links))+' local links checked.\n')
    manifest={}
    for path in sorted(ROOT.rglob('*')):
        if not path.is_file() or path.is_relative_to(BASELINE):continue
        relative=path.relative_to(ROOT).as_posix()
        if relative.startswith(('intermediate/cache/','intermediate/logs/')) or path.name=='artifact_manifest.json':continue
        manifest[relative]=sha256(path)
    (ROOT/'results/validation/artifact_manifest.json').write_text(json.dumps(dict(files=manifest,
        excludes='Mutable logs/cache and this manifest; frozen v14 files have their own checked manifest.'),indent=2)+'\n')
    print(json.dumps(dict(passed=result['passed'],checked_links=len(links),active_hashed_files=len(manifest),failures=failures),indent=2))
    return 0 if result['passed'] else 1


if __name__=='__main__':raise SystemExit(main())
