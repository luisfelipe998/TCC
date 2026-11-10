#!/usr/bin/env python3
"""Verify final scope, report provenance, units and local artifact integrity."""
from pathlib import Path
import ast
import hashlib
import json
import re
from urllib.parse import unquote
import pandas as pd

BASE=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    results=BASE/'results'
    from publish_latest_summary import check_publication
    main_reports=check_publication()
    latest_path=BASE/'all_scenarios_adversarial/results/final_verification.json'
    latest=None
    if latest_path.exists():
        latest=json.loads(latest_path.read_text())
        assert latest['status']=='passed'
    audits={}
    for name,root in [('independent_verification.json',BASE.parent),
                      ('independent_regional_verification.json',BASE),
                      ('independent_paper_verification.json',BASE)]:
        audit=json.loads((results/name).read_text())
        assert audit['status']=='passed',name
        for relative,value in audit['input_sha256'].items():
            assert sha(root/relative)==value,('Stale audit',name,relative)
        script=BASE/'scripts'/name.replace('.json','.py')
        assert sha(script)==audit['verifier_sha256'],('Changed verifier',script)
        audits[name]=dict(status='passed',inputs_checked=len(audit['input_sha256']))
    inventory=pd.read_csv(results/'scenario_inventory.csv')
    selected=inventory[inventory.included]
    assert len(inventory)==31 and len(selected)==26 and selected.run.is_unique
    assert (selected.frames==501).all()
    status=pd.read_csv(results/'all_run_status.csv').set_index('run')
    assert set(status.index)==set(selected.run) and status.index.is_unique
    assert not status.gas_liquid_material_pair.any() and not status.superficial_velocity_pair.any()
    assert int(status.passage_comparison.sum())==4
    assert int(status.regional_comparison.sum())==10
    assert int(status.outer_regional_comparison.sum())==12
    assert int(status.core_optical_propagation.sum())==20
    assert status.loc['m04q40_001','regional_reasons'].find('lag_peak_evidence')>=0
    link_count=0
    for path in BASE.rglob('*.md'):
        for target in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)',path.read_text()):
            target=target.strip('<>')
            if target.startswith(('http:','https:','mailto:','#')):continue
            actual=path.parent/unquote(target.split('#')[0])
            assert actual.exists(),('Broken link',path,target)
            link_count+=1
    scripts=list((BASE/'scripts').glob('*.py'))
    for script in scripts:ast.parse(script.read_text())
    json_count=0
    for path in results.glob('*.json'):
        json.loads(path.read_text());json_count+=1
    from run_all import STAGES
    assert all((BASE/'scripts'/stage).exists() for stage in STAGES)
    manifest={}
    for path in sorted(BASE.rglob('*')):
        rel=path.relative_to(BASE)
        if not path.is_file() or path.name in {'artifact_manifest.json','final_verification.json'}:continue
        if any('cache' in part or part in ('__pycache__','logs') for part in rel.parts):continue
        manifest[str(rel)]=dict(sha256=sha(path),bytes=path.stat().st_size)
    (results/'artifact_manifest.json').write_text(json.dumps(manifest,indent=2)+chr(10))
    summary=dict(status='passed',independent_audits=audits,selected_runs=26,
                  passage_pairs=4,inner_region_pairs=10,outer_region_pairs=12,supported_core_regions=20,
                  documented_material_phase_pairs=0,relative_links_checked=link_count,
                  original_result_scope='Original passage and fixed-region analysis; current report coverage is recorded under latest_investigation.',
                  main_reports_match_latest=main_reports,
                  latest_investigation=latest,
                  scripts_parsed=len(scripts),json_files_parsed=json_count,manifest_artifacts=len(manifest),
                  main_units='px/s',conditional_physical_scale_m_px=.001/53,
                  precision_limit='These checks validate numerical consistency and provenance, not physical phase velocities or calibrated measurement accuracy.')
    (results/'final_verification.json').write_text(json.dumps(summary,indent=2)+chr(10))
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    run()
