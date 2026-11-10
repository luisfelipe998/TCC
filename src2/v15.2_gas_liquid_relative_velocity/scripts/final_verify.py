"""Verify scope, quantitative claims, audit manifests, publication links and outputs."""
from pathlib import Path
import hashlib,json,platform,re,sys,os
from urllib.parse import unquote
import numpy as np
import pandas as pd
BASE=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(BASE/'intermediate/mplcache'))
import scipy,matplotlib,PIL

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run():
    inv=pd.read_csv(BASE.parent/'v15.1_gas_liquid_velocity/results/scenario_inventory.csv');inv=inv[inv.included]
    final=pd.read_csv(BASE/'results/final_per_run_status.csv');metrics=json.loads((BASE/'results/report_metrics.json').read_text())
    assert len(final)==26 and final.run.is_unique and set(final.run)==set(inv.run)
    assert not final.pattern.str.contains('smooth|plug-slug|slug-churn').any() and 'm05q05_002' not in set(final.run)
    assert (inv.frames==501).all() and int(inv.frames.sum())==13026
    for col in ('identified_gas_material_velocity_px_s','identified_liquid_material_velocity_px_s','identified_material_velocity_ratio'):
        assert final[col].isna().all()
    assert not final.material_phase_pair.any()
    assert metrics['small_stable']==int(final.small_threshold_stable.sum())==5
    assert metrics['meniscus_stable']==int(final.meniscus_threshold_stable.sum())==3
    assert metrics['common_object_pairs']==int(final.threshold_stable.sum())==2
    pair=pd.read_csv(BASE/'results/paired_object_consensus.csv').set_index('run')
    expected={'m05q05_001':(.934516,120),'m07q05_001':(.963872,133)}
    assert set(pair.index[pair.threshold_stable])==set(expected)
    for name,(ratio,intervals) in expected.items():
        assert abs(pair.loc[name,'ratio']-ratio)<.0000006 and pair.loc[name,'common_intervals']==intervals
    inter=pd.read_csv(BASE/'results/inter_recording_object_ratios.csv')
    assert len(inter)==10 and inter.same_pattern.sum()==3 and not inter.phase_liquid_ratio.any()
    residual=pd.read_csv(BASE/'results/residual_counts.csv')
    assert residual.supported_pairs.tolist()==[1,2,2,0,4,4]
    stable=pd.read_csv(BASE/'results/residual_all_options_stability.csv')
    assert not stable.jointly_supported_four_background_options.any()
    conditional=json.loads((BASE/'results/conditional_piston_validation.json').read_text())
    assert conditional['conditional_tracks']==0 and conditional['conditional_runs']==0
    source_checks=[]
    for name,script in [('object_validation.json','object_motion.py'),('meniscus_validation.json','meniscus_motion.py'),
                        ('robust_meniscus_metadata.json','robust_meniscus.py'),('physical_bounds_controls.json','physical_bounds.py'),
                        ('conditional_piston_validation.json','conditional_bubble_bounds.py')]:
        data=json.loads((BASE/'results'/name).read_text())
        assert data.get('script_sha256',data.get('source_sha256'))==sha(BASE/'scripts'/script),(name,'stale source')
        source_checks.append(name)
    manifest_checks=0
    audits={}
    for name in ('object_independent_verification.json','comparison_independent_verification.json'):
        data=json.loads((BASE/'results'/name).read_text());assert data['status'].startswith('passed')
        for rel,digest in data['manifest'].items():
            path=BASE/rel
            assert path.is_file() and sha(path)==digest,(name,rel,'stale audited input')
            manifest_checks+=1
        audits[name]=dict(status=data['status'],counts=data['counts'])
    for name in ('residual_audit.json','residual_real_numerical_check.json'):
        data=json.loads((BASE/'results'/name).read_text())
        audits[name]=data
    broken=[];links=0
    for doc in [BASE/'README.md',*(BASE/'findings').glob('*.md')]:
        for target in re.findall(r'\]\(([^)]+)\)',doc.read_text()):
            if '://' in target or target.startswith('#'):continue
            dest=(doc.parent/unquote(target.split('#')[0])).resolve()
            if dest==BASE/'results/final_verification.json':continue
            links+=1
            if not dest.exists():broken.append((str(doc),target))
    assert not broken,broken
    for summary in ('EXECUTIVE_SUMMARY.md','EXECUTIVE_SUMMARY_PT_BR.md'):
        text=(BASE/'findings'/summary).read_text()
        assert '**0/26**' in text and '0.9345' in text.replace(',','.') and len(re.findall(r'^\| m\d\dq\d\d_001 \|',text,re.M))==26
    figures=['coverage_all_scenarios_PT_BR','paired_objects_PT_BR','inter_recording_ratios_PT_BR','object_size_speed_PT_BR','moving_object_examples_PT_BR']
    for name in figures:
        for ext in ('png','pdf','svg'):assert (BASE/'figures'/f'{name}.{ext}').stat().st_size>1000
    artifacts={}
    for folder in ('scripts','results','findings','figures'):
        for p in sorted((BASE/folder).glob('*')):
            if p.is_file() and p.name not in ('final_verification.json','artifact_manifest.json'):
                artifacts[str(p.relative_to(BASE))]=sha(p)
    artifacts['README.md']=sha(BASE/'README.md')
    (BASE/'results/artifact_manifest.json').write_text(json.dumps(artifacts,indent=2)+'\n')
    result=dict(status='passed',metrics=metrics,publication_links_checked=links,audited_manifest_entries=manifest_checks,
        source_metadata_checks=source_checks,artifact_files=len(artifacts),audits=audits,
        versions=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,pandas=pd.__version__,
                      matplotlib=matplotlib.__version__,pillow=PIL.__version__,platform=platform.platform()),
        limitations='Numerical, provenance and reporting verification only; no experimental material-phase velocity ground truth.')
    (BASE/'results/final_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('audits','versions')},indent=2))

if __name__=='__main__':run()
