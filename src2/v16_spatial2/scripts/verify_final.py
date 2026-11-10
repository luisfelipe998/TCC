"""Verify coverage, equations in saved maps, provenance and local report links."""
import ast
from datetime import datetime, timezone
import hashlib
import importlib.metadata as metadata
import inspect
import json
from pathlib import Path
import re
import sys
import zipfile
from urllib.parse import unquote
import numpy as np
import pandas as pd
import pywt
import scipy
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        while chunk:=f.read(8*1024*1024):h.update(chunk)
    return h.hexdigest()

def main():
    inv=pd.read_csv(ROOT/'results/input_audit.csv')
    assert len(inv)==31 and inv.included.sum()==26
    assert inv.reconstruction_exact.all() and (inv.xlsx_max_error<1e-10).all()
    for r in inv.itertuples():
        base=ROOT.parent/'data_processed'/r.run
        assert sha(base/'__panorama_v11.png')==r.panorama_sha256
        assert sha(base/'__dados_cropped_v11.xlsx')==r.xlsx_sha256
    p=pd.read_csv(ROOT/'results/pair_methods.csv')
    assert len(p)==261 and p.groupby('method').size().eq(87).all()
    chosen=p[p.method=='v16_biased']
    assert (chosen.max_wco<=1+1e-8).all()
    names=set(inv[inv.included].run)
    assert set(chosen.run_x)|set(chosen.run_y) <= names
    assert not chosen.duplicated(['run_x','run_y']).any()
    run=pd.read_csv(ROOT/'results/per_run_methods.csv')
    assert len(run)==93 and (run.self_wco_zero_error<1e-12).all()
    assert (run[run.method=='v16_biased'].self_wco_max<=1+1e-8).all()
    error=0.
    migration=json.loads((ROOT/'results/migracao_cache_sem_mascaras.json').read_text())
    preserved={r['file']:r['preserved_array_sha256'] for r in migration['records']}
    for row in chosen.itertuples():
        name=f'{row.run_x}__{row.run_y}'
        assert (ROOT/'figures/pairs'/f'{name}.png').exists()
        with np.load(ROOT/'intermediate'/f'{name}_maps.npz') as z:
            zero=np.flatnonzero(z['lags']==0)
            assert len(zero)==1
            den=abs(z['wxx'][:,zero[0]])*abs(z['wyy'][:,zero[0]])
            expected=abs(z['wxy'])**2/den[:,None]
            finite=np.isfinite(z['wco'])
            error=max(error,float(np.max(abs(expected[finite]-z['wco'][finite]))))
            assert not {'supported','coi','mean_x_coi','mean_y_coi'}.intersection(z.files)
            assert all(hashlib.sha256(z[k].tobytes()).hexdigest()==value
                       for k,value in preserved[f'{name}_maps.npz'].items())
    assert error<1e-12
    for path in (ROOT/'scripts').glob('*.py'):ast.parse(path.read_text())
    pngs=list((ROOT/'figures').rglob('*.png'))
    for path in pngs:
        with Image.open(path) as image:image.verify()
    for filename in ['numerical_validation.json','legacy_literal_validation.json']:
        assert json.loads((ROOT/'results'/filename).read_text())['status']=='PASS'
    null=pd.read_csv(ROOT/'results/null_diagnostic.csv')
    assert len(null)==5 and (null.surrogates==199).all()
    assert len(pd.read_csv(ROOT/'results/null_replicates.csv'))==995
    screen=pd.read_csv(ROOT/'results/preprocessing_screen_pairs.csv')
    assert len(screen)==10*325 and screen.groupby('variant').same_pattern.sum().eq(87).all()
    assert np.isfinite(screen.weighted_scale_peak_power05).all()
    confirm=pd.read_csv(ROOT/'results/preprocessing_confirm_pairs.csv')
    assert len(confirm)==len(confirm.variant.unique())*325
    assert json.loads((ROOT/'results/preprocessing_confirm_metadata.json').read_text())['scale_count']==200
    assert (screen.max_wco_all<=1+1e-8).all() and (confirm.max_wco_all<=1+1e-8).all()
    # The primary pair summaries and the preprocessing baseline must agree on
    # full-resolution maxima: the cached image maps cannot establish this.
    legacy=pd.read_csv(ROOT/'results/preprocessing_v5_pairs.csv')
    for frame,method in [(confirm[confirm.variant=='baseline'],'v16_biased'),
                         (legacy,'v5_recomputed')]:
        match=p[p.method==method].merge(frame,on=['run_x','run_y'],validate='one_to_one')
        assert len(match)==87
        assert np.allclose(match.max_wco,match.max_wco_all,rtol=1e-10,atol=1e-12)
    observed=null.merge(confirm[confirm.variant=='baseline'],on=['run_x','run_y'],validate='one_to_one')
    assert len(observed)==5
    assert np.allclose(observed.observed_max,observed.max_wco_all,rtol=1e-10,atol=1e-12)
    replicates=pd.read_csv(ROOT/'results/null_replicates.csv')
    for row in null.itertuples():
        values=replicates[(replicates.run_x==row.run_x)&(replicates.run_y==row.run_y)].maximum_wco
        assert len(values)==199
        assert np.isclose(np.quantile(values,.95),row.null_q95,rtol=1e-12)
        assert np.isclose((1+np.count_nonzero(values>=row.observed_max))/200,row.conditional_monte_carlo_p)
    for name in ['preprocessing_screen_metadata.json','preprocessing_confirm_metadata.json',
                 'preprocessing_v5_metadata.json','null_metadata.json']:
        assert 'sem COI' in json.loads((ROOT/'results'/name).read_text())['domain']
    for csv_path in (ROOT/'results').glob('*.csv'):
        columns=pd.read_csv(csv_path,nrows=0).columns
        assert not any('coi' in c or 'supported' in c for c in columns),csv_path
    # Confirm that the mathematical functions and every retained full-domain
    # value from the earlier run survived this presentation/domain revision.
    origin=json.loads((ROOT/'results/revisao_sem_mascaras_origem.json').read_text())
    tree=ast.parse((ROOT/'scripts/wavelet_core.py').read_text())
    current_functions={f.name:ast.dump(f,include_attributes=False) for f in tree.body
                       if isinstance(f,ast.FunctionDef)}
    assert all(current_functions[name]==value for name,value in origin['core_equation_ast'].items())
    changed=[name for name,old in origin['previous_figure_sha256'].items() if sha(ROOT/name)!=old]
    assert len(changed)==len(origin['previous_figure_sha256'])==136
    with zipfile.ZipFile(ROOT/'results/historico_com_mascaras.zip') as archive:
        for name in ['profiles.csv','per_run_methods.csv','pair_methods.csv']:
            old=pd.read_csv(archive.open('results/'+name))
            new=pd.read_csv(ROOT/'results'/name)
            pd.testing.assert_frame_equal(new,old[new.columns],check_exact=False,rtol=1e-12,atol=1e-12)
    assert (ROOT/'EXECUTIVE_SUMMARY.md').exists() and (ROOT/'PREPROCESSING.md').exists()
    broken=[]
    for path in ROOT.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',path.read_text()):
            if '://' in link or link.startswith('#'):continue
            target=unquote(link.split('#')[0])
            resolved=(path.parent/target).resolve()
            # This report is produced at the end of the current verification.
            planned=(ROOT/'results/final_verification.json').resolve()
            if not resolved.exists() and resolved!=planned:broken.append(f'{path.name}: {link}')
    assert not broken,broken
    versions={pkg:metadata.version(pkg) for pkg in ['numpy','scipy','pandas','PyWavelets','matplotlib','Pillow','openpyxl']}
    source_paths=[ROOT.parent/'v5_spatial/wavelet_v5.py',ROOT.parent/'pipeline_v11.py',ROOT.parent/'pipeline_v13.py',ROOT.parent/'pipeline_v9.py']
    source_paths += [ROOT.parents[1]/'assets'/x for x in ['schuck-bardo-wavelet-coherence.pdf','deoliveira2021-flow-patterns-pressure-drop.pdf','torrence1998-wavelet-analysis-guide.pdf','grinsted2004-wavelet-coherence-geophysical.pdf']]
    source_paths += [ROOT/'literature/maraun2004.pdf',Path(inspect.getfile(pywt.cwt))]
    manifest=dict(created_utc=datetime.now(timezone.utc).isoformat(),python=sys.version,
        distribution_versions=versions,module_versions=dict(numpy=np.__version__,scipy=scipy.__version__,pywt=pywt.__version__),
        pywavelets_version_discrepancy='Distribution metadata 1.9.0; loaded module reports 1.8.0. Core CWT source hashed below.',
        source_sha256={str(path):sha(path) for path in source_paths},
        script_sha256={str(path.relative_to(ROOT)):sha(path) for path in sorted((ROOT/'scripts').glob('*.py'))},
        parameters=dict(omega0=6,spatial_sampling_px=1,wavelength_limits_px=[50,3000],main_scales=200,
                        unbiased='sensitivity only',primary_correlation='biased',denominator_relative_floor=1e-10,
                        padding_sigma=8,domain='Plano completo, sem COI e sem máscara conservadora',
                        lag_mean='Média do módulo em todos os lags'),
        output_sha256={str(path.relative_to(ROOT)):sha(path) for path in sorted(ROOT.rglob('*'))
                       if path.is_file() and path.suffix in ['.md','.csv','.png','.npz','.html','.txt','.json','.zip','.log']
                       and '__pycache__' not in str(path)
                       and path.name not in ['artifact_manifest.json','final_verification.json']})
    (ROOT/'results/artifact_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    result=dict(status='PASS',input_hashes_verified=62,panoramas=31,eligible_runs=26,pairs=87,
                primary_maps_equation_max_error=error,png_files_verified=len(pngs),null_replicates=995,
                preprocessing_screen_pairs_per_variant=325,preprocessing_variants=10,
                preprocessing_confirm_variants=sorted(confirm.variant.unique()),broken_local_links=broken,
                regenerated_figures=len(changed),equation_functions_unchanged=True,
                retained_full_domain_tables_unchanged=True,cached_wavelet_arrays_unchanged=87,
                masked_columns_in_current_csv=0)
    (ROOT/'results/final_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
