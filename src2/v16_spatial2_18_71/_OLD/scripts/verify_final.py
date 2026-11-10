"""Verifica extração 18–71, equações, comparações e preservação da fonte.

Não executa novamente as transformadas: confere as entradas diretamente nos
PNGs e a Eq. 32.17 nos coeficientes salvos, além de cruzar os resumos integrais.
"""
import ast
from datetime import datetime, timezone
import hashlib
import importlib.metadata as metadata
import inspect
from itertools import combinations
import json
from pathlib import Path
import re
import sys
from urllib.parse import unquote
import numpy as np
import pandas as pd
import pywt
import scipy
from PIL import Image
from extraction_config import (CONFIG, ROW_START, ROW_END, ROW_STOP, ROW_COUNT,
    ROW_SLICE, REFERENCE_ROW_SLICE, FIXED_ROW, UPPER_SLICE, LOWER_SLICE)

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent/'v16_spatial2'
EXCLUDED = {'m03q01_001','m03q80_001','m04q05_001','m05q05_002','m07q10_001'}
METRICS = ['weighted_scale_peak_power05','weighted_scale_peak_power01',
           'weighted_scale_peak_power10','weighted_zero_lag_wco','max_wco_all']


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        while chunk := stream.read(8*1024*1024):
            h.update(chunk)
    return h.hexdigest()


def read_json(name):
    return json.loads((ROOT/'results'/name).read_text())


def verify_source():
    frozen = read_json('source_v16_sha256.json')
    current = {str(p.relative_to(SOURCE)) for p in SOURCE.rglob('*')
               if p.is_file() and '__pycache__' not in p.parts}
    assert current == set(frozen), 'Arquivos adicionados ou removidos na análise original.'
    for name, value in frozen.items():
        assert sha(SOURCE/name) == value, f'Fonte original alterada: {name}'
    assert sha(ROOT/'scripts/wavelet_core.py') == frozen['scripts/wavelet_core.py']
    return frozen


def verify_inputs(inv, frozen):
    assert len(inv) == 31 and inv.included.sum() == 26
    assert not inv.run.duplicated().any()
    assert set(inv.loc[~inv.included,'run']) == EXCLUDED
    assert inv.reconstruction_exact.all() and (inv.xlsx_max_error < 1e-10).all()
    assert not inv.reconstruction_executed_this_run.any()
    assert inv.source_hashes_verified.all()
    for field, value in [('extraction_rows_start',ROW_START),
                         ('extraction_rows_end_inclusive',ROW_END),
                         ('extraction_rows_stop_exclusive',ROW_STOP),
                         ('extraction_row_count',ROW_COUNT),
                         ('xlsx_rows_start',20),('xlsx_rows_end_inclusive',59),
                         ('stitching_rows_start',20),('stitching_rows_end_inclusive',59)]:
        assert inv[field].eq(value).all(), field
    largest_error = 0.0
    for row in inv.itertuples():
        base = ROOT.parent/'data_processed'/row.run
        assert sha(base/'__panorama_v11.png') == row.panorama_sha256
        assert sha(base/'__dados_cropped_v11.xlsx') == row.xlsx_sha256
        with Image.open(base/'__panorama_v11.png') as image:
            panorama = np.asarray(image.convert('L'))
        assert panorama.shape == (row.height,row.n) and row.height >= ROW_STOP
        book = pd.read_excel(base/'__dados_cropped_v11.xlsx',sheet_name='Sheet')
        assert np.array_equal(book.position_px.to_numpy(),np.arange(row.n))
        assert np.allclose(book.Voltage_2,panorama[REFERENCE_ROW_SLICE].mean(axis=0),
                           rtol=0,atol=1e-10)
        source_name = f'intermediate/{row.run}_input.npz'
        assert row.source_input_sha256 == frozen[source_name]
        assert row.source_input_file == str((SOURCE/source_name).relative_to(ROOT.parent))
        with np.load(SOURCE/source_name) as previous, np.load(ROOT/source_name) as current:
            raw = current['raw']
            # Independent integer sum avoids relying only on the extraction helper.
            expected = panorama[ROW_START:ROW_END+1].sum(axis=0,dtype=np.uint64)/ROW_COUNT
            assert np.array_equal(raw,expected) and np.isfinite(raw).all()
            largest_error = max(largest_error,float(np.max(abs(raw-expected))))
            assert np.array_equal(current['raw_20_59'],previous['raw'])
            assert np.array_equal(current['raw_20_59'],book.Voltage_2.to_numpy())
            assert np.array_equal(current['center'],panorama[FIXED_ROW])
            assert np.array_equal(current['upper'],panorama[UPPER_SLICE].mean(axis=0))
            assert np.array_equal(current['lower'],panorama[LOWER_SLICE].mean(axis=0))
            for field in ['raw_shifts','shifts','scores','seams']:
                assert np.array_equal(current[field],previous[field],equal_nan=True)
            assert np.array_equal(current['seams'],np.cumsum(current['shifts'])[:-1])
            assert hashlib.sha256(raw.astype('<f8',copy=False).tobytes()).hexdigest() == row.raw_array_sha256
            signal_path = ROOT/row.signal_csv
            assert sha(signal_path) == row.signal_csv_sha256
            signals = pd.read_csv(signal_path,float_precision='round_trip')
            assert np.array_equal(signals.position_px,np.arange(row.n))
            assert np.array_equal(signals.gray_mean_rows_18_71,raw)
            assert np.array_equal(signals.gray_mean_rows_20_59_reference,current['raw_20_59'])
    summary = read_json('input_audit_summary.json')
    assert summary['all_extractions_exact'] and summary['all_xlsx_match']
    assert summary['reconstruction_executed_this_run'] is False
    assert summary['literal_shift_checks'] == 0
    assert summary['extraction_rows_inclusive'] == [18,71]
    return largest_error


def verify_all_pairs(frame, inv, variants):
    selected = inv[inv.included].set_index('run')
    expected = set(combinations(sorted(selected.index),2))
    assert len(expected) == 325 and set(frame.variant.unique()) == set(variants)
    assert len(frame) == len(variants)*325
    for variant, group in frame.groupby('variant'):
        assert not group.duplicated(['run_x','run_y']).any(), variant
        assert set(zip(group.run_x,group.run_y)) == expected, variant
        assert group.same_pattern.sum() == 87 and (~group.same_pattern).sum() == 238
        assert np.isfinite(group.weighted_scale_peak_power05).all()
        for row in group.itertuples():
            assert row.pattern_x == selected.loc[row.run_x,'pattern']
            assert row.pattern_y == selected.loc[row.run_y,'pattern']
            assert row.same_pattern == (row.pattern_x == row.pattern_y)
            assert row.n == min(selected.loc[row.run_x,'n'],selected.loc[row.run_y,'n'])


def verify_maps(inv, pairs):
    methods = {'v5_recomputed','v16_biased','schuck_unbiased'}
    assert len(pairs) == 261 and set(pairs.method) == methods
    assert pairs.groupby('method').size().eq(87).all()
    selected = inv[inv.included]
    expected = {pair for _, group in selected.groupby('pattern')
                for pair in combinations(sorted(group.run),2)}
    for _, group in pairs.groupby('method'):
        assert set(zip(group.run_x,group.run_y)) == expected
        assert not group.duplicated(['run_x','run_y']).any()
    chosen = pairs[pairs.method == 'v16_biased']
    assert (chosen.max_wco <= 1+1e-8).all()
    runs = pd.read_csv(ROOT/'results/per_run_methods.csv')
    assert len(runs) == 93 and runs.groupby('method').size().eq(31).all()
    assert set(runs.method) == {'v5','v16_biased','schuck_unbiased'}
    assert not runs.duplicated(['run','method']).any()
    assert set(runs.run) == set(inv.run)
    assert (runs.self_wco_zero_error < 1e-12).all()
    assert (runs[runs.method == 'v16_biased'].self_wco_max <= 1+1e-8).all()
    profiles = pd.read_csv(ROOT/'results/profiles.csv')
    assert len(profiles) == 93*200
    assert profiles.groupby(['run','method']).size().eq(200).all()
    assert not profiles.duplicated(['run','method','wavelength_px']).any()
    assert np.isfinite(profiles[['zero_lag_wps','lag_mean_all']]).all().all()
    error = 0.0
    expected_maps = {f'{a}__{b}_maps.npz' for a,b in expected}
    assert {p.name for p in (ROOT/'intermediate').glob('*_maps.npz')} == expected_maps
    for row in chosen.itertuples():
        name = f'{row.run_x}__{row.run_y}'
        assert (ROOT/'figures/pairs'/f'{name}.png').exists()
        with np.load(ROOT/'intermediate'/f'{name}_maps.npz') as z:
            assert set(z.files) == {'wavelengths','wxx','wyy','wxy','wco','lags','mean_x_all','mean_y_all'}
            assert np.array_equal(z['wavelengths'],np.geomspace(50,3000,200))
            assert np.all(np.diff(z['lags']) > 0) and np.max(abs(z['lags'])) <= row.n-1
            zero = np.flatnonzero(z['lags'] == 0)
            assert len(zero) == 1
            for field in ['wxx','wyy','wxy','wco']:
                assert z[field].shape == (200,len(z['lags']))
            zx,zy = abs(z['wxx'][:,zero[0]]),abs(z['wyy'][:,zero[0]])
            valid = (zx > 1e-10*zx.max()) & (zy > 1e-10*zy.max())
            expected_wco = np.divide(abs(z['wxy'])**2,(zx*zy)[:,None],
                out=np.full(z['wco'].shape,np.nan),where=valid[:,None])
            finite = np.isfinite(z['wco'])
            assert finite.any() and np.array_equal(finite,np.isfinite(expected_wco))
            error = max(error,float(np.max(abs(expected_wco[finite]-z['wco'][finite]))))
            assert np.nanmax(z['wco']) <= row.max_wco+1e-12
            assert np.nanmax(z['wco']) <= 1+1e-8
            assert all(z[f'mean_{key}_all'].shape == (200,) for key in ['x','y'])
    assert error < 1e-12
    return error


def verify_band_comparison(confirm, legacy):
    combined = pd.read_csv(ROOT/'results/band_comparison_pairs.csv')
    assert len(combined) == 2*4*325
    assert set(combined.band) == {'18_71','20_59'}
    key = ['variant','run_x','run_y']
    for band,root in [('18_71',ROOT),('20_59',SOURCE)]:
        source = (pd.concat([confirm,legacy],ignore_index=True) if band == '18_71'
                  else pd.concat([pd.read_csv(root/'results/preprocessing_confirm_pairs.csv'),
                                  pd.read_csv(root/'results/preprocessing_v5_pairs.csv')],ignore_index=True))
        merged = combined[combined.band == band].merge(source,on=key,validate='one_to_one',suffixes=('_saved',''))
        assert len(merged) == 1300
        for metric in METRICS:
            original_field = 'fixed20_59_'+metric if band == '18_71' and metric != 'max_wco_all' else metric
            assert np.allclose(merged[metric+'_saved'],merged[original_field],
                               rtol=1e-12,atol=1e-12,equal_nan=True)
    old = combined[(combined.band == '20_59') & (combined.variant == 'central_band20')]
    new = combined[(combined.band == '18_71') & (combined.variant == 'central_band20')]
    control = old.merge(new,on=['run_x','run_y'],validate='one_to_one',suffixes=('_old','_new'))
    control_error = max(float(np.max(abs(control[m+'_old']-control[m+'_new']))) for m in METRICS)
    assert control_error < 1e-10
    outcome = read_json('band_comparison.json')
    assert outcome['status'] == 'PASS'
    assert outcome['fixed_control_30_49_score_max_error'] < 1e-10
    changes = pd.read_csv(ROOT/'results/band_change_summary.csv')
    assert len(changes) == 4*len(METRICS)
    assert (changes.same_pairs_improved+changes.same_pairs_worsened+changes.same_pairs_unchanged).eq(87).all()
    for row in changes.itertuples():
        old = combined[(combined.band == '20_59') & (combined.variant == row.variant)]
        new = combined[(combined.band == '18_71') & (combined.variant == row.variant)]
        match = old.merge(new,on=['run_x','run_y'],validate='one_to_one',suffixes=('_old','_new'))
        delta = match[row.metric+'_new']-match[row.metric+'_old']
        same = match.same_pattern_new
        assert np.isclose(delta[same].mean(),row.delta_pooled_same,rtol=1e-12,atol=1e-12)
        assert int((delta[same] > 1e-12).sum()) == row.same_pairs_improved
        assert int((delta[same] < -1e-12).sum()) == row.same_pairs_worsened
    assert len(pd.read_csv(ROOT/'results/band_comparison_summary.csv')) == 2*4*len(METRICS)
    assert len(pd.read_csv(ROOT/'results/band_change_pairs.csv')) == 4*325
    return control_error


def main():
    frozen = verify_source()
    inv = pd.read_csv(ROOT/'results/input_audit.csv')
    extraction_error = verify_inputs(inv,frozen)
    pairs = pd.read_csv(ROOT/'results/pair_methods.csv')
    equation_error = verify_maps(inv,pairs)
    for path in (ROOT/'scripts').glob('*.py'):
        ast.parse(path.read_text())
    expected_figures = {name for name in frozen if name.startswith('figures/') and name.endswith('.png')}
    actual_figures = {str(path.relative_to(ROOT)) for path in (ROOT/'figures').rglob('*.png')}
    assert len(expected_figures) == 136 and expected_figures <= actual_figures
    pngs = list(ROOT.rglob('*.png'))
    for path in pngs:
        with Image.open(path) as image:
            image.verify()
    for filename in ['numerical_validation.json','legacy_literal_validation.json']:
        assert read_json(filename)['status'] == 'PASS'
    screen = pd.read_csv(ROOT/'results/preprocessing_screen_pairs.csv')
    verify_all_pairs(screen,inv,['baseline','demean_only','center_row','central_band20',
        'core_row_median','core_trimmed_mean20','spatial_gaussian2','spatial_median5','rank_uniform','winsor_1pct'])
    confirm = pd.read_csv(ROOT/'results/preprocessing_confirm_pairs.csv')
    verify_all_pairs(confirm,inv,['baseline','central_band20','core_row_median'])
    legacy = pd.read_csv(ROOT/'results/preprocessing_v5_pairs.csv')
    verify_all_pairs(legacy,inv,['v5'])
    assert (screen.max_wco_all <= 1+1e-8).all() and (confirm.max_wco_all <= 1+1e-8).all()
    assert read_json('preprocessing_screen_metadata.json')['scale_count'] == 64
    assert read_json('preprocessing_confirm_metadata.json')['scale_count'] == 200
    for frame,method in [(confirm[confirm.variant == 'baseline'],'v16_biased'),(legacy,'v5_recomputed')]:
        match = pairs[pairs.method == method].merge(frame,on=['run_x','run_y'],validate='one_to_one')
        assert len(match) == 87
        assert np.allclose(match.max_wco,match.max_wco_all,rtol=1e-10,atol=1e-12)
    null = pd.read_csv(ROOT/'results/null_diagnostic.csv')
    replicates = pd.read_csv(ROOT/'results/null_replicates.csv')
    assert len(null) == 5 and (null.surrogates == 199).all() and len(replicates) == 995
    assert not replicates.duplicated(['run_x','run_y','replicate']).any()
    assert np.isfinite(replicates.maximum_wco).all() and (replicates.maximum_wco <= 1+1e-8).all()
    observed = null.merge(confirm[confirm.variant == 'baseline'],on=['run_x','run_y'],validate='one_to_one')
    assert len(observed) == 5
    assert np.allclose(observed.observed_max,observed.max_wco_all,rtol=1e-10,atol=1e-12)
    for row in null.itertuples():
        samples = replicates[(replicates.run_x == row.run_x) & (replicates.run_y == row.run_y)]
        assert len(samples) == 199 and set(samples.replicate) == set(range(199))
        values = samples.maximum_wco
        count = int(np.count_nonzero(values >= row.observed_max))
        assert np.isclose(np.quantile(values,.95),row.null_q95,rtol=1e-12)
        assert np.isclose(np.median(values),row.null_median,rtol=1e-12)
        assert count == row.exceedances
        assert np.isclose((1+count)/200,row.conditional_monte_carlo_p,rtol=1e-12)
    for name in ['preprocessing_screen_metadata.json','preprocessing_confirm_metadata.json',
                 'preprocessing_v5_metadata.json','null_metadata.json']:
        info = read_json(name)
        assert 'sem COI' in info['domain']
        if 'extraction_rows_inclusive' in info:
            assert info['extraction_rows_inclusive'] == [18,71]
    for path in (ROOT/'results').glob('*.csv'):
        assert not any('coi' in c.lower() or 'supported' in c.lower()
                       for c in pd.read_csv(path,nrows=0).columns), path
    control_error = verify_band_comparison(confirm,legacy)
    expected_docs = ['README.md','REPORT.md','EQUATION_AUDIT.md','EXECUTIVE_SUMMARY.md',
                     'PREPROCESSING.md','RESULTS.md']
    assert all((ROOT/name).exists() for name in expected_docs)
    broken = []
    for path in ROOT.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',path.read_text()):
            if '://' in link or link.startswith('#'):
                continue
            target = unquote(link.split('#')[0])
            resolved = (path.parent/target).resolve()
            planned = {(ROOT/'results'/name).resolve()
                       for name in ['final_verification.json','artifact_manifest.json']}
            if not resolved.exists() and resolved not in planned:
                broken.append(f'{path.name}: {link}')
    assert not broken, broken
    versions = {pkg:metadata.version(pkg) for pkg in ['numpy','scipy','pandas','PyWavelets','matplotlib','Pillow','openpyxl']}
    source_paths = [ROOT.parent/'v5_spatial/wavelet_v5.py',ROOT.parent/'pipeline_v11.py',
                    ROOT.parent/'pipeline_v13.py',ROOT.parent/'pipeline_v9.py']
    source_paths += [ROOT.parents[1]/'assets'/name for name in [
        'schuck-bardo-wavelet-coherence.pdf','deoliveira2021-flow-patterns-pressure-drop.pdf',
        'torrence1998-wavelet-analysis-guide.pdf','grinsted2004-wavelet-coherence-geophysical.pdf']]
    source_paths += [ROOT/'literature/maraun2004.pdf',Path(inspect.getfile(pywt.cwt))]
    manifest = dict(created_utc=datetime.now(timezone.utc).isoformat(),python=sys.version,
        distribution_versions=versions,
        module_versions=dict(numpy=np.__version__,scipy=scipy.__version__,pywt=pywt.__version__),
        pywavelets_version_discrepancy=versions['PyWavelets'] != pywt.__version__,
        predecessor_snapshot_sha256=sha(ROOT/'results/source_v16_sha256.json'),
        predecessor_files_verified=len(frozen),predecessor_files_unchanged=True,
        equation_core_identical_to_predecessor=True,extraction_config=CONFIG,
        source_sha256={str(path):sha(path) for path in source_paths},
        script_sha256={str(path.relative_to(ROOT)):sha(path) for path in sorted((ROOT/'scripts').glob('*.py'))},
        parameters=dict(omega0=6,spatial_sampling_px=1,wavelength_limits_px=[50,3000],main_scales=200,
            unbiased='sensitivity only',primary_correlation='biased',denominator_relative_floor=1e-10,
            padding_sigma=8,domain='Plano completo, sem COI e sem máscara conservadora',
            lag_mean='Média do módulo em todos os lags',
            v5_comparator='Implementação v5 aplicada à extração 18–71; referência histórica 20–59 na comparação entre faixas',
            cross_band_weights='Alternative WPS corrigida de referência 20–59, com pesos e pisos fixos'),
        output_sha256={str(path.relative_to(ROOT)):sha(path) for path in sorted(ROOT.rglob('*'))
            if path.is_file() and path.suffix in ['.md','.csv','.png','.npz','.html','.txt','.json','.log']
            and '__pycache__' not in path.parts and path.name not in ['artifact_manifest.json','final_verification.json']})
    (ROOT/'results/artifact_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    result = dict(status='PASS',input_image_and_workbook_hashes_verified=62,panoramas=31,eligible_runs=26,
        extraction_rows_inclusive=[18,71],extraction_row_count=54,extraction_max_error=extraction_error,
        signal_csv_files_verified=31,reconstruction_executed_this_run=False,
        predecessor_files_verified=len(frozen),predecessor_files_unchanged=True,
        equation_core_identical_to_predecessor=True,primary_maps=87,primary_maps_equation_max_error=equation_error,
        main_figures_verified=len(actual_figures),all_png_files_verified=len(pngs),
        null_replicates=995,preprocessing_screen_pairs_per_variant=325,preprocessing_variants=10,
        preprocessing_confirm_variants=sorted(confirm.variant.unique()),
        cross_band_comparison='PASS',fixed_control_30_49_all_scores_max_error=control_error,
        broken_local_links=broken,masked_columns_in_current_csv=0,
        validation_scope='Entradas 18–71 verificadas nos PNGs; Eq. 32.17 nos mapas salvos; resumos integrais cruzados. Transformadas não reexecutadas por este verificador.')
    (ROOT/'results/final_verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__ == '__main__':
    main()
