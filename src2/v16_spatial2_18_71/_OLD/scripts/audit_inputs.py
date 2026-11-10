"""Confere a procedência v11 e extrai sinais novos das linhas 18–71.

A reconstrução pixel a pixel é evidência herdada da análise 20–59, cuja
procedência é validada por hashes. Esta execução não remonta panoramas nem
reestima deslocamentos. A planilha original continua representando 20–59.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from PIL import Image
from extraction_config import (ROW_START,ROW_STOP,ROW_END,ROW_COUNT,ROW_SLICE,
    REFERENCE_ROW_SLICE,FIXED_ROW,UPPER_SLICE,LOWER_SLICE)

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT.parent/'v16_spatial2'
DATA=ROOT.parent/'data_processed'
EXCLUDED={'m03q01_001','m03q80_001','m04q05_001','m05q05_002','m07q10_001'}

def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    out=ROOT/'results';cache=ROOT/'intermediate';signals=out/'signals'
    cache.mkdir(exist_ok=True);signals.mkdir(exist_ok=True)
    frozen=json.loads((out/'source_v16_sha256.json').read_text())
    sources={}
    def checked_source(relative):
        path=SOURCE/relative
        value=sha256(path)
        assert value==frozen[relative],f'Fonte modificada desde a criação da pasta: {path}'
        sources[path]=value
        return path
    source_audit=checked_source('results/input_audit.csv')
    source_summary=checked_source('results/input_audit_summary.json')
    inventory=pd.read_csv(source_audit)
    assert len(inventory)==31 and inventory.included.sum()==26
    assert set(inventory.loc[~inventory.included,'run'])==EXCLUDED
    actual={p.parent.name for p in DATA.glob('m*/__panorama_v11.png')}
    assert actual==set(inventory.run),'O conjunto de panoramas mudou.'
    rows=[]
    for previous in inventory.sort_values('run').to_dict(orient='records'):
        name=previous['run'];folder=DATA/name
        path=folder/'__panorama_v11.png';book=folder/'__dados_cropped_v11.xlsx'
        for source,expected_hash in [(path,previous['panorama_sha256']),(book,previous['xlsx_sha256'])]:
            actual_hash=sha256(source)
            assert actual_hash==expected_hash,f'Fonte alterada: {source}'
            sources[source]=actual_hash
        with Image.open(path) as image:pano=np.asarray(image.convert('L'))
        assert pano.ndim==2 and pano.shape==(previous['height'],previous['n'])
        assert pano.shape[0]>=ROW_STOP
        frame=pd.read_excel(book,sheet_name='Sheet')
        reference_workbook=frame['Voltage_2'].to_numpy(float)
        reference_image=pano[REFERENCE_ROW_SLICE].mean(axis=0)
        assert np.array_equal(frame.position_px.to_numpy(),np.arange(pano.shape[1]))
        assert np.isfinite(reference_workbook).all() and len(reference_workbook)==pano.shape[1]
        xlsx_error=float(np.max(np.abs(reference_workbook-reference_image)))
        assert xlsx_error<1e-10,f'Planilha original não corresponde às linhas 20–59: {name}'
        labels=[p.stem.lstrip('_') for p in folder.glob('*.txt') if p.stem.lstrip('_')!='dados']
        assert labels==[previous['pattern']],f'Rótulo diferente ou ambíguo: {name}'
        input_source=checked_source(f'intermediate/{name}_input.npz')
        with np.load(input_source) as old:
            reference=old['raw'].copy()
            inherited={k:old[k].copy() for k in ['raw_shifts','shifts','scores','seams']}
        assert np.array_equal(reference,reference_workbook),f'Cache antigo não reproduz a planilha: {name}'
        assert previous['reconstruction_exact'] and previous['rebuilt_width']==pano.shape[1]
        expected_shifts=np.maximum(np.r_[inherited['raw_shifts'],inherited['raw_shifts'][-1]],1)
        assert np.array_equal(inherited['shifts'],expected_shifts)
        assert np.array_equal(inherited['seams'],np.cumsum(inherited['shifts'])[:-1])
        assert len(inherited['shifts'])==previous['frames']
        # New input: grayscale arithmetic mean over exactly 54 image rows.
        raw=pano[ROW_SLICE].mean(axis=0)
        independent=pano[ROW_START:ROW_END+1].sum(axis=0,dtype=np.uint64)/ROW_COUNT
        assert np.array_equal(raw,independent) and np.isfinite(raw).all()
        np.savez_compressed(cache/f'{name}_input.npz',raw=raw,raw_20_59=reference,
            center=pano[FIXED_ROW].astype(float),upper=pano[UPPER_SLICE].mean(axis=0),
            lower=pano[LOWER_SLICE].mean(axis=0),**inherited)
        signal_path=signals/f'{name}_rows_18_71.csv'
        pd.DataFrame(dict(position_px=np.arange(len(raw)),gray_mean_rows_18_71=raw,
                          gray_mean_rows_20_59_reference=reference)).to_csv(signal_path,index=False,float_format='%.17g')
        reloaded=pd.read_csv(signal_path,float_precision='round_trip')
        assert np.array_equal(reloaded.gray_mean_rows_18_71.to_numpy(),raw)
        row=dict(previous)
        row.update(source_raw_std_20_59=previous['raw_std'],
            source_saturation_fraction_20_59=previous['saturation_fraction'],
            raw_std=float(raw.std()),saturation_fraction=float((pano[ROW_SLICE]==255).mean()),
            xlsx_max_error=xlsx_error,xlsx_rows_start=20,xlsx_rows_end_inclusive=59,
            extraction_rows_start=ROW_START,extraction_rows_end_inclusive=ROW_END,
            extraction_rows_stop_exclusive=ROW_STOP,extraction_row_count=ROW_COUNT,
            extraction_max_error=float(np.max(np.abs(raw-independent))),
            reference_cache_max_error=float(np.max(np.abs(reference-reference_image))),
            source_input_sha256=sources[input_source],source_input_file=str(input_source.relative_to(ROOT.parent)),
            reconstruction_executed_this_run=False,reconstruction_provenance='v16_spatial2/results/input_audit.csv',
            stitching_rows_start=20,stitching_rows_end_inclusive=59,
            signal_csv=str(signal_path.relative_to(ROOT)),signal_csv_sha256=sha256(signal_path),
            raw_array_sha256=hashlib.sha256(raw.astype('<f8',copy=False).tobytes()).hexdigest(),
            source_hashes_verified=True)
        rows.append(row)
        print(f'{name}: sinal 18–71 conferido; planilha 20–59 conferida; panorama original preservado.',flush=True)
    for path,value in sources.items():assert sha256(path)==value,f'Fonte alterada durante a auditoria: {path}'
    pd.DataFrame(rows).to_csv(out/'input_audit.csv',index=False)
    previous_summary=json.loads(source_summary.read_text())
    summary=dict(runs=len(rows),included=sum(r['included'] for r in rows),
        extraction_rows_inclusive=[ROW_START,ROW_END],numpy_slice=[ROW_START,ROW_STOP],row_count=ROW_COUNT,
        all_extractions_exact=all(r['extraction_max_error']==0 for r in rows),
        all_xlsx_match=all(r['xlsx_max_error']<1e-10 for r in rows),xlsx_comparison_rows_inclusive=[20,59],
        all_reconstruction_exact=all(r['reconstruction_exact'] for r in rows),
        reconstruction_executed_this_run=False,reconstruction_evidence='Auditoria anterior conferida por hashes dos mesmos panoramas.',
        literal_shift_checks=0,source_literal_shift_checks=previous_summary['literal_shift_checks'],
        shifts_and_seams='Reutilizados dos caches verificados da análise 20–59; não reestimados.',
        raw_20_59='Valores do cache anterior, idênticos à planilha original; utilizados somente como referência.',
        source_audit_sha256=sources[source_audit],source_files_verified=len(sources),source_files_unchanged=True,
        signal_array_hash_definition='SHA-256 dos bytes float64 little-endian do vetor raw na ordem das posições.',
        signal_csv_files=len(rows))
    (out/'input_audit_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':main()
