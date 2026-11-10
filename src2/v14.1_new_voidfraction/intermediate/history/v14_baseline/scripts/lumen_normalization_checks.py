#!/usr/bin/env python3
"""Verify the user-confirmed D53 normalization independently of morphology."""
import csv
from datetime import datetime, timezone
import hashlib
import json

import numpy as np
from PIL import Image

from v14_paths import ROOT, DATA, artifact_path, write_markdown
from v14_selection import selected_runs, PATTERN_ORDER
from v14_geometry import (LUMEN_ROWS, LUMEN_SLICE, LUMEN_HEIGHT_PX,
                          DIAMETER_PX, FIXED_PROBE_ROW, BOUNDARY_PROVENANCE)

BASELINE = ROOT / 'intermediate/history/before_confirmed_lumen_D53'
OLD_DIAMETER = 60.0
AREA_FIELDS = ('conditional_circle_alpha_spatial_mean', 'conditional_circle_alpha_spatial_median',
               'threshold_060_080_min', 'threshold_060_080_max',
               'reference_variants_min', 'reference_variants_max',
               'five_spatial_block_means_min', 'five_spatial_block_means_max')
UNSCALED_VOID_FIELDS = ('run','label','panorama_width_px','detected_center_gas_column_fraction',
                       'center_gas_silhouette_duty','radial_crop_touch_column_fraction',
                       'core_saturated_fraction','validity')
PHASE_FIELDS = ('run','pattern','panorama_width','candidate_mask_rows','candidate_mask_n_rows',
                'candidate_pixel_flip_count','candidate_pixel_flip_comparisons','candidate_pixel_flip_rate',
                'center_candidate_flip_count','center_candidate_flip_rate','center_candidate_flips_per1000px',
                'center_candidate_gas_entries','center_candidate_gas_exits','center_candidate_gas_duty',
                'paired_shadow_detection_fraction','eligible')
REMOVED_FIELDS = {'conditional_missing_width_alpha_lower','conditional_missing_width_alpha_upper',
                  'conditional_alpha_D56','conditional_alpha_D64'}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table(path):
    with path.open(newline='') as stream:
        rows=list(csv.DictReader(stream))
    result={r['run']:r for r in rows}
    assert len(result)==len(rows), f'Duplicate run in {path}'
    return result


def close(actual, expected, label):
    actual,expected=float(actual),float(expected)
    assert np.isfinite(actual) and np.isfinite(expected), label
    error=abs(actual-expected)
    assert error<1e-13, (label,actual,expected,error)
    return error


def equal_value(actual, prior, label):
    if actual==prior:
        return
    try:
        close(actual,prior,label)
    except (ValueError,TypeError):
        raise AssertionError((label,actual,prior))


def independent_area(mask):
    streak=np.zeros(mask.shape[1],dtype=np.int64)
    numerator=np.zeros(mask.shape[1],dtype=np.int64)
    for row in mask:
        streak=np.where(row,streak+1,0)
        numerator+=np.where(row,2*streak-1,0)
    # Literal user-confirmed width protects against a silently changed config.
    return numerator / 53**2


def verify_normalization(expected_runs=None):
    assert LUMEN_ROWS==(18,71) and LUMEN_HEIGHT_PX==53 and DIAMETER_PX==53.0
    assert FIXED_PROBE_ROW==40
    assert 'User-confirmed' in BOUNDARY_PROVENANCE and '18:71' in BOUNDARY_PROVENANCE
    selected=set(selected_runs()) if expected_runs is None else set(expected_runs)
    old_void=table(BASELINE/'results/void_fraction/panorama_void_per_run.csv')
    new_void=table(artifact_path('panorama_void_per_run.csv'))
    old_phase=table(BASELINE/'results/flip_rate/phase_mask_flip_per_run.csv')
    new_phase=table(artifact_path('phase_mask_flip_per_run.csv'))
    old_inventory=table(BASELINE/'results/validation/data_inventory.csv')
    new_inventory=table(artifact_path('data_inventory.csv'))
    assert set(old_void)==set(new_void)==set(old_phase)==set(new_phase)==selected
    assert set(old_inventory)==set(new_inventory)==selected
    factor=(60/53)**2
    projected_factor=60/53
    errors=[];records=[]
    for run in sorted(selected):
        old_path=BASELINE/'intermediate/phase_masks'/f'{run}_mask.png'
        new_path=artifact_path('phase_masks')/f'{run}_mask.png'
        old_mask=np.asarray(Image.open(old_path))>0
        new_mask=np.asarray(Image.open(new_path))>0
        assert np.array_equal(old_mask,new_mask), f'Normalization changed mask pixels: {run}'
        assert sha256(old_path)==sha256(new_path), f'Normalization changed mask encoding: {run}'
        assert not new_mask[:18].any() and not new_mask[71:].any(), f'Mask outside confirmed lumen: {run}'
        old,new=old_void[run],new_void[run]
        assert REMOVED_FIELDS.isdisjoint(new), f'Obsolete outside-lumen scenarios remain: {run}'
        assert int(new['lumen_start_row'])==18 and int(new['lumen_stop_row_exclusive'])==71
        assert float(new['lumen_diameter_px'])==53 and new['boundary_source']==BOUNDARY_PROVENANCE
        for field in AREA_FIELDS:
            errors.append(close(new[field],float(old[field])*factor,f'{run}:{field}'))
            assert 0<=float(new[field])<=1+1e-13, f'{run}:{field} outside 0..1'
        errors.append(close(new['gas_silhouette_projected_fraction'],
                            float(old['gas_silhouette_projected_fraction'])*projected_factor,
                            f'{run}:projected fraction'))
        for field in UNSCALED_VOID_FIELDS:
            equal_value(new[field],old[field],f'{run}:{field}')
        for field in PHASE_FIELDS:
            equal_value(new_phase[run][field],old_phase[run][field],f'{run}:phase:{field}')
        per_column=independent_area(new_mask)
        close(np.mean(per_column),new['conditional_circle_alpha_spatial_mean'],f'{run}:independent D53 mean')
        close(np.median(per_column),new['conditional_circle_alpha_spatial_median'],f'{run}:independent D53 median')
        close(new_mask.sum()/(53*new_mask.shape[1]),new['gas_silhouette_projected_fraction'],f'{run}:independent projection')
        assert np.all((per_column>=0)&(per_column<=1)), f'Column area outside 0..1: {run}'
        flips=int(np.count_nonzero(new_mask[18:71,1:]!=new_mask[18:71,:-1]))
        assert flips==int(new_phase[run]['candidate_pixel_flip_count'])
        source=DATA/run/'__panorama_v11.png'
        assert sha256(source)==old_inventory[run]['sha256']==new_inventory[run]['sha256']
        records.append(dict(run=run,pattern=new['label'].replace('_',' '),eligible=new_phase[run]['eligible']=='True',
            diameter_before_px=60,diameter_after_px=53,lumen_rows='18:71',
            alpha_before_D60=float(old['conditional_circle_alpha_spatial_mean']),
            alpha_after_D53=float(new['conditional_circle_alpha_spatial_mean']),alpha_scale_factor=factor,
            projected_before_D60=float(old['gas_silhouette_projected_fraction']),
            projected_after_D53=float(new['gas_silhouette_projected_fraction']),projected_scale_factor=projected_factor,
            flip_count_before=int(old_phase[run]['candidate_pixel_flip_count']),flip_count_after=flips,
            flip_rate_before=float(old_phase[run]['candidate_pixel_flip_rate']),
            flip_rate_after=float(new_phase[run]['candidate_pixel_flip_rate']),changed_mask_pixels=0,
            before_mask_sha256=sha256(old_path),after_mask_sha256=sha256(new_path),
            input_panorama_sha256=sha256(source)))
    result=dict(passed=True,completed_utc=datetime.now(timezone.utc).isoformat(),
        selected_acquisitions=len(records),statistical_acquisitions=sum(r['eligible'] for r in records),
        diameter_before_px=60,diameter_after_px=53,confirmed_lumen_rows=[18,71],full_lumen_maximum=1,
        boundary_source=BOUNDARY_PROVENANCE,alpha_scale_factor=factor,projected_scale_factor=projected_factor,
        identical_masks=len(records),unchanged_flip_metrics=len(records),unchanged_input_panorama_hashes=len(records),
        scaled_value_checks=len(errors),max_rescaling_error=max(errors),obsolete_primary_fields_removed=True,
        segmentation_source_sha256=sha256(artifact_path('panorama_void_alternatives.py')),
        geometry_source_sha256=sha256(artifact_path('v14_geometry.py')),failures=[],
        interpretation='Only normalization changes. User-confirmed lumen boundaries are provenance, not independent optical/phase calibration.')
    return result,records


def main():
    result,records=verify_normalization()
    artifact_path('lumen_normalization_checks.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    with artifact_path('lumen_normalization_changes.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    lines=['# Verificação da normalização pelo lúmen confirmado','',
        '**PASS:** o usuário confirmou que a fatia Python 18:71 é todo o interior do tubo. São 53 células de pixel, portanto D=53 pixels. A declaração do usuário fornece a geometria de referência; não calibra refração ou forma tridimensional do gás.','',
        'As máscaras e as métricas de flip-rate ficaram idênticas às anteriores. A mudança de α é exclusivamente a correção do denominador: α₅₃ = α₆₀(60/53)². A fração projetada muda por 60/53.','',
        f'Foram verificados {result["identical_masks"]} arquivos de máscara, {result["unchanged_flip_metrics"]} conjuntos de métricas de flip-rate e {result["unchanged_input_panorama_hashes"]} hashes de imagens originais. Maior erro numérico de reescala: {result["max_rescaling_error"]:.3g}.','',
        '| Padrão | n | α anterior, D=60 | α atual, D=53 | Flip-rate atual / 1.000 |','|---|---:|---:|---:|---:|']
    for pattern in PATTERN_ORDER:
        rows=[r for r in records if r['pattern']==pattern and r['eligible']]
        lines.append(f'| {pattern} | {len(rows)} | {np.median([r["alpha_before_D60"] for r in rows]):.6f} | {np.median([r["alpha_after_D53"] for r in rows]):.6f} | {1000*np.median([r["flip_rate_after"] for r in rows]):.3f} |')
    lines+=['','O máximo matemático atual é 1. O antigo máximo 0,7803 era consequência do denominador D=60 e foi superado. Valores perto de 1 não provam gás puro: a máscara e a reconstrução circular continuam sujeitas a erro.','',
        'Os antigos intervalos de extensão para fora do recorte e as alternativas D=56/64 foram retirados dos resultados principais. A região externa à fatia confirmada não integra o lúmen. Recortes internos menores servem apenas para medir o efeito de descartar linhas conhecidas do interior.','',
        'A correção anterior de borda continua separada na [revisão morfológica](deep_method_review.md), com ambos os lados daquela comparação em D=60.','',
        'Detalhes: [comparação por aquisição](lumen_normalization_changes.csv), [verificações](lumen_normalization_checks.json) e [código](lumen_normalization_checks.py).']
    write_markdown(artifact_path('lumen_normalization_review.md'),'\n'.join(lines)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
