#!/usr/bin/env python3
"""Independent morphology, geometry and before/after checks; run after outputs."""
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json

import numpy as np
from PIL import Image
from scipy.ndimage import binary_fill_holes, gaussian_filter1d, label

from v14_paths import ROOT, artifact_path, write_markdown
from v14_selection import iter_runs, is_eligible, PATTERN_ORDER
from v14_geometry import (LUMEN_ROWS, LUMEN_SLICE, LUMEN_HEIGHT_PX, DIAMETER_PX,
                          FIXED_PROBE_ROW, BOUNDARY_PROVENANCE)
import panorama_void_alternatives as estimator

HISTORY = ROOT / 'intermediate/history/before_deep_method_review'
MORPHOLOGY_COMPARISON_DIAMETER = 60.0


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def independent_closing(crop):
    """3x3 dilation then erosion on an infinite zero-extended binary plane."""
    h, w = crop.shape
    padded = np.pad(np.asarray(crop, dtype=bool), 2)
    dilated = np.logical_or.reduce([padded[dy:dy+h+2, dx:dx+w+2]
                                  for dy in range(3) for dx in range(3)])
    return np.logical_and.reduce([dilated[dy:dy+h, dx:dx+w]
                                  for dy in range(3) for dx in range(3)])


def direct_loop_closing(crop):
    h, w = crop.shape
    padded = np.pad(np.asarray(crop, dtype=bool), 2)
    dilated = np.zeros((h+2, w+2), dtype=bool)
    for y in range(h+2):
        for x in range(w+2):
            dilated[y, x] = any(bool(padded[y+dy, x+dx]) for dy in range(3) for dx in range(3))
    closed = np.zeros((h,w), dtype=bool)
    for y in range(h):
        for x in range(w):
            closed[y, x] = all(bool(dilated[y+dy, x+dx]) for dy in range(3) for dx in range(3))
    return closed


def isolated_shadow_candidates(closed):
    components, count = label(closed)
    radial_labels = np.unique(np.r_[components[0], components[-1]])
    sizes = np.bincount(components.ravel(), minlength=count+1)
    retained = (~np.isin(components, radial_labels)) & (sizes[components] >= 4)
    return binary_fill_holes(retained)


def independent_mask(im, reference, threshold=.70, rows=(18,71), min_open_interval=20):
    normalized = im / np.maximum(np.percentile(im[30:50],75,axis=0),1)[None,:]
    ratio = gaussian_filter1d(normalized,1,axis=0) / np.maximum(reference[:,None],.01)
    cols = np.arange(im.shape[1])
    upper = 19 + np.argmin(ratio[19:38],axis=0)
    lower = 42 + np.argmin(ratio[42:71],axis=0)
    accepted = (ratio[upper,cols]<threshold) & (ratio[lower,cols]<threshold)
    shadow = ratio<threshold
    closed = independent_closing(shadow[rows[0]:rows[1]])
    components,_ = label(closed)
    top = components[np.clip(upper-rows[0],0,closed.shape[0]-1),cols]
    bottom = components[np.clip(lower-rows[0],0,closed.shape[0]-1),cols]
    same = (top>0) & (top==bottom)
    intervals,_ = label(accepted)
    sizes = np.bincount(intervals)
    accepted &= same | (sizes[intervals]>=min_open_interval)
    outer_top, outer_bottom = upper.copy(),lower.copy()
    for _ in range(32):
        go = (outer_top>rows[0]) & shadow[np.maximum(outer_top-1,0),cols]
        outer_top[go]-=1
        go = (outer_bottom<rows[1]-1) & shadow[np.minimum(outer_bottom+1,79),cols]
        outer_bottom[go]+=1
    mask = np.zeros_like(im,dtype=bool)
    for y in range(rows[0],rows[1]):
        mask[y]=accepted & (y>=outer_top) & (y<=outer_bottom)
    mask[rows[0]:rows[1]] |= isolated_shadow_candidates(closed)
    return mask


def area_by_intervals(mask, diameter=DIAMETER_PX):
    """Independent interval endpoints, rather than the estimator accumulator."""
    transitions = np.diff(np.pad(mask.astype(np.int8),((1,1),(0,0))),axis=0)
    result=np.zeros(mask.shape[1],dtype=float)
    for x in range(mask.shape[1]):
        starts=np.flatnonzero(transitions[:,x]==1)
        stops=np.flatnonzero(transitions[:,x]==-1)
        result[x]=np.sum((stops-starts)**2)/diameter**2
    return result


def flips(arr):
    return int(np.count_nonzero(arr[...,1:]!=arr[...,:-1]))


def analytic_controls():
    assert LUMEN_ROWS == (18,71) and LUMEN_HEIGHT_PX == 53 and DIAMETER_PX == 53.0
    assert estimator.DIAMETER_PX == 53.0 and FIXED_PROBE_ROW == 40
    cases=[np.zeros((9,12),bool),np.ones((9,12),bool)]
    top=np.zeros((9,12),bool);top[0:5,4:9]=True;cases.append(top)
    left=np.zeros((9,12),bool);left[3:7,0:4]=True;cases.append(left)
    ring=np.zeros((11,13),bool);ring[2:9,3:10]=True;ring[3:8,4:9]=False;cases.append(ring)
    rng=np.random.default_rng(20261023)
    cases.extend(rng.random((9,12))<p for p in (.05,.2,.5,.8,.95))
    for case in cases:
        expected=direct_loop_closing(case)
        assert np.array_equal(independent_closing(case),expected)
        assert np.array_equal(estimator.close_shadow_crop(case),expected)
        assert np.all(~case | expected), 'Zero-extended closing must be extensive'
    assert not isolated_shadow_candidates(estimator.close_shadow_crop(top)).any()
    assert np.array_equal(isolated_shadow_candidates(estimator.close_shadow_crop(left)),left)
    assert isolated_shadow_candidates(estimator.close_shadow_crop(left))[:,0].any()
    assert not isolated_shadow_candidates(estimator.close_shadow_crop(np.ones((9,12),bool))).any()
    ring_filled=isolated_shadow_candidates(estimator.close_shadow_crop(ring))
    assert ring_filled[5,6] and int(ring_filled.sum())==49
    for level in (50,100,200):
        assert not estimator.silhouette(np.full((80,120),level,dtype=float),np.ones(80))[0].any()
    noise=100+rng.uniform(-1,1,size=(80,120))
    assert not estimator.silhouette(noise,np.ones(80))[0].any()
    reference,_=estimator.make_reference()
    assert not estimator.silhouette(np.tile(100*reference[:,None],(1,120)),reference)[0].any()
    body=np.full((80,140),100.,dtype=float)
    body[24:28,40:100]=25;body[62:66,40:100]=25
    body_mask=estimator.silhouette(body,np.ones(80))[0]
    assert body_mask[40,70] and not body_mask[:,10].any() and flips(body_mask[40])==2
    empty=np.zeros((80,10),bool)
    full=empty.copy();full[18:71]=True
    separate=empty.copy();separate[20:30]=True;separate[50:55]=True
    for mask in (empty,full,separate):
        assert np.allclose(estimator.geometry(mask),area_by_intervals(mask),atol=1e-15,rtol=0)
    assert np.allclose(estimator.geometry(full),1.0)
    assert np.allclose(estimator.geometry(separate),(10**2+5**2)/53**2)
    assert np.allclose(area_by_intervals(full,60),(53/60)**2)
    return dict(morphology_reference_cases=len(cases),direct_loop_and_boolean_reduction_agree=True,
                radial_contact_rejected=True,axial_contact_preserved=True,
                interior_hole_filled=True,flat_intensity_controls=3,
                bounded_noise_control=True,liquid_profile_control=True,bright_gas_interior_control=True,
                separated_interval_geometry_control=True,full_lumen_maximum=1.0,
                confirmed_lumen_rows=list(LUMEN_ROWS),confirmed_diameter_px=DIAMETER_PX,
                historical_D60_full_mask_value=(53/60)**2,
                boundary_source=BOUNDARY_PROVENANCE)


def table(path):
    with path.open(newline='') as fp:return list(csv.DictReader(fp))


def main():
    controls=analytic_controls()
    old_source=HISTORY/'scripts/panorama_void_alternatives.py'
    spec=importlib.util.spec_from_file_location('pre_review_estimator',old_source)
    old_estimator=importlib.util.module_from_spec(spec);spec.loader.exec_module(old_estimator)
    reference,_=estimator.make_reference()
    old_reference,_=old_estimator.make_reference()
    assert np.array_equal(reference,old_reference), 'Reference construction unexpectedly changed'
    current_void={r['run']:r for r in table(artifact_path('panorama_void_per_run.csv'))}
    current_flips={r['run']:r for r in table(artifact_path('phase_mask_flip_per_run.csv'))}
    measurements=[]
    for folder in iter_runs():
        run=folder.name;im=estimator.read(run)
        prior_path=HISTORY/'intermediate/phase_masks'/f'{run}_mask.png'
        current_path=artifact_path('phase_masks')/f'{run}_mask.png'
        prior=np.asarray(Image.open(prior_path))>0
        current=np.asarray(Image.open(current_path))>0
        old=old_estimator.silhouette(im,old_reference)[0]
        independent=independent_mask(im,reference)
        assert np.array_equal(prior,old), f'{run}: archived mask differs from archived code'
        assert np.array_equal(current,independent), f'{run}: saved mask differs from independent correction'
        assert np.array_equal(estimator.silhouette(im,reference)[0],independent), f'{run}: corrected estimator disagrees'
        # Compare BOTH masks at D60 to isolate the historical morphology change.
        alpha_old=float(area_by_intervals(prior,MORPHOLOGY_COMPARISON_DIAMETER).mean())
        alpha_new=float(area_by_intervals(current,MORPHOLOGY_COMPARISON_DIAMETER).mean())
        current_alpha=float(area_by_intervals(current,DIAMETER_PX).mean())
        denominator=53*(im.shape[1]-1)
        count_old=flips(prior[18:71]);count_new=flips(current[18:71])
        assert abs(current_alpha-float(current_void[run]['conditional_circle_alpha_spatial_mean']))<1e-14
        assert abs(current_alpha-alpha_new*(60/53)**2)<1e-14
        assert count_new==int(current_flips[run]['candidate_pixel_flip_count'])
        pattern=next(folder.glob('__*.txt')).stem.lstrip('_').replace('_',' ')
        measurements.append(dict(run=run,pattern=pattern,eligible=is_eligible(run,pattern,float((im[20:60]==255).mean())),
            alpha_before=alpha_old,alpha_after=alpha_new,alpha_change=alpha_new-alpha_old,
            morphology_comparison_diameter_px=MORPHOLOGY_COMPARISON_DIAMETER,
            current_alpha_D53=current_alpha,current_normalization_diameter_px=DIAMETER_PX,
            flip_count_before=count_old,flip_count_after=count_new,flip_comparisons=denominator,
            flips_per1000_before=1000*count_old/denominator,flips_per1000_after=1000*count_new/denominator,
            changed_pixels=int((prior!=current).sum()),added_pixels=int((~prior & current).sum()),
            removed_pixels=int((prior & ~current).sum()),
            fixed_row40_flips_per1000=1000*flips(current[40])/(im.shape[1]-1),
            lumen_midpoint_row44_flips_per1000=1000*flips(current[44])/(im.shape[1]-1),
            input_panorama_sha256=sha256(folder/'__panorama_v11.png'),
            before_mask_sha256=sha256(prior_path),after_mask_sha256=sha256(current_path)))
    with artifact_path('deep_method_changes.csv').open('w',newline='') as fp:
        writer=csv.DictWriter(fp,fieldnames=list(measurements[0]));writer.writeheader();writer.writerows(measurements)
    groups=[]
    for pattern in PATTERN_ORDER:
        rows=[r for r in measurements if r['pattern']==pattern and r['eligible']]
        groups.append(dict(pattern=pattern,n=len(rows),
            alpha_before=float(np.median([r['alpha_before'] for r in rows])),
            alpha_after=float(np.median([r['alpha_after'] for r in rows])),
            flips_per1000_before=float(np.median([r['flips_per1000_before'] for r in rows])),
            flips_per1000_after=float(np.median([r['flips_per1000_after'] for r in rows]))))
    result=dict(passed=True,completed_utc=datetime.now(timezone.utc).isoformat(),controls=controls,
        selected_acquisitions=len(measurements),statistical_acquisitions=sum(r['eligible'] for r in measurements),
        changed_masks=sum(r['changed_pixels']>0 for r in measurements),
        max_abs_alpha_change=max(abs(r['alpha_change']) for r in measurements),
        morphology_comparison_diameter_px=MORPHOLOGY_COMPARISON_DIAMETER,
        current_normalization_diameter_px=DIAMETER_PX,boundary_source=BOUNDARY_PROVENANCE,
        source_before_sha256=sha256(old_source),source_after_sha256=sha256(artifact_path('panorama_void_alternatives.py')),
        segmentation_source_sha256=sha256(artifact_path('panorama_void_alternatives.py')),
        group_comparisons=groups,measurements=measurements,failures=[],
        interpretation='Morphology ledger compares both masks at historical D60; current D53 is independently verified. User-supplied boundaries do not calibrate optical refraction or gas shape.')
    artifact_path('deep_method_checks.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    write_report(result)
    print(json.dumps({k:result[k] for k in ('passed','selected_acquisitions','statistical_acquisitions','changed_masks','max_abs_alpha_change')},indent=2))


def write_report(result):
    c=result['controls']
    lines=['# Revisão adversarial da segmentação e das equações','',
        '**Resultado: erro real de borda corrigido e verificado.** Os valores anteriores estão preservados no histórico. Os controles verificam o algoritmo; não fornecem calibração física de gás/líquido.','',
        '## Erro de fechamento morfológico','',
        'O fechamento 3×3 era aplicado diretamente ao recorte. Na erosão, o valor zero fora da matriz apagava a primeira e a última linha. A etapa seguinte tentava rejeitar componentes nessas linhas já apagadas; assim, nunca rejeitava um componente de primeiro plano por contato radial.','',
        'Correção: margem de um pixel com zeros → dilatação → erosão → retirada da margem → identificação de componentes → rejeição dos que tocam uma borda radial. Contatos axiais continuam permitidos. O tamanho mínimo de quatro pixels e o preenchimento de buracos foram mantidos.','',
        f'Verificação independente: {c["morphology_reference_cases"]} matrizes compararam laços explícitos, operações booleanas e a implementação corrigida. Controles adicionais cobrem contato radial/axial, buracos, iluminação uniforme, ruído limitado, perfil líquido e interior claro de gás. Todos passaram.','',
        f'{result["changed_masks"]} de {result["selected_acquisitions"]} máscaras mudaram na correção morfológica. Maior variação absoluta de α nessa correção, mantendo D=60 em ambos os lados: {result["max_abs_alpha_change"]:.6f}. A tabela histórica abaixo isola a morfologia e usa as {result["statistical_acquisitions"]} aquisições elegíveis. A atualização atual para D=53 é verificada separadamente.','',
        '| Padrão | n | α antes, D=60 | α após morfologia, D=60 | Trocas/1.000 antes | Trocas/1.000 corrigido |',
        '|---|---:|---:|---:|---:|---:|']
    for g in result['group_comparisons']:
        lines.append(f'| {g["pattern"]} | {g["n"]} | {g["alpha_before"]:.6f} | {g["alpha_after"]:.6f} | {g["flips_per1000_before"]:.3f} | {g["flips_per1000_after"]:.3f} |')
    lines += ['', 'São medianas por padrão. [Todas as alterações por aquisição](deep_method_changes.csv).','',
        '**Limite:** cumprir o critério de borda não prova que os componentes restantes sejam gás. Sombras do líquido podem permanecer; bolhas reais conectadas à borda podem ser rejeitadas. Os componentes aceitos não precisam conter um buraco fechado.','',
        '## Equações e limites','',
        'Reconstrução atual: α̂ = (1/W) Σₓ Σᵢ (dᵢ(x)/D)², com D = 53 pixels. O usuário confirmou que a fatia Python 18:71 contém todo o interior do tubo: linhas 18 a 70, total de 53 células de pixel. O quadrado vem da razão entre áreas circulares. Cada intervalo radial separado é tratado como um disco. Máscaras vazia, completa e com dois intervalos separados verificaram essa soma dos quadrados.','',
        'Flip-rate: F = ΣᵧΣₓ₌₂ᵂ 1[M(y,x)≠M(y,x−1)] / [53(W−1)]. Multiplica-se por 1.000 apenas para apresentação. São transições espaciais. A forma tridimensional e a identidade das fases continuam sendo hipóteses da máscara.','',
        'O valor máximo atual é (53/53)² = 1. O antigo valor 0,7803 resultava de dividir o mesmo interior de 53 linhas por D=60 e foi superado pela informação do usuário. Para máscaras idênticas, α₅₃ = α₆₀(60/53)²; os flip-rates permanecem iguais. Uma máscara cheia ainda não prova gás puro: segmentação, refração e geometria circular continuam condicionais.','',
        'Não se aplicam extensões de gás para fora da fatia confirmada: essa região fica fora do lúmen definido pelo usuário. Os antigos campos de extensão e as alternativas D=56/64 foram retirados dos resultados principais. A normalização está documentada na [verificação do lúmen](lumen_normalization_review.md).','',
        '## Sonda e independência','',
        'A sonda permanece na linha fixa y=40. O ponto médio das linhas confirmadas 18 a 70 é y=44; os mínimos ópticos antigos em y=15/75 não definem o lúmen. A tabela inclui a sensibilidade entre y=40 e y=44. A confirmação geométrica do usuário não calibra a refração.','',
        'A imagem de bolhas ficou fora da construção da referência, mas foi inspecionada no desenvolvimento do detector. É um controle interno. Os trechos líquidos inferidos no mesmo conjunto também não são verdade física independente.','',
        'Detalhes: [verificações completas](deep_method_checks.json) e [código](deep_method_checks.py).']
    write_markdown(artifact_path('deep_method_review.md'),'\n'.join(lines)+'\n')


if __name__=='__main__':
    main()
