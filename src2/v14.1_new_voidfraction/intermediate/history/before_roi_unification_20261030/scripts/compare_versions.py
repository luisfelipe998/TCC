"""Quantify the seam correction and report all existing threshold settings."""
import csv
from itertools import combinations
import json
import os

import numpy as np
import pandas as pd
from PIL import Image
from scipy.optimize import linprog
from scipy.spatial import ConvexHull

from v14_paths import ROOT, BASELINE, CACHE, artifact_path, original_panorama_path, write_markdown
from v14_selection import PATTERN_ORDER
import panorama_void_alternatives as estimator

os.environ.setdefault('MPLCONFIGDIR', str(CACHE/'matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from summary_figures import COLORS, PATTERN_LABELS_PT


def read(root, relative):
    return pd.read_csv(root/relative, float_precision='round_trip')


def hull_intersections(frame):
    """Half-space feasibility, independent of the summary's vertex combination test."""
    polygons={}
    for pattern,group in frame.groupby('pattern'):
        points=group[['alpha','F1000']].to_numpy()
        if len(points)>=3:
            polygons[pattern]=ConvexHull(points).equations
    pairs=[]
    for a,b in combinations(PATTERN_ORDER,2):
        if a not in polygons or b not in polygons:continue
        bounds=np.vstack([polygons[a],polygons[b]])
        result=linprog([0.,0.],A_ub=bounds[:,:2],b_ub=-bounds[:,2],bounds=[(None,None)]*2,method='highs')
        if result.success:pairs.append([a,b])
        elif result.status!=2:raise RuntimeError(result.message)
    return pairs


def threshold_coordinates(root, condition):
    void=read(root,'intermediate/sensitivity/panorama_void_sensitivity.csv')
    flip=read(root,'intermediate/sensitivity/phase_mask_flip_sensitivity.csv')
    void=void[void.variant.str.startswith('shadow_ratio_')].copy()
    void['shadow_threshold']=void.variant.str.removeprefix('shadow_ratio_').astype(float)
    merged=void.merge(flip,on=['run','shadow_threshold'],validate='one_to_one')
    assert len(merged)==130 and set(merged.shadow_threshold)=={.6,.65,.7,.75,.8}
    merged=merged[merged.eligible].copy()
    assert len(merged)==125
    merged['condition']=condition
    merged['alpha']=merged.alpha_conditional
    merged['F1000']=1000*merged.candidate_pixel_flip_rate
    return merged[['condition','run','pattern','shadow_threshold','alpha','F1000']]


def main():
    old=read(BASELINE,'results/void_fraction/panorama_void_per_run.csv').set_index('run')
    new=read(ROOT,'results/void_fraction/panorama_void_per_run.csv').set_index('run')
    old_f=read(BASELINE,'results/flip_rate/phase_mask_flip_per_run.csv').set_index('run')
    new_f=read(ROOT,'results/flip_rate/phase_mask_flip_per_run.csv').set_index('run')
    assert old.index.equals(new.index) and old.index.equals(new_f.index) and old.index.equals(old_f.index)
    assert old_f.eligible.equals(new_f.eligible), 'Report any eligibility change; do not silently change the comparison population'
    old_reference=np.array(json.loads((BASELINE/'intermediate/features/panorama_liquid_reference.json').read_text())['radial_reference'])
    records=[]
    for run in new.index:
        image=estimator.read(run)
        fixed,_,_=estimator.silhouette(image,old_reference)
        with Image.open(original_panorama_path(run)) as source:
            original=np.asarray(source,dtype=float)
        old_mask,_,_=estimator.silhouette(original,old_reference)
        with Image.open(BASELINE/'intermediate/phase_masks'/f'{run}_mask.png') as encoded:
            assert np.array_equal(old_mask,np.asarray(encoded)==255), 'Baseline is not reproducible under the preserved method'
        row=dict(run=run,pattern=new_f.loc[run,'pattern'],eligible=bool(new_f.loc[run,'eligible']),
                 alpha_v14=float(old.loc[run,'conditional_circle_alpha_spatial_mean']),
                 alpha_v141=float(new.loc[run,'conditional_circle_alpha_spatial_mean']),
                 F1000_v14=float(1000*old_f.loc[run,'candidate_pixel_flip_rate']),
                 F1000_v141=float(1000*new_f.loc[run,'candidate_pixel_flip_rate']),
                 alpha_corrected_frozen_reference=float(estimator.geometry(fixed).mean()),
                 F1000_corrected_frozen_reference=float(1000*np.count_nonzero(fixed[18:71,1:]!=fixed[18:71,:-1])/(53*(fixed.shape[1]-1))))
        row.update(alpha_change=row['alpha_v141']-row['alpha_v14'],F1000_change=row['F1000_v141']-row['F1000_v14'])
        row['flip_relative_change_percent']=100*row['F1000_change']/row['F1000_v14'] if row['F1000_v14'] else None
        records.append(row)
    comparison=pd.DataFrame(records)
    comparison.to_csv(artifact_path('fix_comparison_per_run.csv'),index=False)
    eligible=comparison[comparison.eligible]
    grouped=[]
    for pattern in PATTERN_ORDER:
        group=eligible[eligible.pattern==pattern]
        grouped.append(dict(pattern=pattern,n=len(group),**{key:float(group[key].median()) for key in
            ['alpha_v14','alpha_v141','F1000_v14','F1000_v141']}))
    pd.DataFrame(grouped).to_csv(artifact_path('fix_comparison_by_pattern.csv'),index=False)
    coordinates=pd.concat([threshold_coordinates(BASELINE,'v14'),threshold_coordinates(ROOT,'v14.1')],ignore_index=True)
    coordinates.to_csv(artifact_path('threshold_plot_data.csv'),index=False)
    sensitivity=[]
    for (condition,threshold),frame in coordinates.groupby(['condition','shadow_threshold']):
        pairs=hull_intersections(frame)
        sensitivity.append(dict(condition=condition,threshold=float(threshold),
            overlap_pairs=json.dumps(pairs),plug_slug_intersect=['plug','slug'] in pairs,
            churn_annular_intersect=['churn','wavy annular'] in pairs,n_acquisitions=len(frame)))
    pd.DataFrame(sensitivity).to_csv(artifact_path('threshold_hull_intersections.csv'),index=False)
    # Temporal features are unaffected by a spatial reconstruction correction.
    native_before=read(BASELINE,'results/flip_rate/flip_rate_per_run.csv').set_index('run')
    native_after=read(ROOT,'results/flip_rate/flip_rate_per_run.csv').set_index('run')
    native_fields=[c for c in native_before if c.startswith('raw_') or c in
        ['n_native_frames','native_camera_fps','frame_index_step','effective_native_sample_hz','elapsed_s','sample_exposure_window_s']]
    pd.testing.assert_frame_equal(native_before[native_fields],native_after[native_fields],check_exact=True)
    classifier_rows=[]
    for family,filename in [('morphology','pattern_identification.json'),('silhouette','silhouette_pattern_results.json')]:
        prior=json.loads((BASELINE/'results/patterns'/filename).read_text())['results']
        current=json.loads((ROOT/'results/patterns'/filename).read_text())['results']
        key=lambda r:(r['model'],r['validation'],r.get('variant','baseline'))
        mapping={key(r):r for r in prior}
        assert len(mapping)==len(prior)==len(current)
        for row in current:
            before=mapping[key(row)]
            classifier_rows.append(dict(family=family,model=row['model'],validation=row['validation'],variant=row.get('variant','baseline'),
                accuracy_v14=before['accuracy'],accuracy_v141=row['accuracy'],
                balanced_accuracy_v14=before['balanced_accuracy'],balanced_accuracy_v141=row['balanced_accuracy'],
                delta_balanced_accuracy=row['balanced_accuracy']-before['balanced_accuracy']))
    classifiers=pd.DataFrame(classifier_rows)
    classifiers.to_csv(artifact_path('classifier_comparison.csv'),index=False)
    fixed=classifiers[(classifiers.variant=='baseline')&classifiers.validation.isin(['leave_one_run_out','leave_one_mass_flux_out'])]
    result=dict(passed=True,n_selected=len(comparison),n_eligible=len(eligible),
        primary_threshold=.70,lumen_rows=[18,71],diameter_px=53,
        changed_alpha_runs=int((comparison.alpha_change.abs()>1e-14).sum()),
        changed_flip_runs=int((comparison.F1000_change.abs()>1e-12).sum()),
        max_absolute_alpha_change=float(eligible.alpha_change.abs().max()),
        max_absolute_F1000_change=float(eligible.F1000_change.abs().max()),
        min_relative_flip_change_percent=float(eligible.flip_relative_change_percent.min()),
        max_relative_flip_change_percent=float(eligible.flip_relative_change_percent.max()),
        max_absolute_alpha_change_frozen_reference=float((eligible.alpha_corrected_frozen_reference-eligible.alpha_v14).abs().max()),
        native_temporal_fields_unchanged=native_fields,group_medians=grouped,threshold_hull_intersections=sensitivity,
        fixed_classifier_comparisons=len(fixed),
        fixed_classifiers_improved=int((fixed.delta_balanced_accuracy>1e-12).sum()),
        fixed_classifiers_worsened=int((fixed.delta_balanced_accuracy< -1e-12).sum()),
        fixed_classifiers_unchanged=int((fixed.delta_balanced_accuracy.abs()<=1e-12).sum()),
        effect='The column-order defect is removed. Numerical changes do not establish better segmentation, physical void accuracy or classification.',
        threshold_sensitivity='All five existing thresholds are reported; 0.70 remains primary. No setting is selected for better separation.')
    artifact_path('fix_effectiveness.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    plots(eligible,coordinates)
    lines=['# Efeito da correção das emendas na v14.1','',
        'A v14.1 preservou as tiras e os deslocamentos da v14 e retirou a sobrescrita de colunas nas emendas. '
        'A referência líquida foi reconstruída com a mesma regra. O controle com a referência antiga fixa isolou a alteração da entrada.','',
        f'Nas {len(eligible)} aquisições elegíveis, a maior mudança absoluta de α foi {result["max_absolute_alpha_change"]:.6f}; '
        f'a mudança relativa de flip-rate variou de {result["min_relative_flip_change_percent"]:.2f}% a {result["max_relative_flip_change_percent"]:.2f}%.','',
        'A correção foi efetiva para remover a inversão/duplicação local de colunas. A exatidão física da máscara e do volume de gás continua sem calibração. '
        'O valor maior ou menor de flip-rate não é, por si só, um critério de melhora.','',
        '| Padrão | n | α v14 | α v14.1 | F×1000 v14 | F×1000 v14.1 |','|---|---:|---:|---:|---:|---:|']
    lines += [f'| {r["pattern"]} | {r["n"]} | {r["alpha_v14"]:.6f} | {r["alpha_v141"]:.6f} | {r["F1000_v14"]:.4f} | {r["F1000_v141"]:.4f} |' for r in grouped]
    lines += ['', 'Valores por padrão são medianas, com peso igual por aquisição.','',
        '![Mudanças por aquisição](fix_comparison.png)','',
        '## Dependência do limiar','',
        'As envoltórias observadas foram recalculadas em todos os limiares já previstos. Elas não são fronteiras de classificação.','',
        '| Versão | Limiar | Interseções |','|---|---:|---|']
    lines += [f'| {r["condition"]} | {r["threshold"]:.2f} | {r["overlap_pairs"]} |' for r in sensitivity]
    lines += ['', '![Sensibilidade do gráfico ao limiar](threshold_sensitivity.png)','',
        f'Nos {len(fixed)} pares fixos de modelo/validação, a acurácia balanceada aumentou em {result["fixed_classifiers_improved"]}, '
        f'diminuiu em {result["fixed_classifiers_worsened"]} e ficou igual em {result["fixed_classifiers_unchanged"]}. '
        'Todos os modelos e variantes constam na tabela; não houve seleção do resultado mais favorável. '
        'A classe de bolhas dispersas continua com uma aquisição e sem exemplo de treino quando ela é retirada.','',
        '[Resultados estruturados](fix_effectiveness.json) · [Por aquisição](fix_comparison_per_run.csv) · '
        '[Classificadores](classifier_comparison.csv) · [Interseções](threshold_hull_intersections.csv)\n']
    write_markdown(artifact_path('fix_effectiveness.md'),'\n'.join(lines))
    print(json.dumps(result,indent=2))


def plots(comparison,coordinates):
    fig,axes=plt.subplots(2,1,figsize=(13,7),sharex=True,layout='constrained')
    colors=[COLORS[p] for p in comparison.pattern]
    x=np.arange(len(comparison))
    for ax,column,label in zip(axes,['alpha_change','F1000_change'],['Δα (v14.1 − v14)','Δ transições / 1.000 pixels']):
        ax.bar(x,comparison[column],color=colors);ax.axhline(0,color='black',lw=.8);ax.set_ylabel(label);ax.grid(axis='y',alpha=.2)
    axes[-1].set_xticks(x,comparison.run,rotation=70,ha='right',fontsize=8)
    fig.suptitle('Efeito da correção da entrada e da reconstrução da referência; mesmos parâmetros')
    fig.savefig(artifact_path('fix_comparison.png'),dpi=180);plt.close(fig)
    current=coordinates[coordinates.condition=='v14.1']
    fig,axes=plt.subplots(1,5,figsize=(18,4.5),sharex=True,sharey=True,layout='constrained')
    for ax,(threshold,group) in zip(axes,current.groupby('shadow_threshold')):
        for pattern in PATTERN_ORDER:
            points=group.loc[group.pattern==pattern,['alpha','F1000']].to_numpy()
            ax.scatter(points[:,0],points[:,1],s=25,color=COLORS[pattern],label=PATTERN_LABELS_PT[pattern])
            if len(points)>=3:
                polygon=points[ConvexHull(points).vertices]
                ax.fill(polygon[:,0],polygon[:,1],alpha=.15,color=COLORS[pattern])
        ax.set_title(f'Limiar {threshold:.2f}'+(' (principal)' if threshold==.7 else ''))
        ax.set(xlim=(-.02,1.02),xlabel='α condicional');ax.grid(alpha=.2)
    axes[0].set_ylabel('Transições / 1.000 pixels')
    handles,labels=axes[-1].get_legend_handles_labels()
    fig.legend(handles,labels,fontsize=9,loc='outside lower center',ncol=len(PATTERN_ORDER))
    fig.suptitle('v14.1: todos os limiares previstos; envoltórias da amostra, sem validação de classificação')
    fig.savefig(artifact_path('threshold_sensitivity.png'),dpi=180);plt.close(fig)


if __name__=='__main__':main()
