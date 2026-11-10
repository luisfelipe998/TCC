"""Bilingual summaries rendered from the same checked numerical tables."""
PATTERN_PT = {'dispersed bubbles': 'Bolhas dispersas', 'plug': 'Plug', 'slug': 'Slug',
              'churn': 'Churn', 'wavy annular': 'Anular ondulado'}


def number(value, places, pt):
    value = f'{float(value):.{places}f}'
    return value.replace('.', ',') if pt else value


def hull_overlap_pairs(chart_checks, plot_rows):
    """Test intersection of observed convex hulls by convex-combination feasibility."""
    import numpy as np
    from scipy.optimize import linprog
    lookup = {r['run']: r for r in plot_rows}
    polygons = {}
    for pattern, hull in chart_checks['primary_chart_hulls'].items():
        if hull['shaded']:
            polygons[pattern] = np.array([[float(lookup[run]['conditional_circle_alpha_spatial_mean']),
                                          float(lookup[run]['candidate_pixel_flips_per1000px'])]
                                         for run in hull['polygon_vertex_runs']])
    pairs = []
    patterns = list(polygons)
    for i, first in enumerate(patterns):
        for second in patterns[i+1:]:
            a, b = polygons[first], polygons[second]
            equality = np.vstack([np.hstack([a.T, -b.T]),
                                  np.r_[np.ones(len(a)), np.zeros(len(b))],
                                  np.r_[np.zeros(len(a)), np.ones(len(b))]])
            result = linprog(np.zeros(len(a)+len(b)), A_eq=equality, b_eq=[0,0,1,1],
                             bounds=(0,None), method='highs')
            if result.success:
                pairs.append((first, second))
            elif result.status != 2:
                raise RuntimeError('Hull intersection could not be established: ' + result.message)
    return pairs


def build_summary(*, pt, groups, phase_groups, n_runs, n_rows, diameter,
                  crop_values, effect, overlap_pairs):
    n = sum(int(r['n']) for r in groups)
    assert n_rows == diameter == 53
    table = ['| Padrão | n | α estimado | Trocas / 1.000 pares de pixels |' if pt else
             '| Pattern | n | Estimated α | Changes / 1,000 pixel pairs |', '|---|---:|---:|---:|']
    for group in groups:
        pattern = group['pattern']
        label = PATTERN_PT[pattern] if pt else pattern
        table.append(f"| {label} | {group['n']} | {number(group['conditional_alpha_median'],4,pt)} | "
                     f"{number(1000*float(phase_groups[pattern]['candidate_pixel_flip_rate_median']),3,pt)} |")
    table_text = '\n'.join(table)
    pairs = '; '.join((PATTERN_PT[a]+' × '+PATTERN_PT[b]) if pt else (a+' × '+b) for a,b in overlap_pairs)
    overlap = (f'Interseções entre as envoltórias observadas: **{pairs}**.' if pairs else 'Não há interseções entre as envoltórias calculadas.') if pt else (
               f'Observed hull intersections: **{pairs}**.' if pairs else 'The observed hulls do not intersect.')
    delta = number(effect['max_absolute_alpha_change'],6,pt)
    flip_low = number(effect['min_relative_flip_change_percent'],2,pt)
    flip_high = number(effect['max_relative_flip_change_percent'],2,pt)
    crop_low,crop_high = [number(v,4,pt) for v in (min(crop_values),max(crop_values))]
    better,worse,same = [effect[k] for k in ('fixed_classifiers_improved','fixed_classifiers_worsened','fixed_classifiers_unchanged')]
    thresholds = ['| Limiar | Interseções na v14.1 |' if pt else '| Threshold | Intersections in v14.1 |','|---:|---|']
    import json
    for row in effect['threshold_hull_intersections']:
        if row['condition'] != 'v14.1':
            continue
        intersections = json.loads(row['overlap_pairs'])
        names = '; '.join((PATTERN_PT[a]+' × '+PATTERN_PT[b]) if pt else (a+' × '+b) for a,b in intersections)
        thresholds.append(f"| {number(row['threshold'],2,pt)} | {names or ('nenhuma' if pt else 'none')} |")
    threshold_text = '\n'.join(thresholds)
    if pt:
        return f'''# Resumo executivo — v14.1: fração de vazio e flip-rate

[English version](EXECUTIVE_SUMMARY.md) · [Efeito das correções](fix_effectiveness.md)

**A inversão/duplicação de colunas nas emendas foi removida. A fração volumétrica real não está calibrada.** A análise estima gás visível sob hipótestes ópticas e geométricas. Foram preservados seleção, fórmulas, limiares, modelos e a geometria confirmada pelo usuário: **linhas Python 18:71, D = 53 pixels**.

## As correções foram efetivas?

A montagem antiga podia transformar `[1,2,3,4 | 5,6,7,8]` em `[1,2,3,6,5,6,7,8]`, gerando três transições em vez de uma ao binarizar por intensidade > 5. A v14.1 concatena as mesmas tiras na ordem original, sem essa sobrescrita. Os {n_runs} panoramas antigos foram reproduzidos pixel a pixel como controle; as entradas novas conservam os deslocamentos e as larguras. A referência líquida foi reconstruída com a mesma regra. Os arquivos originais permanecem preservados.

Nas {n} aquisições elegíveis, **máxima mudança absoluta de α: {delta}**; **mudança relativa de flip-rate: {flip_low}% a {flip_high}%**. A correção removeu o defeito de montagem, mas não demonstra maior exatidão física. A tabela [por aquisição](fix_comparison_per_run.csv) também mantém a referência antiga fixa para isolar o efeito direto das colunas.

Nos {effect['fixed_classifier_comparisons']} pares fixos de modelo/validação, a acurácia balanceada aumentou em **{better}**, diminuiu em **{worse}** e ficou igual em **{same}**. Não houve ajuste para favorecer o gráfico ou os classificadores. Menor flip-rate não significa automaticamente melhor resultado.

## Fração de vazio × flip-rate espacial

![Fração de vazio estimada versus flip-rate espacial](void_fraction_vs_flip_rate.png)

**α [−] é adimensional:** a soma das áreas circulares inferidas de cada intervalo de gás, dividida pela área circular do tubo, seguida da média espacial. Multiplicar por 100 expressa a estimativa em porcentagem. **O eixo Y é 1.000 F**, transições entre pixels vizinhos na direção longitudinal, com média das {n_rows} linhas; **não é uma frequência em Hz**.

O gráfico usa **{n} aquisições elegíveis**, em cinco padrões. A aquisição saturada m03q40_001 fica fora do gráfico e das estatísticas, mas permanece na rastreabilidade dos {n_runs} registros selecionados. As regiões coloridas são envoltórias convexas da amostra, não fronteiras de classificação. **Estrelas: médias por padrão, com peso igual por aquisição. Tabela: medianas.**

{table_text}

{overlap} Esse resultado refere-se ao limiar principal **0,70**. Bolhas dispersas têm **uma aquisição**, m07q01_001, sem informação sobre a variabilidade da classe. Os dois eixos vêm da mesma máscara; sua associação não valida a fração volumétrica real.

## Dependência do limiar

{threshold_text}

![Todos os limiares já previstos, com os mesmos limites de eixos](threshold_sensitivity.png)

A separação entre envoltórias muda com o limiar. Portanto, a separação observada no gráfico principal não é uma conclusão robusta para todos os parâmetros. Foram mantidos 0,70 como principal e os cinco limiares já previstos, sem escolher o mais favorável. As envoltórias descrevem os pontos disponíveis; não constituem teste de classificação nem intervalo de confiança.

## Etapas, métodos e equações

1. **Entrada, seleção e saturação.** Mantivemos os quadros, deslocamentos e tiras centrais; retiramos apenas a sobrescrita defeituosa nas emendas. Usamos os rótulos originais e a [seleção anterior](dataset_selection.json). Mistos, anular liso e m05q05_002 continuam excluídos. A fração de pixels 255 nas linhas 20:60 define S; **S > 0,20** exclui a aquisição das estatísticas. A montagem continua dependendo da hipótese de advecção e dos deslocamentos medidos.
2. **Normalização e referência líquida.** Dividimos cada coluna pelo P75 das linhas 30:50, limitado inferiormente a 1. Nos mesmos três panoramas plug, selecionamos as colunas com variação radial nas linhas 20:68 no percentil 5 ou abaixo; usamos medianas e filtro gaussiano σ = 1 pixel. São trechos interpretados como líquido, sem certificação de fase pura.
3. **Segmentação.** **R = Gσ(I normalizada) / max(L, 0,01)**; **R < 0,70** identifica sombras. Conservamos pares superior/inferior, fechamento 3×3 com borda preservada, componentes sem contato radial com o limite e preenchimento de buracos. Pares desconectados exigem ao menos 20 colunas consecutivas; componentes isolados exigem ao menos 4 pixels. A implementação da segmentação permaneceu a mesma.

![Panoramas e máscaras de gás em azul](panorama_phase_masks.png)

4. **Área circular condicional.** Para intervalos separados de largura dᵢⱼ na coluna j: **αⱼ = Σᵢ(dᵢⱼ/D)²**; **ᾱ = (1/W) Σⱼ αⱼ**, D = 53. A máscara preenchendo todo o interior produz α = 1. Uma projeção lateral não determina a forma tridimensional; circularidade, refração e espessura de filmes permanecem limitações.
5. **Flip-rate espacial.** **F = Σᵣ Σⱼ |Bᵣ,ⱼ₊₁ − Bᵣ,ⱼ| / [H(W−1)]**, H = 53, B = máscara binária. O gráfico apresenta 1.000 F. Larguras das tiras, emendas e resolução ainda influenciam F. A [linha fixa y=40](void_fraction_vs_center_flips.png) é uma comparação, não um centro físico calibrado.
6. **Sensibilidade e controles adversariais.** Repetimos todos os limiares, referências e recortes anteriores. Controles independentes verificaram a ordem das colunas, máscaras, áreas, denominadores e resultados escritos. Não alteramos D para compensar os efeitos da entrada. Para h linhas internas mantidas, o teto matemático permanece (h/53)²; em m06q80_001, os recortes produzem α de {crop_low} a {crop_high}. Isso não é intervalo de confiança.

![Sensibilidade à retirada de linhas internas](void_crop_sensitivity.png)

7. **Classificação e verificação.** Mantivemos modelos, validação por aquisição/grupo de fluxo mássico e referências ajustadas somente no treino. Reconstruímos matrizes de confusão, médias, medianas e coordenadas. A classe de bolhas não possui exemplo no treino quando sua única aquisição é retirada; esse teste continua nos escores. O detector foi desenvolvido nestes dados, portanto a validação é interna. A [literatura](source_audit.md) orienta os métodos, sem fornecer calibração volumétrica destas imagens.

[Auditoria adversarial da correção](seam_adversarial_review.md) · [Verificação após a regeneração](post_regeneration_adversarial_review.md) · [Resultados estruturados](fix_effectiveness.json)
'''
    return f'''# Executive summary — v14.1: void fraction and flip-rate

[Versão em português](EXECUTIVE_SUMMARY_PT_BR.md) · [Correction effectiveness](fix_effectiveness.md)

**The seam column-order defect is removed. Physical volumetric void fraction remains uncalibrated.** The analysis estimates visible gas under optical and geometric assumptions. Selection, formulas, thresholds, models and the user-confirmed geometry remain unchanged: **Python rows 18:71, D = 53 pixels**.

## Were the corrections effective?

The old assembly could turn `[1,2,3,4 | 5,6,7,8]` into `[1,2,3,6,5,6,7,8]`, producing three transitions when binarizing at intensity > 5 instead of one. v14.1 concatenates the same strips in their original order without that overwrite. All {n_runs} old panoramas were reproduced pixel for pixel as a control; corrected inputs preserve shifts and widths. The liquid reference was rebuilt with the same rule. Original files are preserved.

Across {n} eligible acquisitions, **maximum absolute α change: {delta}**; **relative flip-rate change: {flip_low}% to {flip_high}%**. This fixes the assembly defect, without proving better physical accuracy. The [per-run comparison](fix_comparison_per_run.csv) also freezes the old reference to isolate the direct input effect.

Of {effect['fixed_classifier_comparisons']} fixed model/validation pairs, balanced accuracy improved in **{better}**, worsened in **{worse}**, and stayed equal in **{same}**. Settings were not tuned to favor the plot or classifiers. Lower flip-rate is not automatically a better result.

## Void fraction versus spatial flip-rate

![Estimated void fraction versus spatial flip-rate](void_fraction_vs_flip_rate.png)

**α [−] is dimensionless:** inferred circular gas areas divided by circular tube area, then averaged spatially. Multiplying by 100 expresses the estimate as a percentage. **The Y axis is 1,000 F**, adjacent-pixel changes along the panorama, averaged over {n_rows} rows. This is **not a time frequency**, and is not per second.

The chart uses **{n} eligible acquisitions** in five patterns. Saturated acquisition m03q40_001 stays out of the chart and statistics but remains in the audit of {n_runs} selected records. Shaded regions are observed convex hulls, not classification boundaries. **Stars are arithmetic means with equal acquisition weights; the table contains medians.**

{table_text}

{overlap} This result applies to the primary **0.70** threshold. Dispersed bubbles have only **one acquisition**, m07q01_001, so within-class variation is unknown. Both axes derive from the same mask; their association does not validate physical void fraction.

## Threshold dependence

{threshold_text}

![Every existing threshold, using the same axis limits](threshold_sensitivity.png)

Hull separation changes with the threshold. Separation in the primary chart therefore is not robust across all tested settings. The original five thresholds and primary 0.70 setting were retained without selecting the most favorable outcome. Hulls describe available observations; they are neither classification tests nor confidence intervals.

## Steps, methods and equations

1. **Input, selection and saturation.** Preserve frames, shifts and center strips; remove only the defective seam overwrite. Original labels and the [previous selection](dataset_selection.json) remain. Mixed patterns, smooth annular and m05q05_002 stay excluded. The fraction S of intensity-255 pixels in rows 20:60 excludes acquisitions from statistics when **S > 0.20**. Reconstruction still assumes advection and depends on measured shifts.
2. **Normalization and liquid reference.** Divide each column by P75 of rows 30:50, bounded below by 1. In the same three plug panoramas, select columns at/below the fifth percentile of radial variability in rows 20:68, then apply medians and Gaussian σ = 1 pixel. These visually liquid-like regions are not certified pure-liquid calibration.
3. **Segmentation.** **R = Gσ(normalized I) / max(L, 0.01)**; **R < 0.70** detects shadows. Preserve upper/lower pairs, padded 3×3 closing, components without radial boundary contact and hole filling. Disconnected pairs need at least 20 consecutive columns; isolated components need at least 4 pixels. Segmentation implementation is unchanged.

![Panoramas with inferred gas masks in blue](panorama_phase_masks.png)

4. **Conditional circular area.** For separate radial intervals dᵢⱼ in column j: **αⱼ = Σᵢ(dᵢⱼ/D)²**; **ᾱ = (1/W) Σⱼ αⱼ**, D = 53. A full-lumen mask yields α = 1. A side view does not determine three-dimensional gas shape; circularity, refraction and unresolved films remain limitations.
5. **Spatial flip-rate.** **F = Σᵣ Σⱼ |Bᵣ,ⱼ₊₁ − Bᵣ,ⱼ| / [H(W−1)]**, H = 53, B = binary gas mask. The chart displays 1,000 F. Strip widths, seams and resolution still affect F. The [fixed y=40 probe](void_fraction_vs_center_flips.png) is a comparison, not a calibrated physical center.
6. **Sensitivity and adversarial controls.** Repeat all previous thresholds, references and internal crops. Independent controls check column order, masks, areas, denominators and written results. D is not changed to compensate for input effects. Keeping h interior rows imposes the mathematical ceiling (h/53)²; for m06q80_001 the crop results span α = {crop_low}–{crop_high}. This is not a confidence interval.

![Sensitivity to removal of known interior rows](void_crop_sensitivity.png)

7. **Classification and verification.** Preserve models, leave-one-run/group-out validation and references fitted only on training data. Rebuild confusion matrices, means, medians and plotted coordinates. Holding out the only bubble acquisition leaves no bubble training example; that test remains in the scores. Detector development used these data, so validation is internal. The [literature audit](source_audit.md) supports method discussion, not physical volume calibration of these images.

[Adversarial correction audit](seam_adversarial_review.md) · [Post-regeneration verification](post_regeneration_adversarial_review.md) · [Structured results](fix_effectiveness.json)
'''
