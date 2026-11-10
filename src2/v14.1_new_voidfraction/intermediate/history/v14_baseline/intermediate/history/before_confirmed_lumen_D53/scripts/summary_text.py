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
                  crop_values, changed_masks, max_alpha_change, overlap_pairs):
    n = sum(int(r['n']) for r in groups)
    ceiling = number((n_rows / diameter) ** 2, 4, pt)
    low, high = [number(v, 4, pt) for v in (min(crop_values), max(crop_values))]
    change = number(max_alpha_change, 6, pt)
    table = ['| Padrão | n | α estimado | Trocas / 1.000 pares de pixels |'] if pt else [
        '| Pattern | n | Estimated α | Changes / 1,000 pixel pairs |']
    table.append('|---|---:|---:|---:|')
    for group in groups:
        p = group['pattern']
        label = PATTERN_PT[p] if pt else p
        alpha = number(group['conditional_alpha_median'], 4, pt)
        flips = number(1000 * float(phase_groups[p]['candidate_pixel_flip_rate_median']), 3, pt)
        table.append(f"| {label} | {group['n']} | {alpha} | {flips} |")
    table_text = '\n'.join(table)
    pair_text = '; '.join((PATTERN_PT[a] + ' × ' + PATTERN_PT[b]) if pt else (a + ' × ' + b)
                          for a, b in overlap_pairs)
    overlap_text = ('Interseções entre as envoltórias observadas: **' + pair_text + '**.' if pair_text else
                    'Não há interseções entre as envoltórias calculadas nesta amostra.') if pt else (
                   'Observed hull intersections: **' + pair_text + '**.' if pair_text else
                   'The calculated hulls do not intersect in this sample.')
    if pt:
        opening = f'''# Resumo executivo — fração de vazio e flip-rate

[English version](EXECUTIVE_SUMMARY.md)

**Os panoramas permitem estimar o gás visível e contar alternâncias espaciais. A fração volumétrica real ainda não está calibrada.** Os valores dependem da segmentação e da hipótese de seções circulares.

## Fração de vazio × flip-rate espacial

![Fração de vazio estimada versus flip-rate espacial](void_fraction_vs_flip_rate.png)

**O que significa α [−]?** α é a fração de vazio estimada a partir das larguras de gás visíveis, assumindo geometria circular. **[−] significa adimensional**, sem unidade; multiplicar α por 100 expressa a estimativa em porcentagem. Isso não elimina as incertezas ópticas e geométricas.

**Eixo Y:** transições por 1.000 pixels na direção longitudinal do panorama, com média das {n_rows} linhas. Cada pixel é comparado com o da coluna seguinte; W colunas fornecem W−1 intervalos de um pixel. Assim, a normalização por 1.000 pares de pixels vizinhos equivale à escala indicada no eixo. Não é uma frequência em Hz.

O gráfico e a tabela usam **{n} aquisições elegíveis pelo critério de saturação**, em cinco padrões. A aquisição saturada m03q40_001 foi retirada do gráfico e das estatísticas; permanece na rastreabilidade dos {n_runs} registros selecionados. As áreas coloridas são envoltórias da dispersão observada, não fronteiras de classificação.

{table_text}

Valores: medianas por padrão. A classe de bolhas dispersas contém somente m07q01_001; uma aquisição não representa a variabilidade da classe. {overlap_text} A separação aparente de outros grupos nesta amostra não valida fronteiras para novas aquisições. Os dois eixos vêm da mesma máscara, portanto seu comportamento conjunto não valida a fração de vazio real.

## Etapas, métodos e equações
'''
        steps = [
            '1. **Seleção e controle de saturação.** Usamos os rótulos originais e a [política de seleção](dataset_selection.json). Classificações mistas, anular liso e m05q05_002 ficam excluídos. Método: fração S de pixels com intensidade 255 nas linhas 20:60; registros com **S > 0,20** não entram nas estatísticas. **Limite:** S ≤ 0,20 ainda permite alguns pixels saturados. É uma regra de qualidade herdada, não uma calibração.',
            '2. **Normalização por percentil e referência robusta de líquido.** Em três panoramas plug escolhidos visualmente, selecionamos colunas com desvio-padrão radial no percentil 5 ou abaixo. Pipeline: **intensidade / P75 das linhas 30:50 → seleção pela variação nas linhas 20:68 → medianas → filtro gaussiano (σ = 1 pixel)**. O divisor é limitado inferiormente a 1. **Limite:** os trechos parecem líquidos, mas não são uma referência certificada de líquido puro.',
            '3. **Segmentação por sombras e morfologia binária.** Calculamos **R = Gσ(I normalizada) / max(L, 0,01)**, com perfil de referência L. Pipeline: **R < 0,70 → pares de sombras superior/inferior → preenchimento entre elas**, somado a **fechamento 3×3 → componentes sem contato radial com o recorte → preenchimento de buracos**. Pares desconectados exigem intervalos de pelo menos 20 pixels; componentes pequenos exigem ao menos 4 pixels. **Correção:** o fechamento agora preserva os pixels das bordas antes de testar o contato. **Limite:** sombras não são contornos de gás comprovados; ruído e líquido dentro de envelopes podem ser classificados incorretamente.',
            '![Panoramas e máscaras de gás estimado em azul](panorama_phase_masks.png)',
            f'4. **Reconstrução circular e média espacial.** Cada intervalo de gás separado na coluna j tem largura dᵢⱼ. Equações: **αⱼ = Σᵢ(dᵢⱼ/D)²** e **ᾱ = (1/W) Σⱼ αⱼ**, com D = {diameter:g} pixels e W = número de colunas. É a razão entre áreas circulares πd²/4 e πD²/4. **Limite:** D é um diâmetro aparente; refração, formas não circulares e sobreposição em profundidade não estão resolvidas.',
            f'5. **Densidade espacial de transições binárias — flip-rate.** B = 1 indica gás estimado; B = 0 indica região não marcada. Equação: **F = Σᵣ Σⱼ |Bᵣ,ⱼ₊₁ − Bᵣ,ⱼ| / [H(W−1)]**, H = {n_rows}. O gráfico usa **1.000 F**. Um segmento de gás isolado e completo normalmente gera duas trocas. **Limite:** emendas, espaçamento em pixels e recorte afetam o resultado. A [comparação na linha fixa y=40](void_fraction_vs_center_flips.png) não representa um centro físico validado.',
            f'6. **Sensibilidade e controles adversariais.** Variamos limiares de 0,60 a 0,80, retiramos uma referência por vez e variamos recorte, diâmetro aparente e regra de intervalos. Controles sintéticos verificam bordas, áreas e contagens. A correção das bordas alterou {changed_masks} máscaras; a maior mudança absoluta de α foi {change}. **Limite:** são testes internos; o detector foi desenvolvido com estes mesmos dados.',
            '![Sensibilidade da fração de vazio ao recorte radial](void_crop_sensitivity.png)',
            f'**Por que aparece um teto de {ceiling}?** A máscara principal só ocupa {n_rows} linhas. Com D = {diameter:g}, a fórmula impõe **α ≤ (53/60)² ≈ {ceiling}**. Mesmo com intervalos separados, a soma das larguras não ultrapassa 53. **Esse teto é do cálculo, não do escoamento.** Em m06q80_001, mudar somente o recorte produz α entre {low} e {high}. Essa faixa é sensibilidade, não intervalo de confiança. Recortes maiores podem incluir reflexos da parede.',
            '7. **Verificação independente da aritmética e confronto com a literatura.** Recalculamos áreas pelo início/fim dos intervalos, trocas, denominadores, medianas e coordenadas dos gráficos. Fontes locais: Schepperle et al. (2024), §2.3, para segmentação e referência líquida; de Oliveira et al. (2021) e Dias et al. (2025), para distinguir suas estimativas por correlações das nossas estimativas ópticas. **Limite:** reproduzir contas e citar métodos não comprova o volume real de gás destes registros.',
        ]
        ending = '''## Conclusões e próximos passos

O anular ondulado tende a apresentar uma região de gás ampla e poucas trocas na máscara. A proximidade de alguns pontos ao teto do recorte impede interpretar o agrupamento como alta precisão. A separação automática de todos os padrões continua sem validação externa.

- Anotar contornos em trechos representativos de todos os padrões e medir erros da máscara.
- Verificar bordas internas do tubo e distorção óptica antes de afirmar frações volumétricas absolutas.

O intervalo de extensão nas tabelas cobre somente um cenário para gás já detectado tocando as bordas. Ele exige ausência de gás totalmente escondido fora do recorte; não é um limite geral da fração real.

Detalhes: [revisão dos métodos](deep_method_review.md), [gráficos para a discussão do TCC](GRAFICOS_PARA_DISCUSSAO_PT_BR.md), [valores por aquisição](panorama_void_per_run.csv), [contagens de flip-rate](phase_mask_flip_per_run.csv), [fontes](source_audit.md), [auditoria após o resumo](post_regeneration_adversarial_review.md).
'''
    else:
        opening = f'''# Executive summary — void fraction and flip-rate

[Versão em português](EXECUTIVE_SUMMARY_PT_BR.md)

**The panoramas support estimates of visible gas and spatial phase alternation. Absolute gas-volume fractions remain uncalibrated.** Values depend on segmentation and assumed circular cross-sections.

## Void fraction × spatial flip-rate

![Estimated void fraction versus spatial flip-rate](void_fraction_vs_flip_rate.png)

**What does α [−] mean?** α is the void estimate reconstructed from visible gas widths under circular geometry. **[−] means dimensionless**, with no unit; multiplying α by 100 expresses the estimate as a percentage. Optical and geometric uncertainty still applies.

**Y-axis:** transitions per 1,000 pixels along the panorama, averaged over {n_rows} rows. Compare each pixel with the next column in the same row; W columns provide W−1 one-pixel intervals. Scaling the count per 1,000 adjacent pixel pairs therefore gives the stated spatial rate, not per second.

The plot and table use **{n} acquisitions passing the saturation criterion** in five patterns. Saturated m03q40_001 is omitted from the plot and statistics but retained in the provenance of the {n_runs} selected records. Shaded areas are observed convex hulls, not classification boundaries.

{table_text}

Values are group medians. Only m07q01_001 represents dispersed bubbles; one acquisition cannot establish within-pattern variability. {overlap_text} Apparent separation of other groups in this sample does not validate boundaries for new acquisitions. Both axes use the same mask, so their joint behavior does not validate true void fraction.

## Steps, methods and equations
'''
        steps = [
            '1. **Selection and saturation control.** Original labels and the [selection policy](dataset_selection.json) exclude mixed classifications, smooth annular and m05q05_002. Method: fraction S of intensity-255 pixels in rows 20:60; **S > 0.20** excludes an acquisition from statistics. **Limit:** S ≤ 0.20 still allows some saturated pixels. This is an inherited quality rule, not a calibration.',
            '2. **Percentile normalization and robust liquid reference.** In three visually selected plug panoramas, choose columns at or below the 5th percentile of radial standard deviation. Pipeline: **intensity / P75 in rows 30:50 → select by variation in rows 20:68 → medians → Gaussian filter (σ = 1 pixel)**. The divisor is floored at 1. **Limit:** liquid-looking gaps are not certified liquid-only references.',
            '3. **Shadow segmentation and binary morphology.** Compute **R = Gσ(normalized image) / max(L, 0.01)**, with reference profile L. Pipeline: **R < 0.70 → paired upper/lower shadows → fill between them**, plus **3×3 closing → components without radial crop contact → hole filling**. Disconnected pairs require at least 20 pixels of accepted interval; small components require at least 4 pixels. **Correction:** closing now preserves border pixels before testing crop contact. **Limit:** shadows are not proven gas outlines; noise and liquid inside envelopes can be mislabeled.',
            '![Original panoramas and inferred gas masks in blue](panorama_phase_masks.png)',
            f'4. **Circular reconstruction and spatial averaging.** Each separate gas interval in column j has width dᵢⱼ. Equations: **αⱼ = Σᵢ(dᵢⱼ/D)²** and **ᾱ = (1/W) Σⱼ αⱼ**, with apparent D = {diameter:g} pixels and panorama width W. This divides circular areas πd²/4 by πD²/4. **Limit:** tube edges are unvalidated; refraction, noncircular shapes and depth overlap remain unresolved.',
            f'5. **Binary spatial transition density — flip-rate.** B = 1 denotes inferred gas and B = 0 unmarked pixels. Equation: **F = Σᵣ Σⱼ |Bᵣ,ⱼ₊₁ − Bᵣ,ⱼ| / [H(W−1)]**, H = {n_rows}. The chart uses **1,000 F**. A complete isolated gas segment normally contributes two changes. **Limit:** seams, pixel spacing and crop affect this measure. The [fixed-row y=40 comparison](void_fraction_vs_center_flips.png) is not a validated physical tube center.',
            f'6. **Sensitivity analysis and adversarial controls.** Vary thresholds from 0.60 to 0.80, leave one reference out, and perturb crop, apparent diameter and interval rules. Synthetic controls test borders, areas and counts. Border correction changed {changed_masks} masks; the largest absolute α change was {change}. **Limit:** these are internal checks on data used to develop the detector.',
            '![Void sensitivity to radial crop](void_crop_sensitivity.png)',
            f'**Why is there a ceiling of {ceiling}?** The main mask occupies only {n_rows} rows. With D = {diameter:g}, the formula imposes **α ≤ (53/60)² ≈ {ceiling}**. Even separate intervals cannot have total width above 53. **This is a processing ceiling, not a flow limit.** For m06q80_001, changing only the crop produces α from {low} to {high}. This is sensitivity, not a confidence interval. Wider crops can include wall reflections.',
            '7. **Independent arithmetic verification and literature comparison.** Recompute areas from interval starts/ends, flip counts, denominators, medians and chart coordinates. Local sources: Schepperle et al. (2024), §2.3, for segmentation and liquid references; de Oliveira et al. (2021) and Dias et al. (2025), to distinguish their correlation estimates from these optical estimates. **Limit:** reproducible arithmetic and source methods do not validate the true gas volume of these records.',
        ]
        ending = '''## Findings and next steps

Wavy annular images tend to show broad gas regions and few mask changes. Clustering near the crop ceiling cannot establish precision. Automatic separation of all patterns still lacks external validation.

- Annotate representative outlines across all patterns and measure mask errors.
- Validate internal tube boundaries and optical distortion before claiming absolute gas-volume fractions.

The tabulated extension interval describes only a scenario for already-detected gas touching crop edges. It requires no wholly unseen gas outside the crop; it is not a general physical bound.

Details: [method review](deep_method_review.md), [thesis discussion graphics](GRAFICOS_PARA_DISCUSSAO_PT_BR.md), [per-run void](panorama_void_per_run.csv), [flip counts](phase_mask_flip_per_run.csv), [literature](source_audit.md), [post-summary audit](post_regeneration_adversarial_review.md).
'''
    return opening + '\n' + '\n\n'.join(steps) + '\n\n' + ending
