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
                  crop_values, historical_crop_values, historical_crop_ceiling,
                  changed_masks, max_alpha_change, overlap_pairs):
    n = sum(int(r['n']) for r in groups)
    assert n_rows == diameter == 53, 'The confirmed full lumen must contain 53 rows'
    low, high = [number(v, 4, pt) for v in (min(crop_values), max(crop_values))]
    historical_low, historical_high = [number(v, 4, pt) for v in
                                      (min(historical_crop_values), max(historical_crop_values))]
    historical_ceiling = number(historical_crop_ceiling, 4, pt)
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

**Correção da escala:** o usuário confirmou que o interior completo ocupa as linhas Python **18:71**, ou seja, 53 pixels. Usar D = 60 era incompatível com esse limite. Agora **D = 53**, e **α₅₃ = α₆₀ × (60/53)²** para a mesma máscara. As máscaras e o flip-rate não mudam com essa correção.

## Fração de vazio × flip-rate espacial

![Fração de vazio estimada versus flip-rate espacial](void_fraction_vs_flip_rate.png)

**O que significa α [−]?** α é a fração de vazio estimada a partir das larguras de gás visíveis, assumindo geometria circular. **[−] significa adimensional**, sem unidade; multiplicar α por 100 expressa a estimativa em porcentagem. Isso não elimina as incertezas ópticas e geométricas.

**Eixo Y:** transições por 1.000 pixels na direção longitudinal do panorama, com média das {n_rows} linhas. Cada pixel é comparado com o da coluna seguinte; W colunas fornecem W−1 intervalos de um pixel. Assim, a normalização por 1.000 pares de pixels vizinhos equivale à escala indicada no eixo. Não é uma frequência em Hz.

O gráfico e a tabela usam **{n} aquisições elegíveis pelo critério de saturação**, em cinco padrões. A aquisição saturada m03q40_001 foi retirada do gráfico e das estatísticas; permanece na rastreabilidade dos {n_runs} registros selecionados. As áreas coloridas são envoltórias da dispersão observada, não fronteiras de classificação. **As estrelas indicam o centroide de cada padrão: média de α e do flip-rate, com peso igual por aquisição.** A tabela abaixo apresenta medianas. Para bolhas dispersas, a estrela coincide com a única aquisição.

{table_text}

Valores: medianas por padrão. A classe de bolhas dispersas contém somente m07q01_001; uma aquisição não representa a variabilidade da classe. {overlap_text} A separação aparente de outros grupos nesta amostra não valida fronteiras para novas aquisições. Os dois eixos vêm da mesma máscara, portanto seu comportamento conjunto não valida a fração de vazio real.

## Etapas, métodos e equações
'''
        steps = [
            '1. **Seleção e controle de saturação.** Usamos os rótulos originais e a [política de seleção](dataset_selection.json). Classificações mistas, anular liso e m05q05_002 ficam excluídos. Método: fração S de pixels com intensidade 255 nas linhas 20:60; registros com **S > 0,20** não entram nas estatísticas. **Limite:** S ≤ 0,20 ainda permite alguns pixels saturados. É uma regra de qualidade herdada, não uma calibração.',
            '2. **Normalização por percentil e referência robusta de líquido.** Em três panoramas plug escolhidos visualmente, selecionamos colunas com desvio-padrão radial no percentil 5 ou abaixo. Pipeline: **intensidade / P75 das linhas 30:50 → seleção pela variação nas linhas 20:68 → medianas → filtro gaussiano (σ = 1 pixel)**. O divisor é limitado inferiormente a 1. **Limite:** os trechos parecem líquidos, mas não são uma referência certificada de líquido puro.',
            '3. **Segmentação por sombras e morfologia binária.** Calculamos **R = Gσ(I normalizada) / max(L, 0,01)**, com perfil de referência L. Pipeline: **R < 0,70 → pares de sombras superior/inferior → preenchimento entre elas**, somado a **fechamento 3×3 → componentes sem contato radial com o limite interno → preenchimento de buracos**. Pares desconectados exigem pelo menos 20 colunas consecutivas aceitas; componentes pequenos exigem ao menos 4 pixels. **Correção:** o fechamento preserva os pixels das bordas antes de testar o contato. **Limite:** identificar fase pela imagem exige validação.',
            '![Panoramas e máscaras de gás estimado em azul](panorama_phase_masks.png)',
            f'4. **Reconstrução circular e média espacial.** Cada intervalo de gás separado na coluna j tem largura dᵢⱼ. Equações: **αⱼ = Σᵢ(dᵢⱼ/D)²** e **ᾱ = (1/W) Σⱼ αⱼ**, com D = {diameter:g} pixels e W = número de colunas. É a razão entre áreas circulares πd²/4 e πD²/4. Uma máscara que preenche todo o interior produz **α = 1**. **Limite:** preencher a máscara não comprova ausência de líquido; a geometria continua sendo uma hipótese.',
            f'5. **Densidade espacial de transições binárias — flip-rate.** B = 1 indica gás estimado; B = 0 indica região não marcada. Equação: **F = Σᵣ Σⱼ |Bᵣ,ⱼ₊₁ − Bᵣ,ⱼ| / [H(W−1)]**, H = {n_rows}. O gráfico usa **1.000 F**. Um segmento de gás isolado e completo normalmente gera duas trocas. **Limite:** emendas, espaçamento em pixels e recorte afetam o resultado. A [comparação na linha fixa y=40](void_fraction_vs_center_flips.png) não representa um centro físico validado.',
            f'6. **Sensibilidade e controles adversariais.** Variamos limiares de 0,60 a 0,80, retiramos uma referência por vez e variamos recortes internos e regras de intervalos. Controles sintéticos verificam bordas, áreas e contagens. A correção anterior da morfologia alterou {changed_masks} máscaras; a maior mudança absoluta de α foi {change}, **calculada historicamente com D = 60**. Essa mudança é separada da correção atual para D = 53. **Limite:** são testes internos; o detector foi desenvolvido com estes mesmos dados.',
            '![Sensibilidade da fração de vazio ao recorte radial](void_crop_sensitivity.png)',
            f'**Leitura do teste de recorte:** mantemos 50, 51, 52 ou 53 linhas do interior, sempre com **D = 53**. Para h linhas mantidas, o máximo matemático é **(h/53)²**; o interior completo permite **α = 1**. Em m06q80_001, esses recortes produzem α entre {low} e {high}. **Limite:** os recortes menores removem linhas conhecidas do interior; não representam outros diâmetros nem um intervalo de confiança. Regiões externas ao interior confirmado não entram no cálculo.',
            '7. **Verificação independente da aritmética e confronto com a literatura.** Recalculamos áreas pelo início/fim dos intervalos, trocas, denominadores, medianas e coordenadas dos gráficos. Fontes locais: Schepperle et al. (2024), §2.3, para segmentação e referência líquida; de Oliveira et al. (2021) e Dias et al. (2025), para distinguir suas estimativas por correlações das nossas estimativas ópticas. **Limite:** reproduzir contas e citar métodos não comprova o volume real de gás destes registros.',
        ]
        limitations = f'''## Limitações científicas para a discussão

- **Identificação do gás:** a câmera registra luz e sombras. Preencher contornos pode contar reflexos ou líquido como gás. A conta correta não comprova que cada pixel recebeu a fase correta.
- **Máscara binária e resolução:** os tons de cinza ajudam a detectar contornos, mas a saída é **1 = gás estimado; 0 = não marcado como gás**. Não estimamos ocupação parcial de cada pixel nem atribuímos uma probabilidade à sua classificação. Filmes líquidos e bolhas abaixo da resolução podem não ser resolvidos, afetando α e o flip-rate. **Não temos calibração que converta intensidade de cinza em proporção de gás/líquido.**
- **Reconstrução da área:** para um intervalo circular, **α = (d/D)²**; para vários, somam-se as áreas estimadas. Um tubo circular não garante gás circular, sobretudo em churn. É como olhar um objeto de lado: a largura visível não revela necessariamente a profundidade. A parede curva e os fluidos também desviam a luz por refração e podem deslocar as bordas aparentes.
- **Histórico do recorte, já corrigido:** usar 53 linhas com D = 60 impunha **(53/60)² ≈ {historical_ceiling}**. No [experimento antigo de m06q80_001](intermediate/history/before_confirmed_lumen_D53/results/validation/adversarial_checks.json), recortes de 50/53/57/60 linhas, todos com D = 60, produziram **{historical_low}–{historical_high}**. Os recortes mais largos ultrapassavam o interior agora confirmado, **18:71, D = 53**. Essa faixa antiga não representa a incerteza atual nem indica gás ausente fora do tubo. O gráfico acima testa somente remoção de linhas internas. A [auditoria da normalização](lumen_normalization_review.md) separa as duas versões.

**Os resultados são reproduzíveis sob essas hipóteses; sua correspondência com o volume real de gás ainda precisa de validação.**
'''
        ending = '''## Conclusões e próximos passos

O anular ondulado tende a apresentar uma região de gás ampla e poucas trocas na máscara. A separação automática de todos os padrões continua sem validação externa.

- Anotar contornos em trechos representativos de todos os padrões e medir erros da máscara.
- Calibrar a relação entre a silhueta e o volume real, considerando refração e formas não circulares, antes de afirmar frações volumétricas absolutas.

Detalhes: [revisão dos métodos](deep_method_review.md), [gráficos para a discussão do TCC](GRAFICOS_PARA_DISCUSSAO_PT_BR.md), [valores por aquisição](panorama_void_per_run.csv), [contagens de flip-rate](phase_mask_flip_per_run.csv), [fontes](source_audit.md), [auditoria após o resumo](post_regeneration_adversarial_review.md).
'''
    else:
        opening = f'''# Executive summary — void fraction and flip-rate

[Versão em português](EXECUTIVE_SUMMARY_PT_BR.md)

**The panoramas support estimates of visible gas and spatial phase alternation. Absolute gas-volume fractions remain uncalibrated.** Values depend on segmentation and assumed circular cross-sections.

**Scale correction:** the user confirmed that the complete tube interior occupies Python rows **18:71**, or 53 pixels. D = 60 was inconsistent with that boundary. We now use **D = 53**, and **α₅₃ = α₆₀ × (60/53)²** for the same mask. This correction leaves the masks and flip-rate unchanged.

## Void fraction × spatial flip-rate

![Estimated void fraction versus spatial flip-rate](void_fraction_vs_flip_rate.png)

**What does α [−] mean?** α is the void estimate reconstructed from visible gas widths under circular geometry. **[−] means dimensionless**, with no unit; multiplying α by 100 expresses the estimate as a percentage. Optical and geometric uncertainty still applies.

**Y-axis:** transitions per 1,000 pixels along the panorama, averaged over {n_rows} rows. Compare each pixel with the next column in the same row; W columns provide W−1 one-pixel intervals. Scaling the count per 1,000 adjacent pixel pairs therefore gives the stated spatial rate, not per second.

The plot and table use **{n} acquisitions passing the saturation criterion** in five patterns. Saturated m03q40_001 is omitted from the plot and statistics but retained in the provenance of the {n_runs} selected records. Shaded areas are observed convex hulls, not classification boundaries. **Stars mark each pattern centroid: mean α and mean flip-rate, with equal weight per acquisition.** The table below reports medians. The dispersed-bubble star coincides with its sole acquisition.

{table_text}

Values are group medians. Only m07q01_001 represents dispersed bubbles; one acquisition cannot establish within-pattern variability. {overlap_text} Apparent separation of other groups in this sample does not validate boundaries for new acquisitions. Both axes use the same mask, so their joint behavior does not validate true void fraction.

## Steps, methods and equations
'''
        steps = [
            '1. **Selection and saturation control.** Original labels and the [selection policy](dataset_selection.json) exclude mixed classifications, smooth annular and m05q05_002. Method: fraction S of intensity-255 pixels in rows 20:60; **S > 0.20** excludes an acquisition from statistics. **Limit:** S ≤ 0.20 still allows some saturated pixels. This is an inherited quality rule, not a calibration.',
            '2. **Percentile normalization and robust liquid reference.** In three visually selected plug panoramas, choose columns at or below the 5th percentile of radial standard deviation. Pipeline: **intensity / P75 in rows 30:50 → select by variation in rows 20:68 → medians → Gaussian filter (σ = 1 pixel)**. The divisor is floored at 1. **Limit:** liquid-looking gaps are not certified liquid-only references.',
            '3. **Shadow segmentation and binary morphology.** Compute **R = Gσ(normalized image) / max(L, 0.01)**, with reference profile L. Pipeline: **R < 0.70 → paired upper/lower shadows → fill between them**, plus **3×3 closing → components without radial lumen-edge contact → hole filling**. Disconnected pairs require at least 20 consecutive accepted columns; small components require at least 4 pixels. **Correction:** closing preserves border pixels before testing contact. **Limit:** image-based phase identification requires validation.',
            '![Original panoramas and inferred gas masks in blue](panorama_phase_masks.png)',
            f'4. **Circular reconstruction and spatial averaging.** Each separate gas interval in column j has width dᵢⱼ. Equations: **αⱼ = Σᵢ(dᵢⱼ/D)²** and **ᾱ = (1/W) Σⱼ αⱼ**, with D = {diameter:g} pixels and panorama width W. This divides circular areas πd²/4 by πD²/4. A mask filling the complete lumen gives **α = 1**. **Limit:** a filled mask does not prove that no liquid is present; geometry remains an assumption.',
            f'5. **Binary spatial transition density — flip-rate.** B = 1 denotes inferred gas and B = 0 unmarked pixels. Equation: **F = Σᵣ Σⱼ |Bᵣ,ⱼ₊₁ − Bᵣ,ⱼ| / [H(W−1)]**, H = {n_rows}. The chart uses **1,000 F**. A complete isolated gas segment normally contributes two changes. **Limit:** seams, pixel spacing and crop affect this measure. The [fixed-row y=40 comparison](void_fraction_vs_center_flips.png) is not a validated physical tube center.',
            f'6. **Sensitivity analysis and adversarial controls.** Vary thresholds from 0.60 to 0.80, leave one reference out, and perturb internal crops and interval rules. Synthetic controls test borders, areas and counts. The earlier morphology correction changed {changed_masks} masks; the largest absolute α change was {change}, **historically calculated with D = 60**. This change is separate from the current correction to D = 53. **Limit:** these are internal checks on data used to develop the detector.',
            '![Void sensitivity to radial crop](void_crop_sensitivity.png)',
            f'**Reading the crop test:** retain 50, 51, 52 or 53 rows of the interior, always with **D = 53**. For h retained rows, the mathematical maximum is **(h/53)²**; the full lumen permits **α = 1**. These crops produce α from {low} to {high} for m06q80_001. **Limit:** smaller crops discard known-interior rows; they are not alternative diameters or a confidence interval. Regions outside the confirmed lumen are excluded from the calculation.',
            '7. **Independent arithmetic verification and literature comparison.** Recompute areas from interval starts/ends, flip counts, denominators, medians and chart coordinates. Local sources: Schepperle et al. (2024), §2.3, for segmentation and liquid references; de Oliveira et al. (2021) and Dias et al. (2025), to distinguish their correlation estimates from these optical estimates. **Limit:** reproducible arithmetic and source methods do not validate the true gas volume of these records.',
        ]
        limitations = f'''## Scientific limitations for the discussion

- **Gas identification:** the camera records light and shadows. Filling contours may count reflections or liquid as gas. Correct arithmetic does not prove that each pixel received the correct phase label.
- **Binary mask and resolution:** grayscale helps detect contours, but the output is **1 = inferred gas; 0 = unmarked as gas**. We do not estimate partial occupancy within each pixel or assign a probability to its classification. Liquid films and bubbles below image resolution may remain unresolved, affecting α and flip-rate. **We have no calibration converting grayscale intensity into a gas/liquid proportion.**
- **Area reconstruction:** for one circular interval, **α = (d/D)²**; for several, sum the estimated areas. A circular tube does not ensure circular gas, especially in churn. It is like looking at an object from the side: visible width does not necessarily reveal depth. The curved wall and fluids also refract light and may shift apparent edges.
- **Crop history, now corrected:** using 53 rows with D = 60 imposed **(53/60)² ≈ {historical_ceiling}**. The [old m06q80_001 experiment](intermediate/history/before_confirmed_lumen_D53/results/validation/adversarial_checks.json) used 50/53/57/60-row crops, all with D = 60, and produced **{historical_low}–{historical_high}**. Its wider crops extended beyond the now-confirmed lumen, **18:71, D = 53**. This old range is not current uncertainty or evidence of gas missing outside the tube. The chart above tests only removal of interior rows. The [normalization audit](lumen_normalization_review.md) separates the two versions.

**Results are reproducible under these assumptions; their correspondence with true gas volume still needs validation.**
'''
        ending = '''## Findings and next steps

Wavy annular images tend to show broad gas regions and few mask changes. Automatic separation of all patterns still lacks external validation.

- Annotate representative outlines across all patterns and measure mask errors.
- Calibrate the relation between silhouettes and true volume, including refraction and noncircular shapes, before claiming absolute gas-volume fractions.

Details: [method review](deep_method_review.md), [thesis discussion graphics](GRAFICOS_PARA_DISCUSSAO_PT_BR.md), [per-run void](panorama_void_per_run.csv), [flip counts](phase_mask_flip_per_run.csv), [literature](source_audit.md), [post-summary audit](post_regeneration_adversarial_review.md).
'''
    return opening + '\n' + '\n\n'.join(steps) + '\n\n' + limitations + '\n' + ending
