# Resumo executivo — fração de vazio e flip-rate

[English version](EXECUTIVE_SUMMARY.md)

**Os panoramas permitem estimar o gás visível e contar alternâncias espaciais. A fração volumétrica real ainda não está calibrada.** Os valores dependem da segmentação e da hipótese de seções circulares.

## Fração de vazio × flip-rate espacial

![Fração de vazio estimada versus flip-rate espacial](../figures/void_fraction/void_fraction_vs_flip_rate.png)

**O que significa α [−]?** α é a fração de vazio estimada a partir das larguras de gás visíveis, assumindo geometria circular. **[−] significa adimensional**, sem unidade; multiplicar α por 100 expressa a estimativa em porcentagem. Isso não elimina as incertezas ópticas e geométricas.

**Eixo Y:** transições por 1.000 pixels na direção longitudinal do panorama, com média das 53 linhas. Cada pixel é comparado com o da coluna seguinte; W colunas fornecem W−1 intervalos de um pixel. Assim, a normalização por 1.000 pares de pixels vizinhos equivale à escala indicada no eixo. Não é uma frequência em Hz.

O gráfico e a tabela usam **25 aquisições elegíveis pelo critério de saturação**, em cinco padrões. A aquisição saturada m03q40_001 foi retirada do gráfico e das estatísticas; permanece na rastreabilidade dos 26 registros selecionados. As áreas coloridas são envoltórias da dispersão observada, não fronteiras de classificação.

| Padrão | n | α estimado | Trocas / 1.000 pares de pixels |
|---|---:|---:|---:|
| Bolhas dispersas | 1 | 0,0014 | 2,464 |
| Plug | 3 | 0,4628 | 9,730 |
| Slug | 3 | 0,4731 | 10,419 |
| Churn | 10 | 0,6419 | 11,190 |
| Anular ondulado | 8 | 0,7673 | 1,133 |

Valores: medianas por padrão. A classe de bolhas dispersas contém somente m07q01_001; uma aquisição não representa a variabilidade da classe. Interseções entre as envoltórias observadas: **Churn × Anular ondulado**. A separação aparente de outros grupos nesta amostra não valida fronteiras para novas aquisições. Os dois eixos vêm da mesma máscara, portanto seu comportamento conjunto não valida a fração de vazio real.

## Etapas, métodos e equações

1. **Seleção e controle de saturação.** Usamos os rótulos originais e a [política de seleção](../results/validation/dataset_selection.json). Classificações mistas, anular liso e m05q05_002 ficam excluídos. Método: fração S de pixels com intensidade 255 nas linhas 20:60; registros com **S > 0,20** não entram nas estatísticas. **Limite:** S ≤ 0,20 ainda permite alguns pixels saturados. É uma regra de qualidade herdada, não uma calibração.

2. **Normalização por percentil e referência robusta de líquido.** Em três panoramas plug escolhidos visualmente, selecionamos colunas com desvio-padrão radial no percentil 5 ou abaixo. Pipeline: **intensidade / P75 das linhas 30:50 → seleção pela variação nas linhas 20:68 → medianas → filtro gaussiano (σ = 1 pixel)**. O divisor é limitado inferiormente a 1. **Limite:** os trechos parecem líquidos, mas não são uma referência certificada de líquido puro.

3. **Segmentação por sombras e morfologia binária.** Calculamos **R = Gσ(I normalizada) / max(L, 0,01)**, com perfil de referência L. Pipeline: **R < 0,70 → pares de sombras superior/inferior → preenchimento entre elas**, somado a **fechamento 3×3 → componentes sem contato radial com o recorte → preenchimento de buracos**. Pares desconectados exigem intervalos de pelo menos 20 pixels; componentes pequenos exigem ao menos 4 pixels. **Correção:** o fechamento agora preserva os pixels das bordas antes de testar o contato. **Limite:** sombras não são contornos de gás comprovados; ruído e líquido dentro de envelopes podem ser classificados incorretamente.

![Panoramas e máscaras de gás estimado em azul](../figures/void_fraction/panorama_phase_masks.png)

4. **Reconstrução circular e média espacial.** Cada intervalo de gás separado na coluna j tem largura dᵢⱼ. Equações: **αⱼ = Σᵢ(dᵢⱼ/D)²** e **ᾱ = (1/W) Σⱼ αⱼ**, com D = 60 pixels e W = número de colunas. É a razão entre áreas circulares πd²/4 e πD²/4. **Limite:** D é um diâmetro aparente; refração, formas não circulares e sobreposição em profundidade não estão resolvidas.

5. **Densidade espacial de transições binárias — flip-rate.** B = 1 indica gás estimado; B = 0 indica região não marcada. Equação: **F = Σᵣ Σⱼ |Bᵣ,ⱼ₊₁ − Bᵣ,ⱼ| / [H(W−1)]**, H = 53. O gráfico usa **1.000 F**. Um segmento de gás isolado e completo normalmente gera duas trocas. **Limite:** emendas, espaçamento em pixels e recorte afetam o resultado. A [comparação na linha fixa y=40](../figures/void_fraction/void_fraction_vs_center_flips.png) não representa um centro físico validado.

6. **Sensibilidade e controles adversariais.** Variamos limiares de 0,60 a 0,80, retiramos uma referência por vez e variamos recorte, diâmetro aparente e regra de intervalos. Controles sintéticos verificam bordas, áreas e contagens. A correção das bordas alterou 24 máscaras; a maior mudança absoluta de α foi 0,005417. **Limite:** são testes internos; o detector foi desenvolvido com estes mesmos dados.

![Sensibilidade da fração de vazio ao recorte radial](../figures/void_fraction/void_crop_sensitivity.png)

**Por que aparece um teto de 0,7803?** A máscara principal só ocupa 53 linhas. Com D = 60, a fórmula impõe **α ≤ (53/60)² ≈ 0,7803**. Mesmo com intervalos separados, a soma das larguras não ultrapassa 53. **Esse teto é do cálculo, não do escoamento.** Em m06q80_001, mudar somente o recorte produz α entre 0,6826 e 0,8575. Essa faixa é sensibilidade, não intervalo de confiança. Recortes maiores podem incluir reflexos da parede.

7. **Verificação independente da aritmética e confronto com a literatura.** Recalculamos áreas pelo início/fim dos intervalos, trocas, denominadores, medianas e coordenadas dos gráficos. Fontes locais: Schepperle et al. (2024), §2.3, para segmentação e referência líquida; de Oliveira et al. (2021) e Dias et al. (2025), para distinguir suas estimativas por correlações das nossas estimativas ópticas. **Limite:** reproduzir contas e citar métodos não comprova o volume real de gás destes registros.

## Conclusões e próximos passos

O anular ondulado tende a apresentar uma região de gás ampla e poucas trocas na máscara. A proximidade de alguns pontos ao teto do recorte impede interpretar o agrupamento como alta precisão. A separação automática de todos os padrões continua sem validação externa.

- Anotar contornos em trechos representativos de todos os padrões e medir erros da máscara.
- Verificar bordas internas do tubo e distorção óptica antes de afirmar frações volumétricas absolutas.

O intervalo de extensão nas tabelas cobre somente um cenário para gás já detectado tocando as bordas. Ele exige ausência de gás totalmente escondido fora do recorte; não é um limite geral da fração real.

Detalhes: [revisão dos métodos](deep_method_review.md), [gráficos para a discussão do TCC](GRAFICOS_PARA_DISCUSSAO_PT_BR.md), [valores por aquisição](../results/void_fraction/panorama_void_per_run.csv), [contagens de flip-rate](../results/flip_rate/phase_mask_flip_per_run.csv), [fontes](source_audit.md), [auditoria após o resumo](post_regeneration_adversarial_review.md).
