# Gráficos para a discussão do TCC

**Use três figuras no corpo principal:** fração de vazio × flip-rate, sensibilidade ao recorte e panorama com máscara. Juntas mostram o resultado, a dependência do processamento e o que foi marcado como gás.

Foram revisados os 17 PNG atuais e os três SVG correspondentes. Os arquivos históricos ficam fora dos resultados atuais. O gráfico principal conserva seu caminho original e usa 25 aquisições elegíveis pelo critério de saturação.

**Na legenda da tese:** α é adimensional e condicionado à reconstrução circular, com D=53 pixels no interior confirmado 18:71. O flip-rate conta mudanças entre pixels vizinhos ao longo das colunas, com média de 53 linhas. As áreas coloridas representam dispersão da amostra. A aquisição saturada está omitida. A explicação está no [resumo em PT-BR](EXECUTIVE_SUMMARY_PT_BR.md).

## Figuras existentes

| Figura | Destino | Contribuição / limite |
|---|---|---|
| [Fração de vazio × flip-rate](../figures/void_fraction/void_fraction_vs_flip_rate.png) | Discussão principal | Mostra a posição e a dispersão dos padrões. Os dois eixos vêm da mesma máscara. |
| [Sensibilidade ao recorte](../figures/void_fraction/void_crop_sensitivity.png) | Discussão principal | Mostra o efeito de retirar linhas do interior confirmado, mantendo D=53. Faixas de sensibilidade não são intervalos de confiança. |
| [Panorama e máscara](../figures/void_fraction/panorama_phase_masks.png) | Método e discussão | Mostra o gás estimado em azul. Exemplos não validam todas as bordas. |
| [Linha fixa y=40](../figures/void_fraction/void_fraction_vs_center_flips.png) | Apêndice | Compara uma linha com a média de 53. O ponto médio das linhas internas é y=44; y=40 foi preservado como sonda fixa. |
| [Referência líquida](../figures/void_fraction/panorama_liquid_reference_candidates.png) | Método ou apêndice | Documenta a origem da referência inferida. Não é fotografia certificada de líquido puro. |
| [Conjunto de panoramas](../figures/patterns/pattern_panorama_contact_sheet.png) | Apêndice | Audita rótulos e a imagem saturada. Exibe recortes centrais de 2.400 pixels. |
| [Matrizes de morfologia](../figures/patterns/pattern_confusion_matrices.png) | Apêndice, se discutir classificação | Validação interna; a única aquisição de bolhas fica sem exemplos da classe em seu treino. |
| [Matrizes de silhuetas](../figures/patterns/silhouette_pattern_confusion_matrices.png) | Apêndice, se discutir classificação | Mesma limitação; concordância com rótulos não valida volume de gás. |
| [Trechos de inspeção de fase](../figures/validation/adversarial_phase_anchor_patches.png) | Apêndice | Interiores selecionados manualmente; não medem erro global da máscara. |
| [Intervalo líquido de slug](../figures/validation/adversarial_slug_gap_mask.png) | Apêndice | Exemplo de ponte espúria. Esquerda: original; direita: máscara. Identificar intervalo e aquisição. |
| [Flip óptico espacial e temporal](../figures/flip_rate/flip_rate_spatial_vs_native.png) | Apêndice | Compara bases de amostragem. Usa brilho, não a máscara principal de gás; barras pretas são medianas. |
| [Quimogramas nativos](../figures/flip_rate/flip_rate_native_kymographs.png) | Apêndice | Diagnóstico óptico temporal curto. Picos não equivalem automaticamente à passagem de bolhas. |
| [Quadro bruto × panorama](../intermediate/exploratory/figures/raw_comparison.png) | Apêndice metodológico | Mostra iluminação e reconstrução. As colunas têm extensões espaciais diferentes. |
| [Largura clara: média × mediana](../intermediate/exploratory/figures/bright_width_mean_vs_median.png) | Exploratório | Diagnostica o descritor anterior; não representa α. |
| [Detector antigo de interfaces](../intermediate/exploratory/figures/interface_detector_overlay.png) | Exploratório | Documenta limitações visíveis em pequenas bolhas e interfaces abertas. |
| [Perfis radiais de brilho](../intermediate/exploratory/figures/radial_interface_profiles.png) | Exploratório / método | Diferencia brilho central e sombras. As faixas vermelhas devem ser lidas conforme a legenda do detector. |
| [Quadros brutos representativos](../intermediate/exploratory/figures/representative_raw_images.png) | Exploratório | Redundante com quadro × panorama; não usar ambos no corpo principal. |

As imagens foram redimensionadas para caber nos painéis. Não medir dimensões físicas diretamente nessas figuras. Os SVG dos três gráficos de resultados permitem exportação vetorial.

## Dois gráficos novos relevantes

1. **Sensibilidade de α e flip-rate ao limiar.** Dois painéis alinhados por aquisição, com ponto em 0,70 e extensão mínimo–máximo entre 0,60 e 0,80. Unir aquisição e limiar nas tabelas [de α](../intermediate/sensitivity/panorama_void_sensitivity.csv) e [de flip-rate](../intermediate/sensitivity/phase_mask_flip_sensitivity.csv), usando as 25 elegíveis. Responde quanto a posição no gráfico depende do limiar. A extensão representa sensibilidade, não erro estatístico.

2. **Panorama, α por coluna e trocas locais alinhados.** Mostrar plug, churn e anular ondulado no mesmo intervalo de pixels. Usar as máscaras para calcular **αⱼ = Σᵢ(dᵢⱼ/D)²**, com **D=53 pixels**, e **fⱼ = Σᵣ |Bᵣ,ⱼ − Bᵣ,ⱼ₋₁| / 53**. A média de fⱼ sobre as W−1 comparações gera F. Explica como uma imagem vira um ponto e por que muito gás pode coexistir com poucas trocas. O eixo horizontal é espacial; o preenchimento pode ocultar líquido interno.

Estas duas figuras são sugestões, ainda não geradas. Evitar curvas de ajuste, novas fronteiras de classificação ou intervalos de confiança obtidos tratando pixels como réplicas independentes.
