# Revisão crítica dos gráficos e uso no TCC

Revisão visual de **17 PNG atuais**: 12 em `figures/` e 5 em `intermediate/exploratory/figures/`. Os três SVG correspondentes foram abertos como XML válido; não contêm nomes excluídos. Os três pares PNG/SVG são exportados a partir da mesma figura pelo gerador. Arquivos históricos não representam resultados atuais.

## Correções e cuidados encontrados

- **Linha y = 40 não é um centro físico calibrado.** As posições aparentes de parede 15 e 75 têm ponto médio 45. O gráfico secundário agora identifica a linha amostrada explicitamente como `y = 40`; nomes antigos de arquivos e colunas são mantidos por compatibilidade. Não se alterou a contagem para escolher uma linha que favoreça o resultado.
- **Seleção do gráfico principal:** conforme a nova instrução, os dois diagramas de dispersão usam apenas as 25 aquisições elegíveis. A imagem saturada `m03q40_001` permanece na tabela de procedência, com `plotted=False`. Não aparece como ponto nem participa dos polígonos.
- **Apresentação:** o gerador dos dois diagramas foi simplificado e traduzido para PT-BR. Foram removidos o nome da aquisição de bolhas, as notas internas e a linha/faixa do teto de recorte. A explicação da estimativa e suas limitações devem acompanhar a figura no texto. O teto continua registrado nos dados e no diagnóstico de recorte.
- **Matrizes de confusão:** as contagens e as acurácias balanceadas estão coerentes visualmente, mas as imagens isoladas não explicam que a única aquisição de bolhas fica sem exemplos dessa classe no treinamento do seu teste. A legenda da tese deve registrar essa limitação e `n = 25`. Não usar a linha de bolhas como evidência de desempenho de um classificador treinado nessa classe.
- **Diagnósticos exploratórios corrigidos no gerador:** o gráfico de média/mediana de largura clara agora usa somente as aquisições elegíveis pela saturação do panorama e apresenta pares de pontos por aquisição, sem trajetórias entre nomes. Os valores brutos por aquisição foram preservados. O perfil radial agora identifica as faixas vermelhas como regiões externas à janela nominal de linhas 18:71; elas não são posições de parede calibradas. Ambos continuam sendo diagnósticos de métodos anteriores.
- **Escalas de imagem:** panoramas e quadros brutos aparecem com proporções reescaladas para caber nos painéis. As figuras permitem comparar a aparência, mas não medir comprimento ou geometria física diretamente. O gráfico de comparação já informa que são registros espaciais distintos.

O gerador revisado foi compilado e as três figuras foram pré-visualizadas em `chart_previews/`, usando os valores então salvos. A regeneração final deve usar a segmentação corrigida pela auditoria de método. As prévias não substituem os resultados finais.

## Inventário e destino recomendado

| Figura | Uso recomendado | Motivo / ressalva essencial |
|---|---|---|
| [Fração de vazio × flip-rate espacial](../../figures/void_fraction/void_fraction_vs_flip_rate.png) | Discussão principal | Mostra a sobreposição dos padrões. Polígonos delimitam apenas a dispersão observada; os dois eixos vêm da mesma máscara. |
| [Sensibilidade ao recorte](../../figures/void_fraction/void_crop_sensitivity.png) | Discussão principal | Demonstra a dependência de α em relação ao recorte. São cenários de sensibilidade, não intervalos de confiança. |
| [Panorama e máscara](../../figures/void_fraction/panorama_phase_masks.png) | Método e discussão principal | Torna visível o que o algoritmo marca como gás. É inspeção qualitativa de exemplos, não validação de todas as bordas. |
| [Linha fixa y = 40](../../figures/void_fraction/void_fraction_vs_center_flips.png) | Apêndice | Mostra a dependência da taxa na linha de amostragem; não é uma segunda validação da fração de vazio. |
| [Colunas de referência líquida](../../figures/void_fraction/panorama_liquid_reference_candidates.png) | Método ou apêndice | Documenta a origem da referência inferida. A aparência líquida não equivale a uma calibração adquirida separadamente. |
| [Contato com os 26 panoramas](../../figures/patterns/pattern_panorama_contact_sheet.png) | Apêndice | Permite auditar os rótulos e ver a aquisição saturada. São recortes centrais de 2.400 pixels, não cada registro inteiro. |
| [Matrizes: morfologia](../../figures/patterns/pattern_confusion_matrices.png) | Apêndice | Relevante apenas se a tese discutir identificação de padrões. Identificar a classe de bolhas sem suporte em treino. |
| [Matrizes: silhueta](../../figures/patterns/silhouette_pattern_confusion_matrices.png) | Apêndice | Mesma ressalva; concordância com rótulos não valida volume de gás. |
| [Patches de inspeção de fase](../../figures/validation/adversarial_phase_anchor_patches.png) | Apêndice | Mostra interiores selecionados manualmente; não estima o erro global de segmentação. |
| [Intervalo líquido de slug](../../figures/validation/adversarial_slug_gap_mask.png) | Apêndice | Exemplo específico de ponte espúria; ao exportar, informar aquisição, intervalo axial, esquerda original e direita com máscara. |
| [Flip óptico espacial e temporal](../../figures/flip_rate/flip_rate_spatial_vs_native.png) | Apêndice | Ajuda a distinguir normalizações, mas usa classes de brilho, não a máscara principal de gás. Barras pretas são medianas. |
| [Quimogramas nativos](../../figures/flip_rate/flip_rate_native_kymographs.png) | Apêndice | Sustenta o diagnóstico temporal óptico em um registro curto. Os picos não equivalem automaticamente à passagem de bolhas. |
| [Quadro bruto × panorama](../exploratory/figures/raw_comparison.png) | Exploratório / apêndice metodológico | Mostra diferenças de iluminação e reconstrução. As duas colunas têm extensões espaciais diferentes. |
| [Largura clara: média × mediana](../exploratory/figures/bright_width_mean_vs_median.png) | Exploratório | Explica uma limitação do descritor anterior; não é α. A apresentação corrigida usa pares de pontos das aquisições elegíveis. |
| [Detector antigo de interfaces](../exploratory/figures/interface_detector_overlay.png) | Exploratório | Documenta uma alternativa anterior, com falhas visíveis em pequenas bolhas e interfaces abertas. |
| [Perfis radiais de brilho](../exploratory/figures/radial_interface_profiles.png) | Exploratório / apêndice metodológico | Mostra que brilho central e sombras são observáveis diferentes. A legenda corrigida explica a janela demarcada em vermelho. |
| [Quadros brutos representativos](../exploratory/figures/representative_raw_images.png) | Exploratório | Redundante com a comparação quadro–panorama; não é necessário usar os dois no corpo principal. |

## Dois gráficos novos realmente úteis

1. **Sensibilidade conjunta de α e flip-rate ao limiar.** Dois painéis alinhados por aquisição, com ponto no limiar 0,70 e extensão mínimo–máximo dos limiares 0,60–0,80. Usar `panorama_void_sensitivity.csv` e `phase_mask_flip_sensitivity.csv`, unindo aquisição e limiar; manter apenas as 25 elegíveis. Responde se a posição no diagrama depende da escolha do limiar. A extensão é sensibilidade, nunca barra de erro ou intervalo de confiança.
2. **Panorama, α por coluna e contagem local de trocas alinhados.** Mostrar exemplos de plug, churn e anular ondulado no mesmo intervalo de pixels, usando as máscaras salvas. Calcular `α_j = Σ_i(d_ij/D_app)²` e `f_j = Σ_y |M_yj − M_y,j−1| / 53`; a média de `f_j` gera o flip-rate. Responde como uma imagem vira um ponto e por que ocupação elevada pode coexistir com poucas trocas. Marcar que a posição é espacial e que o preenchimento da silhueta pode ocultar líquido interno.

Evitar curvas de ajuste ou regiões de classificação adicionais, distribuições suavizadas com uma única aquisição de bolhas e barras de confiança baseadas em pixels tratados como réplicas independentes. Nenhum desses elementos é sustentado pelos dados atuais.

## Integração do gerador revisado

`summary_figures.py` conserva os oito destinos já registrados. `void_fraction_flip_plot_data.csv` continua com 26 linhas e acrescenta `plotted` e `plot_exclusion_reason`. `summary_chart_checks.json` usa `schema_version=2`, `figure_language=pt-BR`, `scatter_selection.displayed_runs`, `scatter_selection.excluded_runs` e `scatter_crop_ceiling_drawn=False`. Em `checks`, `pattern_counts_all_shown` foi substituído por `pattern_counts_retained` e `pattern_counts_displayed`; também constam `displayed_acquisitions`, `displayed_runs` e `excluded_from_scatter_runs`. As contagens e as coordenadas são derivadas dos arquivos de resultados, sem valores de α ou flip-rate escritos manualmente no gráfico.

`image_void_alternatives_metadata.json` registra `bright_width_plot_runs`, `bright_width_plot_excluded_runs` e a regra em `bright_width_plot_selection`. A seleção usa a mesma saturação do panorama da análise principal, e não a média de saturação dos quadros brutos.

Observação adicional para a revisão de fontes: o texto gerado por `image_void_alternatives.py` ainda cita Schepperle §2.2 para a segmentação e usa a expressão “axial view” ao descrever a reconstrução. A seção correta é §2.3; para este experimento, a descrição deve especificar vista lateral. Essas correções de texto foram encaminhadas ao responsável pela revisão geral.
