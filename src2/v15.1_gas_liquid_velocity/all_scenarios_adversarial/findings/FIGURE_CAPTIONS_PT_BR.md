# Legendas para publicação — investigação adicional

As duas figuras principais foram exportadas em PT-BR e inglês, com os mesmos valores, em PDF e SVG vetoriais e PNG de 600 dpi. As grandezas são ópticas e os critérios avaliam suporte e sensibilidade; nenhuma figura estabelece velocidades materiais ou superficiais independentes de gás e líquido.

## Figura 1 — Comparação óptica das 26 gravações

[PNG](../figures/all_26_selected_wall_PT_BR.png) · [PDF](../figures/all_26_selected_wall_PT_BR.pdf) · [SVG](../figures/all_26_selected_wall_PT_BR.svg)

Comparação dos módulos das velocidades ópticas estimadas no núcleo e em uma faixa junto à parede selecionada por observabilidade, para as 26 gravações incluídas. Os eixos são expressos em 10³ px/s; os valores com sinal foram negativos, indicando movimento para a esquerda. Cores e símbolos distinguem os padrões fornecidos, e os números identificam cada gravação. A parede foi escolhida em uma metade temporal e avaliada na outra; o fundo temporal foi estimado somente na metade de seleção e mantido fixo na avaliação. As coordenadas são médias das estimativas locais dos blocos conjuntamente aceitos na avaliação, ponderadas por seus números de intervalos, e não médias do vídeo inteiro ou do filme líquido completo.

Pontos preenchidos indicam sete gravações com suporte na referência e nas quatro sensibilidades, além de concordância de ambas as coordenadas. Dezoito pontos vazios passaram na referência e falharam em alguma sensibilidade. O ponto vazio com × identifica m03q05_001, que não satisfez o mínimo de três blocos conjuntos em cada metade: suas coordenadas provêm de apenas um bloco de avaliação, com 49 intervalos, e permanecem candidatas. O símbolo † indica saturação em m03q40_001, junto ao número 4 e à identificação lateral. Nenhum ponto ausente foi substituído por zero, e não foram inferidas fronteiras de padrões.

As sensibilidades compreendem troca das metades de seleção/avaliação, ausência de subtração temporal e filtragem passa-altas axial com escalas gaussianas de 32 e 8 px. Cada alternativa precisou apresentar suporte próprio e concordância das duas coordenadas dentro de max(2 px/quadro, 20% do módulo da coordenada de referência). Esses limites são critérios empíricos; os símbolos não representam intervalos de confiança nem níveis de exatidão experimental. O núcleo e a parede são regiões ópticas, sem identificação independente das velocidades materiais das fases.

[Dados dos pontos](../results/final_all_26_plot_data.csv) · [Sensibilidades](../results/final_sensitivity.csv).

## Figura 2 — Disponibilidade dos descritores por método e gravação

[PNG](../figures/method_coverage_PT_BR.png) · [PDF](../figures/method_coverage_PT_BR.pdf) · [SVG](../figures/method_coverage_PT_BR.svg)

Cobertura dos critérios de suporte dos métodos de movimento óptico nas 26 gravações. As colunas correspondem à referência anterior com duas faixas externas (12 pares), nova correlação em blocos locais com ambas as paredes (20), parede selecionada com fundo estimado na metade de seleção (25), correlação de intensidade dentro de máscaras com erosão de 3 px (7) e 1 px (8), análise espectral com critérios conjuntos completos (2), e parede selecionada com suporte e concordância em todas as sensibilidades (7). O traço indica ausência de suporte suficiente na definição correspondente, não velocidade zero nem ausência física de uma fase.

Os métodos diferem em grandeza, amostragem temporal, critérios e suporte espacial. Assim, maior cobertura não demonstra maior exatidão física, e as colunas não podem ser combinadas seletivamente para produzir velocidades de fase por gravação. As duas erosões mascaradas têm seis pares em comum com velocidades compatíveis; esse consenso não é uma coluna adicional da figura. Nenhuma coluna identifica independentemente duas velocidades materiais.

[Tabela completa de cobertura](../results/final_method_coverage.csv).

## Diagnósticos preservados

O [gráfico espectral](../figures/spectral_all_candidates_PT_BR.png) mantém seus 26 pares candidatos e critérios próprios; não substitui os valores da Figura 1. Os [exemplos de máscaras](../figures/phase_local_masks_examples.png) mostram quadros nativos e regiões candidatas, sem constituírem anotações materiais de referência. Consulte os relatórios [espectral](spectral_review.md), [mascarado](phase_local_review.md) e [adversarial](adversarial_identifiability_PT_BR.md) antes de usar esses diagnósticos em um artigo.
