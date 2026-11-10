# Revisão de literatura e estatística

Escopo: 26 aquisições retidas; 25 sem saturação usadas nas estatísticas. Esta revisão não calculou velocidades e não executou o pipeline completo.

## Correções implementadas

1. **Método de Mann–Whitney registrado.** O script de flip-rate da máscara usava o padrão `auto` do SciPy sem salvar o método efetivo. Nos dados anteriores à revisão, as cinco comparações de 53 linhas eram exatas; as cinco da linha y=40 eram assintóticas devido a empates. Isso era uma omissão de método, não um erro aritmético. A função comum em `scripts/raw_flip_analysis.py:201` agora resolve explicitamente a escolha e registra método solicitado, método usado, empates e correções. `scripts/phase_mask_flip_analysis.py` também registra os dois tamanhos amostrais. As comparações ópticas antigas preservam sua escolha explícita do método assintótico.

2. **Interpretação do teste.** Os relatórios esclarecem que Mann–Whitney compara postos/distribuições, não necessariamente apenas medianas. Os testes são exploratórios e não têm correção para múltiplas comparações. Não demonstram calibração física nem eliminam diferenças entre condições experimentais.

3. **Validação interna.** Os dois relatórios de classificação agora dizem que a extração de características e o detector foram desenvolvidos com inspeção deste conjunto. Excluir uma aquisição do ajuste do classificador ou da referência evita vazamento direto nesse ajuste, mas não desfaz a inspeção prévia. Os escores são validação cruzada interna. O teste por permutação de rótulos também depende de intercambialidade e não controla automaticamente o efeito das condições de aquisição.

4. **Linha fixa y=40.** Textos de flip-rate e classificação deixam de chamar essa linha de centro físico. Os nomes antigos `center_*` permanecem nos CSVs por compatibilidade, com definição explícita. O centro geométrico não foi calibrado. A aritmética desses campos não mudou por esta correção de nome.

## Problema físico comunicado ao responsável pela segmentação

O limite superior baseado somente no contato com a borda do recorte precisa de uma hipótese adicional: ausência de gás inteiramente fora do recorte em colunas sem contato detectado. Contraexemplo: uma bolha pequena totalmente nos pixels removidos produz máscara observada vazia e indicador de contato igual a zero, embora a fração de gás real seja positiva. Portanto, o intervalo salvo é um cenário condicionado ao contato com a borda, não uma cobertura garantida de todo gás oculto. O responsável pelo método trata essa correção na reconstrução principal.

## Equações e fontes verificadas

| Item | Verificação |
|---|---|
| Xu–Fang reproduzida em de Oliveira (2021), Eq. 5, p. 8 | Conferida visualmente: o fator é `(1 + 2 Fr_lo^(-2) alpha_h^3.5)` multiplicando `((1-x)/x)(rho_v/rho_l)`. Todos os termos são adimensionais. A definição de Froude não está explicitada no artigo; a escolha usada na auditoria continua uma hipótese declarada. |
| Rouhani–Axelsson reproduzida em Dias (2025), Eq. 5, p. 5 | Conferida diretamente na página renderizada. `C0=1+0.12(1-x)` multiplica ambos os termos da mistura. O termo de arrasto escrito na auditoria tem unidade m³/kg, compatível com `x/rho_v` e `(1-x)/rho_l`. Não foi encontrado erro de transcrição. |
| Balanço de energia | `Delta x=4 q'' L/(G D h_lv)` é adimensional: numerador e denominador têm unidade W/m. Depende do estado de entrada, das perdas e das propriedades assumidas; não fornece uma calibração das imagens. |
| Reconstrução circular | `alpha_j=sum_i(d_ij/D)^2` decorre da soma das áreas de discos dividida pela área do tubo. É uma identidade somente sob as hipóteses de seções circulares, escala óptica correta e ausência de sobreposição oculta. Deve-se fazer a média das áreas por coluna, não elevar ao quadrado a largura média. |
| Flip-rate espacial | `F=sum_y sum_j abs(B[y,j+1]-B[y,j])/[H(W-1)]`, com H=53. A multiplicação por 1.000 só muda a escala de exibição. Não há unidade de Hz. |
| Schepperle (2024), seção 2.3, Eq. 2–3 | O artigo usa contornos e comparação com foto de líquido puro, com módulo da correlação de Pearson dos histogramas e limiar 0,85. O v14 usa uma referência inferida e outro algoritmo; trata-se de adaptação, não reprodução desse algoritmo nem de sua calibração. |
| Zeghloul e Al-Sarkhi (2023), resumo | Confirmados diâmetros 19–152 mm, viscosidade até 600 cP e erro relativo absoluto médio em torno de 16%. Isso não valida diretamente o tubo de 1 mm do projeto. |

As anotações de x e alpha da Fig. 7 de de Oliveira (2021) foram conferidas visualmente. As inconsistências numéricas descritas em `findings/source_audit.md` reproduzem com as densidades do CoolProp para R600a a 27 °C. A densidade calculada nesta revisão foi 548,131160603092 kg/m³ para o líquido saturado e 9,649731786620363 kg/m³ para o vapor saturado. Esses números servem apenas à auditoria do artigo; não foram usados para substituir os resultados das imagens.

## Verificação executada

- Cinco scripts modificados compilam e seus módulos necessários importam.
- Três controles cobrem amostra pequena sem empates, amostra com empates e grupos maiores.
- As dez comparações de flip-rate da máscara reproduzem os p-valores anteriores quando recebem exatamente os dados arquivados anteriores à revisão.
- As sete comparações ópticas reproduzem os p-valores anteriores, preservando o método assintótico.
- Controle de classe ausente no centróide e controles de comprimento/censura de segmentos passaram.

Evidências: [verificações estatísticas](literature_stats_checks.json) e [página de Dias com a Eq. 5](literature_pages/dias2025_p05.png). Mudanças posteriores no detector podem alterar legitimamente as máscaras e os p-valores; a preservação acima isola apenas as correções de estatística e documentação.

**Limite:** esta revisão não transforma validação computacional em validação do volume real de gás.

