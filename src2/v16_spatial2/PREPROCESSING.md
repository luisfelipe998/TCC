# Pré-processamento mantendo alternative WPS, alternative WCS e alternative WCO

**Estado atual:** a triagem, a comparação com v5 e a avaliação retirando uma aquisição foram recalculadas sem COI nem máscara conservadora. Os máximos por escala são buscados em todos os lags. As equações alternativas, os pesos e os pisos de potência de referência foram mantidos. Os resultados anteriores estão no [arquivo histórico com máscaras](results/historico_com_mascaras.zip).

## Conclusão

**Os ganhos frente ao v5 foram marginais no escore principal; não foi demonstrada melhora robusta e geral da similaridade nem melhor separação dos padrões.** A referência corrigida, com a média original 20–59, aumentou o escore nos 87 pares, mas a diferença média foi pequena: **+0.003827**. O maior ganho na média desses pares entre as três opções confirmadas foi o da mediana entre linhas, **+0.006256**, com **58** aumentos e **29** reduções individuais. A faixa 30–49 ganhou **+0.002577** frente ao v5, mas apresentou diferença de **-0.001250** frente à referência corrigida.

As duas alternativas foram opções do conjunto de confirmação fixado no estudo anterior. Os resultados foram apresentados sem selecionar apenas a opção ou a ponderação mais favorável nesta atualização.

## Protocolo

Foram considerados 26 registros, todos os 325 pares: 87 intrapadrão e 238 entre padrões. As definições de Schuck, o estimador enviesado, omega0=6, a amostragem em pixels e o sentido de lag do método corrigido foram mantidos. Não foram feitas deformações do eixo, buscas por segmentos favoráveis, mudanças de denominador ou clipping de alternative WCO.

Dez opções foram recalculadas em 64 escalas. O conjunto de três opções definido no estudo anterior, referência 20–59, faixa 30–49 e mediana das linhas 20–59, foi conferido em 200 escalas, juntamente com uma recomputação direta do v5 nos 325 pares. Esse conjunto reduzido foi fixado antes de observar a nova triagem. Não se trata de um ensaio prospectivo com protocolo registrado nem de um teste externo independente.

| Opção | Operação, antes das mesmas fórmulas alternativas |
| --- | --- |
| baseline | Média das linhas 20–59; corte comum; remoção da tendência linear; normalização. |
| demean_only | Mesmo sinal, removendo apenas a média, sem ajuste linear. |
| center_row | Apenas linha 40. |
| central_band20 | Média fixa das linhas 30–49; não escolhida separadamente para cada par. |
| core_row_median | Mediana entre as linhas 20–59 em cada coluna. |
| core_trimmed_mean20 | Média após retirar os 20% menores e 20% maiores valores entre as 40 linhas, por coluna. |
| spatial_gaussian2 | Filtro gaussiano com sigma=2 px no segmento de intensidade; bordas refletidas. |
| spatial_median5 | Mediana espacial de 5 pixels; bordas refletidas. |
| rank_uniform | Postos médios da intensidade divididos por N+1, antes da remoção da tendência linear. Empates são preservados. |
| winsor_1pct | Limitação pelos percentis 1 e 99 de cada segmento antes da remoção da tendência linear. |


A transformação por postos tem motivação na discussão de distribuições em Grinsted et al.; aquela discussão é de coerência convencional. Sua utilização antes das equações alternativas aqui é uma experiência nossa, não uma prescrição de Schuck. Postos e limitação de extremos mudam a informação de intensidade. Não há garantia de que valores extremos sejam ruído: podem corresponder às próprias interfaces.

## Como resumimos os mapas sem mudar as equações

O escore principal é uma média ponderada dos máximos de alternative WCO ao longo de **todos os lags** em cada escala. Os pesos brutos são a média geométrica das alternative WPS no lag zero da referência corrigida; entram escalas de 50–3000 px em que cada alternative WPS da referência tem ao menos 5% de seu próprio máximo. **Os pesos brutos e os critérios de potência são fixos entre candidatos de cada par.** A proteção numérica pode excluir escalas com denominador inválido em uma variante e, nesse caso, alterar os pesos normalizados efetivos. Não se afirma que isso ocorreu na triagem, pois essa ocorrência não foi quantificada. Não há exclusão geométrica de bordas ou restrição por sobreposição. A operação resume a alternative WCO já calculada; não é apresentada como uma equação de Schuck nem como uma nova definição de coerência.

O critério de 5% é uma escolha de análise, não um limiar de significância nem uma máscara de borda. Foram conferidos também 1% e 10%, além da alternative WCO em zero e do máximo simples no plano completo. A proteção contra denominadores numericamente nulos foi preservada. Usar o máximo por escala permite atrasos diferentes em escalas diferentes; portanto o escore não comprova que exista um único atraso físico comum. Retirar as restrições permite encontrar máximos em mais posições e escalas. Um escore maior que na versão histórica não demonstra, por si só, melhor similaridade experimental.

Na média balanceada, cada uma das quatro classes com pares intrapadrão tem peso igual; para controles, cada uma das dez combinações entre classes tem peso igual. Isso impede que os padrões agitado e wavy annular dominem toda a análise, mas pode fazer um benefício em apenas três pares do padrão pistonado parecer maior do que na média dos 87 pares. Por isso ambas as agregações são explicitadas.

## Resultados da triagem — 64 escalas

| Opção | Escore intrapadrão balanceado | Escore entre padrões | Variação intra vs referência corrigida | Variação do contraste intra−entre |
| --- | --- | --- | --- | --- |
| Faixa central: linhas 30–49 | 0.235530 | 0.284991 | +0.006634 | +0.012989 |
| Mediana das linhas 20–59 | 0.233378 | 0.290423 | +0.004481 | +0.005404 |
| Transformação uniforme por postos | 0.226321 | 0.284427 | -0.002575 | +0.004343 |
| Média aparada das linhas 20–59 | 0.230309 | 0.288640 | +0.001412 | +0.004118 |
| Referência: linhas 20–59 | 0.228897 | 0.291346 | +0.000000 | +0.000000 |
| Filtro gaussiano espacial (σ = 2 px) | 0.228901 | 0.291355 | +0.000005 | -0.000004 |
| Remoção apenas da média | 0.229542 | 0.292231 | +0.000646 | -0.000240 |
| Mediana espacial (5 px) | 0.228857 | 0.291686 | -0.000040 | -0.000380 |
| Winsorização de 1% | 0.228398 | 0.291813 | -0.000498 | -0.000965 |
| Linha central (40) | 0.235672 | 0.304969 | +0.006775 | -0.006848 |

A regra da triagem completa selecionou **Faixa central: linhas 30–49**. Essa escolha é exploratória e não redefiniu silenciosamente o conjunto de três opções da comparação com 200 escalas. A tabela apresenta também as variações nos controles entre padrões: elevar a média intrapadrão não basta se o mesmo efeito elevar ainda mais as comparações entre padrões.

## Comparação direta com v5 — 200 escalas

| Opção | Escore intra | Escore entre | Intra−entre | Média dos máximos integrais intra | Média dos máximos integrais entre |
| --- | --- | --- | --- | --- | --- |
| v5 | 0.224815 | 0.287206 | -0.062391 | 0.758758 | 2.732929 |
| Referência: linhas 20–59 | 0.228604 | 0.291655 | -0.063051 | 0.748518 | 0.804779 |
| Faixa central: linhas 30–49 | 0.235414 | 0.285280 | -0.049867 | 0.713887 | 0.843732 |
| Mediana das linhas 20–59 | 0.233077 | 0.290720 | -0.057643 | 0.720517 | 0.835689 |

O escore intra balanceado do v5 foi **0.224815**, contra **0.235414** na faixa central. Na média simples dos 87 pares, os valores foram **0.192835 → 0.195411**, diferença de **+0.002577**. Melhoraram **52** pares e pioraram **35**. O contraste entre a média intrapadrão e a média dos controles entre padrões foi **-0.049867** para a faixa central; valores negativos indicam controles com média maior.

Os máximos simples e o escore ponderado são resumos diferentes. A busca no plano completo torna o máximo mais exposto a pontos extremos, inclusive às violações numéricas do v5. Os resultados de todos os cinco resumos estão em [preprocessing_comparison_rankings.csv](results/preprocessing_comparison_rankings.csv). Escolher apenas o resumo mais favorável não seria uma comparação suficiente.

O mesmo domínio de comparação significa os mesmos lags e os mesmos comprimentos de onda nominais. O comparador literal do v5 manteve sua conversão antiga de escalas e o preparo antes do corte, diferenças metodológicas documentadas em [EQUATION_AUDIT.md](EQUATION_AUDIT.md) e [REPORT.md](REPORT.md). A comparação com o v5 não isola apenas o efeito do pré-processamento.

### Peso igual para cada um dos 87 pares

| Opção | Média dos 87 pares | Diferença vs v5 | Pares com aumento | Pares com redução | Média balanceada por padrão |
| --- | --- | --- | --- | --- | --- |
| v5 | 0.192835 | +0.000000 | 0 | 0 | 0.224815 |
| Referência: linhas 20–59 | 0.196662 | +0.003827 | 87 | 0 | 0.228604 |
| Faixa central: linhas 30–49 | 0.195411 | +0.002577 | 52 | 35 | 0.235414 |
| Mediana das linhas 20–59 | 0.199091 | +0.006256 | 58 | 29 | 0.233077 |

A mediana entre linhas apresentou escore médio **0.199091**, com **58** pares acima do v5 e **29** abaixo. A média por par e a média balanceada por padrão respondem a ponderações distintas e foram mantidas juntas para evitar atribuir aos 87 pares um benefício concentrado em classes pequenas.

Frente à referência corrigida do v16, a mediana entre linhas acrescentou somente **+0.002429** à média dos 87 pares. A faixa 30–49 apresentou diferença de **-0.001250**. Portanto, o pequeno ganho desta faixa frente ao v5 não deve ser atribuído integralmente ao pré-processamento.

### Ganhos e perdas por padrão — faixa central contra v5

| Padrão | Pares | Média v5 | Média faixa central | Diferença |
| --- | --- | --- | --- | --- |
| churn | 45 | 0.236160 | 0.236181 | +0.000022 |
| plug | 3 | 0.232501 | 0.254688 | +0.022187 |
| slug | 3 | 0.304536 | 0.321808 | +0.017272 |
| wavy annular | 36 | 0.126065 | 0.128976 | +0.002911 |

Padrões com aumento da média neste escore: **churn, plug, slug, wavy annular**. Padrões com redução: **nenhum**. A comparação permaneceu por aquisições independentes, sem procurar os segmentos que produzissem maior coerência.

### Ganhos e perdas por padrão — mediana entre linhas contra v5

| Padrão | Pares | Média v5 | Média da mediana entre linhas | Diferença |
| --- | --- | --- | --- | --- |
| churn | 45 | 0.236160 | 0.242391 | +0.006232 |
| plug | 3 | 0.232501 | 0.236459 | +0.003957 |
| slug | 3 | 0.304536 | 0.321833 | +0.017297 |
| wavy annular | 36 | 0.126065 | 0.131623 | +0.005558 |

## Escolha sem usar a aquisição avaliada

Em cada rodada, retirou-se uma aquisição e todos os seus 25 pares antes de escolher o pré-processamento. Nesta atualização, a escolha foi realizada **entre as três opções fixadas, usando 200 escalas**, a partir de preprocessing_confirm_heldout.csv. A regra exigiu aumentar a média intrapadrão e melhorar o contraste intra−entre em relação à referência corrigida; caso contrário, manteve a referência. Contagem das escolhas: **Faixa central: linhas 30–49: 24; Referência: linhas 20–59: 2**. A aquisição do padrão dispersed bubble não tem comparador intrapadrão; por isso 25 aquisições entram no resumo intrapadrão.

Usando o v5 como comparador, as mudanças balanceadas nas quatro classes com réplicas foram: intrapadrão **+0.004729**, entre padrões **-0.003527**, contraste **+0.008256**. Os pares e as rodadas compartilham registros; isso não constitui um conjunto de teste independente nem justifica tratar os pares como 325 réplicas estatísticas. Essa avaliação do conjunto reduzido não valida como melhor uma opção que apareça na nova triagem de dez alternativas.

## O que poderia elevar valores sem demonstrar maior similaridade

- Ganho, inversão de sinal e normalização da amplitude cancelam na alternative WCO; não oferecem uma melhora intrínseca.
- Deslocar um sinal pode levar um máximo para lag zero. A alternative WCO já procura relações em lag; deslocar o mapa não cria uma relação mais forte.
- Estreitar excessivamente a banda ou elevar omega0 pode aumentar coerência em sinais não relacionados. No limite de um único bin espectral comum, a Eq. 32.17 tende a 1 porque numerador e denominador ficam proporcionais ao mesmo produto de amplitudes. Isso não seria uma melhoria de similaridade demonstrada.
- Alinhar/deformar o eixo ou selecionar o segmento com maior alternative WCO mudaria a pergunta e exigiria uma validação independente da correspondência. Não foi feito.
- Suavizar numerador/denominador da maneira de Grinsted mudaria a definição. Não foi feito nas equações alternativas.
- Estimar fases/contornos ou corrigir deslocamentos do panorama pode alterar a qualidade do observável, mas precisaria de validação da imagem, e não apenas de um aumento da alternative WCO. O registro m03q40_001 continua contendo 60 deslocamentos truncados; nenhum filtro simples testado resolve essa procedência física.

## Decisão

Manter as equações alternativas matematicamente verificadas é necessário, mas não faz duas aquisições independentes de um mesmo padrão terem as mesmas posições, fases ou espaçamentos. A referência corrigida apresentou ganho pequeno e sistemático no escore principal. A mediana entre linhas apresentou o maior ganho médio entre as opções confirmadas quando os 87 pares receberam peso igual, mas esse ganho permaneceu marginal e incluiu perdas em 29 pares. A faixa central foi mais favorável com peso igual por padrão. **Nenhuma dessas observações demonstrou melhora robusta e geral da similaridade nem separação convincente dos padrões.** A retirada das máscaras simplificou a apresentação e alterou os resumos estatísticos, sem criar correspondência entre eventos.

Dados: preprocessing_comparison_pairs.csv, preprocessing_comparison_rankings.csv, preprocessing_pattern_tradeoffs.csv e preprocessing_heldout_vs_v5.csv. Figura: [comparação dos candidatos](figures/preprocessing/preprocessing_comparison.png). Mapas completos do candidato para o primeiro par de cada classe estão em figures/preprocessing, escolhidos por ordem alfabética, não pelo maior ganho.
