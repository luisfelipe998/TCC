# Pré-processamento mantendo alternative WPS, alternative WCS e alternative WCO

**Estado atual:** a entrada principal passou a ser a média das linhas 18–71. A triagem, a comparação com a implementação v5 e a avaliação retirando uma aquisição foram recalculadas sem COI nem máscara conservadora. Os máximos por escala são buscados em todos os lags. As equações alternativas e as regras dos pesos/pisos foram mantidas; a comparação entre faixas usa pesos fixos de 20–59, enquanto a comparação entre pré-processamentos usa a referência 18–71. A análise anterior está preservada em [v16_spatial2](../v16_spatial2/EXECUTIVE_SUMMARY.md).

## Conclusão

A ampliação da média de 20–59 para 18–71 diminuiu o escore médio de alternative WCO dos 87 pares intrapadrão de 0.196662 para 0.195636, diferença de -0.001026 (-0.52%), usando os mesmos pesos brutos da referência 20–59. Houve 41 aumentos, 46 reduções e 0 empates. Frente ao v5 original com entrada 20–59 (0.192835), a diferença da média 18–71 corrigida foi +0.002801 (+1.45%). 

O contraste balanceado intrapadrão menos entre padrões passou de -0.063051 para -0.059861. Não foi obtida uma melhora da média dos 87 pares ao ampliar a faixa principal. Os padrões com aumento da média foram plug, slug; com redução, churn, wavy annular. Com peso igual por padrão, a média intrapadrão passou de 0.228604 para 0.234430, variação de +0.005826.

Nesta execução, a piora na média dos 87 pares foi pequena. A vantagem numérica remanescente sobre o v5 original foi marginal em magnitude e não foi criada pela ampliação da faixa. O resultado é misto entre padrões; não sustenta uma melhora global da alternative WCO. Isso não implica mudanças pequenas em cada par: o exemplo plug predeterminado apresentou uma redução substancial do máximo integral, detalhada no diagnóstico do resumo executivo.

Na comparação entre implementações, usando os pesos nativos da referência 18–71, o escore médio dos 87 pares foi 0.184176 para a implementação v5 aplicada à nova entrada e 0.187651 para a referência corrigida de Schuck, diferença de +0.003474. O controle fixo 30–49 teve média 0.185980, com 49 aumentos e 38 reduções frente ao v5 aplicado a 18–71. A mediana da faixa 18–71 teve média 0.192501, com 47 aumentos e 40 reduções. A maior média por par entre as três opções confirmadas foi de Mediana das linhas 18–71. Esta comparação entre implementações e controles não mede, sozinha, a vantagem da nova faixa sobre a antiga 20–59.

As duas alternativas foram opções do conjunto de confirmação fixado no estudo anterior. Os resultados foram apresentados sem selecionar apenas a opção ou a ponderação mais favorável nesta atualização.

## Efeito da extração 18–71 frente à análise anterior 20–59

Esta comparação preserva, em cada par, os pesos brutos e os pisos de alternative WPS da referência corrigida **20–59**. A análise anterior e a nova usam todos os lags, sem COI nem máscara conservadora. Assim, a troca da referência de potência causada pela nova extração não é contabilizada como ganho da faixa. Os mapas e as equações alternativas continuam sendo os de cada entrada; somente seu resumo é calculado com a mesma referência. Pesos brutos fixos só garantem pesos normalizados idênticos quando as escalas elegíveis também têm denominadores numericamente válidos em todas as opções comparadas; as contagens da proteção numérica desta execução estão em [PREPROCESSING.md](PREPROCESSING.md).

Na tabela, a faixa principal é 20–59 na coluna anterior e 18–71 na coluna nova. O controle 30–49 permanece exatamente na mesma posição nas duas análises. Ganhos/perdas são diferenças superiores a 1e−12 em módulo; diferenças menores são empates numéricos.

| Operação | Média dos 87 pares, 20–59 | Média dos 87 pares, 18–71 | Variação | Aumentos | Reduções | Empates |
| --- | --- | --- | --- | --- | --- | --- |
| Implementação v5 sobre a faixa principal | 0.192835 | 0.192020 | -0.000815 | 41 | 46 | 0 |
| Média da faixa principal; Schuck corrigido | 0.196662 | 0.195636 | -0.001026 | 41 | 46 | 0 |
| Controle fixo: média 30–49 | 0.195411 | 0.195411 | +0.000000 | 0 | 0 | 87 |
| Mediana entre as linhas da faixa principal | 0.199091 | 0.200907 | +0.001816 | 45 | 42 | 0 |

O escore do controle fixo 30–49 reproduziu a análise anterior, com diferença absoluta máxima de 0. Isso verifica a continuidade desse controle sob os mesmos pesos. A comparação da média corrigida de 18–71 com o v5 original inclui tanto a ampliação da faixa quanto as correções da implementação; a comparação entre as duas médias corrigidas isola a mudança de entrada.

### Efeito por padrão — média corrigida, pesos fixos de 20–59

| Padrão | Pares | Média 20–59 | Média 18–71 | Variação |
| --- | --- | --- | --- | --- |
| churn | 45 | 0.240848 | 0.237582 | -0.003266 |
| plug | 3 | 0.236329 | 0.245193 | +0.008863 |
| slug | 3 | 0.308429 | 0.326802 | +0.018373 |
| wavy annular | 36 | 0.128809 | 0.128142 | -0.000667 |

Padrões com aumento da média: **plug, slug**. Padrões com redução: **churn, wavy annular**. A média conjunta dos 87 pares atribui a cada padrão peso proporcional à sua quantidade de pares. A média balanceada atribui peso igual a cada um dos quatro padrões com pares e passou de **0.228604** para **0.234430**. Por isso as duas agregações podem variar em sentidos opostos. O padrão dispersed bubble possui uma única aquisição elegível e não entra nessas médias intrapadrão. As duas agregações foram mantidas; os pares continuam compartilhando aquisições.

### Sensibilidade ao resumo do mapa — média corrigida

| Resumo | Média dos 87 pares, 20–59 | Média dos 87 pares, 18–71 | Variação intra | Variação do contraste balanceado |
| --- | --- | --- | --- | --- |
| Máximos por escala, piso de 1% | 0.199966 | 0.199223 | -0.000744 | +0.003592 |
| Máximos por escala, piso de 5% | 0.196662 | 0.195636 | -0.001026 | +0.003190 |
| Máximos por escala, piso de 10% | 0.194852 | 0.194019 | -0.000833 | +0.003446 |
| alternative WCO no lag zero | 0.089638 | 0.089826 | +0.000188 | -0.000973 |
| Máximo no plano completo, sem pesos | 0.654949 | 0.658998 | +0.004049 | -0.018473 |

Os escores nas seções de pré-processamento abaixo usam os pesos nativos da referência **18–71**, fixos entre candidatos. Não devem ser subtraídos diretamente dos antigos escores 20–59 para atribuir um efeito à extração. As duas tabelas respondem a comparações distintas e mantêm as mesmas alternative WPS, alternative WCS e alternative WCO.

Dados da comparação: [pares com pesos fixos](results/band_comparison_pairs.csv), [resumos](results/band_comparison_summary.csv), [mudanças por resumo](results/band_change_summary.csv), [mudanças individuais](results/band_change_pairs.csv) e [verificação estruturada](results/band_comparison.json). Os pares compartilham aquisições; estas diferenças não constituem um teste independente de superioridade.

## Protocolo

Foram considerados 26 registros, todos os 325 pares: 87 intrapadrão e 238 entre padrões. As definições de Schuck, o estimador enviesado, omega0=6, a amostragem em pixels e o sentido de lag do método corrigido foram mantidos. Não foram feitas deformações do eixo, buscas por segmentos favoráveis, mudanças de denominador ou clipping de alternative WCO.

Dez opções foram recalculadas em 64 escalas. O conjunto de três opções definido no estudo anterior, referência 18–71, faixa 30–49 e mediana das linhas 18–71, foi conferido em 200 escalas, juntamente com a implementação v5 aplicada à entrada 18–71 nos 325 pares. Esse conjunto reduzido foi fixado antes de observar a nova triagem. Não se trata de um ensaio prospectivo com protocolo registrado nem de um teste externo independente.

| Opção | Operação, antes das mesmas fórmulas alternativas |
| --- | --- |
| baseline | Média das linhas 18–71; corte comum; remoção da tendência linear; normalização. |
| demean_only | Mesmo sinal, removendo apenas a média, sem ajuste linear. |
| center_row | Apenas linha 40. |
| central_band20 | Média fixa das linhas 30–49; não escolhida separadamente para cada par. |
| core_row_median | Mediana entre as linhas 18–71 em cada coluna. |
| core_trimmed_mean20 | Média após retirar 10 valores de cada extremo entre as 54 linhas, por coluna; floor(0,20 × 54), mantendo 34 valores. |
| spatial_gaussian2 | Filtro gaussiano com sigma=2 px no segmento de intensidade; bordas refletidas. |
| spatial_median5 | Mediana espacial de 5 pixels; bordas refletidas. |
| rank_uniform | Postos médios da intensidade divididos por N+1, antes da remoção da tendência linear. Empates são preservados. |
| winsor_1pct | Limitação pelos percentis 1 e 99 de cada segmento antes da remoção da tendência linear. |


A transformação por postos tem motivação na discussão de distribuições em Grinsted et al.; aquela discussão é de coerência convencional. Sua utilização antes das equações alternativas aqui é uma experiência nossa, não uma prescrição de Schuck. Postos e limitação de extremos mudam a informação de intensidade. Não há garantia de que valores extremos sejam ruído: podem corresponder às próprias interfaces.

## Como resumimos os mapas sem mudar as equações

O escore principal é uma média ponderada dos máximos de alternative WCO ao longo de **todos os lags** em cada escala. Os pesos brutos são a média geométrica das alternative WPS no lag zero da referência corrigida; entram escalas de 50–3000 px em que cada alternative WPS da referência tem ao menos 5% de seu próprio máximo. **Os pesos brutos e os critérios de potência são fixos entre candidatos de cada par.** A proteção numérica pode excluir escalas com denominador inválido em uma variante e, nesse caso, alterar os pesos normalizados efetivos. As ocorrências foram contadas e são apresentadas abaixo. Não há exclusão geométrica de bordas ou restrição por sobreposição. A operação resume a alternative WCO já calculada; não é apresentada como uma equação de Schuck nem como uma nova definição de coerência.

O critério de 5% é uma escolha de análise, não um limiar de significância nem uma máscara de borda. Foram conferidos também 1% e 10%, além da alternative WCO em zero e do máximo simples no plano completo. A proteção contra denominadores numericamente nulos foi preservada. Usar o máximo por escala permite atrasos diferentes em escalas diferentes; portanto o escore não comprova que exista um único atraso físico comum. As duas faixas foram comparadas no mesmo domínio integral; diferenças de amplitude e distribuição de potência causadas pela extração também precisam ser consideradas ao interpretar o resumo.

Na média balanceada, cada uma das quatro classes com pares intrapadrão tem peso igual; para controles, cada uma das dez combinações entre classes tem peso igual. Isso reduz a influência da quantidade de pares de churn e wavy annular, mas aumenta o peso relativo dos três pares de plug e dos três de slug frente à média dos 87 pares. Por isso ambas as agregações são explicitadas.

### Proteção numérica do denominador

Nenhuma combinação de par/opção teve escalas descartadas pelo denominador numérico nestes três lotes. O descarte numérico e o piso de potência são critérios diferentes: o primeiro evita divisão por valores numericamente nulos; o segundo define as escalas usadas no resumo ponderado. Nenhum deles restringe geometricamente o mapa em lag. A contagem de escalas é somada sobre as combinações de par/opção, que não são observações independentes.

| Lote | Combinações par/opção | Escalas por combinação | Combinações com descarte numérico | Total de escalas descartadas | Escalas no escore de 5%, mínimo–máximo |
| --- | --- | --- | --- | --- | --- |
| Triagem | 3250 | 64 | 0 | 0 | 37–64 |
| Confirmação Schuck | 975 | 200 | 0 | 0 | 117–200 |
| Comparador v5, entrada 18–71 | 325 | 200 | 0 | 0 | 117–200 |

## Resultados da triagem — 64 escalas

| Opção | Escore intrapadrão balanceado | Escore entre padrões | Variação intra vs referência corrigida | Variação do contraste intra−entre |
| --- | --- | --- | --- | --- |
| Faixa 30–49 (controle fixo) | 0.246900 | 0.298694 | -0.002202 | +0.012242 |
| Transformação uniforme por postos | 0.258377 | 0.313340 | +0.009275 | +0.009073 |
| Mediana das linhas 18–71 | 0.258262 | 0.315769 | +0.009160 | +0.006529 |
| Média aparada das linhas 18–71 | 0.252369 | 0.312632 | +0.003267 | +0.003775 |
| Remoção apenas da média | 0.249918 | 0.313818 | +0.000817 | +0.000137 |
| Referência: linhas 18–71 | 0.249102 | 0.313139 | +0.000000 | +0.000000 |
| Filtro gaussiano espacial (σ = 2 px) | 0.249106 | 0.313146 | +0.000004 | -0.000004 |
| Mediana espacial (5 px) | 0.249092 | 0.313214 | -0.000010 | -0.000085 |
| Winsorização de 1% | 0.248830 | 0.314282 | -0.000272 | -0.001415 |
| Linha 40 (controle fixo) | 0.245063 | 0.322327 | -0.004039 | -0.013227 |

A regra da triagem completa selecionou **Transformação uniforme por postos**. Essa escolha é exploratória e não redefiniu silenciosamente o conjunto de três opções da comparação com 200 escalas. A tabela apresenta também as variações nos controles entre padrões: elevar a média intrapadrão não basta se o mesmo efeito elevar ainda mais as comparações entre padrões.

## Implementação v5 na entrada 18–71 — 200 escalas

| Opção | Escore intra | Escore entre | Intra−entre | Média dos máximos integrais intra | Média dos máximos integrais entre |
| --- | --- | --- | --- | --- | --- |
| v5 aplicado a 18–71 | 0.244365 | 0.308589 | -0.064224 | 0.767470 | 3.062082 |
| Referência: linhas 18–71 | 0.249002 | 0.313324 | -0.064322 | 0.727288 | 0.802022 |
| Faixa 30–49 (controle fixo) | 0.246895 | 0.298755 | -0.051861 | 0.713887 | 0.843732 |
| Mediana das linhas 18–71 | 0.258232 | 0.315989 | -0.057757 | 0.734353 | 0.815905 |

O escore intra balanceado do v5 foi **0.244365**, contra **0.246895** na faixa fixa. Na média simples dos 87 pares, os valores foram **0.184176 → 0.185980**, diferença de **+0.001804**. Melhoraram **49** pares e pioraram **38**. O contraste entre a média intrapadrão e a média dos controles entre padrões foi **-0.051861** para a faixa fixa; valores negativos indicam controles com média maior.

Os máximos simples e o escore ponderado são resumos diferentes. A busca no plano completo torna o máximo mais exposto a pontos extremos, inclusive às violações numéricas do v5. Os resultados de todos os cinco resumos estão em [preprocessing_comparison_rankings.csv](results/preprocessing_comparison_rankings.csv). Escolher apenas o resumo mais favorável não seria uma comparação suficiente.

O mesmo domínio de comparação significa os mesmos lags e os mesmos comprimentos de onda nominais. O comparador literal do v5 manteve sua conversão antiga de escalas e o preparo antes do corte, diferenças metodológicas documentadas em [EQUATION_AUDIT.md](EQUATION_AUDIT.md) e [REPORT.md](REPORT.md). A comparação com o v5 não isola apenas o efeito do pré-processamento.

### Peso igual para cada um dos 87 pares

| Opção | Média dos 87 pares | Diferença vs v5 | Pares com aumento | Pares com redução | Média balanceada por padrão |
| --- | --- | --- | --- | --- | --- |
| v5 aplicado a 18–71 | 0.184176 | +0.000000 | 0 | 0 | 0.244365 |
| Referência: linhas 18–71 | 0.187651 | +0.003474 | 87 | 0 | 0.249002 |
| Faixa 30–49 (controle fixo) | 0.185980 | +0.001804 | 49 | 38 | 0.246895 |
| Mediana das linhas 18–71 | 0.192501 | +0.008325 | 47 | 40 | 0.258232 |

A mediana entre linhas apresentou escore médio **0.192501**, com **47** pares acima do v5 e **40** abaixo. A média por par e a média balanceada por padrão respondem a ponderações distintas e foram mantidas juntas para evitar atribuir aos 87 pares um benefício concentrado em classes pequenas.

Frente à referência corrigida do v16, a mediana entre linhas apresentou diferença de **+0.004850** na média dos 87 pares. A faixa 30–49 apresentou diferença de **-0.001671**. Essas diferenças permitem separar o efeito dos controles do efeito da implementação numérica.

### Ganhos e perdas por padrão — faixa fixa contra v5

| Padrão | Pares | Média v5 | Média faixa fixa | Diferença |
| --- | --- | --- | --- | --- |
| churn | 45 | 0.215767 | 0.216164 | +0.000397 |
| plug | 3 | 0.277326 | 0.286454 | +0.009128 |
| slug | 3 | 0.362284 | 0.359548 | -0.002736 |
| wavy annular | 36 | 0.122082 | 0.125413 | +0.003330 |

Padrões com aumento da média neste escore: **churn, plug, wavy annular**. Padrões com redução: **slug**. A comparação permaneceu por aquisições independentes, sem procurar os segmentos que produzissem maior coerência.

### Ganhos e perdas por padrão — mediana entre linhas contra v5

| Padrão | Pares | Média v5 | Média da mediana entre linhas | Diferença |
| --- | --- | --- | --- | --- |
| churn | 45 | 0.215767 | 0.222833 | +0.007066 |
| plug | 3 | 0.277326 | 0.294531 | +0.017205 |
| slug | 3 | 0.362284 | 0.385570 | +0.023286 |
| wavy annular | 36 | 0.122082 | 0.129994 | +0.007911 |

## Escolha sem usar a aquisição avaliada

Em cada rodada, retirou-se uma aquisição e todos os seus 25 pares antes de escolher o pré-processamento. Nesta atualização, a escolha foi realizada **entre as três opções fixadas, usando 200 escalas**, a partir de preprocessing_confirm_heldout.csv. A regra exigiu aumentar a média intrapadrão e melhorar o contraste intra−entre em relação à referência corrigida; caso contrário, manteve a referência. Contagem das escolhas: **Mediana das linhas 18–71: 23; Faixa 30–49 (controle fixo): 2; Referência: linhas 18–71: 1**. A aquisição do padrão dispersed bubble não tem comparador intrapadrão; por isso 25 aquisições entram no resumo intrapadrão.

Usando o v5 como comparador, as mudanças balanceadas nas quatro classes com réplicas foram: intrapadrão **+0.005884**, entre padrões **+0.001807**, contraste **+0.004077**. Os pares e as rodadas compartilham registros; isso não constitui um conjunto de teste independente nem justifica tratar os pares como 325 réplicas estatísticas. Essa avaliação do conjunto reduzido não valida como melhor uma opção que apareça na nova triagem de dez alternativas.

## O que poderia elevar valores sem demonstrar maior similaridade

- Ganho, inversão de sinal e normalização da amplitude cancelam na alternative WCO; não oferecem uma melhora intrínseca.
- Deslocar um sinal pode levar um máximo para lag zero. A alternative WCO já procura relações em lag; deslocar o mapa não cria uma relação mais forte.
- Estreitar excessivamente a banda ou elevar omega0 pode aumentar coerência em sinais não relacionados. No limite de um único bin espectral comum, a Eq. 32.17 tende a 1 porque numerador e denominador ficam proporcionais ao mesmo produto de amplitudes. Isso não seria uma melhoria de similaridade demonstrada.
- Alinhar/deformar o eixo ou selecionar o segmento com maior alternative WCO mudaria a pergunta e exigiria uma validação independente da correspondência. Não foi feito.
- Suavizar numerador/denominador da maneira de Grinsted mudaria a definição. Não foi feito nas equações alternativas.
- Estimar fases/contornos ou corrigir deslocamentos do panorama pode alterar a qualidade do observável, mas precisaria de validação da imagem, e não apenas de um aumento da alternative WCO. O registro m03q40_001 continua contendo 60 deslocamentos truncados; nenhum filtro simples testado resolve essa procedência física.

## Decisão

A ampliação da média de 20–59 para 18–71 diminuiu o escore médio de alternative WCO dos 87 pares intrapadrão de 0.196662 para 0.195636, diferença de -0.001026 (-0.52%), usando os mesmos pesos brutos da referência 20–59. Houve 41 aumentos, 46 reduções e 0 empates. Frente ao v5 original com entrada 20–59 (0.192835), a diferença da média 18–71 corrigida foi +0.002801 (+1.45%). 

O contraste balanceado intrapadrão menos entre padrões passou de -0.063051 para -0.059861. Não foi obtida uma melhora da média dos 87 pares ao ampliar a faixa principal. Os padrões com aumento da média foram plug, slug; com redução, churn, wavy annular. Com peso igual por padrão, a média intrapadrão passou de 0.228604 para 0.234430, variação de +0.005826.

Nesta execução, a piora na média dos 87 pares foi pequena. A vantagem numérica remanescente sobre o v5 original foi marginal em magnitude e não foi criada pela ampliação da faixa. O resultado é misto entre padrões; não sustenta uma melhora global da alternative WCO. Isso não implica mudanças pequenas em cada par: o exemplo plug predeterminado apresentou uma redução substancial do máximo integral, detalhada no diagnóstico do resumo executivo.

Na comparação entre implementações, usando os pesos nativos da referência 18–71, o escore médio dos 87 pares foi 0.184176 para a implementação v5 aplicada à nova entrada e 0.187651 para a referência corrigida de Schuck, diferença de +0.003474. O controle fixo 30–49 teve média 0.185980, com 49 aumentos e 38 reduções frente ao v5 aplicado a 18–71. A mediana da faixa 18–71 teve média 0.192501, com 47 aumentos e 40 reduções. A maior média por par entre as três opções confirmadas foi de Mediana das linhas 18–71. Esta comparação entre implementações e controles não mede, sozinha, a vantagem da nova faixa sobre a antiga 20–59.

Todas as opções confirmadas apresentaram média entre padrões maior ou igual à intrapadrão neste escore. Os pares compartilham aquisições, e os resultados permanecem uma avaliação interna. Manter a formulação de Schuck não cria correspondência entre eventos de registros independentes.

Dados: preprocessing_comparison_pairs.csv, preprocessing_comparison_rankings.csv, preprocessing_pattern_tradeoffs.csv e preprocessing_heldout_vs_v5.csv. Figura: [comparação dos candidatos](figures/preprocessing/preprocessing_comparison.png). Mapas completos do candidato para o primeiro par de cada classe estão em figures/preprocessing, escolhidos por ordem alfabética, não pelo maior ganho.
