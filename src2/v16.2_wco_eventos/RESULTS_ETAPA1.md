# Resultados detalhados da comparação de eventos

**Não conseguimos melhorar de forma convincente a alternative WCO com este protocolo em todo o escopo.** Para estruturas alongadas em pistonado (plug), o ganho sobre eventos de extensão distinta foi **0,0275**, uma **melhora marginal**. Em golfadas (slug), a diferença foi **0,1520**, mas dependeu de **um único par de eventos de duas aquisições** e se inverteu para **-0,0427** na faixa 30–49. O estudo fornece evidência exploratória restrita; não demonstra melhora geral, classificação dos padrões ou identidade física entre estruturas.

## Domínio e cobertura

| Padrão | Aquisições elegíveis | Eventos aceitos | Aquisições com equivalentes intrapadrão | Pares de eventos equivalentes |
| --- | --- | --- | --- | --- |
| Agitado (churn) | 10 | 0 | 0 | 0 |
| Bolhas dispersas | 1 | 0 | 0 | 0 |
| Pistonado (plug) | 3 | 20 | 3 | 18 |
| Golfadas (slug) | 3 | 6 | 2 | 1 |
| Anular ondulado | 9 | 0 | 0 | 0 |

N=2.048; amostragem 1 px; 200 comprimentos de onda de 50–1.024 px; lags principais −32…+32. Peso igual por escala logarítmica e por lag. A média conservadora é outro domínio e aparece separadamente nas tabelas. Os três eventos sem qualquer parceiro equivalente continuam catalogados, inclusive quando contribuem como controles distintos.

## Contrastes principais

| Padrão | Pares de aquisições | Equivalentes | Distintos | Diferença | Diferença frente a recortes fixos |
| --- | --- | --- | --- | --- | --- |
| Pistonado (plug) | 3 | 0,2387 | 0,2112 | 0,0275 | 0,0289 |
| Golfadas (slug) | 1 | 0,3656 | 0,2136 | 0,1520 | 0,2192 |

Os equivalentes intrapadrão ocorreram em quatro pares de aquisições: três de plug e um de slug. Há 19 pares de eventos equivalentes, 46 distintos intrapadrão, 20 equivalentes entre padrões e 80 distintos entre padrões. O limite de amostragem foi aplicado sem métricas wavelet. Os pares intermediários e não selecionados permanecem em pares_candidatos.csv.

A versão balanceada intrapadrão atribui 50% do peso a plug e 50% a slug: diferença de 0,0897 frente a eventos distintos e 0,1240 frente a recortes fixos. Esse balanceamento dá ao único par de slug o mesmo peso dos três pares de plug; por isso não fortalece a evidência nem substitui os resultados separados por padrão.

## Todas as variantes predefinidas

| Variante | Padrão | Equivalentes − distintos | Equivalentes − fixos |
| --- | --- | --- | --- |
| principal | Pistonado (plug) | 0,0275 | 0,0289 |
| principal | Golfadas (slug) | 0,1520 | 0,2192 |
| v5_mesmos_segmentos | Pistonado (plug) | 0,0264 | 0,0270 |
| v5_mesmos_segmentos | Golfadas (slug) | 0,1581 | 0,2253 |
| faixa30_49 | Pistonado (plug) | 0,0140 | 0,0297 |
| faixa30_49 | Golfadas (slug) | -0,0427 | 0,0079 |
| somente_media | Pistonado (plug) | 0,0276 | 0,0260 |
| somente_media | Golfadas (slug) | 0,1425 | 0,2119 |
| bordas_menos16 | Pistonado (plug) | 0,0275 | 0,0273 |
| bordas_menos16 | Golfadas (slug) | 0,1574 | 0,2214 |
| bordas_mais16 | Pistonado (plug) | 0,0270 | 0,0302 |
| bordas_mais16 | Golfadas (slug) | 0,1367 | 0,2120 |
| ancoras_opostas_a | Pistonado (plug) | 0,0167 | 0,0204 |
| ancoras_opostas_a | Golfadas (slug) | 0,1525 | 0,2170 |
| ancoras_opostas_b | Pistonado (plug) | 0,0281 | 0,0298 |
| ancoras_opostas_b | Golfadas (slug) | 0,1426 | 0,2157 |

As variantes de centros usam os mesmos controles fixos principais; as demais usam controles recalculados com o mesmo núcleo/faixa/preparo. Todos os pares de eventos foram recalculados em cada variante; não se escolheu a variante de maior valor por par. A coluna diferenca_conservador de [contrastes_por_padrao.csv](results/contrastes_por_padrao.csv) usa somente as escalas com suporte conservador e não substitui a banda principal.

## Separação do núcleo numérico

| Padrão | v5 nos mesmos recortes | Morlet corrigida | Diferença |
| --- | --- | --- | --- |
| Pistonado (plug) | 0,2329 | 0,2387 | 0,0058 |
| Golfadas (slug) | 0,3644 | 0,3656 | 0,0013 |

A diferença acima usa a média dos equivalentes nos mesmos pares de aquisições e a mesma entrada 18–71. O v5 preserva sua wavelet cmor, escala e orientação corrigida apenas no eixo reportado. A validação literal compara suas correlações extraídas por AST e a CWT de convolução com o cálculo usado. Os sinais são preparados dentro de cada segmento de mesmo comprimento; não se reaproveita tendência estimada do panorama inteiro.

## Comparação histórica de panoramas

A base 18–71 apresentou escore médio nativo de aproximadamente **0,187651**, contra **0,184176** do núcleo v5 na mesma faixa. Esses números resumem máximos por escala, em todos os lags e ondas de 50–3.000 px, com pesos de WPS. **Não é válido subtraí-los da média de eventos deste estudo e chamar o resultado de ganho.** O [registro histórico estruturado](results/referencia_historica.json) preserva a fonte; os controles fixos de 2.048 px fornecem a comparação com o mesmo domínio desta análise.

## Dependência, suporte e controles negativos

Os resultados de retirada de aquisição estão em [retirada_aquisicao.csv](results/retirada_aquisicao.csv), com ausência explícita quando não sobra contraste. Não são intervalos de confiança. Não houve descarte de escalas por denominador nas comparações experimentais Morlet; o domínio principal teve suporte geométrico completo. O maior valor experimental de WCO Morlet foi **0,999659**; no v5 foi **0,998794**. Não houve clipping.

Self em lag zero foi 1 por identidade algébrica, não como evidência de equivalência física. Os controles com eventos distintos e fundo comum estão em [controles_sinteticos.csv](results/controles_sinteticos.csv); comparações dentro da mesma aquisição estão em [controles_mesma_aquisicao.csv](results/controles_mesma_aquisicao.csv), fora dos resultados entre aquisições.

## Arquivos numéricos

- [Resultados dos 165 pares em oito versões](results/resultados_eventos.csv).
- [Resultados dos controles fixos](results/resultados_fixos.csv).
- [Agregação por par de aquisições](results/agregados_eventos.csv).
- [Contrastes com controles do mesmo par de aquisições](results/contrastes_pareados.csv).
- [Contrastes por padrão](results/contrastes_por_padrao.csv) e [balanceados](results/contrastes_balanceados.csv).
- [Distribuições](results/distribuicoes.csv), [cobertura](results/cobertura.csv) e [conclusão estruturada](results/conclusao.json).

As matrizes em maps são integrais nos exemplos preservados; os máximos foram calculados antes de reduzir a resolução somente para renderização. Os perfis descritivos dos 165 pares principais estão em results/perfis_eventos.npz. O catálogo, os CSVs e os scripts permitem recalcular qualquer par sem usar mapas de baixa resolução.
