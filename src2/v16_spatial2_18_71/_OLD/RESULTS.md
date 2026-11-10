# Resultados calculados — v16_spatial2_18_71

Esta síntese foi gerada diretamente das tabelas atualizadas. As figuras e as estatísticas atuais usam o plano completo, sem COI nem máscara conservadora. As médias da alternative WPS incluem todos os lags; o denominador da alternative WCO permanece no lag zero. As decisões e demonstrações estão em [REPORT.md](REPORT.md). A referência anterior, com entrada 20–59 e também sem máscaras, foi preservada em [v16_spatial2](../v16_spatial2/README.md).

## Cobertura e reprodução

- 31 panoramas existentes auditados por hash; a reconstrução pixel a pixel foi verificada na análise anterior e não foi repetida. As 31 planilhas originais foram conferidas contra a faixa 20–59; elas não contêm o novo sinal.
- A nova entrada foi extraída diretamente de cada PNG: média de panorama[18:72, :], isto é, 54 linhas, de 18 a 71 inclusive. Os panoramas e os deslocamentos de reconstrução permaneceram os mesmos.
- 26 registros no escopo confirmado; 87 pares intrapadrão; três versões por par.
- Os cinco registros excluídos continuam na auditoria individual, mas não nas tabelas por padrão abaixo.
- 200 comprimentos de onda logaritmicamente espaçados, de 50 a 3000 px; cálculos feitos em todos os lags, sem redução da resolução do sinal.

## Autocoerência: 26 registros elegíveis

O valor em lag zero foi 1 em todos os métodos, por identidade algébrica. As violações abaixo são máximos buscados sobre outros lags e escalas; não são erros dessa identidade.

| Método | Mediana do máximo no plano completo | Maior máximo no plano completo |
| --- | --- | --- |
| Schuck não enviesado, 18–71 | 27.1164 | 407539.1429 |
| v16 enviesado, 18–71 | 1.0000 | 1.0000 |
| v5 aplicado a 18–71 | 1.1887 | 813.8921 |

A Morlet analítica com correlação enviesada respeitou o limite no plano completo, sem limitar artificialmente os resultados a 1. Os máximos das demais implementações estão na tabela; eventuais valores acima de 1 não foram recortados. A análise inteira desta pasta já foi executada sem máscaras, assim como a referência imediata 20–59.

## alternative WCO dos 87 pares

| Método | Mediana do máximo no plano completo | Maior máximo | Pares com máximo > 1 |
| --- | --- | --- | --- |
| Schuck não enviesado, 18–71 | 53.8762 | 261165.2731 | 87 |
| v16 enviesado, 18–71 | 0.6888 | 0.9917 | 0 |
| v5 aplicado a 18–71 | 0.7078 | 1.1315 | 4 |

O máximo sobre todos os lags e escalas e o valor no lag zero descrevem aspectos diferentes. Valores altos em uma banda ou em um par não demonstram alta alternative WCO generalizada entre todas as aquisições de um padrão. A comparação detalhada encontra-se em [pair_methods.csv](results/pair_methods.csv). Nesta pasta, o comparador v5 é sua implementação aplicada à nova entrada 18–71: conserva a conversão antiga de escala e o pré-processamento antes do corte. Não corresponde ao v5 original com entrada 20–59. Seu lag é reportado na convenção do artigo após a transformação. A comparação com a entrada anterior está em [PREPROCESSING.md](PREPROCESSING.md).

## Escalas dos perfis

| Padrão | Registros | Faixa do pico alternative WPS(0,a), px | Mediana desse pico, px | Mediana do pico da média em lag, px |
| --- | --- | --- | --- | --- |
| churn | 10 | 407.7–1238.5 | 527.4 | 609.8 |
| dispersed bubble | 1 | 139.9–139.9 | 139.9 | 928.5 |
| plug | 3 | 375.5–566.7 | 521.9 | 566.7 |
| slug | 3 | 641.2–1213.3 | 682.0 | 682.0 |
| wavy annular | 9 | 175.4–480.7 | 433.7 | 480.7 |

Os picos de alternative WPS em zero e da média em lag podem ser diferentes no mesmo sinal. Não são duas estimativas intercambiáveis do comprimento de uma bolha. Não houve ajuste das bandas para separar classes, e [pattern_profiles.png](figures/pattern_profiles.png) mostra sobreposição e variabilidade. Os picos são descrições do intervalo investigado, sem teste de significância de cada máximo.

## Sensibilidade da entrada e da CWT

Os erros relativos abaixo são normas L2 entre perfis de alternative WPS no lag zero; não são erros percentuais de cada pixel nem estimativas de incerteza física. A comparação de precisão usa escalas idênticas e corrige a constante de normalização da cmor.

| Comparação | Mediana do erro relativo | Maior erro relativo |
| --- | --- | --- |
| Linha 40 fixa versus média das linhas 18–71 | 26.922% | 93.973% |
| Remoção da média versus tendência linear | 0.258% | 4.935% |
| Linhas 45–71 versus média das linhas 18–71 | 11.398% | 31.134% |
| PyWavelets, precisão 12, mesmas escalas | 0.099% | 0.774% |
| PyWavelets, precisão 16, mesmas escalas | 0.012% | 0.084% |
| Linhas 18–44 versus média das linhas 18–71 | 18.329% | 38.622% |

No teste de interpolação de vizinhanças das emendas, 25 dos 26 registros permitiram o teste com menos de metade dos pixels alterada. O maior erro relativo do perfil foi 2.989%. Isso não limita o erro causado por um deslocamento incorreto: interpolar uma emenda não recompõe uma estrutura perdida ou duplicada.

O controle [window_sensitivity.png](figures/diagnostics/window_sensitivity.png) compara começo, centro e fim com N=3171 px, preservando a amostragem. As três janelas de um registro não são réplicas independentes.

## Diagnóstico nulo condicional

Foram recalculados 199 pares de fases aleatórias por exemplo. A estatística foi o máximo em todos os lags e nas 200 escalas de 50–3000 px, aplicando o mesmo domínio ao par observado e a cada realização aleatória. O diagnóstico da entrada 20–59 permaneceu na pasta anterior; aqui os sinais e todas as realizações aleatórias foram calculados para a entrada 18–71.

| Exemplo | Par | Máximo observado | Quantil 95% do modelo | Valor-p condicional |
| --- | --- | --- | --- | --- |
| churn | m03q05_001 / m04q10_001 | 0.6990 | 0.9174 | 0.675 |
| plug | m05q01_001 / m06q01_001 | 0.7284 | 0.8839 | 0.625 |
| slug | m04q01_001 / m05q05_001 | 0.9093 | 0.9865 | 0.440 |
| wavy annular | m03q10_001 / m03q20_001 | 0.3701 | 0.6461 | 0.915 |
| Controle entre padrões | m04q01_001 / m05q01_001 | 0.8068 | 0.9878 | 0.820 |

Esses valores dependem do modelo estacionário de fases aleatórias. Não demonstram sincronização entre aquisições nem identificam causalidade. Não foram feitas decisões de rejeição corrigidas entre os cinco exemplos. A resolução Monte Carlo é 0,005 e nenhum contorno desse diagnóstico foi transferido para os mapas principais. Um máximo incomum sob esse modelo, quando observado, não estabelece ligação física entre registros independentes.

## Verificações e conclusão

Os testes numéricos, a comparação com funções extraídas do próprio v5, os hashes e a cobertura dos artefatos estão nos JSONs de validação em results. O controle de três realizações recuperou −500 amostras de lag e demonstrou que média de módulo e módulo da média não são intercambiáveis.

**Recomendação:** usar a correlação enviesada com a Morlet e a conversão de escala de Schuck, conservar a versão não enviesada como controle e interpretar alternative WPS/alternative WCS/alternative WCO como descritores de lag–escala. A avaliação do efeito da nova extração, com ganhos e perdas frente à entrada anterior, está em [EXECUTIVE_SUMMARY.md](EXECUTIVE_SUMMARY.md) e [PREPROCESSING.md](PREPROCESSING.md). Correção matemática, maior escore de similaridade e validação física são questões distintas.
