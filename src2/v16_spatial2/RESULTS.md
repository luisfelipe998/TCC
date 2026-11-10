# Resultados calculados — v16

Esta síntese foi gerada diretamente das tabelas atualizadas. As figuras e as estatísticas atuais usam o plano completo, sem COI nem máscara conservadora. As médias da alternative WPS incluem todos os lags; o denominador da alternative WCO permanece no lag zero. As decisões e demonstrações estão em [REPORT.md](REPORT.md). A versão anterior está no [arquivo histórico com máscaras](results/historico_com_mascaras.zip).

## Cobertura e reprodução

- 31 panoramas reconstruídos e comparados pixel a pixel; 31 planilhas conferidas contra a média das linhas 20–59.
- 26 registros no escopo confirmado; 87 pares intrapadrão; três versões por par.
- Os cinco registros excluídos continuam na auditoria individual, mas não nas tabelas por padrão abaixo.
- 200 comprimentos de onda logaritmicamente espaçados, de 50 a 3000 px; cálculos feitos em todos os lags, sem redução da resolução do sinal.

## Autocoerência: 26 registros elegíveis

O valor em lag zero foi 1 em todos os métodos, por identidade algébrica. As violações abaixo são máximos buscados sobre outros lags e escalas; não são erros dessa identidade.

| Método | Mediana do máximo no plano completo | Maior máximo no plano completo |
| --- | --- | --- |
| Schuck não enviesado | 28.2971 | 6022.3682 |
| v16 enviesado | 1.0000 | 1.0000 |
| v5 | 1.2220 | 295.7177 |

A Morlet analítica com correlação enviesada respeitou o limite no plano completo, sem limitar artificialmente os resultados a 1. O PyWavelets do v5 e o estimador não enviesado apresentaram violações. Retirar as máscaras não corrigiu esses extremos nem modificou as equações; apenas deixou o domínio integral explícito nas figuras e estatísticas.

## alternative WCO dos 87 pares

| Método | Mediana do máximo no plano completo | Maior máximo | Pares com máximo > 1 |
| --- | --- | --- | --- |
| Schuck não enviesado | 83.0232 | 2423.8677 | 87 |
| v16 enviesado | 0.6567 | 0.9832 | 0 |
| v5 recalculado | 0.6610 | 1.0212 | 1 |

As medianas de alternative WCO no lag zero ao longo das 200 escalas são baixas; não se confirmou a descrição geral de alternative WCO 0,7–1,0 entre todas as aquisições do mesmo padrão nas notas antigas. Isso não exclui valores altos em bandas específicas de alguns pares. A comparação detalhada encontra-se em [pair_methods.csv](results/pair_methods.csv). O comparador v5 conserva o pré-processamento antes do corte; seu lag é reportado na convenção do artigo após a transformação.

## Escalas dos perfis

| Padrão | Registros | Faixa do pico alternative WPS(0,a), px | Mediana desse pico, px | Mediana do pico da média em lag, px |
| --- | --- | --- | --- | --- |
| churn | 10 | 416.2–1029.2 | 544.9 | 603.3 |
| dispersed bubble | 1 | 139.9–139.9 | 139.9 | 928.5 |
| plug | 3 | 281.5–500.9 | 360.4 | 555.2 |
| slug | 3 | 299.5–682.0 | 628.1 | 696.2 |
| wavy annular | 9 | 171.8–480.7 | 433.7 | 480.7 |

Os picos de alternative WPS em zero e da média em lag podem ser diferentes no mesmo sinal. Não são duas estimativas intercambiáveis do comprimento de uma bolha. Não houve ajuste das bandas para separar classes, e [pattern_profiles.png](figures/pattern_profiles.png) mostra sobreposição e variabilidade. Os picos são descrições do intervalo investigado, sem teste de significância de cada máximo.

## Sensibilidade da entrada e da CWT

Os erros relativos abaixo são normas L2 entre perfis de alternative WPS no lag zero; não são erros percentuais de cada pixel nem estimativas de incerteza física. A comparação de precisão usa escalas idênticas e corrige a constante de normalização da cmor.

| Comparação | Mediana do erro relativo | Maior erro relativo |
| --- | --- | --- |
| Linha central versus média das linhas 20–59 | 14.377% | 71.970% |
| Remoção da média versus tendência linear | 0.302% | 7.069% |
| Faixa inferior versus média das linhas 20–59 | 11.625% | 64.268% |
| PyWavelets, precisão 12, mesmas escalas | 0.081% | 0.826% |
| PyWavelets, precisão 16, mesmas escalas | 0.011% | 0.082% |
| Faixa superior versus média das linhas 20–59 | 47.767% | 69.125% |

No teste de interpolação de vizinhanças das emendas, 25 dos 26 registros permitiram o teste com menos de metade dos pixels alterada. O maior erro relativo do perfil foi 2.621%. Isso não limita o erro causado por um deslocamento incorreto: interpolar uma emenda não recompõe uma estrutura perdida ou duplicada.

O controle [window_sensitivity.png](figures/diagnostics/window_sensitivity.png) compara começo, centro e fim com N=3171 px, preservando a amostragem. As três janelas de um registro não são réplicas independentes.

## Diagnóstico nulo condicional

Foram recalculados 199 pares de fases aleatórias por exemplo. A estatística foi o máximo em todos os lags e nas 200 escalas de 50–3000 px, aplicando o mesmo domínio ao par observado e a cada realização aleatória. Os valores anteriores, calculados com máscara, foram substituídos.

| Exemplo | Par | Máximo observado | Quantil 95% do modelo | Valor-p condicional |
| --- | --- | --- | --- | --- |
| churn | m03q05_001 / m04q10_001 | 0.7370 | 0.8889 | 0.515 |
| plug | m05q01_001 / m06q01_001 | 0.9832 | 0.9283 | 0.005 |
| slug | m04q01_001 / m05q05_001 | 0.9363 | 0.9842 | 0.275 |
| wavy annular | m03q10_001 / m03q20_001 | 0.3793 | 0.6464 | 0.855 |
| Controle entre padrões | m04q01_001 / m05q01_001 | 0.9238 | 0.9911 | 0.440 |

Esses valores dependem do modelo estacionário de fases aleatórias. Não demonstram sincronização entre aquisições nem identificam causalidade. Não foram feitas decisões de rejeição corrigidas entre os cinco exemplos. A resolução Monte Carlo é 0,005 e nenhum contorno desse diagnóstico foi transferido para os mapas principais. Um máximo incomum sob esse modelo, quando observado, não estabelece ligação física entre registros independentes.

## Verificações e conclusão

Os testes numéricos, a comparação com funções extraídas do próprio v5, os hashes e a cobertura dos artefatos estão nos JSONs de validação em results. O controle de três realizações recuperou −500 amostras de lag e demonstrou que média de módulo e módulo da média não são intercambiáveis.

**Recomendação:** usar a correlação enviesada com a Morlet e a conversão de escala de Schuck, conservar a versão não enviesada como controle e interpretar alternative WPS/alternative WCS/alternative WCO como descritores de lag–escala. A melhoria demonstrada é matemática e numérica; não foi demonstrado ganho de classificação ou reconstrução física do escoamento.
