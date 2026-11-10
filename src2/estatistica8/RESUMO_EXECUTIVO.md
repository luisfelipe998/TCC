# Resumo executivo: avaliação `estatistica8`

Síntese dos resultados registrados em 26/10/2026. Base analisada: `src2/v16.3_ancoras`.

**A largura espectral da WPS foi o resultado mais promissor para distinguir os padrões *churn* e *wavy annular*.** Combinada com Pearson, atingiu **95,0% de concordância balanceada** nos ensaios avaliados. O resultado justificou continuar investigando essa combinação, mas permaneceu exploratório: nenhuma hipótese passou pelo critério estatístico final após a correção para múltiplos testes. A análise também não encontrou uma solução consistente para identificar os quatro padrões simultaneamente.

**O estudo avaliou se medidas adicionais a Pearson ajudavam a reconhecer padrões de escoamento.** Esses padrões descrevem como líquido e vapor se organizam dentro de um tubo durante a ebulição. A análise utilizou sinais de brilho extraídos de panoramas reconstruídos a partir das imagens dos ensaios. Um sinal, nesse contexto, é uma sequência de números que representa as variações de claro e escuro ao longo do panorama.

Foram utilizados **25 ensaios**, distribuídos da seguinte forma:

| Padrão | Explicação simplificada no contexto do estudo | Ensaios |
| --- | --- | ---: |
| *Churn* | Escoamento agitado e irregular, com forte desorganização da interface entre líquido e vapor. | 10 |
| *Wavy annular*, ou anular ondulado | Vapor no centro do tubo e uma camada de líquido junto à parede, com ondas na interface. | 9 |
| *Plug* | Bolhas alongadas separadas por trechos de líquido, sem bolhas pequenas dispersas nesses trechos. | 3 |
| *Slug* | Alternância entre bolsões de vapor e trechos de líquido, com interface mais instável e pequenas bolhas dispersas. | 3 |

Cada ensaio foi representado por uma **janela**, isto é, um trecho de 2.048 pixels do panorama, centrado em uma posição de referência chamada **âncora**. Foram examinadas 200 escalas entre 50 e 1.024 pixels. Essas escalas representam comprimentos de variação no sinal da imagem. Os 300 pares comparados reutilizaram os mesmos ensaios; portanto, a amostra continuou sendo de 25 aquisições.

**As medidas avaliadas descreveram aspectos diferentes dos sinais.** Foram comparadas 36 representações válidas, nas seis combinações de dois padrões e no conjunto dos quatro padrões.

| Termo técnico | Significado e papel na análise |
| --- | --- |
| **Pearson residual** | Medida de semelhança entre os formatos dos perfis espectrais após retirar sua tendência geral em escala logarítmica. Foi a referência de comparação. |
| **Wavelet e WPS** | Wavelet é uma técnica matemática que examina variações em diferentes escalas. A WPS, ou espectro de potência wavelet, descreve como a intensidade dessas variações se distribui entre escalas. |
| **Largura espectral** | Resume quanto essa distribuição se espalha entre escalas. Uma largura maior indica uma faixa central mais ampla. |
| **WCO e autocoerência** | WCO é a coerência wavelet, que mede semelhanças entre dois sinais em diferentes escalas e deslocamentos. A autocoerência compara um sinal com versões deslocadas dele mesmo. |
| **Hellinger e Wasserstein** | Distâncias matemáticas usadas para comparar distribuições espectrais. Hellinger mede diferenças nas proporções; Wasserstein considera quanto seria necessário deslocar a distribuição entre escalas para aproximar uma da outra. |

Também foram examinados pico, inclinação e combinações de descritores, isto é, números que resumem características dos sinais. Dois **kernels**, matrizes de semelhança que seriam convertidas em distâncias, foram excluídos por não atenderem às condições matemáticas dessa conversão.

**O principal ganho apareceu na comparação entre os 19 ensaios de churn e anular ondulado.** A mediana da largura WPS foi **1,032 em churn** e **0,728 em anular ondulado**. Mediana é o valor central de uma lista ordenada. As faixas se sobrepuseram parcialmente, portanto a largura isolada não separou todos os casos.

A largura foi calculada por `L = ln(λ75 / λ25)`: `λ25` e `λ75` são as escalas em que a distribuição normalizada da WPS, ponderada em escala logarítmica, acumula 25% e 75% do total; `ln` é o logaritmo natural. O índice descreve a dispersão relativa do espectro, sem unidade física. Ele não corresponde diretamente ao tamanho de bolhas em milímetros nem à energia física do escoamento.

| Abordagem | Concordância balanceada entre churn e anular ondulado |
| --- | ---: |
| Pearson residual | 58,3% |
| Largura WPS | 81,1% |
| Pearson + largura WPS, com pesos iguais | **95,0%** |

A **vizinhança balanceada**, ou concordância balanceada por vizinhos, foi calculada comparando cada ensaio com o mais semelhante entre os demais. A classe desse vizinho foi usada como previsão. O procedimento retirou, a cada vez, o próprio ensaio avaliado da base de comparação, prática chamada **validação por exclusão de um ensaio**. Na combinação, as distâncias foram ajustadas usando apenas os ensaios disponíveis para comparação.

A combinação acertou **9 dos 10 casos de churn e os 9 casos de anular ondulado**. Os 95,0% correspondem à média entre 90% e 100%, dando o mesmo peso às duas classes. A contagem total foi de 18 acertos em 19 ensaios. Como a abordagem promissora foi escolhida após explorar essa mesma amostra, esse resultado ainda não estima de forma imparcial o desempenho em novos dados. Não houve conjunto de teste independente.

A **silhueta macro**, medida de separação que dá o mesmo peso às classes, passou de **0,038 com Pearson para 0,435 com a largura**. Valores positivos indicam maior proximidade entre ensaios da mesma classe; valores negativos indicam mistura com outras classes. Entretanto, o **agrupamento não supervisionado**, que forma grupos sem usar os rótulos conhecidos, não recuperou dois grupos bem definidos. Seu **ARI**, índice de concordância com os rótulos ajustado pelo acaso, foi de apenas **0,040** para a largura; 1 representaria concordância perfeita. Assim, o bom resultado por vizinhos não se repetiu na formação automática de grupos.

**A evidência estatística ainda foi insuficiente para confirmar o achado.** O teste de permutação redistribuiu os rótulos entre os ensaios, preservando a quantidade de cada classe, e verificou quantas redistribuições produziram resultados tão fortes quanto o observado. Para churn e anular ondulado, foram examinadas todas as 92.378 distribuições possíveis.

O **valor-p** da silhueta da largura foi **0,000184 antes da correção** e **0,0926 após a correção de Holm**. Essa correção considera que testar muitas alternativas aumenta a possibilidade de encontrar um resultado favorável por acaso. A família principal reuniu 504 testes: 36 representações × 7 comparações × 2 medidas de desempenho.

O valor corrigido ficou acima do limite adotado de **0,05**. O ganho da combinação sobre Pearson também não passou pela correção, aplicada em uma família própria de testes. Isso deixou o resultado como uma evidência preliminar que merece investigação, sem demonstrar sua generalização, isto é, sua capacidade de funcionar em novos ensaios.

**Os controles mostraram estabilidade parcial e limitações relevantes.** Ao deslocar a janela em 64 pixels para cada lado, a combinação manteve concordância balanceada entre **94,4% e 95,0%**. A largura isolada variou mais, entre **66,7% e 81,1%** nas três posições. Foram reanálises dos mesmos ensaios, sem aumento da amostra.

A largura também apresentou associação com o código operacional `q`, usado para agrupar ensaios pela condição de operação, e com o comprimento original do panorama. No pequeno subconjunto com o mesmo `q=40`, formado por dois ensaios de churn e três de anular ondulado, a concordância da largura caiu para **33,3%**. Esses resultados indicam **possível confusão operacional**: parte da diferença observada pode estar ligada às condições dos ensaios, além do padrão de escoamento. Os controles realizados não eliminaram essa possibilidade.

**Os demais resultados limitaram o alcance da proposta:**

- **Quatro padrões simultaneamente:** todas as 36 representações apresentaram silhueta macro negativa na janela original. A melhor combinação por vizinhos chegou a **42,2%**, contra **26,7%** de Pearson. O ganho parcial não resolveu a identificação das quatro classes.
- **Plug versus slug:** o pico WPS, escala em que o perfil atinge seu maior valor, obteve apenas **33,3%** de concordância balanceada. Havia somente três ensaios de cada padrão.
- **WCO e autocoerência:** os resumos testados não apresentaram separação útil e estável nos objetivos prioritários. Alguns resultados favoráveis em outras comparações foram frágeis aos controles.
- **Máximos de WCO:** simulações com sinais aleatórios independentes também produziram valores altos. Um máximo elevado, isoladamente, mostrou-se insuficiente como evidência de identificação. Esse controle foi um diagnóstico do método.

**A recomendação é manter a combinação Pearson + largura WPS como candidata para churn versus anular ondulado.** A próxima etapa deve avaliá-la em novos ensaios, com o método definido antes da avaliação e melhor distribuição das condições operacionais entre as classes. Para identificar os quatro padrões ou separar plug de slug, os resultados atuais não justificam priorizar os resumos adicionais de WCO e o pico WPS.

A auditoria registrada em `validacao_final.json` teve resultado **PASS**, indicando aprovação nas verificações de cálculo e preservação dos 479 arquivos-fonte utilizados. Essa verificação sustenta a integridade do processamento; a capacidade de reconhecer novos ensaios continua dependendo de validação independente.

Fontes para consulta: [relatório técnico](RELATORIO.txt), [protocolo](protocolo.json), [métricas completas](dados/metricas.csv), [comparação com Pearson](dados/ganho_pearson.csv), [testes estatísticos](dados/permutacoes.csv), [controles de sensibilidade](dados/ganho_sensibilidades.csv) e [auditoria registrada](dados/validacao_final.json). As descrições simplificadas dos padrões seguem a [fundamentação do TCC](../../TCC/fundamentacao.tex). As instruções de reprodução estão no [README](README.txt).
