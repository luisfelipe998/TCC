# Análise adversarial da metodologia v16

A auditoria identificou problemas de qualidade na entrada de **m03q40_001**, correspondente ao antigo índice **13**, e limitações na interpretação dos resumos. Não foi identificado erro de implementação nas definições alternativas de WPS, WCS e WCO para a adaptação declarada, nem no cálculo de Pearson. As evidências não demonstraram que o desvio desse ensaio tivesse uma causa única.

O conjunto original de 26 aquisições foi mantido como referência, com ressalva explícita para m03q40. A sensibilidade à retirada do registro foi calculada e apresentada. A seleção não foi alterada para elevar a correlação. Os nomes mXXqYY foram incluídos nos dois eixos das matrizes de Pearson e de escore de WCO. As versões anteriores dos documentos e das figuras foram preservadas em `_OLD/revisao_adversarial_20261027`.

## Escopo e evidência independente

Foram relidos os **15.531 quadros JPEG** dos 31 registros, reestimados os **15.500 deslocamentos** entre quadros e reconstruídos os 31 panoramas em memória. Todos os panoramas coincidiram pixel a pixel com os PNGs utilizados. A extração das 54 linhas 18–71 coincidiu exatamente com as entradas da análise, e as planilhas foram conferidas contra sua faixa original 20–59.

O estimador vetorizado de deslocamento, implementado com produtos por FFT e somas por sobreposição, foi confrontado com o laço literal de Pearson em três intervalos de cada registro e em todos os 60 intervalos substituídos de m03q40. Também foram conferidos os rótulos dos arquivos locais contra a tabela manual: os 26 ensaios incluídos concordaram. O registro já excluído m05q05_002 apresentou rótulo local dispersed bubbles e rótulo slug na tabela da condição m05q05; essa divergência foi preservada como ressalva de procedência, sem contaminar o conjunto analisado.

Reproduzir um panorama demonstrou a execução correta do algoritmo existente. Não demonstrou que as larguras das tiras correspondessem ao deslocamento físico correto das interfaces. Os mapas de 325 pares foram novamente submetidos à verificação final de integridade, denominadores, escores, cobertura e coerência com as entradas.

Fontes: [auditoria das 31 entradas](../04_estatistica/dados/auditoria_entradas.csv), [deslocamentos reestimados](../intermediario/auditoria/adversarial/deslocamentos_reestimados.csv) e [identificação dos quadros por SHA-256](../intermediario/auditoria/adversarial/quadros_sha256.json).

## Metadados de aquisição

As configurações `_config.xsv` e `_config.mii` concordaram quanto à identificação, dimensões, quantidade de quadros, taxa e tempo de exposição. Entre os 26 ensaios incluídos, foram encontrados **62, 127 e 132 µs** de exposição e **4.000, 5.000 e 7.000 quadros/s**. O código `Brightness` variou entre 9 e 20; não foi interpretado como unidade fotométrica nem convertido em uma função de correção. A profundidade nominal foi 8 bits. Os códigos `SensorGain`, `GammaCorrection` e `Contrast` foram preservados na tabela, sem supor que representassem fatores lineares.

As diferenças de aquisição constituíram fatores de confusão entre registros. Normalizar média e desvio padrão não desfez eventual borramento por exposição, saturação ou transformação não linear anterior. Entretanto, **m03q40 e m03q20 tiveram as mesmas configurações nos campos verificados**, incluindo 132 µs e 4.000 quadros/s. Assim, não foi identificado um ajuste de câmera exclusivo que explicasse o ensaio 13. A fonte da iluminação, a cena e a resposta óptica não foram calibradas por esses metadados.

Fonte: [parâmetros nominais e hashes dos 62 arquivos de configuração](../04_estatistica/dados/metadados_aquisicao.csv).

## Conferência das fórmulas de Schuck e Bodmann

Foram conferidas visualmente as páginas impressas **411, 413 e 414**, páginas 3, 5 e 6 do [PDF original](../../../assets/schuck-bardo-wavelet-coherence.pdf), e confrontadas com o [núcleo implementado](../scripts/wavelet_core.py).

| Elemento | Referência | Resultado da conferência |
| --- | --- | --- |
| CWT e Morlet | Eqs. 32.4 e 32.6 | Normalização temporal a⁻¹ᐟ² e Morlet com ω₀ = 6 corretas; quadratura complexa independente confirmou módulo e fase. |
| Escala e comprimento de onda | Eq. 32.8 | λ = 4πa/(6 + √38), com amostragem de 1 px. |
| WPS alternativa | Eqs. 32.14–32.15 | Módulo da CWT da autocorrelação; elevar novamente ao quadrado alteraria a definição. |
| WCS alternativa | Eq. 32.16 | CWT complexa da correlação cruzada. |
| WCO alternativa | Eq. 32.17 | Módulo quadrático da WCS dividido pelas WPS no lag zero, na mesma escala. |
| Estimador de correlação | Eqs. 32.21–32.23 | O PDF utilizou N−|m|; a análise utilizou N. Houve uma adaptação explícita, não reprodução literal desses estimadores. |
| Número de realizações | Eqs. 32.18–32.20 | K = 1 por aquisição. A Eq. 32.20 imprimiu média de módulos, enquanto 32.16 definiu média complexa; a diferença não alterou a WCO para K = 1. |

O estimador com divisor N correspondeu à correlação de registros com extensão nula. Sua autocorrelação apresentou espectro não negativo, e o multiplicador de Fourier da Morlet utilizado também foi não negativo. Essa estrutura permitiu normalizar a correlação filtrada sem abandonar o limite de coerência em 1. A troca automática para N−|m| não constituiu uma correção neutra: modificaria o estimador e suas propriedades. O núcleo numérico foi preservado.

Os controles independentes obtiveram erro máximo de **3.61e-15** na correlação e **2.67e-17** na CWT temporal complexa. A autocoerência no lag zero foi unitária. Um atraso conhecido de 100 amostras foi recuperado como lag −100, conforme a convenção x[n+m]·conjugado(y[n]). A inversão da ordem dos sinais refletiu o mapa em lag. Os 26 perfis experimentais foram recalculados por quadratura temporal, com diferença máxima de **3.55e-15** frente à implementação espectral.

Foi corrigida a descrição do sinal de deslocamento em [pipeline_v11.py](../../pipeline_v11.py): um deslocamento de B para a direita produziu lag negativo no código. A correção foi documental; não foram invertidos sinais nem reconstruídas entradas com uma nova convenção.

Fonte: [controles numéricos das equações](../intermediario/auditoria/adversarial/equacoes.json).

## Diagnóstico do ensaio m03q40

O rótulo **wavy annular** concordou com a classificação manual. Não foi encontrada troca de identificação, de arquivo, de ordem da matriz ou de faixa na extração. O Pearson médio de m03q40 com os outros oito ensaios do padrão foi **0,7424**; entre esses outros oito, a média foi **0,9508**. Seu pico de WPS no lag zero ocorreu em **175,4 px**, enquanto os demais picos ficaram entre **353,0 e 480,7 px**.

Foram encontradas duas evidências de qualidade relevantes:

- **Substituição de deslocamentos:** 48 estimativas nulas e 12 negativas foram convertidas em largura de 1 px, totalizando **60/500 intervalos, ou 12%**. Entre os demais 25 ensaios incluídos, somente m04q01 apresentou uma substituição, 1/500. Os 60 pixels resultantes ocuparam apenas 0,289% do panorama, mas isso não limitou o erro espacial a essa proporção: a extensão omitida ou duplicada não foi conhecida. No controle histórico de deslocamento com verificação em outras linhas e consistência temporal, nenhum dos 60 intervalos recebeu suporte; esse controle foi citado como evidência anterior, sem reexecução nesta auditoria.
- **Acúmulo no teto de intensidade:** **25,43%** dos pixels das linhas 18–71 do panorama atingiram 255. Nos quadros JPEG originais completos, a mesma faixa apresentou **15,71%** no teto. Os denominadores dessas porcentagens foram diferentes, pois o panorama utilizou tiras centrais. A evidência antecedeu a normalização do sinal, mas não identificou a origem do teto entre sensor, processamento anterior e cena. Informações eventualmente perdidas por saturação não seriam recuperadas por normalização posterior.

A primeira substituição ocorreu no intervalo dos quadros **280→281**, após **14.008 px**. Recalculou-se Pearson com esse mesmo prefixo para os nove ensaios wavy annular. Também foram executados controles de faixa espectral, remoção de tendência, faixa vertical e equalização do histograma. Cada controle utilizou a mesma regra para todos os nove registros.

| Controle exploratório | m03q40 versus outros oito | Entre os outros oito |
| --- | ---: | ---: |
| Referência integral, linhas 18–71 | 0,7424 | 0,9508 |
| Prefixo comum de 14.008 px, anterior à primeira substituição | 0,8491 | 0,9350 |
| Panoramas completos, λ de 50–1.024 px | 0,8675 | 0,9465 |
| Centralização e normalização, sem remover a tendência linear do sinal | 0,7405 | 0,9508 |
| Média das linhas 20–59 | 0,8754 | 0,9508 |
| Equalização global do histograma, estimada por aquisição | 0,8763 | 0,9249 |

O prefixo anterior à primeira substituição reduziu a discrepância, mas não a eliminou. Nesse trecho, 24,99% dos pixels já atingiram 255. A mudança para as linhas 20–59 e a equalização também reduziram o desvio; a retirada da tendência linear teve efeito pequeno. Os controles mostraram sensibilidade à faixa óptica, ao intervalo observado e ao processamento, sem separar suas contribuições da variabilidade física do escoamento. Não foi demonstrado que as substituições explicassem sozinhas o resultado, nem que houvesse erro no rótulo físico do padrão. Os controles não foram adotados como novo método nem escolhidos pelo maior Pearson.

Fontes: [controles de m03q40](../04_estatistica/dados/controles_m03q40.csv), [panorama e sinal](../01_panoramas/wavy_annular_m03q40_001.png), [controle histórico do deslocamento](../_OLD/displacement/README.md).

## Sensibilidade à exclusão

A retirada de m03q40 eliminou oito pares intrapadrão e 17 entre padrões, deixando **25 ensaios, 79 pares intrapadrão e 221 interpadrões**. Os descritores dos ensaios restantes não foram recalibrados. As médias por par e por grupo foram recalculadas com as mesmas definições da referência.

| Medida | Conjunto | Intra, por par | Inter, por par | Intra − inter | Intra − inter, balanceado |
| --- | --- | ---: | ---: | ---: | ---: |
| Pearson | 26 ensaios | 0,8697 | 0,8038 | 0,0659 | 0,1453 |
| Pearson | 25, sem m03q40 | 0,8826 | 0,8277 | 0,0549 | 0,1403 |
| Escore de WCO | 26 ensaios | 0,1877 | 0,2635 | -0,0759 | -0,0643 |
| Escore de WCO | 25, sem m03q40 | 0,1956 | 0,2687 | -0,0730 | -0,0660 |

O Pearson dentro de wavy annular aumentou de **0,9045 para 0,9508**. Entretanto, o contraste global de Pearson diminuiu nas duas ponderações, porque as correlações entre padrões também aumentaram. O contraste balanceado do escore de WCO continuou negativo e diminuiu de −0,0643 para −0,0660. Assim, retirar o ensaio deixou o bloco anular mais homogêneo, mas não demonstrou melhora global de discriminação.

As falhas de procedência espacial justificaram uma ressalva de qualidade e a análise com e sem o registro. Não forneceram um deslocamento corrigido para os intervalos problemáticos nem estabeleceram uma regra geral de exclusão validada antes da inspeção dos resultados. Manteve-se o conjunto original como referência auditável; a variante de 25 ensaios permaneceu explicitamente exploratória. Uma futura filtragem por qualidade deverá aplicar critérios de aquisição ou reconstrução a todos os registros, sem usar o valor de Pearson como critério.

Fonte: [tabela completa de sensibilidade](../04_estatistica/dados/sensibilidade_exclusao_m03q40.csv).

## Limitações do método principal

1. **Amostragem espacial reconstruída:** o eixo em pixels dependeu de translação óptica estimada, montagem de tiras variáveis e observações feitas em tempos distintos. A reprodução exata da montagem não validou uma fotografia espacial instantânea nem uma velocidade física. A média de intensidade não correspondeu diretamente à fração de fase.
2. **Faixa e comprimento:** λ de 50–3.000 px foi aplicada a todos os ensaios. O menor panorama comportou somente 1,057 comprimentos de onda de 3.000 px, enquanto o maior comportou 7,693. Isso limitou a comparação nas grandes escalas. m03q40 teve 20.743 px e seu pico esteve próximo de 175 px; o seu desvio não correspondeu simplesmente ao registro mais curto.
3. **Maximização da WCO:** o escore procurou um máximo em todos os lags para cada escala, admitindo um lag diferente a cada vez. Esse procedimento não mediu um único alinhamento físico e incluiu a seleção de extremos. O peso e o piso de 5% foram resumos próprios do estudo, não equações de Schuck nem limiares de significância. A ausência de COI não removeu efeitos de registro finito.
4. **Alcance de Pearson:** a comparação utilizou resíduos de log(WPS) após remover uma reta em log(λ). A amplitude, a inclinação removida e a fase não foram avaliadas por esse coeficiente. Um valor alto indicou semelhança desse descritor e não comprovou identidade do padrão, correspondência de eventos ou acurácia de classificação.
5. **Dependência e composição dos grupos:** os pares compartilharam aquisições e as 200 escalas foram correlacionadas. Não houve 325 observações independentes. Dispersed bubble teve uma única aquisição, impossibilitando estimar sua repetibilidade. As estatísticas permaneceram descritivas.
6. **Aquisição heterogênea:** exposição, taxa de quadros e código de brilho variaram entre registros. A análise em pixels não eliminou possíveis diferenças de resposta óptica, borramento e saturação. Esses efeitos não foram isolados experimentalmente do padrão de escoamento.

## Correções e rastreabilidade

Foram corrigidos os identificadores nos eixos das duas matrizes, a descrição invertida do sinal do deslocamento no pipeline de origem e a omissão das ressalvas de entrada nos documentos ativos. A metodologia explicitou novamente a adaptação do estimador de Schuck e as limitações do suporte nas grandes escalas. Os cálculos existentes de WPS, WCS, WCO e Pearson não precisaram de alteração.

Os controles foram implementados em [auditoria_adversarial.py](../scripts/auditoria_adversarial.py), com [resultado estruturado](../intermediario/auditoria/adversarial/resultado.json). A reprodução foi incluída no [roteiro de execução](../scripts/README.txt). As evidências científicas e o sucesso dos testes numéricos foram mantidos distintos: os testes validaram as operações implementadas, enquanto as ressalvas de entrada e interpretação permaneceram abertas.
