# Correspondência entre metodologia e implementação

Esta versão reuniu o processamento descrito no capítulo 3 do TCC, para os cinco padrões puros e os 26 ensaios de `config/ensaios.json`. As operações numéricas, a seleção dos exemplos e a produção das figuras ficaram no próprio pacote. A única entrada experimental externa foi `data_processed`.

| Trecho do capítulo 3 | Implementação local | Produto |
|---|---|---|
| 3.1 — etapas | `scripts/metodologia_fluxograma.py` | Fluxograma dos três caminhos |
| 3.2 — dados e montagem | `scripts/reconstruir.py`, `config/ensaios.json` | 26 panoramas espectrais e 26 ópticos; deslocamentos de todos os quadros |
| 3.3 — preparo dos sinais | `scripts/wavelet_core.py`, `scripts/recalcular_wavelet.py` | Sinais completos e trechos comuns preparados |
| 3.4 — correlações, CWT e resumos dos campos | Os dois scripts anteriores e `scripts/resumir_pares.py` | 26 perfis individuais e mínimo, média, mediana e máximo da coerência de cada par |
| 3.5 — comparação dos ensaios e exemplos | `04_analise_simplificada/scripts/analisar.py`, `scripts/resumir_pares.py` e `scripts/figuras_resumos.py` | Matriz, cinco mapas de pares, estatísticas por grupo e distribuições de Pearson e da coerência média |
| 3.6 — máscara e medidas ópticas | `scripts/optical_model.py`, `scripts/process_optical.py` | Referência de líquido, 26 máscaras e medidas; 25 ensaios elegíveis |
| 3.6 — gráficos e comparação conjunta | `scripts/gerar_figuras_opticas.py`, `scripts/finalize_optical.py` e análise da etapa 04 | Médias, medianas, envoltórias e medidas dos pares discutidos |
| Figuras de apoio | `scripts/gerar_panoramas.py`, `scripts/gerar_metodologia.py` | Exemplos reais, fluxograma e ilustrações sintéticas identificadas |

## Dados e montagem

Foram utilizadas 501 imagens `ImgA000000.jpg` a `ImgA000500.jpg`, com 1016 × 80 pixels, de cada ensaio. Os arquivos `_config.xsv` forneceram as taxas de aquisição e exposições; os arquivos `__<padrão>.txt` forneceram os rótulos. Incluíram-se dez ensaios de churn, nove de anular ondulado, três de plug, três de slug e um de bolhas dispersas. Rótulos compostos e anular liso ficaram fora do conjunto.

Em todas as operações que dependiam da faixa radial, adotaram-se as linhas de índices 18 a 71, inclusive: o recorte Python `18:72`, com 54 linhas. Na montagem, calculou-se a média dessas linhas por coluna. Compararam-se perfis de quadros consecutivos por Pearson nas posições sobrepostas, para deslocamentos inteiros de −200 a 200 pixels. Escolheu-se o maior coeficiente; em empate, prevaleceu o primeiro deslocamento da busca crescente. O cálculo vetorizado foi acompanhado de verificações por soma direta e de busca literal nos máximos numericamente ambíguos.

Retiraram-se tiras centrais com largura igual ao deslocamento escolhido, adotando-se um pixel quando ele foi nulo ou negativo. O último quadro repetiu a largura anterior. Os panoramas espectrais substituíram a última coluna anterior à emenda pela segunda coluna seguinte somente quando ambas as tiras tinham pelo menos quatro colunas. Os panoramas ópticos justapuseram as mesmas tiras sem substituir pixels. Os arquivos conservaram as 80 linhas; o processamento utilizou as 54 linhas definidas acima.

O sinal de intensidade **I[n]** foi a média dessas 54 linhas na coluna de índice **n**. O panorama constituiu uma composição de instantes distintos; seu eixo permaneceu em pixels, sem conversão em distância ou velocidade do escoamento.

## Preparo, correlações e CWT

Para cada perfil individual, preparou-se o sinal completo. Para cada par, cortaram-se primeiro os dois sinais brutos no menor comprimento e prepararam-se os trechos resultantes. Aplicaram-se inversão da intensidade (255 − I), remoção da tendência linear por mínimos quadrados e padronização pelo desvio padrão populacional. Não se alinharam eventos antes do corte.

As autocorrelações e a correlação cruzada foram calculadas com `scipy.signal.correlate`, modo completo, e divididas pelo número de amostras **N**, mantendo esse divisor em todos os atrasos. Assim, os mapas cobriram os atrasos de −(N − 1) a N − 1.

Utilizou-se a Morlet complexa com frequência angular central **ω₀ = 6**. A grade reuniu 200 comprimentos de onda **λ** geometricamente espaçados entre 50 e 3000 pixels. Cada comprimento de onda correspondeu a uma escala da CWT e a uma linha do mapa. A escala **a** foi obtida dividindo-se o comprimento de onda **λ** pelo fator de Fourier **4π/(6 + √38)**, com espaçamento de um pixel entre amostras. A CWT foi aplicada às correlações, com prolongamento por zeros e cálculo por FFT; a implementação reservou uma margem de oito desvios da envoltória gaussiana da maior escala.

A potência individual foi o módulo da CWT da autocorrelação. O espectro cruzado foi a CWT complexa da correlação cruzada. A coerência wavelet alternativa mediu a associação entre os dois sinais de cada par de ensaios. Seu cálculo dividiu o módulo ao quadrado desse espectro pelo produto das duas potências individuais no atraso zero, em cada escala. O código verificou o denominador antes da divisão, e todas as 200 escalas foram válidas nos 325 pares. Os resumos utilizaram todos os atrasos, sem excluir as extremidades por cone de influência.

O perfil de potência no atraso zero **P** reuniu 200 valores extraídos **após transformar a autocorrelação completa**. Portanto, não se transformou apenas o valor da autocorrelação em zero. O perfil de potência média **M** reuniu a média dos módulos sobre todos os atrasos. O perfil **P** dos sinais completos forneceu o descritor de Pearson; o perfil **M** dos trechos comuns forneceu as curvas laterais dos mapas dos pares.

## Comparação espectral

Calculou-se Pearson diretamente entre os 200 valores de potência no atraso zero de dois ensaios completos. Não se aplicou logaritmo nem se removeu tendência do perfil espectral. A remoção de tendência ocorreu apenas no preparo do sinal de intensidade.

As comparações reuniram 325 pares distintos: 87 intrapadrão e 238 interpadrão, sem comparações de um ensaio consigo mesmo. Cada par recebeu o mesmo peso nos resumos. Selecionou-se o maior Pearson de cada um dos quatro padrões com pelo menos dois ensaios e o maior pico de coerência entre os 238 pares interpadrão. A seleção ocorreu antes da extração dos mapas destinados às figuras.

Os valores mínimos, médios, medianos e máximos da coerência e as coordenadas dos picos vieram dos campos completos. Cada mapa forneceu um valor de cada medida, considerando todas as células válidas com igual peso. As 200 escalas foram válidas nos 325 pares. A mediana foi calculada sobre todos os valores do mapa, não apenas sobre a coluna de atraso zero. Para armazenar e desenhar os cinco mapas selecionados, reduziu-se apenas a quantidade de colunas de exibição, preservando o atraso zero e as extremidades. Essa redução não alterou as métricas.

A etapa complementar `06_resumos` reuniu os quatro resumos de cada mapa e o coeficiente de Pearson correspondente. Para Pearson, calcularam-se mínimo, média, mediana e máximo entre os coeficientes dos pares de cada padrão e combinação de padrões. Na comparação da coerência entre grupos, essas quatro estatísticas foram calculadas entre as médias dos mapas, dando igual peso a cada par. Não se juntaram todas as células de mapas de tamanhos diferentes. Os mesmos resumos foram produzidos separadamente para os 87 pares intrapadrão e 238 interpadrões. O único ensaio de bolhas dispersas não forneceu pares intrapadrão.

A coerência wavelet média do mapa de cada par é um único valor. Ela não corresponde à potência média por escala das curvas laterais, que permanece um perfil com 200 valores. A tabela completa dos 325 pares ficou nos arquivos CSV, com sínteses nas tabelas do TCC.

A avaliação permaneceu descritiva. Os pares compartilharam ensaios e as escalas compartilharam informação do mesmo sinal; essas quantidades não foram usadas como réplicas independentes para significância ou intervalos de confiança.

## Máscaras e medidas ópticas

Todas as operações abaixo foram executadas somente nas 54 linhas de índices 18 a 71, antes de calcular percentis, desvios ou filtros:

1. Dividiu-se cada coluna pelo percentil 75 da intensidade, com divisor mínimo de 1.
2. Nos três ensaios de plug (`m05q01_001`, `m06q01_001`, `m07q05_001`), selecionaram-se as colunas com desvio padrão radial populacional até o percentil 5, incluindo empates. Calculou-se a mediana por linha em cada ensaio e, depois, a mediana dos três perfis.
3. Suavizaram-se a imagem normalizada e a referência com filtro gaussiano radial de desvio padrão de um pixel (`gaussian_filter1d`, demais parâmetros padrão). Dividiu-se a imagem pela referência, com divisor mínimo de 0,01. Razões inferiores a 0,70 identificaram sombras.
4. Procuraram-se os mínimos nas metades superior e inferior da faixa de 54 linhas. Aplicou-se fechamento morfológico de 3 × 3 pixels, com margem externa de zeros. Usou-se a conectividade padrão bidimensional do `scipy.ndimage.label`, pelos quatro vizinhos.
5. Preencheu-se o intervalo entre os limites externos das sombras quando os dois mínimos ficaram abaixo de 0,70 e pertenceram ao mesmo componente conectado, ou quando as duas sombras ocorreram simultaneamente por pelo menos 20 colunas consecutivas. Incluíram-se também componentes com pelo menos quatro pixels, sem contato com as bordas superior e inferior, preenchendo seus interiores.

A máscara binária **B** marcou gás inferido com 1 e os demais pixels com 0. Para visualização, os PNGs codificaram esses valores lógicos como 255 e 0, respectivamente. Os arquivos armazenaram 80 linhas, com zeros fora da faixa de 54 linhas. A triagem contou pixels exatamente iguais a 255 no panorama original dessa faixa. Excluiu-se dos resumos ópticos o ensaio com mais de 20% de saturação (`m03q40_001`); ele permaneceu na análise espectral. As máscaras e medidas foram conservadas para os 26 ensaios, com indicação explícita de elegibilidade.

Considerou-se o tubo circular de diâmetro interno de 1 mm, representado pelo diâmetro na imagem **D = 54 pixels**. Em cada coluna de índice **j**, cada intervalo vertical contínuo de gás teve sua largura **dᵢⱼ** tratada como diâmetro de uma seção circular, sem sobreposição em profundidade. A fração de vazio estimada na coluna **αⱼ** foi a soma das razões de áreas: **αⱼ = Σᵢ(dᵢⱼ/D)²**. A fração de vazio estimada do panorama **α** foi a média das frações de suas **W** colunas. O valor dependeu dessa representação geométrica; não correspondeu simplesmente à fração de pixels marcados.

A densidade de transições **F** foi **Ntrans/[H(W − 1)]**. O número de transições **Ntrans** contou mudanças entre pixels horizontalmente adjacentes de cada linha, o número de linhas **H** foi 54 e o número de colunas **W** foi a largura do panorama. Apresentou-se **1000F**, em transições por mil pares de pixels, sem conversão em frequência temporal.

As tabelas utilizaram medianas por padrão; as estrelas dos gráficos utilizaram médias com igual peso por ensaio. Envoltórias convexas delimitaram grupos com pelo menos três ensaios elegíveis. Verificaram-se as interseções entre envoltórias por dois procedimentos independentes. Para os pares de churn e anular ondulado discutidos no TCC, reuniram-se Pearson, coerência e medidas ópticas dos mesmos ensaios.

## Reprodutibilidade e proteção dos resultados

`executar.py` executou a sequência inteira, incluindo figuras, em um diretório novo. `scripts/comparar_resultados.py` comparou as saídas com referências locais congeladas antes da consolidação. Essas referências serviram apenas para validação; não forneceram sinais ou máscaras aos cálculos. Diferenças de pixels, seleção, elegibilidade ou resultados científicos interromperam a execução. Diferenças de ponto flutuante foram limitadas à tolerância registrada no relatório, com verificação adicional dos valores arredondados usados no texto.

`scripts/verificar_portabilidade.py` repetiu a execução a partir de uma cópia contendo somente código e configuração, com a saída vazia e acesso a arquivos do projeto restrito ao pacote copiado e a `data_processed`. Os registros de proveniência acompanharam entradas, código e produtos por SHA-256. Versões anteriores e arquivos do manuscrito não participaram do cálculo.
