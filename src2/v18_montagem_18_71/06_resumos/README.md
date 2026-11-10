# Resumos de Pearson e coerência wavelet

Esta etapa complementa as análises anteriores sem mudar os sinais, o preparo, a CWT, os coeficientes de Pearson ou a análise óptica. Não calcula significância estatística, intervalos de confiança ou análise de sensibilidade.

## Reprodução

Sobre as saídas existentes da v18, execute a partir da raiz do projeto:

```sh
python3 -B src2/v18_montagem_18_71/scripts/resumir_pares.py --workers 2
```

Para uma saída alternativa que já contenha as etapas 02 e 04:

```sh
TCC_V18_OUTPUT=/caminho/saida python3 -B /caminho/v18/scripts/resumir_pares.py --workers 2
```

O comando completo `executar.py` inclui esta etapa por padrão. A opção `--sem-resumos` permite reproduzir somente as etapas anteriores. Os scripts e as configurações são locais ao pacote e as fontes experimentais não são acessadas por esta etapa.

## Definições

Cada par fornece um coeficiente de Pearson dos perfis de potência no atraso zero dos sinais completos, preservado da etapa 04. O mapa de coerência usa os trechos iniciais de comprimento comum, cortados antes do preparo, conforme a etapa 02.

Para obter a mediana exata, os 325 campos de coerência são recalculados com o núcleo local inalterado. As 200 escalas cobrem comprimentos de onda geometricamente espaçados de 50 a 3.000 pixels. Para um comprimento comum N, o campo contém todos os 2N−1 atrasos, de −(N−1) a N−1. As quatro estatísticas são mínimo, média aritmética, mediana e máximo de todas as células. Cada célula recebe o mesmo peso na grade logarítmica de comprimentos de onda, sem ponderação pela potência. Não se usam imagens nem mapas com colunas reduzidas para calcular estatísticas.

A mediana do campo completo difere da variável antiga `wco_mediana_zero`, que resume somente as 200 escalas no atraso zero. A coerência média do campo é um escalar por par e também difere das curvas laterais de potência média, calculadas separadamente em cada escala.

A máscara de validade exige ambas as potências de normalização acima de 10⁻¹⁰ dos seus respectivos máximos. As 200 escalas devem permanecer válidas em todos os 325 pares. O domínio inclui as extremidades e mantém a correlação com divisor N. Não se aplica um novo corte, cone de influência ou limiar de coerência.

Nos resumos entre pares, cada par recebe o mesmo peso. Não se concatenam células de mapas de comprimentos diferentes. As tabelas apresentam mínimo, média, mediana, máximo e os quartis necessários à figura para Pearson e para a coerência média de cada par. Os quartis usam interpolação linear. Os 325 pares compartilham 26 aquisições e não são tratados como 325 réplicas independentes.

## Produtos e verificação

- `dados/pares_pearson_coerencia.csv`: 325 pares com identificação, comprimentos, domínio, Pearson e quatro estatísticas da WCO.
- `dados/resumo_por_combinacao.csv`: 14 combinações de padrões, com uma linha para cada uma das duas métricas.
- `dados/resumo_intra_inter.csv`: resumos globais dentro e entre padrões.
- `dados/pares_ilustrados.csv`: os cinco pares selecionados na etapa 04, sem nova seleção.
- `resumo_tcc.json`: as mesmas tabelas e as diferenças entre médias intrapadrão e interpadrões.
- `intermediario/pares/`: registros individuais com estatísticas, assinaturas e conferência contra a etapa 02.
- `verificacao.json`: comparação das métricas anteriores, preservação de Pearson e hashes das etapas 02/03/04 e configuração.
- `manifesto.json`: fontes, parâmetros, versões, hashes e estimativas de memória.
- `figuras/distribuicoes_pearson_coerencia.pdf/.png`: distribuições dos pares. Caixas mostram os percentis 25 e 75, linhas mostram medianas, losangos mostram médias e hastes abrangem mínimo e máximo.

Mínimo, média, máximo, estatísticas no atraso zero e posição do pico são conferidos contra todos os 325 registros da etapa 02. A mediana utiliza `numpy.median` sobre o campo completo. O maior campo WCO tem 9.013.000 valores, aproximadamente 68,8 MiB. Os três campos complexos mais a WCO ocupam aproximadamente 481,3 MiB por processo, antes dos temporários. Os campos completos são descartados após os resumos. Reexecuções podem reutilizar registros da etapa 06 com assinaturas idênticas.

Para regenerar apenas a figura a partir da tabela verificada:

```sh
python3 -B src2/v18_montagem_18_71/scripts/figuras_resumos.py
```

Após mudanças na figura, execute novamente `resumir_pares.py` para renovar o manifesto, reutilizando os registros numéricos quando as fontes de cálculo permanecerem idênticas.
