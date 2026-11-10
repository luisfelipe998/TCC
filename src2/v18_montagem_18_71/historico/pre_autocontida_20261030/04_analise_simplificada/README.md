# Análise simplificada: Pearson direto da potência no atraso zero

Esta é a análise estatística vigente após a simplificação autorizada pelo autor em 29/10/2026. O coeficiente de Pearson foi calculado diretamente entre os 26 perfis de potência wavelet no atraso zero, sem logaritmo nem retirada de uma reta do perfil espectral. O preparo dos sinais de intensidade e os cálculos wavelet permaneceram iguais.

A execução anterior, com comparação dos resíduos, permanece integralmente preservada em `../02_wavelet/`. Seus arquivos não foram sobrescritos. A simplificação de Pearson não recalculou a óptica. A revisão posterior da região óptica para `18:72` tem uma atualização separada da associação entre os resultados, descrita abaixo. Os testes exploratórios entre alternativas ficaram na auditoria e não foram convertidos em resultados de classificação do TCC.

## Reprodução

Executar na raiz do repositório:

```bash
python3 -B src2/v18_montagem_18_71/04_analise_simplificada/scripts/analisar.py
```

Para reproduzir somente os números e tabelas, sem redesenhar figuras:

```bash
python3 -B src2/v18_montagem_18_71/04_analise_simplificada/scripts/analisar.py --statistics-only
```

Essa opção grava o manifesto com estado `estatistica_pronta`. Para consolidar novamente o estado `complete` e as origens de todas as figuras, executar o comando completo. São necessários Python, NumPy, SciPy, Matplotlib, Pillow e PyWavelets. O programa lê a execução original e escreve somente nesta pasta e em `TCC/figuras/v18_pearson_direto/`.

### Atualização exclusiva das medidas ópticas

Depois de recalcular `../03_optica/`, atualizar somente a associação com as medidas espectrais existentes:

```bash
python3 -B src2/v18_montagem_18_71/04_analise_simplificada/scripts/atualizar_optica.py
```

Esse comando não recalcula Pearson, wavelets ou mapas, nem redesenha figuras. Atualiza `dados/optica_por_ensaio.csv`, `dados/comparacao_optica_pares_ilustrados.csv`, o campo `discussion_examples` de `resumo_tcc.json` e a rastreabilidade. Os hashes dos descritores, das tabelas espectrais, dos mapas e das oito figuras são conferidos antes e depois. O número de ensaios elegíveis vem da triagem óptica, sem impor uma contagem fixa. O manifesto permanece em estado `complete`.

O relatório está em `intermediario/atualizacao_optica.json`. O estado anterior da associação e a cópia exata do gerador que produziu as figuras preservadas ficam em `historico/antes_roi_optica_18_71_20261030/`. A proveniência das figuras aponta para essa cópia; o script atual de reprodução completa está identificado separadamente no manifesto.

## Procedimento

1. Ler os perfis `wps_zero` dos 26 arquivos individuais de `02_wavelet/intermediario/individuais/`. Cada perfil contém 200 valores, nos mesmos comprimentos de onda de 50 a 3000 pixels, calculados com o sinal completo.
2. Calcular Pearson diretamente entre cada par desses perfis. A centralização e a divisão pela norma apenas implementam a fórmula de Pearson. Não são aplicados logaritmo ou ajuste de tendência espectral.
3. Conferir os 325 coeficientes com `scipy.stats.pearsonr` e a matriz inteira com `numpy.corrcoef`. Excluir autopares dos resumos. Separar os 87 pares intrapadrão e os 238 interpadrões.
4. Calcular médias por padrão, por combinação de padrões e por ensaio. Selecionar o maior coeficiente intrapadrão de cada um dos quatro padrões com repetição. Manter o par interpadrões selecionado pelo maior valor pontual da coerência.
5. Reunir as medidas ópticas já calculadas dos 25 ensaios elegíveis. Para a comparação na Discussão, associar os novos exemplos de churn e anular ondulado às medidas dos respectivos ensaios.
6. Gerar os gráficos estatísticos. Nos perfis exibidos, dividir cada curva pelo próprio máximo somente para visualização. Calcular a mediana, em cada comprimento de onda, dessas curvas normalizadas. Os coeficientes vêm diretamente dos valores originais de potência.

## Resultados centrais

| Resumo | Valor |
|---|---:|
| Média intrapadrão | 0,7642329553866295 |
| Média interpadrões | 0,5814270296142461 |
| Diferença entre médias | 0,18280592577238342 |

| Padrão | Par ilustrado | Comprimento comum (pixels) | Coeficiente de Pearson |
|---|---|---:|---:|
| Churn | m05q40_001 / m06q20_001 | 9873 | 0,9572096997749503 |
| Anular ondulado | m03q20_001 / m04q80_001 | 22533 | 0,9882522803968814 |
| Plug | m05q01_001 / m06q01_001 | 5914 | 0,9400531818172181 |
| Slug | m04q01_001 / m05q05_001 | 5071 | 0,8245855859782485 |

O exemplo interpadrões continua sendo m04q01_001 / m05q80_001, escolhido pela coerência máxima de 0,9936204511. Seu coeficiente de Pearson direto é 0,4513811577. A seleção de churn e anular mudou porque depende de Pearson. A escolha dos quatro panoramas individuais apresentados antes dos pares no TCC não depende desse critério e foi preservada.

## Mapas e figuras

As métricas wavelet dos 325 pares vêm da execução original, pois a simplificação muda apenas a comparação dos perfis individuais. Os mapas de plug, slug e do exemplo interpadrões foram reutilizados com igualdade de hashes. Os mapas dos dois novos exemplos foram recalculados com as mesmas entradas e o mesmo núcleo `wavelet_core.py`. Máximos, mínimos, médias e posições dos picos de coerência foram conferidos contra os registros originais desses pares.

Os mapas utilizam os trechos iniciais de comprimento comum de cada par, preparados depois do corte. Suas médias laterais foram calculadas sobre todos os atrasos. Os dados foram reduzidos em colunas somente depois do cálculo das métricas, para apresentação. Os descritores usados em Pearson permanecem os dos ensaios completos.

As oito figuras e seu manifesto de origem foram copiados para `TCC/figuras/v18_pearson_direto/`. Os quatro panoramas individuais de `v18_panoramas/` e as figuras ópticas de `v18_optica/` permanecem independentes desta geração. Os rótulos dos mapas foram adaptados para “Atraso (pixels)”, “Amplitude padronizada” e “Sinais preparados”, sem alterar os dados.

## Arquivos

- `resumo_tcc.json`: resumos completos, seleção e associação com as medidas ópticas.
- `dados/pares_pearson_wco.csv`: 325 coeficientes e métricas wavelet de referência.
- `dados/perfis_espectrais.csv`: 5200 potências e valores utilizados na identidade de Pearson.
- `dados/resumo_por_padrao.csv`, `resumo_intra_inter.csv` e `medias_por_ensaio.csv`: todos os resumos usados no manuscrito.
- `dados/comparacao_optica_pares_ilustrados.csv`: valores dos dois pares utilizados na Discussão.
- `intermediario/descritores.npz`: potências, médias, valores centralizados, vetores unitários e matriz de Pearson. Não contém resíduos logarítmicos.
- `intermediario/pares/`: cinco mapas compactos, com origem e modo de obtenção.
- `intermediario/verificacao_pearson.json`: conferência numérica dos coeficientes.
- `intermediario/figuras.json`: origem dos dados desenhados e ajustes de rótulos.
- `manifesto.json`: hashes das entradas preservadas, dos scripts e dos resultados.

Os números descrevem o conjunto analisado. Não são acurácia de classificação ou evidência de significância estatística. A comparação das alternativas que motivou esta simplificação está em `auditorias/metodologia/comparacao-descritores/`, separada dos resultados atuais.
