# Análise wavelet após montagem nas linhas 18 a 71

Esta etapa utilizou exclusivamente a variante `wavelet` dos 26 panoramas de `../01_panoramas/manifesto.json`. A faixa usada na montagem e na extração do sinal foi `18:72`, com 54 linhas. A regra histórica das emendas da análise wavelet foi preservada. A variante óptica corrigida não foi utilizada nesta etapa.

O núcleo de `src2/v16_spatial2_18_71/scripts/wavelet_core.py` foi carregado sem alteração, com hash registrado. Mantiveram-se a inversão de intensidade, retirada da tendência linear e padronização, correlação com divisor N, Morlet com frequência central 6, 200 comprimentos de onda entre 50 e 3.000 pixels e todos os atrasos. Não se aplicou cone de influência nem máscara geométrica.

Os 26 perfis individuais utilizaram os sinais completos. O coeficiente de Pearson foi calculado entre os resíduos do logaritmo da potência no atraso zero após ajuste contra o logaritmo do comprimento de onda. Cada um dos 325 campos de pares utilizou os prefixos de comprimento comum, preparados depois do corte. Os escores ponderados retirados do manuscrito não foram recalculados.

Foram selecionados os quatro maiores coeficientes de Pearson intrapadrão, um por grupo com mais de um ensaio, e o maior valor pontual de coerência entre os 238 pares interpadrões. Os quatro pares intrapadrão permaneceram iguais. O exemplo interpadrões passou a ser `m04q01_001 / m05q80_001`, slug e anular ondulado.

Somente esses cinco pares conservaram mapas compactos para apresentação. Todas as métricas e a seleção foram calculadas nos mapas completos, antes da redução do número de colunas. As potências médias laterais também utilizaram todos os atrasos.

## Artefatos

- `resumo_tcc.json`: números atuais, seleção, extensões e resumos usados no manuscrito.
- `comparacao_tcc.json`: comparação dos valores anteriores e atuais, incluindo as interpretações.
- `04_estatistica/dados/`: tabelas dos 325 pares, perfis, matriz, médias por grupo e por ensaio.
- `intermediario/`: entradas, descritores, perfis individuais, métricas dos 325 pares e cinco mapas compactos.
- `intermediario/auditoria/`: verificações numéricas e apresentação das figuras.
- `manifesto.json`: hashes de todos os arquivos desta etapa, fontes, entradas e scripts.
- `../../../TCC/figuras/v18_wavelet/origens.json`: origem das dez figuras do manuscrito.

## Verificação

Os 325 coeficientes foram conferidos com `scipy.stats.pearsonr`, além de `numpy.corrcoef`. A retirada da tendência foi conferida por uma fórmula independente da regressão linear. A potência no atraso zero foi comparada à soma direta da Morlet nos 200 comprimentos de onda de cada ensaio e em seis comprimentos de onda por sinal de cada par. As correlações cruzadas foram conferidas por produtos escalares em sete atrasos, incluindo ambos os extremos e o atraso zero.

Todas as 200 escalas foram válidas nos 325 pares. A coerência permaneceu finita e entre zero e um, sem recorte numérico. Os maiores erros encontrados foram 1,11×10⁻¹⁵ no coeficiente de Pearson, 2,50×10⁻¹⁶ na correlação e 3,55×10⁻¹⁵ na potência dos pares. A verificação independente da equipe ficou registrada separadamente na execução global.

As figuras reutilizaram as funções históricas de apresentação. Os gráficos de perfis e distribuição conservaram as revisões de rótulos já aprovadas. Os dados desenhados foram conferidos contra os CSVs e NPZs. Os arquivos anteriores permaneceram preservados.

## Reprodução

Executar na raiz do repositório, após a conclusão do manifesto da montagem:

```bash
python3 src2/v18_montagem_18_71/scripts/recalcular_wavelet.py --workers 2
python3 src2/v18_montagem_18_71/scripts/gerar_figuras_wavelet.py
python3 src2/v18_montagem_18_71/scripts/finalizar_wavelet.py
```

Os scripts leem as fontes históricas e escrevem somente na nova execução e em `TCC/figuras/v18_wavelet/`.
