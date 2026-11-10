# v16.4 — Âncoras como entrada dos perfis WPS

**A combinação faz sentido metodológico, mas uma única âncora de 2.048 px não melhorou os resultados neste conjunto.** A pasta foi criada para uma análise exploratória controlada; não representa uma nova abordagem validada.

Abra o [relatório visual, com imagens incorporadas](RESUMO_EXECUTIVO.html) ou o [relatório em Markdown](RELATORIO.md).

| Entrada · Hellinger em 50–1.024 px | A por par | A balanceada |
| --- | --- | --- |
| Âncora · 2.048 px | 55,7% | 57,0% |
| Início · 2.048 px | 62,1% | 65,8% |
| Centro · 2.048 px | 61,8% | 67,9% |
| Fim · 2.048 px | 61,0% | 67,8% |
| Panorama completo | 76,3% | 74,6% |

A é concordância de ordenação, **não acurácia**. Todos os controles usam os mesmos 25 ensaios, com 87 pares intra e 213 inter. Os controles de início/centro/fim têm o mesmo comprimento das âncoras. A desvantagem principal persistiu nas 25 retiradas de aquisição, nas duas ponderações. Na faixa conservadora de 50–348,69 px, as âncoras chegaram a 65,3% por par, ainda abaixo do panorama completo e da janela final na mesma faixa.

A proposta foi aplicada à média por coluna do **recorte ancorado**, seguida de WPS no lag zero, normalização de área e Hellinger da v16.1. Não se utilizou a potência média sobre lags. A seleção da âncora depende do padrão conhecido; o estudo não valida classificação automática.

## Reprodução

Na raiz do projeto, com [requirements.txt](requirements.txt):

```sh
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/v16_4_ancoras_wps_mpl python3 src2/v16.4_ancoras_perfils_wps/scripts/run_analysis.py
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/v16_4_ancoras_wps_mpl python3 src2/v16.4_ancoras_perfils_wps/scripts/render_results.py
PYTHONDONTWRITEBYTECODE=1 python3 src2/v16.4_ancoras_perfils_wps/scripts/verify_results.py
```

As pastas irmãs v16.1, v16.3, v16_spatial2_18_71 e data_processed continuam necessárias; o catálogo da v16.3 também verifica suas fontes da v16.2. A primeira execução congela as entradas; mudanças posteriores são recusadas. Nenhum mapa ou resultado das análises anteriores é reescrito. O protocolo declara a triagem que antecedeu os controles.

- [Protocolo e limitações](PROTOCOLO.md) · [configuração](config.json).
- [Resultados de todos os métodos](results/concordancia.csv) · [resultados por padrão](results/por_padrao.csv).
- [Todos os pares](results/pares.csv) · [retiradas de aquisição](results/influencia_aquisicoes.csv).
- [Verificação numérica](results/verificacao.json) · [fontes congeladas](results/entradas_congeladas.json).
- [Manifesto de artefatos](results/artefatos.json) · [verificação final](results/verificacao_final.json).

Recomendação: manter o panorama completo como referência. Múltiplas âncoras/janelas por aquisição são uma hipótese para outro experimento, ainda não executado.
