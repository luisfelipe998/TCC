# estatistica5 — Pearson direto, todos contra todos

**Versão atual:** 26 ensaios, 325 pares únicos, escala de Pearson de −1 a 1, sem células não calculadas. Mediana r intrapadrão: 0.8301; entre padrões: 0.7546. As distribuições apresentam sobreposição; correlação alta de perfil não identifica isoladamente o padrão.

- [Resumo executivo com gráficos](RESUMO_EXECUTIVO.md) · [PDF para compartilhar](RESUMO_EXECUTIVO.pdf)
- [Passo a passo e fórmulas](METODOLOGIA_TODOS_CONTRA_TODOS.md)
- [Matriz principal PNG](todos_contra_todos/figures/matriz_pearson_todos.png) · [PDF vetorial](todos_contra_todos/figures/matriz_pearson_todos.pdf)
- [325 coeficientes](todos_contra_todos/results/pearson_325_pares.csv) · [Matriz CSV](todos_contra_todos/results/matriz_pearson.csv)
- [Índice dos ensaios](todos_contra_todos/results/indice_ensaios.csv) · [Descritores](todos_contra_todos/results/descritores_por_ensaio.npz)
- [Validação](todos_contra_todos/results/validation.json) · [Manifesto](todos_contra_todos/results/manifest.json)

## Reproduzir a versão atual

Na raiz do repositório:

    PYTHONDONTWRITEBYTECODE=1 python3 src2/estatistica5/todos_contra_todos/run_all_pairs.py

A execução confere as fontes e gera todos os resultados, cinco gráficos em PNG/PDF, o resumo executivo em Markdown/PDF e o manifesto. A condição principal usa os primeiros 3.171 px de cada registro; um descritor fixo por ensaio permite a comparação integral. O controle com registros completos é apresentado separadamente.

## Versão anterior

A [análise dos 87 pares](README_87_PARES.md) foi preservada no [arquivo histórico completo](historico_87_pares.zip), incluindo o README original e seu manifesto conferido. Os arquivos antigos em results/ e figures/ e o script run_pearson.py também permanecem disponíveis. Nela havia cinza porque as comparações entre padrões não eram incluídas. Ela usa distância normalizada e trechos dependentes do par; não é a versão atual. O script histórico regrava o README antigo se executado; use o comando acima para reproduzir a comparação completa.

Nenhum ensaio foi reclassificado. O único dispersed bubbles agora participa das comparações com as outras classes, mantendo seu rótulo original. A [explicação da ECDF anterior](EXPLICACAO_ECDF.md) permanece histórica e não descreve os novos coeficientes.
