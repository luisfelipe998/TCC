# Comparação de eventos por WCO alternativa

Entrega local da execução do plano v16.2, com base em **v16_spatial2_18_71**. A pasta foi renomeada para **v16.2_wco_eventos**, conforme o nome solicitado pelo usuário. O atalho v16.2 foi removido. Os dados e as análises anteriores permaneceram intactos.

Comece pelo [resumo executivo com gráficos](EXECUTIVE_SUMMARY.md) ou pelo [PDF compacto](RESUMO_EXECUTIVO.pdf). Veja o [passo a passo](PASSO_A_PASSO.md), os [resultados](RESULTS.md), o [protocolo congelado](PROTOCOL.md), as [decisões](DECISIONS.md) e o [plano original](PLAN_ORIGINAL.md).

A [nota sobre âncoras e autocorrelação](ANCORAS_E_AUTOCORRELACAO.md) explica os limites de x/y, o domínio de lag, a aplicação da CWT às correlações e as condições para uma possível extensão a churn e anular.

## Organização

- scripts: auditoria de entrada, detecção por imagem, pareamento sem wavelets, núcleo idêntico à base, validações, execução, agregação, gráficos, relatórios e verificador final.
- results: catálogo com rejeições, coordenadas e procedência; pares candidatos/selecionados; tabelas numéricas; congelamento SHA-256; validações, logs e manifesto.
- figures/observabilidade: atlas das 26 aquisições, cada um com três quadros originais fixos.
- figures/catalogo: todas as propostas das seis aquisições com morfologia alongada, inclusive rejeitadas.
- figures/eventos: evidência individual dos 26 eventos aceitos, com segmento, quadro da âncora e intensidade.
- maps: coeficientes complexos e WCO integrais dos exemplos, sem redução de lag.

## Reprodução

O Python padrão deste terminal não possuía as dependências. A execução usou **Python 3.13.2** em `/Users/i521862/.pyenv/versions/3.13.2/bin/python3`. As versões da distribuição e dos módulos estão no manifesto. PyWavelets informa distribuição 1.9.0 e módulo 1.8.0; o núcleo Morlet principal usa NumPy/SciPy. requirements.txt registra o ambiente declarado e não garante identidade binária; o manifesto registra hashes do código.

Na raiz do projeto, usando o ambiente indicado:

```sh
export MPLCONFIGDIR=/tmp/wco-eventos-mpl
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/inspect_inputs.py
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/catalog_events.py
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/validate_and_benchmark.py
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/run_analysis.py
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/summarize_results.py
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/render_figures.py
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/write_reports.py
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/verify_final.py
```

A revisão visual em results/revisao_visual.csv é um insumo congelado. Reproduzir o cálculo não produz uma nova anotação independente. Se o catálogo ou as regras forem alterados, o congelamento impedirá a execução silenciosa; uma análise nova deve ter outra versão e registrar que os resultados desta já são conhecidos. Os scripts não executam lotes da pasta de origem.

Para reconstruir somente figuras e relatórios, use render_figures.py, write_reports.py e verify_final.py com os CSVs existentes. O lote completo principal/sensibilidades levou aproximadamente sete minutos nesta máquina; use results/execucao.json como tempo observado e results/estimativa.json como orçamento anterior ao lote. Redação, inspeção visual, renderização e verificação são etapas separadas.

## Interpretação e conservação

Todos os valores são descritivos de intensidade de imagem, lag e escala em pixels. A conclusão é restrita à categoria óptica e ao contexto selecionados. Não há calibração espacial física, identidade entre aquisições, validação externa ou inferência por réplicas independentes de eventos.

O .gitignore existente ignora src2. Nenhum arquivo desta entrega foi presumido versionado e nenhum commit foi criado. Para conservar a análise, preserve esta pasta junto às entradas e à base indicadas no manifesto.
