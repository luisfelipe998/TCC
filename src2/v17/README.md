# Auditoria v17 das análises de padrões de escoamento

Leia primeiro o [relatório consolidado](findings/RELATORIO_CONSOLIDADO.md). Ele reúne as três auditorias independentes, as contraprovas e os caminhos recomendados sem aprendizado de máquina.

- **[findings/](findings/)**: interpretações, inconsistências, limites, evidências e recomendações. Inclui um relatório por pasta, confronto com o TCC e três contraprovas cruzadas.
- **[resultados/](resultados/)**: CSV, JSON, NPZ, gráficos PNG/PDF e registros de procedência. Resultados experimentais separados dos findings.
- **[scripts/](scripts/)**: cálculos reproduzíveis, com saídas restritas a v17.
- **[logs/](logs/)**: registros de execução.
- **[PROTOCOLO.txt](PROTOCOLO.txt)**: escopo, parâmetros fixados e adendo exploratório posterior aos primeiros resultados.

As fontes são `src2/estatistica8`, `src2/v14.1_new_voidfraction`, `src2/v16_spatial2_18_71`, suas entradas/dependências e o texto atual em `TCC/`. TCC_1 não é a versão confrontada. Os resultados principais reproduziram o código e as tabelas do TCC; os achados mais relevantes são metodológicos. O novo descritor de largura Welch mostrou contraste promissor churn/anular, ainda exploratório e sensível ao estimador. Não houve treinamento de classificador; o resultado de vizinhos já existente foi somente auditado.

## Navegação

| Pergunta | Arquivo |
| --- | --- |
| O que foi confirmado e o que precisa mudar na interpretação? | [Síntese](findings/RELATORIO_CONSOLIDADO.md) |
| Os 95% e a largura WPS de estatistica8 se sustentam? | [Auditoria estatistica8](findings/auditoria_estatistica8.md) |
| A máscara mede vazio físico? Flips são independentes da ocupação? | [Auditoria v14.1](findings/auditoria_voidfraction.md) |
| As fórmulas wavelet estão corretas? Como comprimento e janela afetam as métricas? | [Auditoria v16](findings/auditoria_wavelets.md) |
| O texto corresponde aos dados e qual alternativa foi testada? | [Confronto e exploração](findings/confronto_tcc_e_exploracao.md) |
| A separação Welch resistiu a uma segunda implementação? | [Contraprova Welch](findings/contraprova_exploracao.md) |
| O singleton realmente pesa tanto? | [Contraprova da composição](findings/contraprova_composicao.md) |
| A escolha de âncoras usa regras diferentes entre churn e anular? | [Contraprova das âncoras](findings/contraprova_estatistica8.md) |
| Os arquivos originais foram preservados? | [Verificação final](resultados/validacao_entrega.json) |

## Reprodução

Execute da raiz do repositório, com o ambiente científico já disponível. Python 3.13.2, NumPy 2.3.4, SciPy 1.18.0, pandas 2.3.3 e Matplotlib 3.10.7 foram usados; Pillow, PyWavelets e openpyxl são dependências de leituras/núcleos herdados. As versões também estão nos JSONs de execução. Os scripts pressupõem que as fontes locais permaneçam nos caminhos originais.

```sh
export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
python3 -B src2/v17/scripts/audit_estatistica8.py
python3 -B src2/v17/scripts/audit_voidfraction.py
python3 -B src2/v17/scripts/plot_voidfraction_audit.py
python3 -B src2/v17/scripts/audit_wavelets.py
python3 -B src2/v17/scripts/plot_wavelet_audit.py
python3 -B src2/v17/scripts/exploracao_espectral.py
python3 -B src2/v17/scripts/audit_comparacao_tcc.py
python3 -B src2/v17/scripts/contraprova_composicao.py
python3 -B src2/v17/scripts/contraprova_exploracao.py
python3 -B src2/v17/scripts/counterreview_estatistica8.py
python3 -B src2/v17/scripts/verificar_entrega.py
```

As etapas podem regravar seus próprios resultados em v17; os relatórios interpretativos não são reescritos automaticamente. Não execute os `main` ou `run_all` das pastas originais para reproduzir esta entrega: eles possuem destinos próprios e alguns incluem classificadores. As contraprovas de v17 não substituem novas aquisições experimentais.

O [manifesto inicial](resultados/fontes_antes.json) preserva hashes de 2.710 arquivos do TCC e das três pastas auditadas. O [estado dos outros arquivos](resultados/estado_outros_arquivos_antes.json) registra tamanho/data de 46.379 arquivos fora de v17 em src, src2, TCC e TCC_1. Esses registros permitem detectar mudanças posteriores; não são gerados novamente ao reproduzir os cálculos.

`src2` já é ignorado pelo `.gitignore` do projeto. A v17 está salva no disco, sem alterar essa configuração e sem criar commits. A validação final registra eventuais mudanças externas de estado do Git separadamente da preservação das fontes auditadas.

Durante a execução, houve alterações concorrentes no TCC e recompilação por atividade externa a esta tarefa. As três pastas de análise permaneceram idênticas ao início. O TCC final foi reconferido; veja a [reconciliação](findings/reconciliacao_texto_atual.md) e o [primeiro registro que detectou a alteração](resultados/validacao_entrega_primeira.json). Por isso a validação distingue preservação das análises, verificações computacionais e mudança do snapshot do TCC. As citações de linha dos relatórios iniciais devem ser lidas com a reconciliação final.
