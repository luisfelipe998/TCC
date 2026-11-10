# Âncoras e wavelets alternativas — v16.2_wco_eventos

A pasta contém a etapa 1 de plug/slug e a extensão executada para **agitado (churn)** e **anular ondulado (wavy annular)**, ambas baseadas em **v16_spatial2_18_71**. Os arquivos numéricos originais foram preservados; os novos ficam em subpastas extensao_churn_anular.

Comece pelo [resumo executivo com gráficos](EXECUTIVE_SUMMARY.md), [PDF integrado](RESUMO_EXECUTIVO.pdf), [passo a passo com limitações](PASSO_A_PASSO.md) e [explicação de âncoras, autocorrelação e limites de x/y](ANCORAS_E_AUTOCORRELACAO.md). O [resumo original](EXECUTIVE_SUMMARY_ETAPA1.md) e a [reprodução da etapa 1](README_ETAPA1.md) continuam disponíveis.

## Arquivos principais

- [Protocolo congelado da extensão](PROTOCOL_CHURN_ANULAR.md) e [parâmetros](config_ancoras_textura.json).
- [Catálogo completo](results/extensao_churn_anular/catalogo_ancoras.csv): candidatos, seleção, centro, largura óptica, proeminência, estrato, padrão, recorte, costuras, procedência e coordenadas alternativas.
- [Pares selecionados](results/extensao_churn_anular/pares_ancoras.csv) e [todos os pares candidatos](results/extensao_churn_anular/pares_candidatos.csv).
- [Resultados detalhados](RESULTS.md), [resultados ancorados](results/extensao_churn_anular/resultados_ancoras.csv), [controles fixos](results/extensao_churn_anular/resultados_fixos.csv) e [contrastes](results/extensao_churn_anular/contrastes_por_padrao.csv).
- [Pranchas das 19 aquisições](figures/extensao_churn_anular/ancoras), [figuras principais](figures/extensao_churn_anular) e [coeficientes integrais dos exemplos](maps/extensao_churn_anular).
- [Congelamento](results/extensao_churn_anular/congelamento.json), [validação numérica](results/extensao_churn_anular/validacao_numerica.json), [verificação da extensão](results/extensao_churn_anular/verificacao_final.json), [verificação integrada](results/verificacao_final.json) e [manifesto de hashes](results/manifesto.json).

## Reprodução dos relatórios e conferência da entrega

Use Python 3.13.2 em `/Users/i521862/.pyenv/versions/3.13.2/bin/python3`, com as dependências de requirements.txt. O manifesto registra versões, caminhos e hashes. A distribuição PyWavelets declara 1.9.0, enquanto o módulo declara 1.8.0; essa discrepância foi preservada. A Morlet principal usa NumPy/SciPy. Dados e base anteriores precisam estar nos caminhos relativos originais.

Para refazer as agregações, figuras e relatórios a partir dos coeficientes já calculados, use a sequência abaixo na raiz do projeto TCC. Ela conserva os lotes numéricos e o histórico da pausa:

```sh
export MPLCONFIGDIR=/tmp/wco-eventos-mpl
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/summarize_texture_analysis.py
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/render_texture_results.py
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/write_texture_reports.py
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/verify_final.py
```

O verificador integrado chama a verificação da extensão. Para somente conferir a entrega, basta o último comando: ele inclui reprodução temporária do catálogo e recálculos amostrais nas oito variantes. O manifesto inclui os arquivos de log; ao redirecionar a saída desse comando, use um arquivo fora desta pasta para que o log em gravação não invalide seu próprio hash. A geração do PDF integrado usa o utilitário local pdfunite; a verificação textual usa pdftotext, ambos disponíveis nesta execução.

O lote numérico integral está implementado em [catalog_texture_anchors.py](scripts/catalog_texture_anchors.py), [validate_texture_analysis.py](scripts/validate_texture_analysis.py) e [run_texture_analysis.py](scripts/run_texture_analysis.py), nessa ordem, seguido de [audit_texture_legacy.py](scripts/audit_texture_legacy.py). Uma execução integral deve usar uma cópia separada da pasta, mantendo os caminhos relativos da base e dos dados. Ela gera metadados e tempos próprios; não é a retomada histórica. Os verificadores finais desta entrega também exigem os registros e hashes daquela pausa/retomada, portanto uma nova execução precisa de registro e auditoria próprios, sem substituir nem flexibilizar a conferência histórica. A sequência acima é a reprodução direta suportada para a entrega atual.

As etapas concluídas somaram 36,4 minutos de cálculo. Esse total exclui a pausa e o tempo parcial dos diagnósticos interrompidos, que não foi cronometrado separadamente. O pico de memória do primeiro processo não foi salvo antes da pausa; o processo de retomada dos diagnósticos teve pico de 400 MiB. [benchmark.json](results/extensao_churn_anular/benchmark.json) preserva o orçamento anterior. O histórico da pausa e retomada está em [RETOMADA.md](RETOMADA.md); o script resume_texture_analysis.py concluiu somente os diagnósticos, mantendo os metadados das etapas anteriores. O cache é limitado a 16 autocorrelações transformadas, conferidas contra o núcleo original. A revisão visual e o congelamento são insumos históricos: reproduzir o cálculo não cria nova validação independente. Não regenere o congelamento para contornar uma divergência; uma alteração metodológica exige outra versão e declaração de que os resultados atuais já são conhecidos.

## Conservação e escopo

Os relatórios de etapa 1 são arquivos históricos, com o escopo anterior à extensão. Os scripts atuais de relatórios mantêm a integração quando os resultados da extensão existem. O núcleo matemático, as entradas e os resultados originais permanecem preservados. Nenhum lote da pasta de origem é executado por estes scripts.

Todos os valores são descritivos de intensidade óptica, lag e escala em pixels. Não há identidade física comprovada, calibração para mm/tempo, validação externa nem classificação validada. Pares compartilham aquisições e contextos. src2 continua ignorado pelo Git; a entrega é local, sem commit ou sincronização presumidos.
