# v16_spatial2 — análise espacial de Schuck

**Estado atual:** todas as figuras com nomes dos padrões em inglês e demais textos em PT-BR, sem COI nem máscara conservadora. As estatísticas dependentes dessas restrições foram recalculadas no domínio integral; as equações alternativas foram mantidas.

**Entrada:** panorama v11, média das linhas 20–59. **Escopo confirmado:** 26 registros para comparações; auditoria dos 31. **Método principal:** CWT Morlet analítica da correlação enviesada, alternative WPS, alternative WCS, alternative WCO e média da alternative WPS em todos os lags.

## Comece aqui

Veja primeiro a [conclusão executiva honesta](EXECUTIVE_SUMMARY.md), a [verificação das equações com o PDF](EQUATION_AUDIT.md) e o [experimento de pré-processamento](PREPROCESSING.md).

1. [Relatório: equações, fontes, justificativas e lacunas](REPORT.md)
2. [Resultados calculados e comparação com v5](RESULTS.md)
3. [Galeria dos 87 pares](figures/index.html)
4. [Perfis por padrão](figures/pattern_profiles.png)
5. [Diagnóstico de autocorrelação](figures/diagnostics/autocorrelation_estimators.png)
6. [Sensibilidade ao comprimento/posição da janela](figures/diagnostics/window_sensitivity.png)
7. [Controle das definições para múltiplas realizações](figures/diagnostics/ensemble_definition_control.png)

## Panorama, sinal e alternative WPS — exemplos predeterminados

- [churn: m03q05_001](figures/representatives/m03q05_001.png)
- [dispersed bubble: m07q01_001](figures/representatives/m07q01_001.png)
- [plug: m05q01_001](figures/representatives/m05q01_001.png)
- [slug: m04q01_001](figures/representatives/m04q01_001.png)
- [wavy annular: m03q10_001](figures/representatives/m03q10_001.png)

O padrão dispersed bubble tem somente um registro elegível. Sua alternative WPS é apresentada; não há alternative WCS/alternative WCO intrapadrão independente para esse padrão.

## Saídas atuais

- [EXECUTIVE_SUMMARY.md](EXECUTIVE_SUMMARY.md): conclusão atual, cinco exemplos de padrões e quatro comparações de alternative WCS/alternative WCO.
- [RESULTS.md](RESULTS.md): máximos integrais, perfis e diagnóstico nulo recalculado.
- [PREPROCESSING.md](PREPROCESSING.md): triagem, comparação com v5 e avaliação retirando uma aquisição, sem máscaras.
- [Perfis numéricos](results/profiles.csv), [resultados por registro](results/per_run_methods.csv) e [resultados por par](results/pair_methods.csv): campos integrais e no lag zero; as antigas colunas de COI/sobreposição foram retiradas.
- [Comparação numérica com v5](results/preprocessing_comparison_pairs.csv), [resumos dos candidatos](results/preprocessing_comparison_rankings.csv), [médias dos 87 pares](results/preprocessing_pooled_comparison.csv) e [conclusão estruturada](results/preprocessing_conclusion.json): versão atual com busca em todos os lags.
- [Diagnóstico nulo](results/null_diagnostic.csv) e [metadados](results/null_metadata.json): mesmo domínio completo para observações e sinais aleatórios.
- [Galeria](figures/index.html) e figuras em figures/: nomes dos padrões em inglês e demais textos em PT-BR, sem máscaras.
- [Verificação final](results/final_verification.json) e [manifesto](results/artifact_manifest.json): cobertura, equações e procedência.
- [Tempos do recálculo](results/recalculo_sem_mascaras.json), [migração das tabelas](results/migracao_tabelas_integrais.json) e [preservação dos coeficientes armazenados](results/migracao_cache_sem_mascaras.json): registro desta revisão.

## Organização e reprodução

- scripts/: núcleo matemático, auditoria, execução, controles e verificação final.
- results/: tabelas, parâmetros, testes, logs e manifestos.
- figures/: 87 mapas de pares, 31 comparações individuais e diagnósticos.
- [results/historico_com_mascaras.zip](results/historico_com_mascaras.zip): documentos e saídas anteriores à retirada das máscaras; somente para rastreabilidade.
- literature/: extrações dos PDFs, páginas escaneadas de Schuck e artigo adicional.
- intermediate/: sinais conferidos e matrizes reduzidas apenas para renderização. As estatísticas são calculadas na resolução integral; os mapas armazenados têm aproximadamente 1800 posições de lag. Não use esses mapas reduzidos para recalcular máximos integrais.

Para regenerar somente as figuras e os relatórios a partir dos resultados existentes, sem repetir a triagem ou o diagnóstico nulo:

~~~sh
MPLCONFIGDIR=/tmp/v16-mpl python3 src2/v16_spatial2/scripts/regenerate_figures.py
MPLCONFIGDIR=/tmp/v16-mpl python3 src2/v16_spatial2/scripts/verify_final.py
~~~

Os nomes exibidos dos padrões são dispersed bubble, plug, slug, churn e wavy annular; os demais textos estão em PT-BR. Os identificadores internos das tabelas foram preservados. A [verificação dos rótulos](results/revisao_nomes_padroes.json) registra a regeneração e a preservação dos dados numéricos.

Dependências testadas em requirements.txt. Python 3.13.2 foi usado. Para reproduzir também os cálculos completos, execute a partir da raiz do projeto:

~~~sh
export MPLCONFIGDIR=/tmp/v16-mpl
python3 src2/v16_spatial2/scripts/audit_inputs.py
python3 src2/v16_spatial2/scripts/validate_wavelets.py
python3 src2/v16_spatial2/scripts/run_analysis.py --stage all
python3 src2/v16_spatial2/scripts/verify_legacy_baseline.py
python3 src2/v16_spatial2/scripts/additional_diagnostics.py
python3 src2/v16_spatial2/scripts/precision_and_windows.py
python3 src2/v16_spatial2/scripts/null_diagnostic.py --surrogates 199
python3 src2/v16_spatial2/scripts/preprocessing_experiment.py --scales 64 --tag screen
python3 src2/v16_spatial2/scripts/evaluate_preprocessing.py --tag screen
python3 src2/v16_spatial2/scripts/preprocessing_experiment.py --scales 200 --variants baseline central_band20 core_row_median --tag confirm
python3 src2/v16_spatial2/scripts/evaluate_preprocessing.py --tag confirm
python3 src2/v16_spatial2/scripts/compare_v5.py
python3 src2/v16_spatial2/scripts/report_preprocessing.py
python3 src2/v16_spatial2/scripts/finish_outputs.py
python3 src2/v16_spatial2/scripts/verify_final.py
~~~

Os tempos efetivos desta atualização estão nos logs e metadados em results. O diagnóstico nulo tem duração registrada em null_metadata.json; as etapas de pré-processamento, em preprocessing_*_metadata.json. Tempos históricos não são estimativas de uma nova execução e variam com carga e hardware. A extração/OCR dos PDFs é material de pesquisa preservado, não uma dependência para recalcular as wavelets.

O .gitignore existente ignora src2; os arquivos estão disponíveis no workspace, mas não foram adicionados ao Git. Nenhum commit foi feito.
