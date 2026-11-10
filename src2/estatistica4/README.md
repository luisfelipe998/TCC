# estatistica4 — Pearson nos dados 18–71

**Resultado:** separação fraca entre churn e wavy annular. ANOSIM R = 0.0837; p exato por q = 0.085; p ajustado por Holm = 0.090. Não houve evidência convincente de separação robusta pelo descritor principal.

- [Relatório completo](RELATORIO.md)
- [Métodos, fórmulas e plano anterior ao cálculo](METODOLOGIA.md)
- [Resumo numérico](results/summary.json)
- [Verificação das fórmulas](results/formula_validation.json) e [das fontes](results/source_validation.json)
- [Matriz de Pearson](results/pearson_r_matrix.csv), [distâncias por par](results/pair_distances.csv) e [200 permutações exatas](results/exact_q_permutations.csv)
- [Sensibilidades](results/robustness.csv), [influência](results/leave_one_run_out.csv), [dispersão](results/dispersion_summary.csv) e [26 ensaios](results/all26_descriptive_summary.csv)
- [Manifesto, versões e hashes](results/manifest.json)

## Reproduzir

Na raiz do repositório:

    PYTHONDONTWRITEBYTECODE=1 python3 src2/estatistica4/run_analysis.py

O comando verifica as fórmulas antes dos dados, confere as fontes e gera descritores, tabelas, figuras e relatório. Escreve em estatistica4 e usa /tmp/estatistica4-matplotlib para o cache gráfico. As análises anteriores não são alteradas. As dependências estão em requirements.txt e as versões utilizadas, no manifesto.

Somente controles de fórmulas:

    PYTHONDONTWRITEBYTECODE=1 python3 src2/estatistica4/validate_analysis.py

O segundo comando atualiza apenas o JSON de controles. Para renovar o manifesto após mudanças, executar novamente o comando completo. As permutações experimentais são enumeradas integralmente; somente os controles sintéticos usam semente fixa. Campos vazios em comparações não aplicáveis ou dispersão não estimável não representam zero.
