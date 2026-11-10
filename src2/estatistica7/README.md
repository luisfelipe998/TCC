# Exploração estatística da WCO

**Resultado: útil como diagnóstico metodológico; não foi obtido agrupamento estável dos padrões.** Sinais independentes também produziram máximos próximos de 1. Os resumos testados não mantiveram compactação de todas as classes nas posições avaliadas, inclusive plug/slug.

- [Análises selecionadas para o TCC — gráficos, explicação e equações](ANALISES_RELEVANTES_TCC.md)
- [Resumo executivo em PDF — 4 páginas](RESUMO_EXECUTIVO.pdf)
- [Resumo em português com gráficos](RESUMO_EXECUTIVO.md)
- [Método, fórmulas, limites e reprodução](METODOLOGIA.md)
- [Todas as métricas de separação](results/metrics_all.csv)
- [Referência estatística de fases](results/phase_null_summary.csv)
- [Controles com sinais independentes](results/independent_noise_summary.csv)
- [Auditoria numérica](validation/final.json)

Foram utilizados 26 ensaios da `v16_spatial2_18_71`, os mesmos elegíveis nas análises anteriores. A WCO foi recalculada com comprimento comum dentro de cada comparação, em nove condições de comprimento/posição. `estatistica5`, `estatistica6` e as fontes v16 permaneceram preservadas.

## Figuras

Cada figura está disponível em PNG e PDF na pasta [figures](figures).

| Figura | O que permite avaliar |
| --- | --- |
| [01 — Máximos observados](figures/01_maximos_observados.pdf) | Sobreposição intra/inter, concentração na maior escala e efeito de permitir lags diferentes por escala. |
| [02 — Controles independentes](figures/02_controles_independentes.pdf) | Máximos altos entre sinais independentes e sensibilidade ao comprimento. |
| [03 — Triagem das 22 representações](figures/03_triagem_wco.pdf) | Silhueta da classe menos favorecida; evita esconder uma classe ruim na média. |
| [04 — Estabilidade](figures/04_estabilidade.pdf) | Mesmos ensaios, começo/centro/fim, incluindo plug/slug. |
| [05 — Contrastes por escala](figures/05_contrastes_por_escala.pdf) | Perfis completos e influência da retirada de um ensaio, sem selecionar escalas favoráveis. |
| [06 — Modelo de fases](figures/06_modelo_de_fases.pdf) | Seis contrastes de coorte, com correção de Holm. |
| [07 — Autocoerência](figures/07_autocoerencia.pdf) | Descritores próprios dos ensaios: área e meia largura por escala. |
| [08 — Mapas de exemplo](figures/08_mapas_exemplo.pdf) | Seis pares definidos antes da inspeção, com escala de cor comum. |
| [09 — Condições operacionais](figures/09_condicoes_operacionais.pdf) | Cobertura de q e associação entre padrão e comprimento original. |

## Como reproduzir

No diretório do TCC, com Python e as versões registradas no [manifesto inicial](results/manifest_started.json):

```sh
export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
python3 src2/estatistica7/scripts/maps.py
python3 src2/estatistica7/scripts/analyze.py
python3 src2/estatistica7/scripts/phase_null.py
python3 src2/estatistica7/scripts/diagnostics.py
python3 src2/estatistica7/scripts/validate.py
python3 src2/estatistica7/scripts/report.py
python3 src2/estatistica7/scripts/finalize.py
```

Para apenas reconstruir figuras e relatórios, execute `report.py` e `finalize.py`. A configuração arquivada fixa 199 réplicas de fases e as sementes. O `benchmark.py` foi executado antes dos mapas para conferir o motor e definir o custo; não deve ser reexecutado para reproduzir esta execução, pois escreve uma nova configuração de tempo e um novo manifesto inicial.

Os perfis NPZ conservam os resumos por escala calculados no domínio integral. Apenas os mapas de exemplo foram reduzidos a 301 lags para exibição. Os mapas integrais podem ser regenerados pelo código; não é necessário armazenar todas as matrizes volumosas para reproduzir os resumos.

O documento de análises selecionadas reúne resultados de estatistica5/6/7. Sua figura exclusiva do pico WPS pode ser regenerada com `python3 src2/estatistica7/scripts/figura_pico_selecionada.py`, usando somente a tabela arquivada da estatistica6. As outras quatro figuras são referenciadas diretamente nos arquivos originais.
