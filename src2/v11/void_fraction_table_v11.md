# Fração de vazio por padrão de escoamento (v11)

Fração de vazio **cordal média no tempo** α por limiar de Otsu no núcleo do panorama (linhas 20–59, paredes excluídas; gás = `valor <= t`). Para cada coluna (instante) α_j = fração cordal de gás; α = média temporal de α_j. Por ser **cordal** (uma linha, não a seção), α carrega o viés cordal>areal conhecido para núcleo gasoso central (Hewitt; Kataoka–Ishii–Serizawa); em pixels, sem calibração px→mm. A comparação ENTRE padrões é legítima; o valor absoluto **não** é a fração de área da seção transversal.

Padrões: plug, slug, churn, wavy annular — 26 ensaios (25 na estatística). Transições (plug-slug, slug-churn), smooth annular, dispersed bubbles e synthetic foram desconsiderados.

## Estatística por padrão

| padrão | n | α mín | α máx | α média | α mediana |
|---|---|---|---|---|---|
| plug | 4 | 0.344 | 0.402 | 0.375 | 0.377 |
| slug | 3 | 0.386 | 0.448 | 0.422 | 0.431 |
| churn | 10 | 0.483 | 0.582 | 0.546 | 0.555 |
| annular | 8 | 0.426 | 0.566 | 0.509 | 0.524 |

> **Amostra pequena:** slug (n=3) e plug (n=4) têm poucos ensaios — a mediana é praticamente o valor central e `α mín`/`α máx` são os próprios extremos amostrais, não estimativas. Ler com cautela; sem teste estatístico.

> **m05q05 (G=640, q=5):** única condição (m, q) presente em DUAS classes — `m05q05_001` = slug e `m05q05_002` = plug, conforme os marcadores `__*.txt` (fonte de verdade; `_002` é outlier de plug confirmado pelo dono). Ambos mantidos.

> **Descartados da estatística (mantidos na tabela por ensaio):** m03q40_001 (saturado (sat=34%) — descartado da estatística). Saturação do panorama (núcleo recortado em 255) impede o Otsu de separar gás/líquido no lado claro — α fica artificialmente baixo. Correções de de-trend/saturação/limiar local foram testadas e descartadas (perda de informação por recorte).

## Fração de vazio por ensaio (transparência)

Ordenado por α crescente. `sat %` = fração do núcleo saturada em 255; `flag` sinaliza ensaio fora da estatística.

| scenario | pattern | G [kg/m2 s] | q [kW/m2] | alpha | alpha_std | sat % | N [px] | flag |
|---|---|---|---|---|---|---|---|---|
| m03q40_001 | wavy annular | 380 | 40 | 0.275 | 0.077 | 34.3 | 20743 | saturado (sat=34%) — descartado da estatística |
| m07q05_001 | plug | 890 | 5 | 0.344 | 0.225 | 0.0 | 7358 |  |
| m06q01_001 | plug | 760 | 1 | 0.37 | 0.195 | 0.0 | 7798 |  |
| m05q01_001 | plug | 640 | 1 | 0.384 | 0.239 | 0.0 | 5940 |  |
| m05q05_001 | slug | 640 | 5 | 0.386 | 0.249 | 0.0 | 6323 |  |
| m05q05_002 | plug | 640 | 5 | 0.402 | 0.291 | 0.0 | 4510 |  |
| m03q20_001 | wavy annular | 380 | 20 | 0.426 | 0.123 | 14.2 | 23079 |  |
| m06q05_001 | slug | 760 | 5 | 0.431 | 0.199 | 0.0 | 7787 |  |
| m04q01_001 | slug | 510 | 1 | 0.448 | 0.198 | 0.0 | 5071 |  |
| m03q10_001 | wavy annular | 380 | 10 | 0.462 | 0.148 | 12.4 | 14224 |  |
| m04q80_001 | wavy annular | 510 | 80 | 0.473 | 0.159 | 7.2 | 22636 |  |
| m03q05_001 | churn | 380 | 5 | 0.483 | 0.151 | 12.5 | 9586 |  |
| m04q40_001 | wavy annular | 510 | 40 | 0.515 | 0.183 | 8.0 | 14960 |  |
| m04q20_001 | churn | 510 | 20 | 0.519 | 0.163 | 9.8 | 9344 |  |
| m07q20_001 | churn | 890 | 20 | 0.528 | 0.178 | 0.0 | 8255 |  |
| m06q10_001 | churn | 760 | 10 | 0.532 | 0.182 | 0.0 | 7584 |  |
| m05q80_001 | wavy annular | 640 | 80 | 0.533 | 0.207 | 0.0 | 17524 |  |
| m06q80_001 | wavy annular | 760 | 80 | 0.537 | 0.232 | 0.0 | 18665 |  |
| m04q10_001 | churn | 510 | 10 | 0.551 | 0.162 | 4.0 | 8829 |  |
| m06q20_001 | churn | 760 | 20 | 0.559 | 0.197 | 0.0 | 10044 |  |
| m05q40_001 | churn | 640 | 40 | 0.564 | 0.187 | 0.0 | 9937 |  |
| m07q80_001 | wavy annular | 890 | 80 | 0.564 | 0.238 | 0.0 | 16658 |  |
| m06q40_001 | wavy annular | 760 | 40 | 0.566 | 0.228 | 0.0 | 14772 |  |
| m05q20_001 | churn | 640 | 20 | 0.569 | 0.163 | 0.0 | 7391 |  |
| m07q40_001 | churn | 890 | 40 | 0.573 | 0.206 | 0.0 | 10349 |  |
| m05q10_001 | churn | 640 | 10 | 0.582 | 0.168 | 0.0 | 7001 |  |
