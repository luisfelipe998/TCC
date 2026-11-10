# Fração de vazio × taxa de alternação por padrão (v11.1)

Revisão do v11. **Eixo Y = taxa de alternação gás↔líquido por 1000 px** (definição física de escoamento intermitente; costigan-whalley, jones-zuber), substituindo o `std(α_j)` do v11 que invertia a ordem anular×churn. A fração de vazio cordal média **α** permanece no eixo X (cordal, em pixels, sem calibração px→mm; comparação ENTRE padrões é legítima, o valor absoluto não é a fração de área da seção).

A taxa usa um gatilho de Schmitt sobre α_j: estado GÁS quando α_j>0.75, estado LÍQUIDO quando α_j<0.25 (banda morta mantém o estado). Limiar em α_j (fração), não em brilho ⇒ imune ao viés de saturação/q que contaminava o `std`.

Padrões: dispersed bubbles, plug, slug, churn, wavy annular — 26 ensaios (25 na estatística). Transições (plug-slug, slug-churn), smooth annular e synthetic desconsiderados.

## Estatística por padrão

| padrão | n | taxa mín | taxa máx | taxa média | taxa mediana | α mediana |
|---|---|---|---|---|---|---|
| dispersed bubbles | 1 | 0.0 | 0.0 | 0.0 | 0.0 | 0.017 |
| plug | 3 | 3.367 | 4.36 | 3.889 | 3.941 | 0.37 |
| slug | 3 | 2.169 | 5.061 | 3.266 | 2.568 | 0.431 |
| churn | 10 | 0.214 | 2.899 | 1.294 | 0.898 | 0.555 |
| annular | 8 | 0.0 | 1.201 | 0.429 | 0.206 | 0.524 |

**Ordem recuperada (taxa mediana):** plug 3.94 > slug 2.57 > churn 0.90 > annular 0.21 — anular é o MENOS intermitente, como esperado fisicamente.

## Sensibilidade à banda de histerese

Taxa mediana por padrão para quatro bandas [lo, hi]. A ORDEM plug>slug>churn>annular é estável para toda banda larga; só degrada ao estreitar a banda em direção a um corte na mediana.

| padrão | n | taxa 0.20/0.80 | taxa 0.25/0.75 | taxa 0.30/0.70 | taxa 0.35/0.65 |
|---|---|---|---|---|---|
| dispersed bubbles | 1 | 0.0 | 0.0 | 0.0 | 0.0 |
| plug | 3 | 3.126 | 3.941 | 3.941 | 4.757 |
| slug | 3 | 2.312 | 2.568 | 3.082 | 5.137 |
| churn | 10 | 0.3 | 0.898 | 2.616 | 5.327 |
| annular | 8 | 0.062 | 0.206 | 1.534 | 4.576 |

> **dispersed bubbles (n=1):** taxa = 0 — a corda central fica quase sempre líquida (α≈0,02; bolhas pequenas e dispersas nunca enchem a corda de gás), logo nunca atinge o estado gás. Fica no canto inferior-ESQUERDO (pouco vazio, pouca alternação), por motivo DISTINTO do anular (vazio alto e estável). Único ensaio — sem comparação.

> **Amostra pequena:** slug (n=3) e plug (n=3) têm poucos ensaios; mediana ≈ valor central. Ler com cautela.

> **m05q05 (G=640, q=5):** única condição (m, q) em DUAS classes — `_001`=slug, `_002`=plug (outlier confirmado pelo dono). Ambos mantidos.

> **Descartado da estatística (mantido na tabela):** m03q40_001 (saturado (sat=34%) — descartado da estatística). A saturação afeta o eixo X (α); o eixo Y (taxa) é imune (limiar em α_j), mas por coerência com o v11 o ensaio é excluído da estatística.

## Por ensaio (transparência)

Ordenado por taxa de alternação decrescente. `rate` = alternações gás↔líquido por 1000 px; `alpha_std(v11)` repetido para comparação com o descritor antigo.

| scenario | pattern | G [kg/m2 s] | q [kW/m2] | rate [/kpx] | alpha | alpha_std(v11) | sat % | N [px] | flag |
|---|---|---|---|---|---|---|---|---|---|
| m05q05_001 | slug | 640 | 5 | 5.061 | 0.386 | 0.249 | 0.0 | 6323 |  |
| m06q01_001 | plug | 760 | 1 | 4.36 | 0.37 | 0.195 | 0.0 | 7798 |  |
| m07q05_001 | plug | 890 | 5 | 3.941 | 0.344 | 0.225 | 0.0 | 7358 |  |
| m05q01_001 | plug | 640 | 1 | 3.367 | 0.384 | 0.239 | 0.0 | 5940 |  |
| m07q40_001 | churn | 890 | 40 | 2.899 | 0.573 | 0.206 | 0.0 | 10349 |  |
| m07q20_001 | churn | 890 | 20 | 2.786 | 0.528 | 0.178 | 0.0 | 8255 |  |
| m06q05_001 | slug | 760 | 5 | 2.568 | 0.431 | 0.199 | 0.0 | 7787 |  |
| m04q01_001 | slug | 510 | 1 | 2.169 | 0.448 | 0.198 | 0.0 | 5071 |  |
| m06q10_001 | churn | 760 | 10 | 1.714 | 0.532 | 0.182 | 0.0 | 7584 |  |
| m06q20_001 | churn | 760 | 20 | 1.593 | 0.559 | 0.197 | 0.0 | 10044 |  |
| m07q80_001 | wavy annular | 890 | 80 | 1.201 | 0.564 | 0.238 | 0.0 | 16658 |  |
| m03q05_001 | churn | 380 | 5 | 0.939 | 0.483 | 0.151 | 12.5 | 9586 |  |
| m06q80_001 | wavy annular | 760 | 80 | 0.911 | 0.537 | 0.232 | 0.0 | 18665 |  |
| m05q10_001 | churn | 640 | 10 | 0.857 | 0.582 | 0.168 | 0.0 | 7001 |  |
| m06q40_001 | wavy annular | 760 | 40 | 0.812 | 0.566 | 0.228 | 0.0 | 14772 |  |
| m05q40_001 | churn | 640 | 40 | 0.805 | 0.564 | 0.187 | 0.0 | 9937 |  |
| m05q20_001 | churn | 640 | 20 | 0.676 | 0.569 | 0.163 | 0.0 | 7391 |  |
| m03q40_001 | wavy annular | 380 | 40 | 0.482 | 0.275 | 0.077 | 34.3 | 20743 | saturado (sat=34%) — descartado da estatística |
| m04q10_001 | churn | 510 | 10 | 0.453 | 0.551 | 0.162 | 4.0 | 8829 |  |
| m04q20_001 | churn | 510 | 20 | 0.214 | 0.519 | 0.163 | 9.8 | 9344 |  |
| m03q10_001 | wavy annular | 380 | 10 | 0.211 | 0.462 | 0.148 | 12.4 | 14224 |  |
| m04q40_001 | wavy annular | 510 | 40 | 0.201 | 0.515 | 0.183 | 8.0 | 14960 |  |
| m05q80_001 | wavy annular | 640 | 80 | 0.057 | 0.533 | 0.207 | 0.0 | 17524 |  |
| m03q20_001 | wavy annular | 380 | 20 | 0.043 | 0.426 | 0.123 | 14.2 | 23079 |  |
| m04q80_001 | wavy annular | 510 | 80 | 0.0 | 0.473 | 0.159 | 7.2 | 22636 |  |
| m07q01_001 | dispersed bubbles | 890 | 1 | 0.0 | 0.017 | 0.067 | 0.0 | 3171 |  |
