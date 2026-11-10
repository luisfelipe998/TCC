# Vazio AREAL flat-field × flip-rate (v11.3 B — SUPLEMENTAR)

Análise **suplementar** do par v11.3 (a principal é o bright-core, script A). Eixo X = fração de vazio **AREAL 2-D** (área de gás / área da seção, schepperle2024) sobre o quadro **flat-field por baseline de líquido** (percentil 20 por pixel; gás = excesso claro Fn>1,3). Eixo Y = taxa de troca temporal (flip-rate, v11.2 M4).

**Dois propósitos:** (i) PROVAR que a normalização remove o reflexo óptico W/M (o perfil radial colapsa ~100–300×, ver `flatfield_radial_proof_v11_3.png` e a coluna std radial); (ii) dar a métrica de ÁREA que o bright-core cordal não é.

> **LIMITE FÍSICO (não é bug):** mesmo normalizado, o limiar de área clara **inverte churn/annular** — o annular NÃO fica no topo. As bandas de interface ESCURAS do wavy annular (linhas ~20/60) persistem e subcontam o annular. A ordenação correta anular>churn só vem do bright-core (script A).

Padrões: dispersed bubbles, plug, slug, churn, wavy annular — 27 ensaios (26 na estatística).

## Mediana por padrão (saturados fora)

| padrão | n | vazio AREAL X | flip-rate Y (M4) |
|---|---|---|---|
| dispersed bubbles | 2 | 0.245 | 0.0233 |
| plug | 3 | 0.323 | 0.0177 |
| slug | 3 | 0.367 | 0.0174 |
| churn | 10 | 0.453 | 0.0300 |
| annular | 8 | 0.385 | 0.0177 |

**Ordem recuperada (baixo→alto):** dispersed bubbles 0.245 < plug 0.323 < slug 0.367 < annular 0.385 < churn 0.453. Correta até slug; o annular cai abaixo do churn (limite óptico documentado acima).

## Prova da remoção do reflexo W/M

Desvio-padrão radial do perfil de intensidade (col-mean por linha, média sobre quadros) na janela interior, BRUTO vs FLAT-FIELD. O W/M bruto tem std radial grande; após o flat-field colapsa ~100–300×. Ver `flatfield_radial_proof_v11_3.png`.

| scenario | padrão | std radial BRUTO | std radial FLAT | razão |
|---|---|---|---|---|
| m03q10_001 | wavy annular | 69.7 | 0.29 | 242× |
| m03q20_001 | wavy annular | 72.4 | 0.20 | 357× |
| m03q40_001 | wavy annular | 72.1 | 0.14 | 533× |
| m04q40_001 | wavy annular | 60.9 | 0.32 | 188× |
| m04q80_001 | wavy annular | 64.5 | 0.17 | 378× |
| m05q80_001 | wavy annular | 21.3 | 0.08 | 276× |
| m06q40_001 | wavy annular | 22.2 | 0.11 | 198× |
| m06q80_001 | wavy annular | 22.4 | 0.09 | 240× |
| m07q80_001 | wavy annular | 21.1 | 0.12 | 180× |
| m05q05_002 | dispersed_bubbles | 12.3 | 0.65 | 19× |
| m07q01_001 | dispersed bubbles | 12.1 | 0.01 | 1696× |
| m03q05_001 | churn | 67.2 | 0.58 | 115× |
| m04q10_001 | churn | 52.6 | 1.88 | 28× |
| m04q20_001 | churn | 59.0 | 0.94 | 63× |
| m05q10_001 | churn | 19.7 | 0.23 | 86× |
| m05q20_001 | churn | 20.4 | 0.20 | 100× |
| m05q40_001 | churn | 21.2 | 0.12 | 172× |
| m06q10_001 | churn | 22.4 | 0.28 | 80× |
| m06q20_001 | churn | 22.0 | 0.17 | 128× |
| m07q20_001 | churn | 19.4 | 0.31 | 62× |
| m07q40_001 | churn | 20.1 | 0.20 | 103× |
| m05q01_001 | plug | 20.5 | 0.56 | 37× |
| m06q01_001 | plug | 20.1 | 0.54 | 37× |
| m07q05_001 | plug | 16.6 | 0.60 | 28× |
| m04q01_001 | slug | 24.5 | 0.33 | 73× |
| m05q05_001 | slug | 17.0 | 0.62 | 27× |
| m06q05_001 | slug | 20.2 | 0.44 | 46× |

## Por ensaio (transparência)

| scenario | pattern | G | q | vazio areal | flip-rate | sat % | flag |
|---|---|---|---|---|---|---|---|
| m04q10_001 | churn | 510 | 10 | 0.578 | 0.0176 | 4 |  |
| m04q20_001 | churn | 510 | 20 | 0.565 | 0.0196 | 10 |  |
| m04q40_001 | wavy annular | 510 | 40 | 0.541 | 0.0159 | 8 |  |
| m03q05_001 | churn | 380 | 5 | 0.508 | 0.0192 | 12 |  |
| m05q05_002 | dispersed_bubbles | 640 | 5 | 0.490 | 0.0398 | 0 |  |
| m07q20_001 | churn | 890 | 20 | 0.482 | 0.0333 | 0 |  |
| m06q10_001 | churn | 760 | 10 | 0.460 | 0.0296 | 0 |  |
| m07q40_001 | churn | 890 | 40 | 0.446 | 0.0332 | 0 |  |
| m05q10_001 | churn | 640 | 10 | 0.443 | 0.0308 | 0 |  |
| m05q20_001 | churn | 640 | 20 | 0.441 | 0.0304 | 0 |  |
| m03q10_001 | wavy annular | 380 | 10 | 0.438 | 0.0154 | 12 |  |
| m06q20_001 | churn | 760 | 20 | 0.432 | 0.0323 | 0 |  |
| m05q40_001 | churn | 640 | 40 | 0.403 | 0.0286 | 0 |  |
| m04q80_001 | wavy annular | 510 | 80 | 0.393 | 0.0110 | 7 |  |
| m06q05_001 | slug | 760 | 5 | 0.388 | 0.0254 | 0 |  |
| m07q80_001 | wavy annular | 890 | 80 | 0.385 | 0.0221 | 0 |  |
| m06q40_001 | wavy annular | 760 | 40 | 0.385 | 0.0226 | 0 |  |
| m05q05_001 | slug | 640 | 5 | 0.367 | 0.0174 | 0 |  |
| m07q05_001 | plug | 890 | 5 | 0.357 | 0.0192 | 0 |  |
| m06q80_001 | wavy annular | 760 | 80 | 0.336 | 0.0203 | 0 |  |
| m06q01_001 | plug | 760 | 1 | 0.323 | 0.0177 | 0 |  |
| m03q20_001 | wavy annular | 380 | 20 | 0.315 | 0.0121 | 14 |  |
| m05q80_001 | wavy annular | 640 | 80 | 0.293 | 0.0194 | 0 |  |
| m05q01_001 | plug | 640 | 1 | 0.267 | 0.0124 | 0 |  |
| m04q01_001 | slug | 510 | 1 | 0.228 | 0.0156 | 0 |  |
| m03q40_001 | wavy annular | 380 | 40 | 0.139 | 0.0107 | 34 | saturado (sat=34%) |
| m07q01_001 | dispersed bubbles | 890 | 1 | 0.000 | 0.0067 | 0 |  |

## Ressalvas honestas

- **Proxy óptico RELATIVO não-calibrado** (sem px→mm, vista única) — não é a fração de área absoluta da seção.
- **annular ABAIXO de churn** (não fica no topo): limite óptico da vista única — as bandas de interface escuras do annular derrotam o limiar de área clara. A ordenação anular>churn só é recuperável pelo bright-core (script A).
- **dispersed bubbles instável (n=2):** os dois ensaios discordam muito; não entra em comparação estatística.
- **Dependência de iluminação residual** (corr ~0,3 com brilho/saturação); q-confound não testado por subconjunto q-match.
