# v12.1 — Fração de vazio α, velocidades reais de fase (uv, ul) e escorregamento S

Extensão do v12 com o modelo de **escorregamento** (Rouhani–Axelsson). x, Jg e Jl são **idênticos ao v12** (velocidades superficiais são invariantes ao escorregamento); o que o v12.1 acrescenta é α e as velocidades **reais** de fase — as grandezas que a literatura tabula e que de fato mudam com o escorregamento.

## Propriedades do R600a (isobutano) a Tsat = 27 °C
- Fonte: **CoolProp 8.0.0** (`PropsSI`, fluido `R600a`)
- ρ_l = **548.13 kg/m³**  ·  ρ_v = **9.650 kg/m³**
- i_lv = **326.71 kJ/kg**  ·  σ = **9.782 mN/m**  ·  P_sat = 3.716 bar
- Geometria: D = 1,0 mm · L_aquecido = 0,265 m · A = 7.854e-07 m²

## Modelo acrescentado (sobre o v12)
```
# x, Jg=jv, Jl=jl  : idênticos ao v12 (não mudam com escorregamento)
C0    = 1 + 0.12·(1 − x)
drift = 1.18·(1 − x)·[g·σ·(ρl − ρv)]^0.25 / (G·ρl^0.5)
α     = (x/ρv) / ( C0·(x/ρv + (1−x)/ρl) + drift )   # Rouhani–Axelsson (1970)
uv    = jv / α            ul = jl / (1 − α)          # velocidades reais
S     = (x·ρl·(1−α)) / (α·ρv·(1−x))                  # escorregamento = uv/ul
```
> **Grounding:** α — dias2025knn Eq. (5); S — deoliveira2021 Eq. (17). uv=jv/α e ul=jl/(1−α) decorrem da definição de α (implícito em deoliveira2021 via Eq. 17), não citação literal.

## Por ensaio
| Ensaio | Padrão | G [kg/m²s] | q″ [kW/m²] | x [-] | Jg=jv [m/s] | Jl=jl [m/s] | α [-] | uv [m/s] | ul [m/s] | S [-] |
|--------|--------|-----------:|-----------:|------:|------------:|-----------:|------:|---------:|---------:|------:|
| m07q01 | bolha dispersa | 890 | 1 | 0.0036 | 0.336 | 1.6178 | 0.1448 | 2.322 | 1.8917 | 1.228 |
| m05q01 | plug | 640 | 1 | 0.0051 | 0.336 | 1.1617 | 0.1856 | 1.811 | 1.4265 | 1.270 |
| m06q01 | plug | 760 | 1 | 0.0043 | 0.336 | 1.3806 | 0.1635 | 2.057 | 1.6504 | 1.246 |
| m07q05 | plug | 890 | 5 | 0.0182 | 1.681 | 1.5941 | 0.4431 | 3.794 | 2.8626 | 1.325 |
| m04q01 | slug | 510 | 1 | 0.0064 | 0.336 | 0.9245 | 0.2176 | 1.545 | 1.1816 | 1.308 |
| m05q05 | slug | 640 | 5 | 0.0253 | 1.681 | 1.1380 | 0.5125 | 3.281 | 2.3341 | 1.405 |
| m05q05_002 | slug | 640 | 5 | 0.0253 | 1.681 | 1.1380 | 0.5125 | 3.281 | 2.3341 | 1.405 |
| m06q05 | slug | 760 | 5 | 0.0213 | 1.681 | 1.3569 | 0.4766 | 3.527 | 2.5927 | 1.360 |
| m03q05 | churn | 380 | 5 | 0.0427 | 1.681 | 0.6637 | 0.6128 | 2.744 | 1.7139 | 1.601 |
| m04q10 | churn | 510 | 10 | 0.0636 | 3.362 | 0.8712 | 0.6953 | 4.836 | 2.8593 | 1.691 |
| m04q20 | churn | 510 | 20 | 0.1272 | 6.725 | 0.8121 | 0.7964 | 8.444 | 3.9881 | 2.117 |
| m05q10 | churn | 640 | 10 | 0.0507 | 3.362 | 1.1084 | 0.6582 | 5.108 | 3.2429 | 1.575 |
| m05q20 | churn | 640 | 20 | 0.1014 | 6.725 | 1.0492 | 0.7700 | 8.733 | 4.5613 | 1.915 |
| m05q40 | churn | 640 | 40 | 0.2028 | 13.449 | 0.9308 | 0.8478 | 15.863 | 6.1163 | 2.594 |
| m06q10 | churn | 760 | 10 | 0.0427 | 3.362 | 1.3273 | 0.6276 | 5.358 | 3.5639 | 1.503 |
| m06q20 | churn | 760 | 20 | 0.0854 | 6.725 | 1.2681 | 0.7477 | 8.993 | 5.0265 | 1.789 |
| m07q20 | churn | 890 | 20 | 0.0729 | 6.725 | 1.5053 | 0.7254 | 9.271 | 5.4809 | 1.691 |
| m07q40 | churn | 890 | 40 | 0.1458 | 13.449 | 1.3869 | 0.8165 | 16.472 | 7.5571 | 2.180 |
| m03q10 | anular ondulado | 380 | 10 | 0.0854 | 3.362 | 0.6341 | 0.7376 | 4.559 | 2.4162 | 1.887 |
| m03q20 | anular ondulado | 380 | 20 | 0.1708 | 6.725 | 0.5749 | 0.8263 | 8.138 | 3.3102 | 2.458 |
| m03q40 | anular ondulado | 380 | 40 | 0.3415 | 13.449 | 0.4565 | 0.8911 | 15.093 | 4.1904 | 3.602 |
| m03q80 | anular ondulado | 380 | 80 | 0.6831 | 26.898 | 0.2197 | 0.9541 | 28.192 | 4.7874 | 5.889 |
| m04q40 | anular ondulado | 510 | 40 | 0.2545 | 13.449 | 0.6937 | 0.8672 | 15.509 | 5.2231 | 2.969 |
| m04q80 | anular ondulado | 510 | 80 | 0.5089 | 26.898 | 0.4569 | 0.9265 | 29.033 | 6.2126 | 4.673 |
| m05q80 | anular ondulado | 640 | 80 | 0.4056 | 26.898 | 0.6941 | 0.9075 | 29.641 | 7.5011 | 3.952 |
| m06q40 | anular ondulado | 760 | 40 | 0.1708 | 13.449 | 1.1498 | 0.8321 | 16.164 | 6.8462 | 2.361 |
| m06q80 | anular ondulado | 760 | 80 | 0.3415 | 26.898 | 0.9130 | 0.8937 | 30.098 | 8.5885 | 3.504 |
| m07q80 | anular ondulado | 890 | 80 | 0.2916 | 26.898 | 1.1502 | 0.8811 | 30.528 | 9.6726 | 3.156 |

## Resumo por padrão (min / máx / média / mediana)
| Padrão | n | α min | α máx | α média | α med | S min | S máx | S média | S med | uv média | ul média |
|--------|--:|------:|------:|--------:|------:|------:|------:|--------:|------:|---------:|---------:|
| bolha dispersa | 1 | 0.145 | 0.145 | 0.145 | 0.145 | 1.228 | 1.228 | 1.228 | 1.228 | 2.322 | 1.8917 |
| plug | 3 | 0.163 | 0.443 | 0.264 | 0.186 | 1.246 | 1.325 | 1.280 | 1.270 | 2.554 | 1.9799 |
| slug | 4 | 0.218 | 0.512 | 0.430 | 0.495 | 1.308 | 1.405 | 1.370 | 1.383 | 2.908 | 2.1106 |
| churn | 10 | 0.613 | 0.848 | 0.730 | 0.737 | 1.503 | 2.594 | 1.866 | 1.740 | 8.582 | 4.4110 |
| anular ondulado | 10 | 0.738 | 0.954 | 0.872 | 0.886 | 1.887 | 5.889 | 3.445 | 3.330 | 20.695 | 5.8748 |

## Validação de α por padrão × literatura local
| Padrão | α v12.1 (min–máx) | Literatura | Citação | Veredito |
|--------|-------------------|-----------|---------|----------|
| plug | 0.163–0.443 | α 0,042–0,29 (faixa) | deoliveira2021 §4.2 | ordem de grandeza compatível (plug = baixo vazio) |
| slug | 0.218–0.512 | α 0,578–0,714 (faixa, G=520–1220); 0,73 (ponto, G=516) | Dias (2024) Figs. 22d/23b/24a; dias2025knn Tab. 2 | v12.1 **abaixo** da faixa (x_in=0 subestima; ver nota) |
| churn | 0.613–0.848 | α 0,61–0,85 (faixa) | deoliveira2021 §4.4 | faixa sobreposta; v12.1 alcança o topo da faixa |
| anular ondulado | 0.738–0.954 | α 0,86–0,99 (faixa) | deoliveira2021 §4.5 | **bate** (dentro da faixa) |

> **Nota (slug):** de Oliveira (2021) **não publica faixa de α para slug** (confirmado lendo a §4.3 inteira — só dá descritores de caos/Hurst/LLE; faixa de α existe lá só para plug §4.2, churn §4.4 e anular §4.5). A faixa de slug mostrada vem da **dissertação Dias (2024)**, mesmo arranjo LETEF (R-600a, mini-canal 1 mm), três pontos **confirmados visualmente** nas páginas renderizadas (Figs. 22d/23b/24a: α=0,7023/0,7135/0,5782 em G=520/740/1220; o PDF tem OCR corrompido — ver `assets/paper-log.md`). O ponto único α=0,73 de `dias2025knn` Tab. 2 (G=516) **corrobora** a faixa (marcador × na figura). O v12.1 dá α_slug 0.22–0.51, **abaixo** da faixa 0,58–0,71 de Dias (2024) — mesmo viés do plug: o x de equilíbrio (x_in=0) subestima o vazio nos padrões de baixo-a-médio título. A faixa da literatura é mostrada como **referência de onde o slug real vive**, não como concordância de magnitude do v12.1.

## Limitação herdada do v12 (honesta)
α, uv, ul e S são calculados sobre o **x de equilíbrio** (x_in = 0). Não há potência do pré-aquecedor nem vazão mássica medida nos dados brutos (`src/data_processed/<run>/dados.xlsx` contém só `source_file`, `Voltage_0`, `Voltage_1`, `Voltage_2`) para reconstruir o x de entrada real. Logo estas grandezas herdam a mesma natureza de **estimativa de modelo** do v12 — declarado, não escondido. Ver `src2/v12.1/ANALISE_ADVERSARIAL_v12.md`.

## Verificações internas executadas
- herdado v12: x_out(m03q80)=0.683; ensaios={'bubble': 1, 'plug': 3, 'slug': 4, 'churn': 10, 'wavy annular': 10} (28)  [ok]
- 0 < α < 1 em todos os 28 ensaios  [ok]
- α monotônico ↑ com x (G fixo)  [ok]
- S ≥ 1 em todos (28 finitos; máx=5.89)  [ok]
- uv ≥ jv e ul ≥ jl em todos  [ok]
- α anular ∈ [0.738; 0.954], mediana=0.886 (lit. deoliveira2021 0,86–0,99; piso baixo = x_in=0 herdado)  [ok]
