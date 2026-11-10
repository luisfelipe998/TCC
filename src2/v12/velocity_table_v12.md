# v12 — Velocidades superficiais (Jg, Jl) por padrão de escoamento

Mapa de padrões de escoamento no plano velocidade superficial de gás (Jg) × velocidade superficial de líquido (Jl), com clusters por padrão. Padrões incluídos: **bolha dispersa, plug, slug, churn e anular ondulado**.

## Propriedades do R600a (isobutano) a Tsat = 27 °C
- Fonte: **CoolProp 8.0.0** (`PropsSI`, fluido `R600a`)
- ρ_l = **548.13 kg/m³**  ·  ρ_v = **9.650 kg/m³**
- i_lv = **326.71 kJ/kg**  ·  P_sat = 3.716 bar
- Geometria: D = 1,0 mm · L_aquecido = 0,265 m · A = 7.854e-07 m²

## Modelo
```
x_out = 4·q″·L / (G·D·i_lv)     (título de equilíbrio; entrada líquido saturado)
Jg    = G·x / ρ_v              Jl = G·(1−x) / ρ_l
```

## Por ensaio
| Ensaio | Padrão | G [kg/m²s] | q″ [kW/m²] | x [-] | Jg [m/s] | Jl [m/s] |
|--------|--------|-----------:|-----------:|------:|---------:|---------:|
| m07q01 | bolha dispersa | 890 | 1 | 0.0036 | 0.336 | 1.6178 |
| m05q01 | plug | 640 | 1 | 0.0051 | 0.336 | 1.1617 |
| m06q01 | plug | 760 | 1 | 0.0043 | 0.336 | 1.3806 |
| m07q05 | plug | 890 | 5 | 0.0182 | 1.681 | 1.5941 |
| m04q01 | slug | 510 | 1 | 0.0064 | 0.336 | 0.9245 |
| m05q05 | slug | 640 | 5 | 0.0253 | 1.681 | 1.1380 |
| m05q05_002 | slug | 640 | 5 | 0.0253 | 1.681 | 1.1380 |
| m06q05 | slug | 760 | 5 | 0.0213 | 1.681 | 1.3569 |
| m03q05 | churn | 380 | 5 | 0.0427 | 1.681 | 0.6637 |
| m04q10 | churn | 510 | 10 | 0.0636 | 3.362 | 0.8712 |
| m04q20 | churn | 510 | 20 | 0.1272 | 6.725 | 0.8121 |
| m05q10 | churn | 640 | 10 | 0.0507 | 3.362 | 1.1084 |
| m05q20 | churn | 640 | 20 | 0.1014 | 6.725 | 1.0492 |
| m05q40 | churn | 640 | 40 | 0.2028 | 13.449 | 0.9308 |
| m06q10 | churn | 760 | 10 | 0.0427 | 3.362 | 1.3273 |
| m06q20 | churn | 760 | 20 | 0.0854 | 6.725 | 1.2681 |
| m07q20 | churn | 890 | 20 | 0.0729 | 6.725 | 1.5053 |
| m07q40 | churn | 890 | 40 | 0.1458 | 13.449 | 1.3869 |
| m03q10 | anular ondulado | 380 | 10 | 0.0854 | 3.362 | 0.6341 |
| m03q20 | anular ondulado | 380 | 20 | 0.1708 | 6.725 | 0.5749 |
| m03q40 | anular ondulado | 380 | 40 | 0.3415 | 13.449 | 0.4565 |
| m03q80 | anular ondulado | 380 | 80 | 0.6831 | 26.898 | 0.2197 |
| m04q40 | anular ondulado | 510 | 40 | 0.2545 | 13.449 | 0.6937 |
| m04q80 | anular ondulado | 510 | 80 | 0.5089 | 26.898 | 0.4569 |
| m05q80 | anular ondulado | 640 | 80 | 0.4056 | 26.898 | 0.6941 |
| m06q40 | anular ondulado | 760 | 40 | 0.1708 | 13.449 | 1.1498 |
| m06q80 | anular ondulado | 760 | 80 | 0.3415 | 26.898 | 0.9130 |
| m07q80 | anular ondulado | 890 | 80 | 0.2916 | 26.898 | 1.1502 |

## Resumo por padrão (min / máx / média / mediana)
| Padrão | n | Jg min | Jg máx | Jg média | Jg mediana | Jl min | Jl máx | Jl média | Jl mediana |
|--------|--:|-------:|-------:|---------:|-----------:|-------:|-------:|---------:|-----------:|
| bolha dispersa | 1 | 0.336 | 0.336 | 0.336 | 0.336 | 1.6178 | 1.6178 | 1.6178 | 1.6178 |
| plug | 3 | 0.336 | 1.681 | 0.785 | 0.336 | 1.1617 | 1.5941 | 1.3788 | 1.3806 |
| slug | 4 | 0.336 | 1.681 | 1.345 | 1.681 | 0.9245 | 1.3569 | 1.1394 | 1.1380 |
| churn | 10 | 1.681 | 13.449 | 6.556 | 6.725 | 0.6637 | 1.5053 | 1.0923 | 1.0788 |
| anular ondulado | 10 | 3.362 | 26.898 | 18.492 | 20.174 | 0.2197 | 1.1502 | 0.6943 | 0.6639 |

## Ensaios descartados (transições com dois padrões simultâneos)
- `m03q01` — plug-slug
- `m04q05` — slug-churn
- `m07q10` — slug-churn

> `m05q05_002`: o arquivo-marcador diz `__plug.txt`, mas a classificação manual (fonte de verdade) rotula m05q05 como **slug**; adotou-se **slug** (G e q″ idênticos a m05q05_001, logo mesma velocidade).

## Limitações (honestas)
1. **x_in = 0 (entrada líquido saturado):** ignora o sub-resfriamento do pré-aquecedor. x é estimativa de **equilíbrio (limite superior)**, não medida; não há temperatura de entrada nos dados brutos para corrigir. Consequência: ensaios anulares saem com x≈0,3–0,68 — abaixo do esperado (alto vazio) para o rótulo.
2. **G se cancela em Jg:** como x ∝ q″/G, Jg = (4L/(D·i_lv·ρ_v))·q″ depende **só de q″**. No plano Jg×Jl os ensaios de mesmo q″ alinham-se numa vertical; a informação de G está em Jl. Ensaio m07q01 (bolha) tem q″=1 kW/m², logo Jg=0,336 m/s — mesma posição horizontal que m05q01 e m06q01 (plug), separados somente por Jl (G diferente).
3. **x = Bo·4L/D:** o título é um reescalonamento do número de ebulição do trabalho-base (de Oliveira et al. 2021); o **escorregamento é ignorado** (modelo homogêneo). Logo Jg/Jl são ilustrativos do modelo, não velocidades medidas.
4. Ensaios q01 (q″=1 kW/m²) dão x≈0,004–0,006 (≈ monofásico); mantidos com piso em Jg para o eixo log.

## Verificações internas executadas
- A=7.854e-07 m²; mdot(m05)=0.503 g/s  [ok]
- x_out: m03q80=0.683 (máx), m07q01=0.0036 (mín)  [ok]
- x monotônico ↑ com q″ (G fixo)  [ok]
- Jg idêntico entre runs de mesmo q″ (cancelamento de G)  [ok]
- ensaios: {'bubble': 1, 'plug': 3, 'slug': 4, 'churn': 10, 'wavy annular': 10} (total 28); descartados 3  [ok]
