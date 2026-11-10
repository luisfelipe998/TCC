# v13 — Fração de vazio × velocidade superficial de líquido

Combina α cordal (v11, Otsu sobre panorama) com Jl modelo de equilíbrio (v12).

## Por ensaio
| Ensaio | Padrão | G [kg/m²s] | q″ [kW/m²] | α [-] | Jl [m/s] | Jg [m/s] | sat |
|--------|--------|----------:|-----------:|------:|---------:|---------:|-----|
| m05q01_001 | plug | 640 | 1 | 0.384 | 1.1617 | 0.336 |  |
| m06q01_001 | plug | 760 | 1 | 0.370 | 1.3806 | 0.336 |  |
| m07q05_001 | plug | 890 | 5 | 0.344 | 1.5941 | 1.681 |  |
| m04q01_001 | slug | 510 | 1 | 0.448 | 0.9245 | 0.336 |  |
| m05q05_001 | slug | 640 | 5 | 0.386 | 1.1380 | 1.681 |  |
| m06q05_001 | slug | 760 | 5 | 0.431 | 1.3569 | 1.681 |  |
| m03q05_001 | churn | 380 | 5 | 0.483 | 0.6637 | 1.681 |  |
| m04q10_001 | churn | 510 | 10 | 0.551 | 0.8712 | 3.362 |  |
| m04q20_001 | churn | 510 | 20 | 0.519 | 0.8121 | 6.725 |  |
| m05q10_001 | churn | 640 | 10 | 0.582 | 1.1084 | 3.362 |  |
| m05q20_001 | churn | 640 | 20 | 0.569 | 1.0492 | 6.725 |  |
| m05q40_001 | churn | 640 | 40 | 0.564 | 0.9308 | 13.449 |  |
| m06q10_001 | churn | 760 | 10 | 0.532 | 1.3273 | 3.362 |  |
| m06q20_001 | churn | 760 | 20 | 0.559 | 1.2681 | 6.725 |  |
| m07q20_001 | churn | 890 | 20 | 0.528 | 1.5053 | 6.725 |  |
| m07q40_001 | churn | 890 | 40 | 0.573 | 1.3869 | 13.449 |  |
| m03q10_001 | anular ondulado | 380 | 10 | 0.462 | 0.6341 | 3.362 |  |
| m03q20_001 | anular ondulado | 380 | 20 | 0.426 | 0.5749 | 6.725 |  |
| m04q40_001 | anular ondulado | 510 | 40 | 0.515 | 0.6937 | 13.449 |  |
| m04q80_001 | anular ondulado | 510 | 80 | 0.473 | 0.4569 | 26.898 |  |
| m05q80_001 | anular ondulado | 640 | 80 | 0.533 | 0.6941 | 26.898 |  |
| m06q40_001 | anular ondulado | 760 | 40 | 0.566 | 1.1498 | 13.449 |  |
| m06q80_001 | anular ondulado | 760 | 80 | 0.537 | 0.9130 | 26.898 |  |
| m07q80_001 | anular ondulado | 890 | 80 | 0.564 | 1.1502 | 26.898 |  |

## Resumo por padrão (α e Jl)
| Padrão | n | α mín | α máx | α média | Jl mín | Jl máx | Jl média |
|--------|--:|------:|------:|--------:|-------:|-------:|---------:|
| plug | 3 | 0.344 | 0.384 | 0.366 | 1.1617 | 1.5941 | 1.3788 |
| slug | 3 | 0.386 | 0.448 | 0.422 | 0.9245 | 1.3569 | 1.1398 |
| churn | 10 | 0.483 | 0.582 | 0.546 | 0.6637 | 1.5053 | 1.0923 |
| anular ondulado | 8 | 0.426 | 0.566 | 0.509 | 0.4569 | 1.1502 | 0.7833 |

## Notas
- **α cordal**: fração de vazio calculada por Otsu no panorama espaço-tempo (v11); reflete uma corda vertical central, não a seção transversal. Viés cordal > areal para padrões com núcleo gasoso central.
- **Jl modelo**: calculado via balanço de energia (modelo homogêneo, x_in=0, Tsat=27 °C, R600a) — estimativa de equilíbrio, não medido.
- m03q40_001 excluído da estatística (saturação 34% do núcleo → α artificialmente baixo).