# Características dos Padrões de Escoamento

Fontes: de Oliveira et al. (2019, 2021) — R-600a/R-290/R-1270, tubo 1 mm, CWT Morlet.
O paper de 2021 analisa **SO (sensor óptico)** e **PD (queda de pressão)** separadamente;
frequências são **bandas** (não picos únicos, salvo indicado). Este trabalho usa **SO = Voltage_2**.

## Sensor óptico (SO)

| Padrão | Característica SO | Freq. dominante SO (Hz) | Característica WPS |
|---|---|---|---|
| **Bolhas dispersas** | Rajadas intermitentes de energia ampla; sinal irregular | 4–35 (picos FFT); CWT em baixa freq. (p.14) | Ilhas de energia espalhadas em 1–128 Hz; sem banda persistente |
| **Plug** | Intermitente; dois centros de energia distintos | 4–16 (G=510); 2–16 (G=890) (p.15) | Duas regiões de energia em 2–16 Hz; picos isolados de alta freq. |
| **Slug** | Maior energia espectral que plug; bolhas desprendidas da cauda | ≤0,7 + **pico 2,8**; + 8/12/23/32 (G=510) (p.16) | Banda horizontal proeminente ~2,8 Hz; distribuição mais ampla que plug |
| **Churn** | Espectro amplo; instabilidades Kelvin-Helmholtz dominantes | sobe a **15–28** com fluxo mássico (p.16) | Banda de alta energia em 15–28 Hz; centro do retrato de fase denso |
| **Anular ondulado** | Frentes de onda na interface; instabilidades rápidas sobrepostas | 2–16; frentes de onda ≈ **20** (p.17) | Energia em 2–16 Hz; retrato de fase com órbitas espalhadas |
| **Anular liso** | Similar ao anular ondulado, menor amplitude; anti-persistente | não reportada (só Hurst/LLE) (p.17) | Sem banda persistente na TVO; mesma FFT de baixa freq. do anular ondulado |

## Queda de pressão (PD)

| Padrão | Característica PD | Freq. dominante PD (Hz) |
|---|---|---|
| **Bolhas dispersas** | Vários picos de energia; CWT persistente em baixa freq. | 6–100 (picos FFT); CWT em baixa freq. (p.14) |
| **Plug** | Faixa curta de baixa freq. por longos períodos | 0,5–1 (p.16) |
| **Slug** | Flutuações regidas por freq. muito baixas | < 1 (p.16) |
| **Churn** | Freq. dominantes de baixa energia | ≤ 8 (p.16) |
| **Anular ondulado** | Sem freq. característica na CWT; FFT ≈ churn | não observada por CWT; FFT ≈ churn (≤8) (p.17) |
| **Anular liso** | Sem freq. característica na CWT (só Hurst/LLE) | não reportada (p.17) |

Ressalvas: únicos picos únicos claros = slug SO (2,8 Hz) e anular ondulado SO (~20 Hz).
"df=1,9531 Hz" na Fig. 4(b) é resolução da FFT, **não** freq. dominante. Legendas das
Figs. 13/14 parecem trocadas; canais atribuídos pelo texto do corpo.

## Estimador de correlação (Schuck & Bodmann, 2019, cap. 32) — bottom line

O artigo define os estimadores de auto/crosscorrelação na forma **não-enviesada** (divisor `1/(N−m)`, eqs. 32.21–32.23), mas **não argumenta** essa escolha: apenas a atribui à referência Bendat & Piersol (*Random Data*). A robustez a ruído do método vem de outro mecanismo — a **média sobre K realizações** (`1/K Σ`, eqs. 32.18–32.20 e 32.24–32.26), com K=3 nos experimentos do artigo. O estimador não-enviesado é correto *na média*, então combina bem com essa média sobre realizações; além disso, o método só usa lags pequenos (transiente em τ=500 num sinal N=2048) e normaliza a WCO no lag 0 (`R_xx(0) ≥ |R_xx(τ)|`), permanecendo na região de baixa variância do estimador.

Consequência para a análise deste trabalho (regime **K=1**, uma realização por sinal, sem média): perde-se o mecanismo que doma a cauda de lag alto. O estimador não-enviesado, sem média, expõe a variância explosiva da cauda (poucos termos divididos por `N−m → 1`), que na CWT contamina a banda de baixa frequência e satura a WCO — o ruído de baixa freq observado nos gráficos. Fiel ao espírito do artigo, o remédio canônico seria **K>1 (blocos/janelas + média)**. As alternativas exploradas aqui (v3 = truncar a cauda em `m ≤ N/2`; v4 = estimador enviesado `1/N`, cujo divisor aplica uma janela triangular/Bartlett que amortece a cauda) são consertos válidos para K=1, mas **extensões deste trabalho, não propostas do artigo** — o paper não trunca nem usa `1/N`, ele media sobre realizações.
