# Conclusões — Estatística mínima sobre a assinatura WPS (v8/v8_1)

**Pasta:** `src2/estatistica/` · **Scripts:** `extract_descriptors.py`, `run_stats.py`
**Saídas:** `descriptors_{cropped,center}.npz`, `stats_results.txt`, `a2_ecdf.png`
**Base:** 19 ensaios — 10 *churn* + 9 *wavy annular* — reconstrução de imagem sob motor
Schuck & Bodmann (`wavelet_v8.py`, intocado). Data: 2026-09-30.

## Pergunta

A **forma do WPS auto de Schuck (b=0, whitened, resumido por RUN)** é mais parecida entre
réplicas do **mesmo** padrão do que entre padrões **diferentes**, depois de casar o fluxo de
calor q — reconhecida como **propriedade do descritor reconstruído, não do escoamento físico**?

## Método (núcleo mínimo, sem redundância)

- **Descritor:** perfil de escala do WPS auto na coluna de lag-zero (b=0, índice N−1),
  truncado a 7001 px, whitened por escala (resíduo de `log|WPS|` contra `log(escala)`, remove o
  piso 1/f). Um vetor de 200 escalas por run.
- **A1 (resposta):** ANOSIM R de Clarke sobre `D = 1 − Pearson`, com p por **enumeração exata**
  dos 200 arranjos de rótulo q-estratificados (único teste; PERMANOVA-F e PERMDISP descartados
  por redundância — corr 0,98 com R).
- **A2 (baseline):** ECDF das distâncias intra vs inter (inter separado em q-vizinho/distante).
- **A3 (robustez):** o mesmo R em cropped/center × banda-de-costura completa/mascarada.

## Resultados

| Métrica | Valor |
|---|---|
| ANOSIM R (cropped, principal) | **0,177** |
| p estratificado por q (exato, 200 arranjos) | **0,020** (p mínimo possível = 0,005) |
| p não-estratificado (5000 MC, contexto) | 0,007 |
| ANOSIM R (center) | 0,179 |
| ECDF mediana — intra | **0,098** |
| ECDF mediana — inter q-vizinho | 0,159 |
| ECDF mediana — inter q-distante | 0,135 |

Confounds: churn m̄=5,20 q̄=19,5 n̄=8832; annular m̄=4,56 q̄=52,2 n̄=18140.
**corr(shift_mean, n_samples) = 1,0000** (mesma variável de velocidade).

## Conclusão

1. **Existe uma assinatura fraca, repetível e fracamente separável** no descritor. R=0,177
   (>0) com direção correta: distâncias intra (mediana 0,098) são menores que inter (0,135–0,159).
   O efeito sobrevive à estratificação por q (p exato = 0,020), i.e. **não é só efeito de q**.
2. **É fraco e descritivo, não inferencial forte.** A estratificação por q admite só 200
   arranjos → p mínimo 0,005; nunca reportar como significância forte. Reenquadrar sempre como
   *"o pipeline produz uma assinatura repetível e fracamente separável"*, jamais como
   propriedade física do escoamento.
3. **A separação não é efeito trivial de distância em q:** inter q-*vizinho* (0,159) é até
   MAIOR que inter q-*distante* (0,135). Se q dominasse, esperar-se-ia o contrário. Reforça (1).
4. **Estável entre cropped e center** (R = 0,177 vs 0,179) — o resumo 1D não muda a conclusão.

## O que a análise NÃO responde (declarado)

"Qual é o padrão?" (sem classificação); escala espacial física (eixo x sintético); atribuição
causal ao padrão (m não controlado; q colinear — q≤5 só churn, q≥80 só annular); acoplamento
entre padrões (WCO cruzado some sob shift, só figura).

## ⚠️ Inconsistência encontrada na execução — A3 é inaplicável como especificado

A banda de costura λ∈[0.7, 1.4]×shift_mean cai em **[18,5; 37,0] px** (média), mas o eixo de
análise da CWT tem **piso de 50 px** (`FREQ_MAX = 1/50` em `wavelet_v8.py`). Logo a máscara
cobre **0 escalas** e "banda-completa" ≡ "banda-mascarada" (ver `stats_results.txt`). Só 4 dos
19 runs têm a *borda superior* da banda (1,4×shift) passando de 50 px (máx. 64,5 px), e ainda
assim só o sliver 50–64 px.

**Implicação honesta:** o artefato de costura previsto (bump em λ≈shift_mean) está **fora da
banda analisada** — o descritor WPS b=0, por construção (piso 50 px), nunca o vê. Isto é uma
*boa* notícia para a validade de A1 (o R não pode ser inflado por um bump que não está no
descritor), mas significa que **A3 não testa o que pretendia**: não confirma ausência de
artefato, apenas mostra que o artefato mora abaixo do piso. Ver plano de follow-up.

> Correção aplicada durante a execução: a máscara comparava `scales` cru (46,9–2812) em vez de
> comprimento de onda `scales/center_freq` (50–3000 px). Corrigido; o resultado (0 escalas)
> mantém-se — é substantivo (banda abaixo do piso), não bug de unidade.
