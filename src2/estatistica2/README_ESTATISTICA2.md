# Follow-up estatístico — `src2/estatistica2/`

**Base:** 19 ensaios (10 *churn* + 9 *wavy annular*), reconstrução de imagem sob o motor
Schuck & Bodmann (`../v8_spatial_cropped_and_center/wavelet_v8.py`, **intocado**).
**Scripts:** `extract_descriptors.py`, `run_stats.py`, `plot_out_of_band.py`.
**Saídas:** `descriptors_{cropped,center}_{piso50,piso10}.npz`, `stats_results2.txt`,
`a2_ecdf.png`, `a3_out_of_band.png`. Executado 2026-10-03.

## Por que este follow-up existe

A análise de `../estatistica/` respondeu à pergunta central com ANOSIM R=0,177 (p
q-estratificado=0,020). A reverificação contra código e dados confirmou o núcleo A1, mas achou
três problemas que este follow-up corrige — **sem reabrir** a pasta anterior (congelada como
registro) e **sem alterar o motor**:

1. **A3 original era um no-op.** A banda de costura [18,5;37,0] px fica abaixo do piso de 50 px
   da CWT → a máscara cobria 0 escalas. Não testava nada.
2. **Conclusão frágil em A2.** "inter-vizinho > inter-distante reforça o sinal" valia só para a
   mediana; **inverte na média**. Cherry-pick.
3. **Confound de vazão (m) subenfatizado.** A1 só estratificava por q; o efeito pairwise de m
   sobre a distância (+0,16) rivaliza com o de padrão (−0,178).

## O que mudou (evolução, não reescrita)

- `extract_descriptors.py`: **`FREQ_MAX` virou parâmetro**. Gera descritores em dois pisos —
  `piso50` (λ≥50 px, igual ao motor) e `piso10` (λ≥10 px, inclui a região de pequena escala que
  o piso50 exclui). Âncora de regressão confirma `piso50/cropped` ≡ baseline de `../estatistica`
  (max|Δ| = 0,0). Assert de identidade FFT≡laço passa nos dois pisos.
- `run_stats.py`:
  - **A1** agora reporta p por enumeração exata **q-estratificada (200 arranjos)** E
    **m-estratificada (1728 arranjos)**, mais o p livre. Imprime R máximo e rank de R_obs.
  - **A2** reporta **mediana E média** por grupo (expõe a inversão) e destaca as correlações
    pairwise `corr(D,|Δq|)`, `corr(D,|Δm|)`, `corr(D,same_pattern)` como a evidência robusta.
  - **A3′** substitui a máscara no-op por: (b) comparar R em piso50 vs piso10 e (a, figura)
    medir ONDE mora a energia e se há bump na banda de costura.

## Resultados

| Métrica | Valor |
|---|---|
| ANOSIM R (cropped, piso50, principal) | **0,177** |
| p q-estratificado (200 arranjos, p_min 0,005) | 0,020 (rank 3; R_max 0,242) |
| **p m-estratificado (1728 arranjos, p_min 0,00058)** | **0,0156** (rank 26; R_max 0,250) |
| p não-estratificado (5000 MC) | 0,007 |
| ANOSIM R (center, piso50) | 0,179 |
| ECDF intra — mediana / média | 0,098 / 0,136 |
| ECDF inter-vizinho — mediana / média | 0,159 / 0,175 |
| ECDF inter-distante — mediana / média | 0,135 / 0,184 |
| corr(D, \|Δq\|) | **−0,004** |
| corr(D, \|Δm\|) | +0,160 |
| corr(D, mesmo padrão) | −0,178 |
| R piso10 (cropped) vs piso50 | 0,354 vs 0,177 (Δ=+0,177) |
| Energia na banda de costura (cropped, piso10) | **1,1%** do \|WPS b=0\| |

## Vereditos

1. **A assinatura fraca sobrevive ao controle do confound mais bem condicionado.** A
   m-estratificação admite 1728 arranjos (vs 200 da de q, que tem estratos puros q=5/q=80) e
   nenhum estrato puro; R=0,177 dá p=0,0156. Logo a separação não é só efeito de q nem só de m.
   Continua **efeito-tamanho descritivo**, nunca "significância forte".

2. **A separação NÃO é proximidade de ponto operacional em q.** `corr(D,|Δq|)≈0` (robusto). A
   antiga leitura via ordem das medianas de A2 era frágil (inverte na média; as ECDFs inter quase
   se sobrepõem) — abandonada. Em contrapartida, **m é um confound real** (`corr(D,|Δm|)=+0,16`,
   mesma ordem do efeito de padrão −0,178), não controlado — declarado, não consertado.

3. **A costura NÃO contamina o descritor com um bump.** Rebaixar o piso para 10 px dobra R
   (0,177→0,354), MAS a banda de costura carrega só **1,1%** da energia e **não é o pico** (o
   pico fica em λ≈360–950 px; ver `a3_out_of_band.png`). A subida de R vem do **conteúdo difuso
   de pequena escala (10–50 px)** — onde annular (sinal mais longo/fino) e churn diferem, região
   colinear com velocidade/comprimento e sujeita a borda/COI. Como essa região é ambígua, o
   **piso50 (A1 principal) é a leitura CONSERVADORA e defensável**; piso10 é só sonda de robustez,
   nunca número de manchete. O center quase não muda (Δ=+0,019), consistente.

## O que a análise NÃO responde (declarado)

"Qual é o padrão?" (sem classificação/generalização); escala espacial física (eixo x sintético);
atribuição causal ao padrão (m e q colineares com padrão; sem réplicas verdadeiras — 1 run por
célula (m,q), "mesmo padrão" = pontos operacionais distintos com mesmo rótulo de regime);
acoplamento entre padrões (WCO cruzado some sob shift, só figura). **Tudo é propriedade do
descritor reconstruído, não do escoamento físico.**

## Deliberadamente fora (não é rigor, é cerimônia com 19 pontos / 2 grupos / sem réplicas)

PERMANOVA + PERMDISP (r≈0,98 com ANOSIM R); ML/classificador; IC bootstrap de R (o nulo
enumerado JÁ é a distribuição exata); correção de múltiplas comparações entre A1/A2/A3′
(perguntas distintas, não uma família).

## Reproduzir

```
python extract_descriptors.py   # ~8–12 min: 76 leituras xlsx + 2 asserts de identidade; gera 4 .npz
python run_stats.py             # A1/A2/A3′ → stats_results2.txt + a2_ecdf.png
python plot_out_of_band.py      # A3′-a → a3_out_of_band.png
```
