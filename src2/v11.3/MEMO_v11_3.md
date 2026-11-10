# MEMO v11.3 — fração de vazio × intermitência, com o reflexo W/M normalizado

**Duas análises irmãs**, cada uma separada (pedido do dono), ambas com o mesmo
eixo Y = **taxa de troca temporal gás↔líquido** (flip-rate, v11.2 método 4).
`m05q05_002` tratado como **dispersed bubbles** (o marcador `__*.txt` manda;
n=2 com `m07q01_001`).

| | script | eixo X | papel |
|---|---|---|---|
| **A** | `brightcore_void_v11_3.py` | vazio **bright-core** (núcleo CLARO contíguo, quadro bruto) | **PRINCIPAL** — defensável, recupera a física |
| **B** | `areal_flatfield_void_v11_3.py` | vazio **AREAL** 2-D sobre quadro flat-field (área de gás/área seção) | **SUPLEMENTAR** — prova a remoção do reflexo; dá a métrica de ÁREA |

Reusam `void_common_v11_3.py` (portado verbatim do v11.2). Não alteram v11.1/v11.2.

## Por que normalizar: a óptica W/M

O microtubo de 1 mm retroiluminado NÃO é "líquido-claro/gás-escuro". O perfil
radial é **W/M**: núcleo de gás CLARO no centro (linha ~40), interfaces/meniscos
ESCUROS fora do centro (linhas ~20/60), paredes saturadas (0–7). Logo
"gás = pixel escuro" (α cordal do v11/v11.1) **inverte o sinal** — conta a sombra
de interface como gás e o núcleo de gás como líquido. Ver
`flatfield_radial_proof_v11_3.png` (painel esquerdo = o W/M bruto).

## A — bright-core (PRINCIPAL): ordenação física recuperada

| padrão | n | bright-core X | flip-rate Y |
|---|---|---|---|
| dispersed bubbles | 2 | 0.283 / 0.717 (instável) | — |
| plug | 3 | 0.333 | 0.0177 |
| slug | 3 | 0.317 | 0.0174 |
| churn | 10 | **0.283** | **0.0300** |
| annular | 8 | **0.342** | 0.0177 |

- Eixo X: **anular 0,342 > churn 0,283** — ordenação física RECUPERADA
  (MWU p=0,022 unilateral, AUC 0,79). Sobrevive ao controle de brilho sem nível
  (`bright-core-mid`) e ao q-match [10,40] (verificado no v11.2).
- Eixo Y: **churn 0,030 > anular 0,018** (p=0,006, AUC 0,86) — churn é o mais
  INSTÁVEL. Eixo ortogonal de estabilidade, imune a brilho/q.
- Leitura do mapa: **anular = vazio ALTO + estável** (inferior-direito);
  **churn = vazio alto + caótico** (topo); plug/slug no meio; dispersed bubbles
  espalhado (n=2).
- **Re-verificação adversarial independente (`holds=true`):** re-implementação
  do zero reproduz o baseline a 3 decimais. A ordem anular>churn NUNCA inverte
  sob perturbação (9 variantes) e sobrevive ao q-match — o q-confound joga
  CONTRA (anular amostrado em q maior, mas bright-core correla NEGATIVAMENTE com
  q em ambos, r≈−0,4..−0,5, logo o anular vence apesar da penalização). Ressalva:
  o p<0,05 é de DIREÇÃO robusta mas FRÁGIL nas bordas (p=0,022 em 10:70, ≈0,002
  em 12:68–15:65, sobe a 0,11 em 8:72 e 0,24 em Otsu+10) — reportar como efeito
  borderline p≈0,02 com a varredura, sem correção de múltiplas comparações.

## B — areal flat-field (SUPLEMENTAR): reflexo removido, mas annular não no topo

A normalização usa **flat-field por baseline de líquido** (percentil 20 por pixel
= fundo que NÃO acompanha o núcleo de gás móvel), gás = excesso claro Fn>1,3.
(O flat-field por mediana temporal foi testado e DESCARTADO: autonormaliza tudo
para mediana ≈1 e destrói o sinal.)

- **Remoção do reflexo PROVADA:** std radial interior 33 (bruto) → 0,38
  (flat-field), ~88× em média (até ~300× nos anulares de baixo G). Ver a figura
  de prova e a coluna std radial da tabela.
- **Ordem recuperada:** dispersed 0,245 < plug 0,323 < slug 0,367 <
  **annular 0,385 < churn 0,453** — correta até slug, mas o **annular CAI abaixo
  do churn** (não fica no topo).
- **Por quê (limite óptico, não bug):** em wavy annular as bandas de interface
  ESCURAS (linhas ~20/60) persistem mesmo quando o núcleo claro enche o centro;
  o limiar de área clara SUBCONTA o annular. É limitação da vista única + W/M.
  Só o bright-core (A) coloca o annular acima do churn.

## Expectativa da literatura (para comparação honesta)

Monotônico esperado **dispersed bubbles < plug < slug < churn < annular**, annular
o MAIOR (`deoliveira2021` mesmo laboratório; `dias2025knn` "void increases from
slug to churn", sobreposição churn≈annular; `schepperle2024` define
α = área de vapor/área da seção por visão computacional). Valores absolutos
slip-corrigidos: plug 0,04–0,29, churn 0,61–0,85, annular 0,86–0,99.

## Ressalvas obrigatórias (nas duas saídas)

1. **Proxy óptico RELATIVO não-calibrado** (sem px→mm, vista única) — a ordenação
   ENTRE padrões é legítima; o valor absoluto NÃO é a fração de área da seção
   (bright-core ~0,35 e areal ~0,4 ≠ ~0,9 areal real). Coerente com
   `assets/paper-log.md:11`.
2. **Bright-core (A)** recupera a ORDENAÇÃO annular>churn, não a magnitude.
3. **Areal (B)** inverte churn/annular (annular não no topo) — limite óptico.
4. **dispersed bubbles** instável/espúrio: em A os dois reps são 0,283 e 0,717;
   em B são 0,000 e 0,490. n=2, fora da estatística, anotado no mapa.
5. **Dependência de iluminação residual** em B (corr vazio×saturação +0,31);
   q-confound não testado por subconjunto q-match em B.
6. `m03q80` (smooth annular) excluído; `m03q40` fora da estatística (sat 34%).

## Saídas

`brightcore_void_table_v11_3.md`, `brightcore_void_map_v11_3.png`,
`areal_flatfield_void_table_v11_3.md`, `areal_flatfield_void_map_v11_3.png`,
`flatfield_radial_proof_v11_3.png`, `void_common_v11_3.py`.

**Recomendação ao TCC:** ancorar a conclusão na figura A (bright-core,
annular>churn verificado); usar B como figura de controle que demonstra a remoção
do reflexo e expõe honestamente o teto óptico (annular não recuperável no topo por
limiar de área). Ambas rotuladas como proxy óptico relativo. Verificação
adversarial final (q-match + perturbação de janela/limiar/polaridade no bright-core
portado) recomendada antes de fixar no texto.
