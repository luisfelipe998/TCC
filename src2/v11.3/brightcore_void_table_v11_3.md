# Vazio bright-core × taxa de troca temporal (v11.3 A — PRINCIPAL)

Análise **principal** do par v11.3. Eixo X = fração de vazio **bright-core** (maior núcleo CLARO contíguo pelo centro, quadro bruto) — método do v11.2 verificado adversarialmente que RECUPERA anular > churn (o α cordal 'gás = escuro' do v11/v11.1 invertia por causa da óptica W/M: núcleo de gás CLARO, interfaces ESCURAS). Eixo Y = **taxa de troca temporal** gás↔líquido (flip-rate, v11.2 M4), eixo ORTOGONAL de estabilidade (churn caótico > anular estável).

Padrões: dispersed bubbles, plug, slug, churn, wavy annular — 27 ensaios (26 na estatística).

## Mediana por padrão (saturados fora)

| padrão | n | bright-core X (M1) | flip-rate Y (M4) | continuidade axial | α cordal (v11.1, invertido) |
|---|---|---|---|---|---|
| dispersed bubbles | 2 | 0.500 | 0.0233 | 0.972 | 0.210 |
| plug | 3 | 0.333 | 0.0177 | 0.978 | 0.370 |
| slug | 3 | 0.317 | 0.0174 | 0.978 | 0.431 |
| churn | 10 | 0.283 | 0.0300 | 0.974 | 0.555 |
| annular | 8 | 0.342 | 0.0177 | 0.984 | 0.524 |

**Leitura:** no eixo X (bright-core) **anular > churn** (ordenação física recuperada; o α cordal invertido fica na última coluna só para contraste). No eixo Y (flip-rate) **churn > anular** — churn é o mais instável. Anular = vazio ALTO + estável (canto inferior-direito); churn = vazio alto + caótico (topo).

## Controle de brilho — bright-core sem nível absoluto

`bright-core-mid` usa limiar = (min+max)/2 por coluna, removendo o NÍVEL de cinza. Se anular>churn sobrevive, a separação é de FORMA do perfil, não de brilho.

| padrão | bright-core (Otsu) | bright-core-mid (sem nível) |
|---|---|---|
| dispersed bubbles | 0.500 | 0.558 |
| plug | 0.333 | 0.300 |
| slug | 0.317 | 0.317 |
| churn | 0.283 | 0.283 |
| annular | 0.342 | 0.342 |

## Por ensaio (transparência)

| scenario | pattern | G | q | bright-core | flip-rate | α cordal | sat % | flag |
|---|---|---|---|---|---|---|---|---|
| m07q01_001 | dispersed bubbles | 890 | 1 | 0.717 | 0.0067 | 0.017 | 0 |  |
| m03q40_001 | wavy annular | 380 | 40 | 0.483 | 0.0107 | 0.275 | 34 | saturado (sat=34%) |
| m03q20_001 | wavy annular | 380 | 20 | 0.383 | 0.0121 | 0.426 | 14 |  |
| m03q10_001 | wavy annular | 380 | 10 | 0.367 | 0.0154 | 0.462 | 12 |  |
| m04q80_001 | wavy annular | 510 | 80 | 0.367 | 0.0110 | 0.473 | 7 |  |
| m03q05_001 | churn | 380 | 5 | 0.350 | 0.0192 | 0.483 | 12 |  |
| m05q01_001 | plug | 640 | 1 | 0.350 | 0.0124 | 0.384 | 0 |  |
| m05q80_001 | wavy annular | 640 | 80 | 0.350 | 0.0194 | 0.533 | 0 |  |
| m04q01_001 | slug | 510 | 1 | 0.333 | 0.0156 | 0.448 | 0 |  |
| m04q40_001 | wavy annular | 510 | 40 | 0.333 | 0.0159 | 0.515 | 8 |  |
| m06q01_001 | plug | 760 | 1 | 0.333 | 0.0177 | 0.370 | 0 |  |
| m04q20_001 | churn | 510 | 20 | 0.317 | 0.0196 | 0.519 | 10 |  |
| m05q05_001 | slug | 640 | 5 | 0.317 | 0.0174 | 0.386 | 0 |  |
| m07q05_001 | plug | 890 | 5 | 0.317 | 0.0192 | 0.344 | 0 |  |
| m04q10_001 | churn | 510 | 10 | 0.300 | 0.0176 | 0.551 | 4 |  |
| m05q40_001 | churn | 640 | 40 | 0.300 | 0.0286 | 0.564 | 0 |  |
| m06q05_001 | slug | 760 | 5 | 0.300 | 0.0254 | 0.431 | 0 |  |
| m06q80_001 | wavy annular | 760 | 80 | 0.300 | 0.0203 | 0.537 | 0 |  |
| m05q05_002 | dispersed_bubbles | 640 | 5 | 0.283 | 0.0398 | 0.402 | 0 |  |
| m05q10_001 | churn | 640 | 10 | 0.283 | 0.0308 | 0.582 | 0 |  |
| m05q20_001 | churn | 640 | 20 | 0.283 | 0.0304 | 0.569 | 0 |  |
| m06q10_001 | churn | 760 | 10 | 0.283 | 0.0296 | 0.532 | 0 |  |
| m06q40_001 | wavy annular | 760 | 40 | 0.283 | 0.0226 | 0.566 | 0 |  |
| m07q80_001 | wavy annular | 890 | 80 | 0.267 | 0.0221 | 0.564 | 0 |  |
| m06q20_001 | churn | 760 | 20 | 0.250 | 0.0323 | 0.559 | 0 |  |
| m07q20_001 | churn | 890 | 20 | 0.250 | 0.0333 | 0.528 | 0 |  |
| m07q40_001 | churn | 890 | 40 | 0.233 | 0.0332 | 0.573 | 0 |  |

## Ressalvas honestas

- **Ordenação, não magnitude.** bright-core ~0,35 é cordal em pixels (proxy óptico RELATIVO não-calibrado); longe do areal ~0,9 da literatura. Sem px→mm nem segunda vista, o vazio areal verdadeiro NÃO é recuperável destas imagens.
- **dispersed bubbles (n=2):** sai alto ESPÚRIO (centro coincidentemente claro — a regra 'centro claro' precisa de porta de contraste); anotado no mapa, fora da estatística.
- **Sobreposição e n pequeno** (anular 8 × churn 10); plug empata com anular. A separação é de tendência central.
- **Eixo Y (flip-rate)** é ESTABILIDADE, não magnitude: põe churn no topo (mais caótico), consistente com o anular embaixo (estável).
- **Significância é de DIREÇÃO, não de p forte** (re-verificação adversarial independente, `holds=true`): a ordem anular>churn NUNCA inverte sob perturbação (9 variantes de janela/limiar/polaridade) e sobrevive ao q-match (o q-confound joga CONTRA — anular é amostrado em q maior mas bright-core correla NEGATIVAMENTE com q em ambos os padrões, r≈−0,4..−0,5). MAS o p<0,05 é FRÁGIL nas bordas: p=0,022 na janela 10:70, p≈0,002 em 12:68–15:65, mas sobe a p=0,11 na janela larga 8:72 (dilui nas paredes saturadas) e p=0,24 com Otsu+10. Reportar como efeito borderline p≈0,02 COM a varredura de perturbação, não como p forte. Sem correção de múltiplas comparações.
