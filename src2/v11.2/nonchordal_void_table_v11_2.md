# Fração de vazio NÃO-CORDAL por método × padrão (v11.2)

Estudo de métodos NÃO-cordais para recuperar **anular > churn** no vazio — que a métrica cordal α do v11/v11.1 não consegue (churn ≳ anular, invertido). Causa-raiz óptica: o núcleo de gás aparece **CLARO** e as interfaces junto às paredes aparecem **ESCURAS** (perfil radial em W/M), então contar 'gás = pixel escuro' inverte o sinal. Ver `MEMO_v11_2.md` para a leitura completa.

Quatro métodos, cada um verificado adversarialmente (workflow `v11-nonchordal-void`):

| método | ideia | recupera anular>churn? | veredito | principiado? |
|---|---|---|---|---|
| 1 bright-core | maior núcleo CLARO contíguo pelo centro (quadro bruto) | **sim** (0,35 vs 0,28, p≈0,018) | **CONFIRMADO** | **sim** |
| 2 wall-film | largura de núcleo, limiar por-coluna | sim (sem sobreposição, p=0,0004) | PARCIAL | escopado (ordem global dos 5 erra) |
| 3 radial-shape (D1) | centro−parede normalizado | parcial | **REFUTADO** | não (é brilho re-embalado, corr 0,95) |
| 4 flip-rate | troca gás↔líquido por pixel (tempo) | sim, como *estabilidade* | **CONFIRMADO** | **sim** |

Padrões: dispersed bubbles, plug, slug, churn, wavy annular — 26 ensaios (25 na estatística).

## Mediana por padrão (ensaios saturados fora)

| padrão | n | α cordal (v11.1) | bright-core (M1) | wall-film (M2) | flip-rate (M4) |
|---|---|---|---|---|---|
| dispersed bubbles | 1 | 0.017 | 0.717 | 0.431 | 0.0067 |
| plug | 3 | 0.370 | 0.333 | 0.369 | 0.0177 |
| slug | 3 | 0.431 | 0.317 | 0.373 | 0.0174 |
| churn | 10 | 0.555 | 0.283 | 0.361 | 0.0300 |
| annular | 8 | 0.524 | 0.342 | 0.456 | 0.0177 |

**Leitura:** no α cordal (v11.1) churn ≳ anular (invertido). Os métodos 1 e 2 (núcleo CLARO) recuperam **anular > churn**. No método 4 (flip-rate) **churn > anular** (churn é o mais instável) — eixo ORTOGONAL de estabilidade, não de magnitude.

## Controle de brilho — método 1 sem nível absoluto

`bright-core-mid` usa limiar = (min+max)/2 por coluna, removendo o NÍVEL de cinza. Se anular>churn sobrevive, a separação é de FORMA do perfil, não de brilho.

| padrão | bright-core (Otsu) | bright-core-mid (sem nível) |
|---|---|---|
| dispersed bubbles | 0.717 | 0.833 |
| plug | 0.333 | 0.300 |
| slug | 0.317 | 0.317 |
| churn | 0.283 | 0.283 |
| annular | 0.342 | 0.342 |

## Estrutura radial — continuidade axial do núcleo (M1)

Correlação de continuidade axial do núcleo claro. Anular (filme liso, núcleo contínuo) > churn (núcleo esfarrapado).

| padrão | continuidade axial |
|---|---|
| dispersed bubbles | 0.988 |
| plug | 0.978 |
| slug | 0.978 |
| churn | 0.974 |
| annular | 0.984 |

## Por ensaio (transparência)

| scenario | pattern | G | q | α cordal | bright-core | wall-film | flip-rate | sat % | flag |
|---|---|---|---|---|---|---|---|---|---|
| m07q01_001 | dispersed bubbles | 890 | 1 | 0.017 | 0.717 | 0.431 | 0.0067 | 0 |  |
| m03q40_001 | wavy annular | 380 | 40 | 0.275 | 0.483 | 0.545 | 0.0107 | 34 | saturado (sat=34%) |
| m03q20_001 | wavy annular | 380 | 20 | 0.426 | 0.383 | 0.467 | 0.0121 | 14 |  |
| m03q10_001 | wavy annular | 380 | 10 | 0.462 | 0.367 | 0.449 | 0.0154 | 12 |  |
| m04q80_001 | wavy annular | 510 | 80 | 0.473 | 0.367 | 0.487 | 0.0110 | 7 |  |
| m03q05_001 | churn | 380 | 5 | 0.483 | 0.350 | 0.412 | 0.0192 | 12 |  |
| m05q01_001 | plug | 640 | 1 | 0.384 | 0.350 | 0.412 | 0.0124 | 0 |  |
| m05q80_001 | wavy annular | 640 | 80 | 0.533 | 0.350 | 0.475 | 0.0194 | 0 |  |
| m04q01_001 | slug | 510 | 1 | 0.448 | 0.333 | 0.387 | 0.0156 | 0 |  |
| m04q40_001 | wavy annular | 510 | 40 | 0.515 | 0.333 | 0.443 | 0.0159 | 8 |  |
| m06q01_001 | plug | 760 | 1 | 0.370 | 0.333 | 0.369 | 0.0177 | 0 |  |
| m04q20_001 | churn | 510 | 20 | 0.519 | 0.317 | 0.393 | 0.0196 | 10 |  |
| m05q05_001 | slug | 640 | 5 | 0.386 | 0.317 | 0.373 | 0.0174 | 0 |  |
| m07q05_001 | plug | 890 | 5 | 0.344 | 0.317 | 0.361 | 0.0192 | 0 |  |
| m04q10_001 | churn | 510 | 10 | 0.551 | 0.300 | 0.370 | 0.0176 | 4 |  |
| m05q40_001 | churn | 640 | 40 | 0.564 | 0.300 | 0.398 | 0.0286 | 0 |  |
| m06q05_001 | slug | 760 | 5 | 0.431 | 0.300 | 0.348 | 0.0254 | 0 |  |
| m06q80_001 | wavy annular | 760 | 80 | 0.537 | 0.300 | 0.462 | 0.0203 | 0 |  |
| m05q10_001 | churn | 640 | 10 | 0.582 | 0.283 | 0.339 | 0.0308 | 0 |  |
| m05q20_001 | churn | 640 | 20 | 0.569 | 0.283 | 0.354 | 0.0304 | 0 |  |
| m06q10_001 | churn | 760 | 10 | 0.532 | 0.283 | 0.343 | 0.0296 | 0 |  |
| m06q40_001 | wavy annular | 760 | 40 | 0.566 | 0.283 | 0.429 | 0.0226 | 0 |  |
| m07q80_001 | wavy annular | 890 | 80 | 0.564 | 0.267 | 0.433 | 0.0221 | 0 |  |
| m06q20_001 | churn | 760 | 20 | 0.559 | 0.250 | 0.369 | 0.0323 | 0 |  |
| m07q20_001 | churn | 890 | 20 | 0.528 | 0.250 | 0.321 | 0.0333 | 0 |  |
| m07q40_001 | churn | 890 | 40 | 0.573 | 0.233 | 0.350 | 0.0332 | 0 |  |

## Ressalvas honestas

- **Ordenação, não magnitude.** bright-core ~0,35 é cordal em pixels, longe do areal ~0,9 da literatura. Sem calibração px→mm nem segunda vista, o vazio areal verdadeiro NÃO é recuperável.
- **Sobreposição e n pequeno** (anular 8 × churn 10): a separação é de tendência central.
- **Não é vazio universal dos 5 padrões:** dispersed bubbles sai alto espúrio (centro coincidentemente claro — a regra 'centro claro' precisa de porta de contraste); plug empata com anular. Validado para a pergunta anular×churn, não como vazio universal.
- **Método 3 (D1) REFUTADO:** é brilho re-embalado (corr(D1, brilho)=0,95); empata ao casar brilho.
