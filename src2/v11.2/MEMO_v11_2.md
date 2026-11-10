# MEMO — v11.2: dá para recuperar **anular > churn** no vazio com um método NÃO-cordal?

> **Para revisão posterior.** Este memo documenta os testes do `src2/v11.2/`.
> Nada do v11.1 foi alterado. Resumo em uma frase: **sim, a ORDEM anular > churn
> é recuperável por um método não-cordal principiado (o núcleo de gás aparece
> CLARO, não escuro) — mas é recuperação de ordenação entre padrões, não da
> MAGNITUDE do vazio, que continua não-recuperável destas imagens.**

---

## 1. O problema que motivou o v11.2

No v11.1 (e no v11) a fração de vazio é **cordal**: numa corda vertical no centro
do tubo contamos "gás = pixel escuro" (limiar de Otsu) e tiramos a média. Nessa
métrica:

| padrão | α cordal (mediana) |
|---|---|
| churn | **0,55** |
| anular | **0,52** |

O churn sai **ligeiramente acima** do anular — fisicamente invertido. O anular
ondulado deveria ter o **maior** vazio (a literatura reporta α_anular ~ 0,9
*areal*). O orientador apontou o estranhamento e a sobreposição churn×anular.

A pergunta do v11.2: **uma análise NÃO-cordal — explorando a estrutura radial /
2-D / temporal — recupera anular > churn?**

---

## 2. A causa-raiz (achada no quadro bruto): a óptica é "W/M", não "gás-escuro"

Este foi o achado central, confirmado independentemente na verificação
adversarial e reproduzível em `radial_profile_v11_2.png`.

A óptica deste microtubo de 1 mm sob retroiluminação lateral **não** é
"líquido-claro / gás-escuro". O perfil radial de intensidade (média no tempo) tem
forma de **W/M**:

```
parede próxima (reflexo, claro)
   → BANDA ESCURA  (menisco/interface, linha ~20)
      → NÚCLEO CLARO no centro  (linha ~40, intensidade ~150 no anular / ~78 no churn)
   → BANDA ESCURA  (interface, linha ~60)
→ parede distante (escura)
```

Ou seja: **o núcleo de gás aparece CLARO** (a luz atravessa o núcleo) e as
**interfaces junto às paredes aparecem ESCURAS**. Logo a contagem ingênua
"gás = pixel escuro" **inverte o sinal exatamente onde importa** — conta os anéis
escuros de interface como gás e o núcleo claro de gás como líquido. É por isso
que o churn ≳ anular vazava na métrica cordal, e por que nenhum limiar estático
sobre "escuro" (cordal, areal, global, CLAHE) separa a *magnitude* dos dois.

> No `radial_profile_v11_2.png`: anular e churn têm o **mesmo formato W/M**, mas
> o núcleo central do anular é **nitidamente mais claro** (~150 vs ~78). É essa
> diferença de FORMA/brilho do núcleo que os métodos que funcionam exploram.

---

## 3. Os quatro métodos testados e o veredito de cada um

Cada método foi rodado como uma sonda independente e depois **verificado
adversarialmente** (um segundo agente tentou refutá-lo), no workflow
`v11-nonchordal-void`. Reproduzidos em `nonchordal_void_v11_2.py`.

| # | método | ideia | recupera anular>churn? | veredito | principiado? |
|---|---|---|---|---|---|
| 1 | **bright-core** | maior núcleo CLARO contíguo pelo centro (quadro bruto) | **sim** (0,34 vs 0,28; p≈0,02) | **CONFIRMADO** | **sim** |
| 2 | wall-film | largura de núcleo, limiar por-coluna | sim, sem sobreposição (p=0,0002) | **PARCIAL** | escopado (ordem global dos 5 erra) |
| 3 | radial-shape (D1) | centro−parede normalizado | parcial (9–12/13 pares) | **REFUTADO** | não — é brilho re-embalado (corr 0,95) |
| 4 | **flip-rate** | troca gás↔líquido por pixel (tempo) | sim, como **estabilidade** | **CONFIRMADO** | **sim** |

### Método 1 — central-bright-core (o que funciona para a magnitude)
Por coluna axial do quadro bruto, o núcleo de gás = **maior sequência contígua de
pixels CLAROS que contém o centro**; a largura é o vazio "bright-core". Mediana:

- **anular 0,342 > churn 0,283**, Mann-Whitney unilateral **p ≈ 0,02**, AUC 0,79.

Por que é **principiado, não ajustado para dar a resposta**:
- o sinal (anular−churn > 0) **não troca** em nenhuma janela interior
  (8:72 … 15:65) nem offset de limiar (−10 … +20) — 16/16 combinações;
- **sem o confundimento de brilho:** um limiar por-coluna (min+max)/2, que remove
  o nível absoluto de cinza, ainda dá anular > churn (p≈0,002) → é a FORMA do
  perfil (platô central largo), não o brilho;
- **sem o confundimento de q:** restringindo a q ∈ [10,40] (faixa comum), anular
  0,342 > churn 0,283 se mantém;
- leave-one-out: a direção sobrevive à remoção de qualquer ensaio.

### Método 2 — wall-film core-width (parcial)
Mesma ideia, limiar adaptativo por-coluna `thr = 0,5·(p15+p85)`. Numa config dá
separação **sem sobreposição** (anular 0,456 > churn 0,361; p=0,0002; o máximo do
churn fica abaixo do mínimo do anular). **Mas** a ordem **global dos 5 padrões**
sai errada (o churn cai para o MENOR e bubble fica acima de plug/slug) e o
confundimento de fluxo mássico *m* não está totalmente controlado dentro dos
padrões. Resolve o par anular×churn, não serve como vazio universal.

### Método 3 — radial-shape D1 (REFUTADO)
Descritor "centro − parede" normalizado por ensaio. Parecia dar anular > churn em
9–12 de 13 pares casados por q, **mas é brilho re-embalado**: corr(D1, brilho do
centro) = 0,95. Casando o brilho (ensaios escuros, sat=0) vira **empate**
(churn 0,204 ≈ anular 0,207); e um vazio por-linha genuinamente corrigido de
iluminação **inverte** para churn 0,52 > anular 0,47. Não usar.

### Método 4 — flip-rate temporal (o que funciona para a ESTABILIDADE)
Taxa de troca gás↔líquido **por pixel** no núcleo 2-D (espaço×tempo). Não mede
magnitude de vazio — mede caos:

- **churn 0,030 > anular 0,016** (churn troca ~2×), AUC 0,88, permutação
  run-level **p ≈ 0,001**.

Por que é principiado:
- é genuinamente **temporal** (some sob embaralhamento do tempo, que mata o
  descritor alternativo "persistência");
- **robusto ao limiar** (AUC 0,87–0,99 em Otsu/mediana/fixo-70/90/110);
- **imune a brilho** (dentro de escuros e dentro de saturados, AUC=1,00; resíduo
  após regredir brilho mantém anular < churn, p=0,002);
- **sem confundimento de q** (corr(flip, q) ≈ 0,01; o anular é amostrado em q
  MAIOR, o que joga CONTRA o achado, e mesmo assim se mantém);
- **não é artefato de α≈0,5**: normalizar pelo null de Bernoulli 2α(1−α) preserva
  a separação (anular 0,032 < churn 0,061).

> **Nuance importante:** o flip-rate por-pixel põe o **churn no TOPO** (mais
> caótico), enquanto a taxa de alternação CORDAL do v11.1 (gatilho de Schmitt,
> plug>slug>churn>anular) põe plug/slug no topo. Os dois são eixos de
> intermitência diferentes, mas **ambos põem o anular embaixo** (núcleo estável)
> — consistentes entre si.

---

## 4. Veredito honesto

**(A) A ordem anular > churn É recuperável** por um método não-cordal
principiado (bright-core), porque:
1. o núcleo de gás é **CLARO** (não escuro) — a métrica certa mede o núcleo claro;
2. a estrutura radial do anular (núcleo largo, contínuo, filme liso) difere da do
   churn (núcleo estreito, esfarrapado, com lacunas axiais) — a continuidade
   axial do núcleo também ordena anular > churn.

**(B) MAS é recuperação de ORDENAÇÃO entre padrões, não de MAGNITUDE:**
- os valores (~0,35) são **cordais, em pixels**, longe do areal ~0,9 da
  literatura;
- **não há calibração px→mm nem segunda vista óptica** → o vazio areal verdadeiro
  (anular ~0,9) **não é recuperável** destas imagens por nenhum reprocessamento;
- as distribuições **se sobrepõem** na tendência central e o **n é pequeno**
  (anular 8 × churn 10) — a separação é de mediana;
- **não é um vazio universal dos 5 padrões:** dispersed bubbles sai alto espúrio
  (0,72 — o centro é coincidentemente claro; a regra "centro claro" precisaria de
  uma porta de contraste) e o plug praticamente **empata** com o anular, logo o
  anular não é DISTINTAMENTE o maior entre todos os cinco.

---

## 5. Recomendação para o TCC (a decidir com o orientador)

Duas saídas defensáveis, não excludentes:

- **(1) Separação em DOIS eixos (menor risco).** Manter o α cordal como eixo de
  *magnitude* (reconhecidamente degenerado: churn≈anular) e separar churn de
  anular pelo eixo de *estabilidade* — a taxa de troca temporal (flip-rate), o
  achado mais robusto. É o que a taxa de alternação do v11.1 já faz
  conceitualmente. Reportar churn+anular como "agrupamento de alto vazio" no eixo
  α, com a ressalva do teto de medida.

- **(2) Adicionar o eixo bright-core (honra a expectativa física anular > churn).**
  Usar o vazio bright-core do método 1 como eixo X corrigido, reportado
  **estritamente como ordenação entre padrões** (não vazio areal calibrado), com
  o diagnóstico óptico W/M como justificativa. É uma mudança maior e carrega as
  ressalvas de magnitude/sobreposição/5-padrões do §4(B).

A expectativa física (anular > churn) agora tem **base verificada** que não
existia antes. A escolha entre (1) e (2) é de ênfase e de quanto da análise do
v11.1 se quer mexer — por isso está deixada para revisão, sem alterar o v11.1.

---

## 6. Arquivos deste estudo (`src2/v11.2/`)

- `nonchordal_void_v11_2.py` — implementa os 4 métodos + controles de brilho/q +
  estatística (Mann-Whitney unilateral, AUC). Roda sobre os quadros BRUTOS
  (`ImgA*.jpg`, 80×1016) e o panorama.
- `nonchordal_void_table_v11_2.md` — tabela por método × padrão, controle de
  brilho (midpoint), continuidade axial e por-ensaio.
- `nonchordal_void_map_v11_2.png` — mapa **bright-core (magnitude) × flip-rate
  (estabilidade)**: anular = vazio alto + estável; churn = vazio alto + caótico.
- `radial_profile_v11_2.png` — perfil radial **W/M** (diagnóstico óptico:
  núcleo claro, interfaces escuras).
- `MEMO_v11_2.md` — este documento.

### Números-chave (reproduzidos pelo script, conferem com a verificação adversarial)

| comparação anular × churn | mediana | p (unilateral) | AUC |
|---|---|---|---|
| **bright-core (M1)** | anular 0,342 > churn 0,283 | 0,022 | 0,79 |
| **wall-film (M2)** | anular 0,456 > churn 0,361 | 0,0002 | 1,00 |
| **flip-rate (M4, estabilidade)** | churn 0,030 > anular 0,018 | 0,006 | 0,86 |
| α cordal v11.1 (referência, invertido) | churn 0,55 ≳ anular 0,52 | — | — |

> **Relação com o v11.1:** o v11.2 é um estudo PARALELO de métodos de vazio
> não-cordais. Não altera o v11.1 (que mantém o α cordal + taxa de alternação).
> O eixo de intermitência do v11.1 e o flip-rate do v11.2 são consistentes
> (ambos põem o anular como o mais estável).
