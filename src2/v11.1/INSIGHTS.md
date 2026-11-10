# O que é a fração de vazio e a intermitência dos padrões — explicado de forma simples (v11.1)

> Revisão do `v11`. Mesma linguagem acessível, mesma estrutura — **mudou o eixo de
> intermitência** (ver §4 e §7). Para o detalhe técnico e os números exatos, veja
> `void_fraction_table_v11_1.md` e o gráfico `alpha_regions_map_v11_1.png`.

---

## 0. O que mudou do v11 para o v11.1 (leia primeiro)

No `v11`, o eixo vertical do gráfico media a intermitência pelo **desvio-padrão** da
fração de vazio no tempo (`std(α_j)`). Esse número colocava o **annular ACIMA do
churn** em intermitência — e isso está **fisicamente errado**: o annular é o
escoamento mais *estável* (gás correndo contínuo pelo centro), deveria ser o **menos**
intermitente. O orientador apontou esse estranhamento.

A auditoria mostrou que **não havia erro de programação** — os números do v11
reproduzem exatamente. A inversão vinha da **métrica escolhida** somada a um
**efeito de brilho**:

1. o desvio-padrão lia a **ondulação do filme** do annular (a parede molhada tremula)
   como se fosse "intermitência", mesmo o centro do tubo continuando sempre gás;
2. os ensaios annulares foram, por acaso, filmados com **mais aquecimento** (q maior);
   os annulares de baixo aquecimento saem com a imagem **quase estourada de branco**,
   o que **comprime** o desvio-padrão e bagunça a comparação.

**O v11.1 troca o eixo Y por uma medida que corresponde à definição física de
escoamento "intermitente":** com que **frequência** o escoamento **alterna** entre
"momento de gás" e "momento de líquido". Com isso a ordem física correta aparece:

> **plug > slug > churn > annular** — o annular é o **menos** intermitente.

E, a pedido, o v11.1 **inclui o padrão dispersed bubbles** (bolhas dispersas), que o
v11 tinha deixado de fora.

---

## 1. O problema, em uma frase

Dentro de um tubo bem fininho (1 mm de diâmetro — mais fino que um espaguete) passa
um fluido que está **fervendo**. Parte vira vapor (gás), parte continua líquida. A
forma como gás e líquido se arranjam dentro do tubo recebe nomes — **dispersed
bubbles, plug, slug, churn, annular** — e é isso que chamamos de **padrão de
escoamento**.

Este estudo responde a duas perguntas sobre cada padrão:
1. **quanto** do tubo, em média, está ocupado por gás? (a **fração de vazio**, α);
2. o gás passa de forma **estável** ou **aos trancos**? (a **intermitência**).

---

## 2. O que é "fração de vazio" (α) — eixo X do gráfico

Imagine olhar o tubo de frente e perguntar: *"que fatia está preenchida por gás, e não
por líquido?"*

- α = 0 → só líquido. α = 1 → só gás. α = 0,5 → metade e metade.

É o mesmo que perguntar **quanto de espuma tem no copo de chope**: pouca espuma = α
baixo; copo quase só espuma = α alto.

---

## 3. Como medimos isso sem sensor, só com imagem

Não temos um medidor de vazio dentro do tubo — temos **filmagens**. Transformamos o
vídeo numa única imagem comprida, o **panorama**: empilhamos lado a lado uma coluna de
cada quadro. A **altura** da tira = a seção do tubo; a **largura** = o **tempo**
passando. Nessa tira, o **líquido aparece claro** e o **gás escuro**.

**O passo técnico (limiar de Otsu):** para decidir, pixel a pixel, o que é "escuro o
bastante para ser gás", usamos o método clássico de **Otsu** — ele acha sozinho o
ponto de corte que melhor separa claro de escuro. Depois, para cada **coluna** (= um
instante) contamos a fração de pixels escuros → isso é **α_j**, a fração de vazio
**naquele instante**. A média de α_j no tempo todo é o **α do ensaio** (eixo X).

> **Ressalva honesta:** medimos o vazio ao longo de uma **linha vertical no centro** do
> tubo (uma "corda"), não da área inteira. Dá para **comparar** padrões entre si, mas o
> número absoluto não é a fração de área da seção. Está em pixels, sem conversão para mm.

---

## 4. O eixo Y novo — intermitência por "taxa de alternação"

Aqui está a mudança central do v11.1.

**A ideia física.** Um escoamento é "intermitente" quando ele fica **trocando** entre
um momento **cheio de gás** e um momento **cheio de líquido** — pense num **pistão** de
gás, depois uma **rolha** de líquido, depois outro pistão de gás... Quanto mais vezes
isso se repete por metro de tubo, **mais intermitente** é o escoamento. Essa é a
definição usada na literatura de escoamento bifásico (Costigan & Whalley; Jones &
Zuber).

**Como medimos.** Olhamos o sinal α_j (a fração de gás instante a instante) e marcamos
dois estados:
- **"momento de gás"** quando a corda fica quase toda gás (α_j acima de 0,75);
- **"momento de líquido"** quando a corda fica quase toda líquida (α_j abaixo de 0,25);
- no meio, mantemos o estado anterior (isso é um **gatilho de Schmitt** — evita contar
  tremidinhas pequenas como troca de estado).

Depois **contamos quantas vezes** o escoamento alterna entre esses dois estados, por
1000 px de tubo. Esse número é o **eixo Y**.

**Por que isso conserta o problema do annular?** Para contar como "troca", a corda
central tem de ir de *quase toda gás* a *quase toda líquido*. No annular, o centro é
**sempre gás** — a parede ondula, mas o centro nunca fica líquido. Então o annular
**quase não alterna** → fica lá embaixo, **como deveria**. A ondulação do filme, que
enganava o desvio-padrão, agora é **ignorada**. E como o corte é feito sobre a *fração*
α_j (não sobre o brilho do pixel), o **problema da imagem estourada de branco também
some** desse eixo.

**Analogia do semáforo.** Imagine anotar quantas vezes um semáforo troca de verde para
vermelho numa hora. Um cruzamento movimentado (plug/slug) troca o tempo todo; uma
estrada vazia com sinal sempre verde (annular) quase nunca troca. Contamos as trocas,
não o quão "amarelo piscante" está a luz.

---

## 5. Os cinco padrões, em linguagem de torneira

Em ordem de "pouco gás" para "muito gás":

| Padrão | Como é fisicamente | Analogia |
|---|---|---|
| **dispersed bubbles** | muitas bolhas **pequenas** espalhadas no líquido | refrigerante recém-servido, cheio de borbulhas finas |
| **plug** (tampão) | bolhas alongadas, grandes, separadas por bastante líquido | bolhões de ar subindo num canudo |
| **slug** (pistão) | bolhas quase do tamanho do tubo, "empurrando" o líquido | pistões de gás alternando com rolhas de líquido |
| **churn** (agitado) | caótico, gás e líquido misturados e revoltos | máquina de lavar centrifugando |
| **annular** (anular) | gás corre rápido pelo centro, líquido vira um anel na parede | canudo com ar no meio e paredes molhadas |

---

## 6. O que encontramos (as conclusões)

Analisamos **26 ensaios** (25 na estatística): **dispersed bubbles (1), plug (3),
slug (3), churn (10), annular (8 na estatística — 9 com o ensaio saturado)**.
Descartamos as filmagens com **dois padrões ao mesmo tempo** (transições) e o smooth
annular.

### Conclusão 1 — a fração de vazio **cresce** de bubble → plug → slug → churn (eixo X)
| Padrão | α mediana | Leitura |
|---|---|---|
| dispersed bubbles | **0,02** | quase só líquido |
| plug | **0,37** | pouco gás |
| slug | **0,43** | |
| churn | **0,55** | muito gás |
| annular | **0,52** | muito gás |

Quanto mais ferve, mais vapor → mais gás. **Igual ao v11** (o eixo X não mudou).

### Conclusão 2 — churn e annular têm **quase a mesma quantidade** de gás
churn 0,55 e annular 0,52 — praticamente empatados. **A fração de vazio sozinha não
separa esses dois padrões** — o que os separa é o *arranjo*, não o *volume*. (Também
igual ao v11.)

### Conclusão 3 (REVISADA) — a **intermitência** agora ordena os padrões corretamente
Taxa de alternação (mediana, por 1000 px):

| Padrão | taxa | Leitura |
|---|---|---|
| plug | **3,94** | alterna muito (bolhões bem separados) |
| slug | **2,57** | alterna bastante |
| churn | **0,90** | alterna pouco (caótico, mas o centro raramente esvazia) |
| annular | **0,21** | **quase não alterna** (gás contínuo no centro) |
| dispersed bubbles | **0,00** | nunca alterna (centro sempre líquido) |

**Esta é a correção do v11.1:** no v11, o annular aparecia *acima* do churn; agora o
annular é o **menos intermitente**, como o orientador esperava. A intermitência
**separa bem plug/slug (em cima) de churn/annular (embaixo)**.

### Conclusão 4 — o caso especial do dispersed bubbles
As bolhas dispersas também ficam lá embaixo (taxa 0), mas por um **motivo oposto** ao
annular: o annular está embaixo porque é **gás demais e estável**; o bubble está
embaixo porque é **líquido demais e estável** (bolhas pequenas nunca enchem a corda
central). No gráfico, isso os separa no eixo X: bubble no canto **inferior-esquerdo**
(pouco gás), annular no **inferior-direito** (muito gás). É **1 só ensaio** — mostrado
como referência, não como estatística.

### Conclusão 5 — os grupos no gráfico são **rótulos conhecidos**, não "descoberta"
As manchas coloridas são as **regiões** de cada padrão, desenhadas em volta de pontos
cujo rótulo **já conhecíamos** (inspeção visual das filmagens). **Não** rodamos um
algoritmo de agrupamento. (Igual ao v11.)

---

## 7. Por que o v11 dava a ordem errada (o diagnóstico, sem rodeios)

**Não foi bug.** Os números do v11 reproduzem exatamente. O problema era a **métrica**:

- `std(α_j)` mede "o quanto α_j balança". No annular, a **parede ondula** e faz α_j
  balançar — mas isso é **ondulação de filme**, não passagem de pistões. O desvio-padrão
  contava essa ondulação como "intermitência" e **inflava** o annular.
- Além disso, os annulares de baixo aquecimento vinham com imagem **quase estourada de
  branco** (saturada). Isso distorce o Otsu e **mexe** no desvio-padrão de um jeito que
  não tem a ver com física — um **confundimento de brilho/aquecimento**.

A taxa de alternação resolve os dois: (a) só conta uma troca quando a corda central
**realmente** vai de cheia-de-gás a cheia-de-líquido, então ignora a ondulação do
filme; (b) o corte é sobre a *fração* α_j, não sobre o brilho, então o estouro de
branco não afeta esse eixo.

---

## 8. Até onde confiar nestes números (limites honestos)

- **A banda [0,25 ; 0,75] é uma escolha.** Testamos quatro bandas (ver a tabela de
  sensibilidade): a ordem **plug > slug > churn > annular** se mantém para **qualquer
  banda larga**; ela só "embola" quando a banda fica estreita demais (perto de um corte
  na mediana). Reportamos isso abertamente — não escondemos a sensibilidade.
- **Eixo X em pixels, cordal.** Vale para **comparar** padrões, não como fração de área
  absoluta (mesma ressalva do v11).
- **Sem relógio físico.** O eixo horizontal do panorama é um eixo **espacial**
  reconstruído, não um relógio calibrado. Por isso a taxa é dada **por 1000 px** — e
  **não** convertemos para Hz.
- **Amostras pequenas:** slug e plug têm n=3; dispersed bubbles n=1. A ordenação fina
  entre plug e slug é pouco robusta — as afirmações seguras são **plug/slug > churn >
  annular**.
- **m03q40 descartado:** continua fora da estatística por saturação (34% de pixels
  estourados), por coerência com o v11. Observação: a saturação afeta o **eixo X**
  (α); o **eixo Y** (taxa) seria imune — mas mantivemos a exclusão para não mudar o
  critério no meio do caminho.

---

## 9. Resumo em três frases

1. **Eixo X (α):** mais fervura = mais gás — bubble < plug < slug < churn; churn e
   annular **empatam** em quantidade de gás (o que os separa é o arranjo).
2. **Eixo Y novo (taxa de alternação):** recupera a ordem física correta
   **plug > slug > churn > annular** — o annular é o **menos** intermitente, corrigindo
   a inversão do v11, que era artefato do desvio-padrão + brilho, **não** erro de código.
3. **Honestidade:** reportamos a sensibilidade da banda, a natureza cordal/relativa de
   α, a ausência de relógio físico (taxa por 1000 px, não Hz) e as amostras pequenas.

---

### Arquivos deste estudo (`src2/v11.1/`)
- `void_fraction_v11_1.py` — o programa que calcula tudo.
- `void_fraction_table_v11_1.md` — a tabela de números (por padrão, sensibilidade e por ensaio).
- `alpha_regions_map_v11_1.png` — o gráfico α × taxa de alternação com as regiões.
- `INSIGHTS.md` — este documento.

> **Relação com o v11:** o `../v11/` é mantido como registro do descritor anterior
> (`std(α_j)`). O v11.1 **substitui** o eixo de intermitência e **adiciona** o dispersed
> bubbles; o eixo X (α) e as decisões de rotulagem/qualidade são idênticos.
