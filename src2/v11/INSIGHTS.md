# O que é a fração de vazio dos padrões de escoamento — explicado de forma simples

> Documento de apoio ao `v11`. Linguagem acessível, com analogias. Para o detalhe
> técnico e os números exatos, veja `void_fraction_table_v11.md` e o gráfico
> `alpha_regions_map_v11.png`.

---

## 1. O problema, em uma frase

Dentro de um tubo bem fininho (1 mm de diâmetro — mais fino que um espaguete) passa
um fluido que está **fervendo**. Parte vira vapor (gás), parte continua líquida. A
forma como o gás e o líquido se arranjam dentro do tubo recebe nomes — **plug, slug,
churn, annular** — e é isso que chamamos de **padrão de escoamento**.

A pergunta que este trabalho responde é: **quanto do tubo, em média, está ocupado
por gás em cada um desses padrões?** Essa quantidade tem um nome técnico: **fração de
vazio** (símbolo **α**, lê-se "alfa").

---

## 2. O que é "fração de vazio" (α)

Imagine olhar o tubo de frente e perguntar: *"que fatia da seção está preenchida por
gás, e não por líquido?"*

- α = 0 → só líquido (nenhum vazio/gás).
- α = 1 → só gás.
- α = 0,5 → metade gás, metade líquido.

O nome "vazio" vem da engenharia: para o líquido, a bolha de gás é um "espaço vazio".
É o mesmo que perguntar **quanto de espuma tem no seu copo de chope**: pouca espuma =
α baixo; copo quase só espuma = α alto.

---

## 3. Como medimos isso sem sensor, só com imagem

Não temos um medidor de vazio dentro do tubo. Temos **filmagens** do escoamento. O
truque foi transformar o vídeo em uma única imagem comprida, o **panorama**.

**Analogia do panorama:** pense numa foto de "rastro de estrelas" — a câmera fica
aberta e o céu girando vira riscos. Aqui é parecido: empilhamos lado a lado uma
coluna de cada quadro do vídeo. O resultado é uma tira longa onde:

- a **altura** da tira = a seção do tubo (de uma parede à outra);
- a **largura** = o **tempo** passando (esquerda = começo, direita = fim).

Nessa tira, o **líquido aparece claro** e o **gás aparece escuro** (é assim que a luz
atravessa cada fase). Então medir vazio vira "contar quanto da imagem está escuro".

### O passo técnico: "limiar de Otsu"
Precisamos decidir, pixel a pixel, o que é "escuro o bastante para ser gás". Fazemos
isso automaticamente com um método clássico de processamento de imagem chamado
**limiar de Otsu**: ele olha o histograma de brilho da imagem e acha sozinho o ponto
de corte que melhor separa "claro" de "escuro" — sem ninguém escolher o número na mão.

Depois: para cada **coluna** (= um instante no tempo) contamos a fração de pixels
escuros. A **fração de vazio α do ensaio** é a **média** disso ao longo de todo o tempo.

> **Ressalva honesta (importante):** o que medimos é a fração de vazio ao longo de uma
> **linha vertical no centro** do tubo (uma "corda"), não da área inteira da seção. É
> como estimar quanta gordura tem um bife olhando só o corte do meio: dá para
> **comparar** bifes entre si, mas o número absoluto não é a média exata da peça
> toda. Por isso **comparar os padrões entre si é válido**, mas o valor absoluto de α
> não deve ser lido como "fração de área da seção". Além disso, está em pixels, sem
> conversão para milímetros.

---

## 4. Os dois números do gráfico

O gráfico `alpha_regions_map_v11.png` coloca cada ensaio como um ponto, com dois eixos:

### Eixo X — α (fração de vazio média): *"quanto gás, em média?"*
Já explicado acima. Mais à direita = mais gás.

### Eixo Y — intermitência: *"o gás passa de forma estável ou aos trancos?"*
Tecnicamente é o **desvio-padrão** da fração de vazio ao longo do tempo — ou seja,
**o quanto o valor instantâneo balança** em vez de ficar parado numa média.

**Analogia do chuveiro:** dois chuveiros podem jogar a **mesma quantidade** de água
por minuto (mesmo "α médio"), mas um joga um jato **constante** e o outro fica
**pingando-jorrando-pingando**. A quantidade média é igual; o que muda é o *jeito*.
A intermitência mede exatamente esse "liga-desliga".

- **Intermitência alta** = o gás passa aos trancos (bolhões grandes separados por
  líquido). É a cara do **plug** e do **slug**.
- **Intermitência baixa** = fluxo de gás mais estável e contínuo. Mais típico do
  **churn**.

---

## 5. Os quatro padrões, em linguagem de torneira

Em ordem de "pouco gás" para "muito gás":

| Padrão | Como é fisicamente | Analogia |
|---|---|---|
| **plug** (tampão) | bolhas de gás alongadas, grandes, separadas por bastante líquido | bolhões de ar subindo num canudo de refrigerante |
| **slug** (pistão) | bolhas ainda maiores, quase do tamanho do tubo, "empurrando" o líquido | pistões de gás alternando com rolhas de líquido |
| **churn** (agitado) | caótico, gás e líquido todos misturados e revoltos | máquina de lavar na centrifugação |
| **annular** (anular) | gás corre rápido pelo centro, líquido forma um anel/filme na parede | um canudo com ar no meio e as paredes molhadas |

---

## 6. O que encontramos (as conclusões)

Analisamos **26 ensaios** (filmagens), ficando só com os quatro padrões "puros":
**plug (4), slug (3), churn (10), annular (9)**. Descartamos filmagens de outros
padrões e as que tinham **dois padrões ao mesmo tempo** (transições), porque não dá
para atribuí-las a uma única classe.

### Conclusão 1 — a fração de vazio **cresce** de plug → slug → churn
As medianas de α ficaram:

| Padrão | α mediana | Leitura |
|---|---|---|
| plug | **0,38** | menos gás |
| slug | **0,43** | |
| churn | **0,55** | mais gás |
| annular | **0,52** | |

Isso faz sentido físico: quanto mais a água ferve, mais vapor se forma, e o tubo vai
ficando "mais gasoso". **Plug tem menos gás; churn tem mais.**

### Conclusão 2 — churn e annular têm **quase a mesma quantidade** de gás
Olhe de novo: churn 0,55 e annular 0,52 — praticamente empatados. Ou seja, **o que
diferencia churn de annular não é *quanto* gás tem, e sim *como* ele está arranjado**
(bagunçado no churn; organizado em anel no annular). É um resultado interessante: a
fração de vazio sozinha **não** separa esses dois padrões.

### Conclusão 3 — a **intermitência** separa plug/slug de churn
No gráfico, plug e slug ficam **em cima** (intermitência alta: ~0,20–0,29) e churn
fica **mais embaixo** (~0,15–0,20). Faz sentido: plug e slug são o escoamento "aos
trancos" (bolhões passando), enquanto o churn é mais contínuo (só que bagunçado).
A intermitência **não** separa bem churn de annular — os dois se sobrepõem.

### Conclusão 4 — os grupos no gráfico são **rótulos conhecidos**, não "descoberta"
As manchas coloridas no gráfico são as **regiões** de cada padrão, desenhadas em volta
dos pontos que **já sabíamos** a que padrão pertenciam (a classificação veio da
inspeção visual das filmagens). **Não** rodamos um algoritmo que "descobre" grupos
sozinho. Isso é uma escolha de honestidade: dizer "olha como os padrões conhecidos se
posicionam" é diferente de afirmar "o computador achou os grupos" — e seria circular
dizer a segunda coisa desenhando em cima da resposta que já tínhamos.

---

## 7. Uma medição que tivemos que descartar (e por quê)

Um dos 26 ensaios, o **m03q40**, deu um resultado estranho: a **menor** fração de vazio
de todas (α = 0,275), como se fosse o que tem *menos* gás — justo sendo um annular, que
deveria ter bastante.

Investigando, a causa é um **problema da imagem, não do escoamento**: nesse panorama
**34% dos pixels estão "estourados" de branco** (saturados no valor máximo), e a imagem
toda é muito mais clara que as demais (brilho médio 198, contra ~65 de um annular
normal).

**Analogia da foto contra o sol:** quando você fotografa alguém com o sol atrás, o
fundo "estoura" de branco e você perde o contorno da pessoa. Aqui aconteceu o mesmo: o
excesso de luz **apagou o contraste** entre gás e líquido, e o método de corte (Otsu)
confundiu gás claro com líquido, contando *menos* gás do que realmente há. E não dá
para "consertar" depois: a informação que estourou em branco **foi perdida** — não
existe mais no arquivo. (Tentativas de correção já foram testadas e não funcionam.)

Por isso o m03q40 foi **retirado das estatísticas por padrão**, mas **mantido visível**
na tabela por ensaio (marcado como "saturado"), para não esconder nada. É por isso que o
annular aparece com **n = 8** (e não 9) nas médias.

---

## 8. Até onde confiar nestes números (limites honestos)

- **Amostras pequenas:** slug tem só 3 ensaios e plug só 4. Com tão poucos, a "mediana"
  é quase o valor do meio e os extremos são os próprios pontos — servem para comparar
  tendências, **não** como estatística robusta. Não fizemos teste estatístico.
- **Valor absoluto ≠ área:** como explicado, α aqui é "corda central média no tempo",
  com viés conhecido de superestimar; serve para **comparar padrões**, não como valor
  absoluto de engenharia.
- **Caso curioso m05q05:** uma mesma condição de operação (G=640, q=5) foi filmada duas
  vezes e as duas filmagens receberam rótulos **diferentes** — uma como slug, outra como
  plug. Mantivemos as duas (o dono dos dados confirmou que a segunda é mesmo um plug
  atípico), mas o leitor deve saber que esse ponto é uma fronteira.

---

## 9. Resumo em três frases

1. Medimos, a partir das filmagens, **quanto gás** há em média em cada padrão de
   escoamento (fração de vazio α) e **como esse gás varia no tempo** (intermitência).
2. **Mais fervura = mais gás**: plug < slug < churn; mas **churn e annular empatam em
   quantidade de gás** — o que os separa é o *arranjo*, não o volume.
3. Fomos **transparentes** sobre os limites: descartamos uma imagem estragada por
   excesso de luz, sinalizamos as amostras pequenas, e deixamos claro que as regiões do
   gráfico são rótulos já conhecidos, não grupos "descobertos" por algoritmo.

---

### Arquivos deste estudo (`src2/v11/`)
- `void_fraction_v11.py` — o programa que calcula tudo.
- `void_fraction_table_v11.md` — a tabela de números (por padrão e por ensaio).
- `alpha_regions_map_v11.png` — o gráfico α × intermitência com as regiões.
- `INSIGHTS.md` — este documento.
