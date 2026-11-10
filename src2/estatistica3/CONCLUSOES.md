# O que a análise dos quatro padrões descobriu — explicado de forma simples

Pasta: `src2/estatistica3/` · 2026-10-06

Antes comparávamos **dois** jeitos de ferver. Agora comparamos **quatro**:

- **churn** — fervura bagunçada, tipo espuma agitada.
- **anular ondulado** — um filme de água ondulando na parede, com gás rápido passando no meio.
- **plug** — bolhas grandes e alongadas de gás, como rolhas separadas por água.
- **slug** — golfadas de gás ainda maiores, intermitentes e lentas.

A pergunta continua a mesma:

> **Vídeos do mesmo jeito de ferver são mais parecidos entre si do que vídeos de jeitos diferentes?**

E usamos exatamente a mesma "régua" de antes — tanto que, quando olhamos só os dois padrões antigos,
o resultado bate **casa decimal por casa decimal** com a análise anterior. Só acrescentamos dois
convidados novos (plug e slug).

---

## A resposta curta: sim, e um pouco mais nítido que antes

Com quatro padrões, o "ar de família" ficou **um pouco mais visível** (o número de separação subiu de
0,18 para 0,20). Mas o mais interessante não é esse número único — é **quem se parece com quem**.

---

## A descoberta principal: nem todo par de padrões é igual

Pense em quatro irmãos. Alguns são muito diferentes entre si; outros poderiam ser confundidos. Foi
isso que encontramos:

| Comparando... | Dá para distinguir? |
|---|---|
| anular ondulado × **slug** | **Sim, com folga** — o par mais fácil de separar |
| anular ondulado × plug | Sim, bem |
| churn × slug | Mais ou menos |
| churn × anular ondulado | Fraco (como já sabíamos) |
| plug × slug | Fraco — e com poucos vídeos, pouco confiável |
| **churn × plug** | **Quase não dá** — parecem gêmeos para a nossa régua |

**Em palavras:**

- O **anular ondulado** é o "diferentão" da família: a textura rápida do filme na parede não se parece
  com nada dos regimes lentos. É o que a régua separa melhor.
- **churn e plug são praticamente indistinguíveis** pela nossa régua. A fervura agitada e as bolhas
  alongadas deixam uma "assinatura" muito parecida no nosso tipo de medição.
- **plug e slug**, apesar de serem padrões diferentes na física, também mal se separam — e aqui temos
  só 7 vídeos no total, então é a comparação em que menos devemos confiar.

---

## Um vídeo "esquisito" que vale mencionar

Um dos vídeos de **plug** (`m05q05_002`) é o patinho feio: a assinatura dele é **o oposto** de quase
todos os outros. Na figura da matriz ele aparece como uma faixa escura atravessando tudo. Não é erro —
é um vídeo genuinamente atípico. Ele é parte da razão de o grupo "plug" ser tão bagunçado por dentro e
de "churn parecer plug": o plug, como grupo, está espalhado demais.

---

## O que melhorou em relação à análise de dois padrões

**1. O efeito do padrão ficou mais limpo.** Antes, a "quantidade de água" (fluxo mássico) atrapalhava
quase tanto quanto o próprio padrão. Agora, com quatro padrões, o efeito do padrão **domina** com
clareza, e nem a força do fogo nem a quantidade de água explicam a semelhança. Isso é um ganho.

**2. As curvas de distância se separam mais.** No gráfico que compara "parecença dentro do mesmo
padrão" com "parecença entre padrões diferentes", a diferença agora é visível a olho nu — antes as
duas curvas quase se encostavam.

---

## O que continua valendo (as mesmas honestidades de antes)

- **É um "tamanho de efeito", não uma prova fechada.** Temos 26 vídeos e um vídeo por combinação de
  fogo+água — sem repetições de verdade. Os números de "certeza estatística" são modestos de propósito.
- **plug e slug têm pouquíssimos vídeos** (4 e 3). Qualquer conclusão sobre eles — principalmente
  "plug × slug" — é frágil e está marcada como tal.
- **A assinatura é do nosso jeito de processar a imagem**, não da física da fervura em si.
- **Não é um adivinhador:** a análise não pega um vídeo novo e diz qual é o padrão; ela mede o quanto
  a assinatura se repete e separa.
- **A colagem do panorama não cria o padrão falso que temíamos** (fizemos o teste de novo; continua
  carregando quase nada de energia).

---

## Uma ressalva importante para o texto do TCC

O TCC dizia que plug e slug ficavam **de fora** porque suas oscilações são lentas demais (abaixo de ~3
vezes por segundo) para a câmera, que só resolve coisas acima de ~10–14 vezes por segundo. **Isso
continua verdade — mas para o tempo.** Esta análise não olha o tempo; ela olha o **espaço** (a foto
panorâmica do tubo, medida em pixels). No espaço, a assinatura desses padrões **é** mensurável. Então
não há contradição: o argumento do tempo permanece na seção de limitações, e o que ganhamos foi medir a
repetibilidade **espacial** dos quatro padrões.

---

**Resumo em uma frase:** o "ar de família" é real, ficou um pouco mais nítido com quatro padrões, e
revela que o anular ondulado é o mais distinto de todos enquanto churn e plug são quase gêmeos para a
nossa assinatura espectral — sempre lembrando que plug e slug têm poucos vídeos e que tudo isso é
propriedade da assinatura, não da física.
