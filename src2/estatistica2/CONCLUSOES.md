# O que a análise descobriu — explicado de forma simples

Pasta: `src2/estatistica2/` · 2026-10-03

Imagine que filmamos água fervendo dentro de um tubo. Dependendo da força do fogo e da
quantidade de água, a fervura forma **jeitos diferentes** de borbulhar. Dois deles aparecem
bastante nos nossos vídeos:

- **churn** — uma fervura bagunçada, tipo espuma agitada.
- **wavy annular** — um filme de água ondulando na parede, com gás passando no meio.

A pergunta que queríamos responder é bem simples:

> **Vídeos do mesmo jeito de ferver são mais parecidos entre si do que vídeos de jeitos
> diferentes?**

---

## A resposta: sim, mas de leve

Encontramos um "ar de família": vídeos do mesmo jeito de ferver realmente combinam mais entre si.
Mas é um parecido **fraco**, não gritante.

**Pense assim:** não é como separar foto de cachorro de foto de gato — isso é fácil. É mais como
olhar dois irmãos: dá para perceber que são da mesma família, mas eles também têm bastante coisa
diferente. É esse nível de "parecido" que achamos.

E, importante: esse ar de família **continua aparecendo mesmo quando a gente desconta** a força
do fogo e a quantidade de água. Ou seja, não é só "fogo mais forte deixa o vídeo diferente" — tem
algo do próprio jeito de ferver ali.

---

## Três coisas que estavam erradas antes (e consertamos)

Já existia uma análise mais antiga. A ideia central dela estava certa, mas tinha três falhas.
Esta rodada corrige todas, sem bagunçar o trabalho anterior.

**1. Um teste que não testava nada.**
Havia uma verificação que, na prática, olhava para o lugar errado — como procurar uma mancha num
canto da foto que a câmera nem chega a fotografar. Ela dizia "está tudo bem" sem ter olhado para
nada. Trocamos por uma verificação que de fato olha.

**2. Uma conclusão que dependia de como se fazia a conta.**
Uma afirmação da versão antiga virava do avesso dependendo do tipo de média usada. Isso é sinal
de que era coincidência, não prova. Tiramos e colocamos no lugar uma conta que **não** muda de
resultado conforme o gosto de quem calcula.

**3. Um ingrediente esquecido: a quantidade de água.**
A análise antiga descontava a força do fogo, mas esquecia de descontar a quantidade de água — que
também mexe no vídeo, e quase tanto quanto o próprio jeito de ferver. Passamos a descontar as
duas coisas. E o "ar de família" aguentou mesmo assim.

---

## A descoberta mais importante: a "emenda" dos vídeos

Para "ver o tubo inteiro", o programa cola pedacinhos de cada quadro do vídeo lado a lado — como
quem monta uma foto panorâmica no celular. Nosso maior medo era: **e se essa colagem criar um
desenho falso**, que a gente confunda com o tal "ar de família"?

Fizemos o teste decisivo e a notícia é boa:

- **A colagem NÃO cria o desenho falso que temíamos.** Fomos olhar exatamente onde ele apareceria
  e quase não há nada ali.
- Mas percebemos outra coisa: se a gente também olha os **detalhes bem miudinhos** da imagem, o
  "ar de família" fica mais forte. Só que essa parte miudinha é **confusa de interpretar** (mistura
  textura de verdade com defeitos de borda e com o fato de os vídeos terem durações diferentes).

Diante disso, escolhemos a versão **mais cautelosa**: ignorar essa parte confusa. Fica um resultado
menos chamativo, porém mais honesto — daqueles que a gente consegue defender numa banca sem gaguejar.

**Pense assim:** é como medir a altura de alguém sem contar o salto do sapato. Daria uns
centímetros a mais contando o salto, mas aí não seria a altura **dela** de verdade.

---

## Então o método serve?

**Serve — desde que a gente diga certinho para quê.** O que podemos afirmar com segurança:

> Nosso programa consegue extrair, de cada vídeo, uma "assinatura" que **se repete** e que
> **distingue de leve** os dois jeitos de ferver — e isso não é só reflexo do fogo nem da água.

E o que **não** podemos afirmar (está escrito também nos resultados):

- Não é um "adivinhador": ele não pega um vídeo novo e diz qual é o jeito de ferver.
- Essa assinatura é do **nosso jeito de processar a imagem** (câmera, luz, a colagem) — não é uma
  propriedade da física da fervura em si.
- Não temos repetições de verdade: cada combinação de fogo e água foi filmada **uma vez só**.
  "Mesmo jeito de ferver" quer dizer vídeos diferentes que caem no mesmo tipo — não o mesmo vídeo
  repetido.
- Como temos poucos vídeos, os números de "certeza estatística" são propositalmente modestos.
  Eles mostram **o tamanho do efeito**, não uma prova fechada.

---

## Por que não fizemos "mais análises"

Com só 19 vídeos e 2 grupos, encher de testes extras não é rigor — é enfeite que dá falsa
impressão de solidez. Por isso deixamos de fora testes repetidos, inteligência artificial (que
precisaria de centenas de exemplos) e contas que só fazem sentido com muito mais dados.

**Resumo em uma frase:** existe um "ar de família" fraco, real e que se repete; ele é honesto
justamente porque assumimos onde ele começa e onde ele acaba.
