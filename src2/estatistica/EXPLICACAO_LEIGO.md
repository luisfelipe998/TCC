# Explicando para leigos — o que foi feito aqui, passo a passo

Este arquivo conta, em linguagem simples e com analogias, o que os dois programas desta pasta
(`extract_descriptors.py` e `run_stats.py`) fazem. Não é preciso saber matemática para entender.

---

## 1. O problema, numa frase

Temos vídeos de água fervendo e virando vapor dentro de um tubo. Dependendo da intensidade,
o escoamento forma **padrões** diferentes (bolhas, golfadas, "espuma", filme nas paredes...).
A pergunta é simples de dizer e difícil de responder:

> **Duas filmagens do MESMO padrão "se parecem" mais entre si do que filmagens de padrões
> DIFERENTES?**

Se a resposta for sim, existe uma "assinatura" que caracteriza cada padrão. Trabalhamos com
dois padrões que têm filmagens suficientes: **churn** (um regime caótico, "espumoso") e
**wavy annular** (um filme ondulado de líquido na parede com gás no meio). Ao todo, **19
filmagens**: 10 de churn e 9 de wavy annular.

---

## 2. De vídeo a "linha" — a imagem costurada

A câmera filma só um pedacinho do tubo. Para "ver o tubo inteiro", o programa junta os quadros
do vídeo lado a lado, como quem monta um **panorama de celular**: tira uma tira fininha de cada
quadro e cola uma na outra.

- **Analogia:** imagine fotografar um trem passando pela janela. Cada foto pega um vagão. Se
  você recorta uma faixa de cada foto e cola todas em sequência, reconstrói o trem inteiro.
- Depois, essa imagem comprida é "achatada" numa única **linha de números**: para cada posição
  ao longo do tubo, anota-se o quão claro ou escuro está (líquido x gás). É esse sinal 1D
  (uma curva que sobe e desce) que a análise usa.

Há duas formas de achatar, e usamos as duas como conferência:
- **cropped** = média da faixa central do tubo (linhas 20 a 59);
- **center** = só a linha bem no meio (linha 40).

> ⚠️ **Pegadinha importante (a "costura"):** como colamos tiras de largura parecida, a colagem
> pode criar um **padrão repetitivo artificial** — como as emendas de um papel de parede que se
> repetem a cada X centímetros. Isso NÃO vem da física da fervura; vem do jeito de montar o
> panorama. Guardamos essa suspeita para checar depois (ver passo 7, item A3).

---

## 3. A "impressão digital" de cada filmagem (o descritor)

Precisamos resumir cada filmagem num conjunto de números comparável — uma **impressão digital**.
Isso é feito em etapas:

**a) Autocorrelação.** Pergunta: "o sinal se repete a cada quantos passos?". É como bater palmas
numa sala e ouvir o eco: o padrão do eco revela o ritmo e o tamanho do ambiente. A
autocorrelação mede o quanto a curva "combina consigo mesma" quando deslizada.

**b) Transformada wavelet (CWT).** Uma ferramenta que decompõe a curva em **tamanhos de onda**
(ondulações curtas, médias, longas), dizendo quanta "energia" existe em cada tamanho.
- **Analogia:** é como um equalizador de som que mostra quanto tem de grave, médio e agudo —
  só que aqui os "tons" são tamanhos de ondulação ao longo do tubo (medidos em pixels).
- O resultado dessa etapa é o **WPS** (*Wavelet Power Spectrum*, "espectro de potência
  wavelet"): basicamente **o gráfico do equalizador** — quanto de cada tamanho de onda a
  filmagem tem.

**c) "b=0".** O WPS completo é uma tabela grande (tamanho de onda × deslocamento). Pegamos
apenas a **coluna do deslocamento zero** — o "retrato tirado de frente", sem defasagem. É o
resumo mais limpo. Chamamos de **b=0**.

**d) Whitening ("branqueamento").** Sinais naturais sempre têm "muito grave e pouco agudo"
(ondulações grandes dominam). Se não corrigirmos, essa inclinação natural esconde as diferenças
interessantes.
- **Analogia:** é como ajustar o volume para que grave, médio e agudo saiam no mesmo nível,
  revelando os detalhes finos da música que o grave abafava.
- Na prática, removemos a "tendência geral" e ficamos só com o que **destoa** dela.

**e) Cortar no mesmo tamanho (7001 px).** As filmagens têm comprimentos diferentes. Para
compará-las de forma justa, cortamos todas no tamanho da **menor filmagem de churn (7001
pixels)** — como aparar todas as fotos no mesmo tamanho antes de empilhar.
- **Custo honesto:** as filmagens de wavy annular são muito mais longas (até ~23 mil px), então
  esse corte joga fora 30% a 70% delas. Pode apagar sinal útil — e isso está declarado, não
  escondido.

No fim, cada filmagem vira um **vetor de 200 números** (sua impressão digital). Isso é o que
`extract_descriptors.py` produz e salva nos arquivos `.npz`.

> **Conferência de segurança:** usamos um atalho matemático (FFT) para calcular rápido, mas o
> programa confirma que ele dá **exatamente** o mesmo resultado que a conta lenta e original
> (diferença de 0,0000000000000143 — praticamente zero). Assim ganhamos velocidade sem mudar o
> resultado.

---

## 4. Comparando as impressões digitais (a "tabela de distâncias")

Com 19 impressões digitais, comparamos **todas com todas**, montando uma tabela de
**"quão diferentes são"**. Usamos a correlação (semelhança) e viramos em distância:

> distância = 1 − semelhança → **0 = idênticas, 1 = nada a ver.**

Tudo que vem depois sai dessa **única** tabela 19×19. (No começo o plano tinha quatro análises
diferentes; descobrimos que três eram **a mesma coisa vista de ângulos distintos** e cortamos a
repetição — menos é mais.)

---

## 5. A resposta principal (A1): o número ANOSIM "R"

**R** é um número entre −1 e +1 que responde à pergunta central:

- **R perto de 0** → filmagens do mesmo padrão não são mais parecidas que de padrões diferentes
  (sem assinatura).
- **R positivo** → mesmo-padrão é mais parecido que entre-padrões (existe assinatura).
- Quanto maior o R, mais forte a separação.

**Deu R = 0,177.** Positivo, mas modesto: existe uma assinatura, porém **fraca**.

### O "p-valor" e por que ele é fraco aqui

O **p-valor** responde: "e se esse R tiver aparecido por pura sorte?". Para checar, embaralhamos
os rótulos (churn/annular) e vemos quantas vezes a sorte produziria um R tão bom quanto o real.

Aqui mora uma dificuldade honesta. **A intensidade do aquecimento (chamada "q") quase diz
sozinha qual é o padrão:** aquecimento baixo é sempre churn, aquecimento alto é sempre wavy
annular. Então, para não "trapacear" medindo o aquecimento em vez do padrão, só podemos
embaralhar os rótulos **dentro de grupos de mesmo q**.

- **Analogia:** é como testar se um remédio funciona comparando só pessoas da mesma idade e
  peso — justo, mas sobram poucas pessoas para comparar.
- O problema: sobram tão poucas combinações válidas (**exatamente 200**) que o p-valor **nunca
  pode ser menor que 0,005** (= 1 em 200). Ou seja, por mais real que seja o efeito, a estatística
  não consegue "cravá-lo" com força. Deu **p = 0,020** — sugestivo, mas a gente trata como
  **indício descritivo, nunca como prova forte**.

---

## 6. O gráfico de contexto (A2): o ECDF

Para dar escala ao R, desenhamos a figura `a2_ecdf.png`. Ela compara as distâncias de três
grupos de pares de filmagens:

| Grupo de pares | Distância típica (mediana) | Leitura |
|---|---|---|
| **Mesmo padrão** (intra) | **0,098** | mais parecidas ✅ |
| Padrões diferentes, aquecimento parecido | 0,159 | menos parecidas |
| Padrões diferentes, aquecimento distante | 0,135 | menos parecidas |

- Pares do **mesmo padrão são os mais parecidos** (menor distância). Isso confirma o R positivo.
- Detalhe esperto: padrões diferentes com aquecimento **parecido** (0,159) ficaram até MENOS
  parecidos que com aquecimento **distante** (0,135). Se fosse o aquecimento mandando em tudo,
  esperaríamos o contrário. Isso é um bom sinal de que estamos medindo **padrão**, não só
  aquecimento.

*(ECDF é só "o jeito de desenhar como as distâncias se distribuem, da menor para a maior" —
pense numa fila ordenada do mais parecido ao menos parecido.)*

---

## 7. O teste de robustez (A3) — e a surpresa que apareceu

A ideia era: "e se aquele **padrão artificial da costura** (passo 2) for o responsável pela
assinatura?". Para testar, íamos **tapar** a faixa de tamanhos onde a costura deixa marca e ver
se o R some.

**Surpresa ao rodar:** a marca da costura acontece em ondulações **muito curtas** (≈18 a 37
pixels), mas a nossa análise, por regra fixa do método Schuck, **só enxerga ondulações a partir
de 50 pixels**. Ou seja:

- **Analogia:** é como suspeitar de um ruído agudo estragando a música, mas o seu equalizador
  nem tem a faixa de agudos — o ruído está **fora do alcance do aparelho**.
- Consequência: "tapar a faixa da costura" não tapou **nada** (0 faixas afetadas), então o teste
  A3, do jeito que foi planejado, **não testa o que queria**.

**O lado bom:** se a marca da costura está fora da janela de análise, então ela **não pode ter
inflado** o nosso R = 0,177. A assinatura que encontramos não vem da costura. **O lado a
resolver:** ainda não *provamos* isso diretamente — só mostramos que o problema mora fora do
campo de visão. Deixamos um plano de continuação em **`FOLLOWUP_A3.md`** (baixar o "piso" de 50
para 10 pixels, reprocessar e ver se o R muda).

> Durante a execução também corrigimos um errinho de unidade (o programa comparava a régua
> errada); depois de corrigir, a conclusão continuou a mesma.

---

## 8. Conclusão em uma frase

> **O jeito como processamos as imagens produz uma "assinatura" fraca, mas repetível, que
> distingue um pouquinho churn de wavy annular — e essa assinatura é uma propriedade do nosso
> MÉTODO de reconstrução, não uma afirmação sobre a física da fervura.**

E, com honestidade, o que **não** podemos afirmar: dizer "qual é o padrão" de uma filmagem nova
(não é um classificador), medir tamanhos físicos reais, nem provar que o padrão (e não o
aquecimento) é a causa.

---

## Glossário rápido

| Termo | Em miúdos |
|---|---|
| **churn / wavy annular** | Dois tipos de escoamento: churn = caótico/"espumoso"; wavy annular = filme ondulado na parede com gás no meio. |
| **q** | Intensidade do aquecimento aplicado. Aqui ela quase prediz o padrão sozinha (confunde a análise). |
| **m** | Vazão (quanto fluido passa). Não foi controlada — mais um fator que atrapalha a atribuição causal. |
| **panorama / costura** | Juntar tiras dos quadros do vídeo para reconstruir o tubo; a emenda pode criar padrão artificial. |
| **sinal / descritor** | A "linha de brilho" extraída, resumida numa impressão digital de 200 números. |
| **autocorrelação** | Mede de quanto em quanto o sinal se repete (o "eco" dele). |
| **wavelet / CWT / WPS** | Um "equalizador" que mostra quanta energia há em cada tamanho de ondulação. WPS = o gráfico resultante. |
| **b=0** | Pegar só a "coluna de frente" do WPS (sem defasagem) — o resumo mais limpo. |
| **whitening** | Nivelar o espectro para revelar detalhes finos que a tendência natural abafa. |
| **cropped / center** | Dois jeitos de achatar o panorama em linha (faixa central vs. linha do meio); servem de conferência. |
| **distância 1 − correlação** | 0 = iguais, 1 = totalmente diferentes. |
| **ANOSIM R** | Número de −1 a +1: >0 significa "mesmo padrão é mais parecido". Deu 0,177 (fraco). |
| **p-valor** | Chance de o resultado ser sorte. Aqui limitado a ≥0,005 por falta de combinações; tratado como indício. |
| **estratificar por q** | Comparar só dentro de mesmo aquecimento, para não medir aquecimento no lugar de padrão. |
| **ECDF** | Forma de desenhar a distribuição das distâncias, da menor para a maior. |
