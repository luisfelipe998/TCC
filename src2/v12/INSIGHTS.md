# INSIGHTS — O que é este estudo, em linguagem simples

Este documento explica, **sem jargão**, o que a análise da pasta `v12` faz, por que faz,
e o que ela mostra. Cada termo técnico é explicado na primeira vez que aparece.

---

## 1. O problema, em uma frase

Dentro de um **tubo muito fino** (1 milímetro de diâmetro — mais fino que um espaguete),
um líquido é aquecido até ferver e virar vapor. Dependendo de **quanto líquido** passa e de
**quanto calor** é aplicado, o escoamento dentro do tubo assume formatos (padrões)
diferentes. Queremos descobrir, para cada padrão, **com que velocidade o gás e o líquido
estão viajando** pelo tubo — e desenhar um "mapa" que mostre onde cada padrão vive.

---

## 2. Os cinco padrões de escoamento (a parte visual)

Imagine o tubo como uma estrada e o escoamento como o trânsito dentro dele. Conforme se
aquece mais, o "trânsito" muda de cara:

| Padrão | Analogia | O que acontece no tubo |
|--------|----------|------------------------|
| **bolha dispersa** | Refrigerante com gaseificação fraca | Bolhinhas minúsculas e dispersas no líquido; quase tudo ainda líquido |
| **plug** (tampão) | Bolhas de ar grandes e espaçadas numa mangueira de jardim | Pouco vapor: bolhas compridas e separadas boiando no líquido |
| **slug** (pistão) | Trens de bolhas, uma atrás da outra | Mais vapor: bolhas grandes como "pistões" empurrando o líquido |
| **churn** (agitado) | Água fervendo numa panela — caótico, borbulhante | Vapor e líquido se misturam de forma turbulenta e bagunçada |
| **anular ondulado** | Um cano com um filme de água ondulado na parede e vento forte no meio | Muito vapor: o gás corre pelo centro, o líquido vira um anel fino e ondulado na parede |

A sequência **bolha → plug → slug → churn → anular ondulado** é a ordem natural conforme se
ferve cada vez mais o fluido. Descartam-se somente os ensaios de transição com dois padrões
simultâneos (plug-slug, slug-churn).

---

## 3. Os dois "botões" do experimento

Cada ensaio é identificado por um nome como `m05q20`. Esse nome guarda os dois botões que
o pesquisador ajustava:

- **`m` → fluxo mássico `G`**: o quanto de fluido passa pelo tubo por segundo, por área.
  Pense no **volume de carros entrando numa estrada**. Valores: m03=380, m04=510, m05=640,
  m06=760, m07=890 (em kg por m² por segundo). Quanto maior, mais fluido empurrado.

- **`q` → fluxo de calor `q″`**: o quanto de calor é aplicado na parede do tubo.
  Pense no **fogo embaixo da panela**. Valores: q01=1, q05=5, q10=10, q20=20, q40=40,
  q80=80 (em kW por m²). Quanto maior, mais o fluido ferve.

> O fluido usado é o **R-600a (isobutano)** — o mesmo gás da sua geladeira — mantido a
> 27 °C de temperatura de ebulição.

---

## 4. O que queremos calcular: as "velocidades superficiais"

Precisamos de duas velocidades:

- **`Jg` — velocidade superficial do gás**: a velocidade que o **vapor** teria se ele
  sozinho ocupasse o tubo inteiro.
- **`Jl` — velocidade superficial do líquido**: idem, para o **líquido** sozinho.

> **Por que "superficial"?** É um truque de engenharia. Na vida real, gás e líquido
> dividem o espaço do tubo e se atrapalham. É difícil medir a velocidade de cada um
> separadamente. Então fingimos: "e se só o gás existisse, enchendo o tubo todo — quão
> rápido ele iria?". Essa velocidade imaginária é fácil de calcular e é o padrão que
> engenheiros do mundo todo usam para comparar escoamentos. É como medir o trânsito pela
> **quantidade de carros por hora** em vez de perseguir cada carro — uma medida justa e
> comparável.

---

## 5. Como chegamos nas velocidades (a parte da conta)

Tudo depende de **quanto do fluido já virou vapor**. Isso se chama **título de vapor `x`**:

> **Título de vapor `x`** = a fração do fluido que já é vapor. `x = 0` significa tudo
> líquido; `x = 1` significa tudo vapor; `x = 0,3` significa 30% virou vapor.
> Analogia: numa chaleira, `x` é o quanto da água já saiu como vapor pelo bico.

Calculamos `x` por um **balanço de energia** — a mesma lógica de "o calor que entra tem
que ir para algum lugar". Todo o calor aplicado (`q″`) é usado para evaporar o fluido, e
quanto mais fluido passa (`G`), mais difícil é evaporá-lo todo. A fórmula é:

```
x = (4 × calor × comprimento aquecido) / (fluxo mássico × diâmetro × calor de evaporação)
```

Com `x` em mãos, dividimos o fluxo total em parte gás e parte líquido, e convertemos em
velocidade usando a **densidade** de cada fase (o vapor é ~57× mais leve que o líquido,
por isso ele ocupa muito mais espaço e corre muito mais rápido):

```
Jg = G × x / (densidade do vapor)
Jl = G × (1 − x) / (densidade do líquido)
```

As densidades do R-600a (líquido ≈ 548, vapor ≈ 9,65 kg/m³) e o calor de evaporação
(≈ 327 kJ/kg) vêm de uma biblioteca científica padrão de propriedades de fluidos
(**CoolProp**), não foram "chutados".

---

## 6. Os gráficos (e como lê-los)

### Gráfico principal — `velocity_map_v12.png` (Jg × Jl)
Cada ponto é um ensaio. Posição horizontal = velocidade do **gás**; vertical = velocidade
do **líquido**. As cores/formas separam os cinco padrões, e as manchas coloridas (os
"clusters") mostram a região que cada padrão ocupa.

**O que você vê:** da esquerda para a direita, os padrões aparecem na ordem natural —
**bolha dispersa e plug** (quase nenhum vapor, à esquerda extrema), **slug** (pouco vapor),
depois **churn** (meio), depois **anular ondulado** (muito vapor, à direita). Ou seja:
**quanto mais à direita, mais o fluido ferveu**.

> **Nota:** bolha dispersa (m07q01) e dois ensaios plug (m05q01, m06q01) têm o mesmo
> q″=1 kW/m², portanto o mesmo Jg=0,336 m/s (G cancela em Jg). O que os separa é o Jl:
> bolha tem G=890, logo Jl maior (~1,62 m/s) do que os plugs (Jl~1,16–1,38 m/s).

---

## 7. Principais conclusões

1. **Os cinco padrões se organizam numa sequência clara de "quanto o fluido ferveu".**
   Bolha dispersa e plug vivem onde quase nada virou vapor; slug no início; churn no meio;
   anular ondulado onde boa parte já é vapor. O mapa confirma a física esperada.

2. **O gás viaja MUITO mais rápido que o líquido.** No padrão anular ondulado, o gás chega a ~27
   m/s (quase 100 km/h!) enquanto o líquido se arrasta a menos de 1 m/s. Faz sentido: o
   vapor é ~57× mais leve que o líquido, então precisa correr para "caber".

3. **No padrão anular ondulado o líquido DESACELERA conforme o vapor acelera.** É a assinatura do
   anel: o líquido fica espremido num filme fino na parede enquanto o gás domina o centro.

4. **Bolha dispersa fica no extremo baixo-Jg / alto-Jl.** Com q″=1 kW/m² e G=890 kg/(m²·s)
   (o maior fluxo mássico), quase nada virou vapor (x≈0,004). O alto G empurra muito líquido,
   produzindo o maior Jl de todos os ensaios (~1,62 m/s). Posição consistente com a física.

5. **[RETRATADO] Verificação cruzada com o sensor de vazão.** Uma versão anterior
   deste documento afirmava que o fluxo mássico `G` havia sido conferido contra o
   **sensor de vazão real** (`Voltage_1`) e que os dois "concordavam dentro de 4 a
   13%". **Essa alegação foi retirada:** o programa (`velocity_map_v12.py`) **não
   lê** `Voltage_1` nem qualquer arquivo de dados, não existe **nenhuma calibração
   voltagem→g/s** no repositório, e o `G` usado vem do **nome da pasta** (valor
   nominal de projeto), não de medição. A verificação **não foi implementada** e
   **não é reproduzível** com os dados atuais. Ver
   `src2/v12.1/ANALISE_ADVERSARIAL_v12.md` §4.

6. **Toda a conta foi verificada (5 de 5 checagens aprovadas):**
   propriedades do fluido, cada um dos 28 ensaios, as estatísticas por padrão, as unidades
   (dão metros por segundo, como esperado) e a álgebra da fórmula — tudo confere.

### Tabela-resumo das velocidades (m/s)

| Padrão | Gás (Jg) típico | Líquido (Jl) típico |
|--------|----------------|---------------------|
| bolha dispersa | ~0,3 (único ensaio) | ~1,62 (o mais alto) |
| plug   | ~0,3 (lento)   | ~1,4 |
| slug   | ~1,7           | ~1,1 |
| churn  | ~6,7           | ~1,1 |
| anular ondulado | ~20 (rápido) | ~0,7 (o mais lento) |

(Valores aproximados — a tabela completa, ensaio por ensaio, está em
`velocity_table_v12.md`.)

---

## 8. As três ressalvas honestas (limitações)

Nenhum estudo é perfeito. Para não enganar ninguém, aqui estão os limites deste:

1. **O título de vapor `x` é uma estimativa, não uma medida.** Assumimos que o fluido
   entra na seção aquecida já no ponto de ferver (sem "pré-aquecimento sobrando"). Na
   prática há um pré-aquecedor antes, então o `x` real pode ser um pouco menor. Tratamos
   `x` como um **limite superior** (o máximo possível). *Analogia:* é como estimar quanto
   da chaleira virou vapor assumindo que a água já estava a 100 °C ao começar — se ela
   estava mais fria, menos vapor saiu do que a conta diz.

2. **A velocidade do gás `Jg` depende só do calor, não do fluxo mássico.** Por uma
   coincidência da matemática, quando se calcula `x` desse jeito, o fluxo mássico `G` se
   **cancela** na fórmula de `Jg`. Resultado: no gráfico, ensaios com o mesmo calor ficam
   empilhados numa linha vertical — inclusive bolha dispersa (m07q01) e dois plugs
   (m05q01, m06q01), todos com q″=1 kW/m² e Jg=0,336 m/s. Quem carrega a informação de
   "quanto fluido passa" é o eixo do líquido (`Jl`).

3. **Modelo "sem escorregamento".** Assumimos que gás e líquido viajam na mesma
   velocidade média para fazer a divisão — o que não é 100% verdade (o gás costuma
   escorregar mais rápido). Então as velocidades são **ilustrativas do modelo**, ótimas
   para comparar padrões entre si, mas não devem ser lidas como a velocidade exata medida
   de cada fase.

---

## 9. Por que não usamos as imagens (diferença para o `v10`)

Uma tentativa anterior (`v10`) quis medir velocidade **filmando o escoamento** e vendo o
quanto as estruturas se moviam de um quadro para outro. O problema: não existe uma
**régua** nos dados (ninguém registrou quantos milímetros vale cada pixel da imagem), então
a "velocidade" saía em *pixels por segundo* — um número sem significado físico comparável.

O `v12` abandona as imagens e usa só os **parâmetros físicos do experimento** (fluxo
mássico e calor) com as leis da termodinâmica. Não precisa de régua nenhuma, e o resultado
sai direto em **metros por segundo** — unidades de verdade.

---

### Arquivos desta análise
- `velocity_map_v12.py` — o programa que faz tudo (com verificações automáticas embutidas).
- `velocity_map_v12.png` — gráfico principal: velocidade do gás × do líquido (5 padrões).
- `velocity_table_v12.md` — tabelas completas (ensaio a ensaio e resumo por padrão).
- `INSIGHTS.md` — este documento.
