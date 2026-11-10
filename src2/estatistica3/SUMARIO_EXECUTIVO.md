# Sumário executivo — `estatistica3`

**Pergunta.** Vídeos (ensaios) de um *mesmo* padrão de ebulição têm assinatura espectral mais
parecida entre si do que vídeos de padrões *diferentes*?

**Resposta.** Sim, de forma modesta e desigual. A repetibilidade intra-padrão existe e é mensurável
(**R = 0,20**, p = 0,030), mas a separação entre padrões é concentrada: o **anular ondulado** se
destaca de tudo que é lento, enquanto **churn e plug são quase indistinguíveis** por este descritor.
Não é um classificador — é um tamanho de efeito descritivo.

- **Base:** `src2/estatistica3/` · motor Schuck & Bodmann (`wavelet_v8.py`) intocado
- **Fonte dos números:** `stats_results3.txt` (todos os valores abaixo saem dele)
- **Figuras geradas:** `v8-estat-dmatrix.png`, `v8-estat-ecdf.png`, `v8-estat-piso.png`

---

## 1. O que foi computado — o pipeline

```
  26 ensaios (vídeos)                      extract_descriptors.py
  churn 10 · anular 9 · plug 4 · slug 3  ────────────────────────┐
  (slug-churn n=2 e padrões de 1 ensaio DESCARTADOS)             │
                                                                 ▼
  ┌──────────────────────────────────────────────────────────────────────┐
  │ DESCRITOR por ensaio  (idêntico byte-a-byte ao da estatistica2)        │
  │   1. sinal espacial do panorama  →  autocorrelação                     │
  │   2. CWT de Morlet (cmor2.0-0.9549, 200 escalas)                       │
  │   3. lê coluna de atraso nulo  b=0  do |WPS|   (perfil de potência)    │
  │   4. WHITENING por escala: remove o decaimento 1/f de fundo            │
  │      → guarda a FORMA do espectro (onde está a energia), não a amplitude│
  │   5. trunca a 7001 px (menor ensaio de churn)                          │
  │   Saída: Z (26 × 200)  — uma "assinatura" por ensaio                   │
  └──────────────────────────────────────────────────────────────────────┘
                                   │                         run_stats.py
                                   ▼
  ┌──────────────────────────────────────────────────────────────────────┐
  │ MATRIZ DE DISTÂNCIA  D = 1 − ρ(Pearson)   (26 × 26)   ← eq. (1)        │
  │ TUDO o que segue sai desta única matriz:                               │
  │                                                                        │
  │   A0  heatmap da matriz, ordenado por padrão    → v8-estat-dmatrix.png │
  │   A1  ANOSIM R GLOBAL (R de Clarke) + p         ← eq. (2)              │
  │   A1b R par-a-par (6 pares de padrões)          ← descoberta principal │
  │   A2  ECDF intra vs inter + confound            → v8-estat-ecdf.png    │
  │   A4  por que estratificar por m, não por q (pureza dos estratos)      │
  │   A3' robustez à costura do panorama (piso 50 vs 10) → v8-estat-piso   │
  └──────────────────────────────────────────────────────────────────────┘
```

**Equações centrais (as mesmas levadas ao TCC §4.8):**

```
  (1)  distância        d_ij = 1 − ρ(z_i, z_j)
  (2)  ANOSIM / Clarke   R   = ( r̄_entre − r̄_dentro ) / (M/2),   M = C(n,2)
```

---

## 2. O que foi descoberto — resultado global

| Grandeza | Valor | Leitura |
|---|---:|---|
| **R de ANOSIM global (4 grupos)** | **0,204** | separação fraca a moderada, porém real |
| p — estrat. por fluxo mássico `m` (MC, **referência**) | **0,030** | bem condicionado |
| p — estrat. por fluxo de calor `q` (exato, 18 000) | 0,019 | **mal condicionado** (estratos puros) |
| p — sem estratificação (MC) | 0,017 | contexto, sem controle de confound |
| percentil 95 de R sob rótulos embaralhados (`m`) | ≈ 0,18 | régua do acaso |

> R subiu de **0,177** (2 padrões, estatistica2) para **0,204** (4 padrões) porque plug/slug trazem
> assinaturas distintas do anular, puxando a separação média para cima.

---

## 3. A DESCOBERTA PRINCIPAL — quem separa de quem (R par-a-par)

```
  R par-a-par (0 = indistinguível … 1 = separação perfeita)

  anular × slug   +0,52  ██████████████████████████  separação mais forte
  anular × plug   +0,32  ████████████████            separação clara
  churn  × slug   +0,22  ███████████                 separação moderada
  churn  × anular +0,18  █████████                   fraca  (≡ estatistica2)
  plug   × slug   +0,11  █████▌         (n=21)        fraca · n pequeno, cautela
  churn  × plug   +0,03  █▌                           ~ INDISTINGUÍVEIS
                         └────┬────┬────┬────┬────┬───
                            0,0  0,1  0,2  0,3  0,4  0,5
```

**História física que emerge:**

- **O anular ondulado é o padrão mais "diferentão"** — filme fino ondulando rápido na parede, textura
  espacial muito distinta dos pistões/golfadas lentos. Separa-se bem de tudo que é lento.
- **slug** é o regime intermitente mais "assinável" (0,52 vs anular; 0,22 vs churn).
- **churn ≈ plug (R = 0,03)** — a assinatura espectral **não distingue** a fervura agitada do regime
  de bolhas alongadas.
- **plug × slug (0,11)** repousa em só 21 pares (7 ensaios) → leitura menos confiável; **não** ler
  como "plug e slug são distinguíveis".
- Os 6 pares são **exploratórios**, não uma família de testes com correção de múltiplas comparações.

---

## 4. Visão da matriz de distância (`v8-estat-dmatrix.png`)

```
          churn(10)   anular(9)   plug(4)  slug(3)
        ┌──────────┬──────────┬────────┬───────┐
 churn  │  claro   │          │        │       │  blocos CLAROS na diagonal
 (1-10) │ (intra↓) │          │ ~claro │       │  = semelhança intra-padrão
        ├──────────┼──────────┼────────┼───────┤
 anular │          │  claro   │ escuro │ MUITO │  anular × slug/plug ESCURO
 (11-19)│          │ (intra↓) │        │ escuro│  = mais distantes
        ├──────────┼──────────┼────────┼───────┤
 plug   │ ~claro   │  escuro  │ MISTO  │       │  ← ensaio 21 (plug m05q05_002)
 (20-23)│(≈churn)  │          │ ✚ cruz │       │    CRUZ ESCURA: atípico, anti-
        ├──────────┼──────────┼────────┼───────┤    correlacionado com quase tudo
 slug   │          │  escuro  │        │ claro │    (ρ até −0,24). Real, não erro.
 (24-26)│          │          │        │       │    Mantido e DECLARADO.
        └──────────┴──────────┴────────┴───────┘
   amarelo/claro = parecido (distância pequena) · roxo/escuro = diferente
```

O bloco `churn × plug` ser tão claro quanto a diagonal é a imagem do R = 0,03.

---

## 5. Intra vs inter (`v8-estat-ecdf.png`)

```
  ECDF (fração acumulada de pares com distância ≤ x)
  1.0 ┤                        ╭────intra (azul, n=90)
      │                   ╭────╯   mediana 0,098
  0.5 ┤          ╭───╮╭──╯    ╭─── inter (laranja, n=235)
      │      ╭──╯    ╰╯    ╭─╯     mediana 0,153
  0.0 ┤──╭──╯        ╭────╯
      └──┼─────┬─────┬─────┬─────┬──  distância 1 − Pearson
       0.0   0.1   0.2   0.3   0.4

  intra à ESQUERDA de inter  →  mesmo padrão = distâncias menores.
  Separação visível, bem mais nítida que nos 2 padrões (onde quase coincidiam).
```

| grupo | n | mediana | média |
|---|---:|---:|---:|
| intra (mesmo padrão) | 90 | 0,098 | 0,166 |
| inter (padrões distintos) | 235 | 0,153 | 0,249 |

---

## 6. É só efeito das condições de operação? (confound — A2/A4)

```
  Correlação da distância espectral com cada fator (|·| = quanto importa):

  mesmo padrão   −0,14  ███████████████   ← MAIOR magnitude (efeito do padrão)
  |Δq| f. calor  −0,10  ██████████         pequeno, sinal "errado"
  |Δm| f. massa  −0,05  █████              pequeno
```

**Com 4 padrões o efeito do padrão (−0,14) DOMINA os dois confounds operacionais** (ambos ≈ 0 e com
sinal "errado"). Quadro mais limpo que nos 2 padrões, onde `|Δm|` rivalizava com o padrão.

**Por que estratificar por `m` e não por `q`** (A4):

```
  estratos por q (fluxo de calor)        estratos por m (fluxo mássico)
    q=10,20,80 → quase só anular/churn     m=5,6 → MISTURAM os 4 padrões
    3 de 6 estratos puros/quase-puros      1 de 5 puros/quase-puros
    ⇒ q "ordena" os padrões → MAL          ⇒ há o que permutar → REFERÊNCIA
      condicionado (remove o sinal)           honesta (p = 0,030)
```

---

## 7. A costura do panorama cria padrão falso? (robustez — A3′)

Baixando o piso de 50 px → 10 px, R sobe (0,204 → 0,292). **Mas** a banda da emenda carrega só
**0,9 %** da energia e não é o pico (que fica em λ ≈ 360 px). Logo a subida **não** vem do artefato de
costura; vem de conteúdo difuso de pequena escala (colinear com velocidade). **Leitura conservadora:
mantém piso 50 px.** A costura não fabricou o resultado.

---

## 8. Alcance e limitações (declarado)

- **Tamanho de efeito, não teste forte.** 26 ensaios, 1 por célula (m, q) → sem réplicas verdadeiras.
- **plug (4) e slug (3) são minúsculos** → R par-a-par com eles é grosseiro/instável.
- **Propriedade do descritor WPS b=0 reconstruído** (câmera, luz, costura), não da física do escoamento.
- **Não é classificador** — quantifica repetibilidade de uma assinatura, não atribui ensaio novo.
- **Ressalva temporal:** a análise é **espacial** (panorama em px). A exclusão temporal de plug/slug
  (câmera 10–14 Hz vs <3 Hz dos regimes lentos) segue válida para a *dinâmica temporal*, mas **não**
  impede a repetibilidade *espacial*, que é o que estes scripts mediram.

---

## 9. Robustez do readout (nota)

O descritor lê a coluna `b=0`. A alternativa (potência média sobre todos os lags, à la Torrence)
dá r ≈ 0,98 com `b=0`, R marginalmente menor (0,176 vs 0,204) e **ordenação par-a-par idêntica** →
a conclusão não depende desse detalhe. WCS/WCO cross-ensaio **não** entram (somem sob deslocamento
circular; medem fase acidental, não repetibilidade) — ficam só como figura ilustrativa.

---

### Uma frase

> O perfil de potência espectral espacial **repete dentro de cada padrão** (R = 0,20, p = 0,030,
> robusto a confound e à costura), e essa repetibilidade **separa bem o anular rápido dos regimes
> lentos**, mas **não distingue churn de plug** — um "ar de família" modesto, não um classificador.
</content>
</invoke>
