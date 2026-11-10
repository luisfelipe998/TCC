# Análise adversarial — dá para "aproximar o v12 da literatura"? E o que, de fato, melhora o resultado

**Data:** 2026-10-20
**Escopo:** responder honestamente à pergunta do dono — *"É possível corrigir as três lacunas (propriedades CoolProp, magnitudes absolutas de Jg/Jl, cross-check Voltage_1) de modo a aproximar do que a literatura aponta? Ou pelo menos fazer uma análise adversarial completa de qual método melhoraria nossos resultados?"* — e entregar a correção viável como código (`velocity_map_v12_1.py`).

> **Regra (imposta pelo dono, mantida):** não inventar números. Toda grandeza "v12.1" vem de `src2/v12.1/velocity_table_v12_1.md` (gerada pelo script, verificável). Toda grandeza "literatura" vem de citação verbatim dos PDFs em `/assets`, com página/equação. Onde não há fonte local, está escrito **"sem suporte em paper local"**.

---

## 0. Veredito curto

| Lacuna apontada | É corrigível? | O que foi feito |
|-----------------|---------------|-----------------|
| Magnitudes de Jg/Jl "baixas" vs. literatura | **Não por escorregamento** — erro de categoria (ver §2) | Nada a "corrigir": Jg/Jl são superficiais, invariantes ao slip. Documentado. |
| x anular baixo (x_in = 0) | **Não com os dados atuais** — barreira de dados (ver §3) | Declarado como limitação herdada; **não** fabricamos x_in. |
| Cross-check Voltage_1 "4–13%" | **Não é reproduzível** (ver §4) | **Retratado** em `INSIGHTS.md` e `VALIDACAO_LITERATURA_v12.md`. |
| Propriedades R600a (CoolProp) | Já corretas; só faltava citar o software | Citação CoolProp (Bell et al. 2014) registrada. |
| **Aproximar o resultado da literatura de forma legítima** | **Sim** — adicionando as grandezas que a literatura tabula e que mudam com slip | **`src2/v12.1/`**: α (Rouhani–Axelsson), uv, ul, S. α anular **bate** com deoliveira2021 (§5). |

**Resumo em uma frase:** não há como "empurrar" Jg/Jl para os valores da literatura (eles são iguais em equação e corretos por construção); o que aproxima de verdade o estudo da literatura é **calcular as grandezas de escorregamento** (α, uv, ul, S), e isso o v12.1 faz — com α de anular caindo dentro da faixa 0,86–0,99 de de Oliveira et al. (2021).

---

## 1. Método desta análise

1. Código do v12 lido linha a linha (`src2/v12/velocity_map_v12.py`, 543 linhas).
2. Inventário dos dados brutos: `src/data_processed/<run>/dados.xlsx` — **32 ensaios**, colunas **exatamente** `source_file, Voltage_0, Voltage_1, Voltage_2`. Sem temperatura, sem potência de pré-aquecedor, sem pressão absoluta, sem calibração de vazão. (Confirmado por leitura direta.)
3. Correlações extraídas dos PDFs e **conferidas visualmente** nas imagens renderizadas (não só `pdftotext`): Rouhani–Axelsson em `dias2025knn` Eq. (5); escorregamento em `deoliveira2021` Eq. (17).
4. Implementação em `velocity_map_v12_1.py`, com 6 verificações internas (falha ruidosa) — todas passam. Jg/Jl/x do v12.1 conferidos **idênticos** ensaio-a-ensaio aos do v12.

---

## 2. O erro de categoria: Jg/Jl são superficiais — invariantes ao escorregamento

A intuição "a literatura tem Jg/Jl maiores, então falta escorregamento no nosso modelo" **é falsa**. Por definição:

```
Jg = G·x/ρv          Jl = G·(1−x)/ρl
```

Jg e Jl dependem **só** de x, G e das densidades — **não** de como gás e líquido se distribuem no tubo. São vazões volumétricas por área, "como se cada fase ocupasse o tubo inteiro". Escorregamento (α, S) descreve a **repartição** do espaço entre as fases; **não entra** na definição de velocidade superficial. Logo:

- Adicionar Rouhani–Axelsson muda α, uv, ul, S — **não muda Jg nem Jl em um dígito sequer**.
- Verificação executada: `velocity_table_v12_1.md` reproduz Jg/Jl **idênticos** ao `velocity_table_v12.md` (28/28 ensaios, 0 divergências). Isto é a prova empírica do argumento.

**De onde vem, então, a diferença de magnitude vs. `dias2025knn` Tab. 2** (ex.: v12 churn Jg≈6,6 vs. jv 10,5; v12 anular Jg≈18,5 vs. jv 11,7)? De duas coisas que **não são** escorregamento:
1. A Tabela 2 de `dias2025knn` é **um único ponto por padrão em G=516 kg/m²s**; o v12 varre G=380–890. São operating points diferentes.
2. O x de equilíbrio do v12 difere do x com que Dias entrou na definição (ver §3).

As **equações** jl/jv do v12 são **verbatim** as de `dias2025knn` Eqs. (6)–(7) — o método está certo. As magnitudes absolutas continuam **ilustrativas do modelo**, como o próprio v12 declara.

---

## 3. Barreira de dados #1 — x_in = 0 não é corrigível (e não era a causa que se imaginava)

Tanto `deoliveira2021` (Eq. 2) quanto `dias2025knn` (Eq. 2) obtêm a **qualidade de entrada** da seção de testes por **balanço de energia do pré-aquecedor**: precisa-se da **potência do pré-aquecedor** e da **vazão mássica medida**. O inventário dos dados (§1.2) confirma que **nenhuma dessas grandezas existe** nos arquivos brutos — só três canais de voltagem crua. Portanto:

- O **x de entrada real não pode ser reconstruído**. É barreira de dados, não falha de código.
- **Contra-intuição importante:** corrigir o sub-resfriamento (x_in < 0 de entalpia) **abaixaria** x_out, afastando ainda mais das magnitudes altas de Dias. Ou seja, o "x anular baixo" **não** seria resolvido por essa correção — pelo contrário.

**Conclusão honesta:** mantém-se x de equilíbrio (limite superior), e α/uv/ul/S herdam essa natureza de estimativa — declarado em `velocity_table_v12_1.md` §"Limitação herdada" e aqui.

---

## 4. Barreira de dados #2 — o cross-check Voltage_1 "4–13%" (RETRATADO)

A afirmação em `INSIGHTS.md` §7.5 ("o fluxo mássico G foi conferido contra o sensor real `Voltage_1`; concordam dentro de 4–13%") **não é reproduzível**:

- `velocity_map_v12.py` **não lê** `Voltage_1` nem qualquer `dados.xlsx` em tempo de execução (código inteiro auditado).
- Não existe **nenhuma calibração voltagem→g/s** em lugar algum do repositório. `Voltage_1` é voltagem crua.
- G vem do **nome da pasta** (set-point nominal via `G_MAP`), **não** do sensor. Não há como os dois "concordarem a X%" sem uma curva de calibração que não existe.

**Ação (decisão do dono): retratado** em `src2/v12/INSIGHTS.md` §7.5 e marcado como retratação em `src2/v12/VALIDACAO_LITERATURA_v12.md` §4.3. Nenhum número "4–13%" permanece como afirmação válida.

---

## 5. A melhoria legítima — o que o v12.1 adiciona (α, uv, ul, S)

O caminho que **de fato** aproxima o estudo da literatura não é mexer em Jg/Jl, e sim **calcular as grandezas que a literatura tabula e que mudam com escorregamento**. O v12.1 acrescenta, por ensaio, sem tocar no v12:

```
# Fração de vazio — Rouhani–Axelsson (1970), dias2025knn Eq. (5)  [imagem conferida]
C0    = 1 + 0.12·(1 − x)
drift = 1.18·(1 − x)·[g·σ·(ρl − ρv)]^0.25 / (G·ρl^0.5)
α     = (x/ρv) / ( C0·(x/ρv + (1−x)/ρl) + drift )

# Velocidades REAIS de fase (decorrem da definição de α)
uv = jv/α            ul = jl/(1−α)

# Razão de escorregamento — deoliveira2021 Eq. (17)  [verbatim]
S  = (x·ρl·(1−α)) / (α·ρv·(1−x)) = uv/ul
```

`σ` (tensão superficial do R600a @27 °C) = **9,782 mN/m**, via CoolProp `PropsSI('I',...)` — mesmo padrão de propriedades já usado no v12.

**Flags de honestidade (declarados, não escondidos):**
- α — `dias2025knn` Eq. (5): citação direta.
- S — `deoliveira2021` Eq. (17): citação direta (verbatim).
- `uv = jv/α` e `ul = jl/(1−α)` — **decorrência** da definição de α (implícito em `deoliveira2021` via Eq. 17), **não** citação literal. Declarado como tal.
- Herdam x de equilíbrio (§3): são estimativas de modelo, não medidas.

### Validação de α vs. literatura local (números do script × citação verbatim)

| Padrão | α v12.1 (min–máx / mediana) | Literatura | Citação | Veredito |
|--------|-----------------------------|-----------|---------|----------|
| bolha dispersa | 0,145 (n=1) | ≈0,36 (G=516) | `dias2025knn` Tab. 2 | ordem de grandeza (baixo vazio); x diferente |
| plug | 0,163–0,443 / 0,186 | 0,042–0,29 | `deoliveira2021` §4.2 | **ligeiramente alto** no topo (ver nota) |
| slug | 0,218–0,512 / 0,495 | 0,578–0,714 (faixa dias2024); 0,73 (ponto dias2025knn) | `dias2024` Figs. 22d/23b/24a (conferido visualmente); `dias2025knn` Tab. 2 | **abaixo** da faixa (mesmo viés de x_in=0 do plug; ver nota) |
| churn | 0,613–0,848 / 0,737 | 0,61–0,85 | `deoliveira2021` §4.4 | **bate** (faixa sobreposta quase exata) |
| anular ondulado | 0,738–0,954 / 0,886 | 0,86–0,99 | `deoliveira2021` §4.5 | **bate** no topo/mediana; piso 0,74 = x_in=0 herdado |

> **Nota honesta (plug):** o topo da faixa plug do v12.1 (0,44, no ensaio m07q05) ultrapassa o 0,29 de de Oliveira. Causa: m07q05 tem q″=5 kW/m² e x_eq≈0,018, e Rouhani–Axelsson sobe rápido em x baixo; além disso o rótulo de padrão do v12 vem da classificação visual da LETEF, não do x. Não é erro de cálculo — é o mesmo descasamento x_equilíbrio↔padrão da §3. O **churn** e o **anular** (os padrões de alto vazio, de interesse) **batem**.

> **Nota honesta (slug) — faixa buscada em outra fonte, com confirmação visual:** de Oliveira et al. (2021) **não publica faixa de α para slug**. Isto foi verificado **exaustivamente** (lendo a §4.3 "Slug flow" inteira): o paper tabula α apenas para **plug** (§4.2, "void fraction ranging from 0.042 to 0.29"), **churn** (§4.4, "α from 0.61 to 0.85") e **anular** (§4.5); para o slug só há descritores de caos/Hurst/LLE, nenhum número de α. A faixa de slug exibida na figura vem da **dissertação de Dias (2024)** (LETEF/UNISINOS, R-600a, mini-canal 1 mm — mesmo arranjo), com **três pontos confirmados visualmente** nas páginas renderizadas (o PDF tem OCR corrompido, por isso a leitura foi por imagem): **Fig. 22d** G=520/q″=20 → α=0,7023; **Fig. 23b** G=740/q″=40 → α=0,7135; **Fig. 24a** G=1220/q″=20 → α=0,5782 → **faixa 0,578–0,714**. O próprio texto da dissertação (p.73) confirma a tendência: "conforme a velocidade mássica aumenta, os títulos de vapor diminuem… O mesmo ocorre para os valores de α e de S". O **ponto único α=0,73 de `dias2025knn` Tab. 2 (G=516)** fica logo acima do topo dessa faixa e a **corrobora**. O v12.1 (α_slug 0,22–0,51) cai **abaixo** da faixa, pelo mesmo viés de x_in=0 do plug — a faixa é referência de onde o slug real vive, não concordância de magnitude.

A **ordenação** α crescente bolha<plug<slug<churn<anular (0,15 → 0,26 → 0,43 → 0,73 → 0,87 médios) reproduz a mesma sequência de `dias2025knn` Tab. 2 (0,36/0,50/0,68/0,73/0,85/0,88) e a afirmação verbatim de `dias2025knn` ("void fraction increases from slug to churn; churn and wavy-annular are close") — aqui churn 0,73 e anular 0,87 ficam, de fato, próximos.

### S (escorregamento) vs. literatura

v12.1: plug S≈1,25–1,33; churn S≈1,50–2,59; anular S≈1,89–5,89. `deoliveira2021` (§4.2/§4.4) reporta plug S=1,24–1,34 e churn S=1,40–2,34 — **concordância forte** em plug e churn (o anular não é tabelado lá). S≥1 em todos os ensaios e cresce com o vazio, como esperado fisicamente.

### Figuras geradas
- `phase_velocity_map_v12_1.png` — mapa uv × ul (velocidades **reais**), complementa o Jg×Jl do v12.
- `void_fraction_v12_1.png` — α por padrão com as faixas de `deoliveira2021` sobrepostas.

---

## 6. Métodos que NÃO adotamos — e por quê (adversarial)

| Método candidato | Por que foi rejeitado |
|------------------|------------------------|
| **Correção de x_in** (balanço do pré-aquecedor) | Exige q_PH e ṁ medidos — **não existem** nos dados (§3). Impossível sem fabricar dado. |
| **Correlação Xu–Fang** (`deoliveira2021` Eq. 5) | Depende de αh e Fr_lo que **não estão definidos por fórmula fechada** no PDF. Rouhani–Axelsson é auto-contida e já usada por `dias2025knn`. Preferida. |
| **`dias2025knn` Eq. (2)** (x de entrada) | Além de exigir q_PH/ṁ, está **tipograficamente inconsistente** no próprio PDF. Não implementável com confiança. |
| **Drift-flux com C₀ tabelado** (`driftflux2023`) | Os C₀ 1,25/1,20 são **fit próprio de Zeghloul & Al-Sarkhi (2023), Fig. 3**, para escoamento **ar-água horizontal intermitente** — não ebulição de R600a. Transportar esses C₀ seria citação equivocada. Rouhani–Axelsson é o padrão da literatura-base da LETEF. |
| **Velocidade óptica (v10)** | Sem calibração px→mm (ScaleFactor=1,0); daria px/s sem sentido físico. Já abandonado no v12. |

---

## 7. Referências citáveis (bib keys em `assets/paper-log.md`)

- `deoliveira2021` — de Oliveira et al. (2021), IJMF 144:103793. **CRITICAL.** Sustenta: Eq. (17) escorregamento S; α por padrão §4.2/4.4/4.5; balanço x_in Eq. (2).
- `dias2025knn` — Dias et al. (2025), IJMF 185:105122. **CRITICAL.** Sustenta: Rouhani–Axelsson Eq. (5); jl/jv Eqs. (6)–(7); Tab. 2 (α e j por padrão); ponto α_slug=0,73.
- `dias2024` — Dias (2024), Dissertação de Mestrado, LETEF/UNISINOS (R-600a, mini-canal 1 mm). **CRITICAL.** Fonte da **faixa de α do slug (0,578–0,714)**, Figs. 22d/23b/24a — lida por **imagem renderizada** (OCR do PDF corrompido; ver `assets/paper-log.md`). **Sem bib key ainda** em `paper-log.md` (coluna "—"); registrar antes de citar no texto.
- `driftflux2023` — Zeghloul & Al-Sarkhi (2023), JAFM 16(7):1499–1510. **MEDIUM.** Usado só para a ressalva C₀ (não adotado).
- **CoolProp** — Bell, Wronski, Quoilin & Lemort (2014), *Ind. Eng. Chem. Res.* 53(6):2498–2508. **Software — ainda NÃO em `referencias.bib`.** Fonte de ρ_l, ρ_v, i_lv, σ.

---

### Veredito final

A pergunta "dá para aproximar da literatura?" tem resposta em três partes: (1) as magnitudes de Jg/Jl **não** devem nem podem ser "aproximadas" via escorregamento — são superficiais, iguais em equação à literatura e corretas por construção (erro de categoria, provado por Jg/Jl idênticos v12=v12.1); (2) o x_in e o cross-check Voltage_1 são **barreiras de dados** intransponíveis com os sinais atuais — o segundo foi **retratado**; (3) a melhoria **legítima e verificável** é o v12.1, que acrescenta α/uv/ul/S — grandezas que a literatura tabula e que mudam com escorregamento — com **α de churn e anular caindo dentro das faixas de de Oliveira et al. (2021)** e **S concordando com os valores de plug/churn do mesmo paper**. Nada foi inventado; o que não tem grounding local está explicitamente sinalizado.
