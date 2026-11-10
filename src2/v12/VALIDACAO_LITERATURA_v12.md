# v12 — Validação na literatura local (`/assets`) do mapa de velocidades superficiais Jg × Jl

**Data:** 2026-10-20
**Escopo:** auditar se o que `src2/v12/velocity_map_v12.py` faz está **sustentado pelos papers fisicamente presentes em `/assets`**, comparar valor-a-valor com a literatura (com a citação exata para cada número) e apontar **honestamente** o que **não** tem apoio nos papers locais. Nenhum código do v12 foi alterado.

> **Regra desta auditoria (imposta pelo dono):** NÃO inventar números. Todo valor da coluna "v12" vem de `src2/v12/velocity_table_v12.md` (gerado pelo próprio script); todo valor da coluna "literatura" vem de uma **citação verbatim** extraída do PDF correspondente, com página/equação. Onde não há fonte local, está escrito "**sem suporte em paper local**".

---

## 1. Método de auditoria

1. **Código lido linha a linha** (`src2/v12/velocity_map_v12.py`, 543 linhas). As equações, constantes e contagens abaixo citam o número de linha.
2. **Cinco PDFs reextraídos verbatim** via `pdftotext` (plain e `-layout`), com localização (página / equação / tabela) de cada trecho citado:
   - `assets/deoliveira2021-flow-patterns-pressure-drop.pdf` (paper-base LETEF, mesmo arranjo)
   - `assets/dias2025-knn-flow-pattern.pdf` (trabalho imediatamente anterior, mesmo arranjo)
   - `assets/cheng2023-flow-patterns-microchannels.pdf` (revisão de mapas de padrões em microcanais)
   - `assets/schmid2022-flow-pattern-map-microchannels.pdf` (mapa jL×jV, CO₂ 8 mm)
   - `assets/driftflux2023-void-horizontal-intermittent.pdf` (modelo drift-flux, horizontal intermitente)
3. Os números do v12 foram conferidos contra `velocity_table_v12.md` (arquivo em disco), não recalculados.

---

## 2. Passo a passo do que o v12 faz — método, equação e referência de apoio

Para cada ensaio `mXXqYY`, com `A = πD²/4` e área aquecida `= πD·L`:

| Passo | O que o v12 faz | Equação / local no código | Referência local que sustenta | Veredito de grounding |
|------|------------------|---------------------------|-------------------------------|------------------------|
| **1. Entradas** | Decodifica G do prefixo `m` e q″ do prefixo `q` | `G_MAP={3:380,4:510,5:640,6:760,7:890}` kg/m²s (linha 100); q″∈{1,5,10,20,40,80} kW/m² | `deoliveira2021`, **Tabela 1 (p.4)**: "G [kg/(m² s)] 380, 510, 640, 760 and 890"; "q″TS [kW/m²] 1, 5, 10, 20, 40 and 80" | **Confirmado** (mesmo arranjo, mesmos níveis) |
| **2. Geometria / Tsat** | D=1 mm, L_aquecido=0,265 m, Tsat=27 °C | `D=1e-3`, `L_HEATED=0.265`, `T_SAT=300.15` (linhas 94–95, 93) | `deoliveira2021`, **§2 (p.3)** "1.0 mm internal diameter"; TS "heated length of 265 mm"; **Tabela 1** "Tsat [°C] 27" | **Confirmado** |
| **3. Propriedades do fluido** | ρ_l, ρ_v, i_lv do R600a a 27 °C | `PropsSI(...)` (linhas 169–177) → ρ_l=548,13; ρ_v=9,650 kg/m³; i_lv=326,71 kJ/kg; P_sat=3,716 bar | **CoolProp 8.0.0** (biblioteca de software) | **SEM suporte em paper local** — valor de biblioteca, não de artigo (ver §4, item 1) |
| **4. Título de equilíbrio** | `x_out = 4·q″·L/(G·D·i_lv)`, depois `x=clip(x_out,0,1)` | linha 204 (= `Bo·4L/D`), linha 205 | Número de ebulição `Bo = q″/(G·i_lv)` em `deoliveira2021`, **Eq. (18) (p.13)**. O fator `4L/D` é o balanço `q″·(πD·L)=ṁ·x·i_lv` com `ṁ=G·πD²/4` | **Confirmado o Bo**; o reescalonamento `4L/D` é álgebra do balanço, não um valor de tabela |
| **5. x_in = 0** | Entrada como líquido saturado (ignora sub-resfriamento do pré-aquecedor) | `clip(...,0,1)`, limitação declarada (linhas 38–40, 497–500) | — | **Hipótese do modelo** (limite superior), não um valor de literatura; a literatura (ver §4 item 2) indica que o x real é menor |
| **6. Separação de fases** | `Jg = G·x/ρ_v`; `Jl = G·(1−x)/ρ_l` | linhas 206–207 | `dias2025knn`, **Eq. (6)** `jl=G(1−x)/ρl` e **Eq. (7)** `jv=Gx/ρv` (p.4); e **`schmid2022`** (nomenclatura, verbatim): `jL=G(1−x)/ρ`, `jV=G·x/ρ` | **Confirmado por duas fontes independentes** |
| **7. Mapa Jg×Jl** | Scatter log-log Jg×Jl, clusters por padrão | `plot_velocity_map` (linhas 359–397) | `cheng2023`, **§2 (p.6)**: coordenada "superficial gas and liquid velocities Jg and Jl in (m/s)"; **Fig. 36** mapa em jg×jl; `schmid2022`, **Fig. 5/6** (eixos jL×jV) | **Confirmado como sistema de coordenadas padrão** |
| **8. Ordenação dos padrões** | bolha→plug→slug→churn→anular com Jg crescente | ordem `PATTERN_ORDER` (linha 135) | `cheng2023` (p.11): "transitions from slug flow to churn flow and from churn flow to annular flow move to higher vapor velocities"; `schmid2022` Fig.6: "vapour quality increasing … bubbly, slug, churn, annular, vapour/mist" | **Confirmado** |
| **9. Contagem / descarte** | 28 mantidos (bolha 1 / plug 3 / slug 4 / churn 10 / anular 10); descarta 3 transições | `sanity_checks` (linhas 256–258) | Rótulos de `CLASSIFICACAO_PADRAOES_MANUAL.md` (interno ao projeto) | Interno ao projeto (não é número de literatura) |

> **Observação de método, não de valor:** o paper-base `deoliveira2021` **NÃO constrói um mapa Jg×Jl** — ele classifica por `β = S·Bo·10⁵` (**Eq. 19**, cortes 0,5 / 10 / 25 / 152 / 519). Logo a *metodologia do mapa* do v12 apoia-se em `cheng2023`/`dias2025knn`/`schmid2022`, e não no paper flagship do LETEF (ver §4, item 4).

---

## 3. Comparação valor-a-valor (v12 × literatura)

Números do v12 retirados de `src2/v12/velocity_table_v12.md`. Números da literatura são citações verbatim dos PDFs.

### 3.1 Parâmetros de entrada e propriedades

| Grandeza | v12 (disco) | Literatura | Citação exata | Veredito |
|----------|-------------|------------|---------------|----------|
| Faixa de G | 380–890 kg/m²s | 380, 510, 640, 760, 890 | `deoliveira2021` Tab. 1 (p.4) | **confirma** |
| Faixa de q″ | 1–80 kW/m² | 1, 5, 10, 20, 40, 80 | `deoliveira2021` Tab. 1 (p.4) | **confirma** |
| D, L_aq, Tsat | 1 mm; 0,265 m; 27 °C | 1,0 mm; TS aquecida 265 mm; 27 °C | `deoliveira2021` §2 (p.3) + Tab. 1 | **confirma** |
| ρ_l | 548,13 kg/m³ | — | CoolProp 8.0.0 (software) | **sem grounding em paper** |
| ρ_v | 9,650 kg/m³ | — | CoolProp 8.0.0 (software) | **sem grounding em paper** |
| i_lv | 326,71 kJ/kg | — | CoolProp 8.0.0 (software) | **sem grounding em paper** |

### 3.2 Título de vapor por padrão (modelo homogêneo do v12 × modelo com escorregamento da literatura)

> A literatura reporta x **com escorregamento** (modelos Rouhani–Axelsson / Xu–Fang). O v12 usa x de **equilíbrio sem escorregamento** → valores do v12 ficam **abaixo** nos padrões de alto vazio. Isto é esperado e declarado como limitação; **não é concordância de magnitude**.

| Padrão | x — v12 (faixa em disco) | x — literatura | Citação exata | Veredito |
|--------|--------------------------|----------------|---------------|----------|
| churn | 0,043–0,203 | 0,35–0,62 | `deoliveira2021` §4.4: "x between 0.35 and 0.62" | **diverge em magnitude** (v12 subestima; esperado p/ no-slip) |
| anular | 0,085–0,683 | 0,68–0,95 | `deoliveira2021` §4.5: "x varied between 0.68 and 0.95" | **diverge em magnitude** (v12 subestima; esperado) |
| churn | — (v12 não calcula α) | α 0,61–0,85 | `deoliveira2021` §4.4 | contexto (v12 não produz α) |
| anular | — (v12 não calcula α) | α 0,86–0,99 | `deoliveira2021` §4.5 | contexto |

### 3.3 Velocidades superficiais (Jg, Jl) — v12 × `dias2025knn` Tabela 2

> As **equações** (jl, jv) são idênticas às de `dias2025knn` (Eqs. 6–7) — **método confirmado**. Já as **magnitudes** divergem, porque `dias2025knn` usa α de Rouhani–Axelsson (com escorregamento) num **único G=516 kg/m²s**, enquanto o v12 usa x de equilíbrio sobre G=380–890. **Não alegar que os números batem.**

| Padrão | v12 — Jg [m/s] (média, disco) | v12 — Jl [m/s] (média, disco) | Literatura — jv [m/s] | Literatura — jl [m/s] | Citação | Veredito |
|--------|------------------------------:|------------------------------:|----------------------:|----------------------:|---------|----------|
| bolha dispersa | 0,336 | 1,6178 | 0,7 | 0,93 | `dias2025knn` Tab. 2 (linha "Bubbles" dispersas) | método igual; magnitude diferente |
| plug | 0,785 | 1,3788 | 2,4 | 0,88 | `dias2025knn` Tab. 2 | método igual; magnitude diferente |
| slug | 1,345 | 1,1394 | 3,5 | 0,86 | `dias2025knn` Tab. 2 | método igual; magnitude diferente |
| churn | 6,556 | 1,0923 | 10,5 | 0,76 | `dias2025knn` Tab. 2 | método igual; magnitude diferente |
| anular ondulado | 18,492 | 0,6943 | 11,7 | 0,73 | `dias2025knn` Tab. 2 | método igual; magnitude diferente |

> Nota de leitura: `dias2025knn` Tabela 2 é de **um único ponto por padrão** no outlet, com G=516; o v12 reporta a média sobre vários G e q″. A comparação acima é **qualitativa** (ordem de grandeza e ordenação), não ponto-a-ponto.

### 3.4 Ordenação e sentido físico (confirmados pela literatura)

| Afirmação do v12 | Literatura | Citação | Veredito |
|-------------------|-----------|---------|----------|
| Jg sobe de bolha→plug→slug→churn→anular | "transitions from slug to churn and churn to annular move to higher vapor velocities" | `cheng2023` (p.11) | **confirma** |
| Fração de vazio/vapor sobe de slug→churn e churn≈anular se sobrepõem | "the void fraction also increases from slug to churn pattern, however between churn and wavy-annular the values are close" | `dias2025knn` (após Tab. 2, p.4) | **confirma** |
| No anular, Jl cai enquanto Jg sobe (filme fino na parede) | sequência bubbly→…→annular com qualidade crescente | `schmid2022` Fig. 5/6 | **confirma o sentido** |

### 3.5 Limitação declarada "no-slip subestima x anular" — apoiada pela literatura

| Limitação do v12 (honesta) | Literatura que a sustenta | Citação | Veredito |
|-----------------------------|---------------------------|---------|----------|
| Modelo homogêneo ignora escorregamento → Jg/Jl ilustrativos; x anular sai baixo | "the homogenous model, neglects the slippage between phases, and it doesn't give reliable results" | `driftflux2023` §1 (p.1499) | **confirma** a direção do viés |
| idem | escorregamento real S>1 e crescente: plug S=1,24–1,34; churn S=1,40–2,34 | `deoliveira2021` §4.2 e §4.4 | **confirma** que S≫1 no anular/churn, logo x_equilíbrio subestima |

---

## 4. Claims SEM grounding na literatura local (honestidade)

Os itens abaixo **não** são sustentados por nenhum paper de `/assets`. Estão listados explicitamente para não serem apresentados como validados.

1. **Propriedades do R600a (ρ_l=548,13; ρ_v=9,650; i_lv=326,71 kJ/kg).** Vêm da biblioteca **CoolProp 8.0.0** (`PropsSI`), não de um artigo de `/assets`. Para citar a origem corretamente, usar a referência de software **CoolProp (Bell, Wronski, Quoilin & Lemort, 2014, *Ind. Eng. Chem. Res.* 53(6):2498–2508)** — ainda **não presente** em `referencias.bib`. (Não confundir com validação experimental: nenhum paper local tabula essas propriedades.)

2. **Magnitudes absolutas de Jg/Jl.** Divergem de `dias2025knn` Tab. 2 (ex.: v12 churn Jg≈6,6 vs. jv 10,5; v12 anular Jg≈18,5 vs. jv 11,7 m/s). A causa é de modelo (equilíbrio sem escorregamento × Rouhani–Axelsson com escorregamento; faixa de G diferente). **A literatura sustenta o método e a ordenação, não os valores absolutos.** Reportar os Jg/Jl do v12 como **ilustrativos do modelo homogêneo**, comparáveis entre padrões, não como velocidades medidas.

3. **[RETRATADO] "Verificação cruzada com o Voltage_1 concorda dentro de 4–13%"** (afirmado em `src2/v12/INSIGHTS.md` §7, item 5). **Esta checagem não está implementada em `velocity_map_v12.py`** — o script não lê `Voltage_1` nem qualquer arquivo de dados em tempo de execução; a única verificação de vazão interna é um limite de plausibilidade `mdot(m05)=0,503 g/s ∈ [0,3; 0,7]` (linhas 223–224). Não há calibração voltagem→g/s no repositório e `G` vem do nome da pasta (valor nominal). A alegação foi **retirada** de `INSIGHTS.md` em 2026-10-20 (ver `src2/v12.1/ANALISE_ADVERSARIAL_v12.md` §4); não é reproduzível com os dados atuais.

4. **O mapa Jg×Jl não vem do paper-base.** `deoliveira2021` classifica por `β=S·Bo·10⁵` (**Eq. 19**), **não** por mapa Jg×Jl. Portanto o v12 **não** replica o método do flagship do LETEF neste ponto — ele adota um sistema de coordenadas de `cheng2023`/`dias2025knn`/`schmid2022`. Isto deve ser dito explicitamente ao posicionar o v12 em relação ao trabalho-base.

---

## 5. Ressalvas de atribuição (para não citar errado depois)

- **C₀ (drift-flux) do `driftflux2023`:** os valores plug ≈ 1,25 / slug ≈ 1,20 são **o ajuste próprio do artigo (Fig. 3)**, não de França & Lahey. Na **Tabela 1** do mesmo artigo, França & Lahey (1992) aparecem como **1,0 (plug) / 1,2 (slug)**. Se for citar 1,25/1,20, atribuir a **Zeghloul & Al-Sarkhi (2023), Fig. 3**. (O v12 não usa drift-flux; esta ressalva é só para evitar citação equivocada caso o tema seja mencionado no texto.)
- **`cheng2023`:** **não** traz a equação explícita `Jg=Qg/A`, apenas a descrição ("superficial velocities Jg, Jl in m/s"). Para a **definição por quociente**, citar a **nomenclatura de `schmid2022`** (`jL=G(1−x)/ρ`, `jV=G·x/ρ`, verbatim). Além disso, `cheng2023` **não cita Mandhane (1974)** (só Taitel–Dukler, na lista de referências) — não atribuir Mandhane a Cheng.
- **Critério de microescala:** `cheng2023` (p.4) afirma **Co < 0,5** para microescala (Kew & Cornwell). Há um ponto no texto do TCC que escreve "Co > 0,5" (desigualdade invertida) — questão **pré-existente do texto**, não introduzida pelo v12, mas registrada aqui para correção futura.

---

## 6. Referências citáveis (bib keys já em `assets/paper-log.md`)

- `deoliveira2021` — de Oliveira et al. (2021), IJMF 144:103793. **CRITICAL.** Sustenta: faixas G/q″, geometria, Tsat, Bo (Eq.18), escorregamento S, x por padrão (churn/anular).
- `dias2025knn` — Dias et al. (2025), IJMF 185:105122. **CRITICAL.** Sustenta: equações jl/jv (Eqs. 6–7), ordenação (vazio sobe slug→churn, churn≈anular), Tabela 2.
- `cheng2023` — Cheng & Xia (2023), ETFS. **HIGH.** Sustenta: coordenada Jg×Jl como padrão, ordenação com velocidade crescente.
- `schmid2022` — Schmid et al. (2022), ETFS 131:110526. **HIGH.** Sustenta: definição por quociente de jL/jV (verbatim), eixos jL×jV, ordenação por qualidade crescente. (Advertência: CO₂, vertical up/down, 8 mm, adiabático — geometria/fluido diferentes; usar para o **conceito de coordenadas**, não para valores.)
- `driftflux2023` — Zeghloul & Al-Sarkhi (2023), JAFM 16(7):1499–1510. **MEDIUM.** Sustenta: modelo homogêneo "neglects slippage… doesn't give reliable results".
- **CoolProp** — Bell et al. (2014), *Ind. Eng. Chem. Res.* 53(6):2498–2508. **Software; ainda NÃO em `referencias.bib`.** Fonte das propriedades do R600a.

---

### Veredito final da auditoria

O **método** do v12 (balanço de energia → título de equilíbrio → velocidades superficiais → mapa Jg×Jl → ordenação dos padrões) está **sustentado verbatim** pela literatura local (`deoliveira2021`, `dias2025knn`, `cheng2023`, `schmid2022`, `driftflux2023`). As **magnitudes absolutas de Jg/Jl e de x** **não** são corroboradas pela literatura (diferença de modelo: equilíbrio sem escorregamento), o que o próprio v12 declara como limitação — e a literatura confirma a *direção* desse viés. Três itens **não têm grounding local** e estão sinalizados: as propriedades do R600a (CoolProp, não um paper), as magnitudes absolutas (modelo diferente) e, sobretudo, a alegação de "cross-check Voltage_1 em 4–13%" que o **código não implementa** — esta última foi **retratada** em 2026-10-20 (ver item 3 acima e `src2/v12.1/ANALISE_ADVERSARIAL_v12.md`). A melhoria legítima rumo à literatura — fração de vazio α, velocidades reais de fase e escorregamento S — está implementada em `src2/v12.1/`, com α de churn/anular caindo nas faixas de de Oliveira et al. (2021).
