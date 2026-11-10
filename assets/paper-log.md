# Paper Log — TCC: Classificação de Padrões de Ebulição por Wavelet e ML

**Última atualização:** 2026-10-27
**Legenda:** CRITICAL = citação obrigatória no texto | HIGH = recomendado | MEDIUM = complementar | DISCARD = descartado

---

## Auditoria dos originais e downloads — 2026-10-27

O inventário contém **35 PDFs**. A conferência bibliográfica desta revisão utilizou os PDFs originais, extração por `pdftotext` e inspeção visual das páginas relevantes. Os arquivos existentes foram preservados, inclusive o preprint de Brunton; a versão publicada recebeu um nome distinto.

| Arquivo novo | Origem do download | Tamanho (bytes) | SHA-256 |
|--------------|--------------------|----------------|---------|
| `brunton2020-ml-fluid-mechanics-published.pdf` | https://faculty.washington.edu/sbrunton/ARFM2020.pdf | 2236988 | `0644b5e6e6cbb1e8ecbbd137f5a5304946b1c999cc0be0fc210f075a2c4f84ad` |
| `hurlbert1984-pseudoreplication-experimental-design.pdf` | https://sciences.ucf.edu/biology/d4lab/wp-content/uploads/sites/23/2018/08/Hurlbert-1984.pdf | 3362874 | `bcb8b5046b4910abd853451625b386061a1af67f7c14b64094a369a6c42cb185` |

- **Brunton:** a primeira página da versão publicada confirma autores, título, periódico, ano, volume 52, páginas 477–508 e DOI `10.1146/annurev-fluid-010719-060214`. A revisão trata de aprendizado de máquina em mecânica dos fluidos em geral; não demonstra desempenho específico para os panoramas deste TCC.
- **Hurlbert:** *Pseudoreplication and the Design of Ecological Field Experiments*, Stuart H. Hurlbert, *Ecological Monographs* 54(2), 187–211 (1984), confirmado na capa e na primeira página do artigo (PDF p. 2, impressa p. 187). O artigo fundamenta a distinção entre replicação e amostras não independentes; não fornece uma solução estatística específica para pares de panoramas ou escalas wavelet.
- **Sun:** os autores corretos são **Zhiqiang Sun, Yanping Chen e Hui Gong** (PDF p. 2, impressa p. 1). Os prenomes incorretos no BibTeX foram corrigidos. As taxas de identificação são globais de mapas de padrões, e não acurácias de um classificador de aprendizado de máquina.

---

## Validação das análises v11 (fração de vazio) e v12 (velocidades superficiais) — 2026-10-18

Busca bibliográfica (rede de citações dos assets + busca livre) para validar as análises v11/v12.
**Veredito geral:** v12 fortemente suportada; v11 suportada nas conclusões relativas (ordenação plug<slug<churn; sobreposição churn≈anular; legitimidade do limiar Otsu), mas **os valores absolutos de α divergem** da literatura (correlações reportam α maior pois embutem deslizamento) — reportar α de v11 como métrica **relativa/óptica**, não como fração de área absoluta.

**Já cobertos pelos assets existentes (não baixar):** `deoliveira2021` (αh homogêneo, correlação Xu–Fang, S, Bo eq.18, velocidades de fase — sustenta F1, G2, G4), `dias2025knn` ("void fraction increases from slug to churn"; eqs. jl=G(1−x)/ρl, jv=Gx/ρv; sobreposição churn≈anular — sustenta F1, F2, F3, G2, G3), `cheng2023` (mapas Jg×Jl, def. velocidade superficial e Bo — sustenta G1). Mandhane 1974 NÃO encontrado em OA e NÃO necessário (G1 já coberto por cheng2023 + Revellin–Thome).

### PDFs NOVOS baixados nesta sessão (3 arquivos) — não estavam em assets/

| Arquivo | Bib key | Relevância | Descrição |
|---------|---------|-----------|-----------|
| `otsu1979-threshold-selection.pdf` | `otsu1979` | **CRITICAL (v11)** | Otsu (1979) IEEE Trans. SMC 9(1):62–66 — "A Threshold Selection Method from Gray-Level Histograms". Método canônico de limiarização automática (discriminante) usado em v11 para classificar gás (escuro) × líquido (claro). DOI: 10.1109/TSMC.1979.4310076 |
| `schepperle2024-cv-void-fraction-microchannel.pdf` | `schepperle2024` | **HIGH (v11)** | Schepperle, Junaid & Woias (2024) Sensors 24(11):3363 — determinação de regimes de escoamento e **fração de vazio por visão computacional (limiarização adaptativa→segmentação binária)** em ebulição convectiva em microcanal. Antecedente metodológico para F3/F4; Eq. (3), PDF p. 7: razão entre pixels rotulados como vapor e pixels da subseção da imagem. Trata-se de proporção de área projetada, não de validação volumétrica independente no tubo circular do TCC. Gold OA. DOI: 10.3390/s24113363 |
| `driftflux2023-void-horizontal-intermittent.pdf` | `driftflux2023` | **MEDIUM (v12)** | "Improved Drift Flux Void Fraction Model for Horizontal Gas-liquid Intermittent Flow" (2023) J. Applied Fluid Mechanics 16(7):1499–1510 — revisa homogêneo×drift-flux (C₀, v_gj), define Vsg/Vsl, afirma que o modelo homogêneo "neglects the slippage… doesn't give reliable results" (sustenta G3/G4 e o enquadramento intermitente plug/slug de F5); carrega a origem Zuber & Findlay (1965). DOI: 10.47176/jafm.16.07.1693 |

---

## PDFs disponíveis localmente (assets/) — 35 arquivos em 2026-10-27

| Arquivo | Bib key | Relevância | Descrição |
|---------|---------|-----------|-----------|
| `grinsted2004-wavelet-coherence-geophysical.pdf` | `grinsted2004` | **CRITICAL** | Grinsted, Moore & Jevrejeva (2004) Nonlinear Proc. Geophys. — define formalmente a WCO clássica com operador de suavização S[·] (Eq. 8); fonte correta para eq:wco-classica. Também usa λwt ≈ 1,03s para conversão escala→frequência (Morlet ω₀=6). DOI: 10.5194/npg-11-561-2004 |
| `dias2025-knn-flow-pattern.pdf` | `dias2025knn` | **CRITICAL** | Dias et al. (2025), microtubo horizontal de 1 mm com isobutano: kNN com k=15 e 93 atributos (3 estatísticos, 72 GLCM e 18 derivados de DWT bidimensional Daubechies db2, dois níveis; PDF pp. 5–6). Acurácia média no conjunto de teste de 87% para os dez modelos da validação cruzada (Tabela 3, p. 7). O rótulo de cada vídeo foi aplicado aos seus quadros; o kNN recebeu atributos extraídos das imagens. DOI: 10.1016/j.ijmultiphaseflow.2024.105122 |
| `dias2025-cnn-flow-boiling.pdf` | `dias2025cnn` | **CRITICAL** | Dias et al. (2025), CNN ResNet50 no mesmo conjunto de imagens do estudo com kNN. **Acurácia média de teste de 88,1% nas dez dobras; 93% no melhor modelo, da dobra 3** (PDF p. 3, Tabelas 1–2). A rede recebeu as imagens como entrada. Os percentuais avaliam classificação supervisionada e não são diretamente comparáveis aos coeficientes de semelhança dos panoramas deste TCC. DOI: 10.1016/j.aitf.2025.100015 |
| `deoliveira2021-flow-patterns-pressure-drop.pdf` | `deoliveira2021` | **CRITICAL** | de Oliveira et al. (2021) IJMF — análise de padrões de escoamento e flutuações de queda de pressão durante flow boiling; mesmo grupo UNISINOS. DOI: 10.1016/j.ijmultiphaseflow.2021.103793 |
| `deoliveira2019-flow-boiling-hydrocarbons-nonlinear.pdf` | `deoliveira2019` | **CRITICAL** | de Oliveira et al. (2019) IJMF — análise não-linear de frequências características e comportamento caótico de padrões durante flow boiling de hidrocarbonetos em mini-canal 1mm; R-600a, R-290, R-1270. Mesmo arranjo experimental deste trabalho. **Método: FFT/PSD + expoente de Lyapunov + expoente de Hurst + retratos de fase — NÃO usa CWT.** A CWT foi introduzida no arranjo LETEF por deoliveira2021. DOI: 10.1016/j.ijmultiphaseflow.2019.03.015 |
| `torrence1998-wavelet-analysis-guide.pdf` | `torrence1998` | **CRITICAL** | Torrence & Compo (1998) BAMS — guia prático de análise wavelet: CWT, WPS, cross-wavelet spectrum, wavelet coherence. 1.1MB, texto legível. DOI: 10.1175/1520-0477(1998)079<0061:APGTWA>2.0.CO;2 |
| `schuck-bardo-wavelet-coherence.pdf` | `schuck2019` | **CRITICAL** | Schuck Jr. & Bodmann (2019) — formulação alternativa para WPS, WCS e WCO via CWT aplicada às funções de autocorrelação/correlação cruzada (eqs. 32.14–32.17); WCO alternativa (eq. 32.17) normaliza WCS(b,a) por WPS(0,a)·WPSy(0,a), eliminando a necessidade de suavização. Exemplo sintético com três pares de realizações de pulsos gaussianos modulados `s(t)=e^{-t²}cos(6t)`, 25 amostras/s, 2048 amostras, relação sinal-ruído de 5 dB e atraso de 500 amostras entre eventos correspondentes. A Eq. 32.16 conserva o campo cruzado complexo; o estimador da Eq. 32.20 toma módulos antes da média. Para uma única realização, essa distinção não altera a coerência, mas afeta a fase do campo cruzado. **PDF escaneado: `pdftotext` retorna vazio; OCR via pytesseract funciona (15 páginas, texto legível).** Afiliações: Schuck Jr. (DELET/UFRGS), Bodmann (DEMEC/UFRGS). DOI: 10.1007/978-3-030-16077-7_32 |
| `deoliveira2013-thesis-isobutane-minitube.pdf` | `deoliveira2013` | **CRITICAL** | De Oliveira (2013) Dissertação UNISINOS — estudo experimental de transferência de calor e queda de pressão na ebulição do isobutano em mini-tubo. Documento fundacional do arranjo experimental do LETEF utilizado neste trabalho. Texto legível. |
| `dias2024-thesis-flow-pattern-isobutane.pdf` | — | **CRITICAL** | Dias (2024) Dissertação UNISINOS — caracterização de padrões de escoamento em microcanal isobutano no mesmo arranjo experimental. Trabalho imediatamente anterior a este. **PDF escaneado/OCR corrompido: pdftotext retorna texto ilegível; consultar cópia física ou PDF digitalizado pelo autor.** |
| `schmid2022-flow-pattern-map-microchannels.pdf` | `schmid2022` | **HIGH** | Schmid et al. (2022) ETFS — mapa de padrões para escoamento bifásico adiabático de CO₂ vertical (CERN); referência de mapeamento sistemático. DOI: 10.1016/j.expthermflusci.2021.110526 |
| `brunton2020-ml-fluid-mechanics.pdf` | `brunton2020` | **HIGH** | Preprint arXiv:1905.11075v3, de 04/01/2020, preservado. Sua primeira página contém campos provisórios de periódico e DOI; usar a versão publicada abaixo para conferir os metadados e a paginação da referência final. |
| `brunton2020-ml-fluid-mechanics-published.pdf` | `brunton2020` | **HIGH** | Brunton, Noack e Koumoutsakos (2020), *Machine Learning for Fluid Mechanics*, Annual Review of Fluid Mechanics 52:477–508. Revisão de métodos de aprendizado de máquina para compreender, modelar, otimizar e controlar escoamentos. Metadados conferidos na primeira página. DOI: 10.1146/annurev-fluid-010719-060214. Origem e hash registrados acima. |
| `hurlbert1984-pseudoreplication-experimental-design.pdf` | `hurlbert1984` | **HIGH** | Hurlbert (1984), *Pseudoreplication and the Design of Ecological Field Experiments*, Ecological Monographs 54(2):187–211. A definição no resumo (PDF p. 2, impressa p. 187) inclui o uso de inferência quando tratamentos não foram replicados ou réplicas não são estatisticamente independentes. Apoio conceitual para delimitar unidades independentes; não valida um teste específico para este TCC. Origem e hash registrados acima. |
| `morshed2020-slug-wavelet-flow-regime.pdf` | `morshed2020` | **HIGH** | Morshed et al. (2020), Appl. Sci.: sinais de pressão em escoamento horizontal ar/água e ar/solução não newtoniana de goma xantana. Decomposição por pacotes wavelet db4 em quatro níveis, com 16 conjuntos de coeficientes e entropia de norma (PDF p. 20). Evidenciaram características de bolhas, transição bolhas–bolhas alongadas e slug; o resumo ressalva que o sensor não forneceu dados de alta frequência e alta resolução. Não atribuir acurácia de classificação ou robustez geral ao método. DOI: 10.3390/app10093272 |
| `qiao2023-cnn-attention-flow-pattern.pdf` | `qiao2023` | **HIGH** | Qiao et al. (2023) JMSE — CNN ResNet50 com mecanismo de atenção duplo (CBAM+ECA) para padrões gás-líquido. DOI: 10.3390/jmse11040793 |
| `he2016-resnet-deep-residual.pdf` | `he2016resnet` | **HIGH** | He et al. (2016) arXiv — Deep Residual Learning (ResNet), paper fundamental para CNN/ResNet50. arXiv: 1512.03385 |
| `alzubaidi2021-deep-learning-review.pdf` | `alzubaidi2021` | **HIGH** | Alzubaidi et al. (2021) J. Big Data — review completo de deep learning: conceitos, arquiteturas CNN. DOI: 10.1186/s40537-021-00444-8 |
| `cunningham2022-knn-tutorial.pdf` | `cunningham2020` | **HIGH** | O PDF local é o preprint arXiv:2004.04523v2, de 29/04/2020, de Cunningham & Delany, *k-Nearest Neighbour Classifiers: 2nd Edition (with Python examples)*. A referência foi alinhada à versão efetivamente consultada; o nome histórico do PDF foi preservado. Seção 2.2, Eq. (9), pp. 4–5: coeficiente de Pearson, fórmula centrada e intervalo [−1,1]. |
| `lofstedt2019-haralick-texture-invariant.pdf` | `lofstedt2019` | **HIGH** | Löfstedt et al. (2019) PLOS ONE — features de textura Haralick invariantes ao nível de cinza; complementa Haralick 1973. DOI: 10.1371/journal.pone.0212110 |
| `yang2023-cv-ml-heat-transfer-microchannels.pdf` | `yang2023` | **HIGH** | Yang et al. (2023) Energies — review: visão computacional e ML para transferência de calor em microcanais. DOI: 10.3390/en16031500 |
| `suh2021-deep-learning-boiling-heat-transfer.pdf` | `suh2021` | **HIGH** | Suh et al. (2021) Sci. Reports — deep learning para predição de transferência de calor na ebulição. DOI: 10.1038/s41598-021-85150-4 |
| `xu2020-cnn-gas-oil-ect-sensors.pdf` | `xu2020` | **MEDIUM** | Xu et al. (2020) Sensors — CNN com sensores ECT para padrões gás-óleo; metodologia com sensores não-visuais. DOI: 10.3390/s20041200 |
| `he2021-ann-pump-flow-pattern.pdf` | `he2022` | **MEDIUM** | He et al. (2022) Micromachines — ANN + SMOTE para padrões em bomba centrífuga (geometria diferente). DOI: 10.3390/mi13010002 |
| `masiukiewicz2021-flow-pattern-narrow-channel.pdf` | `masiukiewicz2021` | **MEDIUM** | Masiukiewicz & Anweiler (2021) Energies — estereologia para avaliação de padrões em canal retangular estreito. DOI: 10.3390/en14113180 |
| `ruizdiaz2024-transformer-flow-pattern-vertical.pdf` | `ruizdiaz2024` | **MEDIUM** | Ruiz-Díaz et al. (2024) Inventions — transformer neural network para padrões em tubulação vertical óleo-água. Vertical: diferente do setup. DOI: 10.3390/inventions9010015 |
| `zhang2024-arxiv-flow-pattern-vertical-bp.pdf` | `zhang2024bp` | **MEDIUM** | Zhang et al. (2024) arXiv — BP neural network para padrões em tubo vertical. Metodologia simples; vertical. arXiv: 2410.13906 |
| `cheng2023-flow-patterns-microchannels.pdf` | `cheng2023` | **HIGH** | Cheng & Xia (2023), ETFS 148, 110988. DOI conferido na primeira página: 10.1016/j.expthermflusci.2023.110988. A seção 2.2 (p. 5) contém a frase Co < 0,5 para microescala, mas também associa Co < 0,57 ao domínio gravitacional e macroescala. Por essa inconsistência interna, o TCC não adotou o limiar categórico e conservou a definição do número de confinamento e a discussão da diversidade de critérios. |
| `mallat2009-wavelet-tour-signal-processing.pdf` | `mallat1999` | **CRITICAL** | A cópia local contém o prefácio **Second Edition** (p. 17 do PDF), não a terceira edição de 2009. Referência corrigida para Mallat (1999), *A Wavelet Tour of Signal Processing*, 2ª ed., Academic Press, ISBN 978-0-12-466606-1. Metadados conferidos no [catálogo da editora](https://www.sciencedirect.com/book/9780124666061/a-wavelet-tour-of-signal-processing) e na [página do autor](https://www.di.ens.fr/~mallat/book.html). O PDF tem paginação própria; evidências indicam a página do arquivo local. Nome histórico preservado. |
| `tan2015-oil-water-wavelet-sensors.pdf` | `tan2015` | **HIGH** | Tan et al. (2015), CES: sensores combinados de condutância/capacitância mediram a fração local de água em escoamento horizontal de óleo branco e água, tubo de 50 mm. A CWT forneceu mapas de energia wavelet local (LWE), dos quais foram extraídos energia e escala normalizadas para caracterizar estruturas e transições. A distribuição triangular descreve uma tendência global, não uma correspondência rigorosamente biunívoca, e pode mudar com diâmetro e propriedades dos fluidos (PDF p. 15, impressa p. 167). A descrição atual do TCC já usa CWT; o aviso anterior sobre atribuição de DWT ficou obsoleto. DOI: 10.1016/j.ces.2015.04.046 |
| `sun2012-wavelet-entropy-flow-patterns.pdf` | `sun2012` | **HIGH** | Zhiqiang Sun, Yanping Chen e Hui Gong (2012), MST: ar–água horizontal, tubo de 50 mm; flutuações de pressão diferencial através de um corpo não aerodinâmico, DWT db4 em sete níveis e entropia de norma dos detalhes d4. Mapas combinaram entropia média com vazão mássica total ou fração volumétrica de gás. **Taxas globais de identificação de 98,0% e 95,9%**, respectivamente, em 49 dados de teste; pontos nas fronteiras receberam meio erro. Não são taxas superiores a 95% em todas as classes: plug teve 92,9% e 78,6% (PDF p. 7, impressa p. 6, Tabela 2 e seção 3.4). DOI: 10.1088/0957-0233/23/12/125301 |
| `nguyen2010-wavelet-flow-pattern-discrimination.pdf` | `nguyen2010` | **HIGH** | Nguyen, Euh e Song (2010), IJMF: CWT de sinais de fração de vazio de um Impedance Void Meter, em ar–água vertical ascendente, tubo de 80 mm no KAERI. O método combinou energia wavelet local efetiva, fração média de vazio e razão comprimento efetivo/diâmetro estimada pela escala efetiva e velocidade de propagação; distinguiu bubbly discreto, clustered/cap bubbly, slug e churn (PDF p. 11, impressa p. 765). Não demonstrou por comparação experimental que a FFT falharia especificamente na separação dos dois primeiros. DOI: 10.1016/j.ijmultiphaseflow.2010.04.007 |
| `otsu1979-threshold-selection.pdf` | `otsu1979` | **CRITICAL (v11)** | PDF preservado da busca de 18/10/2026; limiarização automática por histogramas. Registro histórico e DOI na tabela daquela sessão acima. |
| `schepperle2024-cv-void-fraction-microchannel.pdf` | `schepperle2024` | **HIGH (v11)** | PDF preservado da busca de 18/10/2026; visão computacional para regimes e fração de vazio em microcanal. Registro histórico e DOI acima. |
| `driftflux2023-void-horizontal-intermittent.pdf` | `driftflux2023` | **MEDIUM (v12)** | PDF preservado da busca de 18/10/2026; modelo drift-flux para escoamento horizontal intermitente. Registro histórico e DOI acima. |
| `hobold2018-ml-boiling-regimes.pdf` | `hobold2018` | **DISCARD** | PDF mantido no acervo, embora a referência tenha sido removida do argumento central. Trata de ebulição em piscina; ver justificativa na seção de artigos descartados. |

---

## Artigos verificados — sem PDF local (acesso pago)

### Pilar 1: Padrões de escoamento em microcanais

| ID bib | Autores | Ano | Título resumido | Relevância | DOI |
|--------|---------|-----|----------------|-----------|-----|
| `wojtan2005` | Wojtan, Ursenbacher, Thome | 2005 | Investigation of flow boiling in horizontal tubes (Part I & II) | **CRITICAL** — mapa de padrões macro-canal; referência base para comparação | 10.1016/j.ijheatmasstransfer.2004.12.012 |
| `revellin2006` | Revellin et al. | 2006 | Characterization of diabatic two-phase flows in microchannels: R-134a in 0.5mm | **CRITICAL** — primeiro mapeamento sistemático em microcanal; citado por Dias 2025 | 10.1016/j.ijmultiphaseflow.2006.02.016 |
| `saisorn2008` | Saisorn, Wongwises | 2008 | Flow pattern and heat transfer in R-134a horizontal circular mini-channel | **HIGH** — dados experimentais de padrões em mini-canal com fluido refrigerante | 10.1016/j.expthermflusci.2007.09.005 |
| `deoliveira2016` | de Oliveira et al. | 2016 | Flow boiling of R-600a in 2.6mm: heat transfer | **HIGH** — mesmo grupo UNISINOS, R600a em mini-canal | 10.1016/j.ijrefrig.2016.08.001 |
| `deoliveira2017` | de Oliveira et al. | 2017 | R-600a flow boiling: pressure drop | **HIGH** — mesmo grupo UNISINOS | 10.1016/j.ijrefrig.2017.08.004 |
| `seal2021` | Seal et al. | 2021 | Machine learning classification of in-tube condensation flow patterns | **HIGH** — ML para padrões de condensação; metodologia comparável | 10.1016/j.ijmultiphaseflow.2021.103755 |

### Pilar 2: DSP e wavelet para escoamento

| ID bib | Autores | Ano | Título resumido | Relevância | DOI |
|--------|---------|-----|----------------|-----------|-----|
| `wei2020` | Wei et al. | 2020 | Gas-liquid two-phase flow pattern identification from multi-scale entropy | **HIGH** — entropia multiescala como feature para classificação | 10.1016/j.ces.2020.115504 |
| `daubechies1990` | Daubechies | 1990 | The Wavelet Transform, Time-Frequency Localization and Signal Analysis | **HIGH** — artigo fundacional da DWT com wavelets Daubechies | 10.1109/18.57199 |

### Pilar 3: ML para classificação de padrões

| ID bib | Autores | Ano | Título resumido | Relevância | DOI |
|--------|---------|-----|----------------|-----------|-----|
| `zhang2020` | Zhang et al. | 2020 | Two-phase flow regime identification based on the liquid-gas interface | **HIGH** — identificação de regime por análise da interface + SVM | 10.1007/s00348-020-03046-x |
| `haralick1973` | Haralick, Shanmugam, Dinstein | 1973 | Textural Features for Image Classification | **CRITICAL** — artigo original da GLCM; referência obrigatória para as 14 features de Haralick | 10.1109/tsmc.1973.4309314 |

---

## Artigos descartados

| Autores | Ano | Razão do descarte |
|---------|-----|------------------|
| Hobold & da Silva | 2018 | Pool boiling (NOT flow boiling) — contexto termodinâmico distinto da ebulição convectiva em microcanais; citação removida do texto por não acrescentar ao argumento central de classificação de padrões por sinais de sensores |
| Singh, Singh | 2021 | "Two phase learning technique in modular neural network" — sobre aprendizado bifásico em redes neurais (conceito diferente de two-phase flow); título enganoso |
| Liu, Lv, Jiang | 2023 | Escoamento horizontal em tubo com *swirl* (rotação forçada) — geometria diferente do setup UNISINOS |
| Polansky, Wang | 2018 | Análise POD de escoamento anular vertical — vertical, não horizontal; diferente do setup |
| Khan et al. | 2022 | Flow Regime ID in Gas-Liquid Horizontal Pipe by Deep Learning — baixa qualidade (15KB); publicado em periódico de baixo prestígio |

---

## Livros de referência e disponibilidade

| ID bib | Autores | Ano | Título | Relevância | Editora |
|--------|---------|-----|--------|-----------|---------|
| `mallat1999` | Mallat | 1999 | A Wavelet Tour of Signal Processing (2ª ed.) | **CRITICAL** — PDF local disponível sob nome histórico `mallat2009-wavelet-tour-signal-processing.pdf`; ver identificação da edição na tabela do acervo | Academic Press. ISBN: 978-0-12-466606-1 |
| `bishop2006` | Bishop | 2006 | Pattern Recognition and Machine Learning | **HIGH** — referência padrão de ML; cobre kNN, validação cruzada, métricas | Springer. ISBN: 978-0-387-31073-2 |
| `gonzalez2018` | Gonzalez, Woods | 2018 | Digital Image Processing (4ª ed.) | **HIGH** — processamento de imagem digital; histograma, escalas de cinza, GLCM | Pearson. ISBN: 978-0-13-335672-4 |
| `thome2004` | Thome | 2004 | Engineering Data Book III | **HIGH** — dados e correlações de transferência de calor em escoamento bifásico | Wolverine Tube Inc. |

---

## Resumo de cobertura por pilar

| Pilar | CRITICAL (local) | CRITICAL (sem PDF) | HIGH local | HIGH sem PDF | Total |
|-------|-----------------|-------------------|------------|-------------|-------|
| Padrões de escoamento | deoliveira2021, deoliveira2019 | wojtan2005, revellin2006 | schmid2022, brunton2020, cheng2023 | saisorn2008, deoliveira2016/17, seal2021 | 10 |
| DSP / Wavelet | torrence1998, schuck2019, mallat1999 | — | morshed2020, tan2015, sun2012, nguyen2010 | wei2020, daubechies1990 | 9 |
| ML classificação | dias2025-knn, dias2025-cnn, he2016resnet | — | qiao2023, brunton2020, alzubaidi2021, suh2021 | zhang2020, seal2021 | 9 |
| GLCM / textura | lofstedt2019 | haralick1973 | — | — | 2 |
| kNN / classificadores | cunningham2020 | — | yang2023 | — | 2 |
| Delineamento / independência | — | — | hurlbert1984 | — | 1 |
| Setup experimental | deoliveira2021, deoliveira2019, deoliveira2013-thesis, dias2024-thesis, dias2025-knn, dias2025-cnn | — | — | deoliveira2016/17 | 6 |

**Para os artigos sem PDF:** acesse via VPN UNISINOS (Portal Capes) ou solicite ao orientador — todos têm DOI verificado.

---

## Histórico de correções bibliográficas

- `dias2025knn`: autores corrigidos (Isadora → Lucas Souza Dias; Copetti adicionada)
- `dias2025cnn`: autores corrigidos (Cardoso adicionada; Oliveira corrigido)
- `deoliveira2021`: título e autores corrigidos via CrossRef (Indrusiak, de Souza, Lima adicionados)
- `deoliveira2019`: título corrigido — era "Experimental analysis of flow patterns and void fraction for flow boiling of isobutane" (esse paper não existe com esse título); correto é "Nonlinear analysis of characteristic frequencies and chaotic behavior..."
- `tan2015`: título e DOI corrigidos — o paper citado por dias2025knn e deoliveira2021 é "Characterization of oil-water two-phase pipe flow with combined conductivity/capacitance sensor and wavelet analysis"; o DOI 10.1016/j.ces.2015.04.046 está correto
- `hobold2018`: nome do primeiro autor corrigido (Gabriel M. → Gustavo M.) e título corrigido
- `cheng2023`: título e autores completamente corrigidos via CrossRef (era paper de imagem e ML, mas o DOI aponta para review de mapas de padrões)

## Conferência conceitual — título de vapor e fração de vazio (28/10/2026)

Não são sinônimos: `x` descreve participação na vazão mássica; `α`, ocupação volumétrica. Não foi necessário acrescentar fontes ou baixar PDFs.

- **deoliveira2021**, §3.1, PDF p. 8, Eqs. (3)–(5): “The local vapor quality along the test section was calculated according to Eq. (4)” e “The void fraction was evaluated with the correlation proposed by Xu and Fang (2014)”. O título foi calculado a partir da entalpia; a fração de vazio foi estimada por correlação.
- **deoliveira2021**, §4, PDF p. 14, Eq. (17), conferida visualmente: `S = u_v/u_l = x ρ_l (1−α)/[α ρ_v (1−x)]`. Trecho literal: “corresponded to vapor and liquid mean velocity in two-phase flow, respectively”. Isolando α, obtém-se `[1+((1−x)/x)(ρ_v/ρ_l)S]⁻¹`, para `0<x<1`, com parâmetros na mesma seção. Mesmo para `S=1`, densidades diferentes mantêm as grandezas distintas.
- **deoliveira2021**, §4.3, PDF p. 16: “The increase of vapor quality and, consequently, of void fraction provided plug coalescence along the flow”. A associação vale para as condições estudadas, sem identidade entre as frações.
- **deoliveira2021**, §4.4, PDF p. 17: “x between 0.35 and 0.62 with α from 0.61 to 0.85”. Os intervalos de churn foram explicitados separadamente no TCC.
- **Limite local:** a estimativa óptica condicional descreve ocupação geométrica inferida das máscaras. Não fornece as vazões de cada fase ou a razão de deslizamento e não foi convertida em título de vapor. As chamadas e os excertos estão reconciliados em `revisao-citacoes-v2.md` e `revisao-resultados-conclusao-v2.md`.
