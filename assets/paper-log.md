# Paper Log — TCC: Classificação de Padrões de Ebulição por Wavelet e ML

**Última atualização:** 2026-05-28
**Legenda:** CRITICAL = citação obrigatória no texto | HIGH = recomendado | MEDIUM = complementar | DISCARD = descartado

---

## PDFs disponíveis localmente (assets/) — 29 arquivos

| Arquivo | Bib key | Relevância | Descrição |
|---------|---------|-----------|-----------|
| `grinsted2004-wavelet-coherence-geophysical.pdf` | `grinsted2004` | **CRITICAL** | Grinsted, Moore & Jevrejeva (2004) Nonlinear Proc. Geophys. — define formalmente a WCO clássica com operador de suavização S[·] (Eq. 8); fonte correta para eq:wco-classica. Também usa λwt ≈ 1,03s para conversão escala→frequência (Morlet ω₀=6). DOI: 10.5194/npg-11-561-2004 |
| `dias2025-knn-flow-pattern.pdf` | `dias2025knn` | **CRITICAL** | Dias et al. (2025) kNN k=15 + 93 features (3 stat + 72 GLCM + 18 DWT) em microtubo 1mm isobutano UNISINOS, 10-fold CV, 87% acc. Paper-base do TCC. DOI: 10.1016/j.ijmultiphaseflow.2024.105122 |
| `dias2025-cnn-flow-boiling.pdf` | `dias2025cnn` | **CRITICAL** | Dias et al. (2025) CNN ResNet50 para mesmos padrões de escoamento, 93% acc. Comparativo direto com abordagem proposta. DOI: 10.1016/j.aitf.2025.100015 |
| `deoliveira2021-flow-patterns-pressure-drop.pdf` | `deoliveira2021` | **CRITICAL** | de Oliveira et al. (2021) IJMF — análise de padrões de escoamento e flutuações de queda de pressão durante flow boiling; mesmo grupo UNISINOS. DOI: 10.1016/j.ijmultiphaseflow.2021.103793 |
| `deoliveira2019-flow-boiling-hydrocarbons-nonlinear.pdf` | `deoliveira2019` | **CRITICAL** | de Oliveira et al. (2019) IJMF — análise não-linear de frequências características e comportamento caótico de padrões durante flow boiling de hidrocarbonetos em mini-canal 1mm; R-600a, R-290, R-1270. Mesmo arranjo experimental deste trabalho. **Método: FFT/PSD + expoente de Lyapunov + expoente de Hurst + retratos de fase — NÃO usa CWT.** A CWT foi introduzida no arranjo LETEF por deoliveira2021. DOI: 10.1016/j.ijmultiphaseflow.2019.03.015 |
| `torrence1998-wavelet-analysis-guide.pdf` | `torrence1998` | **CRITICAL** | Torrence & Compo (1998) BAMS — guia prático de análise wavelet: CWT, WPS, cross-wavelet spectrum, wavelet coherence. 1.1MB, texto legível. DOI: 10.1175/1520-0477(1998)079<0061:APGTWA>2.0.CO;2 |
| `schuck-bardo-wavelet-coherence.pdf` | `schuck2019` | **CRITICAL** | Schuck Jr. & Bodmann (2019) — formulação alternativa para WPS, WCS e WCO via CWT aplicada às funções de autocorrelação/correlação cruzada (eqs. 32.14–32.17); WCO alternativa (eq. 32.17) normaliza WCS(b,a) por WPS(0,a)·WPSy(0,a), eliminando a necessidade de suavização. Validado com pulso gaussiano modulado `s(t)=e^{-t²}cos(6t)`, 25 amostras/s, 2048 amostras, SNR=5 dB, atraso de 500 amostras entre realizações. **PDF escaneado: `pdftotext` retorna vazio; OCR via pytesseract funciona (15 páginas, texto legível).** Afiliações: Schuck Jr. (DELET/UFRGS), Bodmann (DEMEC/UFRGS). DOI: 10.1007/978-3-030-16077-7_32 |
| `deoliveira2013-thesis-isobutane-minitube.pdf` | `deoliveira2013` | **CRITICAL** | De Oliveira (2013) Dissertação UNISINOS — estudo experimental de transferência de calor e queda de pressão na ebulição do isobutano em mini-tubo. Documento fundacional do arranjo experimental do LETEF utilizado neste trabalho. Texto legível. |
| `dias2024-thesis-flow-pattern-isobutane.pdf` | — | **CRITICAL** | Dias (2024) Dissertação UNISINOS — caracterização de padrões de escoamento em microcanal isobutano no mesmo arranjo experimental. Trabalho imediatamente anterior a este. **PDF escaneado/OCR corrompido: pdftotext retorna texto ilegível; consultar cópia física ou PDF digitalizado pelo autor.** |
| `schmid2022-flow-pattern-map-microchannels.pdf` | `schmid2022` | **HIGH** | Schmid et al. (2022) ETFS — mapa de padrões para escoamento bifásico adiabático de CO₂ vertical (CERN); referência de mapeamento sistemático. DOI: 10.1016/j.expthermflusci.2021.110526 |
| `brunton2020-ml-fluid-mechanics.pdf` | `brunton2020` | **HIGH** | Brunton et al. (2020) Ann. Rev. Fluid Mech. — review de ML para mecânica dos fluidos (9.1MB, arXiv). DOI: 10.1146/annurev-fluid-010719-060214 |
| `morshed2020-slug-wavelet-flow-regime.pdf` | `morshed2020` | **HIGH** | Morshed et al. (2020) Appl. Sci. — análise wavelet de regimes slug em escoamento bifásico ar/fluido. DOI: 10.3390/app10093272 |
| `qiao2023-cnn-attention-flow-pattern.pdf` | `qiao2023` | **HIGH** | Qiao et al. (2023) JMSE — CNN ResNet50 com mecanismo de atenção duplo (CBAM+ECA) para padrões gás-líquido. DOI: 10.3390/jmse11040793 |
| `he2016-resnet-deep-residual.pdf` | `he2016resnet` | **HIGH** | He et al. (2016) arXiv — Deep Residual Learning (ResNet), paper fundamental para CNN/ResNet50. arXiv: 1512.03385 |
| `alzubaidi2021-deep-learning-review.pdf` | `alzubaidi2021` | **HIGH** | Alzubaidi et al. (2021) J. Big Data — review completo de deep learning: conceitos, arquiteturas CNN. DOI: 10.1186/s40537-021-00444-8 |
| `cunningham2022-knn-tutorial.pdf` | `cunningham2022` | **HIGH** | Cunningham & Delany (2022) ACM Comput. Surv. — tutorial de kNN; referência teórica para o classificador k=15 usado em Dias 2025. DOI: 10.1145/3459665 |
| `lofstedt2019-haralick-texture-invariant.pdf` | `lofstedt2019` | **HIGH** | Löfstedt et al. (2019) PLOS ONE — features de textura Haralick invariantes ao nível de cinza; complementa Haralick 1973. DOI: 10.1371/journal.pone.0212110 |
| `yang2023-cv-ml-heat-transfer-microchannels.pdf` | `yang2023` | **HIGH** | Yang et al. (2023) Energies — review: visão computacional e ML para transferência de calor em microcanais. DOI: 10.3390/en16031500 |
| `suh2021-deep-learning-boiling-heat-transfer.pdf` | `suh2021` | **HIGH** | Suh et al. (2021) Sci. Reports — deep learning para predição de transferência de calor na ebulição. DOI: 10.1038/s41598-021-85150-4 |
| `xu2020-cnn-gas-oil-ect-sensors.pdf` | `xu2020` | **MEDIUM** | Xu et al. (2020) Sensors — CNN com sensores ECT para padrões gás-óleo; metodologia com sensores não-visuais. DOI: 10.3390/s20041200 |
| `he2021-ann-pump-flow-pattern.pdf` | `he2022` | **MEDIUM** | He et al. (2022) Micromachines — ANN + SMOTE para padrões em bomba centrífuga (geometria diferente). DOI: 10.3390/mi13010002 |
| `masiukiewicz2021-flow-pattern-narrow-channel.pdf` | `masiukiewicz2021` | **MEDIUM** | Masiukiewicz & Anweiler (2021) Energies — estereologia para avaliação de padrões em canal retangular estreito. DOI: 10.3390/en14113180 |
| `ruizdiaz2024-transformer-flow-pattern-vertical.pdf` | `ruizdiaz2024` | **MEDIUM** | Ruiz-Díaz et al. (2024) Inventions — transformer neural network para padrões em tubulação vertical óleo-água. Vertical: diferente do setup. DOI: 10.3390/inventions9010015 |
| `zhang2024-arxiv-flow-pattern-vertical-bp.pdf` | `zhang2024bp` | **MEDIUM** | Zhang et al. (2024) arXiv — BP neural network para padrões em tubo vertical. Metodologia simples; vertical. arXiv: 2410.13906 |
| `cheng2023-flow-patterns-microchannels.pdf` | `cheng2023` | **HIGH** | Cheng & Xia (2023) ETFS — revisão abrangente de padrões de escoamento bifásico gás-líquido em microcanais; define critério de microescala Co < 0,5 (Kew & Cornwell). **ATENÇÃO: TCC cita Co > 0,5 — desigualdade invertida.** DOI: 10.1016/j.expthermflusci.2022.110726 |
| `mallat2009-wavelet-tour-signal-processing.pdf` | `mallat2009` | **CRITICAL** | Mallat (2009) "A Wavelet Tour of Signal Processing, 3ª ed." — referência para DWT: banco de filtros h[n]/g[n] (Fig. 7.12, Seção 7.3) e complexidade O(N) (Seção 7.3.1). Texto legível via pdftotext. ISBN: 978-0-12-374370-1 |
| `tan2015-oil-water-wavelet-sensors.pdf` | `tan2015` | **HIGH** | Tan et al. (2015) CES — escoamento óleo-água com sensor condutividade/capacitância + análise wavelet. **ATENÇÃO: TCC afirma "DWT" e "energia por nível" — o artigo usa CWT + LWE (coeficientes Wn/an). A "robustez a variações de velocidade" não é afirmada no artigo.** DOI: 10.1016/j.ces.2015.04.046 |
| `sun2012-wavelet-entropy-flow-patterns.pdf` | `sun2012` | **HIGH** | Sun et al. (2012) MST — ar-água, tubo 50 mm (macrocanal), DWT db4 7 níveis, entropia norma coef. d4, 4 padrões, >95% acc. Sensor: flutuações de pressão através de bluff body (vortex). DOI: 10.1088/0957-0233/23/12/125301 |
| `nguyen2010-wavelet-flow-pattern-discrimination.pdf` | `nguyen2010` | **HIGH** | Nguyen, Euh & Song (2010) IJMF — CWT aplicada a sinais de fração de vazio (Impedance Void Meter, IVM) para discriminação objetiva de padrões bifásicos ar-água em tubo vertical de 80 mm (KAERI). Define "local wavelet energy coefficient" e "effective scale" como critérios de identificação. Padrões identificados: bubbly discreto, clustered/cap bubbly, slug, churn. **Antecedente mais direto da abordagem por CWT para discriminação de padrões proposta neste TCC.** Citado por deoliveira2021 e deoliveira2019. PDF legível via pdftotext. DOI: 10.1016/j.ijmultiphaseflow.2010.04.007 |

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
| `cheng2023` | Cheng, Xia | 2023 | Flow patterns and flow pattern maps for adiabatic and diabatic gas-liquid two-phase flow in microchannels | **HIGH** (PDF local disponível) — ver tabela acima | 10.1016/j.expthermflusci.2023.110988 |

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

## Livros de referência (sem PDF necessário)

| ID bib | Autores | Ano | Título | Relevância | Editora |
|--------|---------|-----|--------|-----------|---------|
| `mallat2009` | Mallat | 2009 | A Wavelet Tour of Signal Processing (3ª ed.) | **CRITICAL** — texto de referência para teoria de wavelet; cobre CWT, DWT, WPS | Academic Press. ISBN: 978-0-12-374370-1 |
| `bishop2006` | Bishop | 2006 | Pattern Recognition and Machine Learning | **HIGH** — referência padrão de ML; cobre kNN, validação cruzada, métricas | Springer. ISBN: 978-0-387-31073-2 |
| `gonzalez2018` | Gonzalez, Woods | 2018 | Digital Image Processing (4ª ed.) | **HIGH** — processamento de imagem digital; histograma, escalas de cinza, GLCM | Pearson. ISBN: 978-0-13-335672-4 |
| `thome2004` | Thome | 2004 | Engineering Data Book III | **HIGH** — dados e correlações de transferência de calor em escoamento bifásico | Wolverine Tube Inc. |

---

## Resumo de cobertura por pilar

| Pilar | CRITICAL (local) | CRITICAL (sem PDF) | HIGH local | HIGH sem PDF | Total |
|-------|-----------------|-------------------|------------|-------------|-------|
| Padrões de escoamento | deoliveira2021, deoliveira2019 | wojtan2005, revellin2006 | schmid2022, brunton2020, cheng2023 | saisorn2008, deoliveira2016/17, seal2021 | 10 |
| DSP / Wavelet | torrence1998, schuck2019, mallat2009 | — | morshed2020, tan2015, sun2012, nguyen2010 | wei2020, daubechies1990 | 9 |
| ML classificação | dias2025-knn, dias2025-cnn, he2016resnet | — | qiao2023, brunton2020, alzubaidi2021, suh2021 | zhang2020, seal2021 | 9 |
| GLCM / textura | lofstedt2019 | haralick1973 | — | — | 2 |
| kNN / classificadores | cunningham2022 | — | yang2023 | — | 2 |
| Setup experimental | deoliveira2021, deoliveira2019, deoliveira2013-thesis, dias2024-thesis, dias2025-knn, dias2025-cnn | — | — | deoliveira2016/17 | 6 |

**Para os artigos sem PDF:** acesse via VPN UNISINOS (Portal Capes) ou solicite ao orientador — todos têm DOI verificado.

---

## Erros corrigidos nesta sessão

- `dias2025knn`: autores corrigidos (Isadora → Lucas Souza Dias; Copetti adicionada)
- `dias2025cnn`: autores corrigidos (Cardoso adicionada; Oliveira corrigido)
- `deoliveira2021`: título e autores corrigidos via CrossRef (Indrusiak, de Souza, Lima adicionados)
- `deoliveira2019`: título corrigido — era "Experimental analysis of flow patterns and void fraction for flow boiling of isobutane" (esse paper não existe com esse título); correto é "Nonlinear analysis of characteristic frequencies and chaotic behavior..."
- `tan2015`: título e DOI corrigidos — o paper citado por dias2025knn e deoliveira2021 é "Characterization of oil-water two-phase pipe flow with combined conductivity/capacitance sensor and wavelet analysis"; o DOI 10.1016/j.ces.2015.04.046 está correto
- `hobold2018`: nome do primeiro autor corrigido (Gabriel M. → Gustavo M.) e título corrigido
- `cheng2023`: título e autores completamente corrigidos via CrossRef (era paper de imagem e ML, mas o DOI aponta para review de mapas de padrões)
