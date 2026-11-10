
## Identidade dos canais — CONFIRMADA pelo dono dos dados

Até o v3.2 a identidade das três colunas de `dados.xlsx` era **inferida** pelo
comportamento físico (ver `src/v3.1_classic/ANALISE_V3.1_CLASSIC.md`). O dono dos
dados **confirmou** o significado, encerrando a inferência:

| Coluna (xlsx) | Sinal | Grandeza | Unidade |
|---------------|-------|----------|---------|
| **Voltage_0** | **PD** — queda de pressão na seção de teste | pressão diferencial | kPa |
| **Voltage_1** | **G** — vazão mássica do fluido | vazão mássica | g/s |
| **Voltage_2** | **OS** — sensor óptico na seção de visualização | luminosidade da interface | V |

**Implicação:** ao contrário do que o v3.1/v3.2 assumiam (Voltage_0 = aquecedor/q″,
sem PD), **há sinal de PD** nestes dados (Voltage_0). Os dois sinais que de Oliveira
et al. (2021) analisa — **PD e OS** — estão ambos disponíveis. A análise CWT dos dois
canais é feita em `src/v3.3_classic/` (a limitação de "só OS" registrada no v3.1/v3.2
foi superada por esta confirmação).

Notas adicionais do dono dos dados:
- **Não há mais os dados de CWT da época** do paper; os dados disponíveis são as
  colunas dos arquivos Excel (é o que está em `data_processed`).
- As **classificações manuais** foram feitas na época **observando imagens do banco
  de imagens** da seção de visualização.

---

## Rótulos finais — revisão contra o banco de imagens

Revisão das classificações manuais contra o banco de imagens (repassada pelo dono
dos dados). Esta tabela é a **fonte de verdade única** de rótulos: nomeia as pastas
das análises (lida pelos scripts `wavelet_cwt_*.py`) e é a referência para avaliar
qualquer descritor espectral (centroide, coerência *Wavelet*, etc.) em qualquer run.
Um rótulo por condição; nenhuma outra rotulagem é usada como referência. Correções em
relação ao rótulo manual anterior estão marcadas em **(corrigido)**.

Principais mudanças: vários `anular` foram especificados como **wavy annular** ou
**smooth annular**; `m03q10` passou de churn para **wavy annular**; `m07q01` de
bolha dispersa para **dispersed bubbles** (mesmo regime, nomenclatura do banco);
`m07q10` de slug para **slug-churn**.

| Condição | Rótulo final | Mudou? |
|----------|--------------|--------|
| m03q01 | plug-slug | — |
| m03q05 | churn | — |
| m03q10 | wavy annular | **(corrigido: era churn)** |
| m03q20 | wavy annular | **(corrigido: era anular)** |
| m03q40 | wavy annular | **(corrigido: era anular)** |
| m03q80 | smooth annular | **(corrigido: era anular)** |
| m04q01 | slug | — |
| m04q05 | slug-churn | — |
| m04q10 | churn | — |
| m04q20 | churn | — |
| m04q40 | wavy annular | **(corrigido: era anular)** |
| m04q80 | wavy annular | **(corrigido: era anular)** |
| m05q01 | plug | — |
| m05q05 | slug | — |
| m05q10 | churn | — |
| m05q20 | churn | — |
| m05q40 | churn | — |
| m05q80 | wavy annular | **(corrigido: era anular)** |
| m06q01 | plug | — |
| m06q05 | slug | — |
| m06q10 | churn | — |
| m06q20 | churn | — |
| m06q40 | wavy annular | **(corrigido: era anular)** |
| m06q80 | wavy annular | **(corrigido: era anular)** |
| m07q01 | dispersed bubbles | **(corrigido: era bolha dispersa)** |
| m07q05 | plug | — |
| m07q10 | slug-churn | **(corrigido: era slug)** |
| m07q20 | churn | — |
| m07q40 | churn | — |
| m07q80 | wavy annular | **(corrigido: era anular)** |

---

## Assinaturas espectrais esperadas por regime — de Oliveira et al. (2021)

Referência das frequências características de cada padrão nos dois sinais, extraída do
paper-base do grupo (a análise CWT single-channel dele é replicada em
`src/v3.1_classic/`). É a base para avaliar se um descritor de frequência é coerente
com o regime: **um descritor do OS é coerente quando cai dentro da banda esperada do
regime** (coluna "OS — banda").

> de Oliveira, J. D. et al. (2021). *On the nature of flow patterns and pressure drop
> fluctuations during flow boiling.* Int. J. Multiphase Flow, 144, 103793.

**Montagem:** ebulição convectiva de R-600a (isobutano) em minicanal horizontal de
1,0 mm, aquisição a 1 kHz. **Análises aplicadas:** FFT (PSD por Welch) + **Transformada
Wavelet Contínua (CWT)** single-channel (algoritmo Torrence & Compo 1998, wavelet
Morlet), mais métodos não-lineares (retrato de fase / atrator e maior expoente de
Lyapunov). A FFT sozinha foi considerada insuficiente (as frequências mudam ao longo
do tempo); a CWT foi usada para rastrear a frequência predominante no tempo
(escalograma).

**Os dois sinais (analisados separadamente no paper, sem coerência):**
- **PD (queda de pressão)** — transdutor diferencial na seção de teste. Mede a
  flutuação global de pressão bifásica (comportamento hidrodinâmico/de momento).
- **OS (sensor óptico)** — LED + LDR (refração/brilho). Mede a passagem local da
  interface líquido-vapor num ponto (bolhas, pistões, ondas).

Regra geral do paper: **a flutuação de PD tem frequências características
sistematicamente mais BAIXAS que a de OS.** "—" = sem frequência característica clara.

| Padrão | OS — banda esperada (CWT) | OS — leitura no escalograma (tempo×freq) | OS — perfil de potência média | PD (referência do paper) |
|--------|---------------------------|------------------------------------------|-------------------------------|--------------------------|
| **bolha dispersa** | picos 4–35 Hz; energia persiste, de baixa amplitude | manchas fracas e dispersas por toda a faixa 4–35 Hz, sem banda contínua; energia baixa e granulada no tempo | platô/corcova rasa começando em ~4 Hz e decaindo até ~35 Hz; **sem pico único forte** — perfil largo e chato | picos 6–100 Hz; persiste em baixa freq |
| **plug** | banda média 2–16 Hz | faixa horizontal contínua e estável em 2–16 Hz ao longo do tempo (passagem regular de bolhas alongadas) | corcova média centrada em ~2–16 Hz; um máximo suave, sem picos discretos | 0,5–1 Hz sustentado |
| **slug** | pico dominante **~2,8 Hz** (passagem do pistão) + picos 8/12/23/32 Hz (bolhas entre pistões) | banda forte e recorrente em ~2,8 Hz (pistões) com riscos/pontos em 8/12/23/32 Hz entre eles; **padrão multi-banda** bem marcado | **pico principal ~2,8 Hz** com corcovas secundárias em 8/12/23/32 Hz — perfil **multi-modal** (assinatura mais distinta do OS) | <1 Hz; rajadas ~6 Hz |
| **churn** | alta energia; sobe a **15–28 Hz** com o fluxo mássico (Kelvin-Helmholtz) | energia intensa e espalhada, subindo para 15–28 Hz; sem banda horizontal limpa — aspecto caótico/intermitente | pico **largo deslocado para cima (15–28 Hz)**, espalhado; máximo menos definido que no slug | ≤ 8 Hz, sem freq característica |
| **anular (ondulado/liso)** | 2–16 Hz sustentado; frentes de onda **~20 Hz** | banda sustentada 2–16 Hz mais organizada que o churn, com estrias periódicas ~20 Hz (frentes de onda) | banda sustentada 2–16 Hz com realce em **~20 Hz**; menos concentrado que slug | — (flutuação some da CWT) |

Bandas OS compactas usadas na avaliação (derivadas da coluna acima, agregando os
rótulos compostos e o par ondulado/liso): plug e plug-slug **2–16 Hz**; slug e
slug-churn **2,8–32 Hz** (multi-banda); churn **15–28 Hz**; wavy annular e smooth
annular **2–20 Hz**; dispersed bubbles **4–35 Hz**.

**O que cada sinal melhor discrimina:**
- **PD** — melhor nos regimes **intermitentes** com assinatura de baixa freq nítida:
  **plug (0,5–1 Hz)** e **slug (<1 Hz)**. Fica **mudo** (sem freq característica) em
  **churn** e **anular**, onde a flutuação de pressão se suaviza.
- **OS** — melhor na **morfologia de interface** em toda a sequência: resolve
  bolha/plug/slug (2–35 Hz), instabilidades de churn (15–28 Hz) e frentes de onda do
  anular (~20 Hz). Continua **informativo no anular**, onde o PD falha.

**Implicação para o LETEF:** a banda de interface do OS (2–35 Hz) discrimina bolha,
slug e a transição churn/anular; a assinatura de baixa freq (<1 Hz) do PD isola
plug/slug.

Os resultados de cada run de análise (centroides medidos, contagem de acertos,
figuras) ficam documentados no MD da própria pasta da run (ex.:
`src/v3.4_classic/VEREDITO_V3.4_CLASSIC.md`), não aqui.

---

## Mapeamento das figuras CWT-OS de de Oliveira et al. (2021) → nossos ensaios

Cada exemplo de escalograma do sensor óptico (OS) do paper-base foi casado com o
ensaio nosso que o originou. Método: leitura das 8 figuras CWT-OS (Fig. 5, 9, 14,
18, 21) + bandas do texto, contra os V2 reprocessados com detrend + high-pass 2 Hz
(elimina o artefato de piso ~1 Hz) e descritor de banda robusto (pico da mediana
temporal, massa/persistência por banda). G→m documentado (m03=380 … m07=890).

| Fig OS (paper) | regime | G/q″ paper | → ensaio nosso | evidência | confiança |
|----------------|--------|-----------|----------------|-----------|-----------|
| 5(b) | bolha dispersa | 890 / 1 | m07q01 | texto p.14 explícito (G890, q1); único bolha nosso | alta |
| 9(b) | plug | 510 / — | m04q05 | banda 4–16 Hz (paper) = perfil m04q05 (4–8+8–16 dominam) | média |
| 9(c) | plug | 890 / — | m07q05 | banda 2–16 Hz (paper) = perfil m07q05 (pico 2,6 Hz, massa 2–16) | média |
| 14(b) | slug | 510 / — | m04q01 | pico mediano 2,81 Hz idêntico + comb 8/12/23/32 Hz visível | alta |
| 14(c) | slug | 760 / — | m06q05 (palpite) | não resolve: paper pede energia baixa 0,7–2,8 Hz; nossos m06 slug têm pico alto (13–20 Hz) | baixa |
| 18(a) | churn | 510 / — | m04q10 | maior energia 16–64 Hz dos m04, banda larga espalhada | média |
| 18(b) | churn | 760 / — | m06q10 | pico 19,7 Hz, banda 16–32 domina 46% do tempo = churn 15–28 Hz | alta |
| 21(d) | wavy-annular | 510 / 80 | m04q80 | legenda explícita (G510, q80) | alta |

Notas: Fig. 9 = plug (não slug); Fig. 10, 13 e 17 são de PD (não OS) e ficam de
fora. A legenda da Fig. 14 traz "for PD" por engano — título e texto (p.16) confirmam
OS. As figuras 5(b), 21(d), 14(b) e 18(b) são casamentos firmes; 9(b), 9(c) e 18(a)
casam a banda mas nossa classificação funde plug em slug de baixo q; 14(c) não tem
correspondente (nossos slug G760 têm energia mais alta que o slug G760 do paper).
