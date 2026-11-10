# Follow-up — A3 inaplicável: banda de costura abaixo do piso de análise (50 px)

## Problema (encontrado ao rodar `run_stats.py`)

A robustez A3 deveria testar se o artefato de costura (bump quase-periódico em λ≈shift_mean
injetado por `_build_canvas` em `pipeline_v11.py`) contamina a assinatura A1. Mas:

- Banda de costura λ∈[0.7,1.4]×shift_mean = **[18,5; 37,0] px** (média dos 19 runs).
- Eixo de análise da CWT: **λ ∈ [50, 3000] px** (`FREQ_MAX = 1/50` em `wavelet_v8.py:19`).
- Interseção = **0 escalas**. Só 4 runs têm a borda 1,4×shift acima de 50 px (máx. 64,5 px).

Logo `banda-completa ≡ banda-mascarada` (R idêntico nas 4 condições). **A3 não confirma
ausência de artefato — apenas mostra que o artefato está abaixo da janela observada.** O
descritor A1 é imune ao bump por construção (piso 50 px), o que é bom para a validade de A1,
mas deixa a *pergunta do artefato* sem resposta direta.

Nota: já corrigido um bug secundário de unidade na máscara (comparava `scales` cru em vez de
comprimento de onda `scales/center_freq`); não muda o veredito (0 escalas).

## Opções para de fato testar o artefato de costura

### Opção A (recomendada, mínima) — declarar out-of-band e encerrar A3
O bump de costura (~14–46 px) é fisicamente inatingível pelo descritor (piso 50 px). Em vez de
mascarar, **provar** que está fora: computar o WPS auto de UM run SEM piso (estender
`FREQ_MAX` a 1/10) e mostrar o pico em λ≈shift_mean vivendo em 15–40 px, longe da banda 50–3000
usada em A1. Uma figura + um parágrafo. Substitui A3 por "A3′: artefato localizado fora da
janela de análise". Custo: ~20 linhas, 1 run.

### Opção B — baixar o piso e reavaliar A1 com o artefato DENTRO
Reduzir `FREQ_MAX` (ex.: 1/10) no `extract_descriptors.py` (NÃO no motor imutável — a extração
já é cópia), reextrair descritores e recomputar A1. Se R **sobe** ao incluir 15–50 px → parte
da separação vinha do bump de costura (colinear com velocidade→padrão) → A1 atual (piso 50)
é a versão CONSERVADORA e defensável. Se R **não muda** → bump não separa padrões. Custo:
reextração (~1 min FFT) + rerun. Dá o teste que A3 queria, ao contrário de mascarar.

### Opção C — abandonar A3 formalmente
Se o dono aceitar que o piso de 50 px de Schuck já isola o descritor do bump, remover A3 do
plano e do documento, mantendo só A1+A2. Registrar a justificativa (out-of-band) em
`CONCLUSOES.md` (já feito).

## Recomendação
**A + B combinadas:** A (figura out-of-band, 1 run) documenta ONDE o artefato mora; B (reextração
com piso 1/10) mede SE ele separa padrões. Juntas convertem o no-op de A3 numa evidência real
de que A1 (piso 50 px) é a leitura conservadora. Atualizar depois `PLANO_ESTATISTICA_V8.md` e
`CONCLUSOES.md` com o veredito.

## Arquivos afetados
- `src2/estatistica/extract_descriptors.py` — parametrizar `FREQ_MAX` (default 1/50; variante 1/10).
- `src2/estatistica/run_stats.py` — substituir A3 (máscara) por A3′ (out-of-band + comparação de
  R entre pisos), OU remover A3.
- `PLANO_ESTATISTICA_V8.md` — corrigir A3 (a máscara pressupunha banda dentro da janela; não está).
- `CONCLUSOES.md` — acrescentar veredito do teste de piso.
