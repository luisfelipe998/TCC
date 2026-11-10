# v16.2 — Plano: alternative WCO entre eventos equivalentes

**Estado: planejado, não executado.** Este arquivo fornece contexto para um próximo agente. A autorização atual cobre somente este plano; aguardar pedido de execução antes de implementar ou rodar a análise. Não alterar v5/v16 nem executar automaticamente seus lotes de reprodução.

## Pergunta e hipótese

Investigar se estruturas visualmente comparáveis, selecionadas **antes de observar coerência**, apresentam alternative WCO maior e mais específica que comparações sem correspondência entre eventos. Aplicar as mesmas alternative WPS, alternative WCS e alternative WCO aos segmentos dessas estruturas.

Isso aproxima a pergunta de uma comparação entre estruturas. **Ainda não está demonstrado que eventos equivalentes sejam delimitáveis de modo confiável em todos os padrões.** Um resultado positivo responderá primeiro sobre eventos, não automaticamente sobre panoramas completos, classificação geral ou identidade física entre bolhas de aquisições independentes.

## Contexto herdado: leitura obrigatória

A implementação de referência é a corrigida de v16_spatial2, mas correção matemática não significa melhor desempenho. Ler:

1. [Conclusão executiva](../v16_spatial2/EXECUTIVE_SUMMARY.md), com figuras dos cinco padrões e exemplos WCS/WCO.
2. [Auditoria das equações](../v16_spatial2/EQUATION_AUDIT.md) e [PDF original de Schuck](../../assets/schuck-bardo-wavelet-coherence.pdf).
3. [Relatório](../v16_spatial2/REPORT.md), [resultados](../v16_spatial2/RESULTS.md), [pré-processamentos](../v16_spatial2/PREPROCESSING.md) e [README de reprodução](../v16_spatial2/README.md).
4. [Inventário das entradas](../v16_spatial2/results/input_audit.csv), [conclusão numérica](../v16_spatial2/results/preprocessing_conclusion.json) e [núcleo matemático](../v16_spatial2/scripts/wavelet_core.py).

**Resultado anterior, não desta proposta:** o v16 não demonstrou melhora robusta e geral da alternative WCO frente ao v5. Dez pré-processamentos foram triados em 325 pares com 64 escalas; baseline, média central de 20 linhas e mediana foram examinados com 200 escalas. A média das linhas 30–49 foi o candidato exploratório mais promissor, frente à original 20–59.

Na comparação com v5, o escore dos 87 pares intrapadrão passou de 0,163237 para 0,164622: aproximadamente +0,001385, **ganho marginal**. Com peso igual por padrão, passou de 0,182027 para 0,197219. Melhoraram 55 pares e pioraram 32; plug ganhou mais e churn perdeu na média. Os controles entre padrões ainda superaram os intrapadrão no escore principal.

Esse escore usa pesos de alternative WPS do baseline e máximos por escala sobre lags. Não é média simples de WCO nem outra equação de Schuck. Recuperar sua definição exata nos scripts antes de reutilizá-lo. O v16 corrigiu discretização da Morlet, conversão de escala, sinal do lag, recorte antes de preparação e representação de bordas. Não interpretar essas correções como prova de similaridade maior.

## Dados, proveniência e escopo confirmado

Usar os rótulos e a inclusão do CSV de auditoria, não rótulos antigos de INSIGHTS do v5. O escopo autorizado é de **31 entradas auditadas e 26 aquisições nas comparações**:

| Padrão | Aquisições elegíveis | Pares intrapadrão |
| --- | --- | --- |
| Dispersed bubbles | m07q01_001 | 0 |
| Plug | m05q01_001, m06q01_001, m07q05_001 | 3 |
| Slug | m04q01_001, m05q05_001, m06q05_001 | 3 |
| Churn | m03q05_001, m04q10_001, m04q20_001, m05q10_001, m05q20_001, m05q40_001, m06q10_001, m06q20_001, m07q20_001, m07q40_001 | 45 |
| Wavy annular | m03q10_001, m03q20_001, m03q40_001, m04q40_001, m04q80_001, m05q80_001, m06q40_001, m06q80_001, m07q80_001 | 36 |

São 87 pares intrapadrão, 238 entre padrões, 325 pares únicos. Pares que compartilham aquisição não são réplicas independentes. Permanecem excluídos das comparações: m03q01_001, m03q80_001, m04q05_001, m05q05_002 e m07q10_001. Não reincluir para melhorar resultados. Recortes de bubbles não criam novas aquisições independentes.

As entradas ficam em `src2/data_processed/<run>/`: __panorama_v11.png, __dados_cropped_v11.xlsx e quadros ImgA*.jpg. A coluna Voltage_2 é intensidade média das linhas **20–59 inclusive**, não tensão de sensor. O [pipeline_v11.py](../pipeline_v11.py) estima deslocamentos e costura faixas centrais. [pipeline_v13.py](../pipeline_v13.py) extrai a linha 40; [pipeline_v9.py](../pipeline_v9.py) trata sinal temporal e não gerou esses panoramas. Os scripts estão diretamente em src2, apesar da referência inicial a “pipelines”.

A auditoria anterior reproduziu os 31 panoramas pixel a pixel e verificou sinais/posições das planilhas. Isso confirma reprodução, não validade física da reconstrução. O pipeline substitui deslocamentos não positivos por 1 pixel: entre elegíveis, m03q40_001 teve 60/500 e m04q01_001 teve 1/500. Costuras, mistura de bordas e evolução das estruturas podem afetar o espectro. Preservar esses diagnósticos, sem descartar resultados inconvenientes.

A amostragem é **1 pixel/amostra**. Não inventar mm/pixel, Hz, velocidade ou identidade física entre estruturas de registros distintos. [De Oliveira 2021](../../assets/deoliveira2021-flow-patterns-pressure-drop.pdf) aplica CWT a sensores temporais; seus espectros não são a mesma estatística das alternative WPS de panoramas.

## Contrato matemático obrigatório

Autoridade: Schuck Jr. e Bodmann (2019), capítulo 32, DOI 10.1007/978-3-030-16077-7_32, PDF em assets. Conferir diretamente as páginas impressas 411, 413 e 414, também preservadas em [411](../v16_spatial2/literature/schuck-page-03.png), [413](../v16_spatial2/literature/schuck-page-05.png) e [414](../v16_spatial2/literature/schuck-page-06.png). Não resolver divergências escolhendo a versão que aumenta coerência.

CWT com normalização 1/sqrt(a), Eq. 32.4; Morlet e conversão de escala:

$$\psi(t)=\pi^{-1/4}e^{-t^2/2}e^{i\omega_0t},\qquad\omega_0=6\quad(32.6),$$

$$\lambda=\frac{4\pi a}{\omega_0+\sqrt{\omega_0^2+2}}\quad(32.8).$$

**Alternative WPS**, Eqs. 32.14–32.15:

$$P_x^{alt}(b,a)=E\{|\mathrm{CWT}\{R_{xx}\}(b,a)|\},\quad P_y^{alt}(b,a)=E\{|\mathrm{CWT}\{R_{yy}\}(b,a)|\}.$$

**Alternative WCS**, Eq. 32.16:

$$C_{xy}^{alt}(b,a)=E\{\mathrm{CWT}\{R_{xy}\}(b,a)\}.$$

**Alternative WCO**, Eq. 32.17:

$$\Gamma_{xy}^{alt}(b,a)=\frac{|C_{xy}^{alt}(b,a)|^2}{|P_x^{alt}(0,a)|\,|P_y^{alt}(0,a)|}.$$

Alternative WPS é **módulo, não módulo ao quadrado**. Alternative WCS é complexa; seu módulo pode ser mostrado com legenda explícita. A WCO tem numerador ao quadrado e denominador **no lag zero**, transmitido aos demais lags de cada escala. Não substituir por denominador local, máximo, potência convencional, suavização de Grinsted ou ajuste que eleve o resultado. Não aplicar clipping para esconder valores fora de [0,1].

Convenção cruzada: soma de x[n+m] conjugado(y[n]); y atrasado 100 amostras em relação a x produz pico em **−100**. Começar com **K=1**. Para K>1, a Eq. 32.16 (média complexa) diverge da Eq. 32.20 impressa (média dos módulos). CWT da correlação média tampouco substitui em geral a média dos módulos de 32.18–19. Não misturar eventos ou registros em um ensemble silenciosamente; média descritiva de resultados K=1 deve ser identificada como tal.

**Adaptação explícita já autorizada:** Schuck usa divisor N−|m|, não enviesado, nas Eqs. 32.21–32.23. O usuário confirmou divisor N, enviesado. Manter 1/N como ponto de partida, declarando que não reproduz literalmente esses estimadores do artigo. O não enviesado já foi testado no v16, com instabilidade em grandes lags e WCO >1. A estrutura espectral não negativa da correlação enviesada é relevante à estabilidade; não garante confiabilidade experimental de todo ponto.

O núcleo atual fornece prepare, correlation, morlet_cwt, transform_correlation, pair, safe_coherence, masks e masked_mean. prepare inverte intensidade, remove tendência linear/média e normaliza pelo desvio padrão. **Recortar antes dessa preparação.** correlation exige comprimentos iguais. O motor principal é morlet; legacy reproduz a orientação de lag invertida do v5. Verificar os argumentos e registrar hash do código reutilizado.

## Janelas, suporte, evidências reutilizáveis e ambiente

Janelas adicionais não são exigidas por Schuck. No v16, a COI geométrica exige N−1−|b| >= sqrt(2)*a; o controle conservador opcional exige |b|+3*a <= N/2. São regras de suporte/relato, **não novas equações ou testes de significância**. Mostrar bordas e denominadores inválidos. Self-WCO=1 no lag zero é identidade algébrica em pontos válidos, não evidência experimental de semelhança.

Consultar [validação numérica](../v16_spatial2/results/numerical_validation.json), [baseline literal](../v16_spatial2/results/legacy_literal_validation.json), [validações adicionais](../v16_spatial2/results/additional_validation.json), [verificação final](../v16_spatial2/results/final_verification.json) e [manifesto](../v16_spatial2/results/artifact_manifest.json). As verificações anteriores não substituem testes das novas operações propostas.

No v16, `intermediate/<run>_input.npz` contém sinais, deslocamentos e costuras verificados. Arquivos *_maps.npz têm lags reduzidos para renderização (~1800 posições): **não recalcular máximos integrais a partir deles**. Usar tabelas de resultados ou recalcular apenas o necessário.

Ambiente anterior: Python 3.13.2; dependências em [requirements.txt](../v16_spatial2/requirements.txt). PyWavelets apresentou distribuição 1.9.0 e versão do módulo 1.8.0; registrar se repetir. A Morlet principal usa NumPy/SciPy diretamente. Usar MPLCONFIGDIR gravável para figuras. src2 é ignorado pelo Git atual; não alterar essa regra nem presumir versionamento dos novos arquivos.

## Etapas propostas para a execução futura

### 1. Verificar observabilidade antes de calcular coerência

Inspecionar panoramas e quadros originais dos 26 elegíveis, aproveitando a auditoria dos 31. Verificar estruturas delimitáveis, fronteiras, costuras e se o sinal das linhas 20–59 representa o evento visível. Não começar procurando os maiores mapas de WCO.

Plug/slug podem sugerir interfaces de estruturas alongadas; churn pode ter regiões fragmentadas sem limites repetíveis; annular pode sugerir ondulações; bubbles pode conter pequenas estruturas. **São possibilidades a verificar, não anotações locais já comprovadas.** O rótulo global não certifica a morfologia de cada evento. Não impor uma taxonomia única só para completar tabelas.

Se a interpretação física for ambígua, preparar exemplos e consultar o usuário antes de chamar estruturas de equivalentes. Não inventar fronteiras ou identidade. Se não houver eventos observáveis ou suporte adequado, registrar inviabilidade antes do lote completo.

### 2. Congelar regras de anotação sem ver mapas

Criar futuramente PROTOCOL.md, DECISIONS.md e catálogo de eventos. Definir morfologia, âncora visual, qualidade mínima, fronteiras, contexto, sobreposição, desempate e exclusões por artefatos. Preservar exemplos favoráveis, desfavoráveis, rejeitados, ambíguos e sem parceiro; medir cobertura por aquisição/padrão.

Separar **detecção/anotação**, **pareamento** e **avaliação wavelet**. As duas primeiras não poderão usar WCO/WCS nem busca espectral para maximizar o resultado. Preferir regras mensuráveis na imagem. Na anotação manual, registrar incerteza e, quando viável, repetibilidade em uma segunda leitura sem acesso às métricas.

Catálogo mínimo: event_id, aquisição, rótulo global, categoria morfológica local, coordenadas inicial/final no panorama, âncora, linhas utilizadas, contexto, proximidade de costuras, evidência nos quadros quando rastreável, qualidade, inclusão/exclusão, motivo, versão da regra e modo de anotação. Não atribuir quadro ou identidade que a reconstrução não permita recuperar.

Pilotar viabilidade em subconjunto definido por IDs e cobertura de padrões, sem WCO como critério. Congelar regras antes de avaliar outras aquisições. Se não houver dados para separar desenvolvimento e avaliação, declarar exploração; o conjunto já foi examinado no v16.

### 3. Definir equivalência, pareamento e controles

Comparações principais: eventos de **aquisições diferentes**. Equivalência deve usar critérios visuais/tolerâncias predefinidos, por exemplo categoria local e extensão em pixels, se observáveis. Evitar critérios tão estreitos que preservem apenas estruturas quase idênticas. Registrar sem-parceiros e cobertura perdida.

Usar todos os pares admissíveis ou amostragem determinística com semente e limite por aquisição. Resolver múltiplos candidatos por regra fixa/sorteio, nunca pela maior WCO. Não escolher retrospectivamente o melhor evento de cada registro. Se houver otimização do pareamento, seu objetivo deverá usar apenas atributos visuais previstos.

Controles obrigatórios, com comprimento/contexto/qualidade comparáveis quando possível:

- Eventos considerados equivalentes entre aquisições do mesmo padrão.
- Eventos não equivalentes ou pares aleatórios dentro do mesmo padrão.
- Eventos entre padrões distintos.
- Estruturas localmente semelhantes entre padrões distintos, se existirem: semelhança local não implica mesmo rótulo global.
- Segmentos de referência selecionados por regra fixa sem anotação de evento.
- Self-comparação e atrasos conhecidos como controles numéricos; eventos da mesma aquisição reportados separadamente.

Se algum controle não existir, registrar ausência. Bubbles só tem m07q01_001: seus eventos não demonstram repetibilidade entre aquisições desse padrão.

### 4. Resolver comprimentos e contexto sem fabricar semelhança

correlation exige comprimentos iguais. Definir antes do cálculo um protocolo que preserve o evento e **1 px/amostra**. Opções a verificar: segmentos de comprimento comum ao redor de âncoras visuais, com contexto bilateral; ou estratos de comprimento. Escolher por suporte e preservação da estrutura, não por WCO. Registrar perda de fronteira/contexto.

Não usar dynamic time warping, alongamento/compressão do eixo ou alinhamento escolhido pela coerência. A âncora visual pode definir origem; o lag mantém a convenção de Schuck. Zero padding computacional não cria observações, não aumenta N estatístico e não legitima escalas sem suporte.

Recortar antes de inverter/detrender/normalizar; rejeitar sinais degenerados com motivo. Linhas 20–59 são a referência; 30–49 somente sensibilidade predefinida. Se a projeção média não representar o evento, documentar antes de escolher outra faixa. Não escolher o perfil de maior coerência por evento.

O segmento é o objeto desta proposta, **não janela exigida por Schuck**. Não acrescentar Hann ou outra ponderação por padrão; se necessária por evidência de borda, tratar como experimento distinto e controlado.

Eventos curtos podem não sustentar a faixa anterior de 50–3000 px. Derivar escalas válidas de amostragem, extensão e suporte antes de avaliar rótulos; usar faixa comum nas comparações e relatar cobertura. Não reduzir artificialmente a banda para aumentar WCO. Separar COI geométrica de máscara conservadora opcional; ausência de suporte não significa coerência zero.

### 5. Validar o recorte e definir resumos

Controles pequenos antes do lote: evento idêntico, atraso conhecido, evento com ruído, eventos distintos, fundo comum sem evento equivalente e cortes próximos às bordas. Verificar lag, índice zero, denominador, padding e efeito de detrending em eventos curtos. Recalcular a Eq. 32.17 a partir dos coeficientes, sem clipping.

Eventos curtos e suaves podem ter coerência alta com muitos eventos distintos. Investigar nos controles negativos. Máximos buscados em muitos lags/escalas favorecem seleção: definir previamente domínio de lag pela âncora e sua incerteza visual. Não buscar o melhor deslocamento e depois tratá-lo como predefinido.

Começar com **K=1 por par de segmentos**. Agregação posterior de mapas/escores é descritiva; não a chamar automaticamente de WCO de ensemble. Resolver explicitamente a divergência 32.16/32.20 antes de qualquer extensão K>1.

Definir futuramente um resumo principal de WCO na região válida com pesos explícitos em escala/lag; máximo global fica secundário. Se reutilizar o escore do v16, preservar definição e pesos de referência, declarando mudanças de domínio. Mostrar distribuições, suporte e denominadores inválidos. Não aplicar limiar convencional de Grinsted. Eventuais surrogates precisam representar o nulo e o procedimento de seleção/resumo relevante, sem criar significância por ignorar a busca.

### 6. Separar implementação, recorte e pareamento

| Referência | Pergunta |
| --- | --- |
| v5 literal e v16 corrigido nos mesmos segmentos de eventos | Quanto da diferença vem do núcleo numérico? |
| Segmentos de comprimento/contexto equivalentes, sem pareamento, com núcleo v16 | Quanto vem apenas de encurtar sinais? |
| Eventos pareados por critérios visuais, com núcleo v16 | O pareamento acrescenta especificidade? |
| Panoramas/prefixos comuns anteriores de v5/v16 | Quanto mudou o objeto em relação ao resultado global? |

Não comparar médias em escalas/domínios distintos como prova direta de superioridade. Preferir eventos e controles no mesmo domínio e separar resultados que exijam outra faixa. Uma WCO maior em eventos não substitui a conclusão sobre panoramas completos.

Agregar por par de aquisições e por padrão, com versão balanceada; centenas de eventos de uma aquisição não podem dominar a evidência. Para ajuste/validação, manter todos os eventos de cada aquisição no mesmo grupo. Incerteza deve respeitar dependência de eventos/pares e o número pequeno de aquisições; não apresentar eventos como réplicas experimentais independentes.

## Critérios de relevância e adaptação

Há potencial de relevância se o ganho aparecer em eventos selecionados sem WCO, superar controles de recorte/pareamento, não elevar igualmente coerência entre padrões/mismatches, preservar cobertura útil e resistir a variações justificadas das fronteiras. Definir antes do lote os critérios de efeito/cobertura, com justificativa. Consultar o usuário sobre requisitos práticos se necessários à decisão de substituição.

| Descoberta | Caminho justificável |
| --- | --- |
| Fronteiras não são repetíveis | Relatar inviabilidade ou utilidade restrita aos padrões viáveis. |
| Segmentos não sustentam escalas | Rever contexto com regra independente dos mapas ou marcar não estimável. |
| WCO aumenta igualmente em controles | Interpretar efeito genérico do recorte/fundo, não ganho específico. |
| Só poucos casos selecionados funcionam | Relatar baixa cobertura; não generalizar. |
| Ganho exige ajustar fronteiras após ver mapas | Marcar exploração e avaliar novamente com fronteiras congeladas. |
| Ganho estável, restrito a certos eventos | Relatar explicitamente utilidade restrita. |
| Faltam interpretação física, calibração ou mudança de escopo | Consultar o usuário antes de redefinir. |

Registrar descoberta, motivo da mudança, regra anterior/nova, dados já observados e o que ainda poderá ser avaliado sem ajuste ao resultado. Não migrar para classificação por imagem ou métodos fora de wavelets sem autorização. Se o ganho for pequeno, dizer **“melhora marginal”**; se não houver melhora convincente, dizer **“não conseguimos melhorar de forma convincente a alternative WCO com este protocolo”**.

## Entregáveis futuros, tempo e encerramento

Após pedido de execução: catálogo auditável com coordenadas/imagens/exclusões, protocolo e decisões, scripts/configuração, tabelas de eventos/pares/agregados, testes/controles, manifesto e instruções reproduzíveis. Resumo executivo com panorama anotado, sinais, alternative WPS, módulo de alternative WCS e alternative WCO, incluindo bons e maus casos. Mostrar bubbles e churn e suas limitações, mesmo sem comparação intrapadrão viável.

Antes do lote pesado, medir poucos pares de comprimentos distintos, estimar eventos/pares finais, escalas, controles, tempo e memória. **Informar a estimativa ao usuário antes de iniciar**, atualizando se mudar cobertura ou número de controles. Separar tempo de anotação e computação. Os ~5,8 minutos do lote anterior não estimam esta proposta. Não montar produtos cartesianos grandes antes de dimensionar custo.

Encerrar respondendo se equivalência visual acrescentou similaridade wavelet específica, para quais eventos e com que cobertura. Resultado negativo/inconclusivo é aceitável; não prometer WCO maior por obrigação. A [v16.1](../v16.1/PLAN.md) é independente: seus perfis podem ajudar a interpretar divergências, mas não selecionar retrospectivamente eventos com a melhor coerência.
