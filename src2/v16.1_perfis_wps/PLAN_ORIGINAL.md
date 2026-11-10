# v16.1 — Plano: comparar perfis de alternative WPS

**Estado: planejado, não executado.** Este arquivo fornece contexto para um próximo agente. A autorização atual cobre somente este plano; aguardar pedido de execução antes de implementar ou rodar a análise. Não alterar v5/v16 nem executar automaticamente seus lotes de reprodução.

## Pergunta e hipótese

Investigar se aquisições do mesmo padrão compartilham **distribuições de escalas**, sem exigir correspondência espacial entre seus eventos. Priorizar o perfil de alternative WPS no lag zero, já presente na normalização de Schuck; estudar a média em lag secundariamente. Manter alternative WCS e alternative WCO como análises complementares, sem mudar as fórmulas.

Essa representação pode ser mais adequada para comparar padrões do que a correspondência entre eventos independentes, mas **não há garantia de separação melhor que v5**. Perfis semelhantes não implicam WCO alta nem igualdade de tamanho físico das bolhas. Distâncias novas não poderão ser chamadas alternative WCO.

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

### 1. Recuperar antes de recalcular

Conferir inventário, hashes e configuração. Se as entradas não mudaram, aproveitar a auditoria anterior. Ler [profiles.csv](../v16_spatial2/results/profiles.csv), [per_run_methods.csv](../v16_spatial2/results/per_run_methods.csv) e [equal_length_windows.csv](../v16_spatial2/results/equal_length_windows.csv). Identificar exatamente cada método, sem adivinhar nomes.

profiles.csv contém run, pattern, included, method, wavelength_px, zero_lag_wps, lag_mean_all, lag_mean_coi, lag_mean_supported, coi_count e supported_count. Conferir quais linhas correspondem ao método principal, as 200 escalas e sua cobertura. Primeiro reutilizar perfis existentes; recalcular apenas lacunas. **Perfil de registro inteiro não é o perfil do recorte específico por par usado no baseline v5.**

### 2. Congelar representação e métrica antes de avaliar os rótulos

Criar futuramente PROTOCOL.md, configuração e DECISIONS.md: hipótese, métrica principal, normalização, integração em escala, suporte, exclusões técnicas, controles e avaliação. Mudanças posteriores devem ter motivo e data, distinguindo exploração de confirmação.

Proposta inicial a confirmar por viabilidade, sem olhar separação por rótulos:

- Entrada principal: linhas 20–59, preparação v16, correlação enviesada, perfil de lag zero por aquisição inteira. Linhas 30–49 ficam como sensibilidade já explorada.
- Preservar valores do perfil e versão normalizada por área para comparar forma. A preparação já padroniza a variância do sinal; a normalização adicional remove amplitude do perfil intencionalmente.
- Definir faixa comum por amostragem, comprimento e suporte de todas as aquisições elegíveis. Começar verificando a faixa disponível no v16, sem presumir que toda ela seja adequada. Não escolher bandas pelos rótulos nem excluir bubbles porque seu registro é curto.
- Usar medida de integração explícita, inicialmente d(log lambda). Com pesos positivos de quadratura w_j, definir p_j = w_j P(0,a_j) / soma_k w_k P(0,a_k). São massas descritivas entre bandas, **não frações de energia física**. Área degenerada deve produzir registro inválido, não perfil artificial.
- Distância principal proposta: Hellinger, H(p,q) = sqrt(soma_j (sqrt(p_j)−sqrt(q_j))² / 2). Quanto menor, mais semelhantes as formas. Não atribuir essa fórmula a Schuck. Uma sensibilidade simples, como cosseno com pesos coerentes, pode verificar dependência da métrica; evitar busca extensa pela melhor separação.

Definir domínio comum antes de normalizar. Não comparar distribuições normalizadas em faixas diferentes como se fossem iguais. Verificar quadratura, perfis idênticos, invariância à amplitude após normalização e valores degenerados antes do lote.

### 3. Separar forma espectral de comprimento e posição

Comparar registros inteiros com um controle de comprimento comum definido sem consultar resultados. Aproveitar controles começo/centro/fim existentes quando adequados. Recortes da mesma aquisição não são novas réplicas. Não esticar eixos para coincidir picos.

Se houver forte sensibilidade a comprimento, posição, detrending ou suporte, quantificar antes de atribuir diferenças ao padrão. Estudar média em lag total, sob COI e sob suporte conservador como resultados secundários, com contagens de lags por escala. “Mean power” será média descritiva de alternative WPS, não potência física convencional.

### 4. Comparar com v5 e v16 de maneira justa

Manter referências separadas: v5 literal; v16 corrigido com entrada original; nova comparação de perfis. Para isolar representação e implementação, aplicar o mesmo resumo dos perfis a v5 e v16 em segmentos/eixos compatíveis, recalculando apenas o que faltar. Não misturar perfis de registros completos e recortes por par sem identificar a mudança.

Para comparar métodos com unidades diferentes, usar um mesmo desfecho de ordenação: por exemplo, frequência com que um par intrapadrão é mais semelhante que um entre padrões, com direção da métrica explícita. Não comparar numericamente distância de perfil 0,1 com WCO 0,2. Recuperar o escore WCO anterior sem reotimizá-lo.

Relatar todos os 325 pares, variação por padrão e combinação de padrões, agregados por par e balanceados por padrão. Definir pesos dos controles cruzados explicitamente. Mostrar se churn/annular dominam o resultado e se plug/slug perdem desempenho.

### 5. Avaliar incerteza e controles sem vazamento

A unidade experimental é a aquisição. Se houver ajuste de banda, métrica ou pré-processamento, separar ajuste e avaliação por aquisição; nenhum par contendo aquisição reservada pode entrar no ajuste. Validação deixando uma aquisição de fora é uma opção limitada pelo tamanho do conjunto. Os dados já foram vistos no v16: chamar de avaliação interna/exploratória, não teste externo intocado.

Permutações ou estimativas de incerteza devem respeitar aquisições e declarar pressupostos sobre rótulos/condições experimentais. Se não forem defensáveis, apresentar dispersão descritiva e reconhecer a limitação. Não tratar 325 pares como 325 experimentos independentes. Bubbles permite controles cruzados, mas não estimativa intrapadrão entre aquisições.

Controles sintéticos futuros: perfis idênticos; mudança global de amplitude; componentes em escalas distintas; sinais com espectro semelhante e organização/fases distintas. O último caso ajuda a interpretar perfis semelhantes com baixa alternative WCO. Esses controles não substituem comparações experimentais entre padrões.

### 6. Preservar as análises complementares e os casos desfavoráveis

Apresentar panorama, sinal, perfil de lag zero, média em lag, módulo da alternative WCS e alternative WCO. Selecionar exemplos principais por IDs e cobertura de padrões antes das métricas. Acrescentar casos de discordância — perfis semelhantes/WCO baixa, semelhança entre padrões e perdas — identificando que sua seleção posterior serve à interpretação, não à validação.

Reutilizar figuras/mapas existentes quando adequados, sempre indicando configuração e origem. Não é necessário recalcular todos os mapas apenas para ilustrar os perfis.

## Decisões a partir das descobertas

| Descoberta | Caminho justificável |
| --- | --- |
| Perfis existentes atendem ao protocolo | Reutilizar, conferindo proveniência. |
| Faixa comum curta ou perfis instáveis | Avaliar comprimento/suporte e cobertura; não escolher bandas por separação. |
| Ganho só após muitos ajustes | Relatar exploração e exigir avaliação separada. |
| Intrapadrão e controles seguem sobrepostos | Concluir ausência de melhora convincente. |
| Ganho restrito a um padrão | Relatar utilidade restrita e perdas nos demais. |
| Metadados físicos, rótulos ou escopo divergem | Consultar o usuário antes de redefini-los. |

Antes do lote completo, registrar o critério de relevância e justificar efeito mínimo/estabilidade desejados. Não inventar um limiar científico. Se for necessária uma decisão binária com ganho mínimo de aplicação, confirmar esse requisito com o usuário. Reportar efeitos absolutos, relativos, cobertura e controles. Sem evidência convincente, escrever **“não demonstramos melhora”**; com efeito pequeno, **“melhora marginal”**.

## Entregáveis futuros, tempo e encerramento

Após pedido de execução: scripts/configuração, protocolo, decisões, tabelas por aquisição/par, comparação com baselines, diagnósticos, manifesto de hashes/ambiente e instruções reproduzíveis. O resumo executivo deverá incluir os cinco padrões, bubbles e churn explicitamente, exemplos WCS/WCO e casos desfavoráveis; separar resultados herdados dos novos.

Antes de lote pesado, medir um pequeno benchmark com comprimentos distintos e projetar tempo, memória, variantes e controles. **Informar a estimativa ao usuário antes de iniciar**, atualizando se o plano mudar. O lote principal antigo levou ~5,8 minutos, mas não prevê esta proposta. A reutilização de perfis pode reduzir o custo; não prometer duração sem medir.

Encerrar respondendo sobre similaridade de escalas e sua utilidade frente às referências. Um ganho nessa representação não significa aumento de alternative WCO. Se não houver melhora, registrar o resultado negativo, sem migrar silenciosamente para métodos fora de wavelets. Proposta complementar: [v16.2](../v16.2_wco_eventos/PLAN.md).
