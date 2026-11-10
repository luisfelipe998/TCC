# Auditoria adversarial — dados, fórmulas, método e aplicação

A revisão confirmou a integridade das entradas e a implementação das fórmulas principais nos testes independentes realizados. Encontrou, porém, fragilidades de seleção, controles e interpretação, além de um defeito herdado na costura e imprecisões de documentação. Os resultados principais foram preservados; a conclusão permanece exploratória, sem melhora robusta nem identificação física demonstradas.

## Escopo, preservação e critério de leitura

A revisão abrangeu a etapa 1 (plug/slug), a extensão (churn/anular), a cadeia de montagem v11 e a extração 18–71. Os 264 hashes do manifesto anterior foram conferidos antes das alterações. Relatórios e verificadores anteriores foram copiados para [estado_anterior](results/auditoria_adversarial/estado_anterior); os dados e códigos congelados permaneceram intactos. O [plano dos diagnósticos](results/auditoria_adversarial/plano_diagnosticos.json) foi registrado antes das novas WCO desta revisão, com os resultados originais já conhecidos. Portanto, estes testes são exploratórios posteriores, não pré-registro ou validação externa.

O revisor foi o mesmo assistente de IA. “Independente” abaixo designa uma implementação matemática separada, não outro pesquisador. PASS indica integridade/consistência computacional; não certifica identidade de bolhas/ondas, rótulos físicos, estacionariedade ou generalização.

## Achados e tratamento

| ID | Gravidade / natureza | Achado | Tratamento e situação |
| --- | --- | --- | --- |
| A1 | Alta / interpretação | Correspondência próxima foi apresentada como estabilidade sem distinguir a seleção completa. | Corrigidos texto e gráfico; teste completo acrescentado. A repetibilidade continua limitada. |
| A2 | Alta / validade do método | Seleção por mínimos e largura também pode produzir ganho em sinais independentes. | Controle com procedimento completo acrescentado; ganho de WCO não é evidência suficiente de identidade física. |
| A3 | Moderada / amostragem | Poucos parceiros por par de aquisições tornam pequenas diferenças sensíveis à semente. | Todos os 1.095 pares intrapadrão admissíveis avaliados; 100 sementes alternativas por etapa. Original preservado. |
| A4 | Moderada / implementação herdada | Feather de dois pixels usa [0,1], duplicando uma coluna; o comentário sobre o sentido do deslocamento está invertido. | Reprodução e contrafactual quantificados; fonte v11 preservada. Correção física da montagem exige nova análise. |
| A5 | Moderada / relatório | A referência de cobertura de 2/3 migrou de aquisições para pares sem explicitar a mudança do denominador. | Texto corrigido e ambas as coberturas tabuladas; nenhum critério foi redefinido depois dos resultados. |
| A6 | Moderada / uso da fase | Refletir o lag da WCS legada sem conjugar não alinha sua fase à convenção principal. | Restrição documentada: comparação v5 por magnitudes/escores. Não houve uso experimental dessa fase nas conclusões. |

Não foi encontrado erro nas fórmulas principais, na faixa de extração ou nas agregações verificadas. A ausência de erro computacional não remove os achados acima nem garante ausência de falhas não detectadas.

## 1. Dados e cadeia de transformação

Foram relidos **31 conjuntos**, com **15.531 quadros**, seus rótulos locais e **31 planilhas**. A ordem e a sequência dos 501 nomes por aquisição foram verificadas. Todos os **15.500 deslocamentos** foram reestimados por Pearson em todos os lags previstos, sem chamar o pipeline v11; coincidiram com os caches. Os 31 panoramas foram reproduzidos pixel a pixel em memória. A média 18–71, calculada por soma inteira/divisor 54, coincidiu exatamente com os sinais da análise; as planilhas foram conferidas contra sua faixa original 20–59, sem confundi-las com a entrada atual. Não houve quadros inteiros repetidos por SHA-256. Foram registrados hashes de 15624 fontes.

Isso confirma a procedência computacional. Os rótulos foram conferidos contra os arquivos locais, sem uma nova classificação física independente. Os cinco registros excluídos foram auditados, mas suas exclusões não foram revertidas. JPEG, iluminação, paredes incluídas na média e reconstrução de uma dinâmica por faixas permanecem limitações. O eixo horizontal é **coluna do panorama montado**; interpretá-lo como distância percorrida exige hipóteses de advecção/deformação não validadas aqui.

### A4: costura e sinal do deslocamento

Para feather=2, `linspace(0,1,2)` gera [0,1]. A última coluna à esquerda da costura recebe a segunda coluna da faixa seguinte, que também permanece à direita. Não há média gradual entre os quadros nesses dois pontos. Com feather=1, o peso é zero e nada muda. A reprodução independente confirmou exatamente esse comportamento, inclusive os deslocamentos não positivos substituídos por 1 px.

O teste contrafactual retirou **somente o feather**, mantendo deslocamentos, âncoras, recortes, parceiros e a agregação. Foram recalculados 458 pares ancorados intrapadrão e 261 fixos:

| Padrão | Âncoras − fixos original | Sem feather | Mudança |
| --- | --- | --- | --- |
| Pistonado (plug) | 0,0289 | 0,0286 | -0,000241 |
| Golfadas (slug) | 0,2192 | 0,2187 | -0,000478 |
| Agitado (churn) | 0,0283 | 0,0280 | -0,000312 |
| Anular ondulado | 0,0044 | 0,0046 | 0,000141 |

O efeito foi pequeno nesses contrastes. Isso não demonstra que os deslocamentos, a evolução das interfaces ou toda a reconstrução sejam fisicamente corretos. Tampouco inclui redetecção das âncoras no panorama sem feather. Não se substituiu a entrada principal depois de conhecer os resultados.

Um deslocamento sintético real de **+7 px à direita** retorna **−7** na função de Pearson do v11, ao contrário de seu comentário. Isso é uma inconsistência de convenção documental: não prova que o escoamento real se desloca no sentido errado e não justifica inverter os dados. A regra de limitar deslocamentos não positivos a 1 px continua uma hipótese operacional, particularmente relevante em m03q40_001 (60/500).

Evidências: [dados por aquisição](results/auditoria_adversarial/dados_por_aquisicao.csv), [auditoria dos dados](results/auditoria_adversarial/dados.json), [fontes](results/auditoria_adversarial/fontes_sha256.json), [contrastes sem feather](results/auditoria_adversarial/contrastes_sem_feather.csv), código de origem [pipeline_v11.py](../pipeline_v11.py).

## 2. Fórmulas e implementação

O [capítulo original de Schuck Jr. e Bodmann](../../assets/schuck-bardo-wavelet-coherence.pdf), páginas impressas 411, 413 e 414, foi conferido visualmente. A CWT, a Morlet e a conversão de escala correspondem às Eqs. 32.4, 32.6 e 32.8. As WPS usam módulo e a WCO usa o denominador no lag zero, como 32.14–32.17. O divisor 1/N é uma **adaptação enviesada explícita**; o artigo usa 1/(N−|m|) em 32.21–32.23. Não se afirma reprodução literal desse estimador. A remoção de média/tendência também é estimada no próprio recorte.

Foram selecionados o primeiro e o último ID intrapadrão em cada um dos quatro padrões. Nos oito pares, o preparo foi refeito por regressão independente e a correlação por soma direta. A CWT foi conferida por produto escalar explícito com a Morlet temporal em cinco escalas e sete lags, incluindo extremos. A WCO foi conferida nos 65 lags do indicador, nessas cinco escalas. Também foram verificados troca x/y, preenchimento de zeros e escore completo salvo. Maior diferença da CWT: **5,40e-15**; nenhuma discrepância material. Essa amostragem independente não é recálculo por soma direta de todas as células dos milhares de mapas.

O limite unitário da Morlet principal tem fundamento adicional: com correlação 1/N, o multiplicador espectral H é real e não negativo. A WCS é um produto interno ponderado de X e Y deslocado; Cauchy–Schwarz limita seu módulo ao quadrado pelo produto das energias ponderadas no lag zero. Essa justificativa não se transfere automaticamente ao estimador 1/(N−|m|) nem ao kernel discreto legado. As três ultrapassagens já documentadas do v5 (máximo 1,006157, todas fora de ±32) permanecem preservadas.

### A6: lag alinhado não é fase alinhada

Para sinais reais, trocar Rxy por Ryx reflete o lag e conjuga a CWT. O comparador v5 reflete o mapa depois da CWT para orientar o lag, preservando as magnitudes utilizadas. Essa reflexão sem conjugação **não alinha a fase complexa** à principal. No exemplo testado, o erro relativo frente à CWT legada da Rxy canônica foi 1,4132; ao conjugar caiu para 0,0052, com diferença residual da discretização. Não se deve interpretar essa fase refletida como atraso/fase na convenção principal. A fórmula da Morlet principal não apresenta esse problema.

A ambiguidade do artigo entre a média complexa (32.16) e a média de módulos impressa (32.20) já estava identificada na base e foi confirmada. Com K=1, não muda as magnitudes/WCO; com K>1, pode mudar radicalmente. Duas realizações de fases opostas se cancelam na média complexa e não na média de módulos. As médias descritivas entre pares desta entrega **não são WCO de ensemble**.

Evidências: [conferência das fórmulas](results/auditoria_adversarial/formulas.json), [erros por par](results/auditoria_adversarial/formulas_comparacoes.csv), [auditoria anterior das equações](../v16_spatial2_18_71/EQUATION_AUDIT.md).

## 3. Âncoras, seleção e parceiros

### A1: duas noções diferentes de estabilidade

O teste original procura o candidato do mesmo tipo mais próximo do centro congelado, até 64 px. Ele não repete a escolha de maior proeminência por estrato. Refazer o procedimento completo produz:

| Padrão | Correspondência local σ12 / σ20 | Seleção completa σ12 / σ20 | Seleção completa: linhas 30–49 |
| --- | --- | --- | --- |
| Agitado (churn) | 40/40 em ambos | 37/40; 38/40 | 27/40 |
| Anular ondulado | 36/36 em ambos | 29/36; 30/36 | 28/36 |

![Repetibilidade da seleção](figures/auditoria_adversarial/repetibilidade_selecao.png)

As coordenadas dos candidatos correspondentes podem variar até 11 px e, ainda assim, outro candidato vencer a seleção. Quando isso ocorre, não é correto chamar a distância ao novo selecionado de deslocamento físico da mesma estrutura. O novo teste confirma repetibilidade parcial; não invalida a existência das âncoras originais nem altera seu catálogo. Não foram recalculadas todas as WCO de um catálogo totalmente reselecionado: essa é uma análise futura, com novos parceiros e nova interpretação de cobertura.

### A3: sensibilidade ao limite de parceiros

O hash determinístico evita escolher pares pela maior WCO, mas não elimina a variabilidade causada por usar no máximo 12 pares na etapa 1 e 2 na extensão. Foram avaliados **1.095 pares intrapadrão admissíveis** de compatíveis/equivalentes e distintos, sendo 799 novos; os demais escores foram reutilizados. Mantiveram-se as âncoras. Foram aplicadas a semente original e 100 sementes consecutivas adicionais, preservando pesos por par de aquisições e as interseções de comparação:

| Padrão | Controle | Mínimo / máximo em 100 sementes | Sementes com diferença positiva |
| --- | --- | --- | --- |
| Pistonado (plug) | Fixos | 0,0137 / 0,0395 | 100/100 |
| Pistonado (plug) | Larguras/extensões distintas | 0,0018 / 0,0458 | 100/100 |
| Golfadas (slug) | Fixos | 0,2192 / 0,2192 | 100/100 |
| Golfadas (slug) | Larguras/extensões distintas | 0,1520 / 0,1520 | 100/100 |
| Agitado (churn) | Fixos | 0,0120 / 0,0477 | 100/100 |
| Agitado (churn) | Larguras/extensões distintas | -0,0113 / 0,0432 | 86/100 |
| Anular ondulado | Fixos | -0,0137 / 0,0258 | 80/100 |
| Anular ondulado | Larguras/extensões distintas | -0,0711 / 0,0156 | 4/100 |

![Sensibilidade aos parceiros](figures/auditoria_adversarial/sensibilidade_parceiros.png)

Usar todos os parceiros admissíveis resultou em ganhos frente a fixos de **0,0239 em plug**, **0,2192 em slug**, **0,0298 em churn** e **0,0086 em anular**. Slug continua sustentado por um único par equivalente. A faixa entre sementes mede decisões de amostragem dentro do mesmo catálogo, não incerteza populacional. Não se escolheu uma nova semente ou limite que favorecesse a conclusão.

### A5: denominadores de cobertura

No contraste com mínimos distintos, churn tem **10/10 aquisições**, mas **29/45 pares possíveis**; anular tem **8/9 aquisições**, mas **13/36 pares**. A etapa 1 definiu sua referência de 2/3 por aquisição. O texto da extensão não deve transformá-la em um critério por pares sem declarar a mudança. O resumo foi corrigido; os dois denominadores aparecem em [cobertura_denominadores.csv](results/auditoria_adversarial/cobertura_denominadores.csv). Ter todas as aquisições representadas não elimina a escassez de parceiros comparáveis.

## 4. A2: controle de todo o procedimento

Os controles sintéticos anteriores testavam formas e o núcleo, mas não isolavam o efeito de escolher extremos em sinais independentes. Agora foram gerados **40 pares independentes de AR(1), φ=0,95**, com 8.192 pontos após aquecimento de 1.024, intensidade média 128 e desvio 20. Foram aplicados o detector, a seleção por estrato, a compatibilidade de largura, o limite de dois parceiros e o contraste com três janelas fixas, antes de consultar a WCO. Todos os 40 casos foram calculáveis; não houve seleção dos casos de maior efeito.

A média ancorada foi **0,2309**, contra **0,2140** dos fixos; diferença **0,0169**, positiva em **25/40** repetições. Assim, mínimos escolhidos em séries independentes de estrutura estatística semelhante também podem gerar ganho. O resultado não prova que todo ganho experimental seja artefato, nem constitui um nulo calibrado para esses dados. Não foram calculados valores-p, nem subtraído esse ganho sintético dos contrastes reais.

![Controle sintético completo](figures/auditoria_adversarial/nulo_selecao.png)

## 5. Aplicação, dependência e limites que permanecem

As agregações principais foram conferidas por par de aquisições e interseção dos controles. Não houve agrupamento silencioso de plug/slug com churn/anular, nem conversão das 15 ausências de σ20 em zero. As variantes do detector e o painel comum de 815 pares preservam os parceiros originais; não representam a variabilidade de um novo catálogo completo.

Pares compartilham aquisições, recortes se sobrepõem e a CWT usa correlações de todo o contexto. A WCO não localiza automaticamente o evento central. A compatibilidade de largura não controla completamente iluminação, contraste, periodicidade, outras estruturas ou costuras. Mínimo × máximo não é um negativo obrigatório, pois a WCO é invariável à inversão global. As janelas fixas tampouco são necessariamente livres de eventos.

Os 200 comprimentos de onda e 65 lags não são réplicas independentes. O piso de denominador é numérico; o suporte geométrico e o critério de dois comprimentos de onda por segmento não certificam potência, estacionariedade ou significância. O limiar de ganho 0,05 é operacional e descritivo. Não há validação externa, revisão física independente, calibração espaço/tempo ou desempenho de classificador. Esses limites já existiam; os novos testes reforçam sua relevância.

## 6. Alterações e próximas decisões científicas

Foram atualizados resumo executivo, PDF integrado, passo a passo, resultados, nota de âncoras/autocorrelação, README, registro de decisões e estado de retomada. A figura de estabilidade original passou a explicitar que trata de correspondência local. As tabelas principais, mapas, protocolos congelados e bases anteriores foram preservados. Os relatórios históricos permanecem como registros da etapa correspondente, com o estado atual indicado nos documentos principais.

Uma análise posterior deverá definir claramente o estimando: repetibilidade de marco óptico, semelhança do contexto, identidade física ou classificação são objetivos diferentes. Para afirmar identidade física, será necessária anotação independente e rastreamento/validação da montagem. Para avaliar estabilidade operacional completa, será necessário refazer seleção e parceiros em cada configuração. Uma nova montagem deve corrigir e testar a costura e a convenção de deslocamento antes de produzir outra versão dos dados. Mais combinações dos mesmos registros não substituem novas aquisições.

Os diagnósticos WCO desta revisão consumiram **3,0 minutos**; leitura de dados, conferência matemática e redação são etapas separadas. As tabelas contrafactuais exportam somente os escores recalculados, para não misturar outras métricas herdadas; a [revisão de exportação](results/auditoria_adversarial/revisao_exportacao.json) registra que os escores permaneceram exatamente iguais.

Confira a [verificação adversarial](results/auditoria_adversarial/verificacao_final.json), a [verificação integrada](results/verificacao_final.json) e o [manifesto atualizado](results/manifesto.json). Os [scripts de auditoria dos dados](scripts/adversarial_data_audit.py), [âncoras](scripts/adversarial_anchor_audit.py), [fórmulas](scripts/adversarial_math_audit.py) e [diagnósticos](scripts/adversarial_diagnostics.py) registram as operações executadas.
