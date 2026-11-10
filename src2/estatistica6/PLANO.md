# Plano de análises estatísticas complementares

**Estado: planejamento. Nenhuma das análises propostas foi executada nesta etapa.**

Data: 24 de outubro de 2026. Este documento atende à solicitação de preparar um plano e registrar informações suficientes para continuidade por outro agente. A instrução vigente é **não executar as análises**. Uma solicitação posterior definirá quais blocos serão realizados; a existência deste plano não autoriza sua execução automática.

## 1. Objetivo e evidência disponível

Investigar a estabilidade dos descritores e a informação adicional presente nos resultados de v16_spatial2_18_71. Uma conclusão relevante pode ser uma limitação do método, uma associação operacional ou ausência de evidência de separação. Não há expectativa obrigatória de significância ou de melhora da discriminação.

Os números abaixo foram transcritos dos resultados existentes, sem recálculo:

| Resultado de referência | Valor |
| --- | --- |
| Ensaios elegíveis | 26: 10 churn, 9 wavy annular, 3 plug, 3 slug e 1 dispersed bubbles |
| Pares distintos | 325: 87 intrapadrão e 238 interpadrão |
| Mediana de Pearson intrapadrão | 0,830124 |
| Mediana de Pearson interpadrão | 0,754573 |
| Diferença entre medianas | 0,075551 |
| Pares interpadrão que atingem a mediana intrapadrão | 79 de 238 |
| Alteração absoluta de r entre prefixo e registro completo | Mediana 0,101426; máximo 0,703985 |
| Maior diferença contra SciPy na conferência anterior de Pearson | Aproximadamente 1,11 × 10⁻¹⁵ nas duas condições |

Fontes: [resumo numérico atual](../estatistica5/todos_contra_todos/results/summary.json), [validação existente](../estatistica5/todos_contra_todos/results/validation.json) e [resumo executivo](../estatistica5/RESUMO_EXECUTIVO.md).

O quadro atual indica sobreposição entre padrões e sensibilidade ao trecho. Correção numérica não demonstra classificação, repetibilidade experimental ou causalidade. A [estatistica4](../estatistica4/RELATORIO.md) já contém PCA, associações com covariáveis, controles de janelas, influência por ensaio e descrição de dispersão. As extensões propostas não devem ser apresentadas como se esses controles nunca tivessem sido realizados.

Este é um plano exploratório de dados já examinados, não um pré-registro prospectivo. Registrar escolhas antes dos próximos cálculos reduz decisões oportunistas, mas não elimina o conhecimento prévio da amostra.

## 2. Prioridades e dependências

| ID | Análise proposta | Prioridade | Dependência principal |
| --- | --- | --- | --- |
| A1 | Estabilidade do contraste intra/inter e influência por ensaio | Alta | Matriz atual de Pearson |
| A2 | Estabilidade espacial e sensibilidade ao comprimento | Alta | Perfis de janelas já salvos |
| A3 | Características espectrais interpretáveis | Alta | Perfis e coeficientes atuais |
| A4 | Comparação funcional dos perfis WPS | Média | Curvas por ensaio e justificativa de permutação |
| A5 | Distâncias Wasserstein e Hellinger | Média | Normalização definida em A3 |
| A6 | Associação com condições operacionais e reconstrução | Alta como diagnóstico | Características de A3 e inventário |
| A7 | Dispersão dentro dos padrões | Média | Vetores unitários; auditoria adicional para PERMDISP |
| A8 | Ampliação dos controles nulos da WCO | Baixa inicialmente; custo elevado | Modelo nulo e equações da v16 |

Ordem recomendada quando houver solicitação de execução: conferir fontes; realizar A1 e A2; consolidar A3 e o diagnóstico A6; selecionar os complementos necessários. A8 responde a outra pergunta e não é condição para interpretar Pearson.

ANOSIM permanece fora das novas propostas. Não estão previstos novos classificadores, reconstrução de panoramas ou alteração das equações de Schuck. PERMDISP aparece somente como possibilidade específica de A7, condicionada à auditoria metodológica indicada adiante.

## 3. Delineamento comum

### 3.1 Unidade de análise e participantes

- O **ensaio** é a unidade amostral. Pares, pixels, escalas e janelas não são aquisições independentes.
- A referência descritiva contém os 26 ensaios e um descritor fixo por ensaio, calculado nos primeiros 3.171 px. Não usar um trecho que mude conforme o parceiro do par.
- Quando a inferência for justificável, priorizar o subconjunto de 10 churn e 9 wavy annular, mantendo inicialmente os mesmos descritores de 3.171 px. Seus 171 pares compreendem 81 intrapadrão e 90 interpadrão. Não atribuir a esse subconjunto os resultados agregados dos 26 ensaios.
- O controle com 7.001 px usa os mesmos 19 ensaios. A comparação entre comprimentos não deve mudar simultaneamente os participantes.
- Plug e slug, com três ensaios cada, terão interpretação descritiva ou exploratória limitada. Dispersed bubbles, com um registro, não permite estimar variabilidade intraclasse ou comportamento populacional da classe.
- Permanecem excluídos os cinco registros da origem: m03q01_001, m03q80_001, m04q05_001, m05q05_002 e m07q10_001. Não alterar elegibilidade pelo resultado.

As unidades espaciais são pixels do panorama reconstruído. Não converter os resultados em frequência temporal, milímetros ou comprimento de bolhas sem calibração específica. Os campos m e q são códigos operacionais; sua correspondência e suas unidades físicas precisam ser verificadas antes de interpretação física.

### 3.2 Perfil e descritor preservados

A extração é a média das linhas 18–71, inclusivas, com índices começando em zero:

$$s_i[n]=\frac{1}{54}\sum_{h=18}^{71}I_i[h,n].$$

Selecionar o trecho **antes** de inverter a intensidade, retirar a tendência linear e padronizar pelo desvio-padrão populacional. Para o sinal preparado de comprimento N:

$$\widehat R_{ii}[m]=\frac{1}{N}\sum_{n=0}^{N-1-m}x_i[n+m]x_i[n],\qquad 0\le m<N.$$

Usar extensão bilateral, Morlet com frequência angular central 6 e:

$$\lambda=\frac{4\pi a}{6+\sqrt{38}},\qquad P_i(a)=\left|\operatorname{CWT}\{\widehat R_{ii}\}(0,a)\right|.$$

Manter 200 comprimentos de onda logaritmicamente espaçados entre 50 e 3.000 px. P é a magnitude definida na implementação da alternative WPS; não elevar novamente ao quadrado nem denominá-la automaticamente densidade de energia física.

$$z_{ik}=\log(P_i(a_k)+10^{-12})-\widehat c_i-\widehat d_i\log\lambda_k,$$

$$u_i=\frac{z_i-\overline z_i}{\|z_i-\overline z_i\|_2},\qquad r_{ij}=u_i^\mathsf{T}u_j,\qquad d_{ij}=\|u_i-u_j\|_2=\sqrt{2(1-r_{ij})}.$$

Recusar vetores constantes, degenerados ou não finitos. A retirada da reta log-linear não garante branqueamento ou independência entre escalas. A quantidade 1−r é uma dissimilaridade; para operações euclidianas, usar u ou a distância indicada acima.

A banda integral permanece como referência, sem novo COI ou máscara. Escalas próximas de 3.000 px em trechos de 3.171 px têm suporte finito limitado. Uma banda alternativa só poderá ser um controle identificado, com justificativa e limites fixados antes do cálculo; não substituir a referência pelo maior contraste.

Definições completas: [metodologia atual de Pearson](../estatistica5/METODOLOGIA_TODOS_CONTRA_TODOS.md) e [auditoria das equações da v16](../v16_spatial2_18_71/EQUATION_AUDIT.md).

### 3.3 Dependência, inferência e multiplicidade

1. Não aplicar testes de duas amostras ou bootstrap de linhas às correlações de pares como se fossem independentes. Não calcular valor-p paramétrico de Pearson usando as 200 escalas como tamanho amostral.
2. Quando um teste usar rótulos, permutar o rótulo do ensaio inteiro, mantendo seu vetor intacto. Para churn/wavy annular, a estatistica4 identificou 200 arranjos distintos dentro de q. Conferir participantes e contagens antes de reutilizar o esquema.
3. As contagens históricas churn/wavy annular por q são: q=5: 1/0; q=10: 3/1; q=20: 4/1; q=40: 2/3; q=80: 0/4. Estratos puros não trocam rótulos. Dentro de (m,q), existe apenas um arranjo, sem teste de permutação possível para controlar ambos simultaneamente.
4. Permutação restrita não equivale a residualizar covariáveis. Sua validade depende de permutabilidade dentro dos estratos, não garantida pelo delineamento observacional. Sem justificativa suficiente, entregar efeitos descritivos e declarar o teste não realizado.
5. Para enumeração completa, incluindo o observado, usar p = número de estatísticas pelo menos tão extremas quanto a observada dividido pelo número de arranjos. Para Monte Carlo, usar p=(1+b)/(1+B). Não reportar p=0. Registrar tolerância de empates.
6. Antes de uma execução inferencial, fixar quais hipóteses de A1, A3, A4 e A7 serão testadas. Se integrarem a mesma execução, aplicar Holm conjuntamente aos valores-p globais. A4 também controla internamente a busca entre escalas. Não acrescentar testes após selecionar a sensibilidade mais favorável.
7. A8 forma uma família própria, por testar outro objeto. Conclusões que dependam conjuntamente de famílias diferentes exigem explicitar o controle adicional de multiplicidade.
8. Se houver inferência, usar 0,05 apenas como referência nominal. Relatar efeitos, incerteza e todos os testes, inclusive os não significativos. Correção de multiplicidade não transforma esta exploração em confirmação independente.

## 4. A1 — Estabilidade do contraste intra/inter

**Pergunta:** o contraste observado depende de poucos ensaios ou da ponderação das classes?

**Entrada:** matriz, pares e índice da estatistica5 atual. Não exige nova transformada wavelet.

$$\Delta=\operatorname{mediana}\{r_{ij}:g_i=g_j,i<j\}-\operatorname{mediana}\{r_{ij}:g_i\ne g_j,i<j\}.$$

Nos 26 ensaios, a referência arquivada é 0,07555100435301754. O valor do subconjunto de 19 deverá ser calculado separadamente quando autorizado.

**Procedimento principal:** reproduzir o resumo arquivado; retirar um ensaio por vez, mantendo os descritores e o comprimento fixos; calcular Delta sem esse ensaio e a alteração relativa à referência. Mostrar resultados por combinação de padrões e contagens.

Acrescentar uma agregação balanceada nos quatro padrões com réplicas: média das quatro medianas intrapadrão menos média das seis medianas interpadrão entre essas mesmas quatro classes. Essa ponderação responde a outra pergunta e não substitui a mediana global. Mostrar a remoção do singleton separadamente como sensibilidade à composição, conservando a referência com 26.

**Inferência opcional:** teste direto de Delta nos 19 ensaios, com alternativa unilateral de maior semelhança intrapadrão e enumeração dentro de q, se a permutabilidade for justificável. ANOSIM não será executado.

**Reamostragem opcional:** começar por perturbação de pesos de ensaio como análise de estabilidade. Gerar pesos positivos independentes com distribuição exponencial de média um, normalizá-los dentro de cada uma das duas classes para média um e atribuir a cada par distinto o produto dos pesos de seus ensaios. Calcular medianas ponderadas intra/inter em 5.000 repetições, com semente registrada.

Definir a convenção de mediana ponderada antes do cálculo. Uma opção é a inversa da distribuição acumulada em 0,5; quando a massa acumulada for exatamente 0,5, usar a média dos valores adjacentes, reproduzindo a mediana usual com pesos unitários.

Os quantis dessa perturbação são **intervalos de sensibilidade à ponderação**, não intervalos de confiança frequentistas calibrados. Não chamar esse procedimento de bootstrap comum. Caso sejam solicitados intervalos populacionais, definir e validar um método para pares dependentes antes de usá-lo. Duplicações de ensaios não podem inserir artificialmente r=1 como correlação entre aquisições distintas. O singleton não fornece incerteza populacional de sua classe.

**Validação:** invariância à ordenação; exclusão da diagonal e dos pares invertidos; remoção de todas as arestas incidentes no ensaio omitido; conferência das medianas ponderadas; enumeração sem duplicatas. Faixa leave-one-out não é intervalo de confiança.

**Entregas futuras:** tabelas de contrastes, influência por ensaio e ponderações; distribuição das estatísticas permutadas, se houver teste; gráfico de influência e quadro das diferentes agregações.

## 5. A2 — Estabilidade espacial e comprimento

**Pergunta:** quanto o descritor muda ao observar outra região do mesmo panorama?

**Entrada:** perfis já salvos em equal_length_windows.csv, sinais raw e descritores de 7.001 px da estatistica4.

Para L=3.171 e comprimento original N, usar início 0, piso((N−L)/2) e N−L, com fim exclusivo início+L. Conferir coordenadas e deduplicar janelas idênticas. No único dispersed bubbles, prefixo, centro e sufixo coincidem; não existem três observações distintas.

1. Reutilizar os perfis arquivados após conferência metodológica. Recalcular somente o que faltar ou não tiver proveniência verificável.
2. Calcular resíduos, correlações e distâncias entre janelas distintas do mesmo ensaio, além da variação das características de A3.
3. Construir uma matriz completa entre ensaios para cada posição, sempre com um descritor fixo por ensaio naquela condição. Reportar Delta de todas as posições, inclusive mudanças de sinal.
4. Registrar a fração de sobreposição entre janelas, dada pelo comprimento de sua interseção dividido por L. Janelas sobrepostas não demonstram repetibilidade entre aquisições.
5. Comparar 3.171 e 7.001 px nos mesmos 19 ensaios. Separar esse controle da comparação entre posições com comprimento igual. Registros completos continuam um controle com conteúdo e comprimentos diferentes.

**Bootstrap espacial não integra o bloco inicial.** Reamostrar e emendar blocos pode criar descontinuidades e alterar o espectro. Uma extensão exigiria justificar estacionariedade local, tamanho dos blocos, preservação das escalas relevantes e viés introduzido. Sua conclusão se limitaria à variabilidade condicional do registro.

**Validação:** corte antes do preparo; mesma grade; coordenadas e sobreposição corretas; deduplicação do singleton; conferência entre perfis em cache e uma implementação independente. Não selecionar a melhor janela como nova referência.

**Entregas futuras:** coordenadas, sobreposição, estabilidade por ensaio, contraste por condição; gráficos conectados por ensaio, perfis por posição e comparação entre matrizes.

## 6. A3 — Características espectrais interpretáveis

**Pergunta:** há informação em escala dominante, inclinação ou concentração que Pearson dos resíduos não conserva?

**Entrada:** profile, coefficients e wavelengths_px do NPZ atual. Picos e inclinações já foram descritos anteriormente; a extensão será um conjunto pequeno e definido de características, com efeitos e estabilidade.

### Normalização

Definir o suporte logarítmico \(t_k=\log(\lambda_k/1\,\mathrm{px})\) e pesos trapezoidais h: metade do intervalo adjacente nas extremidades e metade da distância entre os dois vizinhos nos pontos interiores. Então:

$$p_{ik}=\frac{P_{ik}h_k}{\sum_\ell P_{i\ell}h_\ell}.$$

Essa é uma distribuição discreta do **perfil** no eixo logarítmico, não probabilidade de padrão nem fração física de energia. Recusar soma nula ou não finita. Não calcular essas probabilidades dos resíduos, que podem ser negativos.

| Característica | Definição | Papel proposto |
| --- | --- | --- |
| Escala dominante | Comprimento de onda do máximo de P; usar seu log na comparação | Principal |
| Inclinação log-linear | Coeficiente d do ajuste existente | Principal |
| Largura central | Quantil 0,75 menos quantil 0,25 de t, com pesos p | Principal |
| Centro logarítmico | Soma de p_k t_k | Complementar descritiva |
| Entropia normalizada | −soma(p_k log p_k)/log K, com 0 log 0=0 | Complementar descritiva |
| Participação por bandas | Soma dos p_k da banda | Complementar descritiva |

Definir quantis pela inversa da distribuição acumulada, sem interpolação não declarada. Registrar empates no máximo, múltiplos máximos e picos na borda; não atribuir precisão subgrade ao pico.

Bandas candidatas: [50,200), [200,800) e [800,3000] px. São convenções exploratórias deste plano, não fronteiras físicas validadas. Não alegar escolha completamente independente dos dados já examinados. As três participações somam um e não são variáveis independentes.

**Comparação:** mostrar todos os ensaios e medianas por classe. Nos 19 ensaios, reportar diferenças de medianas e, opcionalmente, delta de Cliff:

$$\delta_C=\frac{\#\{f_i>f_j\}-\#\{f_i<f_j\}}{n_Cn_A},\qquad i\in churn,\ j\in wavy\ annular.$$

Os produtos de pares dessa fórmula não se tornam observações independentes. Se houver inferência, permutar rótulos de ensaio dentro de q para o valor absoluto da diferença de medianas das três características principais. Aplicar multiplicidade conforme seção 3.3. As complementares permanecem descritivas nesta versão.

**Validação:** probabilidades somam um; quantis ordenados; entropia entre zero e um; entropia zero para massa pontual e um para probabilidades uniformes; invariância da normalização à multiplicação positiva; exatamente uma banda por ponto da grade. Distribuição uniforme de probabilidades não é necessariamente perfil P constante, devido aos pesos de quadratura.

**Entregas futuras:** características por ensaio, efeitos e testes quando justificados; pontos individuais por padrão e estabilidade das características nas janelas de A2.

## 7. A4 — Comparação funcional dos perfis

**Pergunta:** diferenças localizadas em escala ficam ocultas pela correlação global?

**Condição proposta:** 19 ensaios, primeiros 3.171 px, curva \(f_{ik}=\log(P_{ik}+10^{-12})\), antes de remover sua inclinação. Resíduos poderão ser exibidos como controle descritivo, sem segundo teste automático.

Calcular diferença média entre curvas e estatística global:

$$t_k=\frac{\overline f_{Ck}-\overline f_{Ak}}{\sqrt{s^2_{Ck}/n_C+s^2_{Ak}/n_A}},\qquad T=\max_k|t_k|.$$

Recalcular T nos arranjos dentro de q, preservando cada curva inteira. Trata-se de teste global condicional, não de 200 testes t paramétricos. A distribuição dos máximos contempla a busca entre escalas sob a hipótese nula de permutabilidade.

Se houver valores por escala, usar a distribuição dos máximos como referência e explicitar o alcance do controle. Não prometer controle forte sob alternativas parciais sem verificar as hipóteses adicionais. Um envelope nulo de permutação não é intervalo de confiança da diferença populacional. Não converter o resultado em contornos de significância de WCO.

**Validação:** permutar linhas inteiras; manter grade completa; conferir enumeração; detectar denominadores nulos ou quase nulos. Se houver degeneração, revisar e documentar a estatística antes do teste, sem eliminar escalas pelo resultado.

**Entregas futuras:** diferenças de curvas, estatísticas permutadas e metadados; gráficos com dispersão dos ensaios, diferença por escala e referência nula claramente identificada.

## 8. A5 — Wasserstein e Hellinger

**Pergunta:** perfis com Pearson elevado diferem por deslocamento dos picos ou redistribuição de sua magnitude relativa?

Usar probabilidades p e suporte t definidos em A3. Para a distribuição acumulada F:

$$W_1(i,j)=\int|F_i(t)-F_j(t)|\,dt,$$

$$H(i,j)=\sqrt{\frac12\sum_k(\sqrt{p_{ik}}-\sqrt{p_{jk}})^2}.$$

Wasserstein estará em unidades logarítmicas de comprimento de onda, não pixels. Hellinger varia de zero a um. Ambas descartam amplitude global pela normalização; não recuperam intensidade bruta perdida no preparo.

Calcular matrizes 26×26, resumos por combinação e influência por ensaio. Para distâncias, o contraste favorável à semelhança intrapadrão é mediana inter menos mediana intra, invertendo a subtração usada com r.

Nesta versão, ambas são diagnósticos descritivos. Não selecionar a métrica pelo menor valor-p. Spearman dos resíduos já existe na estatistica4 e não será anunciado como novidade.

**Validação:** simetria, diagonal zero, não negatividade, H≤1 e distância zero para distribuições idênticas. Massas pontuais em t_a e t_b devem ter W1=|t_a−t_b|; H=1 em suportes disjuntos. Conferir W1 com scipy.stats.wasserstein_distance, fornecendo suportes e pesos explicitamente, e H por identidade algébrica independente.

**Entregas futuras:** matrizes, resumos e gráficos de comparação com Pearson. Associação entre entradas de matrizes permanece descritiva; suas entradas não constituem amostras independentes.

## 9. A6 — Condições operacionais e reconstrução

**Pergunta:** características se associam às condições ou à reconstrução, além dos rótulos visuais?

**Entrada:** características de A3, m_code, q_code, comprimento original e recomputed_shift_mean. Consultar as associações já existentes na estatistica4 para identificar o que será ampliado.

1. Conferir cobertura de condições por classe, significado dos códigos, campos ausentes e colinearidade. Comprimento e deslocamento médio tiveram correlação 1,000000 nos 19 ensaios da estatistica4; não incluí-los juntos como preditores independentes.
2. Ampliar gráficos de associação com identificação de ensaio e padrão. Igualar o número de pixels não elimina associação com condições operacionais ou reconstrução.
3. Se a codificação numérica for justificável, usar modelos pequenos e previamente fixados, uma característica por resposta. Exemplo: \(f_i=\beta_0+\beta_m m_i+\beta_q q_i+\epsilon_i\). Linearidade nos códigos é hipótese do modelo, não consequência dos nomes dos arquivos. Se não for justificável, manter estratificação e descrição.
4. Usar modelo separado com deslocamento da reconstrução como diagnóstico. Não acrescentar automaticamente interações, transformações escolhidas por ajuste ou todos os descritores como preditores.
5. Um coeficiente adicional de padrão só será estimado se a matriz de desenho tiver posto suficiente e cobertura comparável. A associação continuará condicional ao modelo; extrapolação entre estratos exclusivos não identifica efeito físico do padrão.
6. Inspecionar resíduos, condicionamento, influência e sensibilidade leave-one-out. Usar mínimos quadrados numericamente estáveis, sem inversão explícita da matriz normal. Priorizar coeficientes, suporte observado e instabilidade, sem seleção automática por significância.

**Limite:** regressão não cria réplicas nem resolve ausência de sobreposição. Não atribuir causalidade, explicação exclusiva ao padrão ou origem exclusivamente óptica às associações.

**Entregas futuras:** cobertura de condições, associações e diagnósticos dos modelos aplicáveis; gráficos por ensaio e restrições de identificação.

## 10. A7 — Dispersão intrapadrão

**Pergunta:** um padrão apresenta perfis mais heterogêneos, mesmo quando as médias se sobrepõem?

**Etapa descritiva:** nos vetores unitários u, calcular centroide \(c_g=n_g^{-1}\sum_{i\in g}u_i\), raios \(a_i=\|u_i-c_g\|_2\) e ajuste descritivo \(a_i^*=\sqrt{n_g/(n_g-1)}a_i\). Mostrar raios individuais, mediana e média. O ajuste não torna cada ensaio independente do centroide estimado nem cria um teste. Dispersed bubbles ficará não estimável.

A estatistica4 já descreveu raios dos 19 ensaios de 7.001 px. O bloco proposto usa 3.171 px; distinguir os comprimentos.

**PERMDISP opcional:** restringir a churn/wavy annular. Antes de implementar, consultar uma implementação de referência de betadisper/permutest.betadisper; fixar centroide versus mediana espacial, ajuste de amostra e algoritmo de permutação; verificar a possibilidade e a validade da restrição por q. Os 200 arranjos de A1 não validam automaticamente seu uso nesse teste.

Não improvisar o teste permutando vetores brutos e recalculando centroides: isso pode misturar diferenças de localização e de dispersão. Usar algoritmo de permutação de resíduos ou procedimento equivalente validado para a hipótese declarada. **A versão formal permanece pendente dessa auditoria metodológica.** Se não houver justificativa para o algoritmo e a restrição, entregar somente a descrição.

**Validação:** comparar com implementação de referência; controles sintéticos com centros diferentes e dispersão igual, centros iguais e dispersões diferentes, tamanhos desiguais e heterogeneidade entre estratos. Avaliar comportamento agregado de repetições, sem exigir decisão específica em uma amostra simulada.

**Entregas futuras:** raios por ensaio, resumo por classe, auditoria do método e teste somente se fundamentado; gráficos de dispersão. Diferença de dispersão não constitui acurácia de classificação.

## 11. A8 — Modelo nulo para máximos de WCO

**Pergunta:** máximos excedem o esperado sob sinais estacionários com fases aleatórias independentes e espectros semelhantes aos observados?

A v16 já avaliou cinco exemplos com 199 replicações, com valores-p arquivados de 0,675; 0,625; 0,440; 0,915 e 0,820. Isso não responde pelos 87 pares nem demonstra ausência de relação física.

**Escopo proposto:** os 87 pares intrapadrão elegíveis, com ordenação determinística. Eventual ampliação aos 238 controles interpadrão será definida antes de observar novos valores-p, com família de testes correspondente. Não selecionar pares pelos maiores picos ou menores valores-p.

Preservar a condição da v16: cortar os sinais do par no menor comprimento disponível antes do preparo. Essa condição é própria da WCO e não coincide com os descritores fixos de 3.171 px de Pearson.

$$\operatorname{WCO}_{ij}(b,a)=\frac{|W_{ij}(b,a)|^2}{|W_{ii}(0,a)|\,|W_{jj}(0,a)|},\qquad T_{ij}=\max_{b,a\,\mathrm{válidos}}\operatorname{WCO}_{ij}(b,a).$$

O máximo percorre todos os lags e as 200 escalas. Manter o critério de denominadores numericamente válidos do núcleo v16, registrando escalas válidas de observados e simulados. Esse critério não é significância, COI ou máscara física. Não limitar artificialmente WCO a um.

**Modelo nulo:** randomizar independentemente fases da rFFT de cada sinal, preservando DC e Nyquist quando presente. Manter amplitudes de Fourier antes do novo detrending e reaplicar detrending e normalização como no diagnóstico arquivado. A preservação espectral não é exata depois dessa etapa e deverá ser descrita. O modelo não reproduz necessariamente costuras, morfologia não gaussiana ou não estacionariedade.

Não transferir esse nulo automaticamente para Pearson da WPS: a conservação do espectro restringe fortemente descritores baseados em autocorrelação e não constitui nulo adequado de ausência de semelhança espectral.

**Monte Carlo:** usar p=(1+b)/(1+B), salvar replicações e aplicar Holm à família completa de pares. O máximo já contempla a busca no plano lag–escala; Holm trata a busca entre pares.

Para 87 hipóteses e alfa=0,05, B=199 tem resolução insuficiente até para atingir o primeiro limiar de Holm: é necessário B+1≥87/0,05=1.740. Proposta: B=4.999 para 87 pares. Se a família incluir 325, propor B=9.999; B=4.999 não alcança o primeiro limiar nesse caso. Resolução suficiente não garante precisão perto do limiar: relatar incerteza Monte Carlo e não interpretar pequenas diferenças de p como diferenças físicas.

Antes de um lote autorizado, dimensionar custo em um par predefinido com poucas replicações. Esse piloto servirá somente a tempo e memória, não à seleção de pares ou a conclusões. Fixar o número final de replicações antes de inspecionar resultados inferenciais. Salvar as estatísticas necessárias sem manter milhares de mapas completos desnecessariamente.

**Validação:** convenção de lag, denominador, fases e casos pares/ímpares, reprodutibilidade por par, valores-p não nulos e Holm. Conferir os exemplos anteriores quando configuração, semente e ordem forem idênticas.

**Entregas futuras:** pares, replicações, configuração, sementes, tempo e metadados; máximos observados versus distribuição nula. Um máximo incomum não demonstra identidade de eventos ou acoplamento físico entre aquisições.

## 12. Mapa de arquivos para continuidade

Raiz do projeto: /Users/i521862/Documents/personal/unisinos/TCC. Os links abaixo são relativos a este documento.

| Arquivo | Uso |
| --- | --- |
| [Descritores atuais](../estatistica5/todos_contra_todos/results/descritores_por_ensaio.npz) | Fonte preferencial dos perfis e vetores de 26 ensaios |
| [Índice de ensaios](../estatistica5/todos_contra_todos/results/indice_ensaios.csv) | Identificadores, classes, ordem e comprimentos |
| [Pares atuais](../estatistica5/todos_contra_todos/results/pearson_325_pares.csv) | Correlações de 3.171 px e registros completos |
| [Matriz atual](../estatistica5/todos_contra_todos/results/matriz_pearson.csv) | A1 e conferência de referência |
| [Manifesto atual](../estatistica5/todos_contra_todos/results/manifest.json) | Hashes, versões e proveniência existente |
| [Inventário v16](../v16_spatial2_18_71/results/input_audit.csv) | Elegibilidade, metadados e hashes de entrada |
| [Configuração v16](../v16_spatial2_18_71/results/analysis_config.json) | Faixa de extração, grade, estimador e domínio |
| [Perfis completos](../v16_spatial2_18_71/results/profiles.csv) | Filtrar método e elegibilidade explicitamente |
| [Perfis de janelas](../v16_spatial2_18_71/results/equal_length_windows.csv) | Prefixo, centro e sufixo de 3.171 px |
| [Inventário enriquecido](../estatistica4/results/inventory.csv) | Inclui m_code e q_code |
| [Descritores de 7.001 px](../estatistica4/results/descriptors_principal_7001_inicio.npz) | Controle dos mesmos 19 ensaios |
| [Características existentes](../estatistica4/results/features_primary.csv) | Pico, inclinação, RMS residual e covariáveis de 7.001 px |
| [Coordenadas de janelas](../estatistica4/results/window_coordinates.csv) | Início e fim exclusivo das condições anteriores |
| [Associações existentes](../estatistica4/results/feature_covariate_correlations.csv) | Evitar repetição apresentada como análise nova |
| [Diagnóstico nulo](../v16_spatial2_18_71/results/null_diagnostic.csv) | Cinco exemplos de WCO |
| [Replicações nulas](../v16_spatial2_18_71/results/null_replicates.csv) | Conferência dos valores-p anteriores |

**Distinção essencial:** a referência atual está em estatistica5/todos_contra_todos. A pasta estatistica5/results e o script estatistica5/run_pearson.py pertencem à análise histórica de 87 pares. O script histórico pode sobrescrever documentos da raiz da estatistica5; não executá-lo para atualizar a análise atual.

### 12.1 Estrutura conferida do NPZ atual

Abrir com numpy.load e allow_pickle=False. Alinhar sempre pelo identificador run, mesmo que a ordenação aparente coincida.

| Chave | Dimensão | Significado |
| --- | --- | --- |
| runs; patterns | (26,) cada | Identificadores e classes |
| wavelengths_px | (200,) | Grade de 50 a 3.000 px |
| n_common_px | Escalar | 3171 |
| n_original | (26,) | Comprimentos originais |
| profile; residual; unit_descriptor | (26,200) cada | Perfil, resíduo e vetor unitário do prefixo |
| coefficients | (26,2) | Intercepto e inclinação, nessa ordem |
| profile_full; residual_full; unit_descriptor_full | (26,200) cada | Condição de registro completo |
| coefficients_full | (26,2) | Intercepto e inclinação completos |

O NPZ de 7.001 px tem 19 ensaios e chaves diferentes: runs, labels, q, m, wavelengths_px, profiles, residuals, log_profiles e coefficients. Não supor esquemas idênticos.

### 12.2 Sinais e tabelas

Os sinais estão em ../v16_spatial2_18_71/intermediate/{run}_input.npz, chave raw. Também existem raw_20_59, center, upper, lower, raw_shifts, shifts, scores e seams. Não substituir raw por outra extração sem declarar nova condição.

O panorama correspondente está em src2/data_processed/{run}/__panorama_v11.png, a partir da raiz do projeto. A conferência de extração usa a média de image.convert('L') nas linhas do slice [18:72]. Os panoramas permanecem reconstruídos com os deslocamentos anteriores; não refazê-los nesta estatística.

Em profiles.csv, selecionar os ensaios elegíveis e method igual a v16_biased para a referência completa. O arquivo mistura métodos e registros excluídos. Em equal_length_windows.csv, as colunas relevantes são run, pattern, window, n, wavelength_px e zero_lag_wps; window assume prefix, center e suffix. Os campos old_preprocessed_window_mean e old_preprocessed_window_std são diagnósticos históricos, não descritores para A2.

Validar booleanos explicitamente; conversão genérica da string False pode produzir verdadeiro. Usar float64 e precisão de leitura adequada nas conferências.

### 12.3 Funções reutilizáveis e efeitos de escrita

| Código | Funções úteis | Cuidado |
| --- | --- | --- |
| [wavelet_core.py](../v16_spatial2_18_71/scripts/wavelet_core.py) | prepare, correlation, transform_correlation, pair | Corte antes de prepare; mapas completos podem exigir muita memória |
| [statistics_core.py](../estatistica4/statistics_core.py) | zero_lag_profile, residual_descriptor, unit_rows | Conferir retorno e tolerâncias; calcular só o perfil quando bastar |
| Mesmo arquivo | binary_arrangements, exact_upper_p, holm | Enumeração binária; não generalizar automaticamente a cinco classes |
| [run_all_pairs.py](../estatistica5/todos_contra_todos/run_all_pairs.py) | residuals, pearson_matrix | Contém caminhos globais de saída na estatistica5 |
| [null_diagnostic.py](../v16_spatial2_18_71/scripts/null_diagnostic.py) | scramble, statistic | main escreve na v16; não chamá-lo na nova análise |
| [precision_and_windows.py](../v16_spatial2_18_71/scripts/precision_and_windows.py) | Referência das coordenadas de janelas | main também escreve na origem |

Preferir funções numéricas sem escrita ou adaptadores explícitos na estatistica6. Não executar programas antigos inteiros para obter um cálculo parcial. As fontes devem permanecer somente leitura.

Versões registradas no manifesto anterior: NumPy 2.3.4; SciPy 1.18.0; pandas 2.3.3; Matplotlib 3.10.7; Pillow 12.0.0; PyWavelets 1.9.0. São proveniência histórica, não confirmação do ambiente da próxima execução. Registrar Python e versões efetivas antes de calcular; não atualizar dependências sem necessidade.

## 13. Sequência de implementação futura

**Esta seção é uma instrução para uma etapa posterior, não uma ordem de execução agora.**

1. Ler o escopo mais recente do usuário e selecionar os blocos solicitados. Conferir instruções locais e alterações existentes no workspace.
2. Registrar hashes das fontes lidas e conferir o manifesto anterior. Identificar legitimamente eventuais revisões; divergência de hash, isoladamente, não demonstra corrupção.
3. Registrar configuração e hash antes dos cálculos: participantes, condições, características, estatísticas, famílias de hipóteses, replicações e sementes. Manter histórico das mudanças.
4. Validar elegibilidade, identificadores, dimensões, grade, finitude, contagens e valores de referência. Fazer joins por run.
5. Implementar controles matemáticos e sintéticos capazes de detectar erros substantivos, separados dos resultados experimentais. A7 formal requer auditoria; A8 requer dimensionamento.
6. Executar somente os blocos solicitados, salvando números antes dos gráficos. Resultados desfavoráveis e condições não estimáveis devem permanecer registrados.
7. Gerar gráficos e relatório a partir dos arquivos numéricos salvos, sem números manuais. Conferir identidades, implementações independentes e multiplicidade aplicável.
8. Inspecionar figuras e PDFs, conferir links, verificar hashes das fontes ao terminar e apresentar conclusões proporcionais ao delineamento.

Organização sugerida para a etapa futura: subpastas scripts, configs, results, figures e validation, com arquivos identificados por a1, a2 etc. **Essas subpastas e scripts não foram criados nesta entrega.** Não há comando executável de análise pronto.

Semente base sugerida: 2026102406. Criar fluxos separados por bloco, condição e identificadores ordenados com numpy.random.SeedSequence e identificadores estáveis. Não usar hash() do Python para sementes, pois varia entre processos. Registrar número solicitado e concluído de replicações.

## 14. Critérios de apresentação e conclusão

- Títulos formais, unidades explícitas, texto em português e nomes dos padrões em inglês.
- Evitar vermelho predominante. Para Pearson, conservar roxo–creme–azul e limites fixos de −1 a 1.
- Mostrar observações por ensaio quando possível. Detalhes extensos de processamento ficarão no método, nas legendas ou em comentários Markdown, não no título.
- Distinguir mínimo–máximo, faixa leave-one-out, intervalo de sensibilidade, envelope nulo e intervalo de confiança.
- Informar condição, participantes, pares, estimando, unidade amostral e limitações em cada tabela relevante.
- Ausência de evidência não demonstra equivalência. Equivalência exigiria margem substantiva justificada e análise própria, fora deste plano.
- Nenhum resultado será convertido em causalidade, repetibilidade ou classificação sem evidência específica.
- Aquisições independentes repetidas em condições comparáveis continuam necessárias para conclusões mais fortes.

## 15. Pendências e ponto de parada

A seleção dos blocos ainda não foi solicitada. Permanecem condicionais: permutabilidade dentro de q; interpretação física dos códigos; eventual intervalo populacional para estatísticas de pares; algoritmo e restrição de permutações de PERMDISP; família e custo computacional de A8.

Essas pendências não impedem os componentes descritivos de uma execução que venha a ser solicitada. Não preencher lacunas metodológicas com resultados presumidos.

**Entrega atual: somente este plano.** Foram consultados documentos, metadados e código existentes para prepará-lo. Não foram produzidos novos resultados estatísticos, scripts de análise, valores-p, intervalos, modelos ou gráficos.
