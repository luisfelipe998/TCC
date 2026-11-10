# Plano e verificação dos métodos

Este plano foi registrado antes do cálculo dos novos resultados. Trata-se de uma análise exploratória de dados já estudados, não de um protocolo prospectivo ou de uma validação externa.

## Pergunta e unidade de análise

A pergunta principal é se a forma do perfil da alternative WPS no lag zero apresenta maior similaridade entre ensaios do mesmo padrão do que entre churn e wavy annular. A unidade é o ensaio: 10 churn e 9 wavy annular, os mesmos 19 identificadores da estatística original. Os 200 comprimentos de onda são coordenadas correlacionadas de um descritor; os 171 pares compartilham ensaios. Nenhum dos dois é tratado como amostra independente.

Os 26 ensaios elegíveis da v16 entram no inventário e em resumos descritivos complementares. Plug e slug têm três registros cada; dispersed bubbles tem um único registro e não permite estimar dispersão intraclasse. Os cinco registros excluídos na v16 permanecem excluídos. Os códigos m e q são extraídos dos identificadores, sem atribuir unidades físicas não confirmadas.

## Descritor principal

1. Usar o sinal óptico calculado pela média das linhas 18–71, inclusive, dos PNGs originais. Verificar o vetor salvo diretamente contra os pixels e os hashes.
2. Selecionar os primeiros 7.001 pixels, o menor comprimento entre os 19 ensaios. Remover tendência linear e normalizar dentro desse trecho, seguindo a ordem da v16. Não reamostrar.
3. Usar a autocorrelação bilateral enviesada: Rxx[m] = soma x[n+m] x[n] / N, nos índices válidos. O divisor N é a adaptação declarada na v16; não se confunde com o divisor N−|m| das equações literais do artigo.
4. Usar Morlet de energia unitária, frequência angular central 6, e a conversão lambda = 4 pi a / (6 + sqrt(38)). Manter os 200 comprimentos de onda logaritmicamente espaçados de 50 a 3.000 px.
5. O perfil é P(a) = |CWT(Rxx)(0,a)|, sem elevar novamente ao quadrado. Calcular o coeficiente em zero pela soma temporal da wavelet conjugada e conferir contra a implementação FFT da v16.
6. Calcular w(a) = log(P(a) + 10^(-12)). Ajustar, separadamente em cada ensaio, w = beta0 + beta1 log(a) + resíduo. O descritor z é esse resíduo. Log(lambda) produz o mesmo resíduo porque lambda é proporcional a a. Esse procedimento retira uma tendência log-linear; não demonstra branqueamento estatístico, eliminação de toda dependência entre escalas ou remoção de um piso físico 1/f.

O pré-processamento retira informação de intensidade absoluta; a remoção da inclinação retira parte da forma espectral. Portanto, a pergunta se refere especificamente ao descritor residual. Não se selecionará faixa de escala ou transformação com base no resultado mais favorável.

## Pearson e ANOSIM

Para dois descritores z_i e z_j, r_ij = soma[(z_i−média_i)(z_j−média_j)] / sqrt(soma[(z_i−média_i)^2] soma[(z_j−média_j)^2]). Usar D_ij = 1−r_ij, no intervalo [0,2]. Recusar vetores constantes ou não finitos. Não calcular o valor-p paramétrico de Pearson usando as 200 escalas como tamanho amostral.

ANOSIM usa os postos de todas as distâncias i<j, com posto médio nos empates. Para n ensaios e M=n(n−1)/2 pares, R = (posto médio inter − posto médio intra)/(M/2). A normalização é n(n−1)/4. Uma verificação independente calcula R = 2 AUC − 1, com meio crédito para empates entre distâncias inter e intra.

Permutar os rótulos de ensaio somente dentro dos níveis de q, preservando suas contagens. Enumerar os 200 arranjos distintos, incluindo o observado. O valor-p unilateral exato é a fração R_perm >= R_obs; a resolução é 1/200. A enumeração completa não recebe a correção de mais um usada em amostragem Monte Carlo. Valores iguais dentro de tolerância numérica de 10^(-12) contam como empates. A validade inferencial é condicional à permutabilidade dentro de q, que o delineamento observacional não garante.

Também informar o espaço de permutações dentro de m e de (m,q), para explicitar o limite de controle simultâneo. Os testes não demonstram causalidade ou repetibilidade de aquisições sob condições idênticas. ECDFs e medianas são descritivas, sem testes que tratem pares como independentes.

## Análise complementar de localização e dispersão

Normalizar cada descritor centrado para norma euclidiana um, obtendo u_i. Então ||u_i−u_j||^2 = 2(1−r_ij). A dissimilaridade 1−r não é, em geral, uma distância métrica; sua transformação sqrt(2D) é euclidiana. PCA e somas de quadrados usam diretamente os vetores u.

Aplicar PERMANOVA de um fator nessa geometria: SS_total = soma ||u_i−média_global||^2; SS_intra = soma_g soma_i_em_g ||u_i−média_g||^2; F = [(SS_total−SS_intra)/(G−1)] / [SS_intra/(n−G)]. Calcular R² = (SS_total−SS_intra)/SS_total e o valor-p pelas mesmas 200 permutações de rótulos dentro de q. O R² é a separação marginal do padrão, não a fração ajustada exclusivamente atribuível ao padrão após q. O teste é restrito por q, mas a estatística não faz residualização das covariáveis.

Os dois testes principais, ANOSIM e PERMANOVA, recebem correção de Holm conjunta, além dos valores-p brutos. Mesmo com essa correção permanecem exploratórios. Ambos podem responder a diferenças de dispersão em grupos desbalanceados. Por isso serão apresentadas as distâncias ao centroide por classe, com ajuste descritivo sqrt(n_g/(n_g−1)), sem chamar esse resumo de teste PERMDISP. Uma correlação alta entre estatísticas não torna localização e dispersão perguntas equivalentes.

## Robustez e complementos definidos antes do cálculo

- Repetir R nos trechos de começo, centro e fim de 7.001 px; eles não são réplicas independentes. Informar a sobreposição.
- Comparar comprimentos completos, linha 40, referência 20–59 e trechos de 3.171 px. A comparação 18–71 versus 20–59 mantém o mesmo motor e pré-processamento, isolando a faixa extraída.
- Comparar Pearson dos resíduos, Spearman dos resíduos e Pearson dos log-perfis sem remoção de inclinação. As sensibilidades não serão usadas para substituir o resultado principal.
- Retirar um ensaio por vez e recalcular R, sem interpretar a faixa resultante como intervalo de confiança.
- Representar PCA na geometria euclidiana de Pearson e informar a variância retida. O gráfico não constitui validação de classificador.
- Resumir picos, inclinações log-lineares e associações descritivas com m, q, comprimento original e deslocamento de reconstrução. Não atribuir p-valores a essas buscas nem interpretar correlação como causa.
- Mostrar os perfis de 3.171 px dos 26 ensaios, sem promover classes pequenas a conclusões inferenciais.
- Conferir a estatística anterior a partir dos descritores arquivados. Diferenças com a v8 não isolam a mudança da faixa, pois o motor e a ordem de pré-processamento também diferem.

A máscara antiga de costura não será reproduzida como evidência de robustez quando não atingir nenhuma escala. Uma frequência nominal de costura abaixo da banda não prova imunidade a emendas, modulação, harmônicos, vazamento espectral ou erros de reconstrução. Nenhuma nova máscara ou COI será aplicada.

## Fontes para as definições

- Clarke, K. R. (1993). Non-parametric multivariate analyses of changes in community structure. Australian Journal of Ecology, 18, 117–143. DOI: https://doi.org/10.1111/j.1442-9993.1993.tb00438.x (ANOSIM).
- Anderson, M. J. (2001). A new method for non-parametric multivariate analysis of variance. Austral Ecology, 26, 32–46. DOI: https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x (PERMANOVA).
- Holm, S. (1979). A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics, 6, 65–70 (correção de multiplicidade).
- Documentação SciPy de Pearson e postos: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.pearsonr.html e https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.rankdata.html.
- Equações wavelet e adaptações: [auditoria da v16](../v16_spatial2_18_71/EQUATION_AUDIT.md), [núcleo executado](../v16_spatial2_18_71/scripts/wavelet_core.py) e [PDF local](../../assets/schuck-bardo-wavelet-coherence.pdf).

As fórmulas são verificadas por identidades e exemplos numéricos independentes antes da execução experimental. Os controles sintéticos servem apenas à validação do código; nenhum deles entra nas tabelas de resultados experimentais.
