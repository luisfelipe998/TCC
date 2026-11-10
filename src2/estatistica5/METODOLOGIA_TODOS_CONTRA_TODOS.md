# Método: Pearson direto, todos contra todos

Este protocolo foi definido antes dos novos cálculos. A análise é descritiva e utiliza os 26 ensaios elegíveis da v16_spatial2_18_71: 10 churn, 9 wavy annular, 3 plug, 3 slug e 1 dispersed bubbles. Os cinco ensaios excluídos no inventário original permanecem excluídos. O conjunto completo tem 26 × 25 / 2 = 325 pares únicos, sendo 87 intrapadrão e 238 entre padrões. A diagonal não conta como par experimental.

## 1. Verificar a extração óptica

Para cada panorama original em tons de cinza, o sinal foi conferido contra:

$$s_i[n]=\frac{1}{54}\sum_{h=18}^{71} I_i[h,n].$$

Os índices de linha começam em zero e ambas as extremidades são inclusivas. Foram verificados os hashes dos PNGs e a identidade dos vetores salvos com a extração. Não foi refeita a reconstrução ou a estimação dos deslocamentos dos panoramas.

## 2. Fixar um trecho por ensaio

Na análise principal foi usado o prefixo comum:

$$N=\min_i N_i=3171\ \mathrm{px},\qquad n=0,\ldots,N-1.$$

Todos os registros possuem esses pixels. Não houve interpolação, reamostragem ou escolha de trecho pelo resultado. Cada ensaio recebe um único descritor, que permanece idêntico em todas as 25 comparações com outros ensaios.

Essa decisão substitui os trechos dependentes do parceiro, usados nos mapas dos 87 pares da versão anterior. Portanto, a revisão não consiste apenas em preencher células cinzas ou mudar a barra de cores; os descritores foram recalculados com comprimento comum.

## 3. Retirar tendência e padronizar dentro do trecho

Mantendo o pré-processamento da v16, primeiro se inverte a intensidade e se ajusta uma reta:

$$v_i[n]=255-s_i[n],$$

$$ (\widehat\alpha_i,\widehat\beta_i)=\arg\min_{\alpha,\beta}\sum_{n=0}^{N-1}(v_i[n]-\alpha-\beta n)^2,$$

$$e_i[n]=v_i[n]-\widehat\alpha_i-\widehat\beta_i n,\qquad x_i[n]=\frac{e_i[n]-\overline e_i}{\sqrt{N^{-1}\sum_{n=0}^{N-1}(e_i[n]-\overline e_i)^2}}.$$

A normalização usa o desvio-padrão populacional, como no núcleo v16. Sinais constantes ou numericamente degenerados são recusados. Essa transformação remove intensidade absoluta e tendência linear; a análise posterior não recupera essas informações.

## 4. Estimar a autocorrelação enviesada

Para lags não negativos:

$$\widehat R_{ii}[m]=\frac{1}{N}\sum_{n=0}^{N-1-m}x_i[n+m]x_i[n],\qquad 0\le m<N.$$

Para sinais reais, a extensão bilateral é \(\widehat R_{ii}[-m]=\widehat R_{ii}[m]\). O divisor N preserva a adaptação enviesada usada como principal na v16; não é o divisor N−|m| da alternativa não enviesada. Todos os lags existentes foram mantidos.

## 5. Obter a alternative WPS no lag zero

A Morlet e a relação entre escala a e comprimento de onda lambda são:

$$\psi(t)=\pi^{-1/4}\exp(-t^2/2)\exp(\mathrm{i}6t),$$

$$\lambda=\frac{4\pi a}{6+\sqrt{38}}.$$

Usaram-se 200 comprimentos de onda logaritmicamente espaçados entre 50 e 3.000 px, com amostragem de 1 px. O coeficiente e o perfil utilizados são:

$$C_i(b,a)=\frac{1}{\sqrt a}\sum_{m=-(N-1)}^{N-1}\widehat R_{ii}[m]\overline{\psi\!\left(\frac{m-b}{a}\right)},\qquad P_i(a)=|C_i(0,a)|.$$

Para autocorrelação real e par, a parte imaginária se cancela no lag zero. O programa calcula esse coeficiente pela soma temporal com o cosseno, e o confere contra a implementação FFT da v16 nas 200 escalas de todos os 26 ensaios. Não se eleva novamente P ao quadrado. Não se calculam WCO cruzada nem correlação dos sinais brutos.

Não foram introduzidas máscaras ou COI. A maior escala é próxima do comprimento do trecho, de modo que a interpretação física nessa região permanece limitada pelo registro finito. A correção numérica da transformada não elimina essa limitação.

## 6. Definir o descritor de forma

Conservando a definição usada na estatística anterior:

$$w_{ik}=\log(P_i(a_k)+10^{-12}),$$

$$ (\widehat c_i,\widehat d_i)=\arg\min_{c,d}\sum_{k=1}^{200}(w_{ik}-c-d\log\lambda_k)^2,$$

$$z_{ik}=w_{ik}-\widehat c_i-\widehat d_i\log\lambda_k.$$

Log(lambda) e log(a) produzem o mesmo resíduo, pois lambda é proporcional a a. O ajuste usa mínimos quadrados sem inverter explicitamente a matriz normal. Essa operação retira uma tendência log-linear; não garante branqueamento estatístico, independência entre escalas ou remoção de um fundo físico 1/f.

Usar Pearson “como é” significa usar seu coeficiente assinado, sem transformar sua escala. O objeto comparado continua sendo o descritor residual acima, não a intensidade bruta da imagem.

## 7. Calcular Pearson direto

$$r_{ij}=\frac{\sum_{k=1}^{200}(z_{ik}-\overline z_i)(z_{jk}-\overline z_j)}{\sqrt{\sum_{k=1}^{200}(z_{ik}-\overline z_i)^2}\sqrt{\sum_{k=1}^{200}(z_{jk}-\overline z_j)^2}}.$$

Os coeficientes foram mantidos no domínio [−1,1]. Não foram usados 1−r, (1−r)/2, valor absoluto, quadrado de r ou normalização pelo mínimo/máximo da amostra. O gráfico usa uma escala divergente fixa de −1 a 1, centrada em zero.

- r próximo de 1: formas semelhantes após centralização e escala positiva.
- r próximo de 0: pouca associação linear entre as coordenadas dos descritores.
- r próximo de −1: formas linearmente opostas.

Pearson não mede igualdade de amplitude, causalidade ou probabilidade de pertencer ao mesmo padrão. r = 0 não demonstra independência. As coordenadas espectrais são correlacionadas; nenhum valor-p paramétrico foi calculado tratando as 200 escalas como 200 aquisições independentes.

Se u_i é o descritor centrado dividido por sua norma, a matriz é R=UUᵀ. Essa identidade permite verificar simetria, diagonal unitária e semidefinição positiva, admitindo arredondamento numérico. Todas as 325 comparações foram também conferidas com scipy.stats.pearsonr, usando somente o coeficiente retornado.

## 8. Resumir sem alterar a unidade de observação

Foram exportados todos os coeficientes, medianas e intervalos observados dos pares dentro e entre padrões, além das medianas por combinação de padrões. A diagonal de identidade foi excluída desses resumos. Os pares compartilham ensaios; os resumos não recebem intervalos de confiança calculados como se fossem independentes.

A classe dispersed bubbles tem somente um ensaio. Suas 25 comparações com outras classes são calculáveis, mas sua similaridade intraclasse não pode ser estimada. Na matriz 26×26 sua diagonal continua valendo 1 por identidade; no resumo por classe o valor intraclasse fica não estimável.

## 9. Conferir a sensibilidade ao comprimento

Como controle descritivo, repetiu-se somente Pearson usando o comprimento completo de cada ensaio. Cada registro continua tendo um único descritor fixo nessa segunda condição. Não há corte dependente do parceiro. Os comprimentos, contudo, diferem entre registros, de modo que esse controle não substitui a condição principal de 3.171 px.

O programa apresenta as duas condições sem selecionar aquela que produz maior correlação ou maior contraste entre padrões. Mudanças entre elas refletem conjuntamente comprimento e conteúdo espacial incluído. Nenhum teste de separação ou classificador foi acrescentado.

## 10. Verificações e rastreabilidade

Os controles verificam casos conhecidos r=−1,0,1, invariância a deslocamento/escala positiva, recusa de vetores constantes, equivalência com NumPy e SciPy, extração de pixels, ordem do pré-processamento, autocorrelação por soma literal, conversão escala–comprimento de onda, soma temporal versus FFT, ortogonalidade dos resíduos, cobertura dos 325 pares, ausência de células não calculadas, simetria e autovalores da matriz.

Os perfis calculados são ainda cruzados com as tabelas já existentes da v16: equal_length_windows.csv para o prefixo de 3.171 px e profiles.csv para o comprimento completo. Os arquivos de origem são somente lidos e seus hashes registrados antes e conferidos após a execução. Os resultados, gráficos e documentos também recebem hashes. Os controles sintéticos não entram nas tabelas experimentais.

Referências de definição: [núcleo v16](../v16_spatial2_18_71/scripts/wavelet_core.py), [auditoria das equações](../v16_spatial2_18_71/EQUATION_AUDIT.md), [documentação de Pearson no SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.pearsonr.html).
