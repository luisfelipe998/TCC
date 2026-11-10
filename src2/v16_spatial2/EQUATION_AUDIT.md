# Verificação das equações com o PDF de Schuck em assets

Fonte: [PDF original](../../assets/schuck-bardo-wavelet-coherence.pdf), capítulo 32, Schuck Jr. e Bodmann, 2019. As páginas renderizadas e conferidas visualmente são [411](literature/schuck-page-03.png), [413](literature/schuck-page-05.png) e [414](literature/schuck-page-06.png).

**Revisão atual:** foram retirados o COI e a máscara conservadora dos gráficos e das estatísticas. O domínio passou a ser integral em lag, mantendo 50–3000 px e as equações abaixo. A retirada foi uma escolha de apresentação e de resumo, não uma mudança da formulação de Schuck. As saídas anteriores com máscaras foram preservadas em [historico_com_mascaras.zip](results/historico_com_mascaras.zip).

## Definições alternativas preservadas

O nome **alternative** distingue estas grandezas das definições convencionais baseadas em produtos ponto a ponto das CWTs dos sinais.

**Alternative WPS — Eqs. 32.14 e 32.15:**

$$P_x^{alt}(b,a)=E\left\{\left|\mathrm{CWT}\{R_{xx}\}(b,a)\right|\right\},$$

$$P_y^{alt}(b,a)=E\left\{\left|\mathrm{CWT}\{R_{yy}\}(b,a)\right|\right\}.$$

**Alternative WCS — Eq. 32.16:**

$$C_{xy}^{alt}(b,a)=E\left\{\mathrm{CWT}\{R_{xy}\}(b,a)\right\}.$$

**Alternative WCO — Eq. 32.17:**

$$\Gamma_{xy}^{alt}(b,a)=
\frac{|C_{xy}^{alt}(b,a)|^2}
{|P_x^{alt}(0,a)|\,|P_y^{alt}(0,a)|}.$$

Na análise principal há uma realização por aquisição: **K=1**. O operador de esperança é estimado por essa realização, não por uma média entre pixels independentes. No código, wxx e wyy armazenam os coeficientes complexos intermediários; **abs(wxx)** e **abs(wyy)** são as alternative WPS. wxy é a alternative WCS. O módulo da alternative WCS é usado no gráfico; seu módulo ao quadrado é usado no numerador da alternative WCO.

O denominador contém o valor de cada alternative WPS no **lag zero**, índice N−1 do vetor bilateral. Ele é transmitido a todos os lags da mesma escala. Não se substituiu esse denominador por valores em b, pelo máximo sobre b, por uma potência convencional ou por um valor ajustado para elevar a coerência. Não houve suavização adicional do mapa de alternative WCO nem clipping dos resultados calculados. A proteção numérica relativa de 10^(−10), aplicada a cada alternative WPS no lag zero em relação ao seu próprio máximo entre escalas, permaneceu para evitar divisões numericamente degeneradas. Essa proteção não é COI, máscara de sobreposição ou limiar de significância.

Implementação: [wavelet_core.py](scripts/wavelet_core.py), funções pair, transform_correlation e safe_coherence. A verificação final recalcula a Eq. 32.17 diretamente dos coeficientes salvos nos 87 mapas.

## CWT, Morlet e escala

**Eq. 32.4:**

$$\mathrm{CWT}\{r\}(b,a)=\frac1{\sqrt a}
\int r(t)\overline{\psi\left(\frac{t-b}{a}\right)}dt.$$

**Eq. 32.6**, com ω0=6:

$$\psi(t)=\pi^{-1/4}e^{-t^2/2}e^{i\omega_0t}.$$

**Eq. 32.8:**

$$\lambda=\frac1f=\frac{4\pi a}{\omega_0+\sqrt{2+\omega_0^2}}.$$

V16 implementa a CWT por multiplicação no domínio de Fourier. O fator √(2πa) em sua implementação resulta da convenção FFT/IFFT, cuja transformada inversa contém 1/(2π) no contínuo. Não é uma mudança da normalização temporal acima: a convolução temporal independente com π^(−1/4)/√a foi comparada numericamente com a implementação espectral, com erro máximo da ordem de 10^(−14).

A própria fonte descreve a Morlet não corrigida por admissibilidade e aceita a aproximação para ω0>5,5. V16 conserva essa escolha. Não substitui a Morlet por uma wavelet de banda mais estreita escolhida para elevar a alternative WCO. As unidades de a e b são pixels na análise experimental, pois a amostragem é 1 px/amostra.

## Estimadores: diferença explícita em relação ao artigo

**Eq. 32.21:**

$$\widehat R_{xx,u}[m]=\frac1{N-m}\sum_{n=0}^{N-1-m}x[n+m]\overline{x[n]},\quad 0\le m\le N-1.$$

**Eq. 32.22:**

$$\widehat R_{yy,u}[m]=\frac1{N-m}\sum_{n=0}^{N-1-m}y[n+m]\overline{y[n]}.$$

**Eq. 32.23:**

$$\widehat R_{xy,u}[m]=\frac1{N-m}\sum_{n=0}^{N-1-m}x[n+m]\overline{y[n]}.$$

A mudança dos índices 1…N para 0…N−1 é apenas indexação Python. Para lags negativos utiliza-se a extensão cruzada/conjugada correspondente; o divisor passa a N−|m|.

O usuário confirmou que a metodologia vinha usando o estimador **enviesado**, com divisor **N**. V16 conserva essa versão como principal e calcula a versão **não enviesada**, com divisor **N−|m|**, como comparação literal dos estimadores do artigo. **Assim, não se afirma que o estimador principal reproduz literalmente Eqs. 32.21–32.23. As definições de alternative WPS/WCS/WCO são preservadas; essa adaptação do estimador é declarada.**

A função correlation usa a convenção x[n+m] conjugado(y[n]), como a Eq. 32.23. Um y deslocado 100 amostras à direita de x produz pico em lag −100. O teste recupera esse resultado. O controle inspirado no exemplo do artigo recupera −500.

## O que não pertence às equações citadas

- A média da alternative WPS sobre todos os 2N−1 lags é um resumo descritivo herdado do v5. É uma média de magnitudes, não de magnitudes ao quadrado. Não é apresentada como uma nova equação de potência física do capítulo.
- O COI adaptado da Morlet de Torrence–Compo e a máscara conservadora de sobreposição pertenciam à versão anterior desta análise. Foram retirados dos gráficos e das estatísticas atuais; não eram fatores das fórmulas alternativas.
- As janelas começo/centro/fim permaneceram como controles opcionais de sensibilidade ao comprimento e à posição. Não são fatores multiplicados nas fórmulas alternativas.
- Os escores resumidos usados no experimento de pré-processamento são resumos dos mapas de alternative WCO, agora sem exclusão por COI ou sobreposição. Seus pesos e limiares de alternative WPS são declarados em [PREPROCESSING.md](PREPROCESSING.md). Não substituem os mapas nem são atribuídos a Schuck.
- Os limiares estatísticos convencionais de Grinsted não foram aplicados à alternative WCO.

Ampliar o domínio de busca pode aumentar os máximos e alterar as médias sem mudar qualquer coeficiente wavelet. Por isso, os resultados atuais foram comparados entre métodos no mesmo domínio, e os números antigos com máscaras não foram apresentados como se tivessem sido calculados sem elas. A remoção das máscaras não demonstra desaparecimento dos efeitos de borda nem maior similaridade física entre aquisições.

## Ambiguidade que existe no próprio PDF

A Eq. 32.16 calcula a esperança da CWT **complexa**. A Eq. 32.20 imprime uma média de **módulos**. Para K>1 essas operações não são equivalentes. O PDF também descreve uma CWT de correlações médias, que não equivale em geral à média dos módulos nas Eqs. 32.18–32.19. A divergência é mostrada no relatório e em um controle sintético. Para K=1, ela não altera as magnitudes/WCO da análise principal. Nenhuma escolha entre essas operações foi usada para elevar artificialmente os resultados.

## Evidência numérica

Consulte [numerical_validation.json](results/numerical_validation.json), [legacy_literal_validation.json](results/legacy_literal_validation.json), [additional_validation.json](results/additional_validation.json) e [final_verification.json](results/final_verification.json). O manifesto registra os hashes do PDF original, dos scripts e da implementação PyWavelets usada como baseline.
