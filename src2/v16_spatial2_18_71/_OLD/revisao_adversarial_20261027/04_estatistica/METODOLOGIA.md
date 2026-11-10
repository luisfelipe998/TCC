# Metodologia da análise espacial integral

## Aquisições e extração

Foram auditados 31 panoramas v11 e as respectivas planilhas originais. Mantiveram-se os 26 ensaios elegíveis da v16: 10 de churn, 9 de wavy annular, 3 de plug, 3 de slug e 1 de dispersed bubble. Foram calculados os 325 pares únicos, dos quais 87 foram intrapadrão e 238 entre padrões. Os cinco registros anteriormente excluídos permaneceram na auditoria, sem ingresso nas estatísticas.

O sinal foi extraído pela média aritmética das linhas 18–71, inclusive, totalizando 54 linhas. Os pixels foram novamente lidos dos PNGs originais e conferidos com os sinais arquivados. A montagem dos panoramas não foi refeita. As planilhas originais representaram a extração anterior 20–59 e serviram à conferência de procedência, sem substituir o sinal 18–71.

Cada panorama completo originou um perfil espectral individual fixo. Para os mapas de um par, empregaram-se os primeiros N = min(Na, Nb) pixels de ambos os sinais. O corte precedeu a inversão 255 − intensidade, a remoção da tendência linear, a centralização e a divisão pelo desvio padrão populacional. Não foram utilizadas âncoras, alinhamento de eventos ou escolha de trecho em função do resultado.

## Correlações e transformadas

Utilizaram-se a amostragem espacial de 1 px, uma realização, Morlet com ω₀ = 6 e 200 comprimentos de onda igualmente espaçados em escala logarítmica entre 50 e 3.000 px. A conversão foi λ = 4πa/(6 + √38). A faixa nativa foi preservada inclusive nos registros mais curtos; não foi substituída pela regra N/2 da análise de âncoras.

As correlações bilaterais usaram o estimador enviesado:

$$R_{xy}[m]=\frac{1}{N}\sum_{n\;\mathrm{válidos}}x[n+m]\,\overline{y[n]}.$$

O divisor N constituiu a adaptação já adotada pela v16. As definições alternativas utilizadas foram:

$$WPS_x(m,a)=|\mathrm{CWT}\{R_{xx}\}(m,a)|,$$

$$WCS_{xy}(m,a)=\mathrm{CWT}\{R_{xy}\}(m,a),$$

$$WCO_{xy}(m,a)=\frac{|WCS_{xy}(m,a)|^2}{WPS_x(0,a)\,WPS_y(0,a)}.$$

A WPS foi o módulo da transformada da autocorrelação; esse módulo não foi novamente elevado ao quadrado. A WCS foi armazenada complexa e apresentada em módulo. O denominador da WCO foi avaliado no lag zero. O piso numérico de 10⁻¹⁰ dos máximos individuais protegeu o denominador; todas as 200 escalas foram válidas em todos os pares. Os valores de WCO não foram truncados artificialmente em 1.

Foram utilizados todos os lags de −(N−1) a +(N−1), sem COI ou máscara geométrica. A potência média ao lado de cada WPS correspondeu à média aritmética do módulo em todos os lags. Os eixos dos sinais representaram a extensão do panorama; os eixos dos mapas, o lag da correlação. Essas coordenadas não foram tratadas como equivalentes.

## Resumo da coerência

Preservou-se o escore principal da v16. Para cada escala, calculou-se o máximo da WCO em todos os lags. Esses máximos foram ponderados por √(WPSa(0) WPSb(0)), retendo as escalas em que ambas as WPS atingiram pelo menos 5% dos respectivos máximos. As sensibilidades aos pisos de 1% e 10%, o máximo integral e a WCO ponderada no lag zero também foram registrados.

O piso de potência limitou somente as escalas empregadas no escore, sem remover regiões dos mapas. O máximo de cada escala pôde ocorrer em um lag diferente; o escore não representou um único alinhamento entre os sinais. A média integral da WCO foi fornecida como medida adicional e não substituiu o escore principal.

## Pearson dos perfis espectrais

Empregou-se a definição estatística da v16.3 sobre a entrada e a faixa espectral nativas da v16. Para cada ensaio, extraiu-se a WPS no lag zero do panorama completo e ajustou-se:

$$z_{ik}=\log(P_{ik}+10^{-12})-\widehat c_i-\widehat d_i\log(\lambda_k),\qquad P_{ik}=WPS_i(0,a_k).$$

O ajuste de c e d foi feito por mínimos quadrados nas mesmas 200 escalas. Calculou-se Pearson entre os resíduos centralizados, preservando o sinal da correlação. Cada descritor ficou fixo para todos os parceiros do ensaio. Não foi usado Pearson diretamente das intensidades, da WCO ou da curva de potência média.

A matriz utilizou 26 ensaios e λ de 50–3.000 px. Esse domínio diferiu da comparação anterior entre versões, que utilizou 25 ensaios comuns e λ de 50–1.024 px; os valores não foram comparados diretamente como efeito isolado de âncoras. A diagonal foi unitária e foi excluída das médias. Cada média geral por ensaio reuniu seus 25 parceiros.

## Agregações e alcance

As médias por par atribuíram o mesmo peso a cada par. As médias balanceadas atribuíram peso igual aos quatro padrões com pares intrapadrão e às dez combinações distintas de padrões. Dispersed bubble, com um ensaio, entrou nos 25 controles entre padrões que o envolveram, mas não teve repetibilidade intrapadrão estimada.

Os pares compartilharam aquisições e as escalas foram correlacionadas. Os resumos foram descritivos; não foram produzidos valores-p ou intervalos de confiança tratando pares ou escalas como réplicas independentes. Não houve treinamento de classificador nem estimativa de acurácia nesta análise. As unidades em pixels não foram convertidas em milímetros ou frequência temporal.

## Apresentação e verificação

As métricas e curvas médias foram calculadas com todos os lags, antes de qualquer redução para apresentação. Os arquivos de mapas em `intermediario/pares` conservaram aproximadamente 1.800 posições de lag, incluindo zero e os extremos, os coeficientes complexos, os sinais preparados e os resumos integrais. Esses mapas reduzidos serviram à reprodução das figuras; os máximos integrais não devem ser recalculados apenas deles.

As figuras foram dimensionadas para inserção a 16 cm de largura, com fonte mínima de 8 pt e altura máxima de 24 cm. Foram usadas as mesmas escalas de cor para as duas WPS de cada par e o intervalo 0–1 para exibição da WCO. Os nomes de arquivo identificaram os pares, sem classificação automática como favoráveis ou desfavoráveis.

Os 26 perfis, as métricas dos 325 pares e os mapas intrapadrão foram confrontados com os dados arquivados. Também foram verificadas a correlação por somas diretas, a WPS no lag zero por integração temporal independente da Morlet, a definição de Pearson contra a v16.3, todos os 325 coeficientes por SciPy, as matrizes por NumPy e a integridade dos arquivos arquivados.

Referências de implementação: [núcleo matemático preservado](../scripts/wavelet_core.py), [auditoria original das equações](../_OLD/EQUATION_AUDIT.md), [auditoria dos cálculos atuais](../intermediario/auditoria/calculo.json) e [auditoria estatística](../intermediario/auditoria/estatistica.json).
