# Auditoria e análise espacial por wavelets — v16, linhas 18–71

## Objetivo e decisões confirmadas

Esta pasta repete a análise de [v16_spatial2](../v16_spatial2/README.md), alterando a extração principal para a média das **linhas 18–71, inclusive: 54 linhas, indexadas a partir de zero, equivalentes a panorama[18:72, :]**. O objetivo é avaliar o efeito de incluir uma faixa maior do interior visível do tubo, mantendo a implementação de **Schuck Jr. e Bodmann (2019)** e os mesmos critérios de análise. A auditoria matemática herdada do v16 permanece documentada abaixo; esta extensão não introduz outra correção das equações.

**Estado atual:** as figuras e as estatísticas usam o plano completo de lags e a faixa de comprimentos de onda de 50–3000 px, sem COI e sem máscara conservadora de sobreposição, como a versão 20–59 imediatamente anterior. Os gráficos apresentam os nomes dos padrões em inglês e os demais textos em PT-BR. A mudança de faixa altera o sinal de entrada e exige recalcular as correlações, as transformadas, os mapas e os escores. A configuração está em [analysis_config.json](results/analysis_config.json).

**Nomenclatura:** as Eqs. 32.14–32.15 definem **alternative WPS**, a Eq. 32.16 define **alternative WCS** e a Eq. 32.17 define **alternative WCO**. Os símbolos WPS, WCS e WCO nas transcrições matemáticas abaixo conservam a grafia do artigo, mas, quando se referem a essas equações ou aos resultados do v5/v16, designam as versões alternativas. As definições convencionais são identificadas explicitamente como convencionais.

Foram auditados os 31 panoramas experimentais. As comparações por padrão usam os 26 registros confirmados pelo usuário: dispersed bubble (1), plug (3), slug (3), churn (10) e wavy annular (9). Permanecem fora dessas comparações: m03q01_001, m03q80_001, m04q05_001, m05q05_002 e m07q10_001. Não foram excluídos outros registros por apresentarem espectros inconvenientes. Os 87 pares intrapadrão resultam de 3 + 3 + 45 + 36; dispersed bubble não possui outro registro elegível para um par independente.

Os panoramas originais permanecem intactos. A entrada principal passa a ser a média 18–71; a média 20–59, idêntica à planilha original, é preservada somente como referência. A inspeção das imagens motivou testar estruturas visíveis abaixo da linha 59, especialmente em plug e slug. Ela não determina, com calibração independente, que 18 e 71 sejam exatamente os limites físicos internos da parede: reflexos, refração e outras contribuições ópticas continuam possíveis. Cobrir mais linhas é uma hipótese de melhor representação do escoamento, não garantia de aumento da alternative WCO.

Os controles da linha **40** e da faixa **30–49** foram mantidos fixos para isolar o efeito da expansão da entrada principal. As metades de diagnóstico passaram a **18–44 e 45–71**. A mediana usa as 54 linhas; a média aparada solicitada de 20% por extremo remove **10 linhas de cada extremo da ordenação de intensidades por coluna**, mantendo **34** (fração efetiva de 10/54 por extremo). Não se removem automaticamente dez linhas geométricas superiores e inferiores.

Consulte [RESULTS.md](RESULTS.md) para os resultados numéricos consolidados e [README.md](README.md) para execução e navegação.

A investigação mantém dez pré-processamentos na triagem de 64 escalas, três variantes na confirmação de 200 escalas e 325 pares por variante, incluindo controles entre padrões. Nesta pasta, a implementação v5 também recebe a nova entrada 18–71; isso deve ser distinguido do v5 original aplicado a 20–59. A comparação histórica entre faixas usa pesos e pisos de alternative WPS fixados na referência corrigida 20–59, disponíveis em [band_comparison.json](results/band_comparison.json) e [band_comparison_pairs.csv](results/band_comparison_pairs.csv). Usar diretamente escores nativos de faixas diferentes confundiria o efeito do novo sinal com a mudança dos pesos. A conclusão de desempenho e os números atuais estão em [EXECUTIVE_SUMMARY.md](EXECUTIVE_SUMMARY.md) e [PREPROCESSING.md](PREPROCESSING.md); ganhos, perdas e valores-p da análise 20–59 não são resultados desta extração.

## 1. Fontes e o que cada uma realmente estabelece

| Fonte primária | Evidência usada | Consequência para esta análise |
|---|---|---|
| Schuck Jr. e Bodmann, capítulo 32, 2019, pp. 409–423, DOI [10.1007/978-3-030-16077-7_32](https://doi.org/10.1007/978-3-030-16077-7_32) | Eqs. 32.4–32.8: CWT/Morlet/conversão de escala; 32.14–32.23: espectros de correlações e estimadores; 32.24–32.26: médias entre realizações | Fonte central. As páginas escaneadas foram lidas por OCR e as equações 32.4–32.23 foram verificadas visualmente. |
| de Oliveira et al., IJMF 144 (2021), 103793, DOI [10.1016/j.ijmultiphaseflow.2021.103793](https://doi.org/10.1016/j.ijmultiphaseflow.2021.103793) | Seção 3.2, Eqs. 8–9; seção 2; resultados de frequência dos sensores | Aplica CWT diretamente aos sinais temporais, com algoritmo adaptado de Torrence–Compo. Não define a versão de WCO do v5 nem valida os panoramas costurados. |
| Torrence e Compo, BAMS 79 (1998), 61–78, DOI [10.1175/1520-0477(1998)079<0061:APGTWA>2.0.CO;2](https://doi.org/10.1175/1520-0477(1998)079%3C0061:APGTWA%3E2.0.CO;2) | Tabela 1; seções 3g–3h, 5a–5b e 6c–6d | Morlet, relação escala/período, cone de influência, espectro global convencional e problema da coerência unitária sem suavização. |
| Grinsted, Moore e Jevrejeva, NPG 11 (2004), 561–566, DOI [10.5194/npg-11-561-2004](https://doi.org/10.5194/npg-11-561-2004) | Eqs. 8–10 e seção 3.4 | A coerência convencional exige suavização e calibração estatística. É outra definição, não uma correção automática da Eq. 32.17 de Schuck. |
| Maraun e Kurths, NPG 11 (2004), 505–514, DOI [10.5194/npg-11-505-2004](https://doi.org/10.5194/npg-11-505-2004) | Seções 3–6 e conclusões | Normalização, dependência entre coeficientes, efeitos de borda, falsos picos e testes múltiplos. Um mapa visualmente organizado não demonstra associação física. |

Os quatro primeiros PDFs já estavam em assets. O artigo de Maraun–Kurths foi obtido no site da revista. Os textos extraídos e as páginas de Schuck estão em [literature](literature/). O OCR de fórmulas não deve ser tratado como transcrição matemática: as imagens são a autoridade.

De Oliveira analisa pressão diferencial e sensor óptico medidos simultaneamente a 1 kHz; a câmera é usada para identificar padrões. No projeto, Voltage_2 é apenas o nome da coluna que contém intensidade de imagem. Ela não é uma medição da tensão do sensor óptico do artigo. Frequências em Hz do artigo não podem ser comparadas diretamente com comprimentos de onda em pixels. Fazer essa conversão exigiria uma relação espaço–tempo validada; v16 não a inventa.

## 2. Rastreabilidade dos panoramas

Não existe uma pasta src2/pipelines neste checkout. Os scripts relevantes são [pipeline_v11.py](../pipeline_v11.py), [pipeline_v13.py](../pipeline_v13.py) e [pipeline_v9.py](../pipeline_v9.py). O v11 produz o panorama e a planilha usados pelo v5. O v13 extrai a linha 40 do mesmo panorama. O v9 produz sinais temporais e não é a fonte desta análise espacial.

A auditoria anterior, preservada em [v16_spatial2/results/input_audit.csv](../v16_spatial2/results/input_audit.csv), fez duas verificações separadas:

1. Recalculou os deslocamentos a partir dos quadros ImgA, refez as tiras centrais e a mistura nas emendas e comparou o resultado com o PNG existente. **Naquela auditoria, os 31 panoramas foram reproduzidos exatamente, pixel a pixel.** A implementação vetorizada do objetivo de Pearson foi comparada com a função original em três pares de quadros por registro, além da comparação integral do panorama.
2. Comparou Voltage_2 com a média exata de panorama[20:60, :] e verificou posição, comprimento e valores finitos. **As 31 planilhas coincidiram**, dentro do erro de representação da planilha. Essa planilha continua representando 20–59; não foi reinterpretada como sinal 18–71.

A auditoria desta pasta conferiu os hashes dos PNGs, das planilhas, dos caches e dos registros anteriores; **não remontou panoramas nem reestimou deslocamentos**. Reextraiu a média 18–71 diretamente dos PNGs e a verificou por soma inteira independente dividida por 54. Exportou 31 arquivos em [results/signals](results/signals/) com posição, intensidade média 18–71 e referência 20–59, além dos caches correspondentes. Os deslocamentos e as emendas são dados herdados, com procedência explícita em [input_audit_summary.json](results/input_audit_summary.json). A montagem original continua usando os perfis 20–59 para estimar deslocamentos; ampliar a faixa de extração não reavalia esse passo.

Isso confirma a procedência computacional. Não comprova que a montagem seja uma fotografia espacial instantânea do escoamento. A imagem reúne estruturas observadas em tempos diferentes e depende de deslocamentos estimados por intensidade.

### Limitações concretas da montagem

- O v11 busca o melhor deslocamento inteiro pela correlação de perfis médios e transforma qualquer resultado menor ou igual a zero em 1 px. Não existe uma rejeição por baixa confiança. Em m03q40_001, **60 de 500** deslocamentos foram assim alterados; em m04q01_001, **1 de 500**. Entre os registros excluídos, m03q01_001 teve 116 e m03q80_001 teve 500. São flags de qualidade, não novas exclusões do escopo.
- A descrição do sinal do deslocamento em _best_shift está invertida em relação à expressão implementada: a[s:] é comparado com b[:-s]. Uma translação de b para a direita produz deslocamento estimado negativo. O sinal do código e o movimento físico precisam ser diferenciados; v16 não reinterpreta os negativos como uma velocidade física validada.
- A última tira usa o deslocamento anterior por não haver um próximo quadro. A origem da largura final é, portanto, uma extrapolação explícita.
- A montagem concatena tiras de largura variável, sem registrar uma malha física independente. Com deslocamentos diferentes entre quadros, escolher cada largura como o deslocamento seguinte não é o mesmo que definir fronteiras a meio caminho entre posições consecutivas. Sobreposição/omissão local e deformação das estruturas continuam possíveis, mesmo quando o PNG é reproduzido exatamente.
- A mistura de até 2 pixels altera intensidades nas emendas e a conversão final para uint8 quantiza o resultado. A sensibilidade a substituir vizinhanças de emendas por interpolação está em seam_sensitivity.csv. É um experimento de sensibilidade, não uma imagem supostamente corrigida. Quando mais da metade dos pixels seria alterada, esse teste não é realizado.
- A inversão 255−I não estabelece, por si, uma relação monotônica entre intensidade e fase líquida/gasosa. Para as autocorrelações e para o módulo da WCS quando ambos os sinais são invertidos e centrados, a inversão de sinal cancela. Não é necessário atribuir a cada nível de cinza uma fase para aplicar as equações de wavelets.

## 3. Auditoria das equações do v5

Referência de código: [wavelet_v5.py](../v5_spatial/wavelet_v5.py). A tabela distingue erro, escolha legítima e diferença de definição.

| Parte | v5 | Avaliação e tratamento em v16 |
|---|---|---|
| Entrada e pré-processamento, linhas 67–83 | Inversão, detrend linear, média zero e z-score | Válido como transformação da intensidade. O motivo físico do detrend não foi comprovado; há sensibilidade com apenas remoção da média. V16 rejeita entradas degeneradas/não finitas. |
| Comprimento do par, linhas 154–163 | Pré-processa o registro inteiro e depois corta o maior | O segmento cortado deixa de ter necessariamente média zero e variância unitária. V16 escolhe o prefixo comum e então pré-processa; o número de pixels descartados é exportado. |
| Autocorrelação, linhas 165–180 | Divide por N | Difere das Eqs. 32.21–32.22, que dividem por N−|m|. É uma adaptação deliberada e defensável, detalhada abaixo; não é reprodução literal. |
| Correlação cruzada, linhas 171–186 | Soma x[n] y[n+m] | Sinal de lag oposto ao da Eq. 32.23, que usa x[n+m] conjugado(y[n]). Não altera os módulos globais por princípio, mas espelha a interpretação da defasagem. V16 usa a convenção do artigo. |
| Correlações bilaterais | 2N−1 pontos, zero em N−1 | Correto para sinais reais. A extensão simétrica da autocorrelação é necessária; fornece suporte bilateral, não novas observações independentes. |
| WPS, linhas 200–201 e 256/275 | CWT da autocorrelação; exibe módulo | Compatível com Eqs. 32.14–32.15 para K=1. **Não deve ser substituído por módulo ao quadrado** se o objetivo é Schuck. O nome da variável guarda coeficientes complexos, mas o módulo é usado na exibição e normalização. |
| WCS, linha 202 | CWT da correlação cruzada | Compatível com Eq. 32.16 para K=1, salvo convenção do lag, estimador e detalhes da wavelet. |
| WCO, linha 205 | Módulo da WCS ao quadrado dividido por WPSx(0,a) WPSy(0,a) | A estrutura da Eq. 32.17 está correta. O denominador deve usar lag zero em todas as colunas. V16 adiciona proteção numérica sem somar um epsilon arbitrário nem recortar valores. |
| Morlet | cmor2.0-0.9549 | Tem envoltória compatível e portadora próxima de 6, porém normalização e implementação discretizada diferem da Eq. 32.6. V16 usa diretamente a expressão de Fourier da Morlet. |
| Escala/eixo | central_frequency × comprimento de onda | Não é a Eq. 32.8. Na versão instalada do PyWavelets, central_frequency=0,9375. A relação do artigo leva a um fator 0,968478 entre comprimento efetivo e rótulo do v5. |
| Cone de influência, linhas 24–54 e 96–142 | Um corte horizontal por comprimento/tipo de padrão | O limite no centro não descreve as outras posições. A validade não pode ser estabelecida a partir de uma auto-WCO próxima de 1. A revisão atual do v16 não aplica COI nem máscara conservadora aos mapas ou às estatísticas. |
| Mean power, linhas 211–213 | Média sobre lags de módulo(CWT(Rxx)) | Pode ser mantida como perfil descritivo do WPS de Schuck. Não é o espectro global de potência da CWT direta e não constitui automaticamente uma decomposição da variância. |
| Exibição de WCO, linha 307 | Escala de cor fixa 0–1 | Pode ocultar violações numéricas ou do estimador. V16 exporta valores acima de 1; o mapa primário usa 0–1 somente após verificar a propriedade. A versão não enviesada usa escala sem esse corte. |

A citação “Torrence–Compo Eq. 19” nos comentários do COI do v5 não é a localização correta da regra de e-folding. A referência relevante é a Tabela 1 e a seção 3g. A Eq. 19 trata de intervalos de confiança.

As notas INSIGHTS do v5 também não são a autoridade atual para os rótulos. Há divergências em m05q05_002, no número de registros anulares e nas exclusões por espectro. V16 usa o inventário atual e as exclusões confirmadas nesta conversa.

## 4. Enviesado versus não enviesado: resposta à questão central

Para sinais reais de comprimento N e m≥0:

$$\widehat R_{xy,u}[m]=\frac{1}{N-m}\sum_{n=0}^{N-1-m}x[n+m]y[n],$$

$$\widehat R_{xy,b}[m]=\frac{1}{N}\sum_{n=0}^{N-1-m}x[n+m]y[n].$$

Para lags negativos troca-se a orientação dos sinais. Assim,

$$\widehat R_b[m]=\left(1-\frac{|m|}{N}\right)\widehat R_u[m].$$

A diferença é uma janela triangular de Bartlett no domínio dos lags. Ela atenua progressivamente os grandes lags. Aplicar a CWT depois dessa multiplicação não equivale a multiplicar o resultado final da CWT pela mesma janela: a operação mistura uma vizinhança cujo tamanho depende da escala.

“Não enviesado” pressupõe um modelo estacionário e a formulação apropriada para a média. Após estimar/remover média e tendência a partir do próprio registro, o rótulo não garante ausência de todo viés. Tampouco elimina a variância muito elevada quando N−|m| é pequeno. Em |m|=N−1 existe apenas um produto no numerador.

### Por que a escolha enviesada é especialmente útil aqui

Para um registro finito estendido por zeros, o espectro da autocorrelação enviesada é proporcional a |X(ω)|²/N e é não negativo. A Morlet usada por Schuck possui transformada gaussiana real e não negativa. Escrevendo esse peso como H_a(ω)≥0, o coeficiente de correlação transformada pode ser expresso, a fatores de convenção comuns, como

$$C_{xy}(b,a)=\int X(\omega)\overline{Y(\omega)}H_a(\omega)e^{i\omega b}\,d\omega.$$

Pela desigualdade de Cauchy–Schwarz,

$$|C_{xy}(b,a)|^2\le C_{xx}(0,a)\,C_{yy}(0,a).$$

Logo a WCO fica entre 0 e 1 quando os mesmos pesos não negativos e a mesma normalização são usados nos três termos. Isso vale para a formulação numérica adotada em v16; não é um teorema para qualquer wavelet arbitrária. Dividir cada lag por N−|m| pode destruir a positividade do espectro amostral. Nesse caso, o denominador continua sendo calculável, mas deixa de garantir a desigualdade. WCO muito maior que 1 não pode ser corrigida simplesmente trocando a escala de cor.

**Conclusão:** manter 1/N é justificável e é a recomendação para o panorama finito. Deve-se declarar que se trata de uma adaptação de Schuck e mostrar a versão 1/(N−|m|) como diagnóstico de sensibilidade. Copiar literalmente o estimador não enviesado não foi, neste conjunto, uma melhoria de confiabilidade.

### Autocoerência não é uma verificação de confiabilidade da borda

Para x=y e denominador não nulo, a Eq. 32.17 produz exatamente WCOxx(0,a)=1: numerador e denominador são a mesma quantidade. Isso ocorre até quando o coeficiente tem pouca potência ou sofre efeito de borda. Portanto, observar 0,98–1,00 no lag zero não valida uma escala. Fora do lag zero a autocoerência pode diminuir, oscilar e apresentar lóbulos; ela não deve ser forçada a 1 em todo o plano.

## 5. Implementação da Morlet e fidelidade numérica

V16 usa ω0=6 e

$$\psi(t)=\pi^{-1/4}e^{-t^2/2}e^{i6t},\qquad
\lambda(a)=\frac{4\pi a}{6+\sqrt{38}}.$$

Para a convenção de FFT/IFFT do SciPy, o multiplicador é

$$H_a(\omega)=\sqrt{2\pi a}\,\pi^{-1/4}
\exp\left[-\frac{(a\omega-6)^2}{2}\right].$$

Foi usada a aproximação Morlet não corrigida por admissibilidade, explicitamente discutida no capítulo, cujo valor DC é proporcional a exp(−18). Não foi inserido um corte abrupto de frequências nem uma correção DC que alterasse silenciosamente a equação. A amostragem é 1 px por amostra. O preenchimento com zeros separa imagens periódicas da FFT por pelo menos 16 escalas máximas; os valores fora da região observada não são contados como dados.

Foram verificados: convolução temporal independente com a expressão analítica; aumento do padding de 8 para 12 sigmas por lado; frequência de uma senoide conhecida; recuperação do sinal e magnitude de um atraso conhecido; simetria após troca de x e y; invariância da WCO à escala e inversão de amplitude; auto-WCO em zero; equivalência FFT/convolução do PyWavelets. Os erros e critérios estão em numerical_validation.json.

Para B=2, a amplitude da função cmor do PyWavelets é 1/(√(2π)), enquanto a Morlet do artigo usa π^(−1/4). A razão é **0,531126**. Esse fator muda WPS, WCS e média de magnitude; cancela na WCO. A comparação em escalas idênticas corrige esse fator antes de medir o erro, evitando atribuir uma diferença de normalização a um defeito do método.

Os testes de precisão 10, 12 e 16 mostram que aumentar a resolução interna do PyWavelets ajuda, mas não assegura a propriedade de positividade do filtro discretizado em todo o plano. Um erro global pequeno nos coeficientes pode ser grande relativamente a uma faixa de potência fraca e ser amplificado pela normalização. V16 resolve a expressão espectral com pesos não negativos, em vez de limitar a posteriori WCO a 1.

Os resultados “v5” desta pasta são **recomputações da implementação v5 sobre a extração 18–71**, com a biblioteca e a precisão registradas no manifesto. O comparador histórico 20–59 é identificado separadamente na comparação entre faixas. Não há manifesto do ambiente que gerou todos os PNGs antigos de v5; não se afirma identidade binária com esses arquivos históricos. Na comparação de pares, o eixo do baseline é orientado pela convenção do artigo, para comparar os métodos; o sentido do eixo antigo é oposto. O perfil de potência e a magnitude não devem ser confundidos com a fase complexa sob essa troca.

Há uma discrepância no ambiente instalado: pywt.__version__ informa 1.8.0, enquanto os metadados da distribuição PyWavelets informam 1.9.0. Ambos são registrados no manifesto, junto com hashes do código da CWT. requirements.txt reproduz a versão declarada pela distribuição, mas isso, isoladamente, não comprova identidade com os binários instalados. A CWT primária de v16 usa NumPy/SciPy e a expressão explícita da Morlet, sem depender dessa implementação PyWavelets.

## 6. Domínio completo, bordas e média dos espectros

**As janelas adicionais não são necessárias para executar as equações de Schuck.** O resultado principal usou o registro disponível inteiro, ou o prefixo comum para um par de comprimentos diferentes. Os segmentos de começo/centro/fim permaneceram como controles opcionais de sensibilidade ao comprimento e à posição. Nenhuma janela de Hann ou outra janela adicional foi multiplicada nos sinais principais. A ponderação triangular do estimador enviesado continuou sendo uma consequência de dividir por N.

Nesta pasta, os mapas foram apresentados sobre todos os lags b em [−(N−1), N−1] e os comprimentos de onda de 50–3000 px. Não houve região cinza de COI, linha de sobreposição conservadora ou exclusão de coeficientes por essas regras. O domínio de resumo é o mesmo da análise 20–59 imediatamente anterior.

O COI adaptado da Morlet e a regra conservadora de sobreposição usados em uma revisão mais antiga eram escolhas adicionais, não partes das Eqs. 32.14–32.17. A retirada dessas regras já havia ocorrido em v16_spatial2. A alteração desta pasta é a extração 18–71: ela pode modificar os coeficientes, preservando a definição de alternative WPS, a alternative WCS complexa e o denominador de lag zero da alternative WCO. O histórico de apresentação permanece na [pasta anterior](../v16_spatial2/README.md).

Os perfis atuais são:

- alternative WPS(0,a), o valor no lag zero que participa do denominador da alternative WCO;
- média da alternative WPS(b,a) em todos os 2N−1 lags, equivalente ao perfil descritivo do v5.

Essa média em lag é uma média aritmética de magnitudes, não de magnitudes ao quadrado. Não é uma fração de energia física nem uma estimativa intercambiável com o perfil de lag zero. Médias sobre domínios diferentes também não são a mesma estatística: substituir uma média restrita pela média integral pode elevar ou reduzir seu valor. Um máximo no domínio ampliado, por sua vez, não pode diminuir se os coeficientes e as escalas comparadas forem os mesmos; seu aumento isolado não comprova maior similaridade.

A retirada das máscaras não eliminou o efeito do tratamento das bordas da CWT nem a redução do número de produtos na correlação em grandes lags. A extensão bilateral para 2N−1 amostras não criou novas observações independentes. Essas limitações permaneceram na interpretação, sem excluir posições dos mapas. A proteção relativa de 10^(−10) nos denominadores permaneceu apenas para evitar divisões numericamente degeneradas; ela não é uma máscara geométrica nem um teste de significância.

Comparar registros de comprimentos diferentes também altera a ponderação triangular e o intervalo de média. Por isso foi mantido o controle com janelas de comprimento comum de 3171 px, tomadas do começo, centro e fim de cada registro elegível. O eixo não foi esticado para igualar comprimentos.

## 7. A inconsistência interna de Schuck para múltiplas realizações

Na [página 413](literature/schuck-page-05.png), Eq. 32.16 define E{CWT(Rxy)}, preservando os coeficientes complexos. Na [página 414](literature/schuck-page-06.png), Eq. 32.20 apresenta a média de |CWT(Rxy)|. As operações são diferentes quando K>1:

$$\left|\frac1K\sum_k C_k\right|\ne\frac1K\sum_k|C_k|.$$

Também não se pode trocar, em geral, média de módulos de CWTs por módulo da CWT da correlação média nas Eqs. 32.18–32.19. A CWT é linear; o módulo não é. A descrição posterior do procedimento usando correlações médias não elimina essa distinção.

Para **K=1**, que é o caso de cada par do v5 e desta análise principal, essa inconsistência não muda a magnitude da WCS nem a WCO exibida. Portanto ela não explica, sozinha, os problemas encontrados no v5. Para um futuro agrupamento de realizações, a escolha precisa ser declarada. O controle sintético com três realizações em ensemble_definition_control.png implementa as duas versões, recupera −500 amostras de lag e quantifica a diferença na ordem das operações. É inspirado na experiência do capítulo; não é apresentado como reprodução exata de um código MATLAB que não foi fornecido.

Linhas diferentes do mesmo panorama não são automaticamente realizações independentes. Usá-las como réplicas para estreitar intervalos de confiança criaria pseudorreplicação. Esta pasta mede a sensibilidade à linha fixa 40 e às metades 18–44/45–71, mantendo a média 18–71 como entrada principal. A faixa fixa 30–49 permanece no experimento de pré-processamento.

## 8. Interpretação de WCS/WCO entre aquisições independentes

O exemplo de Schuck possui eventos relacionados e um atraso conhecido entre canais/realizações. Duas aquisições do mesmo regime de escoamento podem compartilhar uma distribuição de escalas sem conter os mesmos eventos ou uma fase espacial comum. Não há razão matemática para exigir WCO≈1 entre todas elas. Os mapas de lag não identificam a localização absoluta de uma bolha, e “x deslocado à direita de y” descreve a convenção matemática, não uma conclusão de causalidade.

Apenas comparar o máximo da alternative WCO pode premiar escalas com poucas oscilações, pontos de baixo denominador e uma busca sobre muitos lags. Nesta revisão, os máximos foram definidos no plano completo de lags e escalas, considerando valores numericamente válidos. A proteção relativa de 10^(−10) no denominador é estritamente numérica; não é um limiar de potência significativa. No experimento de pré-processamento, os pesos e os limiares de alternative WPS que definem os escores continuaram sendo escolhas descritivas explícitas, documentadas em [PREPROCESSING.md](PREPROCESSING.md), sem COI ou restrição conservadora de sobreposição.

Não se transferiram os limiares de Torrence–Compo/Grinsted para a alternative WCO de Schuck. São estatísticas diferentes. O diagnóstico adicional usa 199 pares de fases independentes aleatórias por exemplo, preservando inicialmente as amplitudes do espectro circular de cada sinal e repetindo a remoção de tendência. Sua estatística é o máximo da alternative WCO sobre todos os lags e as mesmas 200 escalas, tanto para os sinais observados quanto para cada par aleatório. Isso inclui a busca dentro do domínio completo na calibração. Os sinais observados e aleatórios desta pasta derivam de 18–71; os valores-p de 20–59 não são reaproveitados.

O modelo nulo não reproduz a morfologia não gaussiana, a evolução temporal nem os artefatos das emendas. Seus valores-p são **condicionais e exploratórios**, com resolução mínima 1/200, e não são usados para declarar relações físicas. Os cinco exemplos foram predeterminados: primeiro par em ordem alfabética de cada classe com réplicas e um controle slug/plug. Não foram selecionados pelo melhor resultado. Os valores-p entre esses cinco exemplos não foram usados para rejeições múltiplas. Uma validação inferencial definitiva exigiria justificar um modelo nulo para esses registros e dispor de réplicas experimentais adequadas.

## 9. Poderíamos aplicar wavelets de outra forma ao mesmo panorama?

| Alternativa | O que responde | Decisão para v16 |
|---|---|---|
| Correlação enviesada + Morlet analítica de Schuck | Associação em lag e escala, com controle da normalização finita | Implementação principal recomendada. Melhora fidelidade numérica e transparência. |
| Correlação não enviesada + mesma Morlet | Sensibilidade à reprodução literal de Eqs. 32.21–32.23 | Executada e exportada. Não recomendada como substituição automática: pode produzir WCO não limitada. |
| PyWavelets com maior precisão e escalas pela Eq. 32.8 | Reduz parte do erro numérico mantendo uma biblioteca conhecida | Testada. É uma alternativa possível, mas a propriedade de limitação ainda precisa ser verificada. |
| Média de linhas, linha central ou bandas verticais | Observáveis ópticos distintos da mesma imagem | Sensibilidade executada. Não há evidência suficiente para declarar uma linha universalmente superior. |
| Janelas de comprimento comum e posições distintas | Dependência do espectro em relação à duração espacial e à região selecionada | Executadas como controle, sem substituir silenciosamente o registro completo. |
| CWT direta da intensidade; módulo ao quadrado dos coeficientes | Localização de estruturas ao longo do panorama; mais próxima da aplicação de de Oliveira | Matematicamente apropriada para essa pergunta, mas **é outra grandeza**. Não substitui os mapas de correlação solicitados. |
| WCS convencional W_x conjugado(W_y) e WCO suavizada | Coerência local entre observáveis alinhados na mesma posição | Exige alinhamento justificável e suavização. Não fornece automaticamente uma medida de semelhança de regime entre aquisições sem correspondência espacial. |

Na CWT direta, a potência média convencional é média_b |W_x|². A integração de energia sobre escalas requer a medida e a normalização apropriadas; Torrence–Compo usa pesos inversos à escala nas médias por banda. Isso não autoriza aplicar mecanicamente o mesmo procedimento ao módulo de CWT(Rxx), que possui outro filtro espectral e outra interpretação. Também não se deve elevar o WPS de Schuck ao quadrado apenas porque a palavra “power” aparece em ambos os métodos.

A contribuição metodológica do v16 frente ao v5 foi uma implementação mais fiel da Morlet e da conversão de escala de Schuck, com limites numéricos verificados, perfis definidos e efeitos do estimador quantificados. Essa implementação foi preservada nesta pasta. **A correção matemática e a inclusão de mais linhas não demonstram, por si, maior alternative WCO entre aquisições do mesmo padrão nem melhor separação entre padrões.** A avaliação empírica, incluindo ganhos e perdas dos pré-processamentos e da mudança de faixa, está em [PREPROCESSING.md](PREPROCESSING.md) e [EXECUTIVE_SUMMARY.md](EXECUTIVE_SUMMARY.md). A comparação entre faixas usa o mesmo domínio de resumo e pesos de referência 20–59 fixos; a faixa 30–49, que não mudou, funciona como controle de consistência.

## 10. Limites que permanecem

Os panoramas ainda dependem da reconstrução original, há um registro elegível com muitos deslocamentos truncados, o contraste óptico não equivale a uma fração de fase, os registros não são réplicas sincronizadas e não há validação física independente do eixo espacial reconstruído. Retirar o COI e a máscara conservadora não resolveu essas limitações nem eliminou os efeitos do registro finito. A reprodução exata dos PNGs de Schuck exigiria os sinais, o código e as escolhas MATLAB originais, não disponíveis no PDF.

A recomendação operacional é usar a versão enviesada como análise principal, manter o diagnóstico não enviesado e apresentar alternative WPS, alternative WCS, alternative WCO e média da alternative WPS em lag sem máscaras geométricas. Os controles de entrada, comprimento e precisão numérica continuam relevantes para interpretar esses descritores. Um máximo no mapa completo não deve ser confundido com significância estatística, identidade de eventos ou validação física do panorama.
