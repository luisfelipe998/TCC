# Identificação física e limites condicionais de velocidade

## Resultado da revisão

Há uma alternativa física útil para **bolhas alongadas isoladas em plug/slug**: a conservação de volume permite limitar a velocidade média do líquido no pistão e, quando uma célula completa é identificada, a média intrínseca do líquido nessa célula. A razão não exige metros por pixel. Entretanto, exige uma hipótese verificável apenas parcialmente nas imagens: **o filme líquido médio, em uma seção do corpo da bolha, escoa no mesmo sentido e não mais rápido que a bolha**. Os intervalos obtidos são condicionais a essa hipótese e à reconstrução da área de gás. Não constituem medição independente de velocidade líquida.

As duas referências consultadas não validam microbolhas como traçadores passivos do líquido, não identificam velocidade média do filme com celeridade de onda e não oferecem uma calibração óptica transferível para os 26 vídeos. A ausência de calibração métrica deixa de ser o obstáculo principal ao usar razões; a identidade do movimento e o fechamento físico permanecem.

Esta revisão foi feita sem preencher lacunas com propriedades, entalpias, escalas ou condições de entrada presumidas. Os controles executados são matemáticos, em unidades adimensionais, e estão separados dos dados experimentais.

## 1. O que os artigos realmente sustentam

**de Oliveira et al. (2021), International Journal of Multiphase Flow 144, 103793, DOI 10.1016/j.ijmultiphaseflow.2021.103793.**

- As páginas 3–4 do PDF descrevem o tubo horizontal, diâmetro interno de 1,0 mm, seção de visualização com o mesmo diâmetro, R600a e temperatura de saturação nominal de 27 °C. O pré-aquecedor estabelece entrada sub-resfriada na seção de testes.
- As páginas 7–8, equações (2)–(5), obtêm título a partir do balanço de energia e fração de vazio por correlação de Xu–Fang. As propriedades são provenientes do REFPROP. O artigo não fornece todas essas variáveis associadas a cada arquivo de câmera atual.
- A página 14, equação (17), define a razão entre velocidades médias de vapor e líquido por `S = x ρ_l (1−α) / [α ρ_g (1−x)]`. Essa relação depende de título, densidades e fração de vazio; ela não representa duas velocidades extraídas por rastreamento das imagens.

**Dias et al. (2025), International Journal of Multiphase Flow 185, 105122, DOI 10.1016/j.ijmultiphaseflow.2024.105122.**

- A página 4 informa 3000 quadros/s e 13,68 µm/pixel; a página 5 informa 200 imagens de 1016 × 88 pixels por vídeo. São condições diferentes dos 501 quadros de 1016 × 80 e 4000–7000 quadros/s desta investigação.
- A página 5, equações (6)–(7), define as velocidades **superficiais** `j_l = G(1−x)/ρ_l` e `j_g = Gx/ρ_g`. A tabela 2 deriva desse tratamento termodinâmico. Os termos “phase velocity” na legenda não alteram a definição dada nas equações.
- Na página 5, seção 2.4, o artigo discute variações de iluminação entre quadros. Isso documenta uma fonte de variação óptica, sem provar que ela explica os movimentos observados nestes arquivos.

Portanto, os artigos sustentam a necessidade de fechar a redução termodinâmica ou hidrodinâmica, mas não autorizam copiar seus valores para os vídeos atuais. Os arquivos originais e as extrações consultadas estão em [de Oliveira, PDF](../../../assets/deoliveira2021-flow-patterns-pressure-drop.pdf), [Dias, PDF](../../../assets/dias2025-knn-flow-pattern.pdf), [texto de de Oliveira](../../v15.1_gas_liquid_velocity/intermediate/calibration_literature/deoliveira2021.txt) e [texto de Dias](../../v15.1_gas_liquid_velocity/intermediate/calibration_literature/dias2025_knn.txt).

## 2. Bolha isolada: o que frente, traseira e centroide medem

Se `x_f` e `x_r` representam as extremidades de um mesmo objeto, então a velocidade do ponto médio é `(dx_f/dt + dx_r/dt)/2`, enquanto a variação do comprimento depende da diferença das velocidades das extremidades, com sinal determinado pela convenção adotada. Comparar essas velocidades identifica alongamento, encurtamento e erros de correspondência. Concordância entre elas é necessária para uma aproximação de translação rígida; não basta para provar volume constante, ausência de evaporação ou identidade material.

O centroide da máscara plana é ponderado pela **área projetada**. O centroide tridimensional da bolha seria ponderado por sua área transversal ao longo do eixo. Se, adicionalmente, cada seção de gás for circular e seu diâmetro observado `d(x)` corresponder ao diâmetro real com escala constante, pode-se usar peso proporcional a `d(x)²`. A comparação entre esses dois centroides é uma sensibilidade geométrica legítima. A reconstrução circular não se torna verdadeira porque ambas as velocidades concordam.

Para uma bolha material de densidade constante, sem mudança de fase, sem fragmentação/coalescência e com volume corretamente reconstruído, a derivada do centroide volumétrico corresponde à velocidade média volumétrica do gás contido. Com evaporação assimétrica, compressão ou variação de forma fora do plano, esse resultado não se segue do centroide projetado. Uma bolha fechada e quase rígida oferece interpretação física mais forte que uma textura anular, mas não mede a velocidade do líquido ao redor.

## 3. Limite útil para o líquido em uma célula plug/slug

### 3.1 Hipóteses e definição das grandezas

Considere uma célula que se translada à velocidade positiva `c`: uma bolha isolada e seu pistão adjacente de líquido. Se o movimento das imagens for para a esquerda, primeiro se inverte o eixo, preservando os sinais. As hipóteses são:

1. A forma e a distribuição de fração de gás se deslocam aproximadamente como `α(x−ct)`, sem variação relevante durante o intervalo analisado.
2. As densidades são aproximadamente constantes; evaporação/condensação, coalescência e fragmentação são desprezíveis no intervalo.
3. Existe uma seção do pistão verdadeiramente sem gás, inclusive fora do plano visível.
4. `α` representa a **fração de área transversal** ocupada por gás. Uma fração projetada da máscara não satisfaz automaticamente essa definição.
5. Em uma seção identificada do corpo da bolha, a velocidade média do filme satisfaz `0 ≤ U_f ≤ c`. Essa é uma hipótese de ramo de escoamento, não uma consequência da condição de não deslizamento na parede.

Definem-se `α_b` como a fração de gás nessa seção do corpo, `U_s` como a média do líquido sobre a seção integral do pistão e `ᾱ` como a fração volumétrica média de gás sobre **uma célula completa**. `U_l,cell` é a média do líquido ponderada por seu volume nessa mesma célula. `U_s` e `U_l,cell` são grandezas diferentes.

### 3.2 Derivação

No referencial da bolha, as equações de conservação sem fonte dão:

```text
α (U_g − c) = constante = 0,
(1−α) (U_l − c) = constante = U_s − c.
```

A primeira constante é nula porque existe uma seção sem gás; portanto, a média transversal do gás é `U_g=c` onde há gás, sob essas hipóteses. Avaliando a segunda expressão no corpo da bolha:

```text
U_s = α_b c + (1−α_b) U_f.
```

Logo, a hipótese `0 ≤ U_f ≤ c` fornece o intervalo:

```text
α_b ≤ U_s/c ≤ 1.
```

Esse intervalo limita o líquido **no pistão**. Ele não pode ser apresentado como média de todo o líquido do vídeo. Para a célula completa, o fluxo total é `U_s`, o fluxo médio de gás é `c ᾱ` e:

```text
U_l,cell/c = (U_s/c − ᾱ)/(1−ᾱ).
max[0, (α_b − ᾱ)/(1−ᾱ)] ≤ U_l,cell/c ≤ 1.
```

Usar uma seção de maior fração de gás fornece um limite inferior mais forte, desde que a hipótese sobre o filme seja válida nessa seção. Se `α_b` for o máximo verdadeiro do perfil, então `α_b ≥ ᾱ`. Usar apenas uma seção central não prova que ela seja o máximo.

Uma aproximação `ᾱ = fração do comprimento ocupada pela bolha × α_b` ignora as pontas da bolha; não deve substituir a integração da área de gás quando esta está disponível. Tampouco é correto usar a fração média de todo o vídeo no lugar da fração da célula rastreada.

### 3.3 Geometria incerta, comparações e casos degenerados

Se `a_L` for um limite inferior válido para `α_b` e `b_U<1` um limite superior válido para `ᾱ`, o limite conservador é:

```text
r_L = max[0, (a_L − b_U)/(1−b_U)],
r_L ≤ U_l,cell/c ≤ 1.
```

Para `r_L>0`, segue `1 ≤ U_g/U_l,cell ≤ 1/r_L`. Se `r_L=0`, a razão não possui limite superior finito sob essas hipóteses. A dispersão entre limiares de segmentação é uma **sensibilidade**, não um intervalo de confiança que garanta conter a geometria verdadeira.

Em duas células/gravações elegíveis, com razão física de velocidades de bolha `C=c_i/c_j` e limites líquidos `r_Li`, `r_Lj`, a conservação fornece:

```text
C r_Li ≤ U_l,i/U_l,j ≤ C/r_Lj.
```

O extremo superior é ilimitado se `r_Lj=0`. Para `C` conhecido apenas no intervalo `[C_L,C_U]`, usam-se `C_L r_Li` e `C_U/r_Lj`. A razão de velocidades em pixels/s entre gravações equivale à razão física somente com escala axial relativa conhecida ou com hipótese explícita de magnificação igual. Sem isso, a comparação continua em coordenadas de imagem.

Essas relações podem permitir uma ordenação condicional entre células quando todo o intervalo de uma razão estiver acima ou abaixo de 1. Intervalos sobrepostos não sustentam uma ordenação forçada. Um único vídeo por condição e poucas células não estabelecem uma distribuição populacional de cada padrão.

### 3.4 Como o limite pode falhar

- Filme médio em contracorrente, `U_f<0`, permite `U_s/c<α_b`. A parede imóvel não impede refluxo médio do filme.
- Filme médio mais rápido que a bolha, `U_f>c`, permite `U_s/c>1`. Não existe limite superior igual a 1 sem excluir esse ramo.
- Pontas que avançam de modo diferente, crescimento volumétrico, bolhas dentro do pistão e troca entre bolhas invalidam a célula rígida simples.
- Churn e anular não apresentam necessariamente uma bolha isolada conectada a um pistão sem gás. Aplicar a mesma fórmula para completar esses padrões seria inadequado.
- A ausência de gás na imagem bidimensional não demonstra ausência de gás em toda a seção. Películas abaixo da resolução tampouco podem ser fixadas em espessura zero.

Consequentemente, esta alternativa deve aparecer como **“intervalo de velocidade líquida condicionado ao modelo de célula”**, ao lado da velocidade de translação da bolha. Mesmo quando estreito, o intervalo não é uma medição experimental independente do líquido.

## 4. Microbolhas como traçadores naturais

O movimento de uma microbolha é um candidato a acompanhar o líquido, mas não é automaticamente um traçador passivo. Ser pequena na imagem não demonstra pequeno número de Stokes, pequeno deslizamento nem neutralidade de flutuação. Os fatores relevantes incluem tamanho tridimensional, viscosidade, posição radial/profundidade, mobilidade da interface, gradiente de pressão, aceleração local, confinamento e mudança de fase.

O modelo de Stokes para uma esfera com interface imobilizada, isolada em líquido quiescente e em regime de Reynolds pequeno, daria a escala de velocidade terminal de flutuação `w = (ρ_l−ρ_g) g d²/(18 μ_l)`. A hipótese de bolha limpa altera o coeficiente. Essa expressão não é um limite superior universal no escoamento presente. Em um tubo horizontal, a gravidade é transversal; sua velocidade terminal não pode ser apresentada como erro axial do rastreamento.

As imagens não fornecem todas as grandezas necessárias para limitar esse deslizamento. Além disso, a exposição registrada corresponde a 43,4–66,0% do intervalo entre quadros; o deslocamento integrado durante a exposição pode ter ordem semelhante à dimensão das menores estruturas. Segmentar uma estrutura de poucos pixels não recupera seu diâmetro real sem considerar desfoque, refração e resolução. A [auditoria de calibração](../../v15.1_gas_liquid_velocity/findings/calibration_review.md) e o [inventário da câmera](../../v15.1_gas_liquid_velocity/results/calibration_inventory.csv) documentam os valores por gravação.

Uma análise descritiva adicional possível é investigar velocidade versus tamanho, separando posição radial e usando trajetórias inteiras como unidades de amostragem. Extrapolar uma regressão para tamanho zero poderia gerar um candidato a velocidade líquida **somente sob um modelo de deslizamento previamente especificado**. Sem objetos próximos de tamanho zero, esse intercepto é extrapolação; covariação entre tamanho e posição, profundidade desconhecida e seleção dos objetos rastreáveis podem produzir um intercepto artificial. Bom ajuste ou validação entre trajetórias não validaria a identidade material do intercepto.

Não foi usado valor presumido de viscosidade, diâmetro físico de cada objeto ou erro de deslizamento para declarar essas microbolhas traçadores líquidos.

## 5. Onda anular e conservação mássica

No modelo local de filme plano, laminar, com cisalhamento `τ` e força volumétrica/gradiente de pressão `K` independentes da espessura `h`, tem-se `q = τh²/(2μ)+Kh³/(3μ)`, `U=q/h` e celeridade cinemática `c_w=dq/dh`. **Se** ambos os termos de condução forem não negativos, resulta `c_w/3 ≤ U ≤ c_w/2`. O intervalo exige pequena perturbação cinemática, ausência das contribuições relevantes de inércia/curvatura e constância dos termos de condução ao variar `h`. Não há evidência suficiente aqui para impor esse modelo às ondas anulares observadas. Um limite deduzido de um modelo inaplicável não é um limite físico do experimento.

Com fração volumétrica média `ᾱ`, densidades e velocidades médias **ponderadas pela fase** bem definidas, a conservação mássica fornece `G = ρ_g ᾱ U_g + ρ_l (1−ᾱ) U_l`. Se `U_g` fosse identificado independentemente, essa equação poderia inferir `U_l`; o resultado dependeria dos dados de vazão, densidades e geometria. A auditoria anterior encontrou vazão em g/s por condição, mas não comprovou sincronização dos registros LabVIEW com cada vídeo. Não é correto associar suas amostras aos quadros por igualdade de contagem ou ordem.

Sem uma velocidade de fase independente ou outro fechamento, essa equação define uma família de pares. Exigir apenas velocidades positivas permite razões de gás/líquido arbitrariamente pequenas ou grandes ao aproximar os extremos admissíveis. Relativizar não elimina essa não unicidade. Também não se deve usar `média(α) × média(U)` no lugar de `média(αU)` sem justificar a ponderação e a covariância.

## 6. Controles executados e recomendação operacional

O script [physical_bounds.py](../scripts/physical_bounds.py) implementou as fórmulas de limite, sem atribuí-las a gravações. Foram executados:

- 25 células sintéticas com perfil de gás de pontas suaves, corpo de fração constante e pistão sem gás; conservação e médias ponderadas foram verificadas.
- 2000 controles de monotonicidade dos limites sob intervalos geométricos.
- 2000 controles da aritmética dos intervalos de comparação entre células.
- 101 decomposições do modelo plano de filme entre contribuição de cisalhamento e pressão.
- Dois contraexemplos de ramo que violam os limites quando o filme escoa para trás ou mais rápido que a bolha.

O maior resíduo de conservação foi `5,03×10⁻¹⁶`; o maior resíduo da média da célula foi `8,88×10⁻¹⁶`. Esses números medem a consistência aritmética dos controles adimensionais. **Não são incertezas experimentais nem validação hidrodinâmica das gravações.** Os parâmetros dos controles não são propriedades de R600a. A [saída completa](../results/physical_bounds_controls.json) registra as hipóteses e o hash do script.

Para os resultados reais, recomenda-se:

1. Priorizar bolhas alongadas isoladas de plug/slug, conferindo correspondência, frente/traseira, comprimento, volume reconstruído e contexto adjacente nos quadros nativos.
2. Publicar a translação de bolhas elegíveis como resultado de objeto, distinguindo centroide projetado de centroide volumétrico condicionado à geometria.
3. Se a seção do corpo e o pistão puderem ser sustentados, calcular `U_s/c ∈ [α_b,1]` com sensibilidade geométrica explícita e rotular o intervalo como modelo condicionado ao filme.
4. Calcular `U_l,cell/c` somente quando a célula inteira e sua fração volumétrica correspondente estiverem disponíveis. Sem a célula, manter esse resultado ausente.
5. Não estender essas fórmulas a churn, bolhas dispersas ou anular apenas para completar pontos. Nessas situações, documentar os observáveis efetivamente recuperados e as hipóteses ainda não testadas.

Esta alternativa acrescenta informação física relativa para subconjuntos elegíveis, mas não torna identificável, por si só, um par de velocidades materiais em todos os padrões.
