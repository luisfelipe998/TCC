# Método e limites da exploração WCO

## Pergunta e relação com as análises anteriores

A pergunta foi se a alternative WCO permitia aproximar ensaios do mesmo padrão e afastar os de outros padrões. A `estatistica5` usou correlação de forma espectral da WPS, não correlação entre mapas WCO. A `estatistica6` testou WPS/mean power e três resumos WCO arquivados; estes últimos usavam comprimento mínimo de cada par. Aqui, a WCO foi recalculada com **um trecho fixo por ensaio em cada condição**, válido em todos os seus pares.

A [configuração principal](configs/wco_v1.json) foi congelada antes dos novos mapas, após conhecimento da estatistica5/6. Não se trata de pré-registro independente. O [adendo](configs/diagnostics_addendum.json), explicitamente posterior à primeira inspeção, introduziu normalização pela autocoerência e controles de ruído. Não acrescentou testes à família inferencial de seis contrastes.

## 1. Unidade experimental e dados

Foram mantidos 26 ensaios: churn 10, wavy annular 9, plug 3, slug 3 e dispersed bubbles 1. As cinco exclusões da origem foram preservadas: `m03q01_001`, `m03q80_001`, `m04q05_001`, `m05q05_002`, `m07q10_001`. O singleton participou da descrição de pares, mas não das medidas de compactação intraclasse.

Usou-se a média das linhas 18–71 inclusive dos panoramas, armazenada nos arquivos `_input.npz` da v16. A extração dos pixels já havia sido auditada nas etapas anteriores. Nesta execução, os hashes das fontes foram conferidos; não houve reconstrução dos panoramas. Os campos históricos de proveniência do inventário foram conservados, enquanto os caminhos efetivamente lidos constam no manifesto.

| N (px) | Participantes em cada posição | Pares por posição | Posições |
| --- | --- | --- | --- |
| 3.171 | Todos os 26 | 325 | Começo, centro, fim |
| 5.071 | Os mesmos 25 de classes replicadas | 300 | Começo, centro, fim |
| 7.001 | Os mesmos 19 churn/wavy annular | 171 | Começo, centro, fim |

Total: 2.388 mapas cruzados e 210 mapas de autocoerência. Cortou-se o sinal bruto **antes** da inversão `255 − intensidade`, retirada da tendência linear e padronização pelo desvio-padrão populacional. O centro iniciou em `floor((N_original − N)/2)`. As coordenadas estão em [window_coordinates.csv](results/window_coordinates.csv).

As comparações entre comprimentos de cada escopo conservaram os mesmos ensaios. As distribuições globais de 26, 25 e 19 ensaios não foram tratadas como um experimento puro de comprimento. Janelas sobrepostas ou idênticas são análises de sensibilidade, não novas aquisições. Para o slug `m04q01_001`, os três trechos de 5.071 px são idênticos.

## 2. Estimador e domínio preservados

Foi mantida a definição da [auditoria de equações da v16](../v16_spatial2_18_71/EQUATION_AUDIT.md), com correlações enviesadas, divisor N, Morlet não corrigida de frequência angular 6 e 200 comprimentos de onda geométricos entre 50 e 3.000 px:

\[
P_i(a)=|\operatorname{CWT}\{\widehat R_{ii}\}(0,a)|,\qquad
\Gamma_{ij}(b,a)=\frac{|\operatorname{CWT}\{\widehat R_{ij}\}(b,a)|^2}{P_i(a)P_j(a)}.
\]

O denominador usa **lag zero**, não o valor em cada lag. A alternative WPS é a magnitude da CWT da autocorrelação, sem um segundo quadrado. Isso difere da coerência wavelet convencional baseada nas CWTs diretas dos sinais. Há uma realização por aquisição, K=1.

Foram utilizados todos os `2N−1` lags, sem COI, clipping, alinhamento ou seleção de regiões dos mapas. O piso numérico exigiu cada P acima de `10⁻¹⁰` do próprio máximo; as 200 escalas passaram em todos os 2.388 mapas observados. O piso não define significância. Escalas e lags estão em pixels.

O motor de convolução foi acelerado calculando somente as posições necessárias, com suporte suficiente para reproduzir a convolução original. Foram conferidos CWT, reversão dos pares, convolução temporal independente e quatro mapas integrais de pares reais. O maior erro absoluto WCO frente ao núcleo original foi **7,55×10⁻¹⁵**. A identidade de autocoerência em lag zero foi verificada como propriedade algébrica, não como validação experimental.

## 3. Representações testadas

Sejam `h_a` os pesos trapezoidais da grade em log(λ), `u_a = h_a/Σh_a` e `w_a ∝ h_a√(P_iP_j)`, normalizados entre escalas válidas. As bandas foram [50,200), [200,800) e [800,3000] px, usando a restrição desses pesos à banda.

| Família | Definição dos resumos |
| --- | --- |
| Máximo global | `max_(a,b) Γ`. |
| Lag zero | `Σ_a w_a Γ(0,a)`. |
| Máximos por escala | `Σ_a w_a max_b Γ(b,a)`. Permite um lag diferente por escala. |
| Máximos com piso de potência | Pesos `√(P_iP_j)`, sem quadratura, exigindo P de ambos ≥5% do próprio máximo; reproduz o resumo anterior. Pisos 1%/10% foram salvos como diagnóstico. |
| Máximo comum | `max_b Σ_a w_a Γ(b,a)`, exigindo um único lag na média das escalas. |
| Médias integrais | Média sobre todos os lags, seguida de média com `u_a` ou `w_a`. |
| Percentil e área | Percentil 90 sobre lags, ponderado por `w_a`; fração ponderada do domínio com Γ≥0,5. O limiar 0,5 é descritivo. |
| Bandas | Média integral não ponderada por potência e máximo comum ponderado, em cada uma das três bandas. |
| Autocoerência | `A_i(a)=média_b Γ_ii`; primeira meia largura positiva de Γ_ii<0,5, dividida por λ. Não houve censura em nenhum perfil calculado. |

Foram 15 resumos cruzados e três representações próprias: distância L2 entre log-perfis de área; L2 entre log-perfis de meia largura; distância entre três áreas logarítmicas de banda, padronizadas por mediana/IQR. A escala desse último vetor foi ajustada nos 25 prefixos e mantida nas demais condições; não é validação fora da amostra.

Para os resumos cruzados usou-se a dissimilaridade `d_ij=1−S_ij`, com diagonal zero. Não se assumiu que todas fossem métricas euclidianas. Como essa convenção pode comprimir silhuetas quando o escore tem autossimilaridade menor que 1, o adendo verificou uma normalização pela diagonal observada:

\[
K_{ij}=\sum_a u_a\,\overline\Gamma_{ij}(a),\quad
C_{ij}=\frac{K_{ij}}{\sqrt{K_{ii}K_{jj}}},\quad
d_{ij}=\sqrt{2(1-C_{ij})}.
\]

Aplicou-se à área integral e às três bandas; o quinto candidato foi o kernel de Γ em lag zero, com pesos comuns e diagonal 1. Exigiu-se matriz empírica semidefinida positiva e C≤1, tolerância `10⁻¹⁰`, antes de usar essa geometria. **Seis das 45 matrizes foram recusadas:** área normalizada longa, em todas as posições de 3.171/5.071 px. As demais 39 foram aceitas. Não houve correção de autovalores para forçar aceitação. Em 7.001 px a banda longa foi admissível, mas também não mostrou separação estável. A rejeição desse kernel não invalida os mapas WCO nem impede a análise descritiva da média longa original.

## 4. Critérios de utilidade e estabilidade

Foram avaliados todos os seis contrastes possíveis entre as quatro classes replicadas e o conjunto de quatro classes. Priorizaram-se churn/wavy annular, plug/slug e quatro padrões.

A silhueta por ensaio comparou a distância média à própria classe, excluindo o próprio ensaio, com a menor distância média a outra classe. A macro foi a média das médias de classe, dando peso igual aos padrões. Registrou-se também a menor média de classe. O vizinho mais próximo usou crédito proporcional nos empates e média balanceada por classe. Sua referência aleatória foi `média_g[(n_g−1)/(n−1)]`, e não 50% por convenção; são 47,2% para churn/wavy annular e 40% para plug/slug.

A triagem exigiu macro≥0,10, cada média de classe>0 e vizinhança≥referência+0,10. Para continuidade, exigiu-se cada classe positiva e vantagem de vizinhança em todas as posições de 3.171/5.071 px, também 7.001 px em churn/wavy annular, além de macro positiva em todas as retiradas de um ensaio. Esses são critérios exploratórios de trabalho, não limiares de significância. Nenhuma representação passou a triagem inicial. **Mesmo sem o limite macro de 0,10, nenhuma manteve todas as classes positivas nas seis condições comuns.**

O agrupamento usou ligação média com k=2 ou k=4 fixado, avaliado por ARI, sem busca de algoritmo/k. Não foram treinados classificadores. Silhueta, ARI e concordância de vizinhos nesta amostra não estimam acurácia futura.

Os perfis funcionais compararam, em todas as escalas, a média das médias intraclasse com a média interclasse de Γ médio, Γ máximo e Γ em lag zero. O envelope ao retirar um ensaio é sensibilidade à influência, não intervalo de confiança. As 200 escalas não foram usadas como 200 réplicas.

### Sinais condicionais conservados

Em 5.071 px, o máximo comum na banda 200–800 px teve ambas as classes positivas nas três posições de wavy annular/slug, com menor média de classe entre 0,0018 e 0,0239; falhou no centro em 3.171 px. O kernel de lag zero também teve médias positivas em churn/slug nas três posições de 5.071 px, porém macros pequenas (0,0107–0,0181), e falhou nas três posições de 3.171 px. Plug/slug tiveram resultado favorável para o máximo global somente no centro/fim de 5.071 px. Esses achados estão preservados nas tabelas; não sustentaram estabilidade espacial e entre comprimentos.

## 5. Controles e referência estatística

**Independência conhecida:** foram simulados 99 pares independentes em cada combinação de N=3.171/5.071/7.001 e ruído branco gaussiano/AR(1) gaussiano com φ=0,95, total 594 pares. O AR(1) descartou 1.000 amostras iniciais. Cada par recebeu sementes independentes e o mesmo preparo/cálculo dos dados. Os intervalos mostrados são percentis das realizações, não IC populacional. O máximo restrito a λ≤800 px foi salvo apenas como diagnóstico de suporte; não foi escolhido como novo classificador. O controle demonstra que máximos altos são possíveis sem relação entre sinais, sem identificar a causa exclusiva dos valores reais.

**Fases aleatórias:** os três resumos da família principal foram máximo global, máximo comum e média integral, em duas coortes: 19 churn/wavy annular e seis plug/slug. O contraste foi:

\[
T=\tfrac12(\overline S_{AA}+\overline S_{BB})-\overline S_{AB}.
\]

Cada réplica gerou um sinal substituto por ensaio, reutilizado em todos os seus pares. Os 186 pares não foram tratados como 186 testes independentes. Foram 199 réplicas, fixadas por custo antes dos mapas observados; não houve parada por valor-p. As sementes derivaram de `2026102407`, réplica e índice do ensaio.

Randomizaram-se independentemente as fases rFFT, preservando amplitudes e DC/Nyquist. Depois reaplicaram-se tendência linear e padronização. Antes dessa reaplicação, o erro de preservação das amplitudes foi <4,3×10⁻¹⁶. Depois, a mudança relativa L2 do espectro de potência teve mediana **4,27%**, P95 **16,45%** e máximo **47,20%**, em 4.975 sinais substitutos. Portanto, a referência é **aproximada**, com espectros estimados nos próprios dados, e não preserva necessariamente costuras, morfologia e não estacionariedade.

Calculou-se `p_model=(1+b)/(199+1)` para a alternativa T observado maior, com empates a `10⁻¹²`. Aplicou-se Holm conjuntamente aos seis contrastes. O menor p bruto foi 0,115 e o menor corrigido 0,69. Os intervalos binomiais salvos medem apenas erro Monte Carlo da cauda. Não são IC de diferenças de padrões. Esses valores-p de referência não são testes de permutação de rótulos, causalidade, generalização ou equivalência; a interpretação depende das aproximações do modelo.

**Operação e reconstrução:** verificaram-se associações descritivas com diferenças de m, q, comprimento original, deslocamento de reconstrução e distância Hellinger da WPS, sem p para linhas de pares. Só q=40 permitiu ambas as classes com ≥2 ensaios: churn 2, wavy annular 3. Nenhum dos três resumos principais teve ambas as médias de classe positivas nesse controle. Não há intercâmbio de rótulos dentro de (m,q), e comprimento/deslocamento original se confundem com padrão. Fixar N não elimina essa confusão.

## 6. Auditoria e reprodução

- [engine.json](validation/engine.json): equivalência de motores e controles algébricos.
- [maps.csv](validation/maps.csv): domínio e limites dos 2.388 mapas observados.
- [kernel_geometry.csv](validation/kernel_geometry.csv): normalizações aceitas e recusadas.
- [final.json](validation/final.json): resumos recalculados por via independente, silhuetas, vizinhos, Holm, coortes, contagens e fontes.
- [manifest_started.json](results/manifest_started.json) e [artifact_manifest.json](results/artifact_manifest.json): fontes, versões e hashes dos arquivos entregues.
- [README](README.md): comandos de reprodução e índice dos nove gráficos.

Nenhum arquivo de origem ou das análises 5/6 foi alterado. A precisão numérica verificada sustenta a execução das fórmulas; não demonstra representatividade ou validade física. A conclusão é restrita aos resumos, estimador, domínio e ensaios avaliados: **há um diagnóstico relevante dos máximos WCO, mas não evidência de um agrupamento de padrões estável nesta exploração.**
