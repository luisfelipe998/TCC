# Protocolo fixado para a comparação de perfis

Versão 1, registrada em 24/10/2026 antes de calcular a separação por rótulos desta análise. Avaliação interna e exploratória; dados e resultados anteriores já conhecidos. A instrução posterior do usuário definiu **v16_spatial2_18_71** como base e substituiu a entrada 20–59 proposta no plano original.

## Fontes e escopo

Entrada principal: média das linhas **18–71 inclusive**, indexação iniciada em zero, diretamente dos panoramas v11. A costura original continuou baseada em 20–59. O CSV auditado definiu 31 entradas, 26 aquisições incluídas, 87 pares intrapadrão e 238 entre padrões. Não houve nova exclusão por resultado. Uma aquisição de bolhas dispersas participou dos controles cruzados, sem réplica intrapadrão.

O núcleo da base foi importado sem modificação, com hash registrado. As páginas impressas 411, 413 e 414 de Schuck Jr. e Bodmann (2019), DOI 10.1007/978-3-030-16077-7_32, foram conferidas visualmente. A fonte primária e o contrato do plano prevaleceram sobre fórmulas convencionais inconsistentes presentes no material auxiliar da habilidade de redação.

## Representação principal e suporte

Cada aquisição inteira foi preparada com inversão 255−I, remoção de tendência linear, média zero e desvio padrão unitário. Utilizou-se K=1, correlação enviesada 1/N, CWT Morlet com ω₀=6 e a conversão λ=4πa/(6+√38). O divisor 1/N é adaptação declarada dos estimadores não enviesados do artigo. A WPS alternativa é o módulo da CWT da autocorrelação, sem quadrado. A WCS alternativa é complexa; a WCO alternativa manteve o numerador ao quadrado e o denominador de WPS no lag zero.

A grade herdada contém 200 pontos logarítmicos entre 50 e 3.000 px. A faixa principal foi definida somente pela grade e pelo menor comprimento incluído, N=3.171: conservar pontos com N−1≥√2a. Isso resultou em **187 pontos entre 50 e 2.295,940454 px**, comuns a todas as aquisições. Trata-se de regra geométrica de suporte no lag zero, não de significância nem de exigência adicional de Schuck. Não se selecionaram bandas pela separação de padrões. A margem no menor registro não garante estacionariedade nem elimina todos os efeitos de borda.

A sensibilidade conservadora conservou |b|+3a≤N/2 em b=0: **117 pontos, 50–543,851438 px**. A faixa integral 50–3.000 px permaneceu como sensibilidade, explicitando o suporte insuficiente no registro mais curto. As três faixas foram normalizadas separadamente; seus valores de distância não foram tratados como a mesma estatística.

Para cada faixa, pesos trapezoidais positivos w na coordenada log(λ) definiram pⱼ=wⱼP(0,aⱼ)/ΣwₖP(0,aₖ). Foram preservados P, a área e p. Essas massas descrevem a forma espectral e **não são frações de energia física**. A área degenerada invalidaria o perfil; não se inseririam distribuições artificiais.

A distância principal foi Hellinger, H(p,q)=√[Σ(√pⱼ−√qⱼ)²/2], de 0 a 1, menor indicando formas mais próximas. A distância não recebeu o nome de WCO nem foi atribuída a Schuck. O cosseno foi a única sensibilidade de métrica, com produto interno ΣwⱼPⱼQⱼ, usando os mesmos pesos de integração.

## Controles fixados

1. Registros inteiros; segmentos de início, centro e fim com N=3.171, recortados antes do preparo. Reutilizar os perfis auditados quando compatíveis. Segmentos sobrepostos não são réplicas.
2. Extração 30–49 como sensibilidade já explorada e 20–59 como referência histórica; nenhuma seleção posterior da melhor faixa.
3. Remoção somente da média para quantificar sensibilidade à tendência linear.
4. Médias descritivas da WPS alternativa sobre todos os lags e dentro de N−1−|b|≥√2a e |b|+3a≤N/2, com contagens por aquisição/escala. A média conservadora foi avaliada apenas na faixa comum conservadora.
5. Perfis v5 e v16 em registros completos e no mesmo eixo de escalas a para isolar a discretização; v5 nominal em separado. Para todos os pares, perfis do prefixo comum foram calculados com v16 e com o preparo literal do v5 (preparo inteiro antes do corte). A comparação desses perfis com o escore WCO herdado usou a faixa nominal integral 50–3.000 px. A diferença de conversão de escala do v5 permaneceu identificada.
6. Controles sintéticos: perfis idênticos, amplitude global, componentes em escalas distintas, espectro semelhante com fases diferentes. São verificações de interpretação, não validação experimental.

## Comparação e relevância

O desfecho comum foi a concordância de ordenação A: probabilidade descritiva de que um par intrapadrão seja mais semelhante que um par entre padrões, atribuindo meio crédito a empates. Para Hellinger/cosseno, semelhança=−distância; para o escore WCO, maior é mais semelhante. A=0,5 é ausência de preferência de ordenação. A não é acurácia de um classificador.

Foram fixadas duas ponderações: (i) peso igual por par dentro dos grupos intra/entre; (ii) peso igual entre os quatro padrões com pares intrapadrão e entre as dez combinações cruzadas. Dentro de cada grupo, pares tiveram pesos iguais. Por padrão, comparar seus pares intrapadrão aos controles que o envolvem, dando pesos iguais às quatro outras classes. Bolhas dispersas não tiveram A intrapadrão.

O escore herdado foi mantido: máximo de WCO por escala em todos os lags, ponderado por √[Pₓ(0,a)Pᵧ(0,a)] da referência, com ambas as WPS ≥5% de seus picos e denominadores numericamente válidos. Foram separados v5 original 20–59, v5 aplicado a 18–71 e v16 corrigido 18–71. Não se compararam diretamente os valores absolutos de distância e WCO.

Não foi inventado um efeito mínimo científico ou de aplicação. Foram relatados A, diferenças absolutas e relativas, cada padrão, casos desfavoráveis e estabilidade. Uma melhora geral exigiria direção consistente nas duas ponderações, no controle de comprimento e ausência de perdas por padrão; o atendimento desses critérios descritivos ainda não constituiria validação externa. Efeito restrito ou instável seria relatado como tal. Sem evidência convincente, a conclusão seria “não demonstramos melhora robusta e geral”.

## Dependência e incerteza

A unidade experimental foi a aquisição. O estudo não ajustou parâmetros aos rótulos. A análise de influência retirou uma aquisição por vez e todos os pares que a contivessem, recalculando pesos e A. As faixas mínimo–máximo são dispersão descritiva de influência, **não intervalos de confiança nem validação preditiva deixando uma aquisição de fora**. A remoção de bolhas dispersas deixa seis combinações cruzadas entre os quatro padrões remanescentes, o que foi explicitado nas contagens de cada repetição.

Não foram usados testes sobre 325 pares independentes. Não foram calculados valores-p por permutação dos rótulos, pois a permutabilidade frente às condições experimentais e aos comprimentos não foi demonstrada. Não foi realizada busca de hiperparâmetros, divisão aleatória por par ou avaliação externa.

## Exemplos e limites

Antes das métricas, fixaram-se os primeiros IDs em ordem alfabética por padrão: m07q01_001, m05q01_001, m04q01_001, m03q05_001 e m03q10_001. Exemplos WCS/WCO: primeiro par por padrão com réplicas; controle m04q01_001 × m05q01_001. Casos extremos de concordância/discordância foram acrescentados depois, explicitamente como seleção interpretativa.

Os mapas reduzidos herdados servem somente à ilustração. Máximos integrais vieram das tabelas integrais ou de novos cálculos integrais. Não houve clipping de WCO, deformação do eixo, conversão inventada para mm/Hz ou inferência de tamanho físico das bolhas. As emendas, os deslocamentos não positivos substituídos por 1 px, os rótulos sem validação independente e o pequeno número de aquisições permaneceram limitações.
