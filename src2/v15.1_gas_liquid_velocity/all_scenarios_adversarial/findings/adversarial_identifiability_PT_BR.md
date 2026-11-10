# Auditoria adversarial: cobertura, médias e identificação das fases

A investigação nova encontrou uma limitação corrigível na análise anterior: exigir uma única celeridade estável durante toda a gravação pode rejeitar movimento óptico real que varia no tempo. Portanto, a cobertura anterior não demonstrava impossibilidade de ampliar o gráfico. Os métodos por blocos ampliaram essa cobertura. Entretanto, obter dois números, identificar duas populações ópticas e medir duas velocidades materiais são resultados distintos. Os controles abaixo foram executados, e os respectivos dados foram preservados.

## Resultados novos conferidos independentemente

| Método | Cobertura conferida | Significado |
|---|---:|---|
| NCC regional em 10 blocos | 20/26 pares com suporte | Média ponderada pelos intervalos dos blocos aceitos simultaneamente no núcleo e nas duas faixas externas |
| NCC distribuída dentro de máscaras, erosão de 3 pixels | 7/26 pares com suporte; 20/26 pares finitos | Movimento de intensidade em pixels candidatos de cada fase na origem e no destino |
| Mesma NCC, erosão de 1 pixel | 8/26 pares com suporte; 23/26 pares finitos | Sensibilidade a uma definição menos conservadora do interior das máscaras |
| Fourier por blocos, verificações simultâneas completas | 2/26 pares com suporte | Propagação dominante com consistência temporal, entre bandas e entre metades axiais |

Seis gravações passaram na NCC mascarada com ambas as erosões. A união das duas listas não deve substituir uma regra fixa escolhida previamente. O método regional de 10 blocos recuperou cinco gravações anulares, mas suas coordenadas continuam sendo núcleo/parede, sem atribuição material. Nos 20 pares aceitos, os intervalos dos blocos usados correspondem a 68,8–98,2% dos 500 intervalos originais. A média deve ser descrita como média sobre os blocos aceitos; ela não recupera o movimento dos intervalos omitidos.

A auditoria reconstruiu os 156 resultados regionais das seis variantes a partir dos blocos: diferença máxima de 8,74×10⁻¹¹ px/s em relação às coordenadas salvas. Conferiu também os 26 pares espectrais, com diferença máxima de 2,92×10⁻¹¹ px/s, e recalculou os critérios de aceitação das 104 estimativas mascaradas. Essas tolerâncias expressam concordância computacional, não erro experimental.

## 1. As regras anteriores realmente podiam rejeitar movimento correto

Foram produzidas três sequências sintéticas de 501 quadros, com translação conhecida e sem identificação de fases. A correspondência entre pares foi calculada separadamente do pico da correlação média.

| Controle sintético | Média temporal verdadeira, px/quadro | Média dos deslocamentos estimados | Pico da correlação média | Rejeição por estabilidade entre blocos/intervalos |
|---|---:|---:|---:|---|
| Velocidade constante | −16,0000 | −16,0000 | −16,0000 | Não |
| 300 intervalos a −8 e 200 a −24 | −14,4000 | −14,4001 | −7,9877 | Sim |
| Aceleração contínua de −4 a −28 | −16,0000 | −16,0001 | −17,1948 | Sim |

O erro absoluto médio por par foi inferior a 0,012 px/quadro nos três controles. No segundo, o pico agregado errou a média temporal por 6,4123 px/quadro: ele privilegiou o modo mais frequente. A regra que exigia concordância de pelo menos quatro dos cinco blocos rejeitou esse movimento conhecido porque somente três blocos concordavam com o modo dominante. A aceleração também alterou os picos dos intervalos de dois e quatro quadros.

Nas tabelas reais anteriores, 31 estimativas regionais eram rejeitadas somente por estabilidade entre blocos/intervalos e passavam nas outras verificações ampliadas. Isso justificou a reanálise local no tempo. Não demonstra que os 31 resultados sejam corretos, nem dispensa verificações locais.

## 2. Correlacionar máscara de gás e seu complemento duplica a informação

Em 26 gravações, três limiares e dez pares adjacentes fixos por gravação/limiar, foram examinadas 780 comparações com máscaras reais salvas. Após calcular a fração radial da máscara dentro do lúmen, aplicou-se NCC diretamente com o mesmo suporte espacial à máscara e ao complemento.

Das 780 comparações, 690 tinham correlações definidas e 90 eram constantes ou degeneradas. A maior diferença entre as curvas da máscara e do complemento foi 4,00×10⁻¹⁵. Em 668 comparações, o índice do máximo foi idêntico; nas outras 22, índices diferentes resolveram empates numéricos entre máximos equivalentes. Houve 636 comparações com máximo único à tolerância de 10⁻¹⁰.

A identidade tem origem algébrica: ao subtrair a média, `1 − máscara` torna-se o negativo da máscara centrada; os dois sinais negativos se cancelam na correlação. O movimento de uma fronteira aparece nos dois lados. A identidade **não** se aplica automaticamente à intensidade grayscale dentro de regiões disjuntas. Por isso, o novo método mascarado, que usa intensidade na origem e no destino, constitui um teste adicional legítimo.

## 3. Ponderar um único movimento pela ocupação produz 26 pontos, mas não duas medições

Foi testada a alternativa de usar um movimento global por par de quadros e ponderá-lo pela fração aparente de gás `w` e pelo complemento `1−w`. Foram executadas 156 combinações: 26 gravações, três limiares e dois filtros de correspondência.

Na referência de limiar 0,70 usando todos os candidatos finitos, obtiveram-se dois números para **26/26** gravações, com razões entre 0,6174 e 0,9994. Esses pontos seriam fáceis de desenhar. Entretanto, por construção,

`média_w(v) − média_(1−w)(v) = Cov(w,v) / [média(w) × (1−média(w))]`.

A identidade foi conferida nos dados, com resíduo máximo inferior a 3,80×10⁻⁸ px/s. Portanto, a diferença mede a associação entre aparência e **um mesmo sinal global de movimento**; ela não fornece duas correspondências materiais independentes. Quando o movimento foi substituído por uma constante, todas as razões ficaram iguais a um. Também foram executados 199 embaralhamentos das ocupações por combinação para remover a associação temporal. Seus percentis são diagnósticos de permutação, não intervalos de confiança experimentais.

O teste que usou todos os candidatos incluiu correspondências que falharam nos filtros anteriores; sua cobertura não deve ser apresentada como 26 resultados aceitos.

## 4. Mesmo a NCC dentro das fases pode acompanhar iluminação compartilhada

Um controle adicional usou máscaras de fase exatas e disjuntas, interiores uniformes e uma textura de iluminação que se deslocava a +12 px/quadro sobre ambas as regiões. Nenhum traçador material foi incluído. O método mascarado encontrou +12,00005 px/quadro em ambas, com correlação de aproximadamente 0,999906; os controles temporais e de embaralhamento passaram.

O estimador recuperou corretamente a propagação óptica apresentada. Isso não é uma falha de correspondência. O exemplo mostra por que uma máscara correta, forte correlação e reversão temporal não demonstram que a textura acompanha o fluido. Ele também não demonstra que a iluminação compartilhada seja a causa dos resultados reais.

O método NCC por blocos foi desafiado ainda com 18 sequências periódicas ambíguas, usando seis períodos e três velocidades. As 18 foram rejeitadas pelas verificações combinadas. Esse resultado sustenta a utilidade dos controles testados, mas não constitui garantia geral contra correspondências ambíguas.

## 5. Onda de filme: converter celeridade em velocidade média exige uma hipótese adicional

Foi construído um exemplo explicitamente sintético de filme laminar localmente plano, sem deslizamento na parede, viscosidade `μ`, cisalhamento interfacial constante `τ` e força de pressão/gravidade `K`. Nesse modelo,

`q(h) = τ h²/(2μ) + K h³/(3μ)`,

`c = dq/dh`, enquanto `u_médio = q/h`.

Para `h=100 μm`, `μ=0,001 Pa·s` e a mesma celeridade de pequena perturbação `c=0,6 m/s`, cinco repartições entre cisalhamento e pressão produziram velocidades médias entre **0,2 e 0,3 m/s**. No caso exclusivamente movido por pressão, `u=c/3`; no caso exclusivamente movido por cisalhamento, `u=c/2`. A espessura e a celeridade foram mantidas iguais.

São exemplos de modelo, não valores de R600a nem limites experimentalmente válidos para os anulares ondulados. Curvatura, inércia, variação do cisalhamento, entranhamento e mudança de fase exigem outras equações. Um fator universal `u=c/2`, `u=c/3` ou `u=c` preencheria coordenadas sem resolver essa escolha.

## 6. Continuidade, fração de vazio e vazão não selecionam sozinhas um ponto

Uma possibilidade legítima é utilizar toda a evolução de uma fração de vazio volumétrica `α(x,t)` e as equações de continuidade para reconstruir **variações de fluxo**. Contudo, mesmo supondo `α` medida perfeitamente e mudança de fase desprezível, integrar `∂t α + ∂x jg=0` exige um fluxo na fronteira.

O controle sintético de onda viajante usou `α=0,5+0,1 sen(x−2t)`. Dois pares de fluxos `jg=2α+Cg` e `jl=2(1−α)+Cl` satisfizeram exatamente as duas equações de continuidade e a mesma vazão mássica média total. As médias intrínsecas correspondentes foram **(4,0; 3,6)** e **(8,0; 3,2)**, em unidades adimensionais. São soluções de continuidade; não se alegou que satisfizessem uma mesma equação de momento ou fechamento hidrodinâmico.

Adicionalmente, mesmo concedendo os valores condicionais de densidade, geometria e vazão previamente usados para cada gravação real, foram construídos cinco pares positivos que satisfazem a conservação por gravação: **130 pares distintos para 26 condições**, com resíduo máximo de 2,28×10⁻¹³ kg/(m²·s). Os valores extremos podem ser fisicamente inadequados e não são medições. O exercício mostra apenas que a conservação isolada restringe uma curva, sem escolher um ponto.

O conjunto fechado ainda permite estimativas condicionais com modelos hidrodinâmicos, inversão de perda de pressão, ou relações de velocidade de bolhas aplicáveis a um regime específico. Isso exigiria declarar o fechamento e verificar geometria, pressão, propriedades e condições de entrada. Nenhuma dessas alternativas transforma automaticamente os 26 pontos em medições de fase. Tampouco justifica transferir uma relação de plug/slug para churn e anular para completar o gráfico.

## Recomendação

Publicar a ampliação de cobertura como um resultado de **movimento óptico**, com o método fixo e o domínio temporal da média explícitos. O método regional por blocos é promissor para reunir mais gravações; a NCC dentro de máscaras é mais próxima da proposta de separar aparências, mas sua pureza de fase continua por demonstrar. Para qualquer gráfico com 26 candidatas, os casos sem suporte precisam continuar identificáveis, e não se devem substituir ausências por zero nem misturar métodos caso a caso.

Não se demonstrou uma impossibilidade absoluta de extrair informação adicional das imagens. Demonstrou-se que algumas maneiras simples de obter cobertura completa — correlacionar complementos, ponderar uma só velocidade, converter ondas por um fator universal ou impor conservação sem fechamento — não identificam as duas velocidades materiais solicitadas.

## Reprodução e evidência

- Executar `python3 -B src2/v15.1_gas_liquid_velocity/all_scenarios_adversarial/scripts/adversarial_identifiability.py`.
- Depois de gerar os três métodos novos, executar `python3 -B src2/v15.1_gas_liquid_velocity/all_scenarios_adversarial/scripts/adversarial_new_methods.py`.
- [Controles e proveniência](../results/adversarial_identifiability.json).
- [Auditoria dos métodos novos](../results/adversarial_new_methods.json).
- [Máscaras reais e complemento](../results/adversarial_real_mask_complement.csv).
- [Ponderações do movimento global](../results/adversarial_weighted_common_motion.csv).
- [Translações sintéticas não estacionárias](../results/adversarial_nonstationary_synthetic.csv).
- [Conversão de onda de filme](../results/adversarial_film_wave_nonuniqueness.csv).
- [Continuidade da fração de vazio](../results/adversarial_void_continuity_nonuniqueness.csv).

As máscaras reais foram reutilizadas dos arquivos salvos, com hashes de entrada; esta auditoria não refez sua segmentação a partir de cada JPEG. A auditoria de agregações não refez todas as correspondências dos três métodos novos. Os exemplos sintéticos testam afirmações específicas e não reproduzem todas as condições ópticas ou hidrodinâmicas experimentais.
