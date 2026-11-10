# Auditoria independente da análise wavelet v16

Data: 27 de outubro de 2026. Escopo: estado ativo de `src2/v16_spatial2_18_71`, confrontado com o texto atual de `TCC/`. `_OLD` foi consultado somente para localizar a fonte bibliográfica e entender verificações anteriores; seus resultados não substituíram os resultados ativos. Os achados abaixo foram formados antes da leitura dos relatórios das outras frentes da v17.

A implementação ativa reproduziu as fórmulas e os números apresentados no TCC. Não se confirmou erro de sinal, quadrado, normalização, fator de Fourier, extração das linhas ou agregação nos cálculos verificados. Os problemas materiais estão no alcance das métricas: a WCO depende fortemente do suporte da realização, a coerência de duas realizações independentes não equivale à semelhança de suas distribuições, e a estabilidade do descritor espectral muda com a faixa e a posição observada.

Isso não invalida a conclusão descritiva do texto atual: o TCC distingue os indicadores, reconhece o domínio reconstruído e não afirma acurácia ou significância. As novas evidências mostram quais controles precisam preceder uma conclusão mais forte sobre separabilidade.

## Evidências produzidas

O script [audit_wavelets.py](../scripts/audit_wavelets.py) escreve somente em `v17/resultados/wavelets`. A reprodução usa uma preparação por mínimos quadrados em forma fechada, quadratura direta da Morlet e uma identidade espectral que evita construir as correlações usadas pelo código original. O núcleo original é importado apenas para confrontos numéricos; nenhum `main` das versões anteriores é executado.

| Verificação | Cobertura | Resultado |
|---|---:|---|
| Leitura dos pixels, média 18–71 e hashes | 31 panoramas | Concordância exata com as entradas arquivadas |
| Preparação do sinal e WPS no atraso zero | 26 aquisições, 200 escalas | Erro máximo da WPS: 1,87 × 10⁻¹⁴ |
| Pearson dos resíduos espectrais | Matriz dos 26 ensaios | Erro máximo: 7,27 × 10⁻¹⁵ |
| Reconstrução do escore a partir das potências e máximos arquivados por escala | 325 pares | Concordância com tolerância 10⁻¹² |
| Recomputação independente dos mapas integrais | 3 pares com N = 3.171, 7.584 e 22.636 | Erro máximo nas posições arquivadas: 8,68 × 10⁻¹⁴ |
| Correlação direta, quadratura complexa, troca de ordem, autocoerência e atraso conhecido | Sinais sintéticos independentes dos dados | Todos os controles passaram |
| Mudança de preenchimento por zeros, de 8 para 16 desvios da maior escala | Sinal sintético | Diferença máxima: 7,63 × 10⁻¹⁷ |
| Suporte e peso das escalas | 325 pares | Quantificados por par, sem testes que tratem pares como independentes |
| Posição e faixa espectral | 8 configurações, 26 aquisições em cada uma | Tabelas e matrizes completas preservadas |
| Ruídos independentes | 48 pares em cada um de 3 comprimentos | Efeito de comprimento isolado experimentalmente |

Os três mapas completos não constituem uma reexecução integral dos 325 mapas. A verificação de todos os escores usou os máximos por escala arquivados; a fórmula de recomputação integral foi confrontada em três comprimentos e nos controles sintéticos. Essa distinção evita atribuir à auditoria uma cobertura que ela não teve.

Arquivos de evidência: [entradas.csv](../resultados/wavelets/entradas.csv), [perfis_reproduzidos.csv](../resultados/wavelets/perfis_reproduzidos.csv), [equacoes.json](../resultados/wavelets/equacoes.json), [pares_reproduzidos.csv](../resultados/wavelets/pares_reproduzidos.csv), [descritores_reproduzidos.json](../resultados/wavelets/descritores_reproduzidos.json) e [execucao.json](../resultados/wavelets/execucao.json).

## Correspondência entre código e TCC

| Ponto | Evidência no código ativo | Evidência no TCC | Julgamento |
|---|---|---|---|
| Média das 54 linhas 18–71 | `scripts/common.py:22`; `scripts/analise.py:46–48` | `metodologia.tex:84–97` | Consistente. O intervalo Python 18:72 contém 54 linhas. |
| Inversão, remoção de reta, centragem e desvio populacional | `scripts/wavelet_core.py:18–28` | `metodologia.tex:111–140` | Consistente. A inversão comum não modifica as métricas baseadas em módulo e autocorrelação. |
| Prefixo comum antes da preparação | `scripts/analise.py:150–153` | `metodologia.tex:111–114` | Consistente; cada parceiro pode mudar o trecho usado por um ensaio. |
| Correlação enviesada e sinal do atraso | `scripts/wavelet_core.py:30–43` | `metodologia.tex:146–172` | Consistente e explicitamente diferente do divisor N−|m| de Schuck e Bodmann. |
| Morlet ω₀ = 6 e conversão λ = 1,03304365a | `scripts/wavelet_core.py:11–12,45–70`; `scripts/common.py:21` | `metodologia.tex:178–223` | Consistente. Quadratura e identidade espectral confirmaram a implementação. |
| WPS = módulo da CWT da autocorrelação | `scripts/analise.py:90–94`; `scripts/wavelet_core.py:87–95` | `fundamentacao.tex:633–657`; `metodologia.tex:229–255` | Consistente. Acrescentar um quadrado à WPS seria mudar a definição alternativa. |
| WCO normalizada pela WPS no atraso nulo | `scripts/wavelet_core.py:76–83` | `fundamentacao.tex:643–646`; `metodologia.tex:233–241` | Consistente. Não falta o operador de suavização da fórmula clássica: a execução usa outra definição. |
| Máximo por escala, peso geométrico e piso de 5% | `scripts/analise.py:170–186` | `metodologia.tex:261–291` | Consistente. O piso não é um limiar de significância. |
| Perfil completo e resíduo de log-potência | `scripts/estatistica.py:12–32,49–62` | `metodologia.tex:298–340` | Consistente; amplitude, inclinação e fase não são conservadas. |
| Médias por par e por grupos | `scripts/estatistica.py:86–115,143–149` | `metodologia.tex:346–366`; `resultados.tex:340–343` | Reproduzidas. A interpretação da ponderação ainda depende da composição dos grupos. |
| Redução dos mapas somente para figura | `scripts/analise.py:170–213` | `metodologia.tex:394–398` | Consistente. Recalcular máximos somente dos mapas reduzidos produziria outra medida. |
| Eixo reconstruído e hipóteses de estacionariedade | `scripts/common.py:54–61`; dados de montagem | `metodologia.tex:99–105,135–140`; `discussao.tex:273–289` | Limitações reconhecidas corretamente, não erros de unidade ocultos. |

## W16-01 Alta relevância metodológica — comprimento e resolução influenciam a WCO

**Estado:** confirmado como mecanismo e associação descritiva; não é um bug de implementação nem prova de que todo o contraste entre padrões seja causado por comprimento.

O escore maximiza a WCO em todos os atrasos, com a mesma faixa de 50–3.000 px para prefixos de comprimentos diferentes (`scripts/analise.py:150–176`). Nos dados, o comprimento comum teve mediana 9.344 px nos pares intrapadrão e 7.358 px nos interpadrões. A correlação de Spearman entre comprimento comum e escore foi **−0,799674**. Entre comprimento e máximo integral, foi −0,458264. São associações entre pares dependentes; não foram calculados valores-p.

Para isolar o mecanismo, geraram-se duas sequências gaussianas independentes, com a mesma distribuição, em cada uma de 48 replicações por comprimento. Aplicaram-se o mesmo preparo, o divisor N e a mesma Morlet. Na escala λ = 3.000 px:

| N | N/λ | Mediana do máximo em atraso | Percentil 5–95% | Fração dos pares com máximo > 0,9 |
|---:|---:|---:|---:|---:|
| 3.171 | 1,057 | 0,968509 | 0,734302–0,998393 | 85,42% |
| 8.192 | 2,731 | 0,695718 | 0,353387–0,933647 | 12,50% |
| 23.079 | 7,693 | 0,330720 | 0,207946–0,609796 | 0% |

Em λ = 1.024 px, as medianas foram respectivamente 0,526753, 0,367015 e 0,204309. Portanto, o efeito não ficou restrito ao último ponto da grade. Mesmo no atraso zero, a WCO média em λ = 3.000 foi 0,872233 no registro curto. A maximização aumenta o resultado, mas a resolução espectral finita também importa antes de maximizar.

As quantidades são resultados de simulação, não probabilidades calibradas dos pares reais. Ruído branco não preserva o espectro, a marginal ou os eventos dos panoramas. A simulação prova que máximos próximos de 1 podem ocorrer sem dependência entre processos sob esse suporte; não demonstra que um par observado específico seja ruído.

**Contraprova à interpretação excessiva:** 112 de 325 pares atingiram o máximo integral em λ = 3.000, mas as escalas com menos de três comprimentos de onda no registro tiveram **peso mediano zero** no escore ponderado. Nenhum par concentrou mais de 50% do peso nessas escalas; o maior peso foi 17,47%. Assim, não é defensável afirmar que “a faixa acima de N/3 explica sozinha o escore”. O peso geométrico protege parcialmente o resumo. O maior peso abaixo de cinco ciclos foi 50,10%, e a mediana foi 2,14%. Os limites de três e cinco ciclos são diagnósticos de suporte, não fronteiras universais de validade.

O TCC já informa o suporte mínimo de 1,057 ciclos (`metodologia.tex:400–407`) e que o piso não é significância (`metodologia.tex:280–286`). A informação nova é a magnitude do efeito nulo e a associação entre o suporte e o indicador nos dados reais.

**Consequência:** a WCO não deve receber a interpretação de evidência de semelhança estatística somente porque seu máximo é alto. Comparações entre padrões com comprimentos distintos precisam controlar suporte e forma espectral.

**Caminho recomendado:** usar comprimento comum entre todos os ensaios em um controle primário; dividir registros mais longos em janelas para avaliar heterogeneidade e resumir por aquisição. Calibrar o escore completo, incluindo a seleção por potência e a maximização, com substitutos que preservem propriedades relevantes de cada sinal. Uma grade em λ/N é útil para estudar o estimador, mas muda a escala em pixels entre ensaios e não substitui automaticamente uma comparação física em λ comum.

Evidências: [suporte_pares.csv](../resultados/wavelets/suporte_pares.csv), [associacao_comprimento.json](../resultados/wavelets/associacao_comprimento.json), [nulos_brancos.csv](../resultados/wavelets/nulos_brancos.csv) e [resumo_nulos_brancos.csv](../resultados/wavelets/resumo_nulos_brancos.csv).

## W16-02 Alta relevância conceitual — coerência entre realizações não é semelhança entre processos

**Estado:** limitação do objetivo respondido pela métrica, não inconsistência da equação implementada.

Para dois processos independentes, de média zero, a esperança da correlação cruzada é zero, mesmo quando ambos têm a mesma autocorrelação, o mesmo espectro e pertencem à mesma população. Na estimativa de uma realização finita, a WCO pode ficar alta por flutuação, suporte reduzido ou por um alinhamento particular, como mostrou W16-01. O indicador mede associação entre as estruturas realizadas dos dois registros, não a igualdade de suas leis estatísticas.

A fonte usada pelo TCC ilustra um problema diferente: detectar o mesmo transiente incorporado a dois sinais, com um atraso conhecido e ruídos distintos (`_OLD/literature/schuck-page-06.txt:31–44`, p. 414 do PDF). Transferir esse objetivo para duas aquisições separadas de um mesmo regime exige uma justificativa adicional. Duas realizações independentes de um regime podem ser semelhantes em distribuição sem compartilhar um evento que a coerência deva recuperar.

O texto atual é prudente: distingue o Pearson espectral do escore e explica que a WCO intrapadrão menor não invalida a semelhança dos resíduos (`discussao.tex:21–35`). A auditoria sustenta essa prudência. A hipótese de que o escore necessariamente deveria ser maior dentro de cada classe é que precisa ser limitada.

**Caminho recomendado:** para semelhança de padrões, priorizar distância entre espectros e distribuições de descritores por aquisição. Reservar a WCO para hipóteses de estrutura compartilhada e deslocamento, incluindo controles com sinais derivados do mesmo registro e eventos efetivamente comparáveis. Distâncias de Hellinger, Jensen–Shannon ou Wasserstein sobre uma representação espectral comum podem ser avaliadas sem treinar um classificador. A normalização deve indicar se se compara amplitude, inclinação ou forma.

## W16-03 Relevância moderada — faixa e posição modificam a estabilidade do Pearson residual

**Estado:** confirmado em controles exploratórios; não demonstra mudança de padrão no interior de uma aquisição e não identifica uma configuração ótima.

O perfil fixo por aquisição evita a ambiguidade de recalcular o descritor para cada parceiro (`scripts/estatistica.py:49–62`; `metodologia.tex:336–340`). Entretanto, “fixo” não significa “estável sob novo trecho observado”. Na v17, foram escolhidas três posições por regra geométrica — prefixo, centro e sufixo —, todas com N = 3.171 px, o menor comprimento disponível. A regra foi aplicada a todos os ensaios, sem escolher janelas pela semelhança.

| Configuração | Média intra | Média inter | Δ por par |
|---|---:|---:|---:|
| Panoramas completos, 50–3.000 | 0,869718 | 0,803828 | 0,065889 |
| Panoramas completos, 50–1.024 | 0,757445 | 0,582220 | 0,175226 |
| Prefixos de 3.171, 50–1.024 | 0,490588 | 0,303344 | 0,187245 |
| Centros de 3.171, 50–1.024 | 0,696998 | 0,426011 | 0,270987 |
| Sufixos de 3.171, 50–1.024 | 0,607685 | 0,460229 | 0,147456 |
| Prefixos de 3.171, 50–3.000 | 0,793236 | 0,699597 | 0,093639 |
| Centros de 3.171, 50–3.000 | 0,766955 | 0,686200 | 0,080755 |
| Sufixos de 3.171, 50–3.000 | 0,756878 | 0,734376 | 0,022502 |

Todos os contrastes agregados permaneceram positivos nessas oito configurações. Isso é uma contraprova contra afirmar que a tendência desaparece inevitavelmente ao controlar comprimento. Entretanto, a magnitude mudou substancialmente, e o contraste agregado não mede separação de cada padrão.

Para o mesmo ensaio `m07q20_001`, o Pearson entre os resíduos do prefixo e do sufixo, na faixa 50–1.024, foi **−0,732192**. A mediana das comparações prefixo–sufixo das 25 aquisições restantes após excluir o menor ensaio, cujas janelas são idênticas, foi 0,734817. Na faixa 50–3.000, a mediana foi 0,793633 e o mínimo 0,389672. As janelas podem se sobrepor; esses números não são estimativas baseadas em 25 réplicas independentes do mesmo processo.

Diferenças entre janelas podem refletir heterogeneidade real, montagem, ruído, não estacionariedade e erro de estimação. Este experimento não separa essas causas. A faixa 50–1.024 aumenta o número mínimo de ciclos, mas também remove informação de escala longa; “melhor suporte” não garante automaticamente melhor descrição física.

**Caminho recomendado:** apresentar a estabilidade dentro da aquisição ao lado da distância entre aquisições; usar várias janelas com a mesma regra e um resumo por ensaio; fixar a faixa por suporte e finalidade física antes de examinar o contraste. Não escolher “centro, 50–1.024” por ter produzido o maior Δ desta tabela. Para discutir separabilidade, examinar padrões individualmente, sobreposição das distribuições e composição dos grupos, incluindo o singleton de bolhas dispersas.

Evidências: [sensibilidade_janelas.csv](../resultados/wavelets/sensibilidade_janelas.csv), [repetibilidade_interna_janelas.csv](../resultados/wavelets/repetibilidade_interna_janelas.csv), os oito arquivos NPZ no mesmo diretório e [janelas_metodo.json](../resultados/wavelets/janelas_metodo.json).

## W16-04 Relevância moderada — comparar os dois indicadores também muda o suporte

**Estado:** limitação já explicitada no TCC, confirmada pelo código.

O Pearson compara panoramas completos e usa um perfil fixo por aquisição; a WCO usa prefixos comuns diferentes para cada par (`scripts/analise.py:150–153` versus `scripts/estatistica.py:49–62`). Para o par do registro mais curto com o mais longo, o trecho do segundo é reduzido de 23.079 para 3.171 px. Como o preparo ocorre depois do corte, a reta, a média e o desvio também mudam.

O TCC declara essa diferença em `resultados.tex:65–68` e `conclusao.tex:22–24`. Portanto, os resultados não estão apresentados como uma comparação experimental isolada entre fórmulas. Qualquer afirmação futura de superioridade de um método, porém, precisa manter a observação comparável.

**Caminho recomendado:** comparar ambas as métricas nas mesmas janelas previamente definidas, com N, extração, banda e ponderação fixos. Manter a análise dos panoramas completos como um resultado complementar que responde a outra pergunta. Quando a WCO é usada com um máximo por escala, relatar separadamente um alinhamento único e a curva de máximos; eles não têm a mesma interpretação.

Os três mapas recomputados também mostram que restringir |lag| ≤ N/4 diminuiu os escores de 0,299186 para 0,234914, de 0,252630 para 0,178669 e de 0,071443 para 0,060687. É uma sensibilidade ao domínio de busca, não uma proposta de limiar universal. Os valores de |lag| ≤ N/2 ficaram arquivados para evitar selecionar somente uma restrição.

## W16-05 Oportunidade de análise — preservar informações removidas pelo resíduo

**Estado:** propriedade matemática deliberada do descritor, não erro.

O ajuste retira o nível e a inclinação de log(P) em função de log(λ) (`scripts/estatistica.py:12–21`; `metodologia.tex:298–318`). Desprezando o piso ε quando P ≫ ε, multiplicar um perfil por Aλᵝ não muda o resíduo. Assim, dois perfis com amplitudes e inclinações diferentes podem produzir Pearson igual a 1.

Isso é adequado quando a hipótese se refere somente à curvatura residual, mas descarta diferenças que também podem carregar informação sobre os padrões. Nos panoramas completos, as inclinações médias foram aproximadamente −0,158 em churn, −0,134 no anular, 0,022 em plug e 0,393 em slug; bolhas dispersas teve −0,759 em um único ensaio. Há sobreposição dos intervalos, e as inclinações mudam com a faixa e o trecho. Esses números não validam um discriminador.

**Caminho recomendado:** conservar lado a lado a inclinação, a amplitude de flutuação antes da padronização e a forma residual, cada uma com sua interpretação e seus controles ópticos. Uma análise de efeitos e distâncias com esses descritores não exige aprendizado de máquina. Não reunir automaticamente todos em uma distância ajustada para maximizar a separação observada.

Evidência exploratória: [inclinacoes_espectrais.csv](../resultados/wavelets/inclinacoes_espectrais.csv), para todas as oito configurações.

## Questões que não se confirmaram como defeitos

- **Ausência de suavização clássica:** não foi omissão da implementação; a CWT é aplicada às correlações e a fórmula alternativa foi declarada.
- **Ausência de COI nos mapas de correlação:** o TCC a informa. Um COI da CWT direta no tempo não deve ser transplantado mecanicamente para o eixo de atrasos. O preenchimento por zeros evita imagens periódicas da FFT, mas não aumenta o suporte de estimação. Controles de comprimento e domínio de busca são mais informativos do que simplesmente chamar toda a figura de inválida.
- **WCO limitada a [0,1]:** para o estimador enviesado e o multiplicador de Fourier não negativo, a desigualdade de Cauchy–Schwarz fornece esse limite. Estar nesse intervalo verifica uma propriedade algébrica, não a validade física da montagem ou a estacionariedade do escoamento.
- **Divisor N em vez de N−|m|:** mudança metodológica explicitada, com vantagem de preservar a estrutura espectral positiva do registro finito. Trocar pelo estimador não enviesado não é uma correção automática e pode produzir coerências fora de [0,1].
- **Inversão 255−I:** modifica o sinal preparado por um sinal global, mas não as autocorrelações, WPS e módulos usados. Não há ganho físico demonstrado apenas por essa inversão; tampouco se encontrou erro numérico causado por ela.
- **Pares e escalas como réplicas independentes:** as estatísticas ativas são descritivas e o TCC nega explicitamente essa independência. Não se confirmou o uso de teste de Pearson com 200 escalas como se fossem amostras independentes para sustentar o resultado.
- **Incerteza óptica ou falha física dos deslocamentos:** a igualdade dos pixels e dos hashes não resolve essas questões. Esta frente não executou nova estimação dos deslocamentos dos 15.531 quadros nem calibrou intensidades ou comprimentos físicos.

## Próximos passos sem aprendizado de máquina

1. Fixar um protocolo de janela e faixa pelo suporte e pela escala de interesse, preservando resultados de sensibilidade sem selecionar a maior separação.
2. Separar semelhança de distribuição espectral de compartilhamento de eventos; usar a WCO apenas quando a segunda hipótese fizer sentido.
3. Comparar PSD em janelas e distâncias de distribuições com o Pearson residual, preservando a aquisição como unidade de agregação.
4. Medir variação dentro de cada aquisição e sua relação com variação entre aquisições, sem transformar janelas sobrepostas em novas réplicas.
5. Incluir controles condicionados às aquisições para máximos de WCO, preservando marginais e estrutura espectral quando esse for o nulo pertinente.
6. Examinar composição dos grupos, condições operacionais e qualidade óptica antes de atribuir efeitos ao rótulo do padrão.

A [figura de diagnóstico](../resultados/wavelets/diagnosticos_wavelets.png) reúne o efeito de suporte nos ruídos e a sensibilidade do contraste. Nenhuma configuração foi promovida a método final, nenhum classificador foi treinado e nenhum arquivo do TCC ou das versões anteriores foi alterado por esta frente.
