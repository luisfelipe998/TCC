# Auditoria independente de estatistica8

Data: 27/10/2026. Escopo: dados, implementação, fórmulas, inferência e pertinência para a separabilidade e a semelhança dos padrões. Os arquivos auditados permaneceram em leitura. Todos os cálculos novos foram gravados em `src2/v17`.

**A largura WPS contém um contraste descritivo entre churn e anular ondulado, mas sua magnitude depende da região do panorama. Não foi confirmado erro aritmético material nos cálculos reconstruídos. Os maiores problemas são de alcance científico: seleção de trechos, confusão operacional e transferência indevida entre versões.** A quantização dos quantis não explicou o contraste. Nas cinco posições uniformes novas, a separação permaneceu menor que na âncora original; a mediana por aquisição conservou um efeito parcial. Estes resultados não demonstram seleção intencional para maximizar a separação.

Esta avaliação foi formada sem ler os relatórios dos outros agentes. Ela confrontou código e dados diretamente; relatórios anteriores foram usados posteriormente para distinguir limitações já reconhecidas de problemas não explicitados.

## O que foi efetivamente conferido

| Verificação independente | Resultado |
| --- | --- |
| Integridade das 479 fontes registradas pela extração | Todos os SHA-256 coincidiram |
| Silhuetas das 36 representações × 7 escopos × 3 posições | 756 configurações; erro absoluto máximo de 1,67 × 10⁻¹⁶ |
| Larguras WPS registradas | Erro máximo de 2,22 × 10⁻¹⁶ |
| Identidade entre distância residual e Pearson | `d² = 2(1-r)`; erro máximo de 1,78 × 10⁻¹⁵ |
| Potência no atraso zero reconstruída a partir dos pixels | 25 ensaios, 200 escalas; erro absoluto máximo de 3,33 × 10⁻¹⁴ |
| Permutação exata largura/churn–anular | 17 alocações tão extremas em 92.378: p = 0,0001840264998159735 |
| Correções Holm das duas famílias de 504 testes | Reproduzidas; erro máximo de 3,69 × 10⁻¹⁴ |
| Combinação Pearson + largura já publicada | Reproduzidos 18 acertos em 19 e concordância balanceada 0,95 |

A reconstrução da potência não importou funções dos programas auditados: releu os panoramas, extraiu as linhas 18–71, refez a remoção de tendência por mínimos quadrados, calculou a autocorrelação linear por `numpy.correlate` e somou diretamente a Morlet no atraso zero. Essa contraprova é distinta da convolução FFT de `core8.py`.

As 756 silhuetas foram reconstruídas a partir das matrizes arquivadas. A verificação da potência cobriu todos os ensaios; ela não equivale a recalcular todos os 900 mapas cruzados das três posições. A auditoria original já conferiu seis mapas completos contra seu motor novo e os 300 resumos centrais contra a fonte. Não se atribui à presente auditoria uma reexecução desses 900 mapas.

Evidências executáveis: [script independente](../scripts/audit_estatistica8.py), [verificação numérica](../resultados/estatistica8/verificacao_independente.json), [registro de execução](../logs/audit_estatistica8.log) e [protocolo das sensibilidades](../resultados/estatistica8/protocolo_sensibilidades.json).

## Fórmulas e processamento

O preparo `x = detrend(255-I)/std` de `core8.py:63–67` está coerente com o sinal usado nas âncoras. Subtrair novamente a média após a remoção de tendência é redundante até arredondamento, mas não é um erro. A autocorrelação usa divisor fixo N (`core8.py:96–98`), preservando a estrutura da estimativa enviesada; não usa o divisor variável N−|lag|.

A Morlet usa ω₀ = 6 e a relação λ/a = 4π/(6+√38), com λ e a em pixels/amostras. O módulo da CWT da autocorrelação no atraso zero, e não seu quadrado, fornece P. O cálculo cruzado `|W_Rxy|²/(Px Py)` é consistente com a formulação alternativa (`core8.py:79–109`). Os 25 perfis P foram recuperados pela soma direta independente.

Os pesos trapezoidais `h` integram em log λ. A distribuição usada na largura é `p_k = P_k h_k / Σ(P_j h_j)`, com L = log λ₇₅ − log λ₂₅. Isso define uma distribuição descritiva da WPS alternativa por log comprimento de onda; não é automaticamente uma distribuição de energia física. O relatório já faz essa distinção (`RELATORIO.txt:46–49`). Não se deve trocar pesos ou interpretação física sem declarar uma nova definição.

Hellinger e Wasserstein foram aplicadas às massas normalizadas. Wasserstein integra a diferença das CDF discretas em log λ, não em λ linear (`core8.py:144–166`). Pearson usa os resíduos de log P após retirar intercepto e inclinação contra log λ; sua distância euclidiana de vetores unitários satisfaz d² = 2(1−r). A igualdade foi conferida numericamente.

Os resumos WCO não são intercambiáveis: máximo integral, média, média central, percentil, área acima de 0,5 e máximos por escala respondem a perguntas diferentes. O piso de 5% do resumo `scale_peak_power05` é um piso relativo de potência, não significância. A diagonal dos kernels de área emprega autocoerência, não uma constante arbitrária. A rejeição dos dois kernels de área com autovalores negativos é coerente com o requisito declarado de geometria euclidiana; não foi encontrada reparação oculta de matriz para favorecer o resultado.

## Achados e limites

### E8-01 — Alta relevância — A análise não corresponde à execução descrita atualmente no TCC

**Evidência:** `estatistica8/protocolo.json`, `scripts/core8.py:12,33`, `scripts/extract.py:16–21`; `TCC/metodologia.tex:3–9,44–62,111–114,178–180,243–255,298–305`.

| Elemento | estatistica8 | TCC atual |
| --- | --- | --- |
| Fonte principal | v16.3_ancoras | v16_spatial2_18_71 |
| Aquisições e classes | 25; quatro classes | 26; cinco classes, com uma bolha dispersa |
| Pares únicos | 300 | 325 |
| Perfil de Pearson | Janela 2048 px centrada na âncora | Panorama completo de cada ensaio |
| Perfil dos pares WCO | Mesma janela 2048 px por aquisição | Prefixos com o menor comprimento de cada par |
| Faixa λ | 50–1024 px | 50–3000 px |
| Estatística | Silhueta, vizinhos, agrupamento, permutações | Descrição de semelhança e contraste intra/inter; sem teste de significância |

**Impacto:** 95% de concordância, a largura e os p-valores de estatistica8 não validam as tabelas atuais do TCC. Nem o controle denominado panorama completo em estatistica8 recompõe a análise do TCC: ele conserva 25 ensaios e 50–1024 px.

**Contraprova/limite:** o TCC identifica explicitamente sua pasta e não incorpora os números de estatistica8. Portanto, a diferença entre versões é confirmada, mas não constitui contradição textual já cometida pelo TCC. É uma fronteira que precisa ser preservada em qualquer conclusão conjunta.

**Recomendação:** comparar versões por uma matriz de fatores explícitos — amostra, posição, comprimento, faixa λ, resumo — alterando um fator por vez. A auditoria atual não modifica o TCC.

### E8-02 — Alta relevância — A seleção herdada consulta a família do padrão antes de gerar os descritores

**Evidência:** `v16.3_ancoras/scripts/select_anchors.py:262–265`; `scripts/README.txt:65–89`; `estatistica8/scripts/extract.py:16–35`; `estatistica8/scripts/analyze.py:8–18,123–148`.

Plug/slug recebem eventos delimitados com aprovação morfológica anterior; churn/anular recebem vales selecionados por proeminência e estabilidade. Depois os descritores resultantes ficam congelados enquanto se permutam rótulos.

**Impacto:** a permutabilidade de rótulos das quatro classes, ou de contrastes que atravessem essas duas famílias, não decorre apenas de permutar aquisições em vez de pares. Sob uma hipótese nula sobre os sinais originais, a extração dependente do grupo pode produzir distribuições diferentes por construção. Os p-valores descrevem a associação dos descritores já selecionados sob a hipótese adicional de sua permutabilidade; não são uma validação integral do processo de seleção.

**Contraprova/limite:** a seleção ocorreu antes dos novos valores WPS/WCO, e não há evidência de escolha das âncoras pelo maior escore de separação. Churn e anular usam a mesma regra geral de vales; a objeção não invalida automaticamente o contraste principal entre eles. Tampouco foi alegada significância corrigida em estatistica8.

**Recomendação:** para estudar separabilidade entre todas as classes, usar extração comum sem consultar rótulos: posições espaciais predefinidas, todas as janelas elegíveis ou resumo por aquisição. Manter as âncoras como análise de morfologia condicionada ao evento.

### E8-03 — Média relevância — Perturbações ±64 px cobrem apenas estabilidade muito local

**Evidência:** `estatistica8/scripts/extract.py:89–110`; `RELATORIO.txt:78–93`; novas tabelas [sensibilidades](../resultados/estatistica8/sensibilidades_descritivas.csv) e [janelas](../resultados/estatistica8/larguras_janelas_distribuidas.csv).

Cada janela deslocada compartilha 1984/2048 = **96,875%** dos pixels com a âncora original. Isso testa pequenas mudanças de fronteira, não representatividade espacial do panorama. O relatório original corretamente chama as janelas de sensibilidades dependentes.

Na contraprova v17, cinco inícios foram fixados por `round(linspace(0,N_panorama−2048,5))`, independentemente do rótulo e do conteúdo. Os quantis foram calculados pela integração do interpolante linear da WPS em log λ. Para churn/anular:

| Condição | Silhueta macro | Menor silhueta média de classe | Delta de Cliff |
| --- | ---: | ---: | ---: |
| Âncora original, quantil discreto | 0,4348 | 0,4151 | 0,9556 |
| Âncora, quantil contínuo | 0,4303 | 0,4079 | 0,9556 |
| Posição uniforme 0, prefixo | 0,1363 | 0,0465 | 0,6000 |
| Posição uniforme 1 | 0,2089 | −0,0559 | 0,5111 |
| Posição uniforme 2, centro | 0,2001 | 0,1925 | 0,6889 |
| Posição uniforme 3 | 0,1445 | 0,0055 | 0,6000 |
| Posição uniforme 4, sufixo | 0,1782 | 0,1380 | 0,7111 |
| Mediana das cinco larguras por aquisição | 0,2648 | 0,1182 | 0,7333 |

**Impacto:** o contraste não desaparece fora da âncora, mas a separação mais forte observada na âncora não caracteriza uniformemente as outras regiões. A mediana das cinco janelas preserva a ordem das medianas de classe e ambas as silhuetas médias positivas, com efeito menor.

**Contraprova/limite:** são cinco posições dos mesmos panoramas; algumas janelas se sobrepõem, especialmente nos registros curtos. Elas não fornecem 125 réplicas, nem garantem amostragem espacial suficiente. O resultado não distingue heterogeneidade física, iluminação e montagem; não demonstra intenção de maximizar o contraste na seleção original.

**Recomendação:** priorizar um descritor agregado por aquisição e apresentar a variabilidade dentro da aquisição. Se a âncora for mantida, interpretar o resultado como propriedade do contexto escolhido, não de todo o padrão.

### E8-04 — Alta relevância — Os controles separados m e q não identificam um efeito independente das condições operacionais

**Evidência:** `scripts/followup.py:13–29,41–63,76–79`; `RELATORIO.txt:84–99,123–126`; [suporte operacional](../resultados/estatistica8/suporte_operacional.csv).

Cada um dos 25 ensaios ocupa uma combinação m×q distinta. Uma permutação exata que preserve simultaneamente os blocos conjuntos m×q teria somente a alocação original. Os controles existentes permutam dentro de m **ou** dentro de q. Eles não preservam simultaneamente ambos, e exclusão de um código no teste de vizinhos é uma operação diferente de controle causal de confusão.

**Impacto:** o contraste largura/padrão continua compatível com influência de condições operacionais e de variáveis não controladas. Usar muitos pares, mais permutações ou uma correção diferente não resolve a falta de sobreposição experimental. O código q é um indicador de condição e não substitui sozinho os valores físicos medidos.

**Contraprova/limite:** o relatório reconhece a confusão e não afirma causalidade. Há 200 alocações dentro de q e 1728 dentro de m para churn/anular; os testes correspondentes estão implementados de forma coerente para esse contraste binário. A inexistência de réplicas m×q não elimina o valor descritivo da largura, mas impede uma conclusão isolada sobre padrão.

**Recomendação:** obter repetições e sobreposição entre padrões nas condições físicas relevantes. Até lá, relatar contraste condicionado à amostra, comparações dentro de condições com suporte e mapas das covariáveis, sem apresentar ajuste estatístico como substituto de desenho experimental.

### E8-05 — Média relevância — O teste denominado ganho usa uma hipótese nula diferente de igualdade de desempenho

**Evidência:** `scripts/analyze.py:129–148`; `RELATORIO.txt:115–121`.

A mesma permutação de rótulos gera a diferença de concordância entre candidato/combinação e Pearson. Isso produz a referência da diferença sob permutabilidade global dos rótulos. Não preserva um possível sinal verdadeiro de Pearson enquanto testa apenas a contribuição adicional da largura; portanto não é, em geral, teste de igualdade de riscos entre dois procedimentos que já predizem um padrão.

**Impacto:** p bruto de 0,0312 não sustenta, por si, superioridade incremental sobre Pearson. O valor corrigido foi 1; não há decisão de superioridade confirmada.

**Contraprova/limite:** o relatório original explicita a referência de permutação, rejeita a leitura confirmatória e não confunde p alto com inexistência do efeito. Não há erro no cálculo da distribuição nula escolhida; a restrição é interpretativa.

**Recomendação:** no caminho sem aprendizado de máquina solicitado agora, abandonar acurácia/ganho de vizinhos como objetivo e comparar dimensões descritivas complementares: forma residual, largura e variabilidade espacial. Uma eventual validação futura de previsão seria outra tarefa.

### E8-06 — Média relevância — “Potência média” possui domínios diferentes em estatistica8 e no TCC

**Evidência:** `estatistica8/scripts/core8.py:93–105`; `TCC/metodologia.tex:237–255`.

Em estatistica8, `mp` é média de |W_Rxx| somente onde `N−1−|lag| ≥ √2 a`. O domínio diminui com a escala. No TCC, M é média de todos os 2N−1 atrasos e os mapas são explicitamente sem máscara geométrica. Os dois perfis podem ter inclinações, quantis e picos diferentes mesmo para o mesmo sinal.

**Impacto:** as sete representações `mean_power_*` não validam diretamente os perfis médios laterais descritos no TCC. O nome abreviado não revela a diferença de domínio.

**Contraprova/limite:** a WPS no atraso zero usada no principal achado de largura não emprega essa média; o problema não invalida sua reconstrução. A máscara é uma definição deliberada herdada da fonte, não erro de índice.

**Recomendação:** qualificar o perfil como média no suporte geométrico e registrar explicitamente o domínio. Comparar a versão integral apenas como sensibilidade identificada, com mesma janela e mesma faixa λ.

### E8-07 — Baixa relevância — A flag robust_pass não exige silhueta positiva em todas as classes após omissão

**Evidência:** `scripts/analyze.py:159–166` e `dados/decisoes.csv`.

`robust_pass` exige o mínimo da silhueta **macro** após omitir um ensaio, mas não `min_class_loo > 0`. Quatro linhas marcadas robustas têm uma classe negativa nessa análise: `wps_slope` e `mean_power_slope` em churn/plug, e as mesmas métricas em churn/slug. Os mínimos de classe são aproximadamente −0,0456, −0,0192, −0,0463 e −0,1283.

**Impacto:** a flag isolada não significa robustez uniforme por classe. A largura churn/anular não sofre esse problema: mínimo macro 0,4058 e mínimo por classe 0,3698.

**Contraprova/limite:** os mínimos estão gravados ao lado da flag, e o protocolo não estabelece explicitamente que todos devam permanecer positivos após a omissão. Trata-se de um critério mais fraco que a palavra “robusto” pode sugerir, não de resultado oculto.

**Recomendação:** usar a expressão “passou o critério operacional declarado” ou impor o critério por classe, se essa for a pergunta científica. Não mudar o critério retrospectivamente para selecionar apenas resultados favoráveis.

## Hipóteses adversariais que não se sustentaram

1. **“O efeito é um artefato dos quantis discretos.”** Não sustentado. Largura contínua manteve delta de Cliff 0,9556 e silhueta 0,4303. A resolução dos quantis originais é aproximadamente 0,01517 em log λ.
2. **“200 escalas são insuficientes para calcular a largura.”** Não sustentado neste intervalo. Passar para 800 escalas produziu silhueta 0,4302; as medianas contínuas ficaram praticamente iguais. Isso não transforma escalas em réplicas independentes.
3. **“Toda a separação vem do limite 1024 px.”** Não sustentado. Os tetos disponíveis de 509,54 e 767,53 px mantiveram silhuetas 0,4497 e 0,5141. Esses valores não autorizam selecionar o melhor teto depois da inspeção; apenas mostram que a borda superior não é explicação suficiente.
4. **“Os 300 pares foram usados como 300 observações independentes.”** Não sustentado. As permutações movem rótulos de aquisições sobre a matriz; a unidade inferencial foi mantida. Os 504 testes por família e a correção Holm foram reconstruídos.
5. **“95% foi um erro de contagem.”** Não sustentado. São 9/10 churn e 9/9 anular, cuja média é 95%; o total é 18/19. A seleção do candidato após exploração continua impedindo tratar esse número como estimativa imparcial de desempenho futuro.
6. **“Nenhum resultado significativo prova ausência de separação.”** Não sustentado e não afirmado no relatório original. A largura possui efeito descritivo; os testes corrigidos não confirmaram uma hipótese. Plug/slug tem apenas 3+3 ensaios, e a simetria das estatísticas de separação impõe p mínimo de 2/20 = 0,10.

## Encaminhamento sem aprendizado de máquina

O caminho mais defensável é descrever aquisições com medidas complementares e examinar sua estabilidade espacial: largura do espectro, forma residual de Pearson e variabilidade dentro da aquisição. A mediana de larguras em janelas distribuídas é candidata concreta, já calculada nesta auditoria, mas apresenta separação apenas parcial e requer avaliação com parâmetros congelados em novos ensaios.

Vizinhos 1-NN com rótulo previsto e agrupamento automático já fazem parte de estatistica8; a ausência de scikit-learn não os torna análises puramente descritivas. Nesta auditoria, o resultado 1-NN existente foi apenas conferido. As alternativas novas usam silhuetas com grupos fixados pelos rótulos experimentais, medidas de efeito e resumos por aquisição; nenhum classificador novo foi ajustado.

Prioridades: (1) usar uma regra comum de extração espacial; (2) separar variação intra-aquisição da variação entre aquisições; (3) harmonizar amostra, faixa λ e definição do perfil antes de comparar versões; (4) manter as condições operacionais visíveis; (5) ampliar replicação antes de alegar generalização. Alterar a correção múltipla para obter p < 0,05 ou multiplicar janelas como réplicas não resolve as limitações observadas.

## Artefatos

- `resultados/estatistica8/verificacao_independente.json`: precisão numérica, hashes, versões e contagens.
- `resultados/estatistica8/reproducao_silhuetas.csv`: 756 conferências.
- `resultados/estatistica8/nulo_exato_largura_independente.npy`: enumeração independente da hipótese original.
- `resultados/estatistica8/wps_soma_direta.npz`: perfis de 200 e 800 escalas, reconstruídos dos pixels.
- `resultados/estatistica8/larguras_por_ensaio.csv`: larguras originais, interpoladas, truncadas e dispersão espacial.
- `resultados/estatistica8/larguras_janelas_distribuidas.csv`: 125 janelas, identificadas como observações dependentes de 25 aquisições.
- `resultados/estatistica8/sensibilidades_descritivas.csv`: resultados completos de todas as configurações novas, inclusive desfavoráveis.
- `resultados/estatistica8/suporte_operacional.csv`: códigos m/q por ensaio.

Reprodução a partir da raiz do repositório: `PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 src2/v17/scripts/audit_estatistica8.py`.
