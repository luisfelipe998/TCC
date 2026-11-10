# Auditoria científica consolidada das análises de escoamento

27 de outubro de 2026. Três avaliações paralelas foram formadas independentemente, uma por pasta, e depois submetidas a contraprovas cruzadas. A comparação com o TCC e a exploração de novos descritores foram conduzidas separadamente. Todas as escritas desta tarefa ficaram em `src2/v17`; as três pastas de análise conferem com os hashes iniciais. O TCC recebeu alterações concorrentes externas a esta tarefa: seu estado final foi reconferido e está identificado na [reconciliação do texto](reconciliacao_texto_atual.md). As citações de linha dos relatórios independentes referem-se às leituras iniciais; a reconciliação fornece as referências atuais.

**Os cálculos principais verificados estão corretos e o TCC atual descreve fielmente a execução v16. O problema predominante é a interpretação e a comparabilidade das medidas, não um erro algébrico geral.** Há evidência descritiva de diferenças entre alguns padrões, especialmente churn e anular, mas os resultados dependem da composição da amostra, dos trechos observados, do suporte espectral e das condições ópticas/operacionais. Nenhuma análise nova treinou classificador.

Uma alternativa concreta, sem aprendizado de máquina, produziu resultado promissor: a largura de uma distribuição espectral de Welch agregada em três janelas fixas apresentou silhueta macro **0,5905** entre os 19 ensaios churn/anular. A reprodução independente confirmou o cálculo. A separação enfraquece e passa a apresentar alguma sobreposição ao mudar a resolução do estimador; continua sendo exploração da mesma amostra. Plug/slug e quatro padrões não apresentaram separação estável nas variantes novas.

## O que os dados permitiram confirmar

| Frente independente | Cobertura nova | Resultado |
| --- | --- | --- |
| estatistica8 | 479 hashes; 25 perfis WPS por soma direta; 756 configurações de silhueta; permutação exata da largura e famílias Holm | Reprodução numérica; p bruto 17/92378 e correção mantidos; 95% herdado conferido, sem interpretação preditiva nova |
| v14.1 | 13.026 JPEGs, 13.000 deslocamentos, 26 panoramas corrigidos, máscaras, área, flips e Otsu | Reconstrução exata e aritmética compatível; uma divergência de hash upstream, sem diferença numérica detectada |
| v16 | 31 entradas, 26 perfis, 325 resumos e três mapas completos; identidades matemáticas e 144 pares sintéticos | Fórmulas e valores reproduzidos; não foi alegada reexecução de todos os 325 mapas |
| TCC | 26 linhas de ensaios, 104 células das outras tabelas principais e 12 figuras | Transcrição e hashes conferem; não houve alteração do texto |
| Exploração v17 | 25 aquisições, duas bandas, três representações espectrais, janelas fixas e âncoras, 360 comparações descritivas | Todas as variantes salvas; nenhuma escolha de limiar, treinamento ou teste preditivo |
| Contraprovas cruzadas | Âncoras sob troca de rótulo, quantis, agregação por classes, Welch manual, Parseval e resolução espectral | Principais números confirmados; hipóteses excessivas foram refutadas ou limitadas |

Reproduzir código não calibra a medida física: a reprodução das máscaras v14.1 usa o detector original e comprova consistência, não uma identificação independente da fase. Os testes de wavelet comprovam operações e identidades, não estacionariedade das imagens ou correção hidrodinâmica dos deslocamentos.

## Achados com maior consequência científica

| Achado | Evidência decisiva e contraponto | Consequência |
| --- | --- | --- |
| **As pastas não são intercambiáveis** | TCC/v16: 26 ensaios e λ 50–3000; estatistica8: 25, âncoras 2048 e λ 50–1024; v14.1: geometria de 53 linhas, panorama corrigido e corte óptico | Harmonizar população, ROI, montagem, domínio e banda antes de comparar métodos. Os dois conjuntos estatísticos de 25 ensaios excluem registros diferentes. |
| **O singleton domina parte do contraste macro** | Sem a única bolha dispersa, Δ Pearson macro cai de 0,145300 para 0,016502; por par, de 0,065889 para 0,028735 | A direção permanece positiva. A queda limita a interpretação global; não é erro da média nem justificativa automática para excluir o registro. |
| **A WCO depende do suporte e da pergunta** | Ruídos independentes com N=3171 e λ=3000 deram mediana do máximo 0,9685. Nos pares reais, Spearman entre N comum e escore foi −0,7997 | Máximo alto não demonstra identidade do padrão. O controle não prova que o escore real inteiro seja artefato: escalas abaixo de três ciclos tiveram peso mediano zero. |
| **Coerência e semelhança estatística respondem a perguntas distintas** | Processos independentes podem ter o mesmo espectro e correlação cruzada populacional zero | Usar distâncias espectrais para semelhança entre realizações; reservar WCO para hipóteses de estrutura compartilhada e deslocamento. |
| **Estabilidade local não assegura representatividade do panorama** | Deslocar 2048 px em ±64 mantém 96,875% do sinal. Largura WPS em cinco posições comuns teve silhueta churn/anular 0,136–0,209, versus 0,435 na âncora | O efeito existe, mas não é uniforme. Quantis contínuos e grade mais densa preservaram o resultado: quantização não o explica. |
| **Extração herdada depende da família do rótulo** | Plug/slug e churn/anular usam regras diferentes de âncora. Trocar churn↔anular preservou as 19 âncoras na contraprova | A permutabilidade das quatro classes exige cautela; não há evidência de vício de seleção entre churn/anular nem de escolha intencional pelo melhor resultado. |
| **Ocupação e flips têm dependência estrutural** | Da mesma máscara, α≤q e F≤2W/(W−1)(1−α). Normalizar F pela esperança de embaralhamento com ocupação fixa inverteu sua ordenação churn/anular | F bruto continua separando descritivamente. Não deve ser interpretado como atividade de interfaces independente da ocupação. O embaralhamento não é um nulo físico validado. |
| **Junções e segmentação afetam os índices ópticos** | Junções: mediana 5,28% das comparações e 17,89% dos flips. Retirá-las preservou a ordenação entre churn/anular | Há influência de montagem, mas ela não explica sozinha a separação. Área circular continua condicional à geometria/óptica; filme líquido não foi calibrado. |
| **Confusão operacional não é resolvida por mais pares** | Cada aquisição ocupa uma combinação m×q diferente; há apenas 3 plug, 3 slug e uma bolha dispersa | Controles separados por m ou q não isolam simultaneamente operação, padrão e aquisição. Janelas e pares não acrescentam ensaios independentes. |
| **Manifesto histórico de código apresenta uma divergência** | 1/13.105 hashes de v14.1 diverge: `src2/pipeline_v11.py`; dados e reconstrução numérica conferem | Registrar nova procedência em execução futura. Não confundir deriva de código com corrupção demonstrada dos resultados. |

Os relatórios especializados distinguem achados novos, limitações já reconhecidas, diferenças de escopo e hipóteses refutadas. Não foi confirmado erro material de fator da Morlet, divisor da autocorrelação, sinal do lag, quadrado da WPS, área circular, denominador dos flips ou correção Holm nos cálculos auditados.

## Resultado novo que merece continuidade

A proposta calculada nesta v17 compara **distribuições espectrais de intensidade**, com extração igual para todos os rótulos. Cada aquisição fornece três janelas de 2048 px — início, centro e fim. Calcula-se a PSD de Welch e transforma-se sua densidade para log comprimento de onda por `f·PSD(f)`, com `f=1/λ`. Cada distribuição é normalizada na banda e as três recebem pesos iguais. A largura é o intervalo interquartil dessa mistura em log λ; não é média das larguras nem energia física do escoamento.

Na configuração definida no protocolo, λ 50–512 e segmentos Welch de 1024, as larguras foram:

| Grupo | Aquisições | Mínimo | Mediana | Máximo |
| --- | ---: | ---: | ---: | ---: |
| Churn | 10 | 0,9803 | 1,0823 | 1,2815 |
| Anular ondulado | 9 | 0,5296 | 0,7266 | 0,9083 |

Isso produziu silhueta macro 0,5905. O controle por retirada de uma aquisição variou entre 0,5701 e 0,6351. Em q=40, dois churn e três anulares preservaram a ordenação; a amostra é insuficiente para eliminar confusão conjunta ou validar generalização. Com segmentos de 2048, a silhueta caiu para 0,4252 e 88/90 comparações preservaram a ordem, com alguma sobreposição. O resultado forte não foi causado por erro de Jacobiano, mistura, quantil ou interpolação, mas depende de resolução/suavização.

Os 90 pares não são 90 réplicas, e a ordenação completa na configuração inicial **não é 100% de acurácia**. A configuração foi examinada junto de muitas alternativas nos mesmos ensaios. Nenhum p-valor confirmatório foi atribuído ao novo resultado. A largura Welch da mistura não separou quatro classes de forma estável: a silhueta foi −0,0421 nessa configuração. Plug/slug apresentou 0,0547 e mudou de sinal nas sensibilidades.

Detalhes, fórmulas, tabelas completas e contrapontos: [confronto com o TCC e exploração](confronto_tcc_e_exploracao.md) e [revisão adversarial de Welch](contraprova_exploracao.md).

## Caminho recomendado sem aprendizado de máquina

1. **Congelar a representação por aquisição.** Usar regra de janela que não dependa do rótulo, banda justificada por suporte e um resumo das distribuições, acompanhando a heterogeneidade entre posições. Manter a largura Welch como candidata exploratória, sem escolher novos parâmetros pelo melhor resultado.
2. **Comparar semelhança com distâncias de distribuições.** Usar largura, Hellinger/Wasserstein e forma residual como perguntas distintas; conservar inclinação e variabilidade quando forem relevantes. Pearson após remover tendência não pode informar sobre a inclinação removida.
3. **Investigar informação morfológica para plug/slug.** Comprimentos de intervalos, intermitência, assimetria radial e evolução temporal são candidatos plausíveis. Validar a máscara e a censura de eventos nas bordas antes de interpretar esses descritores. São próximos experimentos, não resultados demonstrados nesta auditoria.
4. **Usar quadros nativos como controle da montagem.** O controle temporal novo, 70 ms a cada 1 ms, mostrou grande sobreposição churn/anular. Isso delimita a diferença entre transições ópticas temporais e densidade espacial; não converte uma na outra.
5. **Ampliar repetição e condições comparáveis.** Replicar aquisições e documentar exposição, saturação, deslocamentos e condições físicas. Avaliar o protocolo congelado em novos registros. Repetir pares, escalas ou janelas sobrepostas não supre essa necessidade.

Uma conclusão defensável hoje é que **há descritores com contraste descritivo útil para um subconjunto dos padrões, mas ainda não uma demonstração geral de separabilidade ou de medida física de vazio**. Essa conclusão é compatível com o tom prudente do TCC reconferido. Esta tarefa não aplicou recomendações ao texto nem às análises anteriores; as mudanças concorrentes do TCC estão registradas separadamente, sem serem atribuídas a esta auditoria.

## Relatórios e evidências

- [Auditoria de estatistica8](auditoria_estatistica8.md) e [contraprova das âncoras e janelas](contraprova_estatistica8.md).
- [Auditoria óptica e geométrica de v14.1](auditoria_voidfraction.md).
- [Auditoria de fórmulas e suporte da v16](auditoria_wavelets.md) e [contraprova da composição da amostra](contraprova_composicao.md).
- [Confronto com o TCC e exploração sem classificador](confronto_tcc_e_exploracao.md) e [contraprova de Welch](contraprova_exploracao.md).
- [Protocolo, incluindo adendos exploratórios](../PROTOCOLO.txt).
- [Validação da entrega e preservação das fontes](../resultados/validacao_entrega.json).

`PASS` nos arquivos de verificação significa que os critérios computacionais indicados passaram. Não equivale a validação experimental, causal ou preditiva.
