# Reconciliação da auditoria com o texto atual do TCC

A releitura do estado atual de `TCC/` não identificou mudança das fórmulas executadas, dos parâmetros de processamento ou do alcance descritivo que exija retirar os achados da auditoria v17. A correspondência entre a execução v16 e a metodologia permanece válida. As explicitações bibliográficas atuais tornam algumas distinções mais precisas, sem alterar os números auditados.

Esta conferência foi motivada por divergências entre os hashes atuais e os hashes registrados no início da tarefa. A comparação semântica foi feita com os trechos anteriormente lidos na auditoria, não com um snapshot textual completo da versão anterior. Portanto, este registro não é um diff exaustivo. Os hashes identificam os estados; não identificam quem ou qual processo os modificou.

## O que foi relido

Foram relidos integralmente `introducao.tex`, `metodologia.tex`, `fundamentacao.tex`, `discussao.tex` e `referencias.bib`, além da síntese `v17/findings/RELATORIO_CONSOLIDADO.md`. As fórmulas e decisões foram reconfrontadas com o núcleo v16 já auditado e com o extrato local da página 414 de Schuck e Bodmann. `resultados.tex`, `conclusao.tex` e os demais arquivos TeX da raiz de `TCC/` tiveram os hashes novamente registrados; esta etapa não substitui a reconferência de tabelas feita pela frente principal. Arquivos compilados em `dist/` não foram usados como fonte para julgar as fórmulas.

## Pontos reconciliados

| Ponto atual | Localização na leitura atual | Consequência para a auditoria |
|---|---|---|
| Divisor enviesado N e distinção do divisor original N−m, apresentado para atrasos não negativos, com extensão bilateral N−|m| | `metodologia.tex:167–174` | É uma explicitação mais fiel do domínio das equações da fonte. A implementação continua usando N; não muda a reprodução numérica nem a avaliação da adaptação. |
| WPS como módulo, WCS complexa da definição 32.16 e WCO com módulo quadrático no numerador | `fundamentacao.tex:644–677`, rótulos `eq:wps-alt`, `eq:wcs-alt`, `eq:wco-alt` | As três fórmulas permanecem compatíveis com o código. O texto agora distingue explicitamente o estimador 32.20, que toma módulo antes da média, da definição complexa 32.16. Com uma realização, isso não muda a WCO; modifica a retenção da fase na WCS. Essa qualificação deve ser mantida ao falar em fidelidade à fonte. |
| Pseudorreplicação explicitada e apoiada em Hurlbert | `metodologia.tex:401–412`; entrada `hurlbert1984` em `referencias.bib` | Reforça o alcance já reconhecido: 325 pares e 200 escalas não são réplicas independentes. Não foram introduzidos valores-p ou intervalos inferenciais na análise v16. |
| Objetivo de comparação descritiva, sem treinamento e sem equivalência entre correlação alta e identidade do padrão | `introducao.tex:43–60` | Compatível com a auditoria e com a exploração v17. A menção a classificadores na revisão não representa execução nova de aprendizado de máquina. |
| Antecedentes metodológicos distinguem combinações de descritores e variáveis físicas, sensores temporais e classificação de imagens | `discussao.tex:215–237`; revisão bibliográfica em `fundamentacao.tex` | Não altera a extração dos panoramas, a atribuição dos rótulos ou a faixa de análise. A revisão não deve ser interpretada como validação da fração de vazio ou da separabilidade nos dados locais. |
| Dados e processamento: 26 ensaios, 54 linhas 18–71, panoramas completos para Pearson, prefixos comuns para WCO, Morlet ω₀=6, 200 valores λ entre 50 e 3.000, piso principal 5%, remoção de reta em log-potência | `metodologia.tex:44–69,84–140,180–342` | Permanecem os mesmos objetos matemáticos e domínios que foram reproduzidos na auditoria. Nenhum novo cálculo é necessário por causa dessas alterações de texto. |
| Interpretação e limitações: amplitudes e inclinações removidas, lags distintos por escala, condições operacionais, menor registro com cerca de 1,06 ciclos e grupos pequenos | `discussao.tex:12–59,253–300` | Permanecem compatíveis com W16-01 a W16-05 e com a contraprova da composição. A discussão ainda não quantifica a sensibilidade ao singleton calculada na v17. |

A distinção sobre a Equação 32.20 é substantiva para a descrição bibliográfica do campo cruzado, mas não refuta os resultados da v17: nossa verificação usa o campo complexo realmente armazenado e o numerador de WCO é seu módulo quadrático. Para múltiplas realizações, média de módulos e módulo da média não são em geral intercambiáveis; a auditoria e o código examinados limitaram-se a uma realização.

A inclusão de novas linhas deslocou referências por número de linha nos relatórios iniciais. As referências ao estado inicialmente examinado não devem ser lidas como um mapeamento exato para a versão atual. Os rótulos LaTeX e as localizações atuais acima permitem reencontrar as passagens; não foi necessário reescrever os achados científicos.

## Consequência para o relatório consolidado

A síntese científica continua sustentada: os cálculos verificados correspondem à execução descrita, e as restrições principais envolvem suporte, composição da amostra, representatividade das janelas, interpretação óptica e condições operacionais. Os controles novos continuam exploratórios e não comprovam separabilidade geral, predição ou medida física de vazio.

Foi identificada uma correção necessária na afirmação de preservação dos arquivos: expressões como **“os arquivos originais permaneceram como estavam”** eram amplas demais quando os hashes de arquivos do TCC mudaram durante a tarefa. A frente principal corrigiu a síntese, e a releitura confirmou que ela agora distingue as escritas restritas à v17 das alterações concorrentes do TCC. A formulação verificável é: **“Esta tarefa gravou somente em v17. As três pastas de análises auditadas conservaram seus hashes iniciais; houve alterações concorrentes nos arquivos do TCC, registradas e reconciliadas, sem atribuição de autoria.”** A verificação global dos hashes cabe ao manifesto da entrega. Esta reconciliação não alterou o relatório consolidado nem nenhum arquivo do TCC.

## Identificação do estado relido

Horário UTC do registro: **2026-10-27T04:28:31+00:00**. Os hashes a seguir foram calculados após a releitura e reconferidos imediatamente antes de salvar este documento. Para os arquivos de texto especificados, houve consistência entre essas duas leituras. Mudanças posteriores exigem nova comparação com estes hashes.

| Arquivo | SHA-256 |
|---|---|
| `TCC/anexos.tex` | `f1b839f0924914cbff14442b048083a34d878186a0ce389a4a89ccee047edaa6` |
| `TCC/apendices.tex` | `c12438d4f16784f18b374ecc0a5b9c8a9f50d640f0d8b36d25be397dc8a187a3` |
| `TCC/conclusao.tex` | `cf21af9651a24229914950265fe3bc03513605aafe70d207703f2668f8fd5c1d` |
| `TCC/discussao.tex` | `4c8cd684e8d4b6d9c03986ed4cabad1dbb11a03de9b1ce96d103f8d66b45e08f` |
| `TCC/fundamentacao.tex` | `4b94954398d0111e700a21d7e95664f4527b6485e4c3931954125e51c6cc2e13` |
| `TCC/introducao.tex` | `9aaee03aea08dc034c1999031bf4fd6eed6760677df8b44f8150197af6f7cf69` |
| `TCC/main.tex` | `7ac3c7779da8aafd1f1df6e2404fe0c7cf341204176689e9bcaa22299584c59f` |
| `TCC/metodologia.tex` | `a578ad88793d670c31a26db34aad86b0e8afe2216c885256e8e0100781db6b4a` |
| `TCC/pretextual.tex` | `1a62562744bffe1c041a73b212cc3429d05ce16769b878e777d0449404c81dd4` |
| `TCC/referencias.bib` | `b0ff09a79ee8627abbf6967dc267413623da36d8ecbc9001cbca693a9374ecc5` |
| `TCC/resultados.tex` | `587848799cc164cf441dbaf0d85dba05c13114a66ca08f94d257c82022a0e01b` |

A comparação com `v17/resultados/fontes_antes.json` confirmou divergência de hash em `introducao.tex`, `fundamentacao.tex`, `metodologia.tex`, `discussao.tex` e `referencias.bib`. Os demais arquivos TeX da raiz de `TCC/` mantiveram seus hashes nessa comparação. Esse resultado identifica mudança de conteúdo, sem estabelecer sua autoria ou descrever todas as alterações textuais.

A frente principal preservou cópias das 11 fontes finais em [texto_conferido](../resultados/comparacao_tcc/texto_conferido), com o [manifesto do estado final](../resultados/comparacao_tcc/snapshot_texto_final.json). As 11 cópias foram comparadas com os arquivos atuais: todos os hashes coincidiram com os valores desta tabela. Essas cópias documentam o estado final relido; não reconstituem um snapshot textual completo do início da auditoria.

Conclusão: **achados científicos mantidos; atribuição bibliográfica do estimador cruzado melhor delimitada; síntese corrigida para distinguir a atuação desta tarefa das alterações concorrentes.**
