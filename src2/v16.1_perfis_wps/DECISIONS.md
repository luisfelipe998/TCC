# Registro de decisões

- 24/10/2026, antes da avaliação: a execução foi autorizada pelo pedido do usuário; a restrição histórica de somente planejar em PLAN.md deixou de se aplicar.
- 24/10/2026, antes da avaliação: por instrução explícita posterior, a base passou a ser v16_spatial2_18_71. A extração principal foi 18–71. A análise original 20–59 ficou como referência histórica.
- 24/10/2026, antes da avaliação: o v16 atual já havia retirado máscaras e colunas de médias restritas. A proposta antiga citava a versão anterior. Não se misturaram seus escores; médias restritas ausentes seriam recalculadas nesta pasta como sensibilidades explicitamente novas.
- 24/10/2026, antes da avaliação: a faixa principal foi definida pela COI geométrica no lag zero do menor registro, 50–2.295,940454 px. A faixa integral e a regra conservadora foram sensibilidades fixadas. A referência WCO herdada manteve seu domínio integral.
- 24/10/2026, antes da avaliação: concordância de ordenação, duas ponderações e retirada de aquisições foram fixadas no protocolo. Não foram usados valores-p ou intervalos de confiança sem pressupostos experimentais demonstrados.
- 24/10/2026, antes da avaliação: o resumo e todos os novos gráficos usariam PT-BR, inclusive bolhas dispersas, bolhas alongadas (plug), golfadas (slug), agitado (churn) e anular ondulado. Os identificadores de origem seriam preservados nas tabelas de procedência.
- 24/10/2026, antes da avaliação: a coluna central da CWT poderia ser calculada pela soma inversa de Fourier, conservando preenchimento, correlação e Morlet do núcleo. A equivalência numérica seria testada em três comprimentos antes do lote.

## Registro posterior aos resultados

- Nenhuma representação, banda, exclusão ou métrica principal foi alterada após observar a separação. A análise manteve a condição de avaliação interna e exploratória.
- Foram selecionados três casos interpretativos: menor H entre padrões, maior H intrapadrão e menor soma dos postos crescentes de H e do escore WCO entre pares intrapadrão. Os postos foram calculados somente para escolher ilustrações; não alteraram nenhum desfecho agregado. Os IDs e valores estão em results/selected_cases.csv.
- Os segmentos de comprimento fixo foram comparados ao registro inteiro para avaliar estabilidade. A WCO não foi recalculada nesses segmentos fixos; a comparação controlada entre representações usou os prefixos específicos dos pares e a faixa integral comum de rótulos de comprimento de onda.
- A auditoria independente confirmou as 325 distâncias principais, as 21 ordenações, a normalização e a preservação das fontes. Correções de apresentação ajustaram rótulos de eixos, posição de legendas e espaço da classe sem pares; os resultados numéricos não foram alterados.
- A pasta recebeu o nome v16_perfis_wps após a análise. O link simbólico v16.1 preservou as referências históricas sem modificar os documentos das bases. A conferência foi repetida no caminho definitivo.

- Solicitação posterior de renomeação: a pasta passou a v16.1_perfis_wps e o link simbólico v16.1 foi removido. Os nomes antigos não foram mantidos como atalhos. Somente caminhos, documentação e registros de integridade foram atualizados; cálculos e figuras foram preservados.

## Revisão didática de 26/10/2026

- A pedido do usuário, o resumo executivo foi reorganizado como leitura guiada para leitores sem contato prévio com a análise. A sequência passou a explicitar aquisição, sinal, perfil, distância entre pares e concordância A, antes da comparação de versões e das limitações.
- Foram acrescentadas quatro figuras derivadas das tabelas existentes: três exemplos reais de pares, distribuição das distâncias, contagem das ordenações e comparação essencial dos métodos. Os casos ilustrados permaneceram os três já selecionados em `selected_cases.csv`, sem nova escolha para favorecer o resultado.
- O histograma utilizou faixas iguais de largura 0,025 para os dois grupos, com todas as observações e porcentagens calculadas separadamente em cada grupo. As faixas serviram somente à apresentação; não entraram no cálculo de H ou A.
- A contagem explícita de 15.839 comparações favoráveis em 20.706, sem empates, reproduziu o A principal já registrado. Essas comparações reutilizaram aquisições e não constituíram novas réplicas experimentais.
- O guia incorporou um glossário, a leitura dos gráficos técnicos e uma distinção explícita entre distância, escore WCO, ordenação e acurácia. O HTML passou a apresentar navegação por etapas e complementos expansíveis.
- A revisão alterou a apresentação e seus geradores. O protocolo, a configuração, os sinais, os perfis, as distâncias e os resultados da avaliação permaneceram os mesmos. As novas figuras e suas fontes foram registradas em `results/reading_guide.json`.
