# Protocolo da investigação adicional — 24/10/2026

Objetivo: testar se as 26 gravações incluídas podem fornecer duas coordenadas comparáveis em um único gráfico, distinguindo movimento óptico de velocidades materiais de gás e líquido. Preservam-se as exclusões anteriores. Nenhum resultado ausente será substituído por zero, por outro ensaio ou por uma relação de velocidades imposta. Os resultados anteriores não serão sobrescritos.

## Hipóteses a desafiar

1. **A exigência de velocidade constante durante toda a gravação pode ser excessiva.** Testar correlação em janelas temporais curtas e média das velocidades locais. Variabilidade temporal real não equivale, por si só, a erro do estimador. Entretanto, cada janela precisa apresentar evidência de transporte.
2. **Uma janela retangular grande pode perder textura em regiões finas ou descontínuas.** Testar correlação com máscaras de aparência no ponto de origem e no destino, com exclusão de interfaces e sensibilidade ao limiar/erosão. O complemento de uma máscara continua sem identificar líquido independentemente.
3. **Um pico de correlação pode selecionar apenas uma escala óptica.** Comparar com energia direcional no espectro espaço–tempo e testar bandas espectrais, reversão e ordem temporal aleatória.
4. **Ampliar a cobertura pode alterar a grandeza medida.** Examinar médias das duas paredes sem exigir igualdade entre suas velocidades, mantendo cada parede verificável. Essa média tem definição óptica; não será denominada média material do líquido.
5. **Segmentação e equações físicas podem parecer fornecer informação que não está presente.** Testar máscaras complementares, ponderação do mesmo sinal global e exemplos de não identificação cinemática; revisar inferências de conservação e hipóteses de fechamento.

## Ramo de correlação local no tempo

Definição fixa para todas as gravações: núcleo 35:54; faixas externas 18:29 e 60:71; mesma coordenada axial e processamento da referência v15.1. Perfis radiais médios vêm das 501 imagens nativas; será verificada sua correspondência com as imagens. A média temporal de cada posição axial é subtraída.

Referência: dez blocos contíguos, sem pares que atravessem suas fronteiras. Em cada região/bloco: pico da correlação normalizada média, atraso de um quadro. Verificações: pelo menos 60% dos pares com textura >=1 nível de cinza, pico >=0,30, pico maior que nove controles de permutação local em pelo menos 0,05, consistência do atraso de dois quadros, reversão temporal e duas posições axiais (160/320). Picos das verificações devem ser >=0,15 e não atingir o limite de busca. Tolerância de velocidade: max(2 px/quadro, 20% do módulo da estimativa). Estes são critérios empíricos de suporte, não testes de significância ou limites de erro físico.

Um par por gravação requer núcleo e ambas as paredes aceitos simultaneamente em pelo menos 60% dos blocos e distribuição por pelo menos quatro quintos temporais. As coordenadas usam apenas os blocos comuns, ponderados por seu número de intervalos internos; a coordenada de parede é a média das duas faixas de mesma espessura. A concordância entre as velocidades das duas paredes não é exigida: elas podem diferir, mas ambas precisam de suporte.

Sensibilidade predefinida: cinco e vinte blocos, perfis sem subtração temporal e filtragem axial de escala 8/32 px. A faixa entre alternativas não será apresentada como intervalo de confiança. Saturação será explicitamente registrada e testada com controles sintéticos de recorte de intensidade; não basta um número finito para aceitar uma estimativa. Não se alterará um limiar para obter 26 pontos aprovados.

## Validação e decisão

Os controles devem incluir deslocamentos conhecidos, ausência de textura, padrão estacionário, embaralhamento temporal, deformação ou múltiplas velocidades e sensibilidade ao processamento. Parâmetros específicos dos outros ramos serão registrados em seus scripts e relatórios antes dos resultados agregados.

Haverá três avaliações separadas: (a) reprodução aritmética e proveniência; (b) evidência de transporte óptico sob as condições testadas; (c) identificação material da fase. Passar por (a) e (b) não implica (c). Serão apresentados todos os cenários e motivos de falha. Se nenhum método identificar duas velocidades materiais em todos os cenários, essa limitação permanecerá no título/legenda/conclusão, mesmo que se obtenham 26 coordenadas ópticas finitas.

## Extensão exploratória após os primeiros resultados

Os resultados com dez blocos mostraram falhas concentradas em uma das paredes, enquanto a outra ainda continha transporte verificável. Será avaliada uma grandeza diferente: movimento do núcleo e de **uma parede selecionada por observabilidade**. A seleção entre as duas paredes usa somente os cinco primeiros blocos: maior número de blocos com suporte, desempate pelo pico médio de correlação. A escolha é congelada e avaliada nos últimos cinco blocos. Exigir pelo menos três blocos conjuntos núcleo/parede em cada metade. As médias reportadas para essa alternativa usam somente blocos aceitos da segunda metade; não representam a média do vídeo inteiro nem a média de todo o filme líquido.

Esta extensão surgiu após explorar os resultados e não constitui validação em um conjunto experimental novo. Serão preservados o caso sem subtração temporal, a seleção invertendo as metades e a comparação entre a parede escolhida e a outra. Não serão usados padrão, magnitude esperada de velocidade ou distância aos pontos de outro cenário para escolher a parede. A restrição inicial a duas paredes permanece como referência.

A auditoria dessa extensão identificou que o fundo original usava a média dos 501 quadros, compartilhando informação entre as metades. Será executado um controle corrigido: o fundo será estimado somente na metade de seleção e congelado para a metade de avaliação, inclusive nas variantes de filtragem 8/32 px. Os resultados anteriores permanecerão identificados como exploração com fundo compartilhado; nenhuma aprovação será transferida automaticamente à versão corrigida. Mesmo a separação corrigida não equivale a um conjunto experimental externo, pois o desenho do estudo explorou essas gravações.
