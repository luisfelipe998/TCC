# Subtração condicionada à interface: protocolo da investigação

O experimento avaliou se a remoção de uma contribuição óptica associada à geometria instantânea da segmentação permitiria acompanhar texturas interiores diferentes. A variável medida foi deslocamento óptico axial condicionado a máscaras candidatas. A operação não certificou que a textura restante acompanhasse o material líquido ou gasoso.

Foram mantidos os 26 registros selecionados na v15.1 e os 501 quadros nativos de cada registro. O FPS foi lido do inventário previamente conferido contra os arquivos da câmera. O lúmen correspondeu às linhas 18:71. As máscaras candidatas de gás foram a interseção dos limiares 0,60, 0,70 e 0,80; as candidatas de líquido foram o complemento da união, limitado ao lúmen. Essa nomenclatura não implica que todo pixel não segmentado seja líquido.

Três tratamentos foram comparados com os mesmos critérios, máscaras e quadros:

1. Imagem original, sem subtração de fundo.
2. Subtração de um fundo fixo por pixel, calculado exclusivamente na outra metade temporal.
3. Subtração de uma imagem prevista a partir da linha radial e da distância com sinal à fronteira da máscara de limiar 0,70, calculada em cada quadro. A tabela de intensidades médias foi ajustada exclusivamente na outra metade temporal.

O terceiro modelo empregou 25 classes de distância, de largura 2 pixels, com truncamento em ±24 pixels. Em cada linha, classes com menos de 100 pixels de treinamento receberam interpolação das classes povoadas da mesma linha. A tabela foi suavizada somente ao longo da distância, com sigma de uma classe. O preditor não recebeu qualquer estimativa de velocidade. A fração de pixels de avaliação que precisou de classes interpoladas foi registrada. A distância à fronteira da segmentação pode incluir a parede; não é uma distância certificada à interface tridimensional.

Cada metade temporal foi dividida em três blocos; seis blocos no total, sem pares que atravessem suas fronteiras. Dentro de cada bloco, a correlação normalizada usou a interseção da máscara de origem com a máscara de destino deslocada, para atrasos de 1, 2 e 4 quadros. Foram testados deslocamentos de −240 a +240 pixels. O mesmo tratamento foi executado com erosões espaciais de raio 3 pixels e, como sensibilidade, 1 pixel.

Um bloco precisou ter pelo menos 20 quadros com 100 pixels de fase candidata, desvio-padrão de intensidade original e do sinal comparado de pelo menos uma unidade, pico de correlação de pelo menos 0,20, margem de pelo menos 0,02 para concorrentes a mais de 5 pixels e concordância entre atrasos dentro de max(2 pixels/quadro, 20% da magnitude do deslocamento de atraso 1). A sobreposição mínima foi max(1.000 pares de pixels, 25% da maior sobreposição). O pico precisou exceder em 0,05 o maior pico de três controles com ordem dos quadros embaralhada no bloco. Blocos com mais de 20% dos pixels candidatos em intensidade ≥250 foram rejeitados. Esses controles são diagnósticos; não constituem teste de hipótese com valor-p.

O suporte de um registro exigiu pelo menos cinco dos seis blocos. O par gás candidato/líquido candidato exigiu pelo menos cinco blocos aprovados simultaneamente para os dois sinais. As velocidades do par foram médias dos deslocamentos de atraso 1 nesses mesmos blocos, ponderadas pelo número de intervalos comparados e multiplicadas pelo FPS nativo. Assim, duas coordenadas nunca combinaram janelas temporais diferentes. A média representa os intervalos selecionados; não foi identificada como média material de todo o vídeo. Não se exigiu velocidade constante entre os blocos.

A razão com sinal foi calculada apenas para pares aprovados e denominador com magnitude superior a 2 pixels/quadro, para evitar a divisão por um movimento quase nulo em relação à tolerância operacional. A razão foi denominada razão de movimentos ópticos candidatos. Nenhum resultado recebeu unidade m/s, pois a escala axial não foi calibrada independentemente. A razão cancela uma escala comum dentro do registro, mas não certifica a identidade física dos dois sinais.

Os controles sintéticos incluíram texturas interiores independentes com velocidades conhecidas, iluminação móvel compartilhada, halo de interface móvel, interface deformável sem textura material e interiores constantes. Os dois controles de halo foram construídos a partir da mesma família de relação distância/intensidade utilizada pelo preditor e, portanto, favorecem sua especificação; a recuperação nesses controles não demonstra remoção de toda refração experimental. A iluminação móvel compartilhada foi mantida como controle adversarial de interpretação: se ela passar pelos critérios, um resíduo rastreável não pode, por si, ser chamado de velocidade material.

As opções e critérios foram fixados antes da execução dos 26 registros. Resultados rejeitados permaneceram disponíveis nas tabelas; não foram substituídos por outro método para completar um gráfico.
