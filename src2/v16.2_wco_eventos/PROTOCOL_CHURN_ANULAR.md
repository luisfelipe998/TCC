# Extensão: âncoras ópticas em churn e anular ondulado

Este protocolo foi definido antes de calcular novas WCO da extensão. O usuário autorizou continuar na mesma pasta. Os resultados de plug/slug, o protocolo original, o catálogo e seu congelamento permanecem intactos. Esta etapa amplia o objeto para **segmentos centrados em marcos ópticos de textura**, sem exigir fronteiras de um evento físico isolado. A análise continua exploratória: os registros já foram examinados anteriormente.

## Escopo e piloto anterior à WCO

São mantidas as dez aquisições elegíveis de churn e nove de wavy annular, conforme a auditoria da base v16_spatial2_18_71: 45 pares intrapadrão de churn, 36 de anular e 90 entre os dois padrões, totalizando 171 pares de aquisições. Nenhuma exclusão original será revertida. O piloto de imagens usou m03q05_001, m05q40_001, m03q10_001 e m05q80_001, cobrindo imagens claras e escuras. Foram vistos panoramas, recortes e quadros originais desses quatro registros, e depois seus perfis ópticos, sem wavelets novas.

O piloto mostrou depressões locais de intensidade identificáveis mesmo nas texturas contínuas. Elas não foram interpretadas como bolhas isoladas, espessura física de filme, cristas ou identidade de ondas. A palavra âncora nesta extensão significa uma coordenada óptica operacional, cuja estabilidade será quantificada. A revisão posterior das outras aquisições verifica representações e artefatos; não cria um conjunto externo de validação.

## Detecção por imagem e distribuição espacial

A entrada é a média das linhas 18–71 inclusive. Somente para detectar extremos, suaviza-se esse perfil com gaussiana σ=16 px, modo reflect. A análise wavelet usa o perfil original, sem essa suavização adicional. São propostos mínimos de intensidade para a análise principal e máximos como controles de polaridade óptica.

Usa-se find_peaks sobre −I suavizada para mínimos e sobre +I para máximos. A proeminência mínima é max(3 níveis de cinza, 0,25×(P95−P5)), com percentis calculados no perfil suavizado. A proeminência é local, em janela de 1.025 px; a largura à meia proeminência deve estar entre 16 e 512 px. A distância mínima entre candidatos do mesmo tipo é 128 px. Esses parâmetros derivam da viabilidade óptica do piloto, não de WCO. A largura é uma medida do perfil suavizado, não um contorno físico nem o comprimento do recorte.

O centro precisa permitir [c−1024,c+1024), com 16 px adicionais nas duas bordas da aquisição. A região elegível é dividida em quatro estratos espaciais iguais. Em cada estrato e por tipo, seleciona-se a maior proeminência; empate favorece a menor coluna. A seleção é feita da esquerda para a direita, exigindo distância ≥512 px às âncoras já selecionadas do mesmo tipo. Se não houver candidato elegível, o estrato permanece vazio. Assim há no máximo quatro mínimos e quatro máximos por aquisição. Todos os candidatos e motivos de não seleção são preservados.

Registrar: ID, aquisição, padrão global, tipo, centro, proeminência, largura, limites interpolados da meia proeminência, estrato, limites do segmento, costuras, mistura na âncora, deslocamentos não positivos no contexto, quadro e coluna de origem, vizinhos e sobreposição. Os quadros de origem são procedência computacional; não demonstram identidade física. A distância a uma costura não é usada para eliminar casos inconvenientes, pois a montagem é composta de faixas em todo o panorama. O caso m03q40_001 permanece na análise, com seus artefatos assinalados.

## Pareamento e controles congelados

Os pares principais cruzam aquisições diferentes. Antes de calcular WCO, são criados três grupos:

- **Mínimos com larguras compatíveis:** razão maior/menor da largura à meia proeminência ≤1,50. Compatibilidade significa mesmo tipo óptico e esta tolerância de largura, sem equivalência física comprovada.
- **Mínimos com larguras distintas:** razão ≥2,00. A faixa intermediária fica registrada e fora desses contrastes.
- **Polaridades opostas:** mínimo × máximo, um par em cada orientação entre as aquisições quando disponível. Não se exige a mesma largura; esse controle altera polaridade e possivelmente forma, e não isola uma única causa.

Em cada combinação de aquisições/grupo, selecionam-se até dois pares pelos menores hashes SHA-256 da semente, grupo e IDs. Para polaridades opostas, seleciona-se um por orientação. Não se procuram parceiros pela maior coerência. A lista completa de candidatos e a seleção serão congeladas. A presença de diferenças de cobertura entre grupos será explicitada.

Os controles fixos de início, centro e fim têm os mesmos 2.048 pontos. Seus resultados principais para os 171 pares já existem na etapa 1 e serão reutilizados após conferir hashes e parâmetros. Para faixa 30–49, remoção somente da média e núcleo v5, esses controles serão recalculados no mesmo domínio. A comparação usa somente pares de aquisições presentes em ambos os grupos. Não se chama o recorte fixo de negativo físico garantido.

WCO usa módulo ao quadrado e é invariável a uma inversão global do sinal. Portanto, polaridades opostas podem ter WCO alta; esse controle não exige coerência baixa para validar o código. Serão reportados também self, dois mínimos da mesma aquisição quando disponíveis, atraso conhecido e controles sintéticos de mínimos compatíveis, larguras distintas, fundo comum e inversão de sinal. Comparações dentro de uma aquisição ficam fora da evidência entre aquisições.

## Wavelets, resumo e sensibilidades

Mantêm-se N=2.048, recorte antes de preparo, K=1, inversão de intensidade, remoção linear/média e padronização, estimador enviesado 1/N, Morlet ω0=6, 200 comprimentos de onda de 50–1.024 px e média em lags −32…+32. WPS alternativa é módulo; WCS é complexa; WCO usa denominador no lag zero. Não haverá clipping nem otimização de alinhamento. A média conservadora é um diagnóstico separado. Os mapas completos preservam as bordas; não se inferem milímetros, frequência temporal ou velocidade.

Todos os pares selecionados serão avaliados nas oito variantes de config_ancoras_textura.json: principal; faixa 30–49 com os centros originais; apenas remoção da média; deslocamentos opostos −16/+16 e +16/−16 px; redetecção com σ=12 e σ=20; núcleo literal v5 nos mesmos segmentos. As variantes de redetecção procuram o candidato do mesmo tipo mais próximo, até 64 px, com a mesma regra de proeminência/largura e suporte. Sem correspondência, o par fica não estimável naquela variante, sem retorno silencioso ao centro original. Não se refaz o pareamento por largura após redetectar. Contrastes entre variantes serão apresentados também no conjunto comum de pares de eventos, para separar perda de cobertura de mudança de valor.

Além das coordenadas σ12/20, será registrado o extremo local da faixa 30–49 como diagnóstico de representação; a sensibilidade wavelet da faixa conserva o centro principal. Não haverá escolha de uma faixa por evento. Deslocamentos de âncoras são medidos antes da WCO. Mudar σ pode mudar o marco selecionado; a tolerância de 64 px limita esse risco sem demonstrar identidade física.

## Agregação, critérios e preservação

Primeiro se resume cada grupo por par de aquisições; depois se calcula a média por padrão e a versão balanceada. Eventos sobrepostos e pares que compartilham aquisições não serão tratados como réplicas independentes. A retirada de uma aquisição remove todos os pares incidentes; a análise sem m03q40_001 será explicitada como influência de uma aquisição com artefatos, sem mudar a análise principal. Não serão produzidos valores-p ou intervalos artificiais de confiança.

O sucesso operacional pedido pelo usuário será medido por: cobertura de âncoras, rastreabilidade, estabilidade das coordenadas e possibilidade de calcular WPS/WCS/WCO com suporte. Isso é distinto de aumentar a coerência. A margem descritiva de 0,05 e cobertura de 2/3 permanecem referências para um possível ganho, exigindo contraste positivo frente a controles e estabilidade; não são limiares físicos. Ganhos pequenos serão chamados marginais; um resultado sem ganho será preservado.

Antes do lote, conferir recortes e novas operações de detecção/pareamento, validar qualquer cache de coeficientes contra o núcleo original e medir o custo atual em poucos pares. Informar a estimativa ao usuário. Configuração, código de anotação, catálogo, parceiros e revisão de imagens serão registrados por hash antes das WCO experimentais. As figuras deverão explicitar churn e wavy annular e os IDs de x/y. O resumo executivo e o passo a passo serão atualizados com resultados e limitações das duas etapas; a versão anterior do resumo e PDF foi preservada.
