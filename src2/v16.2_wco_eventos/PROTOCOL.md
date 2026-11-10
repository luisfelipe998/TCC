# Protocolo de comparação de eventos por wavelets

Este protocolo foi definido antes de calcular qualquer WCO nova. A análise é exploratória: o conjunto já foi estudado nas versões anteriores e a inspeção visual desta execução participa do desenvolvimento das regras. Não existe conjunto externo de validação. A base indicada pelo usuário é [v16_spatial2_18_71](../v16_spatial2_18_71/README.md); a entrada principal é a média das linhas **18–71 inclusive**, com indexação iniciada em zero.

## Observabilidade e anotação

O piloto inicial foi formado por m05q01_001, m04q01_001, m03q05_001, m03q10_001 e m07q01_001. A inspeção preliminar incluiu depois as outras quatro aquisições plug/slug, ainda sem wavelets, para avaliar extensão e contexto. Por isso nenhum desses seis registros é chamado de avaliação independente. As imagens de início, centro e fim e os quadros 0, 250 e 500 ficam no atlas. Essa amostragem não equivale a inspecionar manualmente todos os 13.026 quadros.

A categoria candidata principal é **estrutura alongada delimitada**: envelope escuro com interior visível, extensão horizontal e duas extremidades discerníveis na imagem. Trata-se de uma categoria óptica operacional, não de identidade física ou classificação automática do escoamento. Junções internas entre estruturas e limites indistintos são marcados como ambíguos, sem inventar uma fronteira. Texturas contínuas/fragmentadas não recebem limites de eventos apenas por oscilações da intensidade. Bolhas pequenas são documentadas separadamente; não são alongamentos artificiais de um evento grande.

A proposta automática de coordenadas utiliza somente a imagem: fração de pixels com intensidade <60 nas linhas 18–71; mediana de 5 colunas; fração ≥0,25; fechamento binário de 13 colunas. O preenchimento da borda do fechamento repete a borda observada, sem criar falsos limites a 6 px. São catalogados intervalos a partir de 96 px, inclusive rejeitados. A revisão visual aceita/rejeita a morfologia; não altera coordenadas para melhorar o sinal. Coordenadas são intervalos **[início, fim)**, com fim exclusivo. A âncora é floor((início+fim)/2). Estes limites por limiar não são uma segmentação física calibrada das interfaces.

Os intervalos precisam ter no máximo 1.536 px e caber inteiros no segmento de 2.048 px centrado na âncora, com pelo menos 256 px de contexto de cada lado. Exige-se margem adicional de 16 px na aquisição para a sensibilidade. Segmentos podem se sobrepor; a sobreposição e estruturas vizinhas são quantificadas. Não se tratam eventos sobrepostos como réplicas independentes. Intervalos tocando as extremidades do panorama, sinais degenerados e morfologia ambígua são rejeitados com motivo. A ausência de eventos sustentáveis é um resultado de observabilidade, não WCO igual a zero.

Costuras e deslocamentos não positivos são registrados por intervalo e por contexto. Não se elimina automaticamente m03q40_001 nem m04q01_001. O quadro que forneceu a faixa da âncora é recuperado das larguras acumuladas do v11; isso é proveniência computacional, não rastreamento da mesma bolha. A mistura de costura é assinalada. A mesma regra de limiar será repetida com intensidades 55 e 65 como diagnóstico de estabilidade das fronteiras, sem escolher a versão de melhor WCO. Não haverá alegação de concordância entre anotadores: a revisão é de um único agente.

## Pareamento e controles

Pares principais usam aquisições diferentes. Dois eventos aceitos da mesma categoria são **equivalentes pelo critério visual operacional** quando max(extensão)/min(extensão) ≤1,25. São **distintos em extensão** quando a razão é ≥1,50. O intervalo intermediário é indeterminado e é registrado, sem virar controle negativo. O contraste isola esta regra de extensão, não todos os atributos possíveis de semelhança visual. Não se usa WPS, WCS, WCO ou frequência para anotar ou parear.

Por par de aquisições e categoria de comparação usam-se todos os pares admissíveis até 12; acima disso, selecionam-se os 12 menores hashes SHA-256 de semente, IDs e categoria. A ordenação não depende de resultados. Pares entre os mesmos eventos são únicos. São separados: equivalentes intrapadrão; distintos intrapadrão; equivalentes entre padrões; distintos entre padrões; intermediários não avaliados. Comparações dentro da mesma aquisição são diagnósticos separados.

Controle de recorte: segmentos de 2.048 px com inícios 0, floor((N−2048)/2) e N−2048, cruzando apenas a mesma posição relativa entre aquisições. São **3 × 325 pares**, incluindo todos os padrões. Eles podem conter eventos; são referências sem anotação, não exemplos garantidamente sem bolhas. Um controle sintético de fundo comum sem evento equivalente testa especificamente essa vulnerabilidade. Não se atribui a esses controles a mesma qualidade morfológica dos eventos: essa diferença é uma limitação.

## Cálculo e domínio

Recorta-se antes de inverter 255−I, remover tendência linear/média e normalizar pelo desvio padrão. A amostragem permanece 1 px/amostra; não se aplica redimensionamento, DTW, Hann ou deslocamento otimizado pela WCO. Usa-se uma realização por par, **K=1**, e divisor N explicitamente enviesado, adaptação das Eqs. 32.21–23. O núcleo é uma cópia idêntica e verificada por hash da base indicada.

Preservam-se Morlet ω0=6, λ=4πa/(6+√38), alternative WPS=|CWT(Rxx)|, alternative WCS=CWT(Rxy) complexa e alternative WCO=|WCS(b,a)|²/(WPSx(0,a)WPSy(0,a)). O índice de lag zero é N−1. O atraso positivo em y corresponde a lag negativo. Não se faz clipping. O piso relativo 10⁻¹⁰ protege apenas denominadores numéricos. Não se agregam coeficientes de eventos como ensemble.

O domínio comum é de **200 comprimentos de onda logarítmicos de 50 a 1.024 px**. O limite inferior conserva a faixa anterior; o superior é N/2, para que o segmento tenha ao menos dois comprimentos de onda. Não se reduz a banda após observar resultados. A média principal usa todos os 65 lags inteiros de **−32 a +32 px**, tolerância operacional para a âncora e interfaces difusas; não é incerteza experimental validada. Pesos uniformes nos lags e nas 200 escalas logarítmicas: média em lag por escala, depois média entre escalas válidas. Sem pesos de potência aprendidos por par. O escore anterior de máximos ponderados **não é reutilizado**.

As condições N−1−|b|≥√2 a (suporte geométrico) e |b|+3a≤N/2 (controle conservador opcional) são diagnósticos de suporte, não novas equações. Todo o domínio principal satisfaz o suporte geométrico. Os mapas integrais exibem bordas; o resumo principal não usa máximos fora do domínio definido. Exportam-se também valor médio em lag zero, máximo integral, média com suporte conservador e contagens de denominadores válidos. A média conservadora cobre outra banda e é identificada separadamente, sem substituir o resultado principal.

## Sensibilidades e separação das causas

O núcleo v5 literal será aplicado **aos mesmos segmentos já recortados**, preservando seu estimador, cmor, conversão de escala e lag refletido após a CWT. Como os comprimentos são iguais, o preparo é feito dentro de cada segmento; não se reintroduz tendência estimada do panorama. Essa comparação separa o núcleo numérico da anotação.

Os mesmos pares congelados serão recalculados para: faixa 30–49; apenas remoção da média; deslocamento de ambos os centros em −16 e +16 px; deslocamento oposto dos centros em −16/+16 e +16/−16 px. Nenhuma variante redefine os parceiros ou seleciona a maior coerência. Deslocar os centros altera as bordas e o contexto, mantendo o evento completo e N. Isso não é uma avaliação exaustiva de todas as possíveis anotações.

Resultados anteriores de panoramas são apresentados em tabela separada, com seu domínio de 50–3.000 px e resumo de máximos ponderados. A diferença numérica direta entre esse escore e a média de eventos não será chamada de ganho. Os controles de mesmo comprimento e domínio fornecem a comparação apropriada de recorte/pareamento.

## Agregação e relevância

Primeiro calcula-se a média por par de aquisições e grupo de comparação; depois média por padrão e versão balanceada com peso igual por padrão. O contraste equivalentes−distintos é **pareado pelo mesmo par de aquisições**, somente onde ambos existem; o mesmo vale para equivalentes−referências fixas. A cobertura de cada contraste é exportada. A sensibilidade por retirada de uma aquisição remove todos os seus eventos e pares; não é intervalo de confiança nem validação externa. Não se usam testes que considerem eventos ou pares como observações independentes.

Como critério descritivo definido antes dos cálculos, uma diferença de pelo menos **0,05** na escala de 0–1 é considerada candidata a relevância; não é um limiar físico validado nem significância estatística. Exige-se cobertura de pelo menos 2/3 das aquisições do padrão e ganho positivo frente a distintos e recortes fixos, conservado nas sensibilidades. Ganhos menores são chamados de melhora marginal. Se pares visualmente semelhantes entre padrões também tiverem WCO alta, isso aponta semelhança local e não separação dos rótulos globais. Sem cobertura, controles ou estabilidade suficientes, a conclusão será restrita/inconclusiva.

## Validações e custo

Antes do lote: correlação por soma direta, convolução temporal independente, índice zero, Eq. 32.17 recalculada, self, atraso conhecido, sinais com ruído, eventos distintos, fundo comum, bordas, denominadores degenerados, padding e reprodução literal do v5. Tempo de poucos pares de comprimentos distintos será medido e informado antes do lote. Anotação, cálculo e renderização terão tempos separados. O manifesto final registrará fontes, código, parâmetros e artefatos. A pasta será renomeada para `wco_eventos` após a execução e novamente verificada.
