# Revisão independente — correlação por blocos e parede selecionada

**A aritmética, os critérios salvos, a seleção e as médias foram reproduzidos sem erro detectado.** A versão final corrigida, com fundo estimado somente na metade de seleção, contém **25/26 pares ópticos** com uma parede selecionada, em cada direção. **22/25** concordam ao inverter as metades e **7/26** passam e concordam em todas as alternativas declaradas. O resultado de referência anterior, usando ambas as paredes, contém20/26 pares. Nenhum desses números estabelece26 pares aceitos ou velocidades materiais de gás/líquido.

## Correção de compartilhamento entre as metades

A auditoria inicial identificou que `selected_wall_per_run.csv` utilizava perfis com fundo estimado pela média das501 imagens. Esse experimento foi preservado como **arquivado, com fundo compartilhado**. A análise final está em `split_background_details.csv` e `split_background_pairs.csv`: na direção inicial→final, o fundo usa somente os quadros0–250 e é congelado para avaliar251–500; na direção oposta, usa somente251–500 e avalia0–250. A seleção da parede também usa exclusivamente sua metade de seleção.

Foi conferido que não existe interseção entre os quadros usados para o fundo e os de avaliação. Os perfis, picos e controles foram recalculados por implementação independente com o fundo correto. **O compartilhamento do fundo foi corrigido na versão final.** Isso ainda não constitui validação experimental externa: o desenho do método explorou essas mesmas gravações. A filtragem8/32 foi repetida com o fundo restrito à seleção, sem transportar automaticamente os resultados anteriores.

## Escopo da conferência

Foram reconstruídas as decisões dos **5070 registros regionais/blocos**, as **156 linhas por gravação/alternativa** e as **208 linhas de seleção de parede** diretamente do CSV detalhado, sem importar os geradores. As regras conferidas incluem textura, picos, controle embaralhado, atraso2, reversão, posições160/320, tolerância, índices de quadros, blocos simultâneos, cobertura temporal, pesos por intervalos, escolha no conjunto de seleção e médias apenas nos blocos aceitos de avaliação.

Adicionalmente, foram reconstruídas **todas as 4680 decisões regionais/blocos e 156 seleções da versão final corrigida**, incluindo as duas direções e três alternativas de processamento. A correção não alterou as contagens25/25 da referência, mas suas coordenadas e verificações foram recalculadas e conferidas; não se presumiu equivalência.

FPS foi comparado com dois campos nativos em cada uma das26 gravações; foram conferidos501 nomes consecutivos e perfis de três imagens nativas por gravação. A média das duas paredes utiliza blocos simultâneos e pesos idênticos: a média das duas médias é matematicamente correta para esse descritor regional.

Recalcularam-se **132 blocos representativos** (60 arquivados e 72 da versão final), incluindo m03q05, m03q40, m04q01, m05q01 e m07q01, com outra implementação de correlação linear e normalização de Pearson. Foram verificadas **15801 avaliações de Pearson por produtos centrados explícitos**, **1848 curvas médias** e **1188 controles embaralhados regenerados**. A maior diferença de Pearson foi 7.87e-13; a maior diferença de velocidade foi 8.88e-13 px/quadro. Esses números verificam implementação, não erro experimental.

## Resultados reproduzidos

| Alternativa de uma parede | Seleção inicial→final | Seleção final→inicial |
|---|---:|---:|
| **Final: fundo somente da seleção**,10 blocos | 25 | 25 |
| Sem subtração temporal | 24 | 23 |
| Fundo somente da seleção + passa-altas32 | 22 | 22 |
| Fundo somente da seleção + passa-altas8 | 8 | 8 |
| Arquivado: fundo compartilhado501 quadros | 25 | 25 |

As alternativas de cinco e vinte blocos forneceram 19 e 18 pares com as duas paredes. A cobertura maior da parede selecionada acompanha uma mudança explícita da grandeza: uma região escolhida por observabilidade substitui a média de ambas.

## Pontos que precisam aparecer no resumo

1. **m03q05 continua rejeitado.** Na seleção inicial→final, a parede inferior foi escolhida, mas somente2/5 blocos conjuntos na seleção e1/5 na avaliação passaram. O ponto finito usa apenas esse bloco de49 intervalos. Ele pode aparecer como candidato rejeitado; não completa26 pontos aceitos.
2. **A versão final separa o fundo e a escolha da parede entre seleção/avaliação, mas não é validação externa.** O vazamento de fundo existia no experimento arquivado e foi corrigido. A escolha do desenho analítico após explorar o conjunto permanece como limitação; não se deve afirmar que os25 resultados foram confirmados em um conjunto experimental independente.
3. **As coordenadas reportadas pertencem à metade de avaliação e a seus blocos aceitos.** Na direção inicial→final, elas usam29,4%–49,0% dos500 intervalos do vídeo. Isso não é a média integral do ensaio. Núcleo e parede usam os mesmos blocos, o que evita uma mistura temporal entre as duas coordenadas.
4. **Inversão das metades e processamento:** na versão final, **22/25** pares concordaram em ambas as coordenadas ao inverter as metades, dentro de max(2px/quadro,20% do módulo da referência). Os casos discordantes foram m04q01_001, m04q10_001, m05q05_001. Exigindo aprovação e concordância também com a variante sem fundo e com as duas filtragens, restaram **7** gravações: m05q01_001, m05q20_001, m06q01_001, m06q05_001, m06q10_001, m07q01_001, m07q05_001. Isso combina variabilidade temporal, mudança de região, escala óptica e mudança dos blocos selecionados; não é um erro experimental calibrado.
5. **m03q40 depende do processamento.** Ele passa na referência de parede selecionada e falha na variante sem subtração e em ambas as variantes passa-altas. A saturação original e essa dependência precisam acompanhar seu ponto. Os controles sintéticos de recorte mostram possibilidade de recuperar movimentos sob saturação, mas não validam automaticamente esta gravação.
6. **As falhas passa-altas não provam velocidades fisicamente falsas.** Alterar a escala óptica pode mudar a estrutura dominante, sua textura e a seleção dos blocos. A evidência limita a representatividade de um único valor e precisa ser mantida como sensibilidade.
7. **Não houve identificação material.** Uma correlação coerente pode acompanhar ondas, interfaces ou brilho compartilhado entre as regiões. A localização núcleo/parede e a escolha da região mais observável não fornecem, por si, gás/líquido material.
8. **A outra parede foi salva como candidata no experimento arquivado.** `other_wall_same_epochs_candidate_px_s` não exige suporte dessa outra parede nos mesmos blocos. Não deve ser apresentada como uma segunda medida de parede validada nem confundida com a parede da versão final.

## Disposição

Os cálculos finais, com fundo restrito à metade de seleção, podem ser utilizados como **descritores de propagação óptica sob o protocolo declarado**, com25 pares na referência exploratória e um candidato explicitamente rejeitado. A sensibilidade ampliada retém apenas7 gravações aprovadas e concordantes em todas as alternativas; os25 pontos não devem ser apresentados como invariantes ao método. O gráfico completo precisa preservar o status do26º ponto, o caráter regional dos eixos, a escolha de uma parede e o período efetivamente amostrado. Nenhuma das conferências converte esses resultados em26 pares de velocidades materiais.

Reprodução da auditoria: `python3 -B src2/v15.1_gas_liquid_velocity/all_scenarios_adversarial/scripts/independent_block_verify.py`. Evidência estruturada: [independent_block_verification.json](../results/independent_block_verification.json).
