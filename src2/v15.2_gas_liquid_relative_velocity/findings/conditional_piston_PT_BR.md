# Triagem experimental de velocidades condicionais do pistão líquido

**Resultado: 0 trajetórias candidatas em 0 gravações.** Foram examinadas 55 trajetórias de objetos grandes nos seis vídeos plug/slug, mantendo na tabela todas as 6602 trajetórias de referência dos 26 vídeos. Os critérios foram registrados no [protocolo](conditional_piston_protocol_PT_BR.md) antes da execução desta etapa.

Os valores publicados são intervalos `U_s/c`, condicionados à reconstrução de área transversal, à interpretação de uma bolha isolada quase rígida com pistão líquido adjacente, à ausência relevante de mudança de fase e à hipótese **0 ≤ velocidade média do filme ≤ velocidade da bolha**. Não são duas velocidades materiais medidas independentemente. `U_s` corresponde somente ao líquido no pistão; a média líquida de uma célula completa ou da gravação não foi calculada.

## Resultados por gravação plug/slug

| Gravação | Padrão | Trajetórias grandes | Aceitas originalmente | Rígidas/volume/centroides | Candidatas condicionais | Intervalos distintos observados |
|---|---|---:|---:|---:|---:|---:|
| m04q01_001 | slug | 0 | 0 | 0 | 0 | 0 |
| m05q01_001 | plug | 1 | 0 | 0 | 0 | 0 |
| m05q05_001 | slug | 27 | 0 | 0 | 0 | 0 |
| m06q01_001 | plug | 0 | 0 | 0 | 0 | 0 |
| m06q05_001 | slug | 6 | 0 | 0 | 0 | 0 |
| m07q05_001 | plug | 21 | 0 | 0 | 0 | 0 |

## Dificuldades observadas

| Critério isolado | Trajetórias grandes que passaram |
|---|---:|
| Aceitação original da trajetória | 0/55 |
| Deformação estrita | 0/55 |
| Variação do volume ≤10% | 31/55 |
| Concordância entre centroides | 43/55 |
| Correspondência simultânea entre limiares ≥80% | 0/55 |
| Contexto sem máscara na mesma lateral ≥80% | 34/55 |

As contagens são marginais; falhas podem se sobrepor. No rastreador original, as rejeições dos objetos grandes foram: `centroid_acceleration`: 47, `front_rear_deformation`: 19, `lag2_correspondence`: 28, `nonlinear_centroid`: 10, `short_track`: 44.

A maior fração de nós com a mesma identificação nos três limiares foi 50.0%, abaixo dos 80% exigidos. O simples achado de um trecho adjacente sem máscara não superou a instabilidade das trajetórias e da identidade entre segmentações. Não se enfraqueceram os filtros para converter esse contexto em velocidade líquida.

## Como interpretar o intervalo

A conservação no referencial da bolha fornece `U_s/c = α_b + (1−α_b) U_f/c`. O valor de `α_b` foi aproximado pelo indicador de área transversal circular do corpo da bolha. Em cada nó correspondido, compararam-se os mesmos objetos nos três limiares. O extremo inferior apresentado é o menor indicador entre esses limiares e nós; o superior é 1, imposto pela hipótese sobre o filme. O extremo inferior combina variação temporal e sensibilidade de segmentação: **não é um limite garantido da geometria real nem um intervalo de confiança**.

A triagem de contexto procurou uma mesma lateral da trajetória com ao menos dez colunas contíguas sem gás em nenhuma das três máscaras, em pelo menos 80% dos seus nós. Isso exige ausência de marcação em toda a altura do lúmen projetado. Ainda assim, máscara vazia pode significar falha de detecção e não exclui gás fora do plano de observação. A triagem não comprova a velocidade do filme.

As trajetórias foram condicionadas ainda à estabilidade de comprimento e de volume reconstruído, concordância dos dois centroides e correspondência única entre segmentações. Os controles de correspondência e contexto passaram, mas verificam apenas o algoritmo. A escolha dos critérios seleciona objetos fáceis de seguir e não garante representatividade do conjunto gasoso. As trajetórias da mesma gravação não são réplicas experimentais independentes.

Os valores em pixels/s usam o módulo da velocidade do centroide volumétrico condicional e o FPS nativo. A razão dentro de uma trajetória cancela uma escala espacial comum. Comparações de velocidade física entre gravações ainda exigem escala axial relativa ou hipótese explícita de magnificação igual.

Nenhuma trajetória satisfez simultaneamente os critérios. Não se produziram intervalos ou gráficos com valores preenchidos por hipótese. A rejeição desta triagem não demonstra inexistência de pistões físicos; indica que a evidência e os critérios aplicados não sustentaram esta redução nas trajetórias disponíveis.

## Arquivos e reprodução

- [Todas as trajetórias, critérios e falhas](../results/conditional_piston_tracks.csv).
- [Contexto e correspondências por nó](../results/conditional_piston_nodes.csv).
- [Sumário dos 26 vídeos](../results/conditional_piston_summary.csv).
- [Controles, critérios e proveniência](../results/conditional_piston_validation.json).
- [Derivação física e limites de interpretação](physical_identifiability_review_PT_BR.md).

Executar `python3 -B src2/v15.2_gas_liquid_relative_velocity/scripts/conditional_bubble_bounds.py` após concluir `object_motion.py`. A etapa lê os catálogos de objetos, as máscaras originais v15 e os quadros nativos usados nos exemplos. Não altera imagens nem resultados anteriores.
