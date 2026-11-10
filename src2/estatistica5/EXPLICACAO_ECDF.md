# Como ler o gráfico ECDF enviado

O gráfico enviado é o da estatistica4: reúne 171 pares entre 19 ensaios (10 churn e 9 wavy annular), com descritores calculados nos primeiros 7.001 px. Ele não corresponde ao conjunto de 87 pares utilizado na matriz da estatistica5. Esta explicação se refere à figura anterior; não foi gerada uma nova análise ECDF.

## Eixos

- **Horizontal:** dissimilaridade d_antiga = 1−r de Pearson. Quanto menor, mais parecida a forma dos dois descritores. Nessa figura, a distância não foi dividida por dois.
- **Vertical:** fração dos pares da curva com distância menor ou igual ao valor horizontal. A definição é F(t) = número de pares com d ≤ t / número de pares naquela curva. Cada curva usa seu próprio denominador.

Por exemplo, x = 0,10 corresponde a r ≥ 0,90. Ler y = 0,47 nesse ponto significa que cerca de 47% dos pares daquela curva têm Pearson de pelo menos 0,90. Não significa 47% de chance de classificar corretamente o padrão.

Com os valores salvos que deram origem à figura:

| Curva | Pares | Pares com d ≤ 0,10 | Fração |
| --- | --- | --- | --- |
| Verde: mesmo padrão | 81 | 38 | 46,9% |
| Laranja: padrões diferentes, códigos q próximos | 36 | 14 | 38,9% |
| Azul: padrões diferentes, códigos q distantes | 54 | 23 | 42,6% |

“Próximos” conserva o critério histórico |Δq| ≤ 20; “distantes”, |Δq| > 20. Trata-se do código q usado na base, sem inferir unidades não documentadas.

Uma curva mais à esquerda concentra pares em distâncias menores. Equivalentemente, em um mesmo limite horizontal, uma curva mais alta contém uma fração maior de pares com distância até esse limite. As curvas sobem em degraus porque cada par adiciona uma observação à distribuição acumulada.

## Por que o gráfico é relevante

Ele mostra se a semelhança encontrada dentro de um padrão é diferente daquela encontrada entre padrões. Se o descritor distinguisse bem as classes neste conjunto, seria desejável que a curva verde se concentrasse claramente em distâncias menores que as curvas entre padrões.

Na figura enviada há forte sobreposição. As medianas de 1−r foram aproximadamente 0,1124 dentro do padrão, 0,1253 entre padrões com q próximo e 0,1174 entre padrões com q distante. Portanto, **perfis de padrões diferentes também podem ser bastante parecidos**. Pearson alto em alguns pares não basta para concluir que o descritor separa os padrões.

Os pares compartilham ensaios: a curva verde não representa 81 aquisições independentes. O gráfico é descritivo; sua sobreposição não prova equivalência entre padrões, e diferenças visuais não estabelecem significância ou causalidade.

## Relação com a nova matriz

Na estatistica5 a escala é (1−r)/2. Uma distância antiga de 0,10 corresponderia a 0,05 nessa nova escala **se os dois descritores fossem exatamente os mesmos**. Dividir por dois apenas muda a escala, preservando a ordem das distâncias.

Entretanto, a nova matriz preserva os trechos específicos de cada par da v16, cujo comprimento é min(N_x,N_y), e usa somente os 87 pares intrapadrão principais. Seus valores não são obtidos apenas dividindo os valores da figura antiga por dois. Sem células entre padrões, a matriz nova descreve semelhança intrapadrão e não responde à pergunta de separação entre padrões que motivava a ECDF anterior.

Fonte numérica: [pair_distances.csv da estatistica4](../estatistica4/results/pair_distances.csv) e [distance_summary.csv](../estatistica4/results/distance_summary.csv).
