# Contraprova independente do efeito de composição dos grupos

Esta revisão foi feita após a avaliação inicial da v16, por solicitação da frente de consolidação. Foram inspecionadas as linhas 130–157 de `v17/scripts/exploracao_espectral.py`, em particular a retirada dos pares incidentes aos registros excluídos, a ordenação canônica dos rótulos e as duas ponderações. Não se identificou erro nesse cálculo.

A reprodução independente não chamou `composition()` e não usou seu CSV de entrada. Recalculou Pearson a partir das potências arquivadas em `v16_spatial2_18_71/intermediario/descritores.npz`, com regressão escalar, e leu os escores de WCO de `intermediario/pares.json`. Percorreu índices de aquisições e reuniu os grupos por rótulo. O CSV produzido pela frente principal foi utilizado somente ao final, para comparação.

O maior erro absoluto entre os oito contrastes reproduzidos e os resultados da frente principal foi **5,55 × 10⁻¹⁶**.

| Conjunto | Pares intra/inter | Δ Pearson por par | Δ Pearson por grupos | Δ WCO por par | Δ WCO por grupos |
|---|---:|---:|---:|---:|---:|
| 26 aquisições | 87 / 238 | 0,065889 | 0,145300 | −0,075888 | −0,064322 |
| 25 sem m07q01, bolhas dispersas | 87 / 213 | 0,028735 | 0,016502 | −0,064628 | −0,023091 |
| 25 sem m03q40, anular | 79 / 221 | 0,054941 | 0,140256 | −0,073040 | −0,065992 |
| 24 sem ambos | 79 / 197 | 0,015936 | 0,007977 | −0,061168 | −0,024219 |

O registro de bolhas dispersas não participa de pares intrapadrão, porque seu grupo tem um único ensaio. Ele participa de quatro das dez combinações interpadrões: bolhas–churn, bolhas–anular, bolhas–plug e bolhas–slug. Portanto, tem participação de **40% da média inter por grupos**, apesar de aparecer em apenas **25/238 = 10,50% dos pares inter**. A ponderação por grupos muda o peso da pergunta, sem criar quatro novas aquisições de bolhas dispersas.

A exclusão desse registro conserva todas as médias intrapadrão e todas as médias dos seis grupos interpadrões restantes. Se B₆ é a média dos seis grupos sem bolhas e B₄ é a média dos quatro grupos com bolhas, então:

`Δ_macro,26 = média_intra − (0,6 B₆ + 0,4 B₄)`

`Δ_macro,25 = média_intra − B₆`

Logo, a redução do contraste é exatamente `0,4(B₆ − B₄)`. Nos dados, ela vale 0,128798, ou 88,64% do contraste macro original. A redução por par é 0,037154, ou 56,39%. Essas porcentagens descrevem a sensibilidade aritmética ao escopo; não são porcentagens de efeito causal “explicadas” por um fator.

O achado é material: o contraste agregado dos cinco padrões oferece informação muito mais favorável do que o contraste restrito aos quatro padrões que permitem comparação intrapadrão. No TCC, a descrição de uma base favorável à continuidade da análise deve ser entendida como uma observação do conjunto particular de 26 ensaios. Não demonstra separabilidade generalizável de todos os regimes, e especialmente não demonstra repetibilidade do padrão de bolhas dispersas.

Há três contrapesos necessários à interpretação:

1. **Não é erro da média balanceada.** O TCC define o peso igual aos quatro grupos intra e dez grupos inter em `metodologia.tex:351–366`. O cálculo responde exatamente a essa definição. Para descrever as cinco categorias observadas, a dissimilaridade do ensaio de bolhas é um resultado real da amostra.
2. **Os contrastes de Pearson continuam positivos.** Nenhuma das quatro exclusões avaliadas inverteu o sinal, e a revisão não autoriza dizer que toda a tendência desapareceu. O que muda fortemente é a magnitude e o suporte à interpretação geral.
3. **Excluir não equivale a corrigir os dados.** Retirar o único ensaio muda a pergunta para quatro padrões e deve ser apresentado como outro escopo. A exclusão não deve ser escolhida para obter um resultado conveniente. O problema de replicação de bolhas é resolvido por novas aquisições, não por uma transformação da média.

Recomendação: relatar conjuntamente os escopos de cinco padrões observados e de quatro padrões com replicação, as médias por combinação e as sensibilidades já calculadas. Caso se deseje inferência, o ensaio deve ser a unidade de permutação ou reamostragem, com as hipóteses de permutabilidade e os fatores operacionais explicitados. Janelas e pares não devem inflar o número de réplicas.

Reprodução: [contraprova_composicao.py](../scripts/contraprova_composicao.py). Resultados: [CSV independente](../resultados/wavelets/contraprova_composicao.csv) e [verificação numérica](../resultados/wavelets/contraprova_composicao.json). Nenhum arquivo anterior à v17 foi alterado.
