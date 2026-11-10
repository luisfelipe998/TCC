# Revisão final dos gráficos e do resumo em PT-BR

Após a regeneração completa, foram inspecionadas visualmente as oito figuras finais alteradas: os dois diagramas de dispersão, recorte, panoramas com máscara, duas matrizes de confusão, média/mediana de largura clara e perfis radiais. As prévias não foram usadas para este parecer.

## Uma correção de redação necessária

Os textos `EXECUTIVE_SUMMARY_PT_BR.md` e `GRAFICOS_PARA_DISCUSSAO_PT_BR.md` dizem “25 aquisições sem saturação”. O critério efetivo é **S ≤ 0,20**, não S = 0. A tabela de procedência mostra 10 das 25 aquisições elegíveis com algum pixel saturado; o maior S elegível é 0,14213570778629922. Portanto, substituir essa expressão por **“25 aquisições elegíveis pelo critério de saturação”** ou **“25 aquisições sem saturação excessiva pelo critério S ≤ 0,20”**. Não é necessário mudar a seleção ou os números.

## Verificações que passaram

- Os dois diagramas exibem somente as 25 aquisições com `plotted=True`; `m03q40_001` está excluída dos pontos e envoltórias. O CSV conserva 26 registros para rastreabilidade.
- O gráfico principal mantém o caminho solicitado, usa PT-BR e separador decimal com vírgula; não contém anotação de teto, referência à aquisição de bolhas ou notas internas desnecessárias.
- O eixo Y descreve trocas entre pixels de colunas vizinhas, com média de 53 linhas. A figura secundária identifica a linha fixa `y = 40`.
- O resumo explica `[−]`, o caráter adimensional de α, o cálculo das áreas circulares, a média espacial e o denominador `H(W−1)` do flip-rate. As limitações aparecem nas respectivas etapas.
- As cinco medianas de α e de `1.000 F` no resumo coincidem, após arredondamento e conversão para vírgula decimal, com as tabelas regeneradas. Os valores de recorte e a alteração das máscaras acompanham os resultados corrigidos.
- A figura de máscaras está em PT-BR. As duas matrizes de confusão agora incluem `n = 25`, validação interna e ausência de exemplo de bolhas no treino do teste dessa aquisição.
- O gráfico exploratório de largura clara contém 25 categorias, pares de pontos e nenhuma trajetória entre aquisições. Os perfis radiais identificam as faixas vermelhas como externas à janela nominal, sem tratá-las como paredes calibradas.
- A revisão de gráficos para o TCC prioriza três figuras e apresenta somente duas sugestões novas, explicitamente ainda não geradas, com fontes de dados e fórmulas adequadas.

Não foi encontrada outra discrepância relevante nesta revisão visual e documental. Ela verifica coerência de apresentação e rastreabilidade; não estabelece calibração da fração volumétrica real.

## Resolução após a revisão

A expressão foi corrigida para aquisições elegíveis pelo critério de saturação nos dois resumos, na revisão de gráficos e nos geradores pertinentes. A seleção continua explicitamente definida por S ≤ 0,20. A pedido do usuário, o eixo Y do gráfico principal passou a se chamar “Flip-rate espacial (transições por 1.000 pixels)”; a escala e a média das linhas estão explicadas no texto. A regeneração final passou pelas 15 etapas e pelas nove categorias de auditoria após os resumos.
