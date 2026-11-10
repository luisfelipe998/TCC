# Protocolo da triagem de pistões líquidos condicionais

Registrado antes de executar esta etapa experimental, após a revisão de conservação e durante a execução do rastreador de objetos. Não se modificam os critérios para obter candidatos.

O resultado pretendido é o intervalo **condicional** `U_s/c ∈ [α_b,1]` para o líquido no pistão adjacente a uma bolha candidata. Não se calcula a média líquida da célula ou da gravação, nem um ponto estimado para a velocidade líquida. As hipóteses físicas estão em [physical_identifiability_review_PT_BR.md](physical_identifiability_review_PT_BR.md).

## Critérios fixados

- Padrão manual plug ou slug; trajetória de objeto grande no limiar 0,70. Todas as trajetórias de referência permanecem na tabela, inclusive fora do domínio desta hipótese.
- Trajetória aceita pelo rastreador e por seu critério `strict_deformation`; razão máximo/mínimo do volume reconstruído ≤1,10; diferença entre velocidades dos centroides volumétrico e projetado ≤`max(1 px/quadro, 10% do módulo da velocidade projetada)`.
- Como condição de domínio da razão, a velocidade projetada deve ter módulo maior que `10⁻¹² px/quadro` e ambos os centroides devem indicar o mesmo sentido. O eixo físico da derivação é orientado nesse sentido.
- Correspondência entre os mesmos objetos nos limiares 0,60, 0,70 e 0,80: classe grande, distância euclidiana de centroides ≤10 pixels, razão de áreas entre 0,5 e 2,0. A correspondência deve ser única em ambas as direções entre candidatos que satisfaçam esses critérios. Exige-se correspondência simultânea com os dois limiares alternativos em pelo menos 80% dos nós da trajetória.
- A partir da união das máscaras originais v15 em 0,60/0,70/0,80, procurar pelo menos dez colunas contíguas sem qualquer pixel marcado como gás nas linhas 18:71. A busca é axialmente externa à caixa do objeto de referência, nas fatias `[xmin−100:xmin−5]` e `[xmax+5:xmax+100]`, limitadas pela imagem. A mesma lateral deve fornecer esse trecho em pelo menos 80% dos nós. Essa regra é mais estrita que alternar de lateral em cada quadro.
- A ausência de gás na máscara é apenas evidência bidimensional. Não comprova pureza líquida da seção, qualidade da segmentação, ausência de gás fora do plano nem a hipótese `0≤U_f≤c`.
- Para os nós correspondidos simultaneamente, preservar mínimo, mediana e máximo de `central_alpha_conditional` de cada limiar. Usar o menor desses valores como extremo inferior do envelope de sensibilidade de `U_s/c`. Esse extremo combina sensibilidade de limiar e variação temporal; **não é um limite geométrico certificado nem intervalo de confiança**.
- Registrar as falhas de cada critério. Não promover ausência de candidato a velocidade zero. Se nenhuma trajetória passar, publicar esse resultado sem relaxar os critérios.

O script [conditional_bubble_bounds.py](../scripts/conditional_bubble_bounds.py) somente executa a leitura experimental após existir `object_validation.json` com as 26 gravações completas. Figuras nativas, quando houver candidatos, exibem caixa do objeto e o trecho sem máscara que sustenta a triagem óptica.
