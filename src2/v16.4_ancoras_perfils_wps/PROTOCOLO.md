# Âncoras como entrada para os perfis WPS

## Pergunta e decisão de realizar o estudo

A proposta foi interpretada como **usar a média por coluna do recorte centrado na âncora como sinal de entrada da v16.1**. Uma curva por aquisição foi preservada. Não se aplicou uma segunda wavelet sobre uma WPS já calculada, nem se construiu uma curva média por classe.

O encadeamento tem fundamento: a v16.3 define o conteúdo espacial observado, enquanto a v16.1 compara a distribuição da WPS entre escalas. As duas usam as mesmas linhas, preparo e núcleo Morlet. A relevância é testar se a seleção óptica de um contexto melhora sua representação espectral; relevância para investigação não implica ganho demonstrado.

**Este protocolo foi escrito depois da triagem dos perfis já disponíveis.** Essa triagem encontrou A por par de 55,73% nas âncoras e 76,26% nos panoramas completos, nos mesmos 25 ensaios e 50–1.024 px. A ausência de ganho inicial motivou controles de comprimento e posição. Não se apresenta o estudo como pré-registro independente dos resultados. Os controles abaixo foram definidos antes de calcular seus resultados nesta pasta.

## Entradas e representação

- Reutilizar as 25 âncoras congeladas da v16.3: 3 plug, 3 slug, 10 churn e 9 wavy annular. A aquisição dispersed bubbles da v16.1 não tem âncora neste catálogo e permanece fora de todas as comparações principais.
- Conservar os centros individuais escolhidos no contexto nativo, usando o perfil `comparavel` de 2.048 px. Não misturar janelas nativas de 2.048 e 4.096 px.
- Extrair a média das linhas 18–71 inclusive; recortar antes do preparo. Usar 255−I, remoção da tendência linear, média zero e desvio padrão populacional unitário.
- Autocorrelação bilateral com divisor N; Morlet ω₀=6; K=1; λ=4πa/(6+√38). A WPS alternativa é |CWT(Rxx)|, sem quadrado adicional.
- Usar **WPS no lag zero**, como a representação principal da v16.1. Não confundir com a curva de potência média sobre lags da v16.3.
- Grade principal fixa: 200 pontos logarítmicos entre 50 e 1.024 px, herdada da v16.3. Todos atendem N−1≥√2a no lag zero para N=2.048. Esse suporte é geométrico, não um teste estatístico.
- Normalizar separadamente cada perfil pela integral trapezoidal em log λ: pⱼ=wⱼPⱼ/ΣwₖPₖ. A densidade desenhada é P/ΣwP; sua área é um. As massas não são energia física.
- H(p,q)=√[Σ(√p−√q)²/2]; menor significa maior semelhança. Importar diretamente as operações da v16.1, conferindo a identidade dos núcleos.

## Controles definidos

1. **Panorama completo:** mesmos 25 ensaios, mesma extração, preparo, métrica e faixa de 50–1.024 px. Comprimento e conteúdo mudam juntos.
2. **Início, centro e fim:** janelas de 2.048 px; início na coluna zero; centro em floor((N−2.048)/2); fim em N−2.048. Mesma grade e preparo após o recorte. Esses controles mantêm o comprimento e variam a localização/conteúdo, permitindo examinar a contribuição da escolha da âncora. Janelas sobrepostas não são réplicas.
3. **Deslocamento de −16/+16 px:** perturbar o centro, mantendo o comprimento. A tolerância da v16.3 motivou a amplitude. Não refazer a seleção; não interpretar esses contextos como novas âncoras morfologicamente aprovadas. Não extrapolar ou preencher sinais; falhar se o recorte sair do panorama.
4. **Faixa conservadora:** normalizar novamente apenas os pontos com 3a≤2.048/2, comuns a todas as condições. Não escolher a faixa que produzir melhor separação.
5. **Cosseno:** única sensibilidade de métrica, com os pesos de integração da v16.1, nas âncoras e nos panoramas completos.
6. **Referências v16.3:** Pearson assinado dos resíduos log-lineares da WPS e WCO média central, usando as tabelas existentes dos mesmos 300 pares. Não recalcular nem alterar os mapas cruzados. A WCO central da v16.3 não é o escore de máximos da WCO antiga comparada na v16.1.
7. **Referência v16.1:** restringir sua representação principal aos mesmos 25 ensaios. Sua faixa de 50–2.295,94 px permanece identificada e não isola o efeito da âncora. Os 76,5%/81,6% publicados com 26 ensaios não serão usados como controle direto.

## Avaliação e interpretação

Calcular todos os 300 pares: 87 intra padrão e 213 inter padrão. A é a fração das 87×213 comparações em que a similaridade intra supera a inter, com meio crédito por empate. A=50% é ausência de preferência de ordenação, não um limiar de significância. Nenhum classificador é treinado.

Relatar A com peso igual por par e com pesos iguais entre os quatro grupos intra e as seis combinações inter. Por padrão, usar seus pares intra e os pares inter que o envolvem, balanceando as três outras classes. Relatar médias/medianas de H, todos os grupos, e o efeito de retirar cada aquisição e seus pares; recalcular os pesos. Comparar diferenças de A sob a mesma retirada. Essas faixas são influência, não intervalos de confiança nem validação deixando uma aquisição de fora.

Não definir um efeito mínimo arbitrário. Só descrever vantagem geral das âncoras se a direção for consistente nas duas ponderações frente aos controles de mesmo comprimento e ao panorama completo, sem esconder perdas por padrão. Caso contrário, relatar ganho restrito ou ausência de melhora geral. Apresentar todos os controles, sem substituir a análise principal pelo melhor resultado.

## Limites que permanecem

As regras de seleção já dependem do padrão conhecido: eventos morfologicamente aprovados para plug/slug e vales para churn/wavy annular. Isso é adequado para um estudo descritivo condicionado aos rótulos, mas não demonstra um extrator independente da classe para classificar novos ensaios. A autocorrelação resume todo o recorte: centralizar uma estrutura não faz a WPS descrever exclusivamente essa estrutura. O benefício possível vem da seleção do conteúdo, não de alinhar o eixo de lag ao centro óptico.

Uma âncora por aquisição pode perder variedade intraensaio. Panoramas completos retêm mais eventos, mas têm comprimentos e condições operacionais desiguais. Costuras, pixels sem calibração física, rótulos sem validação independente, poucas aquisições plug/slug e dependência entre pares permanecem limitações. A análise não separa causalmente esses fatores.

## Preservação e reprodução

Todos os resultados novos ficam nesta pasta. Catálogo, caches, imagens, código e resultados anteriores são somente lidos; nenhuma função de preparação/regeneração é chamada. Os hashes das entradas consumidas são congelados na primeira execução e conferidos nas seguintes. A verificação confronta os perfis com correlação direta e Morlet temporal e confere Hellinger e A por formulações independentes. O relatório registra precisamente quais verificações foram feitas.
