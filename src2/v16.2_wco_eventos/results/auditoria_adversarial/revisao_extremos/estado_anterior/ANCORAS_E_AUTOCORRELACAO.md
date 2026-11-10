# Âncoras, recortes e autocorrelação na análise de eventos

Esta nota reúne a explicação do método executado em **v16.2_wco_eventos**, com base em **v16_spatial2_18_71**. As regras originais permanecem no [protocolo congelado](PROTOCOL.md). A extensão para churn e anular ondulado **foi executada** com um segundo catálogo, descrito abaixo e no [protocolo da extensão](PROTOCOL_CHURN_ANULAR.md). O catálogo original de plug/slug permaneceu intacto.

## Onde a âncora entra

A âncora define a origem do recorte antes do preparo dos sinais e antes das correlações. A CWT é aplicada às correlações desses recortes. O método não procura uma posição que maximize a WCO.

```mermaid
flowchart TD
    A["Imagens: localizar e revisar marcos ópticos sem consultar WCO"] --> B["Definir âncoras e parceiros pelo critério óptico"]
    B --> C["Recortar x e y: 2.048 amostras em cada aquisição"]
    C --> D["Inverter intensidade, remover tendência/média e padronizar"]
    D --> E["Calcular Rxx, Ryy e Rxy dos recortes"]
    E --> F["Aplicar CWT às três correlações"]
    F --> G["Obter WPS, WCS e WCO alternativas"]
```

Os sinais x e y vêm de eventos de aquisições diferentes nas comparações principais. Eles não são as metades esquerda e direita de um evento. Cada autocorrelação usa duas cópias do mesmo sinal; a correlação cruzada compara x com y.

## Evento, âncora e limites dos sinais

A definição de evento por intervalo [L,R) abaixo pertence à **etapa 1: plug/slug**. Na extensão **churn/anular**, c é um mínimo local do perfil óptico suavizado; não se exige um evento físico delimitado. Os limites [c−1024,c+1024) e as equações de correlação/wavelet são os mesmos nas duas etapas.

O sinal bruto de cada aquisição é a média das **54 linhas de 18 a 71 inclusive**, com indexação iniciada em zero:

$$I[p]=\frac{1}{54}\sum_{r=18}^{71}P[r,p].$$

P é a imagem em níveis de cinza e p é a coluna do panorama. A amostragem é **1 pixel por amostra**, sem conversão para milímetros ou tempo.

O detector e a revisão visual delimitam um evento no intervalo **[L, R)**. O início L é incluído; o fim R é exclusivo. Sua extensão é R−L e a âncora é:

$$c=\left\lfloor\frac{L+R}{2}\right\rfloor.$$

Para cada evento, o recorte tem **N=2.048** amostras:

$$[s,e)=[c-1024,\;c+1024),\qquad e-s=2048.$$

Para as aquisições A e B, os sinais ainda sem preparo são:

$$x_0[n]=I_A[c_x-1024+n],\qquad y_0[n]=I_B[c_y-1024+n],\quad n=0,\ldots,2047.$$

As posições globais dos dois recortes podem ser diferentes. Em cada recorte, a âncora ocupa o índice local **1.024**, e as posições relativas a ela vão de **−1.024 a +1.023 px**. O ponto c+1024 não pertence ao recorte.

O catálogo aceita eventos com extensão de 96 a 1.536 px, inteiros no segmento e com pelo menos 256 px de contexto de cada lado. Também exige 16 px adicionais disponíveis nas bordas da aquisição para as sensibilidades de deslocamento dos centros. Um evento sem contexto suficiente é rejeitado; não se preenche a imagem com observações artificiais para fazê-lo caber.

**Exemplo real, usado no caso favorável do resumo:** os dois eventos abaixo pertencem a **pistonado/bolhas alongadas (plug)**. Todas as coordenadas são colunas do respectivo panorama, com fim exclusivo.

| Sinal e evento | Limites do evento [L, R) | Âncora c | Limites do recorte [s, e) | Amostras |
| --- | --- | --- | --- | --- |
| x — m05q01_001_E005 | [2.752, 3.093) | 2.922 | [1.898, 3.946) | 2.048 |
| y — m06q01_001_E011 | [5.423, 5.836) | 5.629 | [4.605, 6.653) | 2.048 |

Esses números vêm do [catálogo de eventos](results/catalogo_eventos.csv). O exemplo favorável foi selecionado após o lote somente para ilustração; não definiu as regras nem os parceiros.

## Preparo dentro de cada recorte

Após recortar, a implementação calcula 255−I, remove a tendência linear, remove a média residual e divide pelo desvio padrão do próprio segmento. Isso produz os sinais preparados x e y. O ajuste de tendência do panorama inteiro não é reaproveitado no segmento. Não há redimensionamento do eixo, deformação temporal, janela de Hann ou alinhamento otimizado por WCO.

As sensibilidades previamente definidas repetem a avaliação com outra faixa vertical, remoção somente da média e deslocamentos de ±16 px dos centros. Os parceiros permanecem os mesmos; nenhuma variante substitui retrospectivamente o resultado principal por apresentar maior coerência.

## Como a autocorrelação é calculada

Para os sinais reais desta análise, a autocorrelação de x compara o recorte com uma cópia deslocada dele próprio:

$$R_{xx}[m]=\frac{1}{N}\sum_{\substack{0\leq n<N\\0\leq n+m<N}}x[n+m]x[n],\qquad N=2048.$$

Ryy é calculada da mesma forma usando somente y. A correlação cruzada é:

$$R_{xy}[m]=\frac{1}{N}\sum_{\substack{0\leq n<N\\0\leq n+m<N}}x[n+m]y[n].$$

Em cada lag m, entram somente os produtos entre amostras que se sobrepõem. A correlação é **linear, sem retorno circular entre as extremidades**. Para os sinais reais, Rxx e Ryy são simétricas em relação ao lag zero; Rxy não precisa ser simétrica.

| Lag m | Produtos sobrepostos | Divisor utilizado |
| --- | --- | --- |
| 0 | 2.048 | 2.048 |
| +100 ou −100 | 1.948 | 2.048 |
| +2.047 ou −2.047 | 1 | 2.048 |

O divisor permanece N, mesmo quando a sobreposição diminui. Esse é o estimador **enviesado 1/N**, adaptação explicitamente adotada no projeto; o artigo de Schuck usa N−|m| nas Eqs. 32.21–23. Os pontos de grandes lags têm pouca sobreposição e não ganham confiabilidade experimental por serem calculáveis.

Com a padronização utilizada, Rxx[0] e Ryy[0] ficam aproximadamente iguais a 1. Isso não obriga a WCO entre eventos diferentes a ser 1. O denominador da WCO usa os **coeficientes wavelet das autocorrelações no lag zero**, que dependem da escala; não usa simplesmente Rxx[0]Ryy[0].

A convenção cruzada é a da soma x[n+m]y[n]: se y estiver atrasado 100 amostras em relação a x, o controle de atraso conhecido produz o pico em **−100 px**.

## O que a CWT recebe e quais são seus eixos

Cada correlação completa tem **2N−1=4.095 amostras**, correspondentes aos lags inteiros de **−2.047 a +2.047 px**. O lag zero ocupa o índice **2.047** do vetor da correlação. Esse índice é diferente do índice **1.024** da âncora no sinal recortado.

Com uma realização por par, K=1, calculam-se:

$$P_x^{alt}(b,a)=|\mathrm{CWT}\{R_{xx}\}(b,a)|,\qquad P_y^{alt}(b,a)=|\mathrm{CWT}\{R_{yy}\}(b,a)|,$$

$$C_{xy}^{alt}(b,a)=\mathrm{CWT}\{R_{xy}\}(b,a),$$

$$\Gamma_{xy}^{alt}(b,a)=\frac{|C_{xy}^{alt}(b,a)|^2}{P_x^{alt}(0,a)P_y^{alt}(0,a)}.$$

A WPS alternativa usa módulo, sem quadrado. A WCS alternativa é complexa; seus gráficos mostram o módulo. A WCO tem numerador ao quadrado e denominador avaliado no lag zero para cada escala, sem clipping. Denominadores numericamente inválidos são marcados como ausentes, não como coerência zero.

**O eixo b dos mapas representa lag da correlação.** Ele não representa uma coluna do panorama nem a posição de uma estrutura dentro do recorte. A Morlet usa ω0=6; a conversão entre escala a e comprimento de onda é λ=4πa/(6+√38).

| Grandeza | Limite/domínio nesta execução | Significado |
| --- | --- | --- |
| Evento visual | [L, R), extensão de 96 a 1.536 px | Intervalo anotado pela regra óptica |
| Sinal x ou y | 2.048 amostras; posição relativa −1.024…+1.023 px | Evento e contexto ao redor da própria âncora |
| Correlação completa | 4.095 lags; −2.047…+2.047 px | Deslocamentos entre cópias dos sinais |
| Comprimentos de onda | 200 valores logarítmicos de 50 a 1.024 px | Escalas examinadas pela CWT |
| Resumo principal da WCO | 65 lags; −32…+32 px | Média por escala e depois entre escalas válidas |

Os **±32 px do resumo são deslocamentos da correlação**. Não significam recortar apenas 65 pixels da imagem ao redor da âncora. Todos os dados do segmento participam das correlações; a CWT também combina lags vizinhos segundo a escala. O resumo principal não procura o máximo sobre todos os deslocamentos.

O domínio principal satisfaz a regra geométrica de suporte declarada no protocolo. A média com suporte conservador é um diagnóstico separado, em outro domínio. Essas regras e o preenchimento computacional usado na FFT não criam amostras observadas, não aumentam o tamanho experimental e não constituem testes de significância.

## O que a âncora muda e as limitações

A âncora determina quais amostras entram no recorte e coloca os eventos em posições relativas próximas antes de calcular a correlação cruzada. Apenas mudar o nome da origem de um mesmo recorte não altera sua autocorrelação. Mover o recorte pode alterá-la, pois muda o conteúdo e o contexto disponíveis.

A análise caracteriza o **segmento contextualizado**. Em 22 dos 26 segmentos aceitos, há pelo menos outro evento aceito dentro do contexto; estruturas menores ou rejeitadas também podem estar presentes. Há sobreposição entre segmentos. Assim, uma WCO alta não pode ser atribuída exclusivamente ao evento central, nem eventos sobrepostos podem ser contados como réplicas independentes. A agregação é feita primeiro por par de aquisições.

Os limites por intensidade são operacionais e não constituem interfaces físicas calibradas. A mudança de limiar de 60 para 55/65 deslocou uma fronteira em até 40 px. As perturbações de centro de ±16 px não cobrem toda essa incerteza. Também permanecem as costuras dos panoramas e a ausência de calibração espacial física.

## Extensão executada: churn e anular ondulado

Foram localizados **76 mínimos principais e 76 máximos de controle nas 19 aquisições elegíveis**. O marco óptico é um extremo local da média das linhas 18–71, suavizada com gaussiana σ=16 px somente para localização. A proeminência e a largura foram limitadas previamente, e quatro estratos espaciais distribuíram os marcos ao longo do panorama. A seleção não consultou WCO. As 19 pranchas do catálogo foram revisadas visualmente; os 76 mínimos mantiveram correspondência local em σ=12/20, com deslocamento máximo de 11 px. Essa correspondência busca o candidato próximo do centro original; não repete a seleção completa por estrato.

**A suavização não alimenta a CWT.** Depois de achar c, o recorte vem da média original das linhas 18–71. A largura à meia proeminência serve para comparar marcos ópticos; não define os limites de x/y e não equivale a comprimento físico de bolha/onda.

| Padrão e âncora de exemplo | Centro c | Recorte [s,e) | Posição relativa |
| --- | --- | --- | --- |
| Agitado (churn), m03q05_001_MIN02256 | 2.256 | [1.232, 3.280) | −1.024…+1.023 px |
| Anular ondulado, m03q10_001_MIN03119 | 3.119 | [2.095, 4.143) | −1.024…+1.023 px |

Esses são exemplos individuais predefinidos de localização. O [catálogo de pares](results/extensao_churn_anular/pares_ancoras.csv) determina quais âncoras foram efetivamente comparadas. Dois marcos podem ter centros globais diferentes: cada um é centralizado na própria janela, sem deslocamento adicional para maximizar WCO.

### Autocorrelação de um contexto ancorado

Rxx continua comparando as 2.048 amostras de x com uma cópia deslocada de x. Ryy faz o mesmo com y; Rxy cruza os dois recortes. **A âncora não substitui o sinal por um impulso e não entra na soma como peso.** Ela só determina o conteúdo da janela. Rebatizar a origem de uma janela idêntica deixa sua autocorrelação igual; mover a janela troca amostras, tendência estimada e contexto, podendo alterar Rxx/Ryy. A correlação é linear, sem emenda circular.

O marco central é apenas uma parte do conteúdo: 120 segmentos contêm outro marco selecionado e 143 sobrepõem outro segmento selecionado. Portanto, a WCO caracteriza a região ancorada inteira. O ponto de maior coerência no mapa não localiza automaticamente a onda/bolha no panorama.

### Polaridade e limites dos resultados

Trocar y por −y inverte Rxy e WCS, mas deixa Ryy e |WCS|² iguais; a WCO permanece igual. Assim, um controle mínimo × máximo pode apresentar coerência alta. Além disso, ele pode alterar largura e forma, não somente polaridade. Esses controles são diagnósticos, sem obrigação matemática de coerência baixa.

Os 830 pares e oito variantes produziram 6625 comparações calculáveis e 15 ausências ópticas explícitas em σ=20. A ancoragem mostrou viabilidade operacional, mas o ganho médio principal sobre janelas fixas foi pequeno em churn e quase nulo no anular. Não foi comprovada identidade física, superioridade geral ou capacidade de classificar os padrões.

Veja o [resumo com gráficos](EXECUTIVE_SUMMARY.md), o [protocolo específico](PROTOCOL_CHURN_ANULAR.md), a [estabilidade das âncoras](results/extensao_churn_anular/estabilidade_ancoras.csv) e as [coordenadas completas](results/extensao_churn_anular/catalogo_ancoras.csv). A [versão anterior desta nota](ANCORAS_E_AUTOCORRELACAO_ETAPA1.md) permanece como registro histórico da proposta.

## Fontes dentro da entrega

- [Protocolo congelado](PROTOCOL.md): regras anteriores aos novos cálculos.
- [Catálogo](results/catalogo_eventos.csv): limites dos eventos e recortes, âncoras, padrão e contexto.
- [Observabilidade da etapa 1](results/observabilidade.csv): inclusões e limitações antes da extensão; para a situação atual, veja [cobertura da extensão](results/extensao_churn_anular/cobertura_ancoras.csv).
- [Detecção e coordenadas](scripts/catalog_events.py): construção de âncora e segmento.
- [Preparo dos recortes e resumo](scripts/analysis_common.py): corte antes de preparar e domínio de lag.
- [Núcleo matemático](scripts/wavelet_core.py): correlações, CWT e denominador da WCO.
- [Passo a passo executado](PASSO_A_PASSO.md) e [resumo executivo](EXECUTIVE_SUMMARY.md): resultados e limitações.

<!-- AUDITORIA_ADVERSARIAL_INICIO -->
## Ressalvas da revisão adversarial

A correspondência local em σ12/20 não equivale à seleção completa: esta recuperou 66/76 e 68/76 mínimos originais em até 64 px. Um marco pode continuar localizável e perder a posição de maior proeminência no estrato.

A correlação 1/N e a Morlet principal com multiplicador espectral real não negativo permitem justificar WCO≤1 por Cauchy–Schwarz. Esse limite não valida fisicamente a janela nem se transfere automaticamente ao estimador não enviesado ou ao kernel legado. As Eqs. 32.16 e 32.20 do artigo divergem para K>1; aqui K=1 e as médias entre pares são descritivas.

**Fase legada:** refletir somente o eixo de lag da WCS v5 conserva as magnitudes, mas não a fase na convenção principal. A comparação legada desta entrega é de magnitudes e escores; não se devem extrair atrasos físicos dessa fase refletida.

**Eixo do sinal:** c é uma coluna da montagem, sem calibração de deslocamento físico. O v11 duplica uma coluna na operação de feather de 2 px e seu comentário sobre o sentido do deslocamento está invertido. Os panoramas foram reproduzidos exatamente e o efeito de retirar o feather foi pequeno nos contrastes testados, mas a validade física da reconstrução permanece em aberto.

O procedimento completo também produziu ganho médio de 0,0169 em sinais AR(1) independentes. Portanto, localizar extremos e obter WCO relativamente maior não demonstra que se reconheceu o mesmo evento. Veja [evidências e limites da revisão](AUDITORIA_ADVERSARIAL.md).
<!-- AUDITORIA_ADVERSARIAL_FIM -->
