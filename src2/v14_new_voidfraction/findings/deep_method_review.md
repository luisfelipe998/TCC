# Revisão adversarial da segmentação e das equações

**Resultado: erro real de borda corrigido e verificado.** Os valores anteriores estão preservados no histórico. Os controles verificam o algoritmo; não fornecem calibração física de gás/líquido.

## Erro de fechamento morfológico

O fechamento 3×3 era aplicado diretamente ao recorte. Na erosão, o valor zero fora da matriz apagava a primeira e a última linha. A etapa seguinte tentava rejeitar componentes nessas linhas já apagadas; assim, nunca rejeitava um componente de primeiro plano por contato radial.

Correção: margem de um pixel com zeros → dilatação → erosão → retirada da margem → identificação de componentes → rejeição dos que tocam uma borda radial. Contatos axiais continuam permitidos. O tamanho mínimo de quatro pixels e o preenchimento de buracos foram mantidos.

Verificação independente: 10 matrizes compararam laços explícitos, operações booleanas e a implementação corrigida. Controles adicionais cobrem contato radial/axial, buracos, iluminação uniforme, ruído limitado, perfil líquido e interior claro de gás. Todos passaram.

24 de 26 máscaras mudaram na correção morfológica. Maior variação absoluta de α nessa correção, mantendo D=60 em ambos os lados: 0.005417. A tabela histórica abaixo isola a morfologia e usa as 25 aquisições elegíveis. A atualização atual para D=53 é verificada separadamente.

| Padrão | n | α antes, D=60 | α após morfologia, D=60 | Trocas/1.000 antes | Trocas/1.000 corrigido |
|---|---:|---:|---:|---:|---:|
| dispersed bubbles | 1 | 0.001836 | 0.001380 | 3.202 | 2.464 |
| plug | 3 | 0.465548 | 0.462812 | 10.045 | 9.730 |
| slug | 3 | 0.475613 | 0.473078 | 9.750 | 10.419 |
| churn | 10 | 0.645083 | 0.641873 | 9.980 | 11.190 |
| wavy annular | 8 | 0.767531 | 0.767346 | 0.963 | 1.133 |

São medianas por padrão. [Todas as alterações por aquisição](../results/validation/deep_method_changes.csv).

**Limite:** cumprir o critério de borda não prova que os componentes restantes sejam gás. Sombras do líquido podem permanecer; bolhas reais conectadas à borda podem ser rejeitadas. Os componentes aceitos não precisam conter um buraco fechado.

## Equações e limites

Reconstrução atual: α̂ = (1/W) Σₓ Σᵢ (dᵢ(x)/D)², com D = 53 pixels. O usuário confirmou que a fatia Python 18:71 contém todo o interior do tubo: linhas 18 a 70, total de 53 células de pixel. O quadrado vem da razão entre áreas circulares. Cada intervalo radial separado é tratado como um disco. Máscaras vazia, completa e com dois intervalos separados verificaram essa soma dos quadrados.

Flip-rate: F = ΣᵧΣₓ₌₂ᵂ 1[M(y,x)≠M(y,x−1)] / [53(W−1)]. Multiplica-se por 1.000 apenas para apresentação. São transições espaciais. A forma tridimensional e a identidade das fases continuam sendo hipóteses da máscara.

O valor máximo atual é (53/53)² = 1. O antigo valor 0,7803 resultava de dividir o mesmo interior de 53 linhas por D=60 e foi superado pela informação do usuário. Para máscaras idênticas, α₅₃ = α₆₀(60/53)²; os flip-rates permanecem iguais. Uma máscara cheia ainda não prova gás puro: segmentação, refração e geometria circular continuam condicionais.

Não se aplicam extensões de gás para fora da fatia confirmada: essa região fica fora do lúmen definido pelo usuário. Os antigos campos de extensão e as alternativas D=56/64 foram retirados dos resultados principais. A normalização está documentada na [verificação do lúmen](lumen_normalization_review.md).

## Sonda e independência

A sonda permanece na linha fixa y=40. O ponto médio das linhas confirmadas 18 a 70 é y=44; os mínimos ópticos antigos em y=15/75 não definem o lúmen. A tabela inclui a sensibilidade entre y=40 e y=44. A confirmação geométrica do usuário não calibra a refração.

A imagem de bolhas ficou fora da construção da referência, mas foi inspecionada no desenvolvimento do detector. É um controle interno. Os trechos líquidos inferidos no mesmo conjunto também não são verdade física independente.

Detalhes: [verificações completas](../results/validation/deep_method_checks.json) e [código](../scripts/deep_method_checks.py).
