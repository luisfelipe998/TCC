# Revisão adversarial da segmentação e das equações

**Resultado: erro real de borda corrigido e verificado.** Os valores anteriores estão preservados no histórico. Os controles verificam o algoritmo; não fornecem calibração física de gás/líquido.

## Erro de fechamento morfológico

O fechamento 3×3 era aplicado diretamente ao recorte. Na erosão, o valor zero fora da matriz apagava a primeira e a última linha. A etapa seguinte tentava rejeitar componentes nessas linhas já apagadas; assim, nunca rejeitava um componente de primeiro plano por contato radial.

Correção: margem de um pixel com zeros → dilatação → erosão → retirada da margem → identificação de componentes → rejeição dos que tocam uma borda radial. Contatos axiais continuam permitidos. O tamanho mínimo de quatro pixels e o preenchimento de buracos foram mantidos.

Verificação independente: 10 matrizes compararam laços explícitos, operações booleanas e a implementação corrigida. Controles adicionais cobrem contato radial/axial, buracos, iluminação uniforme, ruído limitado, perfil líquido e interior claro de gás. Todos passaram.

24 de 26 máscaras mudaram. Maior variação absoluta de α por aquisição: 0.005417. A tabela usa as 25 aquisições elegíveis.

| Padrão | n | α antes | α corrigido | Trocas/1.000 antes | Trocas/1.000 corrigido |
|---|---:|---:|---:|---:|---:|
| dispersed bubbles | 1 | 0.001836 | 0.001380 | 3.202 | 2.464 |
| plug | 3 | 0.465548 | 0.462812 | 10.045 | 9.730 |
| slug | 3 | 0.475613 | 0.473078 | 9.750 | 10.419 |
| churn | 10 | 0.645083 | 0.641873 | 9.980 | 11.190 |
| wavy annular | 8 | 0.767531 | 0.767346 | 0.963 | 1.133 |

São medianas por padrão. [Todas as alterações por aquisição](../results/validation/deep_method_changes.csv).

**Limite:** cumprir o critério de borda não prova que os componentes restantes sejam gás. Sombras do líquido podem permanecer; bolhas reais conectadas à borda podem ser rejeitadas. Os componentes aceitos não precisam conter um buraco fechado.

## Equações e limites

Reconstrução: α̂ = (1/W) Σₓ Σᵢ (dᵢ(x)/Dₐₚₚ)², com Dₐₚₚ = 60 pixels. O quadrado vem da razão entre áreas circulares. Cada intervalo radial separado é tratado como um disco. Máscaras vazia, completa e com dois intervalos separados verificaram essa soma dos quadrados.

Flip-rate: F = ΣᵧΣₓ₌₂ᵂ 1[M(y,x)≠M(y,x−1)] / [53(W−1)]. Multiplica-se por 1.000 apenas para apresentação. São transições espaciais. A forma tridimensional e a identidade das fases continuam sendo hipóteses da máscara.

Para 53 linhas e Dₐₚₚ = 60, o maior valor que a fórmula pode produzir é (53/60)² = 0,7803. É uma consequência matemática da configuração, não um limite físico do escoamento. Retirar essa linha do gráfico não muda o algoritmo.

## Cenário de extensão: não é um limite geral de gás

U = média[α_visível + 1_borda(1−α_visível)] só expande colunas cujo gás detectado toca o recorte. Uma estrutura circular de 3 pixels inteiramente nas linhas omitidas pode ter α = 0,0025, enquanto a máscara visível e U permanecem zero. O controle analítico confirmou esse contraexemplo.

Os campos antigos de limite inferior/superior descrevem um cenário de extensão. Só são limites se não existir gás totalmente oculto, rejeitado ou perdido, além das hipóteses de geometria e fase. Não são intervalos de confiança.

## Sonda e independência

A sonda denominada “centro” era y=40. O ponto médio entre mínimos aparentes y=15 e y=75 é y=45. Mantém-se y=40 por compatibilidade, identificado como sonda fixa. Nenhuma das linhas valida o centro físico. A tabela por aquisição inclui a sensibilidade y=40 versus y=45.

A imagem de bolhas ficou fora da construção da referência, mas foi inspecionada no desenvolvimento do detector. É um controle interno. Os trechos líquidos inferidos no mesmo conjunto também não são verdade física independente.

Detalhes: [verificações completas](../results/validation/deep_method_checks.json) e [código](../scripts/deep_method_checks.py).
