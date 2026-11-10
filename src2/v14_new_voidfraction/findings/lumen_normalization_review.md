# Verificação da normalização pelo lúmen confirmado

**PASS:** o usuário confirmou que a fatia Python 18:71 é todo o interior do tubo. São 53 células de pixel, portanto D=53 pixels. A declaração do usuário fornece a geometria de referência; não calibra refração ou forma tridimensional do gás.

As máscaras e as métricas de flip-rate ficaram idênticas às anteriores. A mudança de α é exclusivamente a correção do denominador: α₅₃ = α₆₀(60/53)². A fração projetada muda por 60/53.

Foram verificados 26 arquivos de máscara, 26 conjuntos de métricas de flip-rate e 26 hashes de imagens originais. Maior erro numérico de reescala: 3.33e-16.

| Padrão | n | α anterior, D=60 | α atual, D=53 | Flip-rate atual / 1.000 |
|---|---:|---:|---:|---:|
| dispersed bubbles | 1 | 0.001380 | 0.001768 | 2.464 |
| plug | 3 | 0.462812 | 0.593138 | 9.730 |
| slug | 3 | 0.473078 | 0.606295 | 10.419 |
| churn | 10 | 0.641873 | 0.822621 | 11.190 |
| wavy annular | 8 | 0.767346 | 0.983426 | 1.133 |

O máximo matemático atual é 1. O antigo máximo 0,7803 era consequência do denominador D=60 e foi superado. Valores perto de 1 não provam gás puro: a máscara e a reconstrução circular continuam sujeitas a erro.

Os antigos intervalos de extensão para fora do recorte e as alternativas D=56/64 foram retirados dos resultados principais. A região externa à fatia confirmada não integra o lúmen. Recortes internos menores servem apenas para medir o efeito de descartar linhas conhecidas do interior.

A correção anterior de borda continua separada na [revisão morfológica](deep_method_review.md), com ambos os lados daquela comparação em D=60.

Detalhes: [comparação por aquisição](../results/validation/lumen_normalization_changes.csv), [verificações](../results/validation/lumen_normalization_checks.json) e [código](../scripts/lumen_normalization_checks.py).
