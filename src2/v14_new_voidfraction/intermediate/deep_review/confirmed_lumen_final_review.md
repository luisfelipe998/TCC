# Revisão final independente — lúmen confirmado

**Resultado:** não foi encontrada divergência numérica ou metodológica nas versões PT-BR/EN regeneradas. Há apenas um detalhe editorial menor identificado abaixo. Esta revisão foi feita após as etapas de resumo e auditoria adversarial, sem executar novamente o pipeline.

## Conferência numérica

- Recalculei diretamente dos intervalos das 26 máscaras a soma de áreas circulares com D=53. Os resultados coincidiram com a tabela por aquisição.
- As 26 máscaras são idênticas, byte a byte, às anteriores à correção D53. Todos os campos comuns dos 26 registros de flip-rate permanecem iguais.
- A reescala alfa53 = alfa60 × (60/53)² e a reescala da fração projetada por 60/53 conferem; maior diferença de ponto flutuante foi 3,33 × 10⁻¹⁶.
- As contagens foram recalculadas da diferença entre pixels vizinhos e os denominadores conferem com 53(W−1).
- As cinco linhas das tabelas dos dois resumos coincidem com as medianas recalculadas, incluindo formatação decimal PT-BR/EN. São 25 aquisições elegíveis para estatísticas, dentro dos 26 registros selecionados.
- Os campos de limites de extensão e D56/D64 não aparecem nos resultados principais. Cada registro contém linhas 18:71, D53 e proveniência da confirmação pelo usuário.

## Interpretação e solicitação mais recente

Os dois resumos agora incluem as três explicações solicitadas: possíveis erros de fase ao converter luz/sombra em gás; soma das áreas sob a hipótese circular, com analogia de visão lateral versus profundidade e refração; e histórico do recorte anterior.

O intervalo antigo de m06q80_001, **0,6826–0,8575**, foi conferido contra o JSON arquivado e é explicitamente identificado como resultado da geometria anterior D60, com recortes 50/53/57/60. Não é apresentado como incerteza atual ou como evidência de gás fora do interior confirmado.

O intervalo atual de remoção de linhas internas, **0,8748–0,9637**, confere com o JSON atual e mantém D53 para recortes de 50/51/52/53 linhas. A confirmação dos limites, a hipótese circular e a calibração de refração estão diferenciadas. O valor máximo atual é 1; máscara preenchida não é apresentado como comprovação de gás puro.

A busca nos textos ativos e na navegação encontrou D60 e 0,7803 somente em comparações explicitamente históricas. O relatório exploratório de alternativas D60 está identificado como legado na abertura. Não foi encontrada alegação atual de gás perdido fora do lúmen completo.

Os arquivos de auditoria pós-resumo, normalização e escopo registravam PASS, sem falhas, no fim desta leitura. A revisão gráfica ficou com o agente principal.

## Pequeno ajuste editorial observado

Em findings/adversarial_review.md, no parágrafo abaixo de “Removal of known-interior rows”, aparece **“Refraction mapping remains uncalibrated.. Every variant”**. O texto de BOUNDARY_PROVENANCE já termina em ponto e o gerador adiciona outro. Basta remover o ponto extra no formato de adversarial_review.py ou usar BOUNDARY_PROVENANCE.rstrip('.'). Isso não afeta dados, métodos ou equações. Não alterei os arquivos revisados nesta etapa de leitura.

