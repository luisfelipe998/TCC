# Recorte interno com diâmetro confirmado de 53 pixels

## Implementado em adversarial_review.py

- Importa a geometria compartilhada de v14_geometry.py: interior completo 18:71, 53 linhas, D=53 e sonda fixa y=40.
- A sensibilidade usa somente (20,70), (19,70), (19,71) e (18,71), mantendo D=53 em todos os casos. Esses recortes removem linhas conhecidas do interior; não representam tubos de diâmetros diferentes.
- As referências líquidas alternativas também calculam alfa explicitamente com D=53.
- Mantém a estrutura crop_variants[str(tuple)] com conditional_alpha, maximum_column_alpha e analytic_crop_ceiling. Acrescenta número de linhas incluídas/removidas e identificação do interior completo.
- O campo analytic_crop_ceiling agora significa somente o máximo matemático após remover linhas conhecidas: (linhas incluídas/53)². Para o interior completo, vale 1.
- Acrescenta proveniência dos limites confirmados e esclarece que a segmentação é recalculada em cada recorte, incluindo as decisões de contorno próximas às bordas artificiais.
- Remove as variantes externas (15,75) e (16,73) do cálculo atual. Nenhuma extensão fora do interior entra nos intervalos ou no resumo.
- O relatório gerado passa a tratar os mínimos ópticos y=15/y=75 como extremos de reflexão, sem usá-los como limites do tubo. O ponto médio das linhas incluídas é y=44; a sonda permanece y=40.
- Mantém os limites físicos: confirmação dos limites pelo usuário não calibra refração, rótulos gás/líquido ou hipótese circular. Máscara preenchida não prova seção fisicamente composta apenas de gás.

## Verificado

O script compila e importa. Foram calculadas 24 reconstruções de recortes internos, em seis aquisições. Todas respeitaram o recorte informado e o máximo matemático com D fixo em 53. Não foi executado o pipeline nem regravado nenhum resultado principal.

## Textos atuais a regenerar ou revisar pelo responsável

Estes eram os pontos ainda presentes nos artefatos gerados antes desta revisão; as mensagens não afirmam que os geradores continuam incorretos após o trabalho paralelo:

- findings/adversarial_review.md:37/46 — minha alteração corrige a tabela e a descrição na próxima execução.
- findings/panorama_void_alternatives.md:11/58/60 — D60, teto 0,7803, y45 e extensão fora do recorte. O gerador pertence ao agente principal.
- findings/EXECUTIVE_SUMMARY.md:37/45/51/56 e findings/EXECUTIVE_SUMMARY_PT_BR.md:37/45/56 — D60, teto e extensão. Geradores atribuídos ao responsável pelos resumos/gráficos.
- findings/silhouette_pattern_review.md:48 — afirma recorte incompleto com D60. O agente principal é responsável pelo gerador nesta revisão.
- findings/deep_method_review.md:29/33/43 — deve distinguir configuração histórica D60 da configuração atual D53 e evitar chamar y45 de centro atual. Responsável: auditoria do método.
- **intermediate/exploratory/image_void_alternatives.md:11** — ainda descreve larguras claras divididas por 60 e reconstruções circulares com D60. Esse arquivo é exploratório e vem de image_void_alternatives.py, fora da minha propriedade. Recomendo identificá-lo explicitamente como diagnóstico legado com normalização histórica, ou revisar seu gerador; não apresentá-lo como a fração de vazio atual com o interior confirmado.

A nota anterior literature_stats_review.md discute a falha do cenário de extensão no contexto anterior à confirmação do interior. Deve ser lida como histórico da revisão; não como proposta atual de gás fora do interior confirmado.

