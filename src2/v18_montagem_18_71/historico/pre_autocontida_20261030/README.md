# Execução com montagem nas linhas 18–71

Esta execução reconstrói os panoramas e recalcula os resultados dependentes com perfis de alinhamento nas 54 linhas de 18 a 71. As versões anteriores permanecem preservadas.

1. `python3 -B src2/v18_montagem_18_71/scripts/reconstruir.py`
2. `python3 -B src2/v18_montagem_18_71/scripts/recalcular_wavelet.py --workers 2`
3. `python3 -B src2/v18_montagem_18_71/scripts/process_optical.py`

4. `python3 -B src2/v18_montagem_18_71/04_analise_simplificada/scripts/analisar.py`

Os geradores das figuras e a finalização dos exemplos ópticos ficam documentados nos manifestos e nos scripts das respectivas etapas. Os diretórios `01_panoramas` e `03_optica` contêm os panoramas e as medidas ópticas vigentes. A análise estatística atual, com Pearson diretamente sobre a potência no atraso zero, está em [`04_analise_simplificada/`](04_analise_simplificada/README.md). A pasta `02_wavelet` preserva os sinais, as potências e os campos wavelet utilizados como entrada, além da comparação anterior em resíduos. Seus coeficientes de Pearson e a seleção antiga são históricos. A análise óptica atual usa exclusivamente `18:72` em normalização, referência, suavização, segmentação, saturação e medidas. A revisão e o estado anterior estão documentados em `03_optica/README.md`.

A simplificação não modificou o preparo dos sinais nem os campos wavelet. Alterou os coeficientes, seus resumos e os exemplos intrapadrão selecionados. A pasta atual registra também a associação dos exemplos às medidas ópticas vigentes. A atualização apenas dessa junção pode ser executada por `04_analise_simplificada/scripts/atualizar_optica.py`, sem recalcular a análise espectral.

Ver a [auditoria da execução](../../auditorias/montagem54/README.md) para critérios, verificações, preservação do histórico e correspondência com o manuscrito.
