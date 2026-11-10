# v15.2 — Investigação de velocidades relativas de gás e líquido

**Comece pelo [EXECUTIVE SUMMARY em PT-BR](findings/EXECUTIVE_SUMMARY_PT_BR.md) ou pelo [English executive summary](findings/EXECUTIVE_SUMMARY.md).** Esses são os resumos canônicos desta versão. A v15.1 foi preservada.

A investigação processou os 26 cenários elegíveis e os 501 quadros nativos de cada um. Recuperou movimentos de objetos e envelopes candidatos a gás; não identificou dois movimentos materiais independentes de líquido e gás. Campos ausentes não foram preenchidos. Não interpretar razões ópticas como escorregamento de fases.

## Resultados para consulta

- [Estado final dos 26 cenários](results/final_per_run_status.csv).
- [Cobertura completa, figura PT-BR](figures/coverage_all_scenarios_PT_BR.png).
- [Comparações nos mesmos intervalos](results/paired_object_consensus.csv) e [figura](figures/paired_objects_PT_BR.png).
- [Comparações entre gravações](results/inter_recording_object_ratios.csv).
- [Legendas para figuras](findings/FIGURE_CAPTIONS_PT_BR.md).
- [Protocolo, critérios e desenvolvimento das alternativas](findings/PROTOCOL.md).
- [Revisão física e literatura](findings/physical_identifiability_review_PT_BR.md), [limites condicionais](findings/conditional_piston_PT_BR.md), [experimento de resíduos](findings/residual_review_PT_BR.md).
- [Auditoria independente](findings/object_independent_review_PT_BR.md) e [verificação final](results/final_verification.json).

## Reprodução

Executar na raiz do repositório:

    python3 -B src2/v15.2_gas_liquid_relative_velocity/scripts/run_all.py

Para conferir os artefatos existentes, seus hashes, escopo, números e links sem refazer a detecção e a correlação:

    python3 -B src2/v15.2_gas_liquid_relative_velocity/scripts/final_verify.py

Dependências: Python, NumPy, SciPy, pandas, Pillow e Matplotlib. As versões utilizadas são registradas na verificação final. Não há dependência de OpenCV, rede, novos dados experimentais ou serviços externos. Figuras principais têm PNG de 300 dpi e versões vetoriais PDF/SVG. As ilustrações diagnósticas de auditoria podem ter outra resolução de exportação.

Entradas preservadas: data_processed, inventários e metadados da v15.1, máscaras e referência radial da v15, e detector v14. Esses arquivos continuam necessários para reproduzir do zero; esta pasta não duplica as 13.026 imagens. Os scripts novos escrevem nesta investigação. Cache de objetos verifica hashes das imagens, detector, referência e funções de extração antes de reutilizar detecções.

A pasta intermediate contém detecções por quadro, nós das trajetórias e logs. Results preserva estimativas rejeitadas e critérios, além das tabelas selecionadas. Exemplos sintéticos permanecem separados dos resultados experimentais.

## Interpretação obrigatória

Velocidade em px/s é deslocamento na imagem dividido pelo tempo nativo. Um traçador material, uma interface, uma onda e uma sombra não são observáveis equivalentes. A conversão métrica não foi usada porque falta calibração axial independente. Razões dentro de uma gravação cancelam uma escala comum; entre gravações, a interpretação física exige magnificação relativa conhecida. Os registros são curtos e as seleções de trajetórias não representam automaticamente todo o escoamento ou uma população experimental de padrões.
