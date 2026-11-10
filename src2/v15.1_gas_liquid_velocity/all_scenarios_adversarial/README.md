# Investigação adicional — todos os cenários no gráfico

Comece pelo [RESUMO EXECUTIVO EM PT-BR](findings/EXECUTIVE_SUMMARY_PT_BR.md) ou pelo [EXECUTIVE SUMMARY](findings/EXECUTIVE_SUMMARY.md). A investigação reuniu 26 pontos ópticos, com 25 pares que passam no método de referência e sete que mantêm suporte e coordenadas compatíveis em todas as sensibilidades. Nenhum método identificou independentemente duas velocidades materiais de gás e líquido.

## Resultados finais

| Arquivo | Conteúdo |
|---|---|
| [Gráfico das 26 gravações — PT-BR](figures/all_26_selected_wall_PT_BR.png) | Núcleo × parede selecionada, com suporte, sensibilidade e saturação visíveis |
| [PDF para publicação](figures/all_26_selected_wall_PT_BR.pdf) | Mesmos pontos em formato vetorial |
| [Cobertura dos métodos](figures/method_coverage_PT_BR.png) | Todas as gravações e critérios específicos de cada alternativa |
| [Dados exatos dos 26 pontos](results/final_all_26_plot_data.csv) | Coordenadas com sinal, seleção, domínio temporal e status |
| [Sensibilidades finais](results/final_sensitivity.csv) | Referência corrigida × troca de metades, perfis brutos e filtragens 32/8 px |
| [Cobertura final](results/final_method_coverage.csv) | Decisões por gravação e método |
| [Métricas finais](results/final_summary_metrics.json) | Contagens e definição da regra conjunta |
| [Legendas em PT-BR](findings/FIGURE_CAPTIONS_PT_BR.md) | Definições e limites para publicação |

As figuras principais têm versões EN e PT_BR, PNG de 600 dpi e PDF/SVG vetoriais. O gráfico mostra módulos em 10³ px/s, preservando os valores com sinal nos CSVs. m03q05_001 permanece um candidato sem cobertura suficiente; m03q40_001 permanece assinalado por saturação. Os pontos vazios não são promovidos a velocidades materiais.

## Escopo e decisões

Foram preservadas as 26 gravações de 501 quadros e as exclusões anteriores: m05q05_002, anular liso e classificações duplas. Os dados originais, v12, v12.1, v14, v15 e resultados anteriores da v15.1 não foram sobrescritos. Esta pasta contém somente a investigação adicional autorizada.

O método principal usa dez blocos temporais, escolhe uma parede em uma metade e avalia a escolha na outra. O fundo é estimado apenas na metade de seleção. As médias usam somente blocos conjuntamente aceitos da metade de avaliação; não representam todos os 500 intervalos. A estabilidade final exige suporte e concordância de ambas as coordenadas nas quatro alternativas, incluindo o passa-altas de 8 px. Os resultados exploratórios anteriores, com fundo compartilhado entre as metades, foram mantidos e não controlam o gráfico final.

[Protocolo e suas extensões documentadas](PROTOCOL.md) · [Revisão adversarial de identificabilidade](findings/adversarial_identifiability_PT_BR.md) · [Revisão regional](findings/block_independent_review.md) · [Método mascarado](findings/phase_local_review.md) · [Método espectral](findings/spectral_review.md).

## Reprodução

Os scripts resolvem seus caminhos a partir da própria localização. As dependências usadas são NumPy, SciPy, Pillow, pandas e Matplotlib, disponíveis no ambiente utilizado. Não há download de dados. Executar a sequência completa a partir da raiz do repositório:

    python3 -B src2/v15.1_gas_liquid_velocity/all_scenarios_adversarial/scripts/run_all.py

Para regenerar somente as figuras a partir dos resultados numéricos existentes:

    python3 -B src2/v15.1_gas_liquid_velocity/all_scenarios_adversarial/scripts/final_figures.py

Para reproduzir a sensibilidade dos critérios de cobertura e tolerância local:

    python3 -B src2/v15.1_gas_liquid_velocity/all_scenarios_adversarial/scripts/gate_sensitivity.py

Com tolerância local de max(2 px/quadro, 20%), exigir 1, 2, 3, 4 ou 5 blocos por metade fornece 26, 25, 25, 22 ou 14 pares, respectivamente. Exigir três blocos com tolerância de max(2 px/quadro, 10%) fornece 19. A opção de um bloco completa a contagem, mas não valida a representatividade do intervalo isolado; a referência permanece em três blocos por metade. [Tabela completa](results/reporting_gate_sensitivity_counts.csv).

Os ramos individuais são reproduzidos por block_ncc.py → selected_wall.py → split_background_ncc.py; phase_local_experiment.py → phase_local_challenge.py → phase_local_review.py; e spectral_transport.py → spectral_finish.py. Os controles adicionais são adversarial_identifiability.py e adversarial_new_methods.py. Os perfis e máscaras anteriores são entradas explícitas, preservadas com hashes; não devem ser removidos ao reproduzir o estudo.

## Limites da validação

O [registro final de verificação](results/final_verification.json) reúne as conferências da versão corrigida e dos artefatos finais. Sua reprodução é feita por scripts/final_verify.py.

As verificações conferem aritmética, correspondência óptica em controles construídos, critérios e proveniência. Os dados experimentais foram reutilizados no desenvolvimento e na avaliação; não há conjunto experimental externo, calibração axial confirmada nem traçadores que validem as duas velocidades materiais. Variação entre filtros ou metades não é intervalo de confiança. Passar por critérios empíricos não estabelece exatidão física.

Evidências estruturadas: [verificação independente dos blocos](results/independent_block_verification.json), [auditoria dos novos métodos](results/adversarial_new_methods.json), [verificação mascarada](results/phase_local_verification.json), [verificação espectral](results/spectral_verification.json) e [controle do fundo por metade](results/split_background_metadata.json). A auditoria regional inicial registrou o compartilhamento do fundo; a figura final utiliza a correção documentada no protocolo.
