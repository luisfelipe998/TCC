# Estado final — extensão para churn e anular ondulado

**CONCLUÍDA após a retomada solicitada pelo usuário. Não há lote pendente nem execução aguardando continuação.**

Pasta: `/Users/i521862/Documents/personal/unisinos/TCC/src2/v16.2_wco_eventos`.
Base: `../v16_spatial2_18_71`, preservada. O escopo permaneceu somente nesta análise; `v16.1_perfis_wps` não foi alterada.

O contexto exato da pausa está preservado em [RETOMADA_CHECKPOINT.md](RETOMADA_CHECKPOINT.md). Suas instruções de espera descrevem o estado anterior e foram superadas pelo pedido explícito de retomada. Este arquivo registra o estado atual.

## O que foi retomado e concluído

1. Foram conferidos os 46 hashes do checkpoint, os congelamentos e o núcleo matemático antes de continuar. As tabelas e matrizes já gravadas estavam intactas.
2. Foram preservadas cópias exatas do contexto pausado e dos metadados de execução. [Conferência da retomada](results/extensao_churn_anular/verificacao_checkpoint_retomada.json) e [execução no instante da pausa](results/extensao_churn_anular/execucao_checkpoint_pausa.json).
3. Somente a etapa interrompida de diagnósticos foi refeita: 152 autocomparações e 19 pares da mesma aquisição. Os lotes ancorado e fixo não foram repetidos. [Execução concluída](results/extensao_churn_anular/execucao.json).
4. Foram consolidados contrastes pareados, sensibilidades, cobertura e retirada de aquisições. O painel das oito variantes usa os mesmos 815 pares locais disponíveis em todas elas; as tabelas completas preservam os 830 pares planejados e suas ausências explícitas.
5. Três pequenas ultrapassagens de WCO=1 do núcleo legado v5 foram investigadas por soma direta das correlações e CWT por convolução temporal. Foram reproduzidas e mantidas, sem clipping; todas ocorreram fora dos lags ±32 usados no indicador. A Morlet principal permaneceu ≤1. [Auditoria específica](results/extensao_churn_anular/auditoria_ultrapassagens_v5.json).
6. Foram produzidos e revisados gráficos, resumo executivo, PDF integrado, explicação de âncoras/autocorrelação e passo a passo com limitações, em PT-BR. Os relatórios anteriores foram preservados como etapa 1.
7. A auditoria integrada passou: integridade das fontes, preservação dos resultados anteriores, catálogo reproduzido byte a byte, coordenadas e procedência, recálculos amostrais das oito variantes, mapas integrais, agregações, links e PDF. [Auditoria da extensão](results/extensao_churn_anular/verificacao_final.json), [auditoria integrada](results/verificacao_final.json) e [manifesto](results/manifesto.json).

## Conteúdo entregue

| Item | Estado final |
| --- | --- |
| Aquisições da extensão | 10 churn + 9 anular ondulado |
| Catálogo óptico | 1.365 candidatos; 152 âncoras selecionadas |
| Âncoras principais e controles | 76 mínimos + 76 máximos |
| Pares de âncoras | 830, congelados antes da nova WCO |
| Oito variantes | 6.640 avaliações previstas: 6.625 calculadas + 15 não estimáveis |
| Motivo das ausências | Três máximos sem correspondência óptica em σ=20; sem substituição de centro |
| Referências fixas | 2.052: 513 principais reutilizadas + 1.539 sensibilidades calculadas |
| Diagnósticos concluídos | 152 autocomparações + 19 pares da mesma aquisição |
| Perfis por escala | 830 pares × 200 comprimentos de onda |
| Mapas completos | 16 NPZ; 12 predefinidos + 4 extremos ilustrativos posteriores |
| Resumo em PDF | 20 páginas; inclui as sete páginas originais após um separador histórico |
| Etapa 1 preservada | 103 tabelas, matrizes e figuras protegidas, além dos relatórios históricos |

## Conclusão e limites

Foi possível localizar os mínimos ópticos e aplicar WPS/WCS/WCO alternativas aos contextos ancorados nas 19 aquisições. Os 76 mínimos principais mantiveram correspondência nos detectores σ=12/20, com deslocamento máximo de 11 px.

O ganho médio de WCO frente a janelas fixas foi **0,0283 em churn** e **0,0044 em anular ondulado**, abaixo da referência descritiva de 0,05. No churn, o ganho caiu para **0,0025** ao usar as linhas 30–49 com os centros principais preservados. Não se demonstrou melhora robusta, identidade física dos marcos nem capacidade de classificar os padrões.

As principais limitações são o caráter óptico/exploratório das âncoras, o desenvolvimento no conjunto já conhecido, contextos sobrepostos com estruturas vizinhas e costuras, cobertura menor dos controles de larguras distintas e ausência de validação física independente. Os sinais têm 2.048 amostras em [c−1024,c+1024), enquanto as correlações têm 4.095 lags. A CWT recebe as correlações desses contextos inteiros; não recebe somente o ponto de âncora.

## Onde continuar a leitura

- [Resumo executivo com gráficos e identificação dos padrões](EXECUTIVE_SUMMARY.md).
- [PDF integrado](RESUMO_EXECUTIVO.pdf).
- [Passo a passo com limitações em cada etapa](PASSO_A_PASSO.md).
- [Âncoras, autocorrelação e limites de x/y](ANCORAS_E_AUTOCORRELACAO.md).
- [Resultados detalhados e tabelas](RESULTS.md).
- [Reprodução dos relatórios e conferência da entrega](README.md).
- [Revisão da apresentação](results/extensao_churn_anular/revisao_entrega.json).
- [Registro de decisões](DECISIONS.md).

## Tempo e memória

As etapas concluídas somaram **2.184,56 s (36,4 min)**: 1.507,85 s de pares ancorados, 660,61 s de referências fixas e 16,10 s de diagnósticos na retomada. Esse total exclui a pausa e o tempo parcial dos diagnósticos interrompidos, não medido separadamente. Figuras, redação e auditoria também são etapas separadas. O pico de memória do processo inicial não foi salvo; **399,75 MiB refere-se somente ao processo de retomada**.

## Para trabalhos futuros

Não é necessário executar novamente o script de retomada nem os lotes principais. A entrega está concluída. Para conferir o estado atual, execute `scripts/verify_final.py` com o Python indicado no README, mantendo qualquer redirecionamento de log fora desta pasta. Para regenerar relatórios a partir dos resultados existentes, siga a sequência de agregação/renderização/redação/auditoria do README.

Uma mudança de detector, parâmetros ou recortes constitui nova análise, com resultados atuais já conhecidos. Deve preservar os registros desta execução, declarar o caráter exploratório e usar congelamento e auditoria próprios. A verificação computacional desta entrega não substitui validação física independente.
