# Retomada — extensão de âncoras para churn e anular ondulado

**Estado: PAUSADO A PEDIDO DO USUÁRIO. Não retomar até o usuário solicitar.**

Pasta desta análise: `/Users/i521862/Documents/personal/unisinos/TCC/src2/v16.2_wco_eventos`.
Base científica/de dados: `../v16_spatial2_18_71`. A pasta da base não foi alterada.

O usuário quer concluir a análise sem perder o trabalho, mas precisou interromper a execução. Foi solicitado identificar o próximo checkpoint seguro, parar e salvar o contexto nesta pasta. **O processo de cálculo foi efetivamente encerrado**, por interrupção solicitada, com código de saída 130. Não ficou um lote aguardando para continuar automaticamente.

## Checkpoint efetivamente alcançado

As etapas **anchors** e **fixed** foram concluídas e seus arquivos foram fechados antes da interrupção. A execução foi interrompida **durante diagnostics**, dentro da CWT de um controle, após a gravação dos resultados principais e das referências. O `KeyboardInterrupt` no fim do log registra essa parada intencional; não indica falha dos resultados já salvos.

| Conteúdo preservado | Quantidade / situação |
| --- | --- |
| Catálogo óptico | 1.365 candidatos; 152 âncoras selecionadas |
| Aquisições | 10 churn + 9 anular ondulado |
| Âncoras principais | 76 mínimos: 40 churn + 36 anular |
| Controles ópticos | 76 máximos: 40 churn + 36 anular |
| Pares locais congelados | 830, selecionados entre 8.208 candidatos |
| Avaliações ancoradas | 6.640 linhas: 6.625 calculadas + 15 não estimáveis |
| Ausências | Somente detector σ=20, envolvendo três máximos sem correspondência óptica |
| Referências fixas | 2.052 linhas: 513 principais reutilizadas + 1.539 novas sensibilidades |
| Perfis principais por escala | 830 pares × 200 comprimentos de onda |
| Mapas predefinidos completos | 12 arquivos NPZ, com WPS/WCS/WCO e sinais |
| Diagnósticos self/mesma aquisição | Interrompidos; não foi gravado o CSV dessa etapa |
| Relatórios integrados da extensão | Scripts preparados; geração e revisão final ainda pendentes |

Os CSVs foram conferidos após a parada: quantidades, unicidade de par/variante, escores finitos onde calculáveis, 200 escalas válidas e 15 ausências explícitas. As quatro variantes fixas têm 513 resultados cada. Os congelamentos da etapa 1 e da extensão foram conferidos por hash. **Não foi perdida nenhuma tabela já gravada.** Os controles diagnósticos que estavam apenas em memória precisarão ser refeitos; não será necessário repetir os pares ancorados nem as referências fixas.

Registro de execução: [execucao.json](results/extensao_churn_anular/execucao.json), agora com estado `PAUSADO_PELO_USUARIO`. O registro conserva os tempos das etapas concluídas: **1.507,848 s em anchors** e **660,615 s em fixed**, total de aproximadamente **36,14 minutos de etapas concluídas**. O tempo parcial de diagnostics não foi cronometrado separadamente no arquivo. Não incluir o período de pausa no futuro tempo de cálculo.

O [checkpoint de hashes](results/extensao_churn_anular/checkpoint_pausa.json) registra os arquivos preservados. A conferência desse checkpoint não substitui a auditoria final científica/computacional, que continua pendente.

## Pedido original que deve ser concluído

Estender o conceito de âncoras para **agitado (churn)** e **anular ondulado (wavy annular)**, aplicando a análise wavelet alternativa na mesma pasta. O foco principal é determinar âncoras e analisar seus contextos. Entregar **resumo executivo com gráficos, tudo em PT-BR, passo a passo e limitações**, preservando o trabalho anterior de plug/slug.

Não trabalhar em `v16.1_perfis_wps`. O nome definitivo desta pasta já é `v16.2_wco_eventos`; não renomear novamente nem apagar a análise. Os antigos caminhos `v16.2` e `wco_eventos` já tinham sido removidos em etapa anterior.

## Método já congelado — não ajustar olhando os resultados

- Entrada: média das linhas **18–71 inclusive**, extraída dos mesmos panoramas auditados da base; 1 px/amostra.
- Detector: gaussiana **σ=16 px somente para localizar extremos**. `find_peaks` em −I para mínimos e +I para máximos; proeminência ≥max(3; 0,25×(P95−P5)); janela de proeminência 1.025 px; largura à meia proeminência entre 16 e 512 px; distância entre candidatos do mesmo tipo ≥128 px.
- Seleção: quatro estratos espaciais; maior proeminência por tipo/estrato, desempate pela menor coluna, separação ≥512 px dos marcos anteriores do mesmo tipo. Até quatro mínimos e quatro máximos por aquisição. Todos os candidatos permanecem catalogados com motivos de não seleção.
- Recorte: **[c−1024,c+1024)**, 2.048 amostras, com 16 px de margem na aquisição para sensibilidades. A wavelet recebe a média original, sem a gaussiana do detector, recortada antes de inverter 255−I, remover tendência/média e padronizar.
- Pareamento: mínimos compatíveis têm razão de larguras ≤1,50; distintos, ≥2,00; intermediários ficam catalogados. Mínimo × máximo é controle óptico. Até dois pares por combinação de aquisições/grupo, escolhidos por hash com semente 20261025; um por orientação no controle de polaridade. Nenhum parceiro foi escolhido pela WCO.
- Núcleo: idêntico à base; Morlet ω0=6, K=1, correlação linear enviesada 1/N, WPS alternativa em módulo, WCS complexa, WCO com numerador ao quadrado e denominador no lag zero. Sem clipping.
- Domínio: 200 comprimentos de onda logarítmicos de 50–1.024 px; correlações com 4.095 lags, −2.047…+2.047; resumo em lags −32…+32, média uniforme por lag e depois por escala válida. A média conservadora é um diagnóstico separado.
- Oito variantes: principal; faixa30_49; somente_media; centros_menos_mais16; centros_mais_menos16; detector_sigma12; detector_sigma20; v5_mesmos_segmentos.
- Detectores σ12/20: buscar o extremo elegível do mesmo tipo mais próximo, até 64 px; sem correspondência, resultado ausente, sem substituição. Não refazer parceiros. A variante wavelet 30–49 mantém os centros principais; sua redetecção foi apenas diagnóstico.
- Agregação: primeiro por par de aquisições, depois por padrão. Controles pareados nas mesmas aquisições; sensibilidades também nos mesmos IDs de pares locais. Retirada de aquisição remove todos os pares incidentes. Sem valores-p ou réplicas independentes artificiais.

Referências: [protocolo da extensão](PROTOCOL_CHURN_ANULAR.md), [configuração](config_ancoras_textura.json), [congelamento](results/extensao_churn_anular/congelamento.json), [catálogo](results/extensao_churn_anular/catalogo_ancoras.csv) e [pares](results/extensao_churn_anular/pares_ancoras.csv).

## O que foi inspecionado e validado

Foram examinadas as quatro pranchas piloto e as 19 pranchas completas, com todas as 152 âncoras, perfis, recortes do panorama e quadros de origem. A revisão é assistida por IA, não uma validação independente por especialista. 142 âncoras fora das zonas de mistura tiveram pixels da coluna conferidos com o quadro de origem; dez ficam em mistura. Os marcos são **ópticos**, não identidade comprovada de ondas/bolhas.

Os testes numéricos passaram antes do lote: detector em sinal constante e vale sintético, centro, translação e inversão; correlação por soma direta; preparo por regressão independente; atraso conhecido; self; cache exatamente equivalente ao núcleo original nas oito variantes. O atraso +100 px de y produziu lag −100 px. A WCO é invariável a y→−y; a WCS muda de sinal. [Validação numérica](results/extensao_churn_anular/validacao_numerica.json).

O cache LRU de `texture_common.py` tem no máximo 16 autocorrelações transformadas. Ele não altera equações. O núcleo v5 usa a avaliação direta herdada. Não mudar os arquivos usados no lote sem registrar a alteração: `execucao.json` guarda seus hashes, e o verificador exige correspondência.

O orçamento inicial comunicado foi aproximadamente 13–30 minutos para o lote. As últimas variantes ficaram mais lentas que o benchmark; a previsão foi revisada quando o usuário perguntou. Depois o usuário solicitou a pausa. [Benchmark](results/extensao_churn_anular/benchmark.json) e [estimativa revisada](results/extensao_churn_anular/estimativa_revisada.json). Não atribuir a lentidão a uma causa não demonstrada.

## Achados já observados — ainda falta auditoria final

O objetivo operacional de localizar e processar marcos foi atendido no lote principal; isso não demonstra melhoria robusta da WCO nem classificação entre padrões.

| Padrão / contraste principal | Diferença média de WCO | Pares de aquisições |
| --- | --- | --- |
| Churn: compatíveis − fixos | +0,028308 | 45 de 45 possíveis |
| Churn: compatíveis − larguras distintas | +0,020218 | 29 de 45 |
| Churn: compatíveis − mínimo × máximo | +0,022939 | 45 |
| Anular: compatíveis − fixos | +0,004441 | 34 de 36 possíveis |
| Anular: compatíveis − larguras distintas | −0,023809 | 13 de 36 |
| Anular: compatíveis − mínimo × máximo | +0,009088 | 34 |
| Churn × anular: compatíveis − fixos | +0,032085 | 89 de 90 possíveis |

Na faixa 30–49, o contraste frente a fixos foi aproximadamente **+0,002468 em churn** e **+0,003896 em anular**, mantendo os centros principais. O ganho de churn depende da representação vertical. As médias principais estão abaixo da referência descritiva de 0,05; a cobertura do controle de larguras distintas fica abaixo de 2/3 em ambos os padrões. Não transformar esses limiares em testes estatísticos ou físicos.

Os **76 mínimos principais** têm correspondência em σ12 e σ20, com deslocamento máximo de **11 px**. As três ausências em σ20 são máximos: m04q80_001_MAX02378, m06q20_001_MAX07196 e m06q40_001_MAX06103, afetando 15 comparações. Redetectar na faixa 30–49 mantém todos os mínimos, mas desloca mínimos de churn em até 47 px.

Dos 152 segmentos, **143 sobrepõem outro segmento** e **120 contêm outro marco selecionado**. Considerando apenas os 76 mínimos: 70 têm sobreposição, 60 contêm outro marco e três caem em mistura. Cada contexto cruza 38–162 costuras. m03q40_001 mantém 60 deslocamentos não positivos herdados; dois contextos selecionados incluem deslocamentos não positivos. Ainda deve ser consolidada a influência de retirar essa aquisição.

Controle sintético: mínimos de larguras distintas tiveram média 0,4465; com o mesmo fundo periódico, 0,9026. Não atribuir a WCO exclusivamente ao ponto central. Sem calibração para mm/Hz/velocidade, revisão independente, validação externa ou classificador validado.

## Arquivos novos já preparados

Em `scripts/`:

- `pilot_texture_anchors.py`: executado; piloto óptico.
- `catalog_texture_anchors.py`: executado e congelado; catálogo, parceiros e pranchas.
- `texture_common.py`: executado; recortes e cache validado.
- `validate_texture_analysis.py`: executado; validações, sintéticos e benchmark.
- `run_texture_analysis.py`: **anchors e fixed concluídos; diagnostics interrompido**.
- `summarize_texture_analysis.py`: preparado, ainda não executado no conjunto completo.
- `render_texture_results.py`: preparado. Somente cobertura/estabilidade e exemplos individuais de âncoras foram gerados preliminarmente e inspecionados. Regerar as figuras finais depois da agregação.
- `write_texture_reports.py`: preparado, ainda não executado; integra resumo, passo a passo, resultados, nota matemática, README e PDF.
- `verify_texture_extension.py`: preparado, ainda não executado; auditoria da extensão.
- `verify_final.py`: adaptado para chamar a auditoria da extensão quando a execução estiver `PASS`, conferir o PDF original arquivado e gerar o manifesto integrado.
- `write_reports.py`: adaptado para reaplicar a integração quando houver resultados completos da extensão; não é necessário reexecutá-lo nesta retomada.

Os scripts novos passaram por compilação sintática durante o desenvolvimento, mas **a sequência completa de pós-processamento ainda não foi exercitada**. Corrigir eventuais problemas de apresentação ou integração sem alterar o catálogo congelado ou os resultados numéricos.

As figuras finais mostrarão padrão e IDs de x/y, sinais, correlações, WPS/WCS/WCO, contrastes, sensibilidades e limitações. O renderizador amostra um lag a cada 4 px fora de ±64 apenas para exibir os PNGs; os cálculos e NPZ completos não são reduzidos. O gerador documenta essa distinção.

## Retomada segura — somente depois de novo pedido do usuário

1. Ler este arquivo e o checkpoint. Conferir seus hashes, o estado `PAUSADO_PELO_USUARIO` e os dois congelamentos. Confirmar que não existe processo remanescente; o processo desta sessão foi encerrado com código 130.
2. Usar **`/Users/i521862/.pyenv/versions/3.13.2/bin/python3`**. O Python padrão não tinha as dependências. Cwd: `/Users/i521862/Documents/personal/unisinos/TCC`. Usar `MPLCONFIGDIR=/tmp/wco-eventos-mpl`.
3. **Refazer somente `diagnostics(signals)` de `run_texture_analysis.py`**, que produz 152 selfs e 19 pares da mesma aquisição. `load_signals()` está em `analysis_common.py`. Verificar os congelamentos antes/depois. Não rodar `anchors`, `fixed` ou o lote completo.
4. **Preservar e completar `execucao.json`.** Atenção: executar diretamente `run_texture_analysis.py --stage diagnostics` recria os metadados apenas dessa etapa e perderia o registro de tempos anteriores. Usar um pequeno orquestrador que leia o JSON existente, chame a função diagnostics, acrescente `etapas.diagnostics`, registre retomada/fim e marque `PASS` somente quando essa etapa e os congelamentos passarem. Manter os hashes de scripts e os tempos de anchors/fixed. Somar os tempos ativos das etapas concluídas; documentar que o tempo parcial interrompido não está separado. Não contabilizar a pausa como cálculo. Preservar uma cópia do JSON pausado antes dessa atualização.
5. Executar sequencialmente, verificando cada retorno:

```sh
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/summarize_texture_analysis.py
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/render_texture_results.py
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/write_texture_reports.py
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/verify_final.py
```

6. Inspecionar visualmente as figuras finais e páginas do PDF, corrigir apresentação se necessário e repetir apenas as verificações afetadas. `pdftoppm`, `pdftotext` e `pdfunite` estão disponíveis em `/opt/homebrew/bin`. A geração planejada produz um PDF integrado com a extensão primeiro e o resumo original de sete páginas após um separador histórico.
7. Confirmar a preservação dos resultados/figuras da etapa 1, fontes, congelamentos e núcleo; conferir exemplos, agregações, conjuntos comuns, ausências, recálculos amostrais, links e hashes finais. **Não declarar a entrega completa enquanto essa auditoria estiver pendente.** Atualizar o estado deste arquivo para concluído, preservando o histórico da pausa, e regenerar o manifesto depois da última alteração documental.
8. Entregar ao usuário links para resumo executivo, PDF, passo a passo e explicação matemática atualizados, com a conclusão honesta: viabilidade de âncoras demonstrada no conjunto, ganho de WCO pequeno/não robusto e limitações explícitas.

## Preservação da etapa anterior e cuidados com o estado atual

Os arquivos `EXECUTIVE_SUMMARY_ETAPA1.md` e `RESUMO_EXECUTIVO_ETAPA1.pdf` já foram preservados. O gerador da extensão arquivará também o passo a passo, resultados, README e nota matemática anteriores antes de integrar os novos documentos.

**Neste checkpoint, os documentos principais da raiz ainda descrevem a etapa 1; não são a entrega final da extensão.** A nota `ANCORAS_E_AUTOCORRELACAO.md` ainda contém a proposta não executada em seu texto antigo; o gerador preparado substituirá essa seção pela extensão efetivamente executada. O resumo/PDF novos não foram gerados.

`DECISIONS.md` e o cabeçalho de `PLAN.md` já receberam a decisão de ampliar o escopo. O plano original continua preservado integralmente; não alterar `PLAN_ORIGINAL.md` nem os arquivos congelados.

`results/verificacao_final.json` e `results/manifesto.json` da raiz ainda são da entrega anterior e **não validam o estado completo desta extensão**. O estado anterior está preservado em [estado_anterior.json](results/extensao_churn_anular/estado_anterior.json).

Não executar o verificador integrado enquanto diagnostics/relatórios estiverem pendentes e interpretar um eventual resultado da etapa 1 como conclusão da extensão. Não modificar os congelamentos para contornar divergências. Os dados e análises anteriores devem continuar intactos.

`src2` permanece ignorado pelo Git. Não houve commit ou sincronização. Alterações preexistentes fora desta análise, em assets e revisões de texto da raiz, foram deixadas intactas. Não há solicitação para apagá-las, publicá-las ou enviar mensagens externas.
