# Passo a passo da análise executada

## 1. Leitura do plano e atualização da referência

Foi lido o plano original de v16.1 e o material obrigatório da análise anterior. Depois da instrução do usuário, adotou-se **v16_spatial2_18_71** como base. A entrada principal passou a ser a média das linhas 18–71 inclusive, com índices iniciados em zero. A versão inicial do plano foi preservada em [PLAN_ORIGINAL.md](PLAN_ORIGINAL.md).

**Limitação identificada:** o plano citava resultados antigos com máscaras. A base atual contém resultados integrais, sem essas restrições. Esses números históricos não foram usados como resultados da nova execução.

## 2. Conferência das equações

Foram inspecionadas visualmente as páginas 411, 413 e 414 de Schuck Jr. e Bodmann. Foram mantidos: Morlet com ω₀=6; CWT com 1/√a; WPS alternativa como módulo da CWT da autocorrelação; WCS alternativa complexa; WCO alternativa com numerador ao quadrado e denominador no lag zero. A correlação com divisor N foi mantida como adaptação explícita.

**Limitação:** o divisor N difere dos estimadores não enviesados impressos no artigo. Não se utilizou uma média entre aquisições como se fosse um ensemble. A habilidade tcc2-writer foi aplicada à redação; as equações primárias prevaleceram sobre inconsistências de seu material auxiliar.

## 3. Auditoria das fontes

O script [preflight.py](scripts/preflight.py) verificou 412 entradas do manifesto, registrou os hashes de 583 arquivos da base e conferiu os hashes de 31 PNGs e 31 XLSX. Recalculou a média 18–71 de cada panorama e comparou-a exatamente ao sinal utilizado. Conferiu também a média 20–59, os rótulos e as contagens de deslocamentos não positivos. Resultado: [preflight.json](results/preflight.json) e [input_verification.csv](results/input_verification.csv).

**Limitação:** a reconstrução dos panoramas e os deslocamentos não foram recalculados. A reprodução pixel a pixel foi evidência herdada, com procedência verificada. Nenhuma dessas verificações comprovou uma reconstrução física instantânea do escoamento.

## 4. Congelamento do protocolo

Foram registrados [PROTOCOL.md](PROTOCOL.md), [config.json](config.json) e [DECISIONS.md](DECISIONS.md) antes de avaliar a separação dos rótulos desta execução. Fixaram-se a distância de Hellinger, a quadratura trapezoidal em log do comprimento de onda, os controles, os pesos e as exclusões. O hash do protocolo foi conferido pelo script de execução.

O menor registro, com 3.171 pixels, definiu a faixa principal: 187 pontos entre 50 e 2.295,94 px, com suporte geométrico no lag zero. Também foram fixadas as faixas integral de 200 pontos e conservadora de 117 pontos.

**Limitação:** o conjunto já havia sido analisado no v16; o congelamento desta etapa não transformou os dados em teste externo intocado. Nenhum limiar científico de ganho mínimo foi inventado.

## 5. Teste de desempenho antes do lote

Foram medidos três comprimentos: 3.171, 9.586 e 23.079 pixels. A CWT integral levou aproximadamente 0,13–0,24 s; a avaliação somente no lag zero, 0,04–0,07 s; a transformada v5, 0,13–0,40 s. A estimativa de 2–4 minutos para o lote e menos de 1 GB de memória foi comunicada antes da execução. A estimativa de memória foi um orçamento a partir das dimensões e do processamento sequencial, não uma medição de pico residente.

**Limitação:** os tempos dependem do ambiente e da carga. A geração das figuras e a auditoria final ficaram fora da estimativa do lote numérico. Os tempos medidos constam de [preflight.json](results/preflight.json).

## 6. Validação das novas operações

O script [validate_operations.py](scripts/validate_operations.py) verificou a quadratura contra uma integral analítica, soma das massas, identidade, invariância à amplitude, distância máxima para massas disjuntas, rejeição de perfis degenerados, contagens de suporte, ordenação ponderada e recuperação de atraso conhecido. A coluna central calculada por soma inversa de Fourier foi comparada ao núcleo integral em três comprimentos.

**Limitação:** controles numéricos validam operações e interpretação matemática, não o modelo físico experimental. Os resultados estão em [operations_validation.json](results/operations_validation.json).

## 7. Reaproveitamento e cálculo dos perfis

Os perfis compatíveis de registros inteiros e dos segmentos de início, centro e fim foram reaproveitados. Para todas as 26 aquisições, a WPS no lag zero e sua média integral foram recalculadas e comparadas aos perfis herdados. Foram calculadas as médias sob COI e suporte conservador, com contagens de lags, além dos controles 30–49, 20–59, remoção somente da média e discretização v5 nos mesmos valores de escala.

Para os 325 pares foram calculados os perfis de 350 combinações únicas de aquisição/comprimento. O v16 recortou antes do preparo; o comparador literal v5 preparou antes do corte. Os registros brutos, normalizados e por prefixo estão nas tabelas [profiles_raw.csv](results/profiles_raw.csv), [profiles_normalized.csv](results/profiles_normalized.csv) e [pair_prefix_profiles.csv](results/pair_prefix_profiles.csv).

**Limitação:** perfis de registros completos, segmentos fixos e prefixos específicos são observações distintas. A normalização removeu amplitude do perfil e não produziu frações de energia física.

## 8. Comparação com as referências

Foram avaliados todos os 325 pares em 21 representações/referências, incluindo as sensibilidades já fixadas. Os escores WCO vieram das tabelas integrais, preservando máximos por escala, pesos de WPS da referência e piso de 5%. Os valores de distância não foram comparados diretamente aos de WCO: foi usada a concordância de ordenação A, com direção explícita.

Foram apresentados peso por par, peso igual entre padrões/combinações, resultados locais por padrão e todas as combinações cruzadas. O controle com prefixos comparou as representações no mesmo segmento de cada par e na faixa nominal integral.

**Limitação:** as diferenças entre referências históricas 20–59 e atuais 18–71 também envolvem seus pesos de origem. Essa comparação não isolou a extração. Os resultados estão em [all_pairs.csv](results/all_pairs.csv), [ranking_summary.csv](results/ranking_summary.csv), [ranking_by_pattern.csv](results/ranking_by_pattern.csv) e [pair_group_summary.csv](results/pair_group_summary.csv).

## 9. Avaliação da influência e dos fatores de confusão

Foi retirada uma aquisição por vez, removendo todos os seus pares e recalculando os pesos. Foram quantificadas a variação por posição, a distância entre segmentos e registro inteiro e a associação descritiva com a razão dos comprimentos. A aquisição de bolhas dispersas continuou incluída na análise principal.

**Limitação:** a retirada é uma análise de influência, não validação preditiva. As faixas não são intervalos de confiança. Não foram calculados valores-p sobre pares dependentes ou permutações de rótulos sem pressupostos defensáveis. As aquisições anulares são sistematicamente mais longas neste conjunto. As tabelas estão em [acquisition_deletion.csv](results/acquisition_deletion.csv) e [within_acquisition_sensitivity.csv](results/within_acquisition_sensitivity.csv).

## 10. Exemplos e casos desfavoráveis

Foram apresentados panorama, sinal, WPS no lag zero e média em lag dos cinco padrões, usando IDs predeterminados. Foram produzidos mapas de módulo da WCS alternativa e WCO alternativa para quatro pares internos predeterminados e um controle cruzado. Três casos adicionais foram selecionados após as métricas: menor H cruzado, maior H interno e um par com H e escore WCO baixos.

**Limitação:** os três casos adicionais explicam os resultados, mas não os validam. Mapas reduzidos herdados foram usados somente na apresentação; nenhum máximo integral foi recalculado a partir deles. Os critérios e a procedência estão em [selected_cases.csv](results/selected_cases.csv) e [illustrated_maps.csv](results/illustrated_maps.csv).

## 11. Redação, gráficos e verificação final

Foram preparados o [resumo executivo](EXECUTIVE_SUMMARY.md), o [relatório detalhado](REPORT.md), 17 figuras em PNG/PDF e uma versão HTML do resumo com imagens incorporadas. Os números foram lidos das tabelas executadas. A verificação final recalculou distâncias e ordenações de modo independente, conferiu a normalização, a cobertura, os mapas, os links, os arquivos e a preservação das fontes.

**Limitação:** a aprovação de verificações computacionais não garante ausência de toda limitação experimental. O escopo exato de cada verificação está em [final_verification.json](results/final_verification.json); os hashes e versões estão em [artifact_manifest.json](results/artifact_manifest.json).

## 12. Renomeação e reprodução

Por solicitação posterior do usuário, a pasta final recebeu o nome **v16.1_perfis_wps**. A pasta v16_perfis_wps foi movida para esse destino e o link simbólico v16.1 foi removido. Os documentos históricos das análises de origem permaneceram preservados. Scripts e links internos foram verificados no caminho definitivo. Os comandos de reprodução estão em [README.md](README.md).

**Limitação:** src2 continua ignorado pelo Git do projeto. Os artefatos estão no workspace, sem commit e sem mudança da regra de exclusão. As análises de origem, os panoramas e as planilhas foram preservados.
