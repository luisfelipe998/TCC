# Resultados detalhados — churn e anular ondulado

Foi possível determinar âncoras ópticas e aplicar WPS, WCS e WCO alternativas em **todas as 19 aquisições elegíveis** de churn e anular ondulado. Foram selecionados **76 mínimos principais e 76 máximos de controle**, antes de calcular a nova WCO. Os 76 mínimos mantiveram correspondência local nas duas suavizações alternativas, com deslocamento máximo de 11 px; esse teste não repete a seleção completa por estrato. **O objetivo operacional foi atendido; não se demonstrou ganho robusto de coerência nem separação entre padrões.** Na análise principal, mínimos compatíveis superaram janelas fixas em 0,0283 no churn e 0,0044 no anular, abaixo da referência descritiva de 0,05.

Os [resultados originais de plug/slug](RESULTS_ETAPA1.md) permanecem separados. As definições de âncora e compatibilidade diferem; não foi feita uma média única entre as etapas.

## Cobertura do pareamento

| Padrão | Grupo | Pares locais | Pares de aquisições | Possíveis |
| --- | --- | --- | --- | --- |
| Agitado (churn) | Mínimos de larguras compatíveis | 90 | 45 | 45 |
| Agitado (churn) | Mínimos de larguras distintas | 46 | 29 | 45 |
| Agitado (churn) | Mínimo × máximo | 90 | 45 | 45 |
| Churn × anular ondulado | Mínimos de larguras compatíveis | 178 | 89 | 90 |
| Churn × anular ondulado | Mínimos de larguras distintas | 79 | 47 | 90 |
| Churn × anular ondulado | Mínimo × máximo | 180 | 90 | 90 |
| Anular ondulado (wavy annular) | Mínimos de larguras compatíveis | 66 | 34 | 36 |
| Anular ondulado (wavy annular) | Mínimos de larguras distintas | 29 | 15 | 36 |
| Anular ondulado (wavy annular) | Mínimo × máximo | 72 | 36 | 36 |

## Extremos ilustrados e intensidade das âncoras

**“Maior/menor escore” ordena a WCO média dos pares de recortes.** Ambos os exemplos usam âncoras em **mínimos locais da intensidade suavizada**, por isso tendem a apontar para regiões escuras. O mínimo é local no perfil médio 18–71 após suavização; não precisa ser o pixel mais escuro da imagem. A seleção usa proeminência e estratos, sem consultar WCO. Após a inversão 255−I, essas depressões tendem a aparecer como picos no sinal preparado.

O escore compara os dois contextos de 2.048 amostras por suas correlações e CWT: média de WCO nas 200 escalas logarítmicas de 50–1.024 px e nos 65 lags de −32 a +32 px. Alinhar dois mínimos não garante que o restante dos sinais seja semelhante.

**Os extremos ilustrados são os dos 90 pares selecionados de churn e 66 de anular**, na variante principal e no grupo de mínimos de larguras compatíveis. A auditoria ampliou o universo e encontrou um churn menor (0,0793) e um anular maior (0,4476). As figuras mantêm os exemplos originais, agora com o universo explícito:

| Padrão | Universo | Pares | Menor WCO média | Maior WCO média |
| --- | --- | --- | --- | --- |
| Agitado (churn) | Amostra original das figuras | 90 | 0,0818 | 0,5007 |
| Agitado (churn) | Todos os admissíveis da auditoria | 484 | 0,0793 | 0,5007 |
| Anular ondulado (wavy annular) | Amostra original das figuras | 66 | 0,0945 | 0,4219 |
| Anular ondulado (wavy annular) | Todos os admissíveis da auditoria | 348 | 0,0945 | 0,4476 |

[Conferência dos IDs, recortes, mínimos ópticos e escores recalculados](results/extensao_churn_anular/verificacao_exemplos.json).

## Contrastes principais

| Padrão | Controle | WCO das âncoras | WCO do controle | Diferença | Pares de aquisições | Fração positiva |
| --- | --- | --- | --- | --- | --- | --- |
| Agitado (churn) | Janelas fixas | 0,2403 | 0,2120 | 0,0283 | 45 | 64,4% |
| Agitado (churn) | Mínimos de larguras distintas | 0,2498 | 0,2296 | 0,0202 | 29 | 55,2% |
| Agitado (churn) | Mínimo × máximo | 0,2403 | 0,2174 | 0,0229 | 45 | 57,8% |
| Anular ondulado (wavy annular) | Janelas fixas | 0,2250 | 0,2205 | 0,0044 | 34 | 55,9% |
| Anular ondulado (wavy annular) | Mínimos de larguras distintas | 0,2113 | 0,2352 | -0,0238 | 13 | 46,2% |
| Anular ondulado (wavy annular) | Mínimo × máximo | 0,2250 | 0,2159 | 0,0091 | 34 | 64,7% |

Somente os mesmos pares de aquisições entram nos dois lados de cada contraste. Por isso a média das âncoras pode mudar entre comparadores. A tabela numérica também preserva mínimos/máximos por par, a fração de diferenças positivas e o contraste no domínio conservador.

## Sensibilidades

| Variante | Churn: âncoras − fixos | Anular: âncoras − fixos |
| --- | --- | --- |
| Principal: linhas 18–71 | 0,0283 | 0,0044 |
| Linhas 30–49 | 0,0025 | 0,0039 |
| Remoção apenas da média | 0,0256 | 0,0034 |
| Centros x −16 / y +16 px | 0,0285 | 0,0063 |
| Centros x +16 / y −16 px | 0,0281 | 0,0043 |
| Detector σ=12 px | 0,0280 | 0,0045 |
| Detector σ=20 px | 0,0287 | 0,0044 |
| CWT da v5; mesmos recortes | 0,0285 | 0,0052 |

As variantes de detector e de centros usam os controles fixos principais; faixa, remoção da média e v5 usam referências correspondentes recalculadas. A proteção de denominadores não descartou escalas no lote experimental. Nos pares ancorados, o máximo integral da Morlet foi 0,998678 e o do legado v5 foi 1,006157, sem clipping; o maior erro da equação foi 0.00e+00. As três pequenas ultrapassagens unitárias do legado (dois mapas ancorados e um fixo) ocorreram fora dos lags do resumo e foram reproduzidas por cálculo independente na [auditoria do v5](results/extensao_churn_anular/auditoria_ultrapassagens_v5.json).

## Tabelas auditáveis

- [Resultados de todos os pares e variantes](results/extensao_churn_anular/resultados_ancoras.csv), incluindo ausências explícitas.
- [Controles fixos calculados e reutilizados](results/extensao_churn_anular/resultados_fixos.csv).
- [Agregados por par de aquisições](results/extensao_churn_anular/agregados_ancoras.csv).
- [Contrastes por par](results/extensao_churn_anular/contrastes_pareados.csv) e [por padrão](results/extensao_churn_anular/contrastes_por_padrao.csv).
- [Sensibilidades no conjunto comum de pares](results/extensao_churn_anular/sensibilidade_conjunto_comum.csv), com principal recalculada no mesmo subconjunto, variante e mudança.
- [Painel das oito variantes na interseção dos mesmos 815 pares locais](results/extensao_churn_anular/sensibilidade_conjunto_global.csv), incluindo a linha principal no mesmo subconjunto.
- [Retirada de uma aquisição](results/extensao_churn_anular/retirada_aquisicao.csv), removendo todos os pares incidentes.
- [Média balanceada entre os dois padrões](results/extensao_churn_anular/contrastes_balanceados.csv); não substitui as conclusões por padrão.
- [Estabilidade óptica](results/extensao_churn_anular/estabilidade_ancoras.csv), [controles sintéticos](results/extensao_churn_anular/controles_sinteticos.csv) e [self/mesma aquisição](results/extensao_churn_anular/controles_mesma_aquisicao.csv).
- [Perfis por escala dos 830 pares principais](results/extensao_churn_anular/perfis_ancoras.npz) e [mapas completos dos exemplos](maps/extensao_churn_anular).

## Limitações

- **Âncora óptica, identidade física não validada.** O mínimo marca uma depressão do perfil de intensidade; não mede espessura do filme nem identifica uma crista, bolha ou a mesma onda em aquisições diferentes. Os rótulos churn/anular são globais, herdados do inventário.
- **Seleção exploratória.** Regra desenvolvida em quatro aquisições piloto e aplicada no mesmo conjunto, já conhecido; revisão visual assistida por IA, sem especialista independente, concordância entre anotadores ou conjunto externo.
- **Montagem e contexto.** 10 das 152 âncoras caem nas zonas de mistura de 2 px das costuras; as outras 142 tiveram a coluna de origem conferida pixel a pixel. Cada contexto cruza 38–162 costuras. 143 recortes sobrepõem outro recorte selecionado e 120 contêm outro marco selecionado. Considerando só os 76 mínimos principais: 70 sobrepõem outro segmento, 60 contêm outro marco e 3 estão em mistura. Isso impede atribuir toda a WCO ao ponto central. As contagens de proximidade de costuras não são prova de artefato, pois as faixas de montagem são densas.
- **Aquisição com ressalva.** m03q40_001 conserva 60 deslocamentos não positivos herdados; 2 dos 152 contextos selecionados incluem deslocamentos não positivos. A análise principal mantém essa aquisição, e sua retirada é mostrada separadamente; os panoramas de entrada não foram substituídos.
- **Comparabilidade dos controles.** Mínimos compatíveis têm razão de larguras ≤1,50; distintos, ≥2,00. As larguras medem meia proeminência óptica. Mínimo × máximo pode diferir também em largura e forma, e a WCO é invariável à inversão global de sinal. Janelas fixas podem conter estruturas e não são negativos físicos garantidos.
- **Cobertura e dependência.** Os controles de larguras distintas existem em menos pares de aquisições, especialmente no anular. Cada contraste usa somente sua interseção de aquisições. Pares compartilham registros; não são réplicas independentes. Não foram calculados valores-p nem intervalos de confiança artificiais.
- **Domínio e potência.** N=2.048, K=1, correlação enviesada 1/N, λ de 50 a 1.024 px e lags principais ±32 px. A WPS alternativa usa módulo; WCO não recebe clipping. Validade numérica do denominador não equivale a potência confiável. Não há conversão para mm, Hz ou velocidade, nem validação de classificador.
- **Núcleo legado v5.** Dois mapas ancorados e um fixo apresentaram pequenas ultrapassagens de WCO=1, chegando a 1,006157. Todas ocorreram fora de ±32 px. A soma direta das correlações e a CWT por convolução temporal reproduziram os valores, com discrepância máxima de 3,04×10⁻¹³. Foram preservados sem clipping como limitação da variante discreta; a Morlet principal permaneceu ≤1. [Auditoria específica](results/extensao_churn_anular/auditoria_ultrapassagens_v5.json).

<!-- AUDITORIA_ADVERSARIAL_INICIO -->
## Revisão adversarial dos resultados

A revisão confirmou a integridade das entradas e a implementação das fórmulas principais nos testes independentes realizados. Encontrou, porém, fragilidades de seleção, controles e interpretação, além de um defeito herdado na costura e imprecisões de documentação. Os resultados principais foram preservados; a conclusão permanece exploratória, sem melhora robusta nem identificação física demonstradas.

| Padrão | Correspondência local σ12 / σ20 | Seleção completa σ12 / σ20 | Seleção completa: linhas 30–49 |
| --- | --- | --- | --- |
| Agitado (churn) | 40/40 em ambos | 37/40; 38/40 | 27/40 |
| Anular ondulado | 36/36 em ambos | 29/36; 30/36 | 28/36 |

| Padrão | Controle | Mínimo / máximo em 100 sementes | Sementes com diferença positiva |
| --- | --- | --- | --- |
| Pistonado (plug) | Fixos | 0,0137 / 0,0395 | 100/100 |
| Pistonado (plug) | Larguras/extensões distintas | 0,0018 / 0,0458 | 100/100 |
| Golfadas (slug) | Fixos | 0,2192 / 0,2192 | 100/100 |
| Golfadas (slug) | Larguras/extensões distintas | 0,1520 / 0,1520 | 100/100 |
| Agitado (churn) | Fixos | 0,0120 / 0,0477 | 100/100 |
| Agitado (churn) | Larguras/extensões distintas | -0,0113 / 0,0432 | 86/100 |
| Anular ondulado | Fixos | -0,0137 / 0,0258 | 80/100 |
| Anular ondulado | Larguras/extensões distintas | -0,0711 / 0,0156 | 4/100 |

| Padrão | Âncoras − fixos original | Sem feather | Mudança |
| --- | --- | --- | --- |
| Pistonado (plug) | 0,0289 | 0,0286 | -0,000241 |
| Golfadas (slug) | 0,2192 | 0,2187 | -0,000478 |
| Agitado (churn) | 0,0283 | 0,0280 | -0,000312 |
| Anular ondulado | 0,0044 | 0,0046 | 0,000141 |

Controle com seleção completa: 0,0169 de ganho médio em AR(1) independentes, sem valor-p. Consulte [a auditoria completa](AUDITORIA_ADVERSARIAL.md), [sementes](results/auditoria_adversarial/sensibilidade_semente.csv), [pares admissíveis](results/auditoria_adversarial/todos_pares_admissiveis.csv) e [o controle sintético](results/auditoria_adversarial/nulo_selecao_completa.csv).
<!-- AUDITORIA_ADVERSARIAL_FIM -->
