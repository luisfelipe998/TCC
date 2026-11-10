# v18 — metodologia completa nas linhas 18 a 71

Esta versão reuniu montagem, preparo dos sinais, CWT, Pearson, segmentação óptica, fração de vazio, transições, resumos e figuras dos **26 ensaios dos cinco padrões puros** descritos na metodologia. Os cálculos e geradores são locais. A única entrada experimental externa é `data_processed`; versões anteriores, auditorias e arquivos do TCC não são necessários para executar a análise.

A consolidação inicial de 30/10/2026 foi aprovada na execução integral e na repetição em cópia isolada, com **186 verificações sem divergências** e 23 processos auditados no teste de portabilidade. As saídas ativas foram substituídas somente após essas verificações, e o estado anterior foi arquivado. Evidências dessa consolidação: [comparação completa](validacao/comparacao_completo.json), [portabilidade](validacao/portabilidade.json) e [teste de rejeição de alterações](validacao/teste_protecao.json).

A atualização posterior acrescentou os resumos de Pearson e coerência em `06_resumos/`, correspondentes às seções 3.4.2, 3.5 e 4.4–4.7 do TCC. Os registros anteriores foram preservados. A [validação da atualização](validacao/resumos_pearson_coerencia.json) reúne os cálculos, a verificação independente, a reprodução isolada e a versão do PDF. A [correspondência com a metodologia](validacao/sincronizacao_metodologia.json) foi renovada para o texto atual.

A consolidação editorial seguinte manteve os cálculos e todos os arquivos de `06_resumos/`. No TCC, a tabela por combinação passou de 6 a 5, e a síntese global passou de 7 a 6. Os quatro resumos de cada exemplo saíram da antiga Tabela 5 e passaram às legendas das Figuras 16–19 e 22. Os dados completos continuam em `06_resumos/dados/pares_ilustrados.csv`, sem recálculo ou mudança na metodologia.

A revisão posterior de consistência da Discussão conferiu as interpretações contra os resultados e a metodologia, preservando os dados, o código e as figuras. A [validação editorial](validacao/consistencia_discussao.json) identifica aquele PDF e os três pareceres independentes.

A revisão seguinte acrescentou texto de abertura à seção 3.5 e explicitou em 3.6 o critério óptico de saturação já implementado. A [conferência estrutural](validacao/secoes_com_texto.json) identifica o PDF daquela etapa. Dados, código e resultados permaneceram iguais.

A revisão didática da Discussão simplificou as seções 5.1–5.4 e conciliou a interpretação da coerência com os Resultados e a Conclusão. As estatísticas dos mapas foram preservadas. O contraste entre máximo e mediana deixou de ser usado como argumento sobre a semelhança global dos sinais. A [validação desta revisão](validacao/discussao_didatica.json) identifica o PDF daquela etapa e os pareceres independentes. O processamento, os dados, as figuras e a metodologia permaneceram iguais.

O esclarecimento posterior da seção 5.1 explicitou a participação do único ensaio de bolhas dispersas em 25 pares interpadrões. A [validação do ajuste](validacao/esclarecimento_bolhas.json) identifica o PDF daquela etapa, com todos os cálculos preservados.

A revisão de clareza das seções 5.1–5.4 priorizou a interpretação dos achados e concentrou suas limitações. A [validação editorial](validacao/discussao_clareza.json) identifica o PDF daquela etapa. Os dados, cálculos, figuras e a metodologia foram preservados.

A síntese posterior da Conclusão conectou os achados à questão da Introdução e aos objetivos, com sugestões breves de trabalhos futuros. A [validação da Conclusão](validacao/conclusao_sintese.json) identifica o PDF daquela etapa. A análise e os capítulos anteriores permaneceram inalterados.

A revisão adversarial integral conferiu novamente o texto, as equações, a implementação e os artefatos, com três pareceres independentes. Os ajustes ficaram no resumo e no abstract e em duas pontuações da metodologia. A [validação integral](validacao/revisao_integral_final.json) identifica o PDF daquela etapa e a auditoria. Dados, código, equações, figuras, resultados, discussão e conclusão permaneceram iguais.

O esclarecimento posterior da Conclusão explicitou a média de cada mapa de coerência wavelet e a média entre pares, com igual peso, de 0,0874 entre padrões e 0,0607 dentro do mesmo padrão. A [validação do esclarecimento](validacao/clareza_medias_conclusao.json) identifica o PDF daquela etapa. Os valores foram reconferidos nos 325 pares, sem alteração dos dados ou cálculos.

A revisão final da Conclusão reforçou a ligação aos objetivos, a interpretação conjunta das medidas espectrais e ópticas e os limites das comparações. A [validação da revisão](validacao/conclusao_objetivos_clareza.json) identifica o PDF atual. Somente a redação da Conclusão foi alterada, sem novos cálculos ou mudanças nos resultados científicos.

## Executar tudo

A revisão de formatação do TCC até o capítulo 4 preservou integralmente os dados e o código desta análise. A correspondência daquele PDF, dos fontes e das citações está registrada na [auditoria de formatação](../../auditorias/formatacao-ate-cap4/README.md), com o estado anterior arquivado. As alterações da metodologia limitaram-se à diagramação.

Na raiz do projeto, com as bibliotecas de `requirements.txt` instaladas:

```sh
python3 -B src2/v18_montagem_18_71/executar.py --workers 2
```

O comando cria uma pasta nova em `execucoes/`, executa todas as etapas e compara os resultados com a referência congelada em `config/`. Ele interrompe a sequência se houver diferenças científicas ou falha de execução. Nunca substitui automaticamente os resultados ativos. O resultado final e os logs ficam na pasta informada ao término; `manifesto_execucao.json` precisa indicar `status: PASS` e `complete: true`.

Para mover a v18 ou escolher outro destino, use caminhos explícitos:

```sh
python3 -B /caminho/v18_montagem_18_71/executar.py \
  --data-root /caminho/data_processed \
  --output /caminho/saida_nova --workers 2
```

A saída deve estar vazia. Os ensaios precisam conservar as pastas e os nomes declarados em `config/ensaios.json`, com os 501 JPEGs, `_config.xsv` e o rótulo `__<padrão>.txt`. As referências numéricas servem somente à verificação; todos os sinais, máscaras e medidas são recalculados dos JPEGs.

## Organização

| Caminho | Conteúdo |
|---|---|
| `config/` | Seleção dos 26 ensaios e referências de validação |
| `scripts/` | Núcleos numéricos, montagem, óptica, geradores e verificadores locais |
| `01_panoramas/` | Duas variantes dos 26 panoramas; deslocamentos e hashes dos dados |
| `02_wavelet/` | Sinais, 26 perfis individuais e métricas completas dos 325 pares |
| `03_optica/` | Referência de líquido, máscaras, medidas, elegibilidade, envoltórias e figuras |
| `04_analise_simplificada/` | Código de Pearson direto; matriz, tabelas, cinco mapas e oito figuras |
| `05_figuras/` | Quatro exemplos de panoramas e quatro figuras de metodologia |
| `06_resumos/` | Estatísticas completas da WCO por par, resumos de Pearson e coerência média, figura global |
| `validacao/` | Comparações numéricas e relatório de portabilidade da execução consolidada |
| `historico/` | Estados anteriores preservados, sem participação no processamento |
| `execucoes/` | Execuções novas, com logs e manifestos próprios |

[METODOLOGIA.md](METODOLOGIA.md) relaciona as seções do capítulo 3 aos scripts e descreve os parâmetros usados. A faixa `18:72` corresponde às **54 linhas de índices 18 a 71, inclusive**, desde o alinhamento até as medidas ópticas.

Pearson compara diretamente a potência no atraso zero dos **sinais completos**, com 200 valores por ensaio. As CWTs foram aplicadas às autocorrelações completas antes dessa extração. Os mapas dos pares usam os **trechos iniciais de comprimento comum**, cortados antes do preparo. Não há comparação de resíduos logarítmicos no processamento vigente.

A etapa complementar [06_resumos](06_resumos/README.md) obtém mínimo, média, mediana e máximo do campo completo de coerência de cada par, preservando Pearson e as etapas anteriores. O comando principal a executa após a comparação com a referência congelada. Use `--sem-resumos` somente para repetir o percurso anterior. A etapa também pode ser executada isoladamente sobre as saídas existentes, sem reconstruir os panoramas ou a análise óptica.

## Verificar independência das outras pastas

```sh
python3 -B src2/v18_montagem_18_71/scripts/verificar_portabilidade.py \
  --data-root src2/data_processed --workers 2
```

Esse comando copia somente código e configuração para uma pasta temporária, inicia uma saída vazia e repete toda a análise. Um monitor de acesso a arquivos bloqueia dependências de outras pastas do projeto e o carregamento de código de `data_processed`. Bibliotecas Python e recursos do sistema, como fontes, continuam disponíveis. O relatório é gravado em `validacao/portabilidade.json`; trata-se de verificação de dependências da execução Python, não de isolamento de segurança do sistema operacional.

## Critérios e limites

As comparações exigem igualdade dos pixels dos panoramas e máscaras, dos deslocamentos, da seleção dos pares e da elegibilidade óptica. Os números são comparados com tolerância de ponto flutuante registrada no relatório e também na precisão usada no texto. As referências anteriores ficam preservadas em `historico/pre_autocontida_20261030/`; o verificador usa apenas a cópia portátil em `config/`.

Os resumos ópticos incluem 25 ensaios. `m03q40_001` ultrapassa 20% de pixels saturados na faixa de 54 linhas e permanece apenas na análise espectral. A fração de vazio continua condicionada à máscara e à representação de cada intervalo de gás como seção circular, com diâmetro do tubo de 54 pixels. A avaliação é descritiva, sem tratar os pares ou escalas como réplicas independentes.

O núcleo preserva opções alternativas de compatibilidade, mas a execução define explicitamente Morlet própria, estimador com divisor N e preparo com retirada de tendência linear. PyWavelets é uma dependência de importação desse núcleo; a CWT executada usa o filtro próprio e a FFT da SciPy.
