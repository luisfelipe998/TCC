# Contraprova cruzada da auditoria de estatistica8

Revisão adversarial realizada pelo auditor da v14.1, após concluir sua avaliação independente. O objeto foi [auditoria_estatistica8.md](auditoria_estatistica8.md), sobretudo os achados E8-02 e E8-03. O relatório do primeiro auditor não foi alterado. Código original, texto do TCC e evidências numéricas foram confrontados; os resultados desta contraprova estão em [resultados/contraprova_estatistica8](../resultados/contraprova_estatistica8/).

**Os achados principais resistiram à contraprova, com uma delimitação essencial: a seleção óptica não distingue churn de anular no algoritmo atual.** A crítica de extração dependente de família não pode ser usada para descartar automaticamente o principal contraste largura/churn–anular. A redução do efeito nas janelas distribuídas foi confirmada independentemente e sustenta a restrição de representatividade, sem evidenciar escolha intencional de uma âncora favorável.

## E8-02 Seleção dependente da família

O ramo de código é real: `v16.3_ancoras/scripts/prepare_comparisons.py:54–65` chama `event_candidates` para plug/slug e `texture_candidates` para os demais. Essa é a evidência direta do lote completo, complementar à citação do piloto em `select_anchors.py:262–265`. O código depois preserva a âncora nativa no contexto comparável de 2048 px (`prepare_comparisons.py:20–35,74–80`).

Para churn e anular, o ranking é o mesmo — maior proeminência, estabilidade e coordenada (`select_anchors.py:145–190`) — e a janela nativa preferida é de 4096 px para ambos (`scripts/config.json:8`). A afirmação do primeiro auditor de que “a objeção não invalida automaticamente o contraste principal” foi submetida a uma verificação adicional:

1. Foram relidos os 19 panoramas de churn/anular e extraídas as mesmas 54 linhas.
2. Cada sinal foi processado duas vezes pelas funções de seleção existentes, uma com rótulo churn e outra com rótulo wavy annular.
3. Foram comparados centro escolhido, janela nativa, identificador da âncora e número de candidatos.

**Nos 19 casos, esses resultados foram idênticos sob a troca de rótulo; todos os centros também coincidiram com as âncoras arquivadas.** A execução ocorreu em memória, sem chamar `setup`, geração de arquivos do pipeline ou cálculo wavelet.

Isso refuta uma formulação excessiva de E8-02 — “o contraste churn/anular é fabricado porque cada um recebe um detector diferente”. Os dois grupos recebem o mesmo detector. Não refuta o risco nos contrastes entre famílias e não resolve confusão operacional ou desenvolvimento informado pelos mesmos dados. Uma operação determinística comum preservaria permutabilidade se os sinais de entrada fossem permutáveis; o desenho observacional não fornece essa premissa automaticamente.

**Juízo:** manter E8-02 com seu contra-argumento atual. Não converter “a família é consultada” em acusação de seleção de âncoras pelo maior escore WPS. O código de seleção não consulta as wavelets. Usar também a referência do lote completo acima torna a procedência mais direta.

## E8-03 Estabilidade local e representatividade espacial

Foi confirmada a sobreposição `(2048−64)/2048 = 0,96875` entre cada janela deslocada e a original. Esse percentual se refere aos pixels brutos compartilhados; após remover tendência e normalizar cada janela, não significa que 96,875% dos valores preparados sejam numericamente idênticos. A proximidade de conteúdo continua suficiente para distinguir o teste local de uma exploração do panorama inteiro.

Os números das cinco posições uniformes e da mediana por aquisição foram recalculados com uma segunda implementação de silhueta e delta de Cliff, a partir das larguras por janela arquivadas na v17. O erro máximo foi **1,25 × 10⁻¹⁶**. A inversão dos quantis contínuos foi conferida por busca de raiz delimitada, em vez da expressão algébrica do primeiro auditor; o erro máximo da largura foi **2,66 × 10⁻¹⁵** nos 25 perfis.

| Comparação | Resultado da contraprova |
|---|---|
| Âncora contínua, silhueta macro | 0,4303 |
| Cinco posições uniformes, silhuetas macro | 0,1363; 0,2089; 0,2001; 0,1445; 0,1782 |
| Mediana por aquisição, silhueta macro | 0,2648 |
| Delta de Cliff na âncora e na mediana por aquisição | 0,9556 e 0,7333 |
| Quantis contínuos | Cálculo confirmado por integração/inversão independente |

**Juízo:** a magnitude mais forte na âncora não representa uniformemente as posições amostradas. Contudo, os dados mantêm ordenação parcial nas alternativas; não demonstram ausência de separabilidade nem suprem a necessidade de novas aquisições. As cinco janelas não são réplicas independentes, nem o experimento foi uma amostragem probabilística exaustiva das posições. A versão contínua muda a representação da densidade, mas a comparação âncora contínua versus janelas contínuas controla esse fator.

Não foi refeita nesta contraprova a extração wavelet das 125 janelas: o primeiro auditor realizou essa etapa diretamente dos pixels. A verificação adicional aqui cobre a seleção de âncoras, o resumo estatístico das larguras e a integração dos quantis. Esta distinção evita contar duas leituras dos mesmos CSVs como duas reconstruções completas de sinal.

## Demais pontos examinados

- **E8-01 e E8-06:** diferenças de amostra, recorte, faixa espectral e domínio da potência média são reais. O TCC identifica sua execução e não importa os números de estatistica8. “Transferência indevida” deve permanecer um risco de interpretação conjunta, não uma inconsistência já demonstrada no texto.
- **E8-04:** permutar em m ou em q não conserva conjuntamente m×q. Com um registro por célula conjunta, não há permutação interna não trivial nesse desenho. Isso restringe uma conclusão isolada sobre padrão; não implica que qualquer ajuste sob hipóteses modeladas seja matematicamente impossível.
- **E8-05:** a distribuição da diferença de concordâncias sob permutação global testa ausência de associação dos rótulos sob aquela permutabilidade. Ela não é, em geral, a distribuição sob igualdade de desempenho de dois procedimentos com sinal. A ressalva do relatório está correta, e a tarefa atual não exige desenvolver uma validação preditiva alternativa.
- **E8-07:** `robust_pass` requer silhueta macro positiva após exclusão; não exige o mínimo de classe positivo nesse mesmo controle (`analyze.py:159–166`). Trata-se de diferença entre critério declarado e possível leitura da palavra “robusto”; o primeiro relatório não o exagera como erro oculto. A largura/churn–anular mantém ambos os mínimos positivos.

## Precaução na síntese das três auditorias

Os dois conjuntos de **25** não são iguais. Estatistica8 usa quatro classes, retém m03q40 e exclui a única bolha dispersa. As estatísticas da v14.1 excluem m03q40 pela saturação e retêm a bolha dispersa. O TCC/v16 usa 26. Uma comparação de métodos deve alinhar os identificadores, não apenas as contagens. Tampouco se deve confundir os pares derivados de cada conjunto com ensaios independentes.

Artefatos: [script da contraprova](../scripts/counterreview_estatistica8.py), [troca de rótulos nas âncoras](../resultados/contraprova_estatistica8/label_swap_anchor_checks.csv), [rechecagem dos resumos](../resultados/contraprova_estatistica8/window_statistics_recheck.csv) e [resumo com hash do relatório revisado](../resultados/contraprova_estatistica8/summary.json).

Reprodução, a partir da raiz:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -B src2/v17/scripts/counterreview_estatistica8.py
```
