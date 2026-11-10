"""Relatórios locais em PT-BR gerados dos resultados conferidos."""
from pathlib import Path
import json
import textwrap
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/wco-eventos-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from analysis_common import ROOT,BASE
from inspect_inputs import LABELS

def f(v,n=4):return f'{float(v):.{n}f}'.replace('.',',')
def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rows])
def write(name,text):(ROOT/name).write_text(text.strip()+'\n')

def main():
    out=ROOT/'results';c=json.loads((out/'conclusao.json').read_text());timing=json.loads((out/'execucao.json').read_text())
    cover=pd.read_csv(out/'cobertura.csv');by=pd.read_csv(out/'contrastes_por_padrao.csv');primary=by[by.variante=='principal']
    contrasts=pd.read_csv(out/'contrastes_pareados.csv');ev=pd.read_csv(out/'resultados_eventos.csv');fi=pd.read_csv(out/'resultados_fixos.csv')
    cat=pd.read_csv(out/'catalogo_eventos.csv');accepted=cat[cat.incluida]
    tests=json.loads((out/'validacao_numerica.json').read_text());synthetic=pd.read_csv(out/'controles_sinteticos.csv').set_index('controle')
    examples=json.loads((out/'exemplos.json').read_text());freeze=json.loads((out/'congelamento.json').read_text())
    patterns=['dispersed bubbles','plug','slug','churn','wavy annular']
    inventory=pd.read_csv(BASE/'results/input_audit.csv').query('included').sort_values('run')
    coverage=cover.set_index('padrao')
    pattern_names={**LABELS,'plug':'Pistonado / bolhas alongadas (plug)','wavy annular':'Anular ondulado (wavy annular)',
                   'dispersed bubbles':'Bolhas dispersas (dispersed bubbles)'}
    roles={p:(f'{int(coverage.loc[p,"eventos_aceitos"])} eventos aceitos; comparação ancorada e referências fixas'
              if p in ['plug','slug'] else 'Referências fixas e atlas; nenhum evento ancorado aceito') for p in patterns}
    pattern_table=table(['Rótulo nos arquivos','Nome em português','Aquisições','Participação na análise'],[
        [p,pattern_names[p],int(coverage.loc[p,'aquisicoes']),roles[p]] for p in patterns])
    main_events=ev[ev.variante=='principal'].set_index('pair_id')
    cross_id=ev[(ev.variante=='principal')&(ev.grupo=='entre_equivalentes')].sort_values('pair_id').pair_id.iloc[0]
    case_ids={'favoravel':examples['favoravel'],'desfavoravel':examples['desfavoravel'],'entre_padroes':cross_id}
    case_names={'favoravel':'Favorável','desfavoravel':'Desfavorável','entre_padroes':'Controle entre padrões'}
    case_rows=[];case_labels={};identification=[]
    for kind,pid in case_ids.items():
        row=main_events.loc[pid]
        case_labels[kind]=f'{row.pattern_x} × {row.pattern_y}'
        case_rows.append([case_names[kind],f'{LABELS[row.pattern_x]} — {row.event_x}',f'{LABELS[row.pattern_y]} — {row.event_y}'])
        identification.append(dict(caso=kind,pair_id=pid,padrao_x=row.pattern_x,padrao_y=row.pattern_y,event_x=row.event_x,event_y=row.event_y))
    (out/'identificacao_exemplos.json').write_text(json.dumps(identification,ensure_ascii=False,indent=2)+'\n')
    cases_table=table(['Exemplo','Sinal x: padrão e evento','Sinal y: padrão e evento'],case_rows)
    def get(pat,comp='distintos',var='principal'):
        return by[(by.variante==var)&(by.escopo=='intra')&(by.comparador==comp)&(by.padrao==pat)].iloc[0]
    plug=get('plug');slug=get('slug');plugfix=get('plug','fixos');slugfix=get('slug','fixos')
    eqtable=table(['Padrão','Pares de aquisições','Equivalentes','Distintos','Diferença','Diferença frente a recortes fixos'],[
        [LABELS[pat],int(get(pat).pares_aquisicoes),f(get(pat).equivalentes),f(get(pat).controle),f(get(pat).diferenca),f(get(pat,'fixos').diferenca)] for pat in ['plug','slug']])
    covtable=table(['Padrão','Aquisições elegíveis','Eventos aceitos','Aquisições com equivalentes intrapadrão','Pares de eventos equivalentes'],[
        [LABELS[r.padrao],r.aquisicoes,r.eventos_aceitos,r.aquisicoes_em_equivalentes_intra,r.pares_eventos_equivalentes_intra] for r in cover.itertuples()])
    sensrows=[]
    for pat in ['plug','slug']:
        g=by[(by.padrao==pat)&(by.escopo=='intra')&(by.comparador=='distintos')&(by.variante!='v5_mesmos_segmentos')]
        sensrows.append([LABELS[pat],f(g.diferenca.min()),f(g.diferenca.max()),int(get(pat).pares_aquisicoes)])
    crosse=primary[(primary.escopo=='entre')&(primary.comparador=='distintos')].iloc[0]
    lead=(f'**Não conseguimos melhorar de forma convincente a alternative WCO com este protocolo em todo o escopo.** '
          f'Para estruturas alongadas em pistonado (plug), o ganho sobre eventos de extensão distinta foi **{f(plug.diferenca)}**, uma **melhora marginal**. '
          f'Em golfadas (slug), a diferença foi **{f(slug.diferenca)}**, mas dependeu de **um único par de eventos de duas aquisições** e se inverteu para **{f(get("slug",var="faixa30_49").diferenca)}** na faixa 30–49. '
          'O estudo fornece evidência exploratória restrita; não demonstra melhora geral, classificação dos padrões ou identidade física entre estruturas.')
    executive=f'''# Resumo executivo da comparação de eventos por wavelets

{lead}

A execução usou a base solicitada **v16_spatial2_18_71**, com média das linhas **18–71 inclusive**, e preservou as equações alternativas de Schuck. O objetivo foi verificar se parear estruturas pela imagem, antes de observar a coerência, acrescentaria similaridade específica aos segmentos. A unidade calculada foi um segmento de **2.048 px centrado em um evento**, com contexto; ele pode conter estruturas vizinhas.

A [explicação de âncoras, limites dos sinais e autocorrelação](ANCORAS_E_AUTOCORRELACAO.md) detalha o processamento com um exemplo real e distingue o método executado da possível extensão a churn e anular.

## Identificação dos padrões

{pattern_table}

**Plug e slug são rótulos diferentes nesta base:** plug corresponde ao pistonado/bolhas alongadas; slug corresponde a golfadas. Os nomes acima são os rótulos globais das aquisições, herdados do [inventário auditado](results/auditoria_entradas.csv). A categoria óptica do evento central não redefine o padrão da aquisição.

**Os resultados de eventos ancorados abaixo se referem somente a plug e slug.** Churn, anular ondulado e bolhas dispersas aparecem na cobertura, nas imagens e nas referências fixas. Nos gráficos de eventos, “entre padrões” significa especificamente **plug × slug**. Os controles fixos gerais incluem os cinco padrões.

Nos gráficos de contraste e sensibilidade, azul-petróleo identifica plug e laranja identifica slug. No gráfico de cobertura, as cores identificam grupos de inclusão. Nos mapas de WPS/WCS/WCO, as cores representam o valor da grandeza da barra lateral; o padrão está escrito no título. Os controles sintéticos são sinais artificiais, sem rótulo de escoamento.

## Resultado principal e sua cobertura

{eqtable}

As médias acima atribuem peso igual a cada par de aquisições, após resumir seus pares de eventos. Os controles distintos e fixos pertencem aos **mesmos pares de aquisições** dos equivalentes. Em plug foram 18 pares de eventos equivalentes, distribuídos pelos três pares possíveis de aquisições. Em slug, somente um dos três pares possíveis teve equivalentes; não é possível generalizar seu ganho às outras aquisições. As diferenças são absolutas na escala da WCO, não percentuais de acurácia.

Em plug, **um dos três pares de aquisições piorou** frente aos recortes fixos. A média positiva não representa melhora em todos os pares.

![Contrastes pareados entre eventos equivalentes, distintos e recortes fixos](figures/contrastes.png)

*Cada ponto representa um par de aquisições, com dependência entre pontos que compartilham registros. A linha de 0,05 é uma margem descritiva definida antes do cálculo; não é significância estatística nem requisito físico validado. Não foram produzidos intervalos de confiança artificiais a partir de eventos sobrepostos.*

Foram conferidas **31 entradas**, mantidas as **26 aquisições elegíveis** e as cinco exclusões anteriores. O detector propôs **{len(cat)} intervalos alongados**; **26** foram aceitos em **6 aquisições** após revisão de morfologia e contexto. Apenas **5/26 aquisições ({f(500/26,1)}%)** participaram de pares equivalentes intrapadrão. Essa é cobertura do protocolo, não a fração de todos os eventos físicos detectados.

{covtable}

![Cobertura do protocolo em cada padrão](figures/cobertura.png)

*Zero eventos aceitos não significa ausência física de estruturas. A delimitação de eventos alongados não foi estabelecida para as texturas de churn/anular sob esta regra. Bolhas dispersas contém objetos pequenos e somente uma aquisição elegível, sem repetição intrapadrão independente.*

## Especificidade limitada do pareamento

O ganho em plug permaneceu pequeno frente a eventos distintos em extensão. O ganho de slug é um caso promissor, porém isolado. Entre plug e slug, eventos localmente equivalentes também superaram os distintos: diferença média **{f(crosse.diferenca)}** nos **{int(crosse.pares_aquisicoes)} pares de aquisições** com ambos os grupos. Portanto, o critério pode refletir alguma semelhança local sem separar os rótulos globais.

![Distribuição das médias por par de aquisições dentro e entre padrões](figures/entre_padroes.png)

*Os controles entre padrões também usam a mesma categoria alongada e a mesma regra de extensão. As distribuições são descritivas; conjuntos com coberturas diferentes não foram usados para atribuir um ganho de classificação.*

## Estabilidade e controles negativos

As mesmas coordenadas e os mesmos parceiros foram avaliados com a faixa 30–49, remoção apenas da média e perturbações de ±16 px nos centros. Os ganhos sobre eventos distintos, excluindo o núcleo v5, variaram assim:

{table(['Padrão','Menor diferença','Maior diferença','Pares de aquisições'],sensrows)}

**O ganho de slug não resistiu à mudança de faixa vertical:** com 30–49, equivalentes ficaram abaixo dos distintos em {f(abs(get('slug',var='faixa30_49').diferenca))}. As variações de centro preservaram o sinal do contraste médio, mas não compensam essa dependência do observável óptico.

![Sensibilidades sem escolher a variante de maior coerência](figures/sensibilidades.png)

O controle sintético com **fundo comum e eventos distintos** atingiu média **{f(synthetic.loc['fundo_comum_eventos_distintos','escore'])}**. Eventos suaves com larguras distintas atingiram **{f(synthetic.loc['eventos_suaves_distintos','escore'])}**, e seu máximo integral chegou a **{f(synthetic.loc['eventos_suaves_distintos','maximo_integral'])}**. Assim, WCO alta, especialmente um máximo isolado, não demonstra equivalência de eventos.

![Controles sintéticos e diferença entre média e máximo integral](figures/controles_sinteticos.png)

## Imagens e espectros dos cinco padrões

![Panoramas, sinais e WPS alternativa dos cinco padrões](figures/cinco_padroes.png)

*De cima para baixo: bolhas dispersas — m07q01_001; pistonado (plug) — m05q01_001; golfadas (slug) — m04q01_001; agitado (churn) — m03q05_001; anular ondulado (wavy annular) — m03q10_001. Plug/slug: primeiro evento aceito da primeira aquisição em ordem alfabética. Demais padrões: referência central fixa, explicitamente sem anotação de evento. Os mapas usam o módulo da CWT da autocorrelação, sem elevá-lo ao quadrado. O eixo de posição do sinal não é o eixo de lag dos mapas. A banda de 50–1.024 px não pretende descrever todas as pequenas bolhas visíveis.*

## Casos favorável e desfavorável

Os dois exemplos abaixo são os extremos do escore intrapadrão **selecionados somente para ilustração após o lote**. Eles não definiram o catálogo nem os parceiros, e não representam desempenho típico.

{cases_table}

As letras **x** e **y** identificam os dois sinais do par. Nas curvas, x usa linha contínua e y usa linha tracejada; ambas podem pertencer ao mesmo padrão. Os títulos de cada panorama e de cada WPS identificam o respectivo padrão e evento.

**Caso favorável: {case_labels['favoravel']} — ambos os sinais são pistonado.**

![Exemplo favorável: pistonado (plug) × pistonado (plug), com WPS, WCS e WCO alternativas](figures/pares/favorável.png)

*Par: {examples['favoravel']}. A WCO usa numerador ao quadrado e denominador no lag zero. O máximo integral não é o resumo principal.*

**Caso desfavorável: {case_labels['desfavoravel']} — ambos os sinais também são pistonado.**

![Exemplo desfavorável: pistonado (plug) × pistonado (plug)](figures/pares/desfavorável.png)

*Par: {examples['desfavoravel']}. Compartilhar a categoria visual e extensão semelhante não garantiu alta coerência em toda a banda. As linhas de suporte não são contornos de significância.*

**Controle entre padrões: {case_labels['entre_padroes']} — x é pistonado (plug); y é golfadas (slug).**

![Controle entre padrões: sinal x pistonado (plug) × sinal y golfadas (slug)](figures/pares/entre-padrões.png)

*Controle entre padrões escolhido pelo primeiro ID, sem selecionar sua WCO. A aparência local não certifica o rótulo global.*

## Limitações que condicionam a conclusão

- **Poucas aquisições e seleção exploratória:** só três registros por padrão plug/slug; o conjunto já havia sido analisado. Não houve teste externo, validação por especialista independente ou medida entre anotadores.
- **Equivalência operacional limitada:** a seleção combinou categoria óptica revisada e extensão com razão ≤1,25. Ela não valida igualdade física, forma completa ou identidade de bolhas. **{c['eventos_sem_parceiro']} eventos aceitos ficaram sem parceiro equivalente entre aquisições.**
- **Contexto e dependência:** {c['eventos_com_vizinhos_aceitos']} dos 26 segmentos contêm ao menos outro evento aceito; estruturas rejeitadas ou menores também podem estar no contexto. Os segmentos se sobrepõem. Não é possível atribuir toda a WCO exclusivamente ao evento central.
- **Fronteiras e montagem:** variar o limiar de intensidade de 60 para 55/65 deslocou uma fronteira em até **{f(c['max_mudanca_fronteira'],0)} px**. A perturbação de centros não cobre toda a incerteza morfológica. Os panoramas mantêm costuras e deslocamentos originais, inclusive os 60 não positivos de m03q40_001 e o caso de m04q01_001; não houve reconstrução física nova.
- **Escopo restrito:** churn e anular foram mantidos nos controles fixos e no atlas, com delimitação de eventos não estabelecida. Bolhas dispersas não oferece repetição entre aquisições. A inspeção de quadros foi amostral; não se afirma ausência de eventos em todos os quadros.
- **Grandezas distintas:** K=1, 1 px/amostra, sem mm/pixel, Hz ou velocidade inferidos. O divisor N é adaptação explícita de Schuck. O resumo de eventos é média em lags −32…+32 e ondas de 50–1.024 px; os escores antigos de máximos ponderados em 50–3.000 px não são diretamente comparáveis.
- **Peso das escalas:** pesos uniformes também incluem escalas de pouca potência, desde que numericamente válidas. A proteção do denominador evita divisões degeneradas; ela não estabelece confiabilidade experimental ou significância dessas escalas.

## Verificações e entrega

Passaram os controles de correlação por soma direta, convolução temporal independente da Morlet, padding, denominadores, recorte antes do preparo, self, atraso conhecido e núcleo literal do v5. O atraso de +100 px em y foi recuperado em **−100 px**. As **2.430 comparações** do lote principal/sensibilidades e mais **31 controles self/mesma aquisição** foram concluídos. O maior erro de recálculo da Eq. 32.17 no lote foi **{c['erro_equacao_maximo']:.2e}**; nenhuma escala experimental foi descartada por denominador numérico nas versões Morlet. Esses testes verificam a implementação, não a validade física das fronteiras.

O lote de cálculo levou **{f(timing['segundos_total']/60,1)} minutos**, com pico de memória reportado de **{f(timing['pico_memoria_mib'],0)} MiB**. A estimativa prévia foi de 7–13 minutos. A revisão visual e a preparação do protocolo ocorreram antes desse lote; geração de figuras e verificação final são etapas separadas.

Veja o [passo a passo com limitações](PASSO_A_PASSO.md), o [protocolo congelado](PROTOCOL.md), as [decisões](DECISIONS.md), os [resultados detalhados](RESULTS.md), o [catálogo de eventos](results/catalogo_eventos.csv), as [instruções de reprodução](README.md) e a [verificação final](results/verificacao_final.json). Uma versão compacta para compartilhar está no [PDF do resumo executivo](RESUMO_EXECUTIVO.pdf).
'''
    write('EXECUTIVE_SUMMARY.md',executive)
    steps=f'''# Passo a passo da análise executada

A análise foi executada a partir do plano histórico v16.2, com a atualização do usuário para a base **v16_spatial2_18_71**. O objetivo foi comparar segmentos de estruturas visualmente equivalentes sem usar coerência na seleção.

Para acompanhar as fórmulas e os limites com um exemplo real, veja a [explicação de âncoras e autocorrelação](ANCORAS_E_AUTOCORRELACAO.md).

1. **Conferência das referências e do escopo.** Foram lidos o plano, os relatórios e o núcleo da base; as Eqs. 32.4, 32.6, 32.8 e 32.14–23 foram conferidas nas imagens das páginas impressas 411, 413 e 414 de Schuck. Foram mantidas as 31 entradas auditadas, 26 elegíveis e cinco exclusões. **Limitação:** os números do plano precediam revisões da base; os atuais foram mantidos separados dos novos resultados.

2. **Auditoria de entrada sem reconstruir os panoramas.** Foram verificados os hashes dos 31 PNGs e das 31 planilhas. A média das linhas 18–71 foi extraída dos PNGs e comparada exatamente ao cache da base; as planilhas continuaram sendo referências de 20–59. Deslocamentos e costuras foram conferidos nos caches. **Limitação:** reprodução computacional herdada não garante reconstrução física correta. As fontes e seus hashes estão em [fontes_sha256.json](results/fontes_sha256.json).

3. **Inspeção anterior à WCO.** Foram produzidos e examinados 26 atlas com panorama completo, recortes fixos de início/centro/fim e quadros 0, 250 e 500. O piloto inicial cobriu os cinco padrões; as outras quatro aquisições plug/slug foram vistas antes de concluir as regras. **Limitação:** desenvolvimento exploratório, sem conjunto independente; os 78 quadros fixos não cobrem todos os 13.026 quadros elegíveis. [Observabilidade por aquisição](results/observabilidade.csv).

4. **Definição operacional de evento e rejeições.** Usou-se a fração de pixels <60 nas linhas 18–71, mediana de 5 colunas, limiar de fração 0,25 e fechamento de 13 colunas. Os candidatos tinham ≥96 px. O preenchimento da operação morfológica repetiu a borda, evitando criar falsos limites. A revisão visual rejeitou junções/agrupamentos e contornos ambíguos; o contexto descartou intervalos sem margem ou >1.536 px. Das {len(cat)} propostas, 26 foram aceitas. **Limitação:** os candidatos não formam inventário completo de todos os eventos físicos; só os 67 candidatos das seis aquisições com morfologia alongada passaram por revisão individual de suas imagens. Nas outras aquisições, não se confirmou a categoria nos atlas. [Revisão visual](results/revisao_visual.csv).

5. **Rastreabilidade e estabilidade da anotação.** Foram registrados limites [início, fim), âncora, linhas, contexto, costuras, qualidade, inclusão/motivo, quadro/faixa que forneceu a âncora e colunas correspondentes. Limiar 55/65 avaliou a estabilidade sem escolher o melhor resultado wavelet. A maior mudança de fronteira foi {f(c['max_mudanca_fronteira'],0)} px. **Limitação:** uma segunda parametrização não é uma segunda leitura independente; o quadro de origem da faixa não prova identidade física da estrutura.

6. **Congelamento antes da avaliação.** Regras, revisão, coordenadas e parceiros foram registrados por SHA-256 em {freeze['data_utc']}. Nenhuma WCO experimental nova havia sido calculada. O pareamento exigiu categoria alongada e razão de extensões ≤1,25; razão ≥1,50 definiu controles distintos. A faixa intermediária ficou fora da avaliação. Até 12 pares por combinação de aquisições/grupo foram selecionados por hash com semente 20261024. **Limitação:** a regra compara extensão e categoria geral, sem certificar equivalência física. [Congelamento](results/congelamento.json) e [todos os pares candidatos](results/pares_candidatos.csv).

7. **Recorte e preparação.** Cada segmento teve 2.048 amostras, com pelo menos 256 px de contexto bilateral para o maior evento permitido. Recortou-se antes de inverter a intensidade, remover tendência/média e padronizar. Não houve deformação do eixo, Hann ou ajuste do alinhamento por WCO. **Limitação:** {c['eventos_com_vizinhos_aceitos']} segmentos contêm outros eventos aceitos e há sobreposição; o método caracteriza o segmento contextualizado.

8. **Equações e resumo fixados.** Foram usados K=1, Morlet ω0=6, correlação enviesada 1/N, 200 ondas de 50–1.024 px, WPS alternativa em módulo, WCS complexa e denominador da WCO no lag zero. A média principal teve pesos iguais nas escalas logarítmicas e nos 65 lags −32…+32. Foram exportados suporte, denominadores, média conservadora e máximos. **Limitação:** banda e resumo diferem da análise de panoramas; regras de suporte não são testes de significância.

9. **Testes pequenos e dimensionamento.** Foram executados controles analíticos e sintéticos antes do lote. O erro da convolução independente foi {tests['erro_convolucao_independente']:.2e}; o atraso conhecido foi recuperado com o sinal correto. O benchmark mediu cerca de 0,16–0,19 s por par em N=2.048/200 escalas. A estimativa de 7–13 minutos foi informada antes do lote. **Limitação:** esse orçamento excluiu redação, revisão visual e renderização. [Testes](results/validacao_numerica.json), [benchmark](results/benchmark.csv) e [estimativa](results/estimativa.json).

10. **Lote experimental e controles.** Foram calculados 165 pares de eventos em oito versões (1.320 comparações); 975 controles fixos principais nos 325 pares de aquisições; mais 135 controles fixos para as variantes v5/faixa30–49/remoção somente da média. Também foram feitos 26 selfs e cinco pares dentro da mesma aquisição. O lote levou {f(timing['segundos_total']/60,1)} minutos. **Limitação:** os controles fixos podem conter eventos e não têm a mesma seleção morfológica; não são negativos físicos garantidos. O tempo ativo de anotação não foi cronometrado separadamente; não deve ser inferido do tempo de cálculo. [Tempos e memória](results/execucao.json).

11. **Agregação e interpretação.** Os escores foram agregados primeiro por par de aquisições, depois por padrão; a versão balanceada atribuiu peso igual aos padrões com resultados. Equivalentes foram contrastados com controles dos mesmos pares de aquisições. A retirada de uma aquisição removeu todos os seus pares. **Limitação:** a análise de slug fica sem contraste quando se retira qualquer uma de suas duas aquisições do único par equivalente; valores ausentes foram mantidos como ausentes. Não foram tratados eventos como réplicas nem calculados valores-p sem modelo nulo adequado.

12. **Relatos, gráficos, preservação e renomeação.** Foram preparados resumo executivo, PDF, figuras com os cinco padrões, exemplos favorável/desfavorável, catálogo e scripts. Os extremos foram usados apenas como ilustrações posteriores. A pasta de trabalho foi renomeada para **v16.2_wco_eventos**, conforme solicitado pelo usuário, e o atalho v16.2 foi removido. A verificação final conferiu os hashes das fontes, o congelamento, os cálculos dos mapas, os vínculos e as imagens. Também foram inspecionadas visualmente a capa, os gráficos e a renderização do PDF. **Limitação:** src2 continua ignorado pelo Git; esta entrega é local, sem commit ou sincronização presumidos. [Verificação final](results/verificacao_final.json) e [manifesto](results/manifesto.json).

{lead}
'''
    write('PASSO_A_PASSO.md',steps)
    alltable=table(['Variante','Padrão','Equivalentes − distintos','Equivalentes − fixos'],[
        [v,LABELS[p],f(get(p,'distintos',v).diferenca),f(get(p,'fixos',v).diferenca)] for v in ev.variante.unique() for p in ['plug','slug']])
    legacy=[]
    for p in ['plug','slug']:
        # Mesmo suporte de pares e mesmos segmentos; separa o núcleo numérico.
        legacy.append([LABELS[p],f(get(p,var='v5_mesmos_segmentos').equivalentes),f(get(p).equivalentes),
                       f(get(p).equivalentes-get(p,var='v5_mesmos_segmentos').equivalentes)])
    results=f'''# Resultados detalhados da comparação de eventos

{lead}

## Domínio e cobertura

{covtable}

N=2.048; amostragem 1 px; 200 comprimentos de onda de 50–1.024 px; lags principais −32…+32. Peso igual por escala logarítmica e por lag. A média conservadora é outro domínio e aparece separadamente nas tabelas. Os três eventos sem qualquer parceiro equivalente continuam catalogados, inclusive quando contribuem como controles distintos.

## Contrastes principais

{eqtable}

Os equivalentes intrapadrão ocorreram em quatro pares de aquisições: três de plug e um de slug. Há 19 pares de eventos equivalentes, 46 distintos intrapadrão, 20 equivalentes entre padrões e 80 distintos entre padrões. O limite de amostragem foi aplicado sem métricas wavelet. Os pares intermediários e não selecionados permanecem em pares_candidatos.csv.

A versão balanceada intrapadrão atribui 50% do peso a plug e 50% a slug: diferença de 0,0897 frente a eventos distintos e 0,1240 frente a recortes fixos. Esse balanceamento dá ao único par de slug o mesmo peso dos três pares de plug; por isso não fortalece a evidência nem substitui os resultados separados por padrão.

## Todas as variantes predefinidas

{alltable}

As variantes de centros usam os mesmos controles fixos principais; as demais usam controles recalculados com o mesmo núcleo/faixa/preparo. Todos os pares de eventos foram recalculados em cada variante; não se escolheu a variante de maior valor por par. A coluna diferenca_conservador de [contrastes_por_padrao.csv](results/contrastes_por_padrao.csv) usa somente as escalas com suporte conservador e não substitui a banda principal.

## Separação do núcleo numérico

{table(['Padrão','v5 nos mesmos recortes','Morlet corrigida','Diferença'],legacy)}

A diferença acima usa a média dos equivalentes nos mesmos pares de aquisições e a mesma entrada 18–71. O v5 preserva sua wavelet cmor, escala e orientação corrigida apenas no eixo reportado. A validação literal compara suas correlações extraídas por AST e a CWT de convolução com o cálculo usado. Os sinais são preparados dentro de cada segmento de mesmo comprimento; não se reaproveita tendência estimada do panorama inteiro.

## Comparação histórica de panoramas

A base 18–71 apresentou escore médio nativo de aproximadamente **0,187651**, contra **0,184176** do núcleo v5 na mesma faixa. Esses números resumem máximos por escala, em todos os lags e ondas de 50–3.000 px, com pesos de WPS. **Não é válido subtraí-los da média de eventos deste estudo e chamar o resultado de ganho.** O [registro histórico estruturado](results/referencia_historica.json) preserva a fonte; os controles fixos de 2.048 px fornecem a comparação com o mesmo domínio desta análise.

## Dependência, suporte e controles negativos

Os resultados de retirada de aquisição estão em [retirada_aquisicao.csv](results/retirada_aquisicao.csv), com ausência explícita quando não sobra contraste. Não são intervalos de confiança. Não houve descarte de escalas por denominador nas comparações experimentais Morlet; o domínio principal teve suporte geométrico completo. O maior valor experimental de WCO Morlet foi **{f(c['maior_wco_morlet'],6)}**; no v5 foi **{f(c['maior_wco_v5'],6)}**. Não houve clipping.

Self em lag zero foi 1 por identidade algébrica, não como evidência de equivalência física. Os controles com eventos distintos e fundo comum estão em [controles_sinteticos.csv](results/controles_sinteticos.csv); comparações dentro da mesma aquisição estão em [controles_mesma_aquisicao.csv](results/controles_mesma_aquisicao.csv), fora dos resultados entre aquisições.

## Arquivos numéricos

- [Resultados dos 165 pares em oito versões](results/resultados_eventos.csv).
- [Resultados dos controles fixos](results/resultados_fixos.csv).
- [Agregação por par de aquisições](results/agregados_eventos.csv).
- [Contrastes com controles do mesmo par de aquisições](results/contrastes_pareados.csv).
- [Contrastes por padrão](results/contrastes_por_padrao.csv) e [balanceados](results/contrastes_balanceados.csv).
- [Distribuições](results/distribuicoes.csv), [cobertura](results/cobertura.csv) e [conclusão estruturada](results/conclusao.json).

As matrizes em maps são integrais nos exemplos preservados; os máximos foram calculados antes de reduzir a resolução somente para renderização. Os perfis descritivos dos 165 pares principais estão em results/perfis_eventos.npz. O catálogo, os CSVs e os scripts permitem recalcular qualquer par sem usar mapas de baixa resolução.
'''
    write('RESULTS.md',results)
    readme='''# Comparação de eventos por WCO alternativa

Entrega local da execução do plano v16.2, com base em **v16_spatial2_18_71**. A pasta foi renomeada para **v16.2_wco_eventos**, conforme o nome solicitado pelo usuário. O atalho v16.2 foi removido. Os dados e as análises anteriores permaneceram intactos.

Comece pelo [resumo executivo com gráficos](EXECUTIVE_SUMMARY.md) ou pelo [PDF compacto](RESUMO_EXECUTIVO.pdf). Veja o [passo a passo](PASSO_A_PASSO.md), os [resultados](RESULTS.md), o [protocolo congelado](PROTOCOL.md), as [decisões](DECISIONS.md) e o [plano original](PLAN_ORIGINAL.md).

A [nota sobre âncoras e autocorrelação](ANCORAS_E_AUTOCORRELACAO.md) explica os limites de x/y, o domínio de lag, a aplicação da CWT às correlações e as condições para uma possível extensão a churn e anular.

## Organização

- scripts: auditoria de entrada, detecção por imagem, pareamento sem wavelets, núcleo idêntico à base, validações, execução, agregação, gráficos, relatórios e verificador final.
- results: catálogo com rejeições, coordenadas e procedência; pares candidatos/selecionados; tabelas numéricas; congelamento SHA-256; validações, logs e manifesto.
- figures/observabilidade: atlas das 26 aquisições, cada um com três quadros originais fixos.
- figures/catalogo: todas as propostas das seis aquisições com morfologia alongada, inclusive rejeitadas.
- figures/eventos: evidência individual dos 26 eventos aceitos, com segmento, quadro da âncora e intensidade.
- maps: coeficientes complexos e WCO integrais dos exemplos, sem redução de lag.

## Reprodução

O Python padrão deste terminal não possuía as dependências. A execução usou **Python 3.13.2** em `/Users/i521862/.pyenv/versions/3.13.2/bin/python3`. As versões da distribuição e dos módulos estão no manifesto. PyWavelets informa distribuição 1.9.0 e módulo 1.8.0; o núcleo Morlet principal usa NumPy/SciPy. requirements.txt registra o ambiente declarado e não garante identidade binária; o manifesto registra hashes do código.

Na raiz do projeto, usando o ambiente indicado:

```sh
export MPLCONFIGDIR=/tmp/wco-eventos-mpl
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/inspect_inputs.py
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/catalog_events.py
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/validate_and_benchmark.py
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/run_analysis.py
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/summarize_results.py
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/render_figures.py
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/write_reports.py
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/verify_final.py
```

A revisão visual em results/revisao_visual.csv é um insumo congelado. Reproduzir o cálculo não produz uma nova anotação independente. Se o catálogo ou as regras forem alterados, o congelamento impedirá a execução silenciosa; uma análise nova deve ter outra versão e registrar que os resultados desta já são conhecidos. Os scripts não executam lotes da pasta de origem.

Para reconstruir somente figuras e relatórios, use render_figures.py, write_reports.py e verify_final.py com os CSVs existentes. O lote completo principal/sensibilidades levou aproximadamente sete minutos nesta máquina; use results/execucao.json como tempo observado e results/estimativa.json como orçamento anterior ao lote. Redação, inspeção visual, renderização e verificação são etapas separadas.

## Interpretação e conservação

Todos os valores são descritivos de intensidade de imagem, lag e escala em pixels. A conclusão é restrita à categoria óptica e ao contexto selecionados. Não há calibração espacial física, identidade entre aquisições, validação externa ou inferência por réplicas independentes de eventos.

O .gitignore existente ignora src2. Nenhum arquivo desta entrega foi presumido versionado e nenhum commit foi criado. Para conservar a análise, preserve esta pasta junto às entradas e à base indicadas no manifesto.
'''
    write('README.md',readme)
    # PDF compacto para leitura/compartilhamento; o Markdown contém os mapas completos.
    with PdfPages(ROOT/'RESUMO_EXECUTIVO.pdf',metadata={'Title':'Resumo executivo da comparação de eventos por wavelets','Author':'Análise local do projeto TCC','Subject':'Resultados e limitações em PT-BR'}) as pdf:
        fig=plt.figure(figsize=(8.27,11.69));fig.patch.set_facecolor('white')
        fig.text(.08,.94,'Comparação de eventos\npor WCO alternativa',fontsize=23,weight='bold',va='top',color='#153c58')
        blocks=[('Conclusão',lead.replace('**','')),
                ('Resultados',f'Plug: diferença de {f(plug.diferenca)} frente a eventos distintos, em três pares de aquisições. Slug: {f(slug.diferenca)}, sustentada por somente um par de eventos. Frente a recortes fixos: {f(plugfix.diferenca)} em plug e {f(slugfix.diferenca)} em slug.'),
                ('Cobertura','31 entradas conferidas; 26 aquisições elegíveis; 26 eventos aceitos em seis aquisições; apenas cinco aquisições participaram de equivalentes intrapadrão. Churn/anular tiveram delimitação não estabelecida; bolhas dispersas tem um único registro.'),
                ('Método','Base v16_spatial2_18_71, linhas 18–71. Catálogo congelado antes da WCO; N=2.048; 200 comprimentos de onda de 50–1.024 px; média em lags −32 a +32 px; K=1; estimador enviesado; sem clipping.'),
                ('Limitações decisivas','Seleção exploratória por um único agente, poucos registros, eventos e contextos sobrepostos, sem validação física das fronteiras. Coerência alta em controles de fundo comum impede usar WCO isolada como prova de equivalência. Não houve teste externo ou classificação validada.')]
        y=.82
        for heading,body in blocks:
            fig.text(.08,y,heading,fontsize=13,weight='bold',va='top',color='#147d92');y-=.03
            lines=textwrap.wrap(body,width=87)
            fig.text(.08,y,'\n'.join(lines),fontsize=10.5,va='top',linespacing=1.4);y-=.022*len(lines)+.035
        fig.text(.08,.035,'Detalhes, dados e reprodução: EXECUTIVE_SUMMARY.md, PASSO_A_PASSO.md e README.md',fontsize=8,color='#58636d')
        pdf.savefig(fig);fig.savefig(ROOT/'figures/capa_resumo.png',dpi=130);plt.close(fig)
        fig=plt.figure(figsize=(8.27,11.69));fig.patch.set_facecolor('white')
        fig.text(.08,.94,'Qual padrão é qual',fontsize=23,weight='bold',va='top',color='#153c58')
        fig.text(.08,.885,'Rótulos globais das aquisições; extraídos do inventário auditado.',fontsize=10,va='top')
        y=.83
        for p in patterns:
            fig.text(.08,y,pattern_names[p],fontsize=13,weight='bold',va='top',color='#147d92')
            n=int(coverage.loc[p,'aquisicoes']);rep=inventory[inventory.pattern==p].run.iloc[0]
            body=f'{n} aquisição(ões). Exemplo de imagem: {rep}. {roles[p]}.'
            fig.text(.08,y-.032,'\n'.join(textwrap.wrap(body,width=84)),fontsize=10.5,va='top',linespacing=1.4)
            y-=.117
        fig.text(.08,.23,'Padrões dos exemplos de pares no relatório completo',fontsize=12,weight='bold',va='top',color='#153c58')
        lines=[f'{case_names[k]}: {case_labels[k]}.' for k in case_ids]
        lines+=['x = linha contínua; y = linha tracejada. As duas linhas podem ser do mesmo padrão.',
                'Mapas WPS/WCS/WCO: a cor indica o valor da barra lateral, e o título identifica o padrão.']
        fig.text(.08,.195,'\n'.join('\n'.join(textwrap.wrap(line,width=86)) for line in lines),fontsize=10,va='top',linespacing=1.6)
        pdf.savefig(fig);fig.savefig(ROOT/'figures/guia_padroes.png',dpi=130);plt.close(fig)
        for filename,title in [('cobertura','Cobertura do protocolo'),('contrastes','Resultados pareados'),('sensibilidades','Estabilidade dos resultados'),('controles_sinteticos','Limitações evidenciadas pelos controles'),('cinco_padroes','Imagens, sinais e WPS dos cinco padrões')]:
            fig=plt.figure(figsize=(11.69,8.27));ax=fig.add_axes([.035,.08,.93,.82]);ax.imshow(plt.imread(ROOT/'figures'/f'{filename}.png'));ax.axis('off')
            fig.suptitle(title,fontsize=18,weight='bold',color='#153c58',y=.965)
            fig.text(.04,.035,'Resultados descritivos; os pares compartilham aquisições. Não são réplicas independentes.',fontsize=9)
            pdf.savefig(fig);plt.close(fig)
    print('Resumo executivo, PDF, passo a passo, resultados e README gerados.')
    extension=ROOT/'results/extensao_churn_anular/execucao.json'
    if extension.exists() and json.loads(extension.read_text()).get('status')=='PASS' and (ROOT/'results/extensao_churn_anular/resumo.json').exists():
        from write_texture_reports import main as integrate_texture
        integrate_texture()

if __name__=='__main__':main()
