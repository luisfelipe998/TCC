"""Write the honest performance conclusion, including the direct v5 comparator."""
import json
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from evaluate_preprocessing import summarize,MAIN,METRICS
from preprocessing_experiment import extract
from wavelet_core import pair,prepare
from run_analysis import plot_pair,small_maps,LAM
from finish_outputs import table
from plot_labels import PATTERN_PT, VARIANT_PT

ROOT=Path(__file__).resolve().parents[1]

def null_update():
    meta=json.loads((ROOT/'results/null_metadata.json').read_text())
    assert 'sem COI' in meta.get('domain',''), 'Conclua o diagnóstico nulo integral antes do executivo.'
    current=pd.read_csv(ROOT/'results/null_diagnostic.csv')
    assert len(current)==5 and (current.surrogates==199).all()
    with zipfile.ZipFile(ROOT/'results/historico_com_mascaras.zip') as archive:
        previous=pd.read_csv(archive.open('results/null_diagnostic.csv'))
    comparison=current.merge(previous,on=['run_x','run_y'],suffixes=('_current','_previous'),validate='one_to_one')
    rows=[[PATTERN_PT[r.pattern_current],f'{r.conditional_monte_carlo_p_previous:.3f}'.replace('.',','),
           f'{r.conditional_monte_carlo_p_current:.3f}'.replace('.',',')] for r in comparison.itertuples()]
    text='''## Diagnóstico condicional após a retirada das máscaras

O diagnóstico com 199 pares de fases aleatórias por exemplo também foi recalculado. A estatística atual foi o máximo da alternative WCO em todos os lags e nas 200 escalas, tanto nos pares observados quanto nos simulados. A tabela compara valores-p de **estatísticas com domínios diferentes**; os valores antigos são apenas uma referência histórica.

'''+table(['Exemplo','Valor-p anterior, com máscara','Valor-p atual, plano completo'],rows)+'''

O indício anterior do exemplo churn não se manteve: o valor-p passou de 0,035 para 0,515. No exemplo plug, ocorreu o inverso, de 0,700 para 0,005. O segundo resultado indica um máximo incomum sob o modelo de fases aleatórias aplicado ao domínio completo. Ele não comprova identidade de eventos, ligação física entre aquisições ou separação dos padrões. Os cinco valores-p são condicionais ao modelo, sem decisão de rejeição corrigida entre exemplos; a resolução de Monte Carlo é 0,005. Não foram adicionados contornos de significância aos mapas.

Os máximos observados, os quantis simulados e as limitações do modelo estão em [RESULTS.md](RESULTS.md) e [null_diagnostic.csv](results/null_diagnostic.csv).

'''
    return text,current.to_dict(orient='records')

def main():
    null_text,null_records=null_update()
    screen=pd.read_csv(ROOT/'results/preprocessing_screen_pairs.csv')
    confirm=pd.read_csv(ROOT/'results/preprocessing_confirm_pairs.csv')
    legacy=pd.read_csv(ROOT/'results/preprocessing_v5_pairs.csv')
    assert all('max_wco_all' in frame and 'max_wco_supported' not in frame for frame in [screen,confirm,legacy]), 'Atualize as estatísticas integrais antes do relatório.'
    combined=pd.concat([confirm,legacy],ignore_index=True)
    combined.to_csv(ROOT/'results/preprocessing_comparison_pairs.csv',index=False)
    ranks=[]
    for metric in METRICS:
        s=summarize(combined,metric)
        for key in ['same','different','gap']:s['delta_vs_v5_'+key]=s[key]-s.loc['v5',key]
        ranks.append(s.reset_index().assign(metric=metric))
    pd.concat(ranks,ignore_index=True).to_csv(ROOT/'results/preprocessing_comparison_rankings.csv',index=False)
    main=summarize(combined,MAIN);mx=summarize(combined,'max_wco_all')
    pivot=combined[combined.same_pattern].pivot(index=['run_x','run_y','pair_group'],columns='variant',values=MAIN)
    pivot['delta_vs_v5']=pivot.central_band20-pivot.v5
    grouped=pivot.groupby('pair_group').agg(base=('v5','mean'),candidate=('central_band20','mean'),
        delta=('delta_vs_v5','mean'),pairs=('delta_vs_v5','size'))
    grouped.to_csv(ROOT/'results/preprocessing_pattern_tradeoffs.csv')
    pivot.reset_index().to_csv(ROOT/'results/preprocessing_pair_tradeoffs.csv',index=False)
    improved=int((pivot.delta_vs_v5>0).sum());worse=int((pivot.delta_vs_v5<0).sum())
    poolbase=float(pivot.v5.mean());poolcandidate=float(pivot.central_band20.mean())
    pooled_direction='redução' if poolcandidate<poolbase else ('aumento' if poolcandidate>poolbase else 'ausência de mudança')
    poolmedian=float(pivot.core_row_median.mean())
    poolcorrected=float(pivot.baseline.mean())
    median_delta=pivot.core_row_median-pivot.v5
    median_improved=int((median_delta>0).sum());median_worsened=int((median_delta<0).sum())
    median_direction='redução' if poolmedian<poolbase else ('aumento' if poolmedian>poolbase else 'ausência de mudança')
    pooled_records=[]
    for variant in ['v5','baseline','central_band20','core_row_median']:
        delta=pivot[variant]-pivot.v5
        pooled_records.append(dict(variant=variant,mean_score=float(pivot[variant].mean()),
            delta_vs_v5=float(delta.mean()),pairs_improved=int((delta>0).sum()),
            pairs_worsened=int((delta<0).sum()),pairs_unchanged=int((delta==0).sum()),
            balanced_mean=float(main.loc[variant,'same'])))
    pd.DataFrame(pooled_records).to_csv(ROOT/'results/preprocessing_pooled_comparison.csv',index=False)
    # Evaluate the choices made without the held-out acquisition against v5 too.
    folds=pd.read_csv(ROOT/'results/preprocessing_confirm_heldout.csv');held=[]
    for row in folds.itertuples():
        test=combined[(combined.run_x==row.run)|(combined.run_y==row.run)]
        s=summarize(test,MAIN);b=s.loc['v5'];c=s.loc[row.selected]
        held.append(dict(run=row.run,pattern=row.pattern,selected=row.selected,
            delta_same=float(c.same-b.same),delta_different=float(c.different-b.different),
            delta_gap=float(c.gap-b.gap)))
    held=pd.DataFrame(held);held.to_csv(ROOT/'results/preprocessing_heldout_vs_v5.csv',index=False)
    hd=held[held.delta_same.notna()].groupby('pattern')[['delta_same','delta_different','delta_gap']].mean().mean()
    folder=ROOT/'figures/preprocessing';folder.mkdir(exist_ok=True)
    fig,axes=plt.subplots(1,3,figsize=(19,7),layout='constrained')
    sr=summarize(screen,MAIN);variants=list(sr.index)
    x=np.arange(len(variants))
    axes[0].bar(x-.18,sr.delta_same,width=.36,label='Variação intrapadrão')
    axes[0].bar(x+.18,sr.delta_different,width=.36,label='Variação entre padrões')
    axes[0].set_xticks(x,labels=[VARIANT_PT[v] for v in variants],rotation=55,ha='right');axes[0].axhline(0,color='k',lw=.6)
    axes[0].tick_params(axis='x',labelsize=8)
    axes[0].set(ylabel='Variação em relação à referência corrigida do v16',title='Triagem: 64 escalas, peso igual por grupo');axes[0].legend(fontsize=8)
    upper=pivot[['v5','central_band20','core_row_median']].max().max()*1.05
    for ax,variant,gains,losses,title in [(axes[1],'central_band20',improved,worse,'Faixa central 30–49'),
                                         (axes[2],'core_row_median',median_improved,median_worsened,'Mediana das linhas 20–59')]:
        for pattern,g in pivot.groupby('pair_group'):
            ax.scatter(g.v5,g[variant],label=PATTERN_PT[pattern],s=20,alpha=.8)
        ax.plot([0,upper],[0,upper],color='.5',ls='--')
        ax.set(xlim=(0,upper),ylim=(0,upper),xlabel='Escore do v5',ylabel='Escore do candidato',
               title=f'{title}, 200 escalas\n{gains} ganhos / {losses} perdas intrapadrão')
        ax.legend(fontsize=8)
    fig.suptitle('Escores da alternative WCO: busca em todas as defasagens, sem nova fórmula de coerência')
    fig.savefig(folder/'preprocessing_comparison.png',dpi=140);plt.close(fig)
    # Full alternative WPS/WCS/WCO maps for predetermined first pairs, not best cases.
    for pattern,g in confirm[(confirm.variant=='central_band20')&confirm.same_pattern].groupby('pair_group'):
        row=g.sort_values(['run_x','run_y']).iloc[0]
        a=extract(row.run_x)['central_band20'];b=extract(row.run_y)['central_band20'];n=min(len(a),len(b))
        result=pair(prepare(a[:n]),prepare(b[:n]),LAM)
        plot_pair(small_maps(result),row.run_x,row.run_y,f'{PATTERN_PT[pattern]} | média exploratória das linhas 30–49',
                  folder/f'central_band20_{row.run_x}__{row.run_y}.png')
    rows=[]
    for variant in ['v5','baseline','central_band20','core_row_median']:
        r=main.loc[variant];q=mx.loc[variant]
        rows.append([VARIANT_PT[variant],f'{r.same:.6f}',f'{r.different:.6f}',f'{r.gap:.6f}',f'{q.same:.6f}',f'{q.different:.6f}'])
    ranking=summarize(screen,MAIN)
    screen_rows=[[VARIANT_PT[v],f'{r.same:.6f}',f'{r.different:.6f}',f'{r.delta_same:+.6f}',f'{r.delta_gap:+.6f}'] for v,r in ranking.iterrows()]
    grouprows=[[PATTERN_PT[k],int(r.pairs),f'{r.base:.6f}',f'{r.candidate:.6f}',f'{r.delta:+.6f}'] for k,r in grouped.iterrows()]
    pooledrows=[[VARIANT_PT[r['variant']],f"{r['mean_score']:.6f}",f"{r['delta_vs_v5']:+.6f}",
        r['pairs_improved'],r['pairs_worsened'],f"{r['balanced_mean']:.6f}"] for r in pooled_records]
    median_grouprows=[]
    for pattern,g in pivot.groupby('pair_group'):
        median_grouprows.append([PATTERN_PT[pattern],len(g),f'{g.v5.mean():.6f}',
            f'{g.core_row_median.mean():.6f}',f'{(g.core_row_median-g.v5).mean():+.6f}'])
    selection_counts='; '.join(f'{VARIANT_PT[k]}: {v}' for k,v in folds.selected.value_counts().items())
    screen_evaluation=json.loads((ROOT/'results/preprocessing_screen_evaluation.json').read_text())
    screen_choice=VARIANT_PT[screen_evaluation['choice_on_full_dataset']]
    positive_patterns=', '.join(PATTERN_PT[k] for k in grouped.index[grouped.delta>0]) or 'nenhum'
    negative_patterns=', '.join(PATTERN_PT[k] for k in grouped.index[grouped.delta<0]) or 'nenhum'
    definitions='''| Opção | Operação, antes das mesmas fórmulas alternativas |
| --- | --- |
| baseline | Média das linhas 20–59; corte comum; remoção da tendência linear; normalização. |
| demean_only | Mesmo sinal, removendo apenas a média, sem ajuste linear. |
| center_row | Apenas linha 40. |
| central_band20 | Média fixa das linhas 30–49; não escolhida separadamente para cada par. |
| core_row_median | Mediana entre as linhas 20–59 em cada coluna. |
| core_trimmed_mean20 | Média após retirar os 20% menores e 20% maiores valores entre as 40 linhas, por coluna. |
| spatial_gaussian2 | Filtro gaussiano com sigma=2 px no segmento de intensidade; bordas refletidas. |
| spatial_median5 | Mediana espacial de 5 pixels; bordas refletidas. |
| rank_uniform | Postos médios da intensidade divididos por N+1, antes da remoção da tendência linear. Empates são preservados. |
| winsor_1pct | Limitação pelos percentis 1 e 99 de cada segmento antes da remoção da tendência linear. |
'''
    text=f'''# Pré-processamento mantendo alternative WPS, alternative WCS e alternative WCO

**Estado atual:** a triagem, a comparação com v5 e a avaliação retirando uma aquisição foram recalculadas sem COI nem máscara conservadora. Os máximos por escala são buscados em todos os lags. As equações alternativas, os pesos e os pisos de potência de referência foram mantidos. Os resultados anteriores estão no [arquivo histórico com máscaras](results/historico_com_mascaras.zip).

## Conclusão

**Os ganhos frente ao v5 foram marginais no escore principal; não foi demonstrada melhora robusta e geral da similaridade nem melhor separação dos padrões.** A referência corrigida, com a média original 20–59, aumentou o escore nos 87 pares, mas a diferença média foi pequena: **{poolcorrected-poolbase:+.6f}**. O maior ganho na média desses pares entre as três opções confirmadas foi o da mediana entre linhas, **{poolmedian-poolbase:+.6f}**, com **{median_improved}** aumentos e **{median_worsened}** reduções individuais. A faixa 30–49 ganhou **{poolcandidate-poolbase:+.6f}** frente ao v5, mas apresentou diferença de **{poolcandidate-poolcorrected:+.6f}** frente à referência corrigida.

As duas alternativas foram opções do conjunto de confirmação fixado no estudo anterior. Os resultados foram apresentados sem selecionar apenas a opção ou a ponderação mais favorável nesta atualização.

## Protocolo

Foram considerados 26 registros, todos os 325 pares: 87 intrapadrão e 238 entre padrões. As definições de Schuck, o estimador enviesado, omega0=6, a amostragem em pixels e o sentido de lag do método corrigido foram mantidos. Não foram feitas deformações do eixo, buscas por segmentos favoráveis, mudanças de denominador ou clipping de alternative WCO.

Dez opções foram recalculadas em 64 escalas. O conjunto de três opções definido no estudo anterior, referência 20–59, faixa 30–49 e mediana das linhas 20–59, foi conferido em 200 escalas, juntamente com uma recomputação direta do v5 nos 325 pares. Esse conjunto reduzido foi fixado antes de observar a nova triagem. Não se trata de um ensaio prospectivo com protocolo registrado nem de um teste externo independente.

'''+definitions+'''

A transformação por postos tem motivação na discussão de distribuições em Grinsted et al.; aquela discussão é de coerência convencional. Sua utilização antes das equações alternativas aqui é uma experiência nossa, não uma prescrição de Schuck. Postos e limitação de extremos mudam a informação de intensidade. Não há garantia de que valores extremos sejam ruído: podem corresponder às próprias interfaces.

## Como resumimos os mapas sem mudar as equações

O escore principal é uma média ponderada dos máximos de alternative WCO ao longo de **todos os lags** em cada escala. Os pesos brutos são a média geométrica das alternative WPS no lag zero da referência corrigida; entram escalas de 50–3000 px em que cada alternative WPS da referência tem ao menos 5% de seu próprio máximo. **Os pesos brutos e os critérios de potência são fixos entre candidatos de cada par.** A proteção numérica pode excluir escalas com denominador inválido em uma variante e, nesse caso, alterar os pesos normalizados efetivos. Não se afirma que isso ocorreu na triagem, pois essa ocorrência não foi quantificada. Não há exclusão geométrica de bordas ou restrição por sobreposição. A operação resume a alternative WCO já calculada; não é apresentada como uma equação de Schuck nem como uma nova definição de coerência.

O critério de 5% é uma escolha de análise, não um limiar de significância nem uma máscara de borda. Foram conferidos também 1% e 10%, além da alternative WCO em zero e do máximo simples no plano completo. A proteção contra denominadores numericamente nulos foi preservada. Usar o máximo por escala permite atrasos diferentes em escalas diferentes; portanto o escore não comprova que exista um único atraso físico comum. Retirar as restrições permite encontrar máximos em mais posições e escalas. Um escore maior que na versão histórica não demonstra, por si só, melhor similaridade experimental.

Na média balanceada, cada uma das quatro classes com pares intrapadrão tem peso igual; para controles, cada uma das dez combinações entre classes tem peso igual. Isso impede que os padrões agitado e wavy annular dominem toda a análise, mas pode fazer um benefício em apenas três pares do padrão pistonado parecer maior do que na média dos 87 pares. Por isso ambas as agregações são explicitadas.

## Resultados da triagem — 64 escalas

'''+table(['Opção','Escore intrapadrão balanceado','Escore entre padrões','Variação intra vs referência corrigida','Variação do contraste intra−entre'],screen_rows)+f'''

A regra da triagem completa selecionou **{screen_choice}**. Essa escolha é exploratória e não redefiniu silenciosamente o conjunto de três opções da comparação com 200 escalas. A tabela apresenta também as variações nos controles entre padrões: elevar a média intrapadrão não basta se o mesmo efeito elevar ainda mais as comparações entre padrões.

## Comparação direta com v5 — 200 escalas

'''+table(['Opção','Escore intra','Escore entre','Intra−entre','Média dos máximos integrais intra','Média dos máximos integrais entre'],rows)+f'''

O escore intra balanceado do v5 foi **{main.loc['v5','same']:.6f}**, contra **{main.loc['central_band20','same']:.6f}** na faixa central. Na média simples dos 87 pares, os valores foram **{poolbase:.6f} → {poolcandidate:.6f}**, diferença de **{poolcandidate-poolbase:+.6f}**. Melhoraram **{improved}** pares e pioraram **{worse}**. O contraste entre a média intrapadrão e a média dos controles entre padrões foi **{main.loc['central_band20','gap']:+.6f}** para a faixa central; valores negativos indicam controles com média maior.

Os máximos simples e o escore ponderado são resumos diferentes. A busca no plano completo torna o máximo mais exposto a pontos extremos, inclusive às violações numéricas do v5. Os resultados de todos os cinco resumos estão em [preprocessing_comparison_rankings.csv](results/preprocessing_comparison_rankings.csv). Escolher apenas o resumo mais favorável não seria uma comparação suficiente.

O mesmo domínio de comparação significa os mesmos lags e os mesmos comprimentos de onda nominais. O comparador literal do v5 manteve sua conversão antiga de escalas e o preparo antes do corte, diferenças metodológicas documentadas em [EQUATION_AUDIT.md](EQUATION_AUDIT.md) e [REPORT.md](REPORT.md). A comparação com o v5 não isola apenas o efeito do pré-processamento.

### Peso igual para cada um dos 87 pares

'''+table(['Opção','Média dos 87 pares','Diferença vs v5','Pares com aumento','Pares com redução','Média balanceada por padrão'],pooledrows)+f'''

A mediana entre linhas apresentou escore médio **{poolmedian:.6f}**, com **{median_improved}** pares acima do v5 e **{median_worsened}** abaixo. A média por par e a média balanceada por padrão respondem a ponderações distintas e foram mantidas juntas para evitar atribuir aos 87 pares um benefício concentrado em classes pequenas.

Frente à referência corrigida do v16, a mediana entre linhas acrescentou somente **{poolmedian-poolcorrected:+.6f}** à média dos 87 pares. A faixa 30–49 apresentou diferença de **{poolcandidate-poolcorrected:+.6f}**. Portanto, o pequeno ganho desta faixa frente ao v5 não deve ser atribuído integralmente ao pré-processamento.

### Ganhos e perdas por padrão — faixa central contra v5

'''+table(['Padrão','Pares','Média v5','Média faixa central','Diferença'],grouprows)+f'''

Padrões com aumento da média neste escore: **{positive_patterns}**. Padrões com redução: **{negative_patterns}**. A comparação permaneceu por aquisições independentes, sem procurar os segmentos que produzissem maior coerência.

### Ganhos e perdas por padrão — mediana entre linhas contra v5

'''+table(['Padrão','Pares','Média v5','Média da mediana entre linhas','Diferença'],median_grouprows)+f'''

## Escolha sem usar a aquisição avaliada

Em cada rodada, retirou-se uma aquisição e todos os seus 25 pares antes de escolher o pré-processamento. Nesta atualização, a escolha foi realizada **entre as três opções fixadas, usando 200 escalas**, a partir de preprocessing_confirm_heldout.csv. A regra exigiu aumentar a média intrapadrão e melhorar o contraste intra−entre em relação à referência corrigida; caso contrário, manteve a referência. Contagem das escolhas: **{selection_counts}**. A aquisição do padrão dispersed bubble não tem comparador intrapadrão; por isso 25 aquisições entram no resumo intrapadrão.

Usando o v5 como comparador, as mudanças balanceadas nas quatro classes com réplicas foram: intrapadrão **{hd.delta_same:+.6f}**, entre padrões **{hd.delta_different:+.6f}**, contraste **{hd.delta_gap:+.6f}**. Os pares e as rodadas compartilham registros; isso não constitui um conjunto de teste independente nem justifica tratar os pares como 325 réplicas estatísticas. Essa avaliação do conjunto reduzido não valida como melhor uma opção que apareça na nova triagem de dez alternativas.

## O que poderia elevar valores sem demonstrar maior similaridade

- Ganho, inversão de sinal e normalização da amplitude cancelam na alternative WCO; não oferecem uma melhora intrínseca.
- Deslocar um sinal pode levar um máximo para lag zero. A alternative WCO já procura relações em lag; deslocar o mapa não cria uma relação mais forte.
- Estreitar excessivamente a banda ou elevar omega0 pode aumentar coerência em sinais não relacionados. No limite de um único bin espectral comum, a Eq. 32.17 tende a 1 porque numerador e denominador ficam proporcionais ao mesmo produto de amplitudes. Isso não seria uma melhoria de similaridade demonstrada.
- Alinhar/deformar o eixo ou selecionar o segmento com maior alternative WCO mudaria a pergunta e exigiria uma validação independente da correspondência. Não foi feito.
- Suavizar numerador/denominador da maneira de Grinsted mudaria a definição. Não foi feito nas equações alternativas.
- Estimar fases/contornos ou corrigir deslocamentos do panorama pode alterar a qualidade do observável, mas precisaria de validação da imagem, e não apenas de um aumento da alternative WCO. O registro m03q40_001 continua contendo 60 deslocamentos truncados; nenhum filtro simples testado resolve essa procedência física.

## Decisão

Manter as equações alternativas matematicamente verificadas é necessário, mas não faz duas aquisições independentes de um mesmo padrão terem as mesmas posições, fases ou espaçamentos. A referência corrigida apresentou ganho pequeno e sistemático no escore principal. A mediana entre linhas apresentou o maior ganho médio entre as opções confirmadas quando os 87 pares receberam peso igual, mas esse ganho permaneceu marginal e incluiu perdas em {median_worsened} pares. A faixa central foi mais favorável com peso igual por padrão. **Nenhuma dessas observações demonstrou melhora robusta e geral da similaridade nem separação convincente dos padrões.** A retirada das máscaras simplificou a apresentação e alterou os resumos estatísticos, sem criar correspondência entre eventos.

Dados: preprocessing_comparison_pairs.csv, preprocessing_comparison_rankings.csv, preprocessing_pattern_tradeoffs.csv e preprocessing_heldout_vs_v5.csv. Figura: [comparação dos candidatos](figures/preprocessing/preprocessing_comparison.png). Mapas completos do candidato para o primeiro par de cada classe estão em figures/preprocessing, escolhidos por ordem alfabética, não pelo maior ganho.
'''
    (ROOT/'PREPROCESSING.md').write_text(text)
    executive=f'''# Conclusão executiva — v16_spatial2

**Os ganhos frente ao v5 foram marginais no escore principal. Não foi demonstrada melhora robusta e geral da similaridade nem separação convincente dos padrões.** A referência corrigida do v16 aumentou esse escore nos 87 pares, com ganho médio de apenas **{poolcorrected-poolbase:+.6f}**. Entre as opções confirmadas, a mediana das linhas 20–59 teve o maior ganho médio com peso igual por par: **{poolmedian-poolbase:+.6f}**, mas **{median_worsened} dos 87 pares pioraram**.

A faixa central 30–49, destacada no estudo anterior, ganhou **{poolcandidate-poolbase:+.6f}** frente ao v5, mas ficou **{poolcorrected-poolcandidate:.6f}** abaixo da referência corrigida na média dos 87 pares. Portanto, não é a melhor opção segundo essa ponderação. Todos esses números foram recalculados sem COI nem máscara conservadora.

## Atualização da apresentação e dos cálculos

Todos os gráficos foram refeitos com os nomes dos padrões em inglês e os demais textos em PT-BR, sem o cinza do COI, a linha da máscara conservadora ou as curvas de médias restritas a essas regiões. A média apresentada da alternative WPS foi calculada sobre todos os lags. Os escores de pré-processamento e o diagnóstico com fases aleatórias foram recalculados no plano completo; não houve alteração das equações de alternative WPS, alternative WCS e alternative WCO.

Remover essas restrições permite pesquisar máximos em mais posições e escalas. Um aumento em relação ao antigo escore com máscara, sozinho, não demonstra maior similaridade. Nesta atualização, o v5 e os candidatos foram comparados nos mesmos lags e comprimentos de onda nominais; o v5 manteve sua conversão antiga de escalas e seu preparo antes do corte, diferenças já auditadas. A versão anterior foi preservada no [arquivo histórico com máscaras](results/historico_com_mascaras.zip).

## O que foi efetivamente verificado

Foram lidos os dois PDFs solicitados e três referências metodológicas complementares, conferidas visualmente as equações de Schuck, reproduzidos exatamente os 31 panoramas e verificadas as planilhas. A análise experimental usa os 26 registros confirmados, com 87 pares intrapadrão e 238 controles entre padrões. As equações alternativas permaneceram em uso durante todos os experimentos.

As definições de **alternative WPS (32.14–32.15), alternative WCS (32.16) e alternative WCO (32.17)** foram preservadas. O estimador enviesado com 1/N continua sendo uma adaptação explícita: o artigo usa 1/(N−|m|) em 32.21–32.23. A versão literal não enviesada também foi calculada e apresentou instabilidade em grandes lags. O denominador da alternative WCO permanece no lag zero.

O v16 corrigiu a discretização da Morlet, a conversão de escala, a orientação do lag e o preparo após o corte comum. A versão corrigida respeitou os limites da alternative WCO nos testes e nos pares avaliados em todo o plano, sem limitar artificialmente os resultados a 1. **Essa correção matemática e numérica não comprova maior similaridade física entre registros.** A retirada das máscaras não foi responsável por essa correção.

## Comparação atual com o v5

Dez opções de pré-processamento foram recalculadas em 64 escalas, conservando as fórmulas. A referência 20–59, a faixa 30–49 e a mediana entre linhas foram conferidas em 200 escalas. A faixa 30–49 permaneceu como comparação definida no estudo anterior; não foi escolhida novamente a partir dos resultados sem máscaras.

'''+table(['Opção','Média dos 87 pares','Diferença vs v5','Pares com aumento','Pares com redução','Média balanceada por padrão'],pooledrows)+f'''

O contraste intrapadrão menos entre padrões foi **{main.loc['central_band20','gap']:+.6f}** para a faixa central e **{main.loc['core_row_median','gap']:+.6f}** para a mediana entre linhas. Os dois valores negativos indicam controles entre padrões com média maior, sem separação convincente neste escore.

Os resultados também dependeram do resumo do mapa: por exemplo, a média balanceada dos máximos integrais passou de **{mx.loc['v5','same']:.3f}** no v5 para **{mx.loc['central_band20','same']:.3f}** na faixa central, uma redução. Esse máximo é sensível a pontos extremos, inclusive às violações numéricas do v5. Não deve ser confundido com o escore ponderado da tabela nem interpretado isoladamente como medida da qualidade do método.

O ganho da mediana entre linhas sobre a referência corrigida foi somente **{poolmedian-poolcorrected:+.6f}**. A faixa central permaneceu mais favorável quando cada padrão recebeu o mesmo peso. Essa inversão entre ponderações é parte da limitação do resultado, não evidência de uma opção universalmente superior.

O escore é um resumo dos mapas, **não uma nova fórmula de WCO**: máximos por escala em todos os lags, ponderados pelas alternative WPS da referência no lag zero. O piso principal de 5% do pico dessas alternative WPS foi mantido, com verificações a 1% e 10%; não é um COI nem uma máscara geométrica. Os pesos brutos são fixos por par; a proteção de denominadores pode alterar sua normalização efetiva se uma variante tiver escalas numericamente inválidas. Essa ocorrência não foi quantificada. A triagem selecionou **{screen_choice}** pela regra declarada. A avaliação retirando uma aquisição usou exclusivamente as três opções fixadas e 200 escalas. As classes possuem poucas aquisições e os pares não são réplicas independentes.

## Comprimento dos registros e janelas

**As janelas adicionais não são exigidas pelas equações de Schuck.** A análise principal usou o registro completo disponível, ou o prefixo comum de um par. O controle de segmentos de início, centro e fim permaneceu como diagnóstico separado de sensibilidade ao comprimento e à posição. A retirada do COI e da máscara conservadora não elimina possíveis efeitos de borda nem a menor quantidade de produtos que compõem a correlação em lags extremos.

## Figuras dos padrões e dos resultados

Os exemplos abaixo usam a entrada principal, média das linhas 20–59, e foram escolhidos pela ordem alfabética dos registros de cada padrão. São exemplos do conjunto, não os casos com maior coerência.

### plug — m05q01_001

![Panorama, sinal de intensidade e alternative WPS do registro plug m05q01_001](figures/representatives/m05q01_001.png)

*Figura 1 — Panorama original, sinal analisado, alternative WPS e sua média em lag. O eixo do sinal é posição no panorama; o eixo do mapa wavelet é lag da autocorrelação. Eles não representam a mesma coordenada.*

### slug — m04q01_001

![Panorama, sinal de intensidade e alternative WPS do registro slug m04q01_001](figures/representatives/m04q01_001.png)

*Figura 2 — Exemplo slug com as mesmas definições e processamento. Um pico de alternative WPS descreve uma escala do sinal; não identifica automaticamente o comprimento de uma bolha.*

### wavy annular — m03q10_001

![Panorama, sinal de intensidade e alternative WPS do registro wavy annular m03q10_001](figures/representatives/m03q10_001.png)

*Figura 3 — Exemplo de wavy annular. As escalas observadas e a aparência do panorama não demonstram, por si, alta alternative WCO entre aquisições independentes.*

### dispersed bubble — m07q01_001

![Panorama de dispersed bubble, sinal de intensidade e alternative WPS do registro m07q01_001](figures/representatives/m07q01_001.png)

*Figura 4 — Exemplo de dispersed bubble, com o panorama, o sinal extraído e a alternative WPS. Há somente um registro elegível desse padrão no escopo confirmado; portanto, não existe um segundo registro independente para apresentar uma comparação intrapadrão de alternative WCS e alternative WCO.*

### churn — m03q05_001

![Panorama, sinal de intensidade e alternative WPS do registro churn m03q05_001](figures/representatives/m03q05_001.png)

*Figura 5 — Exemplo churn com a mesma entrada e as mesmas equações alternativas. A média em lag e o valor no lag zero resumem aspectos diferentes da alternative WPS.*

### Comparação dos cinco padrões

![Perfis de alternative WPS e média em lag para os cinco padrões analisados](figures/pattern_profiles.png)

*Figura 6 — Mediana dos perfis normalizados por registro; as faixas mostram o mínimo e o máximo observados, não intervalos de confiança. A sobreposição entre padrões é visível. O padrão dispersed bubble tem somente um registro elegível. O perfil médio em lag inclui todos os lags.*

### Exemplos de alternative WCS e alternative WCO entre registros

As figuras 7–10 mostram o primeiro par em ordem alfabética de cada padrão com pelo menos duas aquisições elegíveis. Usam a implementação corrigida do v16, a correlação enviesada e a média original das linhas 20–59. No topo estão as alternative WPS individuais e suas médias sobre todos os lags; no centro, o módulo da alternative WCS; abaixo, a alternative WCO da Eq. 32.17, com denominador no lag zero. Os mapas exibem o plano completo, sem COI nem máscara conservadora.

#### plug — m05q01_001 × m06q01_001

![Alternative WPS, alternative WCS e alternative WCO entre os registros plug m05q01_001 e m06q01_001](figures/pairs/m05q01_001__m06q01_001.png)

*Figura 7 — Comparação de duas aquisições plug. O mapa de alternative WCS mostra magnitude de correlação filtrada por escala; a alternative WCO normaliza essa magnitude pelas alternative WPS no lag zero.*

#### slug — m04q01_001 × m05q05_001

![Alternative WPS, alternative WCS e alternative WCO entre os registros slug m04q01_001 e m05q05_001](figures/pairs/m04q01_001__m05q05_001.png)

*Figura 8 — Comparação de duas aquisições slug. As regiões de maior alternative WCO indicam associação em lag e escala; não localizam uma mesma bolha nas duas aquisições independentes.*

#### churn — m03q05_001 × m04q10_001

![Alternative WPS, alternative WCS e alternative WCO entre os registros churn m03q05_001 e m04q10_001](figures/pairs/m03q05_001__m04q10_001.png)

*Figura 9 — Comparação de duas aquisições churn. O exemplo apresenta os mapas calculados; a conclusão sobre ganhos e perdas frente ao v5 depende do conjunto de pares, não de um mapa isolado.*

#### wavy annular — m03q10_001 × m03q20_001

![Alternative WPS, alternative WCS e alternative WCO entre os registros anulares ondulados m03q10_001 e m03q20_001](figures/pairs/m03q10_001__m03q20_001.png)

*Figura 10 — Comparação de duas aquisições de wavy annular. Compartilhar o rótulo do padrão não exige que a alternative WCO seja alta em todas as escalas ou defasagens.*

### Ganhos e perdas frente ao v5

![Comparação dos pré-processamentos e dos escores de alternative WCO dos pares frente ao v5](figures/preprocessing/preprocessing_comparison.png)

*Figura 11 — À esquerda, a triagem com 64 escalas compara as opções com a referência corrigida do v16. Ao centro, a faixa 30–49 é comparada ao v5 em 200 escalas: {improved} pares acima da diagonal e {worse} abaixo. À direita, a mediana entre linhas apresenta {median_improved} aumentos e {median_worsened} reduções. Os eixos dos dois diagramas de dispersão têm os mesmos limites. Os escores usam todos os lags, sem COI nem máscara conservadora; não são novas fórmulas de coerência.*

## Decisão honesta

Recomenda-se a implementação corrigida como base matemática, reconhecendo o pequeno ganho sistemático no escore principal da entrada original. A mediana entre linhas e a faixa 30–49 permaneceram alternativas exploratórias: a primeira foi melhor na média por par; a segunda, na média com peso igual por padrão. **Os ganhos foram marginais e não sustentam a alegação de similaridade experimental robustamente maior ou de classificação melhor.**

O principal limite permanece: pertencer ao mesmo padrão não significa que duas aquisições independentes contenham os mesmos eventos ou a mesma organização espacial. Não seria correto forçar a coerência a subir por seleção de segmentos orientada pela própria coerência, deformação do eixo, redução artificial da banda ou mudança do denominador.

Os próximos caminhos continuam documentados apenas como planos: [v16.1](../v16.1/PLAN.md), comparação de perfis da alternative WPS, especialmente no lag zero; e [v16.2](../v16.2/PLAN.md), comparação de eventos equivalentes definidos por critérios de imagem antes de observar a coerência. Nenhuma dessas novas análises foi executada nesta atualização. A segunda responderia sobre eventos selecionados, não sobre a similaridade dos panoramas completos.

Detalhes: [equações conferidas](EQUATION_AUDIT.md), [experimento e comparação direta com v5](PREPROCESSING.md), [relatório metodológico](REPORT.md), [resultados gerais](RESULTS.md) e [galeria](figures/index.html).
'''
    executive=executive.replace('## Comprimento dos registros e janelas',null_text+'## Comprimento dos registros e janelas')
    (ROOT/'EXECUTIVE_SUMMARY.md').write_text(executive)
    summary=dict(conclusion='Ganhos marginais no escore principal frente ao v5, sem demonstração de melhora robusta e geral da similaridade ou da separação dos padrões. Referência corrigida: pequeno ganho nos 87 pares. Mediana entre linhas: maior média por par entre as opções confirmadas, com perdas individuais.',
        domain='Todos os lags; 50–3000 px; pisos de potência da referência de 1%, 5% e 10% preservados.',
        primary_score=MAIN,screen_choice=screen_evaluation['choice_on_full_dataset'],
        heldout_selection_source='preprocessing_confirm_heldout.csv: três opções fixadas, 200 escalas',
        candidate='central_band20',same_pattern_pairs_improved=improved,same_pattern_pairs_worsened=worse,
        pooled_same_pattern_score_v5=poolbase,pooled_same_pattern_score_candidate=poolcandidate,
        pooled_same_pattern_delta_vs_v5=poolcandidate-poolbase,pooled_change=pooled_direction,
        corrected_reference=dict(pooled_mean=poolcorrected,delta_vs_v5=poolcorrected-poolbase,
            pairs_improved=int(((pivot.baseline-pivot.v5)>0).sum()),
            pairs_worsened=int(((pivot.baseline-pivot.v5)<0).sum())),
        median_between_rows=dict(pooled_mean=poolmedian,delta_vs_v5=poolmedian-poolbase,
            pairs_improved=median_improved,pairs_worsened=median_worsened,
            balanced_mean=float(main.loc['core_row_median','same']),
            balanced_same_minus_different=float(main.loc['core_row_median','gap'])),
        balanced_same_pattern_score_v5=float(main.loc['v5','same']),balanced_same_pattern_score_candidate=float(main.loc['central_band20','same']),
        heldout_vs_v5={k:float(v) for k,v in hd.items()})
    summary['conditional_null_full_domain']=null_records
    (ROOT/'results/preprocessing_conclusion.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
