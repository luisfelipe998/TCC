"""Write the honest performance conclusion, including the direct v5 comparator."""
import json
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
    previous=pd.read_csv(ROOT.parent/'v16_spatial2/results/null_diagnostic.csv')
    assert len(current)==5 and (current.surrogates==199).all()
    comparison=current.merge(previous,on=['run_x','run_y'],suffixes=('_current','_previous'),validate='one_to_one')
    assert len(comparison)==5
    rows=[[PATTERN_PT[r.pattern_current],f'{r.conditional_monte_carlo_p_previous:.3f}'.replace('.',','),
           f'{r.conditional_monte_carlo_p_current:.3f}'.replace('.',','),
           f'{r.observed_max_previous:.6f}',f'{r.observed_max_current:.6f}'] for r in comparison.itertuples()]
    plug=comparison[comparison.pattern_current=='plug'].iloc[0]
    changes=pd.read_csv(ROOT/'results/band_change_pairs.csv')
    group=changes[(changes.variant=='baseline')&changes.same_pattern_new&(changes.pair_group_new=='plug')]
    example=group[(group.run_x==plug.run_x)&(group.run_y==plug.run_y)].iloc[0]
    text='''## Diagnóstico condicional nas duas faixas de extração

Foram recalculados 199 pares de fases aleatórias para cada um dos cinco exemplos. A estatística foi o máximo da alternative WCO em todos os lags e nas 200 escalas. Tanto a referência 20–59 quanto a nova extração 18–71 usaram esse domínio completo, sem COI nem máscara conservadora. Os sinais observados e simulados da nova faixa foram recalculados; os resultados antigos permaneceram na pasta de origem.

'''+table(['Exemplo','Valor-p, entrada 20–59','Valor-p, entrada 18–71','Máximo observado, 20–59','Máximo observado, 18–71'],rows)+f'''

**Uma mudança pequena na média dos 87 pares não significa que cada par mudou pouco.** No exemplo plug predeterminado, {plug.run_x} × {plug.run_y}, o máximo no plano completo caiu de **{plug.observed_max_previous:.6f}** para **{plug.observed_max_current:.6f}**, uma redução substancial de **{plug.observed_max_previous-plug.observed_max_current:.6f}**. Seu valor-p condicional passou de **{plug.conditional_monte_carlo_p_previous:.3f}** para **{plug.conditional_monte_carlo_p_current:.3f}**. O par foi escolhido pela ordem dos registros, antes de observar esses resultados.

O escore ponderado desse mesmo par, usando os pesos fixos 20–59, passou de **{example.old_score:.6f}** para **{example.new_score:.6f}**. Já a média dos escores ponderados dos **{len(group)} pares de plug** aumentou de **{group.old_score.mean():.6f}** para **{group.new_score.mean():.6f}**, porque os outros pares compensaram a perda. O escore ponderado por escala, o máximo simples no plano e o valor-p condicional são resumos diferentes; a melhora da média de um padrão não exige melhora em todos os pares nem em todos esses resumos.

Os valores-p são condicionais ao modelo estacionário de fases aleatórias. A diferença entre dois valores-p não é um teste de superioridade de uma faixa de extração. Um máximo incomum nesse modelo não comprova identidade de eventos, ligação física entre aquisições ou separação dos padrões. A resolução de Monte Carlo é 0,005; não foi tomada uma decisão de rejeição corrigida entre os cinco exemplos e não foram adicionados contornos de significância aos mapas.

Os quantis simulados e as limitações estão em [RESULTS.md](RESULTS.md) e [null_diagnostic.csv](results/null_diagnostic.csv).

'''
    return text,current.to_dict(orient='records')

def band_update():
    result=json.loads((ROOT/'results/band_comparison.json').read_text())
    assert result['status']=='PASS' and result['metric']==MAIN
    summary=pd.read_csv(ROOT/'results/band_comparison_summary.csv')
    changes=pd.read_csv(ROOT/'results/band_change_summary.csv')
    individual=pd.read_csv(ROOT/'results/band_change_pairs.csv')
    principal=changes[changes.metric==MAIN].set_index('variant')
    scores=summary[summary.metric==MAIN].set_index(['band','variant'])
    names={'baseline':'Média da faixa principal; Schuck corrigido',
           'v5':'Implementação v5 sobre a faixa principal',
           'central_band20':'Controle fixo: média 30–49',
           'core_row_median':'Mediana entre as linhas da faixa principal'}
    rows=[]
    for variant in ['v5','baseline','central_band20','core_row_median']:
        r=principal.loc[variant]
        rows.append([names[variant],f'{r.old_pooled_same:.6f}',f'{r.new_pooled_same:.6f}',
                     f'{r.delta_pooled_same:+.6f}',int(r.same_pairs_improved),
                     int(r.same_pairs_worsened),int(r.same_pairs_unchanged)])
    r=principal.loc['baseline'];old=scores.loc[('20_59','baseline')];new=scores.loc[('18_71','baseline')]
    relative=100*r.delta_pooled_same/r.old_pooled_same
    direction='aumentou' if r.delta_pooled_same>0 else ('diminuiu' if r.delta_pooled_same<0 else 'não alterou')
    delta_vs_v5=result['current_delta_vs_original_v5']['baseline']
    relative_vs_v5=100*delta_vs_v5/result['original_v5_pooled_same']
    conclusion=(f'A ampliação da média de 20–59 para 18–71 {direction} o escore médio de alternative WCO dos 87 pares '
        f'intrapadrão de {r.old_pooled_same:.6f} para {r.new_pooled_same:.6f}, diferença de '
        f'{r.delta_pooled_same:+.6f} ({relative:+.2f}%), usando os mesmos pesos brutos da referência 20–59. '
        f'Houve {int(r.same_pairs_improved)} aumentos, {int(r.same_pairs_worsened)} reduções e '
        f'{int(r.same_pairs_unchanged)} empates. Frente ao v5 original com entrada 20–59 '
        f'({result["original_v5_pooled_same"]:.6f}), a diferença da média 18–71 corrigida foi '
        f'{delta_vs_v5:+.6f} ({relative_vs_v5:+.2f}%). '
        f'\n\nO contraste balanceado intrapadrão menos entre padrões passou de {old.balanced_gap:+.6f} '
        f'para {new.balanced_gap:+.6f}.')
    if r.delta_pooled_same<=0:
        conclusion+=' Não foi obtida uma melhora da média dos 87 pares ao ampliar a faixa principal.'
    elif new.balanced_gap<=0:
        conclusion+=' O aumento do escore intrapadrão não produziu separação por esse resumo: a média balanceada entre padrões permaneceu maior ou igual à intrapadrão.'
    else:
        conclusion+=' O aumento observado é descritivo e interno; não estabelece superioridade em novos registros.'
    per_pattern=individual[(individual.variant=='baseline')&individual.same_pattern_new].groupby('pair_group_new').agg(
        pairs=('delta','size'),previous=('old_score','mean'),current=('new_score','mean'),delta=('delta','mean'))
    assert per_pattern.pairs.sum()==87
    pattern_rows=[[PATTERN_PT[name],int(q.pairs),f'{q.previous:.6f}',f'{q.current:.6f}',f'{q.delta:+.6f}']
                  for name,q in per_pattern.iterrows()]
    gains=', '.join(PATTERN_PT[name] for name in per_pattern.index[per_pattern.delta>1e-12]) or 'nenhum'
    losses=', '.join(PATTERN_PT[name] for name in per_pattern.index[per_pattern.delta<-1e-12]) or 'nenhum'
    conclusion+=(f' Os padrões com aumento da média foram {gains}; com redução, {losses}. '
                 f'Com peso igual por padrão, a média intrapadrão passou de {old.balanced_same:.6f} '
                 f'para {new.balanced_same:.6f}, variação de {new.balanced_same-old.balanced_same:+.6f}.')
    # Qualitative interpretation reviewed against this frozen 18–71 experiment;
    # the figures above retain the exact magnitudes and both aggregations.
    conclusion+=('\n\nNesta execução, a piora na média dos 87 pares foi pequena. A vantagem numérica remanescente '
                 'sobre o v5 original foi marginal em magnitude e não foi criada pela ampliação da faixa. '
                 'O resultado é misto entre padrões; não sustenta uma melhora global da alternative WCO. '
                 'Isso não implica mudanças pequenas em cada par: o exemplo plug predeterminado apresentou '
                 'uma redução substancial do máximo integral, detalhada no diagnóstico do resumo executivo.')
    sensitivity=[]
    for metric,label in [('weighted_scale_peak_power01','Máximos por escala, piso de 1%'),
                         (MAIN,'Máximos por escala, piso de 5%'),
                         ('weighted_scale_peak_power10','Máximos por escala, piso de 10%'),
                         ('weighted_zero_lag_wco','alternative WCO no lag zero'),
                         ('max_wco_all','Máximo no plano completo, sem pesos')]:
        q=changes[(changes.variant=='baseline')&(changes.metric==metric)].iloc[0]
        sensitivity.append([label,f'{q.old_pooled_same:.6f}',f'{q.new_pooled_same:.6f}',
                            f'{q.delta_pooled_same:+.6f}',f'{q.delta_balanced_gap:+.6f}'])
    text='''## Efeito da extração 18–71 frente à análise anterior 20–59

Esta comparação preserva, em cada par, os pesos brutos e os pisos de alternative WPS da referência corrigida **20–59**. A análise anterior e a nova usam todos os lags, sem COI nem máscara conservadora. Assim, a troca da referência de potência causada pela nova extração não é contabilizada como ganho da faixa. Os mapas e as equações alternativas continuam sendo os de cada entrada; somente seu resumo é calculado com a mesma referência. Pesos brutos fixos só garantem pesos normalizados idênticos quando as escalas elegíveis também têm denominadores numericamente válidos em todas as opções comparadas; as contagens da proteção numérica desta execução estão em [PREPROCESSING.md](PREPROCESSING.md).

Na tabela, a faixa principal é 20–59 na coluna anterior e 18–71 na coluna nova. O controle 30–49 permanece exatamente na mesma posição nas duas análises. Ganhos/perdas são diferenças superiores a 1e−12 em módulo; diferenças menores são empates numéricos.

'''+table(['Operação','Média dos 87 pares, 20–59','Média dos 87 pares, 18–71','Variação','Aumentos','Reduções','Empates'],rows)+f'''

O escore do controle fixo 30–49 reproduziu a análise anterior, com diferença absoluta máxima de {result['fixed_control_30_49_score_max_error']:.3g}. Isso verifica a continuidade desse controle sob os mesmos pesos. A comparação da média corrigida de 18–71 com o v5 original inclui tanto a ampliação da faixa quanto as correções da implementação; a comparação entre as duas médias corrigidas isola a mudança de entrada.

### Efeito por padrão — média corrigida, pesos fixos de 20–59

'''+table(['Padrão','Pares','Média 20–59','Média 18–71','Variação'],pattern_rows)+f'''

Padrões com aumento da média: **{gains}**. Padrões com redução: **{losses}**. A média conjunta dos 87 pares atribui a cada padrão peso proporcional à sua quantidade de pares. A média balanceada atribui peso igual a cada um dos quatro padrões com pares e passou de **{old.balanced_same:.6f}** para **{new.balanced_same:.6f}**. Por isso as duas agregações podem variar em sentidos opostos. O padrão dispersed bubble possui uma única aquisição elegível e não entra nessas médias intrapadrão. As duas agregações foram mantidas; os pares continuam compartilhando aquisições.

### Sensibilidade ao resumo do mapa — média corrigida

'''+table(['Resumo','Média dos 87 pares, 20–59','Média dos 87 pares, 18–71','Variação intra','Variação do contraste balanceado'],sensitivity)+'''

Os escores nas seções de pré-processamento abaixo usam os pesos nativos da referência **18–71**, fixos entre candidatos. Não devem ser subtraídos diretamente dos antigos escores 20–59 para atribuir um efeito à extração. As duas tabelas respondem a comparações distintas e mantêm as mesmas alternative WPS, alternative WCS e alternative WCO.

Dados da comparação: [pares com pesos fixos](results/band_comparison_pairs.csv), [resumos](results/band_comparison_summary.csv), [mudanças por resumo](results/band_change_summary.csv), [mudanças individuais](results/band_change_pairs.csv) e [verificação estruturada](results/band_comparison.json). Os pares compartilham aquisições; estas diferenças não constituem um teste independente de superioridade.

'''
    return text,conclusion,result

def denominator_update(screen,confirm,legacy):
    records=[]
    for label,frame,scales in [('Triagem',screen,64),('Confirmação Schuck',confirm,200),('Comparador v5, entrada 18–71',legacy,200)]:
        valid=frame.denominator_valid_scales.to_numpy()
        assert np.isfinite(valid).all() and ((valid>=0)&(valid<=scales)).all()
        records.append(dict(stage=label,pair_variant_evaluations=len(frame),scales_per_evaluation=scales,
            evaluations_with_invalid_scales=int((valid<scales).sum()),
            total_invalid_scales=int((scales-valid).sum()),
            minimum_valid_scales=int(valid.min()),
            minimum_eligible_scales_power05=int(frame.eligible_scales_power05.min()),
            maximum_eligible_scales_power05=int(frame.eligible_scales_power05.max())))
    count=sum(x['evaluations_with_invalid_scales'] for x in records)
    note=('Nenhuma combinação de par/opção teve escalas descartadas pelo denominador numérico nestes três lotes.'
          if count==0 else f'A proteção numérica descartou ao menos uma escala em {count} combinações de par/opção nos três lotes.')
    rows=[[r['stage'],r['pair_variant_evaluations'],r['scales_per_evaluation'],
           r['evaluations_with_invalid_scales'],r['total_invalid_scales'],
           f"{r['minimum_eligible_scales_power05']}–{r['maximum_eligible_scales_power05']}"] for r in records]
    text='''### Proteção numérica do denominador

'''+note+''' O descarte numérico e o piso de potência são critérios diferentes: o primeiro evita divisão por valores numericamente nulos; o segundo define as escalas usadas no resumo ponderado. Nenhum deles restringe geometricamente o mapa em lag. A contagem de escalas é somada sobre as combinações de par/opção, que não são observações independentes.

'''+table(['Lote','Combinações par/opção','Escalas por combinação','Combinações com descarte numérico','Total de escalas descartadas','Escalas no escore de 5%, mínimo–máximo'],rows)+'\n\n'
    return text,note,records

def main():
    null_text,null_records=null_update()
    screen=pd.read_csv(ROOT/'results/preprocessing_screen_pairs.csv')
    confirm=pd.read_csv(ROOT/'results/preprocessing_confirm_pairs.csv')
    legacy=pd.read_csv(ROOT/'results/preprocessing_v5_pairs.csv')
    assert all('max_wco_all' in frame and 'max_wco_supported' not in frame for frame in [screen,confirm,legacy]), 'Atualize as estatísticas integrais antes do relatório.'
    denominator_text,denominator_note,denominator_records=denominator_update(screen,confirm,legacy)
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
    for ax,variant,gains,losses,title in [(axes[1],'central_band20',improved,worse,'Faixa fixa 30–49'),
                                         (axes[2],'core_row_median',median_improved,median_worsened,'Mediana das linhas 18–71')]:
        for pattern,g in pivot.groupby('pair_group'):
            ax.scatter(g.v5,g[variant],label=PATTERN_PT[pattern],s=20,alpha=.8)
        ax.plot([0,upper],[0,upper],color='.5',ls='--')
        ax.set(xlim=(0,upper),ylim=(0,upper),xlabel='Escore da implementação v5 (18–71)',ylabel='Escore do candidato',
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
    best_variant=max(['baseline','central_band20','core_row_median'],key=lambda v:float(pivot[v].mean()))
    native_conclusion=(f'Na comparação entre implementações, usando os pesos nativos da referência 18–71, '
        f'o escore médio dos 87 pares foi {poolbase:.6f} para a implementação v5 aplicada à nova entrada '
        f'e {poolcorrected:.6f} para a referência corrigida de Schuck, diferença de {poolcorrected-poolbase:+.6f}. '
        f'O controle fixo 30–49 teve média {poolcandidate:.6f}, com {improved} aumentos e {worse} reduções frente ao v5 aplicado a 18–71. '
        f'A mediana da faixa 18–71 teve média {poolmedian:.6f}, com {median_improved} aumentos e {median_worsened} reduções. '
        f'A maior média por par entre as três opções confirmadas foi de {VARIANT_PT[best_variant]}. '
        'Esta comparação entre implementações e controles não mede, sozinha, a vantagem da nova faixa sobre a antiga 20–59.')
    contrast_note=('Todas as opções confirmadas apresentaram média entre padrões maior ou igual à intrapadrão neste escore.'
                   if (main.loc[['baseline','central_band20','core_row_median'],'gap']<=0).all()
                   else 'Ao menos uma opção apresentou média intrapadrão maior que a média entre padrões neste escore; isso não constitui uma validação de classificação.')
    band_text,band_conclusion,band_records=band_update()
    definitions='''| Opção | Operação, antes das mesmas fórmulas alternativas |
| --- | --- |
| baseline | Média das linhas 18–71; corte comum; remoção da tendência linear; normalização. |
| demean_only | Mesmo sinal, removendo apenas a média, sem ajuste linear. |
| center_row | Apenas linha 40. |
| central_band20 | Média fixa das linhas 30–49; não escolhida separadamente para cada par. |
| core_row_median | Mediana entre as linhas 18–71 em cada coluna. |
| core_trimmed_mean20 | Média após retirar 10 valores de cada extremo entre as 54 linhas, por coluna; floor(0,20 × 54), mantendo 34 valores. |
| spatial_gaussian2 | Filtro gaussiano com sigma=2 px no segmento de intensidade; bordas refletidas. |
| spatial_median5 | Mediana espacial de 5 pixels; bordas refletidas. |
| rank_uniform | Postos médios da intensidade divididos por N+1, antes da remoção da tendência linear. Empates são preservados. |
| winsor_1pct | Limitação pelos percentis 1 e 99 de cada segmento antes da remoção da tendência linear. |
'''
    text=f'''# Pré-processamento mantendo alternative WPS, alternative WCS e alternative WCO

**Estado atual:** a entrada principal passou a ser a média das linhas 18–71. A triagem, a comparação com a implementação v5 e a avaliação retirando uma aquisição foram recalculadas sem COI nem máscara conservadora. Os máximos por escala são buscados em todos os lags. As equações alternativas e as regras dos pesos/pisos foram mantidas; a comparação entre faixas usa pesos fixos de 20–59, enquanto a comparação entre pré-processamentos usa a referência 18–71. A análise anterior está preservada em [v16_spatial2](../v16_spatial2/EXECUTIVE_SUMMARY.md).

## Conclusão

{band_conclusion}

{native_conclusion}

As duas alternativas foram opções do conjunto de confirmação fixado no estudo anterior. Os resultados foram apresentados sem selecionar apenas a opção ou a ponderação mais favorável nesta atualização.

## Protocolo

Foram considerados 26 registros, todos os 325 pares: 87 intrapadrão e 238 entre padrões. As definições de Schuck, o estimador enviesado, omega0=6, a amostragem em pixels e o sentido de lag do método corrigido foram mantidos. Não foram feitas deformações do eixo, buscas por segmentos favoráveis, mudanças de denominador ou clipping de alternative WCO.

Dez opções foram recalculadas em 64 escalas. O conjunto de três opções definido no estudo anterior, referência 18–71, faixa 30–49 e mediana das linhas 18–71, foi conferido em 200 escalas, juntamente com a implementação v5 aplicada à entrada 18–71 nos 325 pares. Esse conjunto reduzido foi fixado antes de observar a nova triagem. Não se trata de um ensaio prospectivo com protocolo registrado nem de um teste externo independente.

'''+definitions+'''

A transformação por postos tem motivação na discussão de distribuições em Grinsted et al.; aquela discussão é de coerência convencional. Sua utilização antes das equações alternativas aqui é uma experiência nossa, não uma prescrição de Schuck. Postos e limitação de extremos mudam a informação de intensidade. Não há garantia de que valores extremos sejam ruído: podem corresponder às próprias interfaces.

## Como resumimos os mapas sem mudar as equações

O escore principal é uma média ponderada dos máximos de alternative WCO ao longo de **todos os lags** em cada escala. Os pesos brutos são a média geométrica das alternative WPS no lag zero da referência corrigida; entram escalas de 50–3000 px em que cada alternative WPS da referência tem ao menos 5% de seu próprio máximo. **Os pesos brutos e os critérios de potência são fixos entre candidatos de cada par.** A proteção numérica pode excluir escalas com denominador inválido em uma variante e, nesse caso, alterar os pesos normalizados efetivos. As ocorrências foram contadas e são apresentadas abaixo. Não há exclusão geométrica de bordas ou restrição por sobreposição. A operação resume a alternative WCO já calculada; não é apresentada como uma equação de Schuck nem como uma nova definição de coerência.

O critério de 5% é uma escolha de análise, não um limiar de significância nem uma máscara de borda. Foram conferidos também 1% e 10%, além da alternative WCO em zero e do máximo simples no plano completo. A proteção contra denominadores numericamente nulos foi preservada. Usar o máximo por escala permite atrasos diferentes em escalas diferentes; portanto o escore não comprova que exista um único atraso físico comum. As duas faixas foram comparadas no mesmo domínio integral; diferenças de amplitude e distribuição de potência causadas pela extração também precisam ser consideradas ao interpretar o resumo.

Na média balanceada, cada uma das quatro classes com pares intrapadrão tem peso igual; para controles, cada uma das dez combinações entre classes tem peso igual. Isso reduz a influência da quantidade de pares de churn e wavy annular, mas aumenta o peso relativo dos três pares de plug e dos três de slug frente à média dos 87 pares. Por isso ambas as agregações são explicitadas.

## Resultados da triagem — 64 escalas

'''+table(['Opção','Escore intrapadrão balanceado','Escore entre padrões','Variação intra vs referência corrigida','Variação do contraste intra−entre'],screen_rows)+f'''

A regra da triagem completa selecionou **{screen_choice}**. Essa escolha é exploratória e não redefiniu silenciosamente o conjunto de três opções da comparação com 200 escalas. A tabela apresenta também as variações nos controles entre padrões: elevar a média intrapadrão não basta se o mesmo efeito elevar ainda mais as comparações entre padrões.

## Implementação v5 na entrada 18–71 — 200 escalas

'''+table(['Opção','Escore intra','Escore entre','Intra−entre','Média dos máximos integrais intra','Média dos máximos integrais entre'],rows)+f'''

O escore intra balanceado do v5 foi **{main.loc['v5','same']:.6f}**, contra **{main.loc['central_band20','same']:.6f}** na faixa fixa. Na média simples dos 87 pares, os valores foram **{poolbase:.6f} → {poolcandidate:.6f}**, diferença de **{poolcandidate-poolbase:+.6f}**. Melhoraram **{improved}** pares e pioraram **{worse}**. O contraste entre a média intrapadrão e a média dos controles entre padrões foi **{main.loc['central_band20','gap']:+.6f}** para a faixa fixa; valores negativos indicam controles com média maior.

Os máximos simples e o escore ponderado são resumos diferentes. A busca no plano completo torna o máximo mais exposto a pontos extremos, inclusive às violações numéricas do v5. Os resultados de todos os cinco resumos estão em [preprocessing_comparison_rankings.csv](results/preprocessing_comparison_rankings.csv). Escolher apenas o resumo mais favorável não seria uma comparação suficiente.

O mesmo domínio de comparação significa os mesmos lags e os mesmos comprimentos de onda nominais. O comparador literal do v5 manteve sua conversão antiga de escalas e o preparo antes do corte, diferenças metodológicas documentadas em [EQUATION_AUDIT.md](EQUATION_AUDIT.md) e [REPORT.md](REPORT.md). A comparação com o v5 não isola apenas o efeito do pré-processamento.

### Peso igual para cada um dos 87 pares

'''+table(['Opção','Média dos 87 pares','Diferença vs v5','Pares com aumento','Pares com redução','Média balanceada por padrão'],pooledrows)+f'''

A mediana entre linhas apresentou escore médio **{poolmedian:.6f}**, com **{median_improved}** pares acima do v5 e **{median_worsened}** abaixo. A média por par e a média balanceada por padrão respondem a ponderações distintas e foram mantidas juntas para evitar atribuir aos 87 pares um benefício concentrado em classes pequenas.

Frente à referência corrigida do v16, a mediana entre linhas apresentou diferença de **{poolmedian-poolcorrected:+.6f}** na média dos 87 pares. A faixa 30–49 apresentou diferença de **{poolcandidate-poolcorrected:+.6f}**. Essas diferenças permitem separar o efeito dos controles do efeito da implementação numérica.

### Ganhos e perdas por padrão — faixa fixa contra v5

'''+table(['Padrão','Pares','Média v5','Média faixa fixa','Diferença'],grouprows)+f'''

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

{band_conclusion}

{native_conclusion}

{contrast_note} Os pares compartilham aquisições, e os resultados permanecem uma avaliação interna. Manter a formulação de Schuck não cria correspondência entre eventos de registros independentes.

Dados: preprocessing_comparison_pairs.csv, preprocessing_comparison_rankings.csv, preprocessing_pattern_tradeoffs.csv e preprocessing_heldout_vs_v5.csv. Figura: [comparação dos candidatos](figures/preprocessing/preprocessing_comparison.png). Mapas completos do candidato para o primeiro par de cada classe estão em figures/preprocessing, escolhidos por ordem alfabética, não pelo maior ganho.
'''
    text=text.replace('## Protocolo',band_text+'## Protocolo')
    text=text.replace('## Resultados da triagem — 64 escalas',denominator_text+'## Resultados da triagem — 64 escalas')
    (ROOT/'PREPROCESSING.md').write_text(text)
    executive=f'''# Conclusão executiva — v16_spatial2_18_71

{band_conclusion}

{native_conclusion}

## Entrada e domínio recalculados

Todos os gráficos foram refeitos com os nomes dos padrões em inglês e os demais textos em PT-BR, sem COI nem máscara conservadora. A média apresentada da alternative WPS foi calculada sobre todos os lags. Os escores de pré-processamento e o diagnóstico com fases aleatórias foram recalculados no plano completo; não houve alteração das equações de alternative WPS, alternative WCS e alternative WCO.

A análise anterior com entrada 20–59 está preservada em [v16_spatial2](../v16_spatial2/EXECUTIVE_SUMMARY.md), também sem máscaras. Portanto, a comparação entre faixas não envolve uma retirada adicional de restrições. Nesta pasta, a implementação v5 e os candidatos foram comparados nos mesmos lags e comprimentos de onda nominais, usando a nova entrada; o v5 manteve sua conversão antiga de escalas e seu preparo antes do corte, diferenças já auditadas.

## O que foi efetivamente verificado

Foram mantidas as equações de Schuck conferidas visualmente na análise de origem. A auditoria desta pasta validou os hashes dos 31 panoramas existentes e das planilhas originais de 20–59; herdou a verificação anterior de reconstrução pixel a pixel, sem refazê-la. O novo sinal foi extraído diretamente das 54 linhas 18–71 dos PNGs. As planilhas originais não foram tratadas como dados da nova faixa. A análise experimental usa os 26 registros confirmados, com 87 pares intrapadrão e 238 controles entre padrões. As equações alternativas permaneceram em uso durante todos os experimentos.

A faixa 18–71 amplia a cobertura vertical pretendida do tubo. Seus limites estão ilustrados na [auditoria visual da extração](results/row_band_visual_audit.png), mas não representam uma segmentação fisicamente validada de todo o interior do tubo. Os panoramas continuam usando os deslocamentos estimados originalmente com as linhas 20–59; a nova faixa altera a extração do sinal, sem reconstruir as imagens. O experimento verifica o efeito dessa faixa declarada sobre as análises wavelet.

A ampliação acrescentou **2 linhas acima** (18–19) e **12 abaixo** (60–71). Antes de remover a tendência e normalizar, a média das 54 linhas é exatamente **(40/54) × média de 20–59 + (14/54) × média das linhas acrescentadas**. Portanto, a região anterior contribui com aproximadamente 74,1% e a região adicionada com 25,9% da intensidade bruta. A motivação é incluir mais da região interna visível, especialmente as partes inferiores das estruturas nos exemplos plug e slug; essa escolha não valida a posição exata das paredes nem garante uma melhora da alternative WCO.

Com a remoção da tendência linear e a padronização utilizadas, uma contribuição média das linhas acrescentadas que fosse constante ou afim na posição produziria o mesmo sinal preparado, salvo arredondamento numérico. Portanto, eventuais perdas não podem ser atribuídas simplesmente à diluição por uma intensidade estática de parede. As diferenças decorrem de variação espacial acrescentada, que pode envolver estruturas e efeitos ópticos; sua origem física não foi identificada nesta análise.

As definições de **alternative WPS (32.14–32.15), alternative WCS (32.16) e alternative WCO (32.17)** foram preservadas. O estimador enviesado com 1/N continua sendo uma adaptação explícita: o artigo usa 1/(N−|m|) em 32.21–32.23. A versão literal não enviesada também foi calculada e apresentou instabilidade em grandes lags. O denominador da alternative WCO permanece no lag zero.

O v16 corrigiu a discretização da Morlet, a conversão de escala, a orientação do lag e o preparo após o corte comum. A versão corrigida respeitou os limites da alternative WCO nos testes e nos pares avaliados em todo o plano, sem limitar artificialmente os resultados a 1. **Essa correção matemática e numérica não comprova maior similaridade física entre registros.**

## Comparação entre implementações na entrada 18–71

Dez opções de pré-processamento foram recalculadas em 64 escalas para a entrada 18–71, conservando as fórmulas. A referência 18–71, a faixa fixa 30–49 e a mediana entre as linhas 18–71 foram conferidas em 200 escalas. A linha 40 e a faixa 30–49 foram mantidas em suas posições anteriores para isolar o efeito da ampliação da entrada principal. A faixa 30–49 permaneceu como comparação definida no estudo anterior; não foi escolhida novamente a partir dos resultados sem máscaras.

'''+table(['Opção','Média dos 87 pares','Diferença vs v5','Pares com aumento','Pares com redução','Média balanceada por padrão'],pooledrows)+f'''

O contraste intrapadrão menos entre padrões foi **{main.loc['central_band20','gap']:+.6f}** para a faixa fixa e **{main.loc['core_row_median','gap']:+.6f}** para a mediana entre linhas. {contrast_note}

Os resultados também dependeram do resumo do mapa: por exemplo, a média balanceada dos máximos integrais passou de **{mx.loc['v5','same']:.3f}** no v5 para **{mx.loc['central_band20','same']:.3f}** na faixa fixa, diferença de **{mx.loc['central_band20','same']-mx.loc['v5','same']:+.6f}**. Esse máximo é sensível a pontos extremos, inclusive às violações numéricas do v5. Não deve ser confundido com o escore ponderado da tabela nem interpretado isoladamente como medida da qualidade do método.

A mediana entre linhas apresentou diferença de **{poolmedian-poolcorrected:+.6f}** frente à referência corrigida na média por par. A média balanceada por padrão está apresentada separadamente. Ponderações diferentes podem mudar a ordenação das opções e não comprovam superioridade universal.

O escore é um resumo dos mapas, **não uma nova fórmula de WCO**: máximos por escala em todos os lags, ponderados pelas alternative WPS da referência no lag zero. O piso principal de 5% do pico dessas alternative WPS foi mantido, com verificações a 1% e 10%; não é um COI nem uma máscara geométrica. Os pesos brutos são fixos por par; a proteção de denominadores pode alterar sua normalização efetiva se uma variante tiver escalas numericamente inválidas. {denominator_note} As contagens estão em [PREPROCESSING.md](PREPROCESSING.md). A triagem selecionou **{screen_choice}** pela regra declarada. A avaliação retirando uma aquisição usou exclusivamente as três opções fixadas e 200 escalas. As classes possuem poucas aquisições e os pares não são réplicas independentes.

## Comprimento dos registros e janelas

**As janelas adicionais não são exigidas pelas equações de Schuck.** A análise principal usou o registro completo disponível, ou o prefixo comum de um par. O controle de segmentos de início, centro e fim permaneceu como diagnóstico separado de sensibilidade ao comprimento e à posição. A retirada do COI e da máscara conservadora não elimina possíveis efeitos de borda nem a menor quantidade de produtos que compõem a correlação em lags extremos.

## Figuras dos padrões e dos resultados

Os exemplos abaixo usam a entrada principal, média das linhas 18–71, e foram escolhidos pela ordem alfabética dos registros de cada padrão. São exemplos do conjunto, não os casos com maior coerência.

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

As figuras 7–10 mostram o primeiro par em ordem alfabética de cada padrão com pelo menos duas aquisições elegíveis. Usam a implementação corrigida do v16, a correlação enviesada e a média das linhas 18–71. No topo estão as alternative WPS individuais e suas médias sobre todos os lags; no centro, o módulo da alternative WCS; abaixo, a alternative WCO da Eq. 32.17, com denominador no lag zero. Os mapas exibem o plano completo, sem COI nem máscara conservadora.

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

![Alternative WPS, alternative WCS e alternative WCO entre os registros wavy annular m03q10_001 e m03q20_001](figures/pairs/m03q10_001__m03q20_001.png)

*Figura 10 — Comparação de duas aquisições de wavy annular. Compartilhar o rótulo do padrão não exige que a alternative WCO seja alta em todas as escalas ou defasagens.*

### Ganhos e perdas frente ao v5

![Comparação dos pré-processamentos e dos escores de alternative WCO dos pares frente ao v5](figures/preprocessing/preprocessing_comparison.png)

*Figura 11 — À esquerda, a triagem com 64 escalas compara as opções com a referência corrigida do v16. Ao centro, a faixa 30–49 é comparada ao v5 em 200 escalas: {improved} pares acima da diagonal e {worse} abaixo. À direita, a mediana entre linhas apresenta {median_improved} aumentos e {median_worsened} reduções. Os eixos dos dois diagramas de dispersão têm os mesmos limites. Os escores usam todos os lags, sem COI nem máscara conservadora; não são novas fórmulas de coerência.*

## Decisão honesta

{band_conclusion}

{contrast_note} A escolha da faixa deve considerar simultaneamente a comparação com a extração anterior, as perdas individuais, os controles entre padrões e a integridade do observável. O resultado de uma única ponderação não estabelece uma opção universalmente superior.

O principal limite permanece: pertencer ao mesmo padrão não significa que duas aquisições independentes contenham os mesmos eventos ou a mesma organização espacial. Não seria correto forçar a coerência a subir por seleção de segmentos orientada pela própria coerência, deformação do eixo, redução artificial da banda ou mudança do denominador.

Os próximos caminhos continuam documentados apenas como planos: [v16.1](../v16.1/PLAN.md), comparação de perfis da alternative WPS, especialmente no lag zero; e [v16.2](../v16.2/PLAN.md), comparação de eventos equivalentes definidos por critérios de imagem antes de observar a coerência. Nenhuma dessas novas análises foi executada nesta atualização. A segunda responderia sobre eventos selecionados, não sobre a similaridade dos panoramas completos.

Detalhes: [equações conferidas](EQUATION_AUDIT.md), [experimento e comparação direta com v5](PREPROCESSING.md), [relatório metodológico](REPORT.md), [resultados gerais](RESULTS.md) e [galeria](figures/index.html).
'''
    executive=executive.replace('## Entrada e domínio recalculados',band_text+'## Entrada e domínio recalculados')
    executive=executive.replace('## Comprimento dos registros e janelas',null_text+'## Comprimento dos registros e janelas')
    (ROOT/'EXECUTIVE_SUMMARY.md').write_text(executive)
    summary=dict(conclusion=band_conclusion+' '+native_conclusion,
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
    summary['band_comparison']=band_records
    summary['denominator_diagnostics']=denominator_records
    (ROOT/'results/preprocessing_conclusion.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
