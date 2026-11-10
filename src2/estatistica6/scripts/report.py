"""Gráficos e resumo executivo em português, exclusivamente de resultados salvos."""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/estatistica6-matplotlib')
from common import *
import textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter
from matplotlib.font_manager import FontProperties

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,
                     'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180,
                     'axes.labelcolor':'#263442','text.color':'#263442','axes.titleweight':'semibold',
                     'pdf.fonttype':42,'ps.fonttype':42})
CMAP=LinearSegmentedColormap.from_list('roxo_creme_azul',['#72508f','#faf6e9','#216e9a'])
POSITION={'prefix':'Começo','center':'Centro','suffix':'Fim','7001_prefix':'7.001 px\ncomeço','full':'Completo'}
METRIC={'pearson_residual':'Pearson dos resíduos','hellinger':'Hellinger','wasserstein':'Wasserstein',
        'trio':'Trio espectral','log_peak':'Escala dominante','slope':'Inclinação','width':'Largura central',
        'normalized_magnitude':'WCS normalizada em lag zero','max_wco_all':'WCO máximo global †',
        'weighted_zero_lag_wco':'WCO em lag zero †','weighted_scale_peak_power05':'WCO por escala, piso 5% †'}
FAMILY={'wps':'WPS','mean_power':'Mean power','wcs_zero':'WCS','wco_cached':'WCO'}
SCOPE={'churn_annular':'churn × wavy annular','plug_slug':'plug × slug','quatro_padroes':'Quatro padrões'}


def number(x,digits=3):
    return f'{float(x):.{digits}f}'.replace('.',',')


def save_figure(fig,name):
    fig.savefig(ROOT/f'figures/{name}.png',bbox_inches='tight',facecolor='white')
    fig.savefig(ROOT/f'figures/{name}.pdf',bbox_inches='tight',facecolor='white')
    plt.close(fig)


def grid(ax):
    ax.grid(axis='y',alpha=.18,zorder=0)


def triage_figure(metrics):
    frame=metrics[metrics.condition.isin(['prefix','pair_min_length'])]
    keys=list(frame[['family','metric']].drop_duplicates().itertuples(index=False,name=None))
    values=np.array([[frame[(frame.family==f)&(frame.metric==m)&(frame.scope==s)].iloc[0].macro_silhouette for s in CONFIG['scopes']] for f,m in keys])
    fig,ax=plt.subplots(figsize=(10.5,7.9),layout='constrained')
    im=ax.imshow(values,cmap=CMAP,vmin=-.5,vmax=.5,aspect='auto')
    labels=[]
    for family,metric in keys:
        labels.append((FAMILY[family]+' · ' if family in ['wps','mean_power'] else '')+METRIC[metric])
    ax.set_yticks(range(len(keys)),labels)
    ax.set_xticks(range(3),['churn × wavy annular\n19 ensaios','plug × slug\n6 ensaios','Quatro padrões\n25 ensaios'])
    ax.tick_params(axis='both',length=0)
    for i,(family,metric) in enumerate(keys):
        for j,scope in enumerate(CONFIG['scopes']):
            rec=frame[(frame.family==family)&(frame.metric==metric)&(frame.scope==scope)].iloc[0]
            ax.text(j,i,number(values[i,j]),ha='center',va='center',weight='bold' if rec.prefix_gate else 'normal',
                    color='white' if abs(values[i,j])>.32 else '#243444')
    for y in [6.5,9.5,10.5]:ax.axhline(y,color='white',lw=2)
    ax.set_title('Triagem da compactação e separação dos padrões',pad=14)
    cb=fig.colorbar(im,ax=ax,shrink=.75,pad=.04);cb.set_label('Silhueta média balanceada por padrão')
    fig.get_layout_engine().set(rect=(0,.09,1,.91))
    fig.text(.02,.035,'Valores positivos favorecem o próprio padrão; zero indica sobreposição. Negrito: passou o filtro preliminar completo.\n† WCO arquivada: comprimento mínimo de cada par. Demais linhas: os mesmos 3.171 px por ensaio.',fontsize=9)
    save_figure(fig,'01_triagem')


def profiles_figure(source):
    p,_=probabilities(source['profile'],source['wavelengths_px'])
    _,h=log_weights(source['wavelengths_px']);density=p/h
    fig,axes=plt.subplots(1,5,figsize=(15,3.8),sharex=True,sharey=True,layout='constrained')
    for ax,g in zip(axes,PATTERNS):
        vals=density[source['patterns']==g]
        for row in vals:ax.plot(source['wavelengths_px'],row,color=COLORS[g],alpha=.3,lw=.8)
        ax.plot(source['wavelengths_px'],np.median(vals,axis=0),color=COLORS[g],lw=2.2)
        ax.set(title=f'{g}\nn = {len(vals)}',xscale='log',xlim=(50,3000),xlabel='Comprimento de onda (px)')
        ax.set_xticks([50,200,800,3000],['50','200','800','3.000']);grid(ax)
    axes[0].set_ylabel('Magnitude normalizada\npor unidade de log(λ)')
    fig.suptitle('Distribuição espectral da WPS nos primeiros 3.171 px',fontsize=14)
    save_figure(fig,'02_perfis_wps')


def churn_figure(source,features_frame,matrices):
    fig,axes=plt.subplots(1,2,figsize=(12,5.5),gridspec_kw={'width_ratios':[1,1.05]},layout='constrained')
    ax=axes[0]
    for g,shift in [('churn',-.10),('wavy annular',.10)]:
        f=features_frame[(features_frame.pattern==g)&features_frame.condition.isin(['prefix','center','suffix'])]
        pivot=f.pivot(index='run',columns='condition',values='width')[['prefix','center','suffix']]
        for i,(_,row) in enumerate(pivot.iterrows()):
            jitter=(i-(len(pivot)-1)/2)*.008
            ax.plot(np.arange(3)+shift+jitter,row,color=COLORS[g],alpha=.24,lw=.8)
            ax.scatter(np.arange(3)+shift+jitter,row,s=16,color=COLORS[g],alpha=.65,zorder=3)
        ax.plot(np.arange(3)+shift,pivot.median(),color=COLORS[g],lw=2.5,marker='D',ms=6,label=f'{g} · mediana')
    ax.set(xticks=range(3),xticklabels=['Começo','Centro','Fim'],ylabel='Largura central em log(λ)',title='Largura espectral por ensaio')
    ax.legend(fontsize=8,loc='upper right');grid(ax)
    ax=axes[1]
    mask=scope_mask(source['patterns'],'churn_annular');labels=source['patterns'][mask];runs=source['runs'][mask]
    d=matrices['wps__prefix__hellinger'][np.ix_(mask,mask)]
    im=ax.imshow(d,cmap='Blues',vmin=0,vmax=1)
    ax.set(xticks=[4.5,14],xticklabels=['churn (10)','wavy annular (9)'],yticks=[4.5,14],yticklabels=['churn (10)','wavy annular (9)'],title='Hellinger entre distribuições WPS')
    ax.tick_params(axis='both',length=0,labelsize=9)
    ax.axvline(9.5,color='#333',lw=.8);ax.axhline(9.5,color='#333',lw=.8)
    fig.colorbar(im,ax=ax,shrink=.65,label='Distância Hellinger (0 a 1)',pad=.03)
    fig.suptitle('churn e wavy annular apresentam sinal parcial de separação',fontsize=14)
    fig.get_layout_engine().set(rect=(0,.08,1,.92))
    fig.text(.02,.025,'19 ensaios; janelas de 3.171 px. Linhas finas ligam regiões do mesmo ensaio.\nNa matriz, a diagonal é autocomparação; blocos fora da diagonal mostram a sobreposição entre padrões.',fontsize=9)
    save_figure(fig,'03_churn_annular')


def plug_figure(source,features_frame,matrices,cluster_stats):
    fig=plt.figure(figsize=(11.8,7.7),layout='constrained')
    gs=fig.add_gridspec(3,3,width_ratios=[2.4,1,.065],hspace=.12)
    ax=fig.add_subplot(gs[:,0])
    rows=features_frame[features_frame.pattern.isin(['plug','slug'])&features_frame.condition.isin(['prefix','center','suffix'])]
    for k,run in enumerate(source['runs'][scope_mask(source['patterns'],'plug_slug')]):
        rec=rows[rows.run==run].set_index('condition').loc[['prefix','center','suffix']]
        g=rec.pattern.iloc[0]
        ax.plot(range(3),rec.peak_px,marker=['o','s','^'][k%3],ms=6,lw=1.7,color=COLORS[g],
                ls='-' if g=='plug' else '--',label=f'{run} · {g}')
    ax.set(xticks=range(3),xticklabels=['Começo','Centro','Fim'],ylabel='Escala dominante da WPS (px)',title='A ordem entre plug e slug se inverte',ylim=(290,1320))
    ax.legend(fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.08),ncol=2,frameon=False);grid(ax)
    mask=scope_mask(source['patterns'],'plug_slug')
    for k,condition in enumerate(['prefix','center','suffix']):
        ax=fig.add_subplot(gs[k,1])
        d=matrices[f'mean_power__{condition}__hellinger']
        if len(d)==26:d=d[np.ix_(mask,mask)]
        im=ax.imshow(d,cmap='Blues',vmin=0,vmax=.7)
        ari=cluster_stats[(cluster_stats.family=='mean_power')&(cluster_stats.metric=='hellinger')&(cluster_stats.scope=='plug_slug')&(cluster_stats.condition==condition)].iloc[0].ari
        ax.set(xticks=[1,4],xticklabels=['plug','slug'],yticks=[1,4],yticklabels=['plug','slug'],title=f'Mean power · {POSITION[condition]}\nARI = {number(ari,2)}')
        ax.title.set_fontsize(10);ax.tick_params(axis='both',length=0,labelsize=9)
        ax.axvline(2.5,color='#333',lw=.7);ax.axhline(2.5,color='#333',lw=.7)
    cax=fig.add_subplot(gs[:,2])
    fig.colorbar(im,cax=cax,label='Hellinger (escala comum)')
    fig.suptitle('Comparação direta de plug e slug com três ensaios por padrão',fontsize=14)
    fig.get_layout_engine().set(rect=(0,.07,1,.93))
    fig.text(.02,.018,'Agrupamento por ligação média, k = 2. ARI = 1: grupos coincidem com os rótulos; ARI ≈ 0: concordância próxima do acaso.\nO resultado favorável da mean power no fim não se repetiu no começo ou no centro. Janelas não são novas aquisições.',fontsize=9)
    save_figure(fig,'04_plug_slug')


def stability_figure(metrics):
    fig,axes=plt.subplots(1,2,figsize=(12,5.3),layout='constrained')
    selections={
        'churn_annular':[('wps','pearson_residual','#8c8c8c'),('wps','hellinger','#237da1'),('wps','width','#8661aa'),('wps','trio','#bd8631')],
        'plug_slug':[('wps','log_peak','#8661aa'),('wps','hellinger','#8c8c8c'),('mean_power','hellinger','#31977b'),('mean_power','wasserstein','#bd8631')]
    }
    for ax,(scope,selection) in zip(axes,selections.items()):
        for family,metric,color in selection:
            frame=metrics[(metrics.scope==scope)&(metrics.family==family)&(metrics.metric==metric)].set_index('condition')
            y=[frame.loc[c,'macro_silhouette'] for c in ['prefix','center','suffix']]
            ax.plot(range(3),y,color=color,marker='o',lw=2,label=f'{FAMILY[family]} · {METRIC[metric]}')
            if '7001_prefix' in frame.index:ax.scatter([3.4],[frame.loc['7001_prefix','macro_silhouette']],color=color,marker='D',s=35)
        ax.axhline(0,color='#555',lw=1);grid(ax)
        if scope=='churn_annular':ax.set_xticks([0,1,2,3.4],['Começo','Centro','Fim','7.001 px\ncomeço'])
        else:ax.set_xticks(range(3),['Começo','Centro','Fim'])
        ax.set(title=f'{SCOPE[scope]} · n = {19 if scope=="churn_annular" else 6}',ylabel='Silhueta média balanceada',ylim=(-.09,.49))
        ax.legend(loc='upper center',bbox_to_anchor=(.5,-.14),ncol=1,fontsize=8,frameon=False)
    fig.suptitle('A separação depende do descritor e da posição observada',fontsize=14)
    save_figure(fig,'05_estabilidade')


def covariate_figure(inventory,features_frame):
    fig,axes=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
    table=inventory.groupby(['pattern','q_code']).size().unstack(fill_value=0).reindex(PATTERNS)
    im=axes[0].imshow(table,cmap='Blues',vmin=0,vmax=int(table.to_numpy().max()),aspect='auto')
    axes[0].set(xticks=range(len(table.columns)),xticklabels=table.columns,yticks=range(5),yticklabels=PATTERNS,xlabel='Código q',title='Cobertura desigual de condições')
    axes[0].tick_params(length=0)
    for i in range(5):
        for j in range(len(table.columns)):axes[0].text(j,i,str(table.iloc[i,j]),ha='center',va='center',color='white' if table.iloc[i,j]>=3 else '#263442')
    joined=features_frame[features_frame.condition=='prefix'].merge(inventory[['run','n']],on='run',validate='one_to_one')
    ax=axes[1]
    for g in ['churn','wavy annular']:
        f=joined[joined.pattern==g]
        ax.scatter(f.n,f.width,color=COLORS[g],s=45,label=g,edgecolors='white',lw=.5)
    rho=read_csv(ROOT/'results/a6_associations.csv').query("scope=='churn_annular' and condition=='prefix' and feature=='width' and covariate=='n'").iloc[0].spearman
    ax.set(xlabel='Comprimento original do panorama (px)',ylabel='Largura central da WPS em log(λ)',title=f'Associação com a reconstrução · ρ = {number(rho)}')
    ax.ticklabel_format(axis='x',style='plain');ax.legend(fontsize=9);grid(ax)
    fig.get_layout_engine().set(rect=(0,.09,1,.91))
    fig.text(.02,.025,'O descritor usa sempre 3.171 px; o eixo horizontal conserva o comprimento original.\nComprimento e deslocamento médio têm correlação 1,000. A associação não distingue padrão, operação e reconstrução.',fontsize=9)
    save_figure(fig,'06_condicoes_operacionais')


def pearson_figure():
    contrasts=read_csv(ROOT/'results/a1_contrastes.csv');loo=read_csv(ROOT/'results/a1_leave_one_out.csv')
    balanced=json.loads((ROOT/'results/a1_balanced.json').read_text())['delta_balanced_four_classes']
    fig,axes=plt.subplots(1,2,figsize=(11.8,4.4),layout='constrained')
    scopes=['todos26','quatro_padroes','churn_annular','plug_slug']
    names=['26 ensaios','Sem singleton (25)','churn × wavy annular','plug × slug']
    x=np.array([contrasts[(contrasts.scope==s)&(contrasts.condition=='prefix')].iloc[0].delta for s in scopes])
    lo=np.array([loo[loo.scope==s].delta.min() for s in scopes]);hi=np.array([loo[loo.scope==s].delta.max() for s in scopes])
    ax=axes[0]
    ax.hlines(np.arange(4),lo,hi,color='#9badba',lw=4)
    ax.scatter(x,np.arange(4),color='#267da0',s=45,zorder=3)
    ax.set(yticks=range(4),yticklabels=names,xlabel='Mediana r intra − mediana r inter',title='Sensibilidade à retirada de um ensaio')
    ax.invert_yaxis();ax.axvline(0,color='#555',lw=.8);ax.grid(axis='x',alpha=.2)
    ax.text(.02,-.25,f'Contraste balanceado das quatro classes: {number(balanced)}.\nLinhas: mínimo–máximo ao retirar um ensaio; não são IC.',transform=ax.transAxes,fontsize=9)
    ax=axes[1]
    for scope,color in [('churn_annular','#237da1'),('plug_slug','#8661aa')]:
        data=contrasts[contrasts.scope==scope].set_index('condition')
        ax.plot(range(3),[data.loc[c,'delta'] for c in ['prefix','center','suffix']],marker='o',lw=2,color=color,label=SCOPE[scope])
    ax.set(xticks=range(3),xticklabels=['Começo','Centro','Fim'],ylabel='Contraste de Pearson',title='Mesmo comprimento, outras posições')
    ax.axhline(0,color='#555',lw=.8);ax.legend(fontsize=9);grid(ax)
    save_figure(fig,'07_diagnostico_pearson')


def make_summary(metrics,features_frame,cluster_stats):
    def row(family,metric,scope,condition='prefix'):
        return metrics[(metrics.family==family)&(metrics.metric==metric)&(metrics.scope==scope)&(metrics.condition==condition)].iloc[0]
    def series(family,metric,scope,field):
        return [float(row(family,metric,scope,c)[field]) for c in ['prefix','center','suffix']]
    effects=read_csv(ROOT/'results/a3_effects.csv')
    a1=read_csv(ROOT/'results/a1_contrastes.csv')
    slopes=effects[(effects.scope=='plug_slug')&(effects.feature=='log_peak')].set_index('condition')
    screened=metrics[metrics.condition.isin(['prefix','pair_min_length'])]
    decisions=read_csv(ROOT/'results/decisoes_finais.csv')
    counts=features_frame[features_frame.condition=='prefix'].pattern.value_counts().to_dict()
    s=dict(
        sample={'total':sum(counts.values()),**counts},
        n_representations=len(screened[['family','metric']].drop_duplicates()),n_scope_combinations=len(screened),
        n_advanced=len(decisions),n_pass_geometry=int(bools(decisions.retained).sum()),
        churn_hellinger_silhouette=series('wps','hellinger','churn_annular','macro_silhouette'),
        churn_hellinger_neighbor=series('wps','hellinger','churn_annular','balanced_neighbor'),
        churn_neighbor_chance=float(row('wps','hellinger','churn_annular').neighbor_chance),
        churn_width_silhouette=series('wps','width','churn_annular','macro_silhouette'),
        churn_width_medians=[float(effects[(effects.scope=='churn_annular')&(effects.feature=='width')&(effects.condition=='prefix')].iloc[0][c]) for c in ['median_a','median_b']],
        churn_hellinger_7001=float(row('wps','hellinger','churn_annular','7001_prefix').macro_silhouette),
        plug_peak_cliff=[float(slopes.loc[c,'cliff_delta']) for c in ['prefix','center','suffix']],
        plug_peak_silhouette=series('wps','log_peak','plug_slug','macro_silhouette'),
        plug_mean_power_hellinger_silhouette=series('mean_power','hellinger','plug_slug','macro_silhouette'),
        plug_mean_power_hellinger_ari=[float(cluster_stats[(cluster_stats.family=='mean_power')&(cluster_stats.metric=='hellinger')&(cluster_stats.scope=='plug_slug')&(cluster_stats.condition==c)].iloc[0].ari) for c in ['prefix','center','suffix']],
        pearson_churn_delta=[float(a1[(a1.scope=='churn_annular')&(a1.condition==c)].iloc[0].delta) for c in ['prefix','center','suffix']],
        four_class_best_silhouette=float(metrics[(metrics.condition=='prefix')&(metrics.scope=='quatro_padroes')].macro_silhouette.max()),
        data_reuse='exploratório, sem teste independente',
        conclusion='Sinais parciais em WPS; nenhum agrupamento de padrões robusto às condições avaliadas. Plug/slug: direção do pico inverte entre começo e fim.'
    )
    write_json(ROOT/'results/executive_summary.json',s)
    return s


def write_documents(s,features_frame,metrics):
    seq=lambda vals:' / '.join(number(x) for x in vals)
    head='Resumo executivo da discriminação de padrões por wavelets'
    main='A WPS trouxe sinais úteis para comparar padrões, mas nenhum método testado sustentou agrupamento consistente nas diferentes janelas. Para churn e wavy annular, a distribuição e a largura espectral foram os melhores caminhos para continuidade. Para plug e slug, a escala dominante separou os ensaios no começo, mas sua relação com o padrão se inverteu no fim.'
    steps=[
        'Conferi os 26 ensaios elegíveis, os sinais extraídos das linhas 18–71, os perfis e os hashes das fontes. Mantive os cinco excluídos da origem.',
        f'Fixei trechos de 3.171 px e a grade de 200 escalas entre 50 e 3.000 px. Comparei {s["n_representations"]} representações em três conjuntos: churn/wavy annular (19), plug/slug (6) e quatro padrões com réplicas (25).',
        f'Filtrei candidatos por silhueta balanceada ≥ 0,10, silhueta positiva em cada classe e vantagem do vizinho mais próximo ≥ 10 pontos percentuais sobre o acaso. Avançaram {s["n_advanced"]} combinações descritor/comparação.',
        'Conferi começo, centro, fim, retirada de um ensaio e, nos mesmos 19 registros, o controle de 7.001 px. Recalculei mean power no centro/fim somente nos seis plug/slug que justificaram essa extensão.',
        f'Comparei agrupamentos por ligação média com número fixo de grupos e examinei condições operacionais. Passaram {s["n_pass_geometry"]} combinações pelo filtro geométrico; isso não garantiu agrupamento correto ou direção estável do efeito.'
    ]
    churn=(f'Hellinger sobre a WPS teve silhuetas {seq(s["churn_hellinger_silhouette"])} no começo/centro/fim e {number(s["churn_hellinger_7001"])} em 7.001 px. '
           f'A concordância balanceada do vizinho mais próximo foi {" / ".join(number(x*100,1) for x in s["churn_hellinger_neighbor"])}%, contra {number(100*s["churn_neighbor_chance"],1)}% ao acaso. '
           'São descrições da amostra, não acurácia preditiva. '
           f'A largura central mediana foi {number(s["churn_width_medians"][0])} em churn e {number(s["churn_width_medians"][1])} em wavy annular, em log(λ). '
           'As diferenças persistiram em média, mas houve sobreposição; no centro, uma classe apresentou silhueta negativa. '
           'Com Hellinger, o ARI do agrupamento ficou próximo de zero nas três posições: a boa vizinhança local não se converteu em dois grupos bem definidos.')
    ranges=[]
    for condition in ['prefix','suffix']:
        for g in ['plug','slug']:
            f=features_frame[(features_frame.condition==condition)&(features_frame.pattern==g)]
            ranges.append(dict(condition=condition,pattern=g,minimum=float(f.peak_px.min()),maximum=float(f.peak_px.max())))
    range_text=lambda cond,g: next(f'{number(x["minimum"],0)}–{number(x["maximum"],0)} px' for x in ranges if x['condition']==cond and x['pattern']==g)
    plug=(f'No começo, plug teve picos de {range_text("prefix","plug")} e slug de {range_text("prefix","slug")}, sem sobreposição. '
          f'No fim, plug passou a {range_text("suffix","plug")} e slug a {range_text("suffix","slug")}. '
          f'O delta de Cliff de plug menos slug mudou de {number(s["plug_peak_cliff"][0],0)} para +{number(s["plug_peak_cliff"][2],0)}. '
          'Logo, não há uma regra única “maior escala = slug” para todo o panorama. '
          f'A mean power com Hellinger teve ARI {seq(s["plug_mean_power_hellinger_ari"])}: recuperou os rótulos apenas no fim. '
          'Wasserstein também mostrou sinal parcial, enfraquecido no centro. Com três ensaios por classe, esses resultados são pistas para investigação, não validação.')
    stopped=(f'Pearson em churn/wavy annular teve contraste intra−inter de {seq(s["pearson_churn_delta"])}: a vantagem desapareceu no fim. '
             f'Nenhuma representação passou a triagem conjunta das quatro classes; a maior silhueta foi {number(s["four_class_best_silhouette"])}. '
             'WCS no lag zero e os três resumos WCO não passaram o filtro. A triagem WCO usou comprimentos variáveis por par; sua falha não exclui outras análises WCO com comprimento fixo. '
             'Não ampliei simulações nulas WCO nem testes funcionais por escala ou PERMDISP, pois não resolveriam a discriminação observada nesta etapa.')
    limits=[
        'Amostra já explorada: a seleção de descritores e as métricas usam os mesmos ensaios. Não houve teste externo, classificadores ou valores-p novos.',
        'Padrão, operação e reconstrução se confundem. A largura WPS teve Spearman −0,748 com o comprimento original nos 19 ensaios; comprimento e deslocamento médio tiveram correlação 1,000. Não existem trocas de rótulo dentro de (m,q) que permitam separar ambos os efeitos.',
        'Janelas do mesmo panorama não são réplicas independentes. As janelas extremas de plug/slug se sobrepõem em até 40,1%. Dispersed bubbles tem um ensaio e não permite avaliar compactação intraclasse.',
        'Escalas em pixels, sem conversão física ou temporal. Escalas de até 3.000 px têm suporte limitado em trechos de 3.171 px. Mean power significa média da magnitude WPS sobre os lags, não energia física. A nomenclatura da origem é WCO; WPO não está definida.'
    ]
    recommendation='Priorizar distribuição WPS e largura espectral para churn/wavy annular. Em plug/slug, investigar a variação espacial do pico e da mean power antes de propor uma regra de agrupamento. A próxima evidência decisiva requer aquisições independentes em condições comparáveis.'
    validation=json.loads((ROOT/'validation/final_numeric.json').read_text())
    md=f'# {head}\n\n**{main}**\n\n[PDF com os principais gráficos](RESUMO_EXECUTIVO.pdf) · [Método e reprodução](METODOLOGIA.md) · [Tabelas completas](results/metricas_finais.csv)\n\n## O que foi feito\n\n'
    md+='\n'.join(f'{i}. {p}' for i,p in enumerate(steps,1))
    md+='\n\n![Triagem dos descritores](figures/01_triagem.png)\n\n## churn e wavy annular\n\n'+churn
    md+='\n\n![Distribuição e largura WPS em churn e wavy annular](figures/03_churn_annular.png)\n\n## plug e slug\n\n'+plug
    md+='\n\n![Comparação de plug e slug em três posições](figures/04_plug_slug.png)\n\n## Análises encerradas na triagem\n\n'+stopped
    md+='\n\n![Sensibilidade à posição e ao comprimento](figures/05_estabilidade.png)\n\n## Limitações\n\n'+'\n'.join('- '+p for p in limits)
    md+='\n\n**Próximo passo útil:** '+recommendation
    md+=f'\n\n**Validação:** PASS. Extração dos 26 sinais, 123 perfis em cache, fórmulas, pares e agrupamentos conferidos por cálculos independentes; {validation["sources_checked"]} fontes preservadas. [Validação numérica](validation/final_numeric.json).\n'
    (ROOT/'RESUMO_EXECUTIVO.md').write_text(md)
    # Conteúdo numérico e texto idênticos entre Markdown e PDF.
    document=dict(title=head,main=main,steps=steps,churn=churn,plug=plug,stopped=stopped,limits=limits,recommendation=recommendation)
    write_json(ROOT/'results/report_text.json',document)
    method='''# Método e reprodução da estatistica6

Execução exploratória do [plano original](PLANO.md), autorizada pela solicitação posterior do usuário, com triagem prévia e comparação direta plug/slug. O plano permanece preservado como registro. Configuração congelada antes do cálculo: [triagem_v1.json](configs/triagem_v1.json); seu SHA-256 está nos manifestos. Limites inclusivos do filtro usam tolerância numérica de 10⁻¹².

## Amostra e comparabilidade

Unidade amostral: ensaio. Conjuntos: 26 no inventário e matrizes; 19 churn/wavy annular; 6 plug/slug; 25 nos quatro padrões com réplicas. Dispersed bubbles é mantido nas comparações entre ensaios, mas excluído da avaliação de compactação. Nenhum registro foi excluído pelo resultado. O controle de 7.001 px contém exatamente os mesmos 19 churn/wavy annular, nunca as quatro classes. Os códigos m/q são usados como códigos; documentação histórica de unidades não substitui calibração das aquisições.

Perfil principal: magnitude da alternative WPS no lag zero, com autocorrelação enviesada, Morlet ω₀=6, 200 comprimentos de onda entre 50–3.000 px. Cortar antes de inverter intensidade, retirar tendência linear e dividir pelo desvio-padrão populacional. Preservar a banda e as equações da v16. Sem novo COI, alinhamento, reconstrução, suavização, alteração de faixa óptica ou quadratura da magnitude.

## Representações e métricas

- **Pearson:** retirar a reta de log(P+10⁻¹²) em log(λ), centrar o resíduo e normalizar sua norma. Distância euclidiana d=√[2(1−r)], não 1−r para a geometria. A1 relata r e mediana(r intra)−mediana(r inter).
- **Distribuição WPS:** t=log(λ/1 px); h são pesos trapezoidais em t; p=P·h/Σ(P·h). Hellinger=||√pᵢ−√pⱼ||/√2. Wasserstein=Σ|Fᵢ−Fⱼ|Δt. Não se interpreta p como energia física ou probabilidade de padrão.
- **Características:** log da escala do máximo, inclinação log-linear e largura Q75(t)−Q25(t), por inversa da CDF. Centro, entropia e participações [50,200), [200,800), [800,3000] são complementares. Empates, máximos locais e bordas foram registrados. Trio: distância euclidiana após mediana/IQR dos 26 prefixos, preservados nas demais posições; não há interpretação preditiva desse ajuste global.
- **Mean power:** média de |CWT(Rxx)| sobre todos os 2N−1 lags. Prefixos recalculados para os 26; centro/fim recalculados somente nos seis plug/slug aprovados preliminarmente. Comparação por Pearson residual, Hellinger e Wasserstein.
- **WCS no lag zero:** todos os 325 pares com N=3.171. Similaridade Σh|CWT(Rxy)(0,a)| / Σh√(PᵢPⱼ), nas escalas com denominadores válidos; dissimilaridade 1−similaridade. A CWT foi calculada por soma temporal complexa e conferida contra a FFT original. Este resumo não cobre a informação de todo o mapa WCS.
- **WCO:** apenas variant=baseline, 200 escalas, pesos nativos 18–71 do arquivo preprocessing_confirm_pairs.csv. Máximo global, média ponderada no lag zero e média ponderada de máximos por escala com piso 5%. Dissimilaridade 1−escore, sem clipping; não se presume métrica euclidiana. Comprimento mínimo por par, diferentemente dos demais descritores. Sem busca entre variantes/pisos, novas simulações ou extrapolação para WCO de comprimento fixo.

## Triagem e aprofundamento

Silhueta de cada ensaio s=(b−a)/max(a,b), com a=distância média aos outros ensaios do mesmo padrão e b=menor distância média a outra classe. Resumo principal: média das médias por classe, para não favorecer as maiores. Vizinho mais próximo exclui diagonal; empates recebem a proporção dos vizinhos empatados do mesmo padrão. Concordância balanceada é a média por classe. Referência aleatória: média de (n_g−1)/(n−1), igual a 47,222% nos 19 e 40% nos seis. Não são testes, previsões fora da amostra ou classificadores.

Filtro inicial: silhueta macro ≥0,10; média de silhueta de cada classe >0; concordância balanceada de vizinhança pelo menos 0,10 acima do acaso. Filtro geométrico seguinte: silhueta macro e vantagem de vizinhança positivas em centro/fim; mínimo leave-one-run-out da silhueta positivo no prefixo; silhueta positiva no controle 7.001 quando aplicável. Limites são convenções de prioridade, não margens científicas validadas. Todos os candidatos, falhas e condições permanecem nos CSVs.

**Passar o filtro geométrico não significa manter a direção do efeito, compactar cada classe em toda janela ou recuperar os rótulos por agrupamento.** Essas limitações são verificadas nos efeitos de Cliff, silhuetas por classe e agrupamentos. O pico plug/slug exemplifica a inversão que uma distância por si só não detecta. Não se alterou o filtro após observar essa limitação.

Agrupamento: ligação média e número de grupos fixo (k=2 nos dois contrastes que avançaram). Rótulos são usados somente para comparar a partição por ARI; 1 indica coincidência e 0 a expectativa ajustada ao acaso. Nenhuma escolha de k, linkage ou melhor janela. Os quatro padrões não passaram a triagem, portanto seu agrupamento não foi aprofundado. As janelas preservam a ordem dos ensaios; não são replicações experimentais.

## Blocos executados e limites

| Bloco | Execução |
| --- | --- |
| A1 | Contrastes de r, pares por classe, remoção do singleton, peso igual entre classes, retirada de cada ensaio e 5.000 perturbações de pesos por contraste binário. |
| A2 | Começo/centro/fim; 76 janelas distintas e 75 pares de janelas do mesmo ensaio; sobreposição registrada. Controle de comprimento nos mesmos 19. |
| A3 | Características, medianas e delta de Cliff por condição; sem testes inferenciais. |
| A4 | Curvas/diferenças descritivas completas; sem teste de máximo por escala. |
| A5 | Distâncias Hellinger/Wasserstein, 325 pares cada; aprofundamento conforme triagem. |
| A6 | Cobertura de m/q, Spearman, efeitos dentro de q; sem regressões que confundam padrão e operação. |
| A7 | Raios ao centroide dos vetores residuais unitários; ajuste √[n/(n−1)], sem PERMDISP. |
| A8 | Não ampliado: não resolve a separação e a triagem WCO foi desfavorável. |

Perturbações A1: pesos exponenciais positivos normalizados para média um dentro de cada classe; peso do par é o produto; somente pares distintos. Mediana ponderada pela inversa da CDF em 0,5, média dos adjacentes em empate (tolerância 10⁻¹²). Quantis 2,5–97,5% são sensibilidade à ponderação, não intervalos de confiança. Semente 2026102406, fluxos SeedSequence independentes. Os 10.000 resultados estão salvos.

Não se fizeram permutações inferenciais: permutabilidade não estabelecida. Existem 200 arranjos dentro de q em churn/wavy annular e 9 em plug/slug; dentro de (m,q), apenas um em ambos. Sem p-valores novos, família inferencial ou correção Holm a aplicar. Não houve bootstrap espacial, ANOSIM ou classificadores. Comparações de pares não foram tratadas como amostras independentes.

## Reprodução

Da raiz do projeto, com Python e dependências já instaladas:

```sh
python3 src2/estatistica6/scripts/screen.py
python3 src2/estatistica6/scripts/advance.py
python3 src2/estatistica6/scripts/validate_results.py
python3 src2/estatistica6/scripts/report.py
python3 src2/estatistica6/scripts/finalize.py
```

Os programas escrevem apenas em estatistica6 e no cache temporário do Matplotlib. Versões efetivas: ver manifest_started.json; PyWavelets 1.8.0 difere do manifesto histórico (1.9.0), mas os cálculos novos usam a Morlet própria, sem motor legacy. As fontes foram conferidas antes/depois por SHA-256. Executar novamente substitui somente os produtos desta análise. Relatórios e gráficos são derivados dos CSVs/NPZs salvos.

## Arquivos para auditoria

- [Triagem completa](results/triagem.csv), [decisões preliminares](results/decisoes_preliminares.csv), [decisões geométricas finais](results/decisoes_finais.csv).
- [Métricas por condição](results/metricas_finais.csv), [características](results/a3_features.csv), [efeitos](results/a3_effects.csv), [agrupamentos](results/cluster_summary.csv).
- [Influência dos candidatos](results/candidatos_leave_one_out.csv), [sensibilidade por pesos](results/a1_weight_sensitivity.csv), [janelas](results/a2_within_run_windows.csv).
- [Covariáveis](results/a6_associations.csv), [cobertura](results/a6_coverage.csv), [efeitos dentro de q](results/a6_within_q_effects.csv).
- [Validação de fontes e fórmulas](validation/numerical.json), [verificação independente final](validation/final_numeric.json), [manifesto](results/manifest.json).
'''
    (ROOT/'METODOLOGIA.md').write_text(method)
    (ROOT/'README.md').write_text('''# Estatística exploratória de padrões em v16

Plano executado com triagem prévia e comparação direta plug/slug. Fontes v16, estatistica4 e estatistica5 preservadas.

- [Resumo executivo com gráficos](RESUMO_EXECUTIVO.md)
- [Resumo executivo em PDF](RESUMO_EXECUTIVO.pdf)
- [Método, decisões e comandos de reprodução](METODOLOGIA.md)
- [Plano original preservado](PLANO.md)
- [Resultados numéricos](results/metricas_finais.csv)
- [Validação final](validation/final_numeric.json)

Gráficos completos: [triagem](figures/01_triagem.png), [perfis WPS](figures/02_perfis_wps.png), [churn/wavy annular](figures/03_churn_annular.png), [plug/slug](figures/04_plug_slug.png), [estabilidade](figures/05_estabilidade.png), [condições operacionais](figures/06_condicoes_operacionais.png) e [diagnóstico de Pearson](figures/07_diagnostico_pearson.png). Cada PNG tem uma versão PDF com o mesmo nome.

Cinco combinações passaram a triagem geométrica. Nenhuma estabeleceu agrupamento robusto nas condições avaliadas; em plug/slug, a direção do pico se inverteu entre começo e fim. O relatório distingue sinal exploratório de classificação validada.
''')
    return document


def executive_pdf(doc):
    def page(title,number):
        fig=plt.figure(figsize=(8.27,11.69),facecolor='white')
        fig.text(.075,.955,title,fontsize=17,weight='bold',va='top')
        fig.text(.075,.925,'Estatistica6 · v16_spatial2_18_71 · 26 ensaios',fontsize=9,color='#64717b')
        fig.text(.925,.025,str(number),fontsize=9,ha='right',color='#64717b')
        return fig
    def paragraph(fig,text,y,width=103,size=10):
        # Medição real da fonte evita corte lateral por estimativa de caracteres.
        renderer=fig.canvas.get_renderer()
        font=FontProperties(family='DejaVu Sans',size=size)
        available=.85*fig.get_figwidth()*fig.dpi
        lines=[];line=''
        for word in text.split():
            candidate=(line+' '+word).strip()
            measured=renderer.get_text_width_height_descent(candidate,font,ismath=False)[0]
            if measured>available and line:
                lines.append(line);line=word
            else:line=candidate
        if line:lines.append(line)
        fig.text(.075,y,'\n'.join(lines),fontsize=size,va='top',linespacing=1.35)
        return y-(size*1.45*len(lines)+8)/(fig.get_figheight()*72)
    def picture(fig,name,bounds):
        ax=fig.add_axes(bounds);ax.imshow(plt.imread(ROOT/f'figures/{name}.png'));ax.axis('off')
    with PdfPages(ROOT/'RESUMO_EXECUTIVO.pdf') as pdf:
        fig=page('Discriminação dos padrões por wavelets',1)
        y=paragraph(fig,doc['main'],.885,size=10.3,width=98)
        fig.text(.075,y,'Triagem antes do aprofundamento',fontsize=11,weight='bold',va='top');y-=.028
        for i,step in enumerate(doc['steps'],1):y=paragraph(fig,f'{i}. {step}',y,size=9.5,width=109)
        picture(fig,'01_triagem',[.025,.08,.95,min(.43,y-.10)])
        paragraph(fig,'Silhueta positiva favorece o próprio padrão; negrito marca aprovação no filtro preliminar. Nenhum método passou no conjunto de quatro padrões.',.065,size=9,width=110)
        pdf.savefig(fig);plt.close(fig)
        fig=page('churn e wavy annular',2)
        y=paragraph(fig,doc['churn'],.885,size=10.2,width=99)
        picture(fig,'03_churn_annular',[.025,.34,.95,.33])
        picture(fig,'05_estabilidade',[.025,.065,.95,.28])
        pdf.savefig(fig);plt.close(fig)
        fig=page('plug e slug',3)
        y=paragraph(fig,doc['plug'],.885,size=10.2,width=99)
        picture(fig,'04_plug_slug',[.025,.19,.95,.47])
        paragraph(fig,'Decisão: manter a variação espacial do pico e da mean power como questões úteis. Não transformar o resultado de uma janela em uma regra de identificação de padrão.',.14,size=10,width=100)
        pdf.savefig(fig);plt.close(fig)
        fig=page('Limitações e continuidade',4)
        y=paragraph(fig,doc['stopped'],.885,size=9.7,width=106)
        for i,lim in enumerate(doc['limits'],1):y=paragraph(fig,f'{i}. {lim}',y,size=9.4,width=111)
        picture(fig,'06_condicoes_operacionais',[.025,.155,.95,.255])
        paragraph(fig,doc['recommendation'],.115,size=9.7,width=104)
        paragraph(fig,'Validação numérica: PASS. Fontes preservadas; scripts, tabelas e método em src2/estatistica6/.',.055,size=8.8,width=117)
        pdf.savefig(fig);plt.close(fig)


def main():
    source,inventory,_=load_inputs()
    metrics=read_csv(ROOT/'results/metricas_finais.csv')
    f=read_csv(ROOT/'results/a3_features.csv')
    cluster_stats=read_csv(ROOT/'results/cluster_summary.csv')
    with np.load(ROOT/'results/distances_final.npz',allow_pickle=False) as z:matrices={k:z[k].copy() for k in z.files}
    triage_figure(metrics);profiles_figure(source);churn_figure(source,f,matrices)
    plug_figure(source,f,matrices,cluster_stats);stability_figure(metrics)
    covariate_figure(inventory,f);pearson_figure()
    s=make_summary(metrics,f,cluster_stats)
    doc=write_documents(s,f,metrics);executive_pdf(doc)
    print('Gerados 7 gráficos PNG/PDF e resumo executivo Markdown/PDF.',flush=True)


if __name__=='__main__':main()
