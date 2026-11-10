"""Relatório e figuras somente a partir de resultados auditados em disco."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/estatistica7-matplotlib')
from core import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import FuncFormatter

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,
    'axes.titleweight':'semibold','axes.spines.top':False,'axes.spines.right':False,
    'text.color':'#243746','axes.labelcolor':'#243746','savefig.dpi':180,
    'pdf.fonttype':42,'ps.fonttype':42})
INK='#243746';BLUE='#247b9f';PURPLE='#8260a1';TEAL='#308c77';ORANGE='#c78328'
CMAP=LinearSegmentedColormap.from_list('signed',['#79558c','#faf6ed','#2b8b86'])
SCOPE={'churn_annular':'churn × wavy annular','plug_slug':'plug × slug','quatro_padroes':'Quatro padrões'}
POSITION={'prefix':'Começo','center':'Centro','suffix':'Fim'}
METRIC={
 'global_max':'Máximo global', 'zero_weighted':'Lag zero · ponderado',
 'scale_peak_weighted':'Máximos por escala · ponderados',
 'scale_peak_power05':'Máximos por escala · piso de potência 5%',
 'common_lag_peak':'Máximo em um lag comum', 'map_mean':'Média integral',
 'map_mean_weighted':'Média integral · ponderada','lag_q90_weighted':'Percentil 90 dos lags · ponderado',
 'area05_weighted':'Fração do mapa com WCO ≥ 0,5',
 'mean_short':'Média · 50–200 px','mean_middle':'Média · 200–800 px','mean_long':'Média · 800–3.000 px',
 'common_peak_short':'Máximo comum · 50–200 px','common_peak_middle':'Máximo comum · 200–800 px',
 'common_peak_long':'Máximo comum · 800–3.000 px',
 'self_area_log_l2':'Autocoerência · perfil de área', 'self_width_log_l2':'Autocoerência · perfil de meia largura',
 'self_area_bands':'Autocoerência · áreas em três bandas',
 'area_normalized_full':'Área normalizada · integral *','area_normalized_short':'Área normalizada · 50–200 px *',
 'area_normalized_middle':'Área normalizada · 200–800 px *','area_normalized_long':'Área normalizada · 800–3.000 px *',
 'zero_uniform_kernel':'Lag zero · distância de kernel *'}


def num(x,d=3):return f'{float(x):.{d}f}'.replace('.',',')
def thousand(n):return f'{int(n):,}'.replace(',','.')
def dec(ax,axis='y'):
    def formatted(x,p):
        # O layout pode mudar o localizador entre desenho e exportação.
        ticks=getattr(ax,'get_'+axis+'ticks')()
        step=float(np.min(np.diff(ticks))) if len(ticks)>1 else .1
        digits=next((d for d in range(2,7) if np.isclose(round(step,d),step,rtol=1e-8,atol=1e-12)),6)
        if abs(x)<.5*10**(-digits):x=0.
        return num(x,digits)
    getattr(ax,axis+'axis').set_major_formatter(FuncFormatter(formatted))
def grid(ax):ax.grid(axis='y',alpha=.18,zorder=0)
def footer(fig,text,bottom=.11):
    fig.get_layout_engine().set(rect=(0,bottom,1,1-bottom))
    fig.text(.015,.025,text,fontsize=9,linespacing=1.45)
def save(fig,name):
    for ext in ['png','pdf']:fig.savefig(ROOT/f'figures/{name}.{ext}',bbox_inches='tight',facecolor='white')
    plt.close(fig)


def observed_figure(scores):
    fig,axes=plt.subplots(1,3,figsize=(13.5,4.5),layout='constrained')
    for same,color,label in [(True,BLUE,'Mesmo padrão · 87 pares'),(False,ORANGE,'Padrões diferentes · 238 pares')]:
        f=scores[scores.same_pattern==same]
        values=np.sort(f.global_max)
        axes[0].step(values,np.arange(1,len(values)+1)/len(values),where='post',color=color,label=label,lw=2)
        axes[1].hist(f.peak_lambda_px,bins=np.geomspace(50,3100,16),weights=np.ones(len(f))/len(f),
                     histtype='step',color=color,lw=2,label=label)
        axes[2].scatter(f.scale_peak_weighted,f.common_lag_peak,s=12,alpha=.5,color=color)
    axes[0].set(title='Máximo alto também entre padrões',xlabel='Máximo global WCO',ylabel='Fração acumulada',xlim=(.4,1),ylim=(0,1.02))
    axes[0].legend(frameon=False,fontsize=8,loc='upper left');dec(axes[0]);dec(axes[0],'x')
    axes[1].set(xscale='log',xlabel='λ do máximo (px)',ylabel='Fração dos pares',title='Máximos no limite da grade')
    axes[1].set_xticks([50,200,800,3000],['50','200','800','3.000']);dec(axes[1])
    axes[1].text(.03,.94,f'{num(100*np.isclose(scores.peak_lambda_px,3000).mean(),1)}% em 3.000 px',transform=axes[1].transAxes,va='top')
    axes[2].plot([0,1],[0,1],color='#8e99a0',ls='--',lw=1)
    axes[2].set(xlim=(0,1),ylim=(0,1),xlabel='Média dos máximos por escala',ylabel='Máximo da média em um lag comum',title='Um lag diferente por escala eleva o resumo')
    dec(axes[2]);dec(axes[2],'x')
    for ax in axes:grid(ax)
    footer(fig,'26 ensaios, começo de 3.171 px, 200 escalas e todos os 6.341 lags. Pares compartilham ensaios: descrição, sem teste por par.\nÀ direita, a diagonal é o limite algébrico: máximo da média ≤ média dos máximos; ambos usam os mesmos pesos de potência.')
    save(fig,'01_maximos_observados')


def noise_figure(scores,noise,summary):
    fig,axes=plt.subplots(1,2,figsize=(12.4,5),layout='constrained')
    subsets=[scores[scores.same_pattern].global_max,scores[~scores.same_pattern].global_max,
             noise[(noise.n==3171)&(noise.model=='gaussian_white')].global_max,
             noise[(noise.n==3171)&(noise.model=='gaussian_ar1_095')].global_max]
    colors=[BLUE,ORANGE,PURPLE,TEAL]
    box=axes[0].boxplot(subsets,positions=range(4),widths=.48,showfliers=False,patch_artist=True,
                        medianprops={'color':INK,'linewidth':1.6})
    rng=np.random.default_rng(7)
    for i,(v,color,patch) in enumerate(zip(subsets,colors,box['boxes'])):
        patch.set_facecolor(color);patch.set_alpha(.2)
        axes[0].scatter(rng.uniform(i-.16,i+.16,len(v)),v,s=6,color=color,alpha=.32,zorder=2)
        axes[0].text(i,1.025,num(v.median()),ha='center',fontsize=10,weight='bold',color=color)
    axes[0].set(xticks=range(4),xticklabels=['Mesmo\npadrão\n87 pares','Outro\npadrão\n238 pares','Ruído\nbranco\n99 pares','AR(1)\nφ = 0,95\n99 pares'],
                ylabel='Máximo global WCO',ylim=(.35,1.07),title='3.171 px: independência permite máximos altos')
    axes[0].text(.02,.06,'Valores acima das caixas: medianas',transform=axes[0].transAxes,fontsize=8)
    dec(axes[0]);grid(axes[0])
    for model,color,offset,label in [('gaussian_white',PURPLE,-45,'Ruído branco'),('gaussian_ar1_095',TEAL,45,'AR(1), φ = 0,95')]:
        f=summary[summary.model==model].sort_values('n')
        axes[1].errorbar(f.n+offset,f.median_max,yerr=np.array([f.median_max-f.q05_max,f.q95_max-f.median_max]),
                         color=color,marker='o',capsize=4,lw=1.5,label=label)
    axes[1].set(xticks=[3171,5071,7001],xticklabels=['3.171','5.071','7.001'],xlabel='Comprimento do sinal (px)',
                ylabel='Máximo global WCO',ylim=(.35,1.04),title='Mesmo modelo, diferentes comprimentos')
    axes[1].legend(frameon=False,loc='lower left');dec(axes[1]);grid(axes[1])
    footer(fig,'Simulações: 99 pares independentes por modelo e comprimento; mesmo preparo, escalas de 50–3.000 px e todos os lags.\nCaixas: quartis; hastes: até 1,5 IQR. À direita: mediana e percentis 5–95 das simulações, sem intervalo de confiança populacional.',.14)
    save(fig,'02_controles_independentes')


def triage_figure(metrics):
    f=metrics[metrics.condition=='3171_prefix']
    names=CFG['cross_candidates']+CFG['self_candidates']+[m for m in METRIC if m not in CFG['cross_candidates']+CFG['self_candidates'] and m in f.metric.unique()]
    scopes=CFG['primary_scopes']
    vals=np.array([[f[(f.metric==m)&(f.scope==s)].iloc[0].min_class_silhouette for s in scopes] for m in names])
    fig,ax=plt.subplots(figsize=(10.5,9.7),layout='constrained')
    im=ax.imshow(vals,cmap=CMAP,norm=TwoSlopeNorm(vmin=-.6,vcenter=0,vmax=.2),aspect='auto')
    ax.set(xticks=range(3),xticklabels=['churn × wavy annular\n19 ensaios','plug × slug\n6 ensaios','Quatro padrões\n25 ensaios'],
           yticks=range(len(names)),yticklabels=[METRIC[m] for m in names])
    ax.tick_params(length=0,labelsize=9)
    for i in range(len(names)):
        for j in range(3):ax.text(j,i,num(vals[i,j]),ha='center',va='center',fontsize=9,color='white' if vals[i,j]<-.35 else INK)
    for y in [14.5,17.5]:ax.axhline(y,color='white',lw=3)
    cb=fig.colorbar(im,ax=ax,shrink=.65,pad=.03);cb.set_label('Menor silhueta média entre as classes');dec(cb.ax)
    ax.set_title('Triagem WCO: uma média global pode esconder uma classe ruim',pad=16)
    footer(fig,'Começo de 3.171 px. Positivo: cada classe ficou, em média, mais próxima de si do que da classe alternativa mais próxima.\nNenhuma das 22 representações passou o filtro completo: macro ≥ 0,10, cada classe > 0 e vizinhança ≥ acaso + 0,10.\n* Adendo após inspeção inicial. A área normalizada longa foi recusada nesta condição por matriz não semidefinida positiva.',.115)
    save(fig,'03_triagem_wco')


def stability_figure(metrics):
    fig,axes=plt.subplots(2,3,figsize=(12.5,7.2),layout='constrained')
    for row,scope in enumerate(['churn_annular','plug_slug']):
        for col,metric in enumerate(['global_max','common_lag_peak','map_mean']):
            ax=axes[row,col]
            for n,color in [(3171,BLUE),(5071,PURPLE),(7001,TEAL)]:
                f=metrics[(metrics.scope==scope)&(metrics.metric==metric)].set_index('condition')
                conditions=[f'{n}_{p}' for p in POSITION]
                if not all(c in f.index for c in conditions):continue
                values=np.array([f.loc[c,'min_class_silhouette'] for c in conditions])
                ax.plot(range(3),values,color=color,marker='o',label=f'{thousand(n)} px',lw=1.8)
            ax.axhline(0,color='#666',lw=1);grid(ax);dec(ax)
            ax.set(xticks=range(3),xticklabels=POSITION.values(),title=METRIC[metric])
            if col==0:ax.set_ylabel(SCOPE[scope]+'\nMenor média de classe')
            if col==2:ax.legend(frameon=False,fontsize=9)
    fig.suptitle('A compactação das duas classes não persistiu em todas as janelas',fontsize=14)
    footer(fig,'Silhueta média da classe menos favorecida: valores > 0 seriam necessários nas duas classes em todas as posições.\nMesmos 19 churn/wavy annular e mesmos 6 plug/slug nas comparações entre comprimentos. Escalas verticais variam por painel.\nJanelas podem se sobrepor e não são réplicas independentes; as linhas descrevem sensibilidade espacial.',.13)
    save(fig,'04_estabilidade')


def functional_figure():
    data=csv_read(ROOT/'results/functional.csv')
    fig,axes=plt.subplots(2,3,figsize=(13,7.3),layout='constrained')
    for row,scope in enumerate(['churn_annular','plug_slug']):
        for col,(observable,label) in enumerate([('mean','Média sobre os lags'),('peak','Máximo sobre os lags'),('zero','Valor em lag zero')]):
            ax=axes[row,col]
            for pos,color in [('prefix',BLUE),('center',PURPLE),('suffix',ORANGE)]:
                f=data[(data.condition=='3171_'+pos)&(data.scope==scope)&(data.observable==observable)]
                ax.plot(f.wavelength_px,f.balanced_contrast,color=color,lw=1.6,label=POSITION[pos])
                if pos=='prefix':ax.fill_between(f.wavelength_px,f.loo_min,f.loo_max,color=BLUE,alpha=.13)
            ax.axhline(0,color='#666',lw=.8);ax.set(xscale='log',xlabel='Comprimento de onda λ (px)',title=label)
            ax.set_xticks([50,200,800,3000],['50','200','800','3.000']);dec(ax);grid(ax)
            if col==0:ax.set_ylabel(SCOPE[scope]+'\nContraste intra − inter')
            if row==0 and col==2:ax.legend(frameon=False,fontsize=9)
    fig.suptitle('Perfis completos: diferenças locais dependem da escala e da posição',fontsize=14)
    footer(fig,'Contraste = média das duas médias intraclasse − média interclasse, em cada escala. Positivo favorece o mesmo padrão.\nSombra azul: mínimo–máximo ao retirar um ensaio no começo; é análise de influência, não intervalo de confiança.\nTodas as 200 escalas são mostradas; não foram selecionadas bandas por significância.',.13)
    save(fig,'05_contrastes_por_escala')


def null_figure(summary,replicates):
    fig,axes=plt.subplots(2,3,figsize=(12.5,6.9),layout='constrained')
    for row,scope in enumerate(['churn_annular','plug_slug']):
        for col,metric in enumerate(['global_max','common_lag_peak','map_mean']):
            ax=axes[row,col]
            f=replicates[(replicates.scope==scope)&(replicates.metric==metric)]
            result=summary[(summary.scope==scope)&(summary.metric==metric)].iloc[0]
            ax.hist(f.balanced_contrast,bins=18,color='#c3d7e1',edgecolor='white',label='199 coortes simuladas')
            ax.axvline(0,color='#8a8a8a',lw=.8,ls=':')
            ax.axvline(result.observed_contrast,color=PURPLE,lw=2.2,label='Observado')
            ax.set_title(METRIC[metric]);ax.set_xlabel('Contraste intra − inter');dec(ax,'x');grid(ax)
            ax.text(.96,.94,'p Holm = '+num(result.p_model_holm,2),transform=ax.transAxes,ha='right',va='top',fontsize=10,
                    bbox={'facecolor':'white','edgecolor':'none','alpha':.9,'pad':1})
            if col==0:ax.set_ylabel(SCOPE[scope]+'\nRéplicas de coorte')
            if row==0 and col==2:ax.legend(frameon=False,fontsize=8,loc='upper left',bbox_to_anchor=(0,.84))
    fig.suptitle('Sem evidência de vantagem intraclasse no modelo de fases',fontsize=14)
    footer(fig,'Contraste balanceado; alternativa positiva; Holm conjunto de seis testes. Um único sinal substituto por ensaio em cada réplica.\nModelo aproximado: a redetrendagem alterou o espectro de potência em mediana 4,3% (P95 16,4%; máximo 47,2%, norma L2 relativa).\nValores-p condicionais ao modelo; não testam generalização, causalidade nem equivalência dos padrões.',.14)
    save(fig,'06_modelo_de_fases')


def self_figure():
    with np.load(ROOT/'results/profiles_3171_prefix.npz',allow_pickle=False) as z:
        labels=z['labels'];area=z['self_area'];width=z['self_width']
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    for ax,values,title in zip(axes,[area,width],['Média da autocoerência sobre os lags','Meia largura da autocoerência / λ']):
        for g in PATTERNS:
            v=values[labels==g]
            for line in v:ax.plot(LAM,line,color=COLORS[g],alpha=.18,lw=.65)
            ax.plot(LAM,np.median(v,axis=0),color=COLORS[g],lw=1.8,label=f'{g} (n={len(v)})')
        ax.set(xscale='log',xlabel='Comprimento de onda λ (px)',title=title)
        ax.set_xticks([50,200,800,3000],['50','200','800','3.000']);dec(ax);grid(ax)
    axes[0].set_ylabel('Área média normalizada pelo número de lags')
    axes[1].set_ylabel('Razão adimensional');axes[1].legend(frameon=False,fontsize=8,loc='upper left')
    footer(fig,'Cada linha fina representa um ensaio de 3.171 px; linhas grossas são medianas de classe. Sem bandas de significância.\nA autocoerência em lag zero é 1 por construção. Suas áreas e larguras são descritores do próprio ensaio, não concordância entre ensaios.')
    save(fig,'07_autocoerencia')


def examples_figure(scores):
    with np.load(ROOT/'results/examples_3171_prefix.npz',allow_pickle=False) as z:ex={k:z[k] for k in z.files}
    fig,axes=plt.subplots(2,3,figsize=(12.8,7.4),layout='constrained')
    for ax,k in zip(axes.ravel(),[0,2,1,3,5,4]):
        i,j=ex['pairs'][k];a,b=ex['runs'][i],ex['runs'][j]
        row=scores[(scores.run_x==a)&(scores.run_y==b)].iloc[0]
        im=ax.pcolormesh(ex['display_lags'],LAM,ex['display_wco'][k],shading='auto',cmap='viridis',vmin=0,vmax=1,rasterized=True)
        ax.set(yscale='log',xlabel='Lag espacial (px)',ylabel='λ (px)',title=f'{ex["labels"][i]} × {ex["labels"][j]}\n{a} / {b}')
        ax.set_yticks([50,200,800,3000],['50','200','800','3.000']);ax.tick_params(labelsize=8)
        ax.text(.03,.06,'Máximo integral: '+num(row.global_max),transform=ax.transAxes,color='white',fontsize=9,
                bbox={'facecolor':'#202a31','alpha':.7,'edgecolor':'none','pad':3})
    fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.8,pad=.02,label='WCO')
    fig.suptitle('Exemplos fixados antes da triagem: faixas altas não são exclusivas do padrão',fontsize=13)
    footer(fig,'Somente a exibição foi reduzida a 301 lags; todos os resumos usaram os 6.341 lags. Escala de cor comum de 0 a 1.\nOs seis pares foram definidos antes de examinar estes resultados; os mapas ilustram o cálculo e não constituem teste estatístico.')
    save(fig,'08_mapas_exemplo')


def covariates_figure():
    inv=csv_read(ROOT/'results/inventory.csv')
    fig,axes=plt.subplots(1,2,figsize=(12.5,4.7),layout='constrained')
    table=inv.groupby(['pattern','q_code']).size().unstack(fill_value=0).reindex(PATTERNS)
    axes[0].imshow(table,aspect='auto',cmap='Blues',vmin=0,vmax=table.to_numpy().max())
    axes[0].set(xticks=range(len(table.columns)),xticklabels=table.columns,yticks=range(5),yticklabels=PATTERNS,
                xlabel='Código q',title='Pouco suporte para comparar sob a mesma condição')
    axes[0].tick_params(length=0)
    for i in range(5):
        for j in range(len(table.columns)):
            axes[0].text(j,i,str(table.iloc[i,j]),ha='center',va='center',color='white' if table.iloc[i,j]>=3 else INK)
    rng=np.random.default_rng(8)
    for i,g in enumerate(PATTERNS):
        values=inv[inv.pattern==g].n
        axes[1].scatter(values,rng.uniform(i-.12,i+.12,len(values)),s=35,color=COLORS[g])
    axes[1].set(yticks=range(5),yticklabels=PATTERNS,xlabel='Comprimento original do panorama (px)',
                title='Fixar o trecho não remove toda a confusão')
    axes[1].invert_yaxis();axes[1].grid(axis='x',alpha=.18)
    footer(fig,'Os descritores foram recalculados com N comum dentro de cada comparação. O comprimento original permanece associado ao padrão.\nSomente churn/wavy annular em q=40 têm ≥2 ensaios de cada classe no mesmo q (2 e 3); não há intercâmbio de rótulos dentro de (m,q).')
    save(fig,'09_condicoes_operacionais')


def make_summary(scores,metrics,noise_summary,null_summary):
    f=metrics[metrics.condition=='3171_prefix']
    def record(scope,metric):return f[(f.scope==scope)&(f.metric==metric)].iloc[0]
    noise3171=noise_summary[noise_summary.n==3171].set_index('model')
    validation=json.loads((ROOT/'validation/final.json').read_text())
    nlow=((scores.global_peak_power_x_fraction<.05)|(scores.global_peak_power_y_fraction<.05)).sum()
    summary=dict(
        sample=csv_read(ROOT/'results/inventory.csv').pattern.value_counts().to_dict(),
        primary_candidates=18,additional_candidates_attempted=5,usable_candidates_at_main_condition=int(f.metric.nunique()),
        scope_combinations_at_main_condition=len(f),initial_pass_count=int(f.screen_pass.sum()),
        positive_classes_all_six_positions=int(validation['stable_positive_classes_without_010_threshold']),
        observed_pairs=len(scores),observed_median_max=float(scores.global_max.median()),
        intra_median_max=float(scores[scores.same_pattern].global_max.median()),inter_median_max=float(scores[~scores.same_pattern].global_max.median()),
        upper_boundary_count=int(np.isclose(scores.peak_lambda_px,3000).sum()),upper_boundary_fraction=float(np.isclose(scores.peak_lambda_px,3000).mean()),
        low_power_peak_count=int(nlow),low_power_peak_fraction=float(nlow/len(scores)),
        median_common_lag_peak=float(scores.common_lag_peak.median()),median_scale_peak_weighted=float(scores.scale_peak_weighted.median()),
        noise_white_median_max=float(noise3171.loc['gaussian_white','median_max']),noise_ar1_median_max=float(noise3171.loc['gaussian_ar1_095','median_max']),
        noise_white_fraction_above09=float(noise3171.loc['gaussian_white','fraction_max_above09']),
        minimum_holm_phase_p=float(null_summary.p_model_holm.min()),
        churn_global_macro=float(record('churn_annular','global_max').macro_silhouette),
        churn_global_minclass=float(record('churn_annular','global_max').min_class_silhouette),
        four_class_best_macro=float(f[f.scope=='quatro_padroes'].macro_silhouette.max()),
        validation_max_wco_map_error=validation['original_engine_max_error'],
        conclusion='Útil como diagnóstico metodológico: máximos altos não evidenciam padrão. Nenhuma representação testada sustentou compactação de todas as classes nas janelas avaliadas; isso não prova ausência de informação na WCO.')
    json_write(ROOT/'results/executive_summary.json',summary)
    csv_save('metrics_all.csv',metrics)
    return summary


def paragraphs(s):
    lead='Vale manter a WCO no TCC como diagnóstico dos limites da similaridade. Nesta exploração, ela não sustentou um método estável para agrupar os padrões. O resultado mais claro foi que sinais independentes também produziram máximos próximos de 1.'
    steps=[
        'Conferi as fontes e recalculei a WCO da v16 para os 26 ensaios: churn 10, wavy annular 9, plug 3, slug 3 e dispersed bubbles 1. Mantive linhas 18–71 e as exclusões da origem.',
        'Usei comprimentos comuns: 3.171 px nos 26, 5.071 px nos mesmos 25 com classes replicadas e 7.001 px nos mesmos 19 churn/wavy annular; começo, centro e fim. Foram 2.388 mapas cruzados e 210 de autocoerência, sempre com todos os lags.',
        'Examinei 18 representações: máximos, lag zero, médias, fração acima de 0,5, três bandas e perfis de autocoerência. No adendo exploratório, testei cinco normalizações; quatro foram utilizáveis em 3.171 e 5.071 px.',
        'Avaliei seis comparações entre duas classes e a comparação conjunta das quatro: silhueta por classe, vizinho mais próximo, agrupamento por ligação média, retirada de um ensaio, perfis por escala e condições operacionais.',
        'Calibrei a interpretação com 594 pares independentes de ruído branco/AR(1) e 199 réplicas de coorte com fases aleatórias. Os seis contrastes do modelo de fases receberam correção de Holm.'
    ]
    maxima=(f'Nos 325 pares reais, a mediana do máximo foi {num(s["observed_median_max"])}: {num(s["intra_median_max"])} dentro do padrão e {num(s["inter_median_max"])} entre padrões. '
            f'{s["upper_boundary_count"]} máximos ({num(100*s["upper_boundary_fraction"],1)}%) caíram em λ = 3.000 px, limite superior da grade. '
            f'Em 3.171 px, as medianas foram {num(s["noise_white_median_max"])} para ruído branco independente e {num(s["noise_ar1_median_max"])} para AR(1) independente. '
            'Isso demonstra que um máximo alto, isoladamente, não é evidência de padrão comum ou acoplamento físico. Não demonstra que os ensaios reais sejam ruído.')
    discrimination=(f'Nenhuma das {s["usable_candidates_at_main_condition"]} representações utilizáveis passou o filtro inicial completo em qualquer uma das sete comparações. '
            'Mesmo retirando o limiar de silhueta macro ≥ 0,10, nenhuma manteve silhueta média positiva em todas as classes nas seis janelas de 3.171/5.071 px. '
            f'Para churn/wavy annular, o máximo global teve silhueta macro {num(s["churn_global_macro"])} no começo, mas a classe menos favorecida teve {num(s["churn_global_minclass"])}. '
            f'Nos quatro padrões, a maior silhueta macro inicial foi apenas {num(s["four_class_best_macro"])}. '
            'Vizinhança local e agrupamento não foram interpretados como acurácia preditiva.')
    plug=('Plug e slug foram comparados diretamente, com três ensaios por classe, em todas as posições de 3.171 e 5.071 px. '
          'O máximo global apresentou compactação das duas classes em algumas janelas de 5.071 px, mas falhou nas três posições de 3.171 px. '
          'As vantagens pontuais não sustentaram uma regra espacialmente estável. Os perfis por escala e as normalizações adicionais também não mudaram essa decisão.')
    phases=(f'Os testes de vantagem intraclasse dos três resumos fixados (máximo global, máximo em um lag comum e média integral) não indicaram contraste maior que o esperado pelo modelo de fases após Holm: p ajustado ≥ {num(s["minimum_holm_phase_p"],2)}. '
            'São testes condicionais a um modelo aproximado, não demonstração de equivalência ou ausência de informação. A redetrendagem alterou o espectro de potência em mediana 4,3% e chegou a 47,2% em uma realização.')
    limits=[
        'Amostra pequena e já explorada; plug/slug têm n=3 cada e dispersed bubbles não permite compactação intraclasse. Sem teste externo ou estimativa de generalização.',
        'Pares compartilham ensaios; escalas e janelas não são aquisições independentes. Em 5.071 px, começo/centro/fim de um slug são exatamente o mesmo trecho.',
        'Padrão, m, q e reconstrução se confundem. Apenas q=40 permite contraste simétrico churn/wavy annular (2 e 3 ensaios); não há separação desses efeitos em (m,q).',
        'Escalas em pixels, sem calibração física. Mantive o estimador sem máscara de cone de influência; os controles medem essa configuração, não todas as definições possíveis de coerência.',
        'A normalização da área na banda longa falhou na condição de matriz semidefinida positiva em 3.171/5.071 px e foi recusada como kernel. Isso não invalida o cálculo WCO original.'
    ]
    recommendation=('Manter no TCC os controles independentes, a concentração dos máximos na maior escala e a ausência de separação estável dos resumos testados. '
                    'Não usar “WCO alta” como prova de mesmo padrão. Para o objetivo de agrupamento, a distribuição e a largura da WPS da estatistica6 continuam sendo pistas mais promissoras, ainda sem validação. '
                    'A estatistica5 descreve semelhança de forma espectral; a estatistica6 mostra sua sensibilidade espacial; esta análise acrescenta a calibração necessária para interpretar a WCO.')
    return dict(lead=lead,steps=steps,maxima=maxima,discrimination=discrimination,plug=plug,phases=phases,limits=limits,recommendation=recommendation)


def markdown_report(s,p):
    text='# WCO: exploração e utilidade para o TCC\n\n**'+p['lead']+'**\n\n'
    text+='[PDF com os gráficos principais](RESUMO_EXECUTIVO.pdf) · [Método e reprodução](METODOLOGIA.md) · [Todas as métricas](results/metrics_all.csv)\n\n'
    text+='## O que foi feito\n\n'+'\n'.join(f'{i}. {v}' for i,v in enumerate(p['steps'],1))+'\n\n'
    for title,body,figure,alt in [
        ('Resultado que vale manter',p['maxima'],'02_controles_independentes','Controles com sinais independentes'),
        ('O objetivo de agrupamento',p['discrimination'],'03_triagem_wco','Triagem das representações WCO'),
        ('Plug e slug',p['plug'],'04_estabilidade','Sensibilidade à posição e ao comprimento'),
        ('Referência estatística de fases',p['phases'],'06_modelo_de_fases','Contrastes observados e modelo de fases')]:
        text+=f'## {title}\n\n{body}\n\n![{alt}](figures/{figure}.png)\n\n'
    text+='## Limitações\n\n'+'\n'.join('- '+v for v in p['limits'])+'\n\n'
    text+='## Decisão para o TCC\n\n'+p['recommendation']+'\n\n'
    text+=f'**Auditoria numérica: PASS.** WCO integral conferida contra o núcleo original, erro máximo {s["validation_max_wco_map_error"]:.2e}; contagens, silhuetas, testes, normalizações e hashes verificados. As 34 fontes diretas permaneceram inalteradas. [Auditoria](validation/final.json). Correção numérica não demonstra validade física ou representatividade.\n'
    (ROOT/'RESUMO_EXECUTIVO.md').write_text(text)


def page(pdf,title,index):
    fig=plt.figure(figsize=(8.27,11.69),facecolor='white')
    fig.text(.06,.958,'ESTATISTICA7  /  WCO  /  ANÁLISE EXPLORATÓRIA',fontsize=8,color='#657888',weight='bold')
    fig.text(.06,.919,title,fontsize=17,weight='bold')
    fig.text(.06,.035,'v16_spatial2_18_71 · 24/10/2026 · escalas espaciais em pixels',fontsize=8,color='#657888')
    fig.text(.94,.035,str(index),fontsize=8,ha='right',color='#657888')
    return fig


def para(fig,text,y,size=10.1,color=INK,bold=False,x=.06,width=.88):
    """Quebra pelo tamanho real da fonte, não por contagem de caracteres."""
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    prop=FontProperties(family='DejaVu Sans',size=size,weight='bold' if bold else 'normal')
    max_px=fig.bbox.width*width;lines=[]
    for paragraph in text.split('\n'):
        current=''
        for word in paragraph.split():
            candidate=(current+' '+word).strip()
            if current and renderer.get_text_width_height_descent(candidate,prop,False)[0]>max_px:
                lines.append(current);current=word
            else:current=candidate
        lines.append(current)
    lineheight=size*1.4/(72*fig.get_figheight())
    for line in lines:
        fig.text(x,y,line,fontproperties=prop,va='top',color=color)
        y-=lineheight
    return y-.009


def embed(fig,name,bottom,top):
    data=plt.imread(ROOT/f'figures/{name}.png')
    ax=fig.add_axes([.04,bottom,.92,top-bottom]);ax.imshow(data);ax.axis('off')


def pdf_report(s,p):
    with PdfPages(ROOT/'RESUMO_EXECUTIVO.pdf',metadata={'Title':'WCO: exploração e utilidade para o TCC','Author':'Análise local reproduzível','Subject':'Limites de similaridade e discriminação de padrões'}) as pdf:
        fig=page(pdf,'O achado útil: máximo alto não identifica o padrão',1)
        y=para(fig,p['lead'],.878,11,bold=True)
        embed(fig,'02_controles_independentes',.42,.79)
        y=para(fig,p['maxima'],.395,10.3)
        y=para(fig,'Decisão para o TCC',y-.012,12,bold=True)
        y=para(fig,p['recommendation'],y,10.1)
        assert y>.06,('page1',y)
        pdf.savefig(fig);plt.close(fig)
        fig=page(pdf,'Como a análise foi conduzida',2)
        y=.876
        for i,t in enumerate(p['steps'],1):y=para(fig,f'{i}. {t}',y,9.6)
        # A figura é alta: usa toda a largura disponível no espaço restante.
        embed(fig,'03_triagem_wco',.075,y+.005)
        pdf.savefig(fig);plt.close(fig)
        fig=page(pdf,'Separação instável, inclusive em plug e slug',3)
        y=para(fig,p['discrimination'],.876,10.2)
        embed(fig,'04_estabilidade',.30,.74)
        y=para(fig,p['plug'],.275,10.2)
        y=para(fig,'Critério exploratório',y-.008,11,bold=True)
        y=para(fig,'Na triagem: silhueta macro ≥ 0,10, média positiva em cada classe e vizinhança ≥ acaso + 0,10. A estabilidade exigiu classes positivas em todas as posições e vantagem de vizinhança. Os limites são critérios de trabalho, não testes de significância. A conclusão desfavorável também persistiu sem o corte macro de 0,10.',y,9.7)
        assert y>.06,('page3',y)
        pdf.savefig(fig);plt.close(fig)
        fig=page(pdf,'Referência estatística e limites da conclusão',4)
        y=para(fig,p['phases'],.876,10.0)
        embed(fig,'06_modelo_de_fases',.425,.76)
        y=para(fig,'Limitações que mudam a interpretação',.402,11.5,bold=True)
        for t in p['limits']:y=para(fig,'• '+t,y,9.25)
        y=para(fig,f'Auditoria: PASS. Erro máximo WCO vs. núcleo original: {s["validation_max_wco_map_error"]:.2e}. Fontes preservadas. Tabelas, figuras extras e fórmulas acompanham este PDF.',y,8.9,color='#516674')
        assert y>.06,('page4',y)
        pdf.savefig(fig);plt.close(fig)


def main():
    scores=csv_read(ROOT/'results/scores_3171_prefix.csv')
    metrics=pd.concat([csv_read(ROOT/'results/metrics.csv'),csv_read(ROOT/'results/metrics_addendum.csv')],ignore_index=True)
    noise=csv_read(ROOT/'results/independent_noise_controls.csv');noise_summary=csv_read(ROOT/'results/independent_noise_summary.csv')
    null_summary=csv_read(ROOT/'results/phase_null_summary.csv');null_reps=csv_read(ROOT/'results/phase_null_replicates.csv')
    observed_figure(scores);noise_figure(scores,noise,noise_summary);triage_figure(metrics)
    stability_figure(metrics);functional_figure();null_figure(null_summary,null_reps)
    self_figure();examples_figure(scores);covariates_figure()
    summary=make_summary(scores,metrics,noise_summary,null_summary);p=paragraphs(summary)
    markdown_report(summary,p);pdf_report(summary,p)
    print('9 figuras PNG/PDF, resumo executivo MD/PDF e tabela consolidada gerados.',flush=True)


if __name__=='__main__':main()
