"""Gráficos em PT-BR, com tabelas de origem e escalas explícitas."""
from pathlib import Path
import json
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/v16-perfis-mpl')
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from profile_core import ROOT, BASE, LAM, wc, distribution, hellinger, quadrature
from run_profiles import PT

OUT=ROOT/'figures'
OUT.mkdir(exist_ok=True)
RESULTS=ROOT/'results'
ORDER=list(PT.values())
COLORS=dict(zip(ORDER,['#208282','#ce8a13','#6d56a4','#cc554b','#4269a5']))
SHORT=dict(zip(ORDER,['Bolhas dispersas','Bolhas alongadas\n(plug)','Golfadas\n(slug)','Agitado\n(churn)','Anular ondulado']))
FMT=FuncFormatter(lambda x,pos: f'{x:g}'.replace('.',','))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,
    'axes.labelsize':10,'axes.spines.top':False,'axes.spines.right':False,
    'figure.facecolor':'white','savefig.facecolor':'white','pdf.fonttype':42})
FIGURES=[]


def save(fig,name):
    fig.savefig(OUT/f'{name}.png',dpi=170,bbox_inches='tight')
    fig.savefig(OUT/f'{name}.pdf',bbox_inches='tight')
    FIGURES.append(name)
    plt.close(fig)


def finish(ax, xlabel=None, ylabel=None, log=False, categorical_y=False):
    if xlabel: ax.set_xlabel(xlabel)
    if ylabel: ax.set_ylabel(ylabel)
    if log:
        ax.set_xscale('log')
        ax.set_xticks([50,100,200,500,1000,2000])
    ax.xaxis.set_major_formatter(FMT)
    if not categorical_y: ax.yaxis.set_major_formatter(FMT)
    ax.grid(alpha=.18)


def main():
    cfg=json.loads((ROOT/'config.json').read_text())
    mask=np.array(cfg['primary_mask']);axis=LAM[mask]
    norm=pd.read_csv(RESULTS/'profiles_normalized.csv')
    pairs=pd.read_csv(RESULTS/'all_pairs.csv')
    ranks=pd.read_csv(RESULTS/'ranking_summary.csv').set_index('metodo')
    pat=pd.read_csv(RESULTS/'ranking_by_pattern.csv')
    infl=pd.read_csv(RESULTS/'acquisition_deletion.csv')
    within=pd.read_csv(RESULTS/'within_acquisition_sensitivity.csv')
    inv=pd.read_csv(BASE/'results/input_audit.csv');inv=inv[inv.included]
    labels=inv.set_index('run').pattern.map(PT).to_dict()
    primary=norm[(norm.variante=='principal')&(norm.dominio=='principal')]

    fig,axes=plt.subplots(1,5,figsize=(17,4.4),sharex=True,sharey=True,layout='constrained')
    for ax,pattern in zip(axes,ORDER):
        g=primary[primary.padrao==pattern]
        a=g.pivot(index='aquisicao',columns='comprimento_onda_px',values='densidade_log')
        xx=a.columns.to_numpy()
        for row in a.to_numpy():ax.plot(xx,row,color=COLORS[pattern],alpha=.3,lw=.8)
        ax.fill_between(xx,a.min(),a.max(),color=COLORS[pattern],alpha=.12)
        ax.plot(xx,a.median(),color=COLORS[pattern],lw=2.2)
        count_label='1 aquisição' if len(a)==1 else f'{len(a)} aquisições'
        ax.set_title(f'{SHORT[pattern]}\n{count_label}')
        finish(ax,'Comprimento de onda (px)',log=True)
    axes[0].set_ylabel('Densidade descritiva em log(λ)\nárea de cada aquisição = 1')
    fig.suptitle('Perfis de WPS alternativa no lag zero • linhas 18–71\nLinhas finas: aquisições; linha espessa: mediana; faixa: mínimo–máximo, sem inferência',fontsize=13)
    save(fig,'01_perfis_cinco_padroes')

    selected=['wps_principal','wps_prefixo_par','v5_prefixo_par','wco_v16_18_71','wco_v5_18_71','wco_v5_20_59']
    names=['Perfil principal • registro inteiro','Perfil v16 • prefixo de cada par','Perfil v5 • prefixo de cada par',
           'Escore WCO v16 • linhas 18–71','Escore WCO v5 • linhas 18–71','Escore WCO v5 • linhas 20–59']
    fig,axes=plt.subplots(1,2,figsize=(13,5.7),layout='constrained',sharey=True)
    for ax,key,title in zip(axes,['A_pares','A_balanceada'],['Peso igual por par','Peso igual por padrão/combinação']):
        y=np.arange(len(selected))
        vals=ranks.loc[selected,key].to_numpy()
        ax.barh(y,vals,color=['#208282']*3+['#6d56a4']*3,alpha=.8,height=.58)
        for j,method in enumerate(selected):
            sub=infl[infl.metodo==method][key]
            ax.plot([sub.min(),sub.max()],[j,j],color='#222',lw=1.5)
            ax.plot([sub.min(),sub.max()],[j,j],'|',color='#222',ms=9)
            ax.text(vals[j]+.01,j+.23,f'{vals[j]:.3f}'.replace('.',','),fontsize=9)
        ax.set_yticks(y,names);ax.set_ylim(len(selected)-.5,-.5);ax.set_xlim(0,1)
        ax.axvline(.5,color='#333',ls='--',lw=1)
        finish(ax,'Concordância de ordenação A (maior = melhor)',categorical_y=True);ax.set_title(title)
    fig.suptitle('Comparação pela mesma pergunta de ordenação\nTraços: mínimo–máximo ao retirar uma aquisição; não são intervalos de confiança',fontsize=13)
    save(fig,'02_comparacao_referencias')

    sorted_runs=[r for pattern in ORDER for r in sorted(k for k,v in labels.items() if v==pattern)]
    matrix=np.zeros((26,26));index={r:i for i,r in enumerate(sorted_runs)}
    mainpairs=pairs[pairs.metodo=='wps_principal']
    for r in mainpairs.itertuples():matrix[index[r.aquisicao_x],index[r.aquisicao_y]]=matrix[index[r.aquisicao_y],index[r.aquisicao_x]]=r.valor
    fig,ax=plt.subplots(figsize=(10.7,9.1),layout='constrained')
    im=ax.imshow(matrix,vmin=0,vmax=1,cmap='viridis_r')
    ax.set_xticks(range(26),sorted_runs,rotation=90,fontsize=7)
    ax.set_yticks(range(26),[f'{r}' for r in sorted_runs],fontsize=7)
    counts=[sum(v==p for v in labels.values()) for p in ORDER]
    for b in np.cumsum(counts)[:-1]:ax.axhline(b-.5,color='white',lw=1.5);ax.axvline(b-.5,color='white',lw=1.5)
    for tick,run in zip(ax.get_yticklabels(),sorted_runs):tick.set_color(COLORS[labels[run]])
    cb=fig.colorbar(im,ax=ax,label='Distância de Hellinger • menor = mais semelhante')
    cb.ax.yaxis.set_major_formatter(FMT)
    ax.set_title('Todos os 325 pares • perfis principais\nBlocos: bolhas dispersas, bolhas alongadas, golfadas, agitado, anular ondulado\nDiagonal zero por identidade; aquisições compartilhadas geram dependência')
    save(fig,'03_matriz_distancias')

    methods=['wps_principal','inicio3171','centro3171','fim3171','wps_cosseno','wps_integral',
             'wps_conservador','media_lags_todos','media_lags_coi','media_lags_conservador',
             'faixa30_49','faixa20_59','somente_media','v5_mesmas_escalas']
    names=['Principal, registro inteiro','Início, N = 3.171','Centro, N = 3.171','Fim, N = 3.171',
           'Cosseno, mesma faixa','Faixa integral: 50–3.000 px','Faixa conservadora: 50–544 px',
           'Média em todos os lags','Média em lags sob COI','Média conservadora: 50–544 px',
           'Extração 30–49','Extração 20–59','Somente remoção da média','v5: mesmos valores de escala a']
    fig,ax=plt.subplots(figsize=(11,7.7),layout='constrained')
    y=np.arange(len(methods))
    ax.plot(ranks.loc[methods,'A_pares'],y,'o',label='Peso por par',color='#208282')
    ax.plot(ranks.loc[methods,'A_balanceada'],y,'s',label='Peso por padrão/combinação',color='#ce8a13')
    for i,m in enumerate(methods):ax.plot(ranks.loc[m,['A_pares','A_balanceada']],[i,i],color='#999',lw=.8)
    ax.set_yticks(y,names);ax.invert_yaxis();ax.set_xlim(0,1);ax.axvline(.5,color='#333',ls='--')
    finish(ax,'Concordância de ordenação A',categorical_y=True);ax.legend(loc='lower left')
    ax.set_title('Sensibilidades fixadas antes da avaliação\nFaixas distintas geram estatísticas distintas; nenhuma variante foi escolhida como nova principal')
    save(fig,'04_sensibilidades')

    fig,ax=plt.subplots(figsize=(12,5.2),layout='constrained')
    ms=['wps_principal','inicio3171','centro3171','fim3171','wco_v16_18_71']
    mnames=['Perfil inteiro','Início 3.171','Centro 3.171','Fim 3.171','Escore WCO v16']
    colors=['#208282','#64b8aa','#b1d8cd','#ce8a13','#6d56a4']
    for k,(m,name,col) in enumerate(zip(ms,mnames,colors)):
        values=pat[pat.metodo==m].set_index('padrao').reindex(ORDER).A_balanceada.to_numpy()
        ax.bar(np.arange(5)+(k-2)*.15,values,width=.14,label=name,color=col)
    ax.set_xticks(range(5),[SHORT[p] for p in ORDER]);ax.set_ylim(0,1);ax.set_xlim(-.6,4.6)
    ax.text(0,.2,'Sem par\nintrapadrão',ha='center')
    ax.axhline(.5,color='#333',ls='--');ax.legend(ncol=3,loc='upper left',bbox_to_anchor=(0,1.24))
    ax.set_ylabel('Concordância A por padrão\ncontroles cruzados balanceados')
    ax.yaxis.set_major_formatter(FMT);ax.grid(axis='y',alpha=.2)
    save(fig,'05_resultados_por_padrao')

    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    for same,col,name in [(False,'#9fa4ad','Entre padrões'),(True,'#208282','Mesmo padrão')]:
        g=mainpairs[mainpairs.mesmo_padrao==same]
        ratio=np.maximum(g.n_x,g.n_y)/np.minimum(g.n_x,g.n_y)
        axes[0].scatter(ratio,g.valor,s=17,alpha=.6,c=col,label=name)
    finish(axes[0],'Razão entre comprimentos do par','Distância de Hellinger')
    axes[0].legend();axes[0].set_title('Comprimento e distância: associação descritiva')
    for j,p in enumerate(ORDER):
        g=inv[inv.pattern.map(PT)==p]
        offsets=np.linspace(-.14,.14,len(g))
        axes[1].scatter(np.full(len(g),j)+offsets,g.n,c=COLORS[p],s=40)
    axes[1].set_xticks(range(5),[SHORT[p] for p in ORDER],fontsize=8)
    axes[1].set_ylabel('Comprimento do panorama (px)');axes[1].grid(axis='y',alpha=.2)
    axes[1].set_title('Comprimentos por padrão • cada ponto é uma aquisição')
    save(fig,'06_comprimento_e_confundimento')

    fig,axes=plt.subplots(5,3,figsize=(16,12),layout='constrained')
    for line,run in enumerate(cfg['representatives']):
        pattern=labels[run]
        with Image.open(ROOT.parent/'data_processed'/run/'__panorama_v11.png') as im:pano=np.asarray(im.convert('L'))
        raw=np.load(BASE/'intermediate'/f'{run}_input.npz')['raw'];n=len(raw)
        axes[line,0].imshow(pano,cmap='gray',vmin=0,vmax=255,aspect='auto',extent=[0,n,79.5,-.5])
        for h in [17.5,71.5]:axes[line,0].axhline(h,color='#d45d43',lw=.8)
        axes[line,0].set_title(f'{pattern}\n{run} • N={n}');axes[line,0].set_ylabel('Linha da imagem')
        axes[line,0].set_xlabel('Posição no panorama (px)')
        axes[line,1].plot(np.arange(n),wc.prepare(raw),lw=.65,color=COLORS[pattern])
        finish(axes[line,1],'Posição no panorama (px)','Intensidade preparada')
        for variant,col,name in [('principal',COLORS[pattern],'Lag zero'),('media_lags_todos','#777','Média em todos os lags')]:
            g=norm[(norm.aquisicao==run)&(norm.variante==variant)&(norm.dominio=='principal')]
            axes[line,2].plot(g.comprimento_onda_px,g.densidade_log,color=col,label=name)
        finish(axes[line,2],'Comprimento de onda (px)','Densidade em log(λ)',log=True)
        axes[line,2].legend(fontsize=8)
    fig.suptitle('Cinco exemplos predeterminados • panorama, sinal e perfis de WPS alternativa\nLinhas vermelhas delimitam 18–71; posição da imagem e lag da correlação são coordenadas distintas',fontsize=13)
    save(fig,'07_panorama_sinal_perfis')

    # Casos desfavoráveis selecionados depois das métricas, por regras registradas.
    pivot=pairs.pivot(index=['aquisicao_x','aquisicao_y'],columns='metodo',values='valor')
    meta=mainpairs.set_index(['aquisicao_x','aquisicao_y'])[['mesmo_padrao','padrao_x','padrao_y']]
    chosen=[]
    for title, idx in [('Controle cruzado mais semelhante',pivot.loc[~meta.mesmo_padrao,'wps_principal'].idxmin()),
                       ('Par intrapadrão menos semelhante',pivot.loc[meta.mesmo_padrao,'wps_principal'].idxmax())]:
        chosen.append(dict(tipo=title,aquisicao_x=idx[0],aquisicao_y=idx[1],
                           H=float(pivot.loc[idx,'wps_principal']),WCO=float(pivot.loc[idx,'wco_v16_18_71'])))
    same=pivot.loc[meta.mesmo_padrao]
    idx=(same.wps_principal.rank(method='min')+same.wco_v16_18_71.rank(method='min')).idxmin()
    chosen.append(dict(tipo='Perfis próximos com escore WCO baixo',aquisicao_x=idx[0],aquisicao_y=idx[1],
                       H=float(pivot.loc[idx,'wps_principal']),WCO=float(pivot.loc[idx,'wco_v16_18_71'])))
    pd.DataFrame(chosen).to_csv(RESULTS/'selected_cases.csv',index=False)
    fig,axes=plt.subplots(1,3,figsize=(15,4.4),layout='constrained')
    for ax,case in zip(axes,chosen):
        for r in [case['aquisicao_x'],case['aquisicao_y']]:
            g=primary[primary.aquisicao==r]
            ax.plot(g.comprimento_onda_px,g.densidade_log,label=f'{r}\n{labels[r]}')
        title=f'{case["tipo"]}\nH={case["H"]:.3f}; escore WCO={case["WCO"]:.3f}'.replace('.',',')
        ax.set_title(title);finish(ax,'Comprimento de onda (px)','Densidade em log(λ)',log=True);ax.legend(fontsize=8)
    fig.suptitle('Casos escolhidos após as métricas • interpretação, sem validação\nH: registros inteiros na faixa principal; WCO: prefixo comum, faixa integral e máximo ponderado',fontsize=13)
    save(fig,'08_casos_desfavoraveis')

    maps=[]
    examples=[(a,b,'Exemplo predeterminado') for a,b in cfg['example_pairs']]
    examples += [(r['aquisicao_x'],r['aquisicao_y'],r['tipo']) for r in chosen]
    scores=pd.read_csv(BASE/'results/preprocessing_comparison_pairs.csv')
    scores=scores[scores.variant=='baseline'].set_index(['run_x','run_y'])
    for a,b,kind in examples:
        a,b=sorted([a,b]);n=int(min(inv.set_index('run').loc[[a,b],'n']))
        source=BASE/'intermediate'/f'{a}__{b}_maps.npz'
        if source.exists():
            with np.load(source) as z:q={k:z[k].copy() for k in z.files}
            origin='Mapa reduzido herdado; máximo integral obtido da tabela da base'
        else:
            ra=np.load(BASE/'intermediate'/f'{a}_input.npz')['raw'][:n]
            rb=np.load(BASE/'intermediate'/f'{b}_input.npz')['raw'][:n]
            full=wc.pair(wc.prepare(ra),wc.prepare(rb),LAM)
            j=np.unique(np.r_[np.arange(0,len(full['lags']),max(1,int(np.ceil(len(full['lags'])/1800)))),n-1])
            q={k:full[k][:,j] for k in ['wxx','wyy','wxy','wco']};q['lags']=full['lags'][j]
            q['wavelengths']=LAM
            assert np.isclose(np.nanmax(full['wco']),scores.loc[(a,b),'max_wco_all'],atol=1e-10)
            origin='Mapa recalculado integralmente e reduzido apenas para apresentação'
            del full
        la=q['wavelengths'];lags=q['lags'];zero=int(np.flatnonzero(lags==0)[0])
        zx,zy=abs(q['wxx'][:,zero]),abs(q['wyy'][:,zero])
        calc=np.abs(q['wxy'])**2/(zx*zy)[:,None]
        valid=np.isfinite(q['wco']);eqerr=float(np.max(np.abs(calc[valid]-q['wco'][valid])))
        assert eqerr<1e-10
        maximum=float(scores.loc[(a,b),'max_wco_all'])
        fig,axes=plt.subplots(3,1,figsize=(11,9),layout='constrained')
        for values,run in [(zx,a),(zy,b)]:
            axes[0].plot(la,values/(quadrature(la)@values),label=f'{run}: {labels[run]}')
        finish(axes[0],'Comprimento de onda (px)','Densidade em log(λ)',log=True);axes[0].legend()
        axes[0].set_title('WPS alternativa do prefixo comum, no lag zero • faixa integral')
        for ax,arr,title,cmap in [(axes[1],abs(q['wxy']),'Módulo da WCS alternativa','magma'),
                                  (axes[2],q['wco'],'WCO alternativa • denominador no lag zero','viridis')]:
            im=ax.pcolormesh(lags,la,arr,shading='auto',cmap=cmap,
                vmin=0,vmax=(max(1.,maximum) if 'WCO' in title else None),rasterized=True)
            ax.set_yscale('log');ax.set_ylabel('Comprimento de onda (px)');ax.set_xlabel('Lag da correlação (px)')
            ax.set_title(title);cb=fig.colorbar(im,ax=ax);cb.ax.yaxis.set_major_formatter(FMT)
            bound=(n-1-np.abs(lags))*wc.FOURIER_FACTOR/np.sqrt(2)
            bound=np.maximum(bound,la[0])
            ax.plot(lags,bound,color='white',ls='--',lw=.85,label='Borda da COI geométrica')
            ax.set_ylim(la[0],la[-1]);ax.set_xlim(lags[0],lags[-1]);ax.legend(fontsize=7,loc='upper right')
            ax.yaxis.set_major_formatter(FMT);ax.xaxis.set_major_formatter(FMT)
        h=float(pivot.loc[(a,b),'wps_principal'])
        fig.suptitle(f'{kind}\n{a} × {b} • N={n} • H inteiro={h:.3f} • máximo WCO integral={maximum:.3f}'.replace('.',','),fontsize=12)
        name=f'par_{a}__{b}';save(fig,name)
        maps.append(dict(aquisicao_x=a,aquisicao_y=b,tipo=kind,figura=f'{name}.png',origem=origin,
            n=n,erro_equacao=eqerr,maximo_integral=maximum,denominadores_invalidos=int((~np.isfinite(q['wco'][:,zero])).sum())))
        print(f'Figura WCS/WCO: {a} × {b}',flush=True)
    pd.DataFrame(maps).to_csv(RESULTS/'illustrated_maps.csv',index=False)

    synth=np.load(RESULTS/'synthetic_controls.npz')
    ctrl=json.loads((RESULTS/'operations_validation.json').read_text())
    fig,axes=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    for field,name in [('seno100','Senoide de 100 px'),('seno500','Senoide de 500 px')]:
        p=synth[field];axes[0].plot(LAM,p/(quadrature(LAM)@p),label=name)
    for field,name in [('phase_x','Realização 1'),('phase_y','Realização 2')]:
        p=synth[field];axes[1].plot(LAM,p/(quadrature(LAM)@p),label=name)
    for ax in axes:finish(ax,'Comprimento de onda (px)','Densidade em log(λ)',log=True);ax.legend()
    axes[0].set_title(f'Escalas distintas • H={ctrl["H_componentes_100_500"]:.3f}'.replace('.',','))
    pc=ctrl['fases_distintas']
    axes[1].set_title(f'Fases distintas • H={pc["H"]:.3f}\nEscore WCO={pc["escore_WCO"]:.3f}'.replace('.',','))
    fig.suptitle('Controles sintéticos • verificações de comportamento, sem validação experimental')
    save(fig,'09_controles_sinteticos')
    (RESULTS/'figure_manifest.json').write_text(json.dumps(dict(figures=sorted(set(FIGURES)),language='pt-BR'),ensure_ascii=False,indent=2)+'\n')
    print(f'{len(set(FIGURES))} figuras únicas, PNG e PDF',flush=True)


if __name__=='__main__':
    main()
