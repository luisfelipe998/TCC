"""Figuras PT-BR da extensão; padrões e sinais explicitamente identificados."""
import json
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/wco-eventos-mpl')
import numpy as np
import pandas as pd
from PIL import Image
from scipy.ndimage import gaussian_filter1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.ticker import ScalarFormatter,FuncFormatter
from analysis_common import ROOT,DATA,LAM,load_signals,crop_prepare,evaluate,save_map
from wavelet_core import correlation
from texture_common import OUT,MAPS,LABELS,GROUPS,VARIANT_LABELS,Evaluator,spec,CFG

FIG=ROOT/'figures/extensao_churn_anular'
COLORS={'churn':'#7545a1','wavy annular':'#147d64','entre':'#597386'}
SHORT={'churn':'Churn','wavy annular':'Anular ondulado','entre':'Churn × anular'}
COMPS={'fixos':'Âncoras − janelas fixas','minimos_distintos':'Compatíveis − larguras distintas',
       'polaridades_opostas':'Compatíveis − mínimo × máximo'}

def f(x,d=3):return f'{x:.{d}f}'.replace('.',',')

def save(fig,name):
    for ax in fig.axes:
        for axis in [ax.xaxis,ax.yaxis]:
            if isinstance(axis.get_major_formatter(),ScalarFormatter):
                axis.set_major_formatter(FuncFormatter(lambda value,pos:format(value,'g').replace('.',',')))
    fig.savefig(FIG/f'{name}.png',dpi=150,bbox_inches='tight');plt.close(fig)

def heat(ax,m,key,title,cmap='viridis'):
    val=m[key] if key=='wco' else abs(m[key])
    # Só apresentação: resolução suficiente para o PNG, preservando todos
    # os lags centrais. Escores e arquivos NPZ usam as 4.095 colunas originais.
    ix=np.unique(np.r_[np.arange(0,len(m['lags']),4),np.where(abs(m['lags'])<=64)[0],len(m['lags'])-1])
    kwargs={'vmin':0,'vmax':1} if key=='wco' else {}
    p=ax.pcolormesh(m['lags'][ix],LAM,val[:,ix],shading='auto',cmap=cmap,rasterized=True,**kwargs)
    ax.set_xlim(-2047,2047);ax.set_ylim(50,1024)
    ax.set_yscale('log');ax.set_yticks([50,100,250,500,1024],labels=['50','100','250','500','1.024'])
    ax.set_ylabel('Comprimento de onda λ (px)');ax.set_xlabel('Lag da correlação (px)');ax.set_title(title,fontsize=10)
    boundary=2047-np.sqrt(2)*m['scales']
    ax.plot(boundary,LAM,'w--',lw=.7);ax.plot(-boundary,LAM,'w--',lw=.7)
    ax.axvspan(-32,32,color='white',alpha=.15);ax.axvline(0,color='white',lw=.4)
    plt.colorbar(p,ax=ax,pad=.015,label='WCO' if key=='wco' else 'Módulo')

def anchor_examples(cat,signals):
    fig,axes=plt.subplots(4,2,figsize=(14,11),layout='constrained')
    for j,pattern in enumerate(CFG['padroes']):
        r=cat[(cat.padrao==pattern)&(cat.tipo=='minimo')].sort_values(['aquisicao','ancora']).iloc[0]
        c=int(r.ancora);lo=int(r.segmento_inicio);hi=int(r.segmento_fim)
        pano=np.asarray(Image.open(DATA/r.aquisicao/'__panorama_v11.png').convert('L'))
        color=COLORS[pattern];ax=axes[0,j]
        ax.imshow(pano[:,lo:hi],cmap='gray',vmin=0,vmax=255,aspect='auto',extent=[-1024.5,1023.5,79.5,-.5])
        ax.axvline(0,color=color,lw=2);ax.axhline(17.5,color='#edc85c',lw=.8);ax.axhline(71.5,color='#edc85c',lw=.8)
        ax.set_title(f'{LABELS[pattern]}\n{r.anchor_id}\nRecorte global [{lo}, {hi}); âncora {c}',fontsize=10)
        ax.set_xlabel('Posição relativa à âncora (px)');ax.set_ylabel('Linha da imagem')
        raw=signals[r.aquisicao]['principal'];sm=gaussian_filter1d(raw,16,mode='reflect');ix=np.arange(c-256,c+256)
        ax=axes[1,j];ax.plot(ix-c,raw[ix],color='#aaa',lw=.6,label='Média original 18–71')
        ax.plot(ix-c,sm[ix],color=color,label='Suavizada só para localizar σ=16')
        ax.axvline(0,color=color);ax.axvspan(r.inicio_meia_prominencia-c,r.fim_meia_prominencia-c,color=color,alpha=.15)
        ax.set_title('Detalhe óptico de 512 px; sombra = largura à meia proeminência',fontsize=9)
        ax.set_xlabel('Posição relativa à âncora (px)');ax.set_ylabel('Intensidade');ax.legend(fontsize=8)
        x=crop_prepare(raw,lo);ax=axes[2,j];ax.plot(np.arange(2048)-1024,x,color=color,lw=.7)
        ax.axvline(0,color='#333',ls='--',lw=.8);ax.set_xlabel('Posição relativa à âncora (px)');ax.set_ylabel('Sinal preparado (adimensional)')
        ax.set_title('As 2.048 amostras originais entram no preparo e na correlação',fontsize=9)
        m,_=evaluate(x,x);heat(axes[3,j],m,'wxx',f'WPS alternativa de {SHORT[pattern]} = |CWT(Rxx)|');del m
    fig.suptitle('Âncoras executadas nos dois padrões\nMínimo óptico central: referência operacional; não certifica uma onda ou bolha física',fontsize=15)
    save(fig,'ancoras_dos_padroes')

def coverage(cat):
    fig,axes=plt.subplots(1,3,figsize=(15,4.8),layout='constrained')
    ax=axes[0]
    for j,pattern in enumerate(CFG['padroes']):
        g=cat[cat.padrao==pattern];n=int((g.tipo=='minimo').sum());k=int((g.tipo=='maximo').sum())
        ax.bar(j,n,color=COLORS[pattern],label='Mínimos' if j==0 else None)
        ax.bar(j,k,bottom=n,color=COLORS[pattern],alpha=.35,label='Máximos' if j==0 else None)
        ax.text(j,n/2,str(n),ha='center',color='white',weight='bold');ax.text(j,n+k/2,str(k),ha='center',weight='bold')
    ax.set_xticks([0,1],['Churn\n10/10 aquisições','Anular ondulado\n9/9 aquisições']);ax.set_ylim(0,90)
    ax.set_ylabel('Âncoras selecionadas');ax.set_title('Quatro mínimos e quatro máximos\npor aquisição');ax.legend(fontsize=8)
    for j,kind in enumerate(['minimo','maximo'],start=1):
        ax=axes[j]
        for pattern in CFG['padroes']:
            for alt,ls in [('sigma12','-'),('sigma20','--')]:
                g=cat[(cat.padrao==pattern)&(cat.tipo==kind)];v=g['deslocamento_'+alt].abs().dropna().sort_values().to_numpy()
                ax.step(np.r_[0,v,64],np.r_[0,np.arange(1,len(v)+1)/len(g),len(v)/len(g)],where='post',color=COLORS[pattern],ls=ls,
                        label=f'{SHORT[pattern]}, {alt.replace("sigma","σ=")} ({len(v)}/{len(g)})')
        ax.axvline(16,color='#666',ls=':',lw=1);ax.set_xlim(-1,65);ax.set_ylim(0,1.04)
        ax.set_xlabel('Deslocamento absoluto da âncora (px)');ax.set_ylabel('Fração de todas as âncoras do tipo')
        ax.set_title('Correspondência local dos '+('mínimos' if kind=='minimo' else 'máximos'));ax.legend(fontsize=8,loc='lower right')
    fig.suptitle('Correspondência condicionada às âncoras originais; não repete a seleção por estrato\nAusências de correspondência não são consideradas deslocamento zero',fontsize=14)
    save(fig,'cobertura_estabilidade')

def contrast_plots():
    cp=pd.read_csv(OUT/'contrastes_pareados.csv');cp=cp[cp.variante=='principal']
    fig,axes=plt.subplots(2,3,figsize=(14,7.5),layout='constrained')
    for i,pattern in enumerate(CFG['padroes']):
        for j,co in enumerate(COMPS):
            ax=axes[i,j];g=cp[(cp.classe==pattern)&(cp.comparador==co)].sort_values(['run_x','run_y']);v=g.diferenca.to_numpy()
            ypos=np.linspace(-.2,.2,len(v))
            ax.scatter(v,ypos,s=23,color=COLORS[pattern],alpha=.65)
            ax.axvline(0,color='#333',lw=.8);ax.axvline(.05,color='#999',ls=':',lw=.8)
            if len(v):
                ax.scatter([v.mean()],[0],marker='D',s=100,color='#202020',zorder=4)
                ax.text(.02,.95,f'Média: {f(v.mean())}\n{len(v)} pares de aquisições; {int((v>0).sum())} positivos',transform=ax.transAxes,va='top',fontsize=10)
            ax.set_ylim(-.5,.5);ax.set_yticks([]);ax.set_xlim(-.4,.4);ax.set_xlabel('Diferença de WCO média')
            ax.set_title(f'{SHORT[pattern]}\n{COMPS[co]}',fontsize=11)
    fig.suptitle('Efeito da ancoragem nas mesmas combinações de aquisições\nCada ponto = média dos pares locais; losango = média entre pares de aquisições',fontsize=14)
    save(fig,'contrastes')
    # Exibe também o conjunto cruzado; não o interpreta como classificação.
    cs=pd.read_csv(OUT/'contrastes_por_padrao.csv');cs=cs[(cs.variante=='principal')&(cs.comparador=='fixos')]
    fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
    for j,cl in enumerate(['churn','wavy annular','entre']):
        r=cs[cs.classe==cl].iloc[0];ax.bar(j-.15,r.ancoras,.3,color=COLORS[cl],label='Mínimos compatíveis' if j==0 else None)
        ax.bar(j+.15,r.controle,.3,color='#a9a9a9',label='Janelas fixas' if j==0 else None)
        ax.text(j,max(r.ancoras,r.controle)+.015,f'Δ={f(r.diferenca)}\nn={r.pares_aquisicoes}',ha='center',fontsize=10)
    ax.set_xticks(range(3),['Churn × churn','Anular × anular','Churn × anular']);ax.set_ylabel('WCO média pareada');ax.set_ylim(0,.85);ax.legend()
    ax.set_title('Mesma banda e mesmos pares de aquisições\nCoerência entre padrões também pode ser alta')
    save(fig,'entre_padroes')

def sensitivity_plots():
    cs=pd.read_csv(OUT/'sensibilidade_conjunto_global.csv')
    cols=[(p,c) for p in CFG['padroes'] for c in COMPS];arr=np.full((8,6),np.nan);counts=np.zeros((8,6),int)
    for i,v in enumerate(CFG['variantes']):
        for j,(p,c) in enumerate(cols):
            g=cs[(cs.variante==v)&(cs.classe==p)&(cs.comparador==c)]
            if len(g):arr[i,j]=g.iloc[0].diferenca;counts[i,j]=g.iloc[0].pares_aquisicoes
    bound=max(.05,float(np.nanmax(abs(arr))))
    fig,ax=plt.subplots(figsize=(14,7),layout='constrained');im=ax.imshow(arr,cmap='RdBu',norm=TwoSlopeNorm(0,-bound,bound),aspect='auto')
    for i in range(8):
        for j in range(6):
            ax.text(j,i,f'{f(arr[i,j])}\n(n={counts[i,j]})',ha='center',va='center',fontsize=9,
                    color='white' if abs(arr[i,j])>bound*.6 else '#222')
    ax.set_yticks(range(8),[VARIANT_LABELS[v] for v in CFG['variantes']]);ax.set_xticks(range(6),[
        f'{SHORT[p]}\nvs '+{'fixos':'fixos','minimos_distintos':'distintos','polaridades_opostas':'mín. × máx.'}[c] for p,c in cols])
    ax.set_title('Sensibilidades: mínimos compatíveis menos o controle\n815 pares locais comuns às oito variantes; n = pares de aquisições',fontsize=13)
    fig.colorbar(im,ax=ax,label='Diferença de WCO média');save(fig,'sensibilidades')
    lo=pd.read_csv(OUT/'retirada_aquisicao.csv');fig,axes=plt.subplots(1,3,figsize=(14,4.7),layout='constrained')
    for j,co in enumerate(COMPS):
        ax=axes[j]
        for i,p in enumerate(CFG['padroes']):
            g=lo[(lo.classe==p)&(lo.comparador==co)];v=g.diferenca.dropna()
            ax.scatter(v,np.full(len(v),i),color=COLORS[p],s=25,alpha=.6)
            if len(v):ax.scatter(g.diferenca_completa.iloc[0],i,marker='D',color='#222',s=65,zorder=4)
            special=g[g.retirada=='m03q40_001']
            if len(special):ax.scatter(special.diferenca,np.full(len(special),i),s=130,facecolor='none',edgecolor='#ca4b28',lw=1.5,zorder=5)
        ax.axvline(0,color='#444',lw=.8);ax.set_yticks([0,1],['Churn','Anular ondulado']);ax.set_ylim(-.7,1.7)
        ax.set_xlabel('Diferença de WCO média');ax.set_title(COMPS[co],fontsize=10)
    fig.suptitle('Retirada de uma aquisição por vez\nLosango = conjunto completo; círculo aberto = sem m03q40_001; pontos não são intervalos de confiança',fontsize=13)
    save(fig,'retirada_aquisicao')

def example_pairs(cat,signals):
    examples=pd.read_csv(OUT/'exemplos.csv');results=pd.read_csv(OUT/'resultados_ancoras.csv').query("variante=='principal'").set_index('pair_id')
    ev=Evaluator(signals);cat=cat.set_index('anchor_id')
    for r in examples.itertuples():
        row=results.loc[r.pair_id];a=cat.loc[row.anchor_x];b=cat.loc[row.anchor_y]
        m,s,x,y=ev.pair(spec(a,'principal'),spec(b,'principal'))
        group=results[(results.pattern_x==r.classe)&(results.pattern_y==r.classe)&(results.grupo=='minimos_compativeis')]
        expected=group.escore.max() if r.tipo=='maior_escore' else group.escore.min() if r.tipo=='menor_escore' else r.escore
        assert np.isclose(s['escore'],expected,rtol=1e-11,atol=1e-11)
        path=MAPS/f'{r.pair_id}.npz'
        if not path.exists():save_map(path,m,x,y,**row.to_dict(),pair_id=r.pair_id,escolha=r.escolha)
        fig,axes=plt.subplots(4,2,figsize=(14,12),layout='constrained')
        color=COLORS[r.classe]
        for j,(anchor,letter) in enumerate([(a,'x'),(b,'y')]):
            pano=np.asarray(Image.open(DATA/anchor.aquisicao/'__panorama_v11.png').convert('L'));lo=int(anchor.segmento_inicio);hi=int(anchor.segmento_fim)
            ax=axes[0,j];ax.imshow(pano[:,lo:hi],cmap='gray',vmin=0,vmax=255,aspect='auto',extent=[-1024.5,1023.5,79.5,-.5])
            ax.axvline(0,color=color,lw=2);ax.set_xlabel('Posição relativa à própria âncora (px)')
            ax.set_title(f'Sinal {letter}: {LABELS[anchor.padrao]}\n{anchor.name}\nRecorte [{lo}, {hi}); centro {int(anchor.ancora)}',fontsize=10)
        ax=axes[1,0];ax.plot(np.arange(2048)-1024,x,color=color,lw=.8,label='x — '+SHORT[a.padrao])
        ax.plot(np.arange(2048)-1024,y,color=color,lw=.8,ls='--',label='y — '+SHORT[b.padrao]);ax.axvline(0,color='#777',lw=.6)
        ax.set_xlabel('Posição relativa à própria âncora (px)');ax.set_ylabel('Sinal preparado');ax.legend(fontsize=8)
        ax=axes[1,1];ax.plot(m['lags'],correlation(x,x),color=color,label='Rxx',lw=.8)
        ax.plot(m['lags'],correlation(y,y),color=color,ls='--',label='Ryy',lw=.8)
        ax.plot(m['lags'],correlation(x,y),color='#555',label='Rxy',lw=.8);ax.set_xlabel('Lag (px)');ax.set_ylabel('Correlação (divisor N)');ax.legend(fontsize=8)
        heat(axes[2,0],m,'wxx',f'WPS alternativa de x — {SHORT[a.padrao]}')
        heat(axes[2,1],m,'wyy',f'WPS alternativa de y — {SHORT[b.padrao]}')
        heat(axes[3,0],m,'wxy','Módulo da WCS alternativa — x × y')
        heat(axes[3,1],m,'wco',f'WCO alternativa — média principal {f(s["escore"],4)}')
        typ={'predeterminado':'Exemplo predefinido (primeiro ID)','maior_escore':'Maior WCO média','menor_escore':'Menor WCO média'}[r.tipo]
        fig.suptitle(f'{LABELS[r.classe]} × {LABELS[r.classe]}\n{typ} entre {len(group)} pares compatíveis selecionados: {f(s["escore"],4)}\n'
                     'Âncoras: mínimos ópticos em ambos os sinais; escore: WCO do contexto, não intensidade',fontsize=12)
        slug='churn' if r.classe=='churn' else 'anular';save(fig,f'par_{slug}_{r.tipo}');del m

def controls():
    c=pd.read_csv(OUT/'controles_sinteticos.csv');fig,ax=plt.subplots(figsize=(10,5),layout='constrained')
    x=np.arange(len(c));ax.bar(x-.18,c.escore,.36,color='#547da0',label='Média principal ±32 px')
    ax.bar(x+.18,c.maximo_integral,.36,color='#b1b8bf',label='Máximo de todo o mapa')
    labels=['Mínimos\ncompatíveis','Larguras\ndistintas','Fundo comum +\nlarguras distintas','Polaridade\ninvertida']
    ax.set_xticks(x,labels);ax.set_ylabel('WCO');ax.set_ylim(0,1.13);ax.legend()
    for j,r in enumerate(c.itertuples()):ax.text(j-.18,r.escore+.02,f(r.escore),ha='center',fontsize=9)
    ax.set_title('Controles sintéticos: WCO alta não demonstra identidade física\nO módulo ao quadrado conserva a WCO quando se troca y por −y')
    save(fig,'controles_sinteticos')

def main():
    FIG.mkdir(exist_ok=True);cat=pd.read_csv(OUT/'catalogo_ancoras.csv').query('selecionada');_,signals=load_signals()
    coverage(cat);anchor_examples(cat,signals);contrast_plots();sensitivity_plots();controls();example_pairs(cat,signals)
    print('Figuras da extensão geradas em PT-BR.')

if __name__=='__main__':main()
