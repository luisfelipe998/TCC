"""Figuras científicas em PT-BR, sem alterar os dados de imagem."""
import json
import os
import time
os.environ.setdefault('MPLCONFIGDIR','/tmp/wco-eventos-mpl')
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter,MaxNLocator
from analysis_common import ROOT,DATA,BASE,CFG,LAM,load_signals,crop_prepare,evaluate,save_map
from inspect_inputs import LABELS

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                    'axes.spines.right':False,'axes.titleweight':'bold','figure.facecolor':'white',
                    'savefig.facecolor':'white'})
COLORS={'plug':'#147d92','slug':'#b75a2a','other':'#687b8e'}
PATTERNS=['dispersed bubbles','plug','slug','churn','wavy annular']
VLABEL={'principal':'Principal','v5_mesmos_segmentos':'Núcleo v5','faixa30_49':'Linhas 30–49',
        'somente_media':'Somente remoção da média','bordas_menos16':'Centros −16 px','bordas_mais16':'Centros +16 px',
        'ancoras_opostas_a':'Centros −16/+16','ancoras_opostas_b':'Centros +16/−16'}
F=ROOT/'figures'
def dec(x,pos=None):return f'{x:.2f}'.replace('.',',')
def save(fig,name):
    fig.savefig(F/f'{name}.png',dpi=150,bbox_inches='tight')
    if '/' not in name:fig.savefig(F/f'{name}.svg',bbox_inches='tight')
    plt.close(fig)

def overview():
    out=ROOT/'results';coverage=pd.read_csv(out/'cobertura.csv').set_index('padrao').loc[PATTERNS]
    fig,ax=plt.subplots(figsize=(11,4.8),layout='constrained');x=np.arange(5);w=.25
    for off,key,label,col in [(-w,'aquisicoes','Elegíveis','#ccd5dc'),(0,'aquisicoes_com_evento','Com evento aceito','#458fa0'),(w,'aquisicoes_em_equivalentes_intra','Com parceiro do mesmo padrão','#153c58')]:
        bars=ax.bar(x+off,coverage[key],w,label=label,color=col)
        ax.bar_label(bars,padding=3,fontsize=9)
    ax.set_xticks(x,[LABELS[k] for k in PATTERNS]);ax.set_ylabel('Número de aquisições')
    ax.yaxis.set_major_locator(MaxNLocator(integer=True));ax.set_ylim(0,12.5)
    ax.set_title('Cobertura: pareamento intrapadrão em 5 das 26 aquisições')
    ax.legend(frameon=False,ncols=3,loc='upper left');save(fig,'cobertura')
    c=pd.read_csv(out/'contrastes_pareados.csv').query('variante=="principal" and escopo=="intra"')
    fig,axs=plt.subplots(1,2,figsize=(12,5),layout='constrained')
    for ax,comp,title in zip(axs,['distintos','fixos'],['Equivalentes − distintos em extensão','Equivalentes − recortes fixos']):
        for i,pattern in enumerate(['plug','slug']):
            g=c[(c.comparador==comp)&(c.padrao==pattern)]
            offsets=np.linspace(-.12,.12,len(g)) if len(g)>1 else np.array([0.])
            ax.scatter(i+offsets,g.diferenca,s=65,color=COLORS[pattern],alpha=.8,zorder=4)
            ax.plot([i-.22,i+.22],[g.diferenca.mean()]*2,color=COLORS[pattern],lw=4)
            ax.text(i,g.diferenca.max()+.022,f'{len(g)} par(es) de aquisições',ha='center',fontsize=9)
        ax.axhline(0,color='#505050',lw=1);ax.axhline(.05,color='#89949d',ls='--',lw=1,label='Margem descritiva de 0,05')
        ax.set_xticks([0,1],['Pistonado (plug)','Golfadas (slug)']);ax.set_xlim(-.5,1.5)
        ax.set_ylim(-.12,.25);ax.yaxis.set_major_formatter(FuncFormatter(dec));ax.set_ylabel('Diferença da média de WCO alternativa')
        ax.set_title(title);ax.legend(frameon=False,loc='lower left',fontsize=8)
    fig.suptitle('Cada ponto representa um par de aquisições; traços = médias descritivas\nOs pares compartilham aquisições e não são réplicas independentes',fontsize=12)
    save(fig,'contrastes')
    s=pd.read_csv(out/'contrastes_por_padrao.csv').query('escopo=="intra"')
    fig,axs=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    keys=list(VLABEL)
    for ax,comp in zip(axs,['distintos','fixos']):
        for pat in ['plug','slug']:
            g=s[(s.comparador==comp)&(s.padrao==pat)].set_index('variante').loc[keys]
            ax.plot(range(len(keys)),g.diferenca,'o-',label=LABELS[pat],color=COLORS[pat])
        ax.set_xticks(range(len(keys)),[VLABEL[k] for k in keys],rotation=45,ha='right',fontsize=8)
        ax.axhline(0,color='#505050',lw=1);ax.axhline(.05,color='#89949d',ls='--',lw=1)
        ax.set_title('Equivalentes − '+('distintos' if comp=='distintos' else 'recortes fixos'))
        ax.set_ylabel('Diferença de WCO alternativa');ax.yaxis.set_major_formatter(FuncFormatter(dec));ax.legend(frameon=False)
    fig.suptitle('Sensibilidade com os mesmos parceiros e o mesmo comprimento de segmento',fontsize=13)
    save(fig,'sensibilidades')
    # O controle entre padrões compara a mesma morfologia com rótulos diferentes.
    a=pd.read_csv(out/'agregados_eventos.csv').query('variante=="principal"')
    fig,ax=plt.subplots(figsize=(10,4.8),layout='constrained')
    groups=['intra_equivalentes','intra_distintos','entre_equivalentes','entre_distintos']
    labels=['Plug × plug / slug × slug\nequivalentes','Plug × plug / slug × slug\ndistintos','Plug × slug\nequivalentes','Plug × slug\ndistintos']
    for i,gname in enumerate(groups):
        g=a[a.grupo==gname].sort_values(['pattern_x','pattern_y','run_x','run_y'])
        for k,r in enumerate(g.itertuples()):
            color=COLORS.get(r.pattern_x,'#687b8e') if gname.startswith('intra') else '#687b8e'
            ax.scatter(i+(k-(len(g)-1)/2)*.035,r.escore,color=color,s=50)
    ax.set_xticks(range(4),labels);ax.set_ylabel('Média de WCO alternativa por par de aquisições')
    ax.yaxis.set_major_formatter(FuncFormatter(dec));ax.set_ylim(0,.5)
    ax.set_title('Semelhança local também pode ocorrer entre padrões diferentes')
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([],[],marker='o',ls='',color=COLORS[p],label=LABELS[p]) for p in ['plug','slug']]
              +[Line2D([],[],marker='o',ls='',color=COLORS['other'],label='Plug × slug (entre padrões)')],frameon=False,fontsize=9)
    save(fig,'entre_padroes')
    synth=pd.read_csv(out/'controles_sinteticos.csv')
    names={'evento_identico':'Evento idêntico','atraso_100':'Atraso conhecido: 100 px','evento_com_ruido':'Evento com ruído',
        'eventos_suaves_distintos':'Eventos suaves distintos','fundo_comum_eventos_distintos':'Fundo comum + eventos distintos',
        'ruidos_independentes':'Ruídos independentes','evento_cortado_na_borda':'Evento cortado na borda'}
    fig,ax=plt.subplots(figsize=(10,5),layout='constrained');y=np.arange(len(synth))
    ax.barh(y-.17,synth.escore,.34,label='Média no domínio principal',color='#147d92')
    ax.barh(y+.17,synth.maximo_integral,.34,label='Máximo em todo o mapa',color='#d9a252')
    ax.set_yticks(y,[names[k] for k in synth.controle]);ax.invert_yaxis();ax.set_xlim(0,1.06)
    ax.set_xlabel('WCO alternativa');ax.xaxis.set_major_formatter(FuncFormatter(dec));ax.legend(frameon=False,loc='lower right')
    ax.set_title('Controles sintéticos: coerência alta não comprova equivalência visual')
    save(fig,'controles_sinteticos')

def event_cards(signals,cat):
    folder=F/'eventos';folder.mkdir(exist_ok=True)
    for r in cat[cat.incluida].itertuples(index=False):
        pano=np.asarray(Image.open(DATA/r.aquisicao/'__panorama_v11.png').convert('L'))
        frame=np.asarray(Image.open(DATA/r.aquisicao/r.quadro_ancora).convert('L'))
        fig,axs=plt.subplots(3,1,figsize=(13,6),layout='constrained')
        s,e=int(r.segmento_inicio),int(r.segmento_fim)
        axs[0].imshow(pano[:,s:e],cmap='gray',vmin=0,vmax=255,aspect='auto',extent=[s,e,80,0])
        axs[0].axvspan(r.inicio,r.fim,color='#1aa2b7',alpha=.12)
        for x in [r.inicio,r.fim]:axs[0].axvline(x,color='#1aa2b7',lw=1.2)
        axs[0].axvline(r.ancora,color='#e3a12f',ls='--');axs[0].set_title('Evento anotado dentro do contexto de 2.048 px')
        axs[0].set_xlabel('Posição no panorama (px)');axs[0].set_ylabel('Linha')
        axs[1].imshow(frame,cmap='gray',vmin=0,vmax=255,aspect='auto',extent=[0,frame.shape[1],80,0])
        axs[1].axvline(r.coluna_quadro_ancora,color='#e3a12f')
        axs[1].set_title(f'Quadro que forneceu a faixa da âncora: {r.quadro_ancora}; não é rastreamento de identidade')
        axs[1].set_xlabel('Coluna do quadro (px)');axs[1].set_ylabel('Linha')
        for ax in axs[:2]:
            ax.axhline(18,color='#1aa2b7',lw=.7);ax.axhline(72,color='#1aa2b7',lw=.7)
        axs[2].plot(np.arange(s,e),signals[r.aquisicao]['principal'][s:e],color='#294e6a',lw=.8)
        axs[2].axvspan(r.inicio,r.fim,color='#1aa2b7',alpha=.12);axs[2].set_xlabel('Posição no panorama (px)')
        axs[2].set_ylabel('Intensidade média\n18–71 (níveis de cinza)')
        fig.suptitle(f'{r.event_id} | extensão operacional: {r.extensao} px | {LABELS[r.padrao]}',fontsize=13)
        save(fig,'eventos/'+r.event_id)

def five_patterns(inv,signals,cat):
    fig,axs=plt.subplots(5,3,figsize=(16,12),layout='constrained',gridspec_kw={'width_ratios':[1.3,1,1]})
    for row,pat in enumerate(PATTERNS):
        run=inv[inv.pattern==pat].run.iloc[0];g=cat[(cat.aquisicao==run)&cat.incluida]
        if len(g):s=int(g.sort_values('inicio').iloc[0].segmento_inicio);kind='primeiro evento aceito'
        else:s=(len(signals[run]['principal'])-2048)//2;kind='referência central sem anotação'
        pano=np.asarray(Image.open(DATA/run/'__panorama_v11.png').convert('L'))
        axs[row,0].imshow(pano[:,s:s+2048],cmap='gray',vmin=0,vmax=255,aspect='auto',extent=[s,s+2048,80,0])
        for r in g.itertuples():
            if s<=r.ancora<s+2048:
                axs[row,0].axvline(r.inicio,color='#1aa2b7',lw=.8);axs[row,0].axvline(r.fim,color='#1aa2b7',lw=.8)
        axs[row,0].set_xlim(s,s+2048)
        display=LABELS[pat] if pat in ['plug','slug','churn'] else f'{LABELS[pat]} ({pat})'
        axs[row,0].set_title(f'{display}\n{run} | {kind}',fontsize=10)
        axs[row,0].set_ylabel('Linha');axs[row,0].set_xlabel('Posição no panorama (px)')
        x=crop_prepare(signals[run]['principal'],s);m,_=evaluate(x,x)
        axs[row,1].plot(np.arange(2048),x,lw=.7,color='#294e6a');axs[row,1].set_xlabel('Posição no segmento (px)')
        axs[row,1].set_ylabel('Intensidade preparada\n(desvios padrão)');axs[row,1].set_title('Sinal após o recorte',fontsize=10)
        im=axs[row,2].pcolormesh(m['lags'][::4],LAM,np.abs(m['wxx'][:,::4]),shading='auto',cmap='magma',rasterized=True)
        axs[row,2].set_yscale('log');axs[row,2].set_xlabel('Lag da autocorrelação (px)')
        axs[row,2].set_ylabel('Comprimento de onda (px)');axs[row,2].set_title('WPS alternativa = |CWT(Rxx)|',fontsize=10)
        fig.colorbar(im,ax=axs[row,2],label='Módulo (unid. arbitrárias)',shrink=.85)
    fig.suptitle('Cinco padrões: imagem, sinal e WPS alternativa\nOs eixos de posição e lag descrevem coordenadas diferentes',fontsize=14)
    save(fig,'cinco_padroes')

def pair_figure(pid,label,signals,cat):
    out=ROOT/'results';rows=pd.read_csv(out/'resultados_eventos.csv');r=rows[(rows.pair_id==pid)&(rows.variante=='principal')].iloc[0]
    c=cat.set_index('event_id');rx,ry=c.loc[r.event_x],c.loc[r.event_y]
    x=crop_prepare(signals[r.run_x]['principal'],int(rx.segmento_inicio));y=crop_prepare(signals[r.run_y]['principal'],int(ry.segmento_inicio))
    path=ROOT/'maps'/f'{pid}.npz'
    if path.exists():
        with np.load(path) as z:m={k:z[k].copy() for k in ['wxx','wyy','wxy','wco','scales','lags','denominator_valid']}
    else:
        m,_=evaluate(x,y);save_map(path,m,x,y,pair_id=pid,uso='ilustracao_posterior',escore=float(r.escore))
    fig,axs=plt.subplots(7,1,figsize=(13,16),layout='constrained',gridspec_kw={'height_ratios':[.65,.65,.9,1,1,1,1.2]})
    for ax,rec,role in zip(axs[:2],[rx,ry],['x','y']):
        s,e=int(rec.segmento_inicio),int(rec.segmento_fim);pano=np.asarray(Image.open(DATA/rec.aquisicao/'__panorama_v11.png').convert('L'))
        ax.imshow(pano[:,s:e],cmap='gray',vmin=0,vmax=255,aspect='auto',extent=[-1024,1024,80,0])
        ax.axvline(rec.inicio-rec.ancora,color='#22b9c5');ax.axvline(rec.fim-rec.ancora,color='#22b9c5')
        ax.set_title(f'Sinal {role} — {LABELS[rec.padrao]} | {rec.name} | extensão {rec.extensao} px',fontsize=10);ax.set_xlabel('Posição relativa à âncora (px)');ax.set_ylabel('Linha')
    axs[2].plot(np.arange(2048)-1024,x,label=f'x: {LABELS[rx.padrao]} — {rx.name} (linha contínua)',lw=.85,color=COLORS[rx.padrao])
    axs[2].plot(np.arange(2048)-1024,y,label=f'y: {LABELS[ry.padrao]} — {ry.name} (linha tracejada)',lw=.95,ls='--',color=COLORS[ry.padrao],alpha=.85)
    axs[2].legend(frameon=False,fontsize=8);axs[2].set_ylabel('Sinal preparado');axs[2].set_xlabel('Posição relativa à âncora (px)')
    values=[np.abs(m['wxx']),np.abs(m['wyy']),np.abs(m['wxy']),m['wco']]
    titles=[f'WPS alternativa de x — {LABELS[rx.padrao]} = |CWT(Rxx)|',f'WPS alternativa de y — {LABELS[ry.padrao]} = |CWT(Ryy)|',
            'Módulo da WCS alternativa = |CWT(Rxy)|','WCO alternativa; denominador no lag zero']
    for i,(ax,arr,title) in enumerate(zip(axs[3:],values,titles)):
        kwargs=dict(vmin=0,vmax=1) if i==3 else {}
        im=ax.pcolormesh(m['lags'][::4],LAM,arr[:,::4],shading='auto',cmap='viridis' if i==3 else 'magma',rasterized=True,**kwargs)
        ax.set_yscale('log');ax.set_ylim(50,1024);ax.set_ylabel('Comprimento de onda (px)');ax.set_xlabel('Lag (px)')
        ax.set_title(title,fontsize=10);edge=2047-np.sqrt(2)*m['scales']
        for sign in [-1,1]:ax.plot(sign*edge,LAM,':',color='white',lw=1)
        if i==3:
            ax.axvline(-32,color='#ffbb55',lw=.8);ax.axvline(32,color='#ffbb55',lw=.8)
        fig.colorbar(im,ax=ax,label='WCO' if i==3 else 'Módulo (unid. arbitrárias)',pad=.015,shrink=.95)
    fig.suptitle(f'{label} | média no domínio principal: {r.escore:.4f}'.replace('.',',')+
                 f'\nPadrões comparados: {LABELS[rx.padrao]} × {LABELS[ry.padrao]}'+
                 '\nPontilhado branco: limite geométrico de suporte; traços laranja: lags −32 a +32 px',fontsize=12)
    save(fig,'pares/'+label.split(' ')[0].lower())

def main():
    t0=time.perf_counter();(F/'pares').mkdir(exist_ok=True)
    inv,signals=load_signals();cat=pd.read_csv(ROOT/'results/catalogo_eventos.csv')
    overview();event_cards(signals,cat);five_patterns(inv,signals,cat)
    ex=json.loads((ROOT/'results/exemplos.json').read_text())
    pair_figure(ex['favoravel'],'Favorável (extremo ilustrativo posterior)',signals,cat)
    pair_figure(ex['desfavoravel'],'Desfavorável (extremo ilustrativo posterior)',signals,cat)
    pairs=pd.read_csv(ROOT/'results/pares_eventos.csv')
    pid=pairs[pairs.grupo=='entre_equivalentes'].sort_values('pair_id').pair_id.iloc[0]
    pair_figure(pid,'Entre-padrões (primeiro par por ID)',signals,cat)
    (ROOT/'results/renderizacao.json').write_text(json.dumps(dict(segundos=time.perf_counter()-t0,figuras_png=len(list(F.rglob('*.png')))),indent=2)+'\n')
    print(f'Figuras concluídas em {time.perf_counter()-t0:.1f} s.')

if __name__=='__main__':main()
