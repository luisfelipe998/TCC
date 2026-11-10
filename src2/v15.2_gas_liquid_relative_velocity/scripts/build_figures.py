"""Publication-sized PT-BR figures. Every plotted coordinate comes from saved CSVs."""
import os,gzip,json
from pathlib import Path
import numpy as np
import pandas as pd
from PIL import Image
BASE=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(BASE/'intermediate/mplcache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
PATTERNS={'dispersed bubbles':'Bolhas dispersas','plug':'Plug','slug':'Slug','churn':'Churn','wavy annular':'Anular ondulado'}
COLORS={'dispersed bubbles':'#6f4c9b','plug':'#187a8c','slug':'#bf5b17','churn':'#537740','wavy annular':'#bb4b73'}

def save(fig,name):
    for ext in ('png','pdf','svg'):fig.savefig(BASE/'figures'/f'{name}.{ext}',dpi=300,bbox_inches='tight',facecolor='white')
    plt.close(fig)

def run():
    d=pd.read_csv(BASE/'results/object_relative_per_run.csv')
    inv=pd.read_csv(BASE.parent/'v15.1_gas_liquid_velocity/results/scenario_inventory.csv');inv=inv[inv.included]
    order=list(PATTERNS);d['order']=d.pattern.map({v:i for i,v in enumerate(order)});d=d.sort_values(['order','run']).reset_index(drop=True)
    fig,axs=plt.subplots(1,2,figsize=(10,10),sharey=True,gridspec_kw={'wspace':.16})
    small=pd.read_csv(BASE/'results/motion_resolved_object_summary.csv');small=small[(small.threshold==.7)&(small.group=='small')].set_index('run')
    large=pd.read_csv(BASE/'results/robust_meniscus_summary.csv');large=large[large.threshold==.7].set_index('run')
    for ax,prefix,title,data in zip(axs,('small','meniscus'),('Pequenos componentes ópticos','Extremidades de envelopes alongados'),(small,large)):
        for i,r in d.iterrows():
            if r[prefix+'_supported']:
                c=COLORS[r.pattern];v=abs(r[prefix+'_px_s'])/1000
                lo,hi=sorted(abs(data.loc[r.run,['q25_px_s','q75_px_s']].astype(float).to_numpy())/1000)
                ax.plot([lo,hi],[i,i],color=c,lw=1.3)
                ax.scatter(v,i,color=c if r[prefix+'_threshold_stable'] else 'white',edgecolor=c,s=38,zorder=3)
            else:ax.text(.97,i,'—',transform=ax.get_yaxis_transform(),ha='right',va='center',color='.6')
        ax.set(title=title,xlabel='Magnitude da velocidade óptica (10³ px/s)',xlim=(0,150))
        ax.grid(axis='x',alpha=.16);ax.set_axisbelow(True)
    axs[0].set_yticks(range(len(d)),[r.run+'  ·  '+PATTERNS[r.pattern] for r in d.itertuples()],fontsize=8)
    axs[0].invert_yaxis()
    fig.suptitle('Cobertura dos 26 cenários: movimento de objetos e interfaces',y=.96,fontsize=13)
    fig.text(.12,.035,'Preenchido: suporte nos três limiares. Aberto: suporte apenas na referência.\nBarras: intervalo interquartil entre trajetórias; não são intervalos de confiança. —: sem suporte.',fontsize=9)
    save(fig,'coverage_all_scenarios_PT_BR')
    paired=pd.read_csv(BASE/'results/paired_object_consensus.csv');paired=paired[paired.threshold_stable]
    fig,ax=plt.subplots(figsize=(6.5,5.2))
    for r in paired.itertuples():
        x=abs(r.meniscus_px_s)/1000;y=abs(r.small_px_s)/1000
        ax.scatter(x,y,color=COLORS[r.pattern],s=65,marker='o' if r.pattern=='plug' else 's',label=PATTERNS[r.pattern])
        ax.annotate(r.run,(x,y),xytext=(5,7),textcoords='offset points',fontsize=9)
    ax.set(xlim=(75,122),ylim=(75,120),xlabel='Translação de envelopes alongados (10³ px/s)',ylabel='Movimento de pequenos componentes (10³ px/s)',title='Comparação nos mesmos intervalos temporais')
    ax.grid(alpha=.16);ax.legend(frameon=False)
    fig.text(.12,.015,'Razão entre movimentos de objetos; não representa escorregamento líquido/gás.',fontsize=8)
    save(fig,'paired_objects_PT_BR')
    stable=d[d.small_threshold_stable].copy();v=abs(stable.small_px_s.to_numpy());matrix=v[:,None]/v[None,:]
    fig,ax=plt.subplots(figsize=(7,6));im=ax.imshow(np.log2(matrix),cmap='RdBu_r',vmin=-1.6,vmax=1.6)
    labels=[r.run+'\n'+PATTERNS[r.pattern] for r in stable.itertuples()]
    ax.set_xticks(range(len(stable)),labels,rotation=40,ha='right',fontsize=8);ax.set_yticks(range(len(stable)),labels,fontsize=8)
    for i in range(len(v)):
        for j in range(len(v)):ax.text(j,i,f'{matrix[i,j]:.2f}'.replace('.',','),ha='center',va='center',color='white' if abs(np.log2(matrix[i,j]))>.9 else 'black')
    ax.set(xlabel='Cenário no denominador',ylabel='Cenário no numerador',title='Razão de velocidades de objetos pequenos (px/s ÷ px/s)')
    cb=fig.colorbar(im,ax=ax,fraction=.046,pad=.04);cb.set_ticks([-1,0,1],labels=['0,5','1','2']);cb.set_label('Razão entre magnitudes')
    save(fig,'inter_recording_ratios_PT_BR')
    tracks=pd.read_csv(BASE/'results/motion_resolved_small_tracks.csv')
    fig,axs=plt.subplots(2,4,figsize=(12,6),sharex=True,sharey=True)
    for ax,(name,g) in zip(axs.flat,[(name,g) for name,g in tracks.groupby('run') if len(g)>=12]):
        ax.scatter(g.height_px,abs(g.signed_px_s)/1000,s=15,color=COLORS[g.pattern.iloc[0]],alpha=.75)
        ax.set_title(name,fontsize=10);ax.grid(alpha=.15)
    axs.flat[-1].axis('off');fig.supxlabel('Altura mediana do objeto (px)');fig.supylabel('Magnitude da velocidade óptica (10³ px/s)')
    fig.suptitle('Tamanho aparente não certifica um traçador do líquido')
    save(fig,'object_size_speed_PT_BR')
    # Direct visual evidence of accepted moving tracks. No generated/synthetic image is mixed with these panels.
    names=['m07q01_001','m05q01_001','m05q05_001','m07q05_001'];fig,axs=plt.subplots(4,3,figsize=(11,7))
    for row,name in enumerate(names):
        with gzip.open(BASE/'intermediate'/f'tracks_{name}_0.7.json.gz','rt') as f:ts=json.load(f)
        ts=[t for t in ts if t['accepted'] and t['object_class']=='small' and abs(t['signed_px_frame'])>2 and t['intervals']>=10]
        target=np.median([t['signed_px_frame'] for t in ts]);t=min(ts,key=lambda t:abs(t['signed_px_frame']-target));nodes=t['nodes']
        selected=[nodes[i] for i in [0,len(nodes)//2,len(nodes)-1]]
        lo=max(0,min(n['xmin'] for n in selected)-25);hi=min(1016,max(n['xmax'] for n in selected)+25)
        for ax,n in zip(axs[row],selected):
            im=np.array(Image.open(BASE.parent/'data_processed'/name/f"ImgA{n['frame']:06d}.jpg").convert('L'))
            ax.imshow(im,cmap='gray',vmin=0,vmax=255,interpolation='none',aspect='equal')
            ax.add_patch(Rectangle((n['xmin']-.5,n['ymin']-.5),n['width'],n['height'],fill=False,edgecolor='#0072b2',lw=1.2))
            ax.set(xlim=(max(0,n['x']-80),min(1016,n['x']+80)),ylim=(71,18),title=f"{name} · quadro {n['frame']} · objeto {t['track_id']}")
            ax.tick_params(labelsize=7);ax.title.set_fontsize(8)
    fig.supxlabel('Coluna da imagem (px)');fig.supylabel('Linha da imagem (px)');fig.tight_layout()
    save(fig,'moving_object_examples_PT_BR')

if __name__=='__main__':run()
