"""Piloto óptico de âncoras; não importa nem calcula wavelets."""
from pathlib import Path
import json
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/wco-eventos-mpl')
import numpy as np
import pandas as pd
from PIL import Image
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT.parent/'v16_spatial2_18_71'
DATA=ROOT.parent/'data_processed'
PILOT=['m03q05_001','m05q40_001','m03q10_001','m05q80_001']

def main():
    out=ROOT/'results/extensao_churn_anular';rows=[]
    inv=pd.read_csv(BASE/'results/input_audit.csv').query('included').sort_values('run')
    for run in PILOT:
        pano=np.asarray(Image.open(DATA/run/'__panorama_v11.png').convert('L'))
        raw=pano[18:72].mean(axis=0);smooth=gaussian_filter1d(raw,16,mode='reflect')
        p05,p95=np.percentile(smooth,[5,95]);threshold=max(3.,.25*(p95-p05))
        candidates={}
        for kind,sign in [('minimo',-1),('maximo',1)]:
            pos,prop=find_peaks(sign*smooth,prominence=threshold,width=(16,512),distance=128,wlen=1025)
            keep=(pos>=1040)&(pos<=len(raw)-1040)
            candidates[kind]=pos[keep]
            for i in np.flatnonzero(keep):
                rows.append(dict(aquisicao=run,padrao=inv.set_index('run').loc[run,'pattern'],tipo=kind,
                    ancora=int(pos[i]),prominencia=float(prop['prominences'][i]),largura_meia_prominencia=float(prop['widths'][i]),
                    limiar_prominencia=threshold,amplitude_robusta=p95-p05))
        fig,axs=plt.subplots(3,1,figsize=(15,8),layout='constrained')
        axs[0].imshow(pano,cmap='gray',vmin=0,vmax=255,aspect='auto',extent=[0,len(raw),80,0])
        axs[0].set_title('Panorama completo; candidatos ópticos sem consulta a WCO')
        for kind,color,marker in [('minimo','#137f8b','v'),('maximo','#c2742c','^')]:
            pos=candidates[kind]
            axs[1].scatter(pos,smooth[pos],c=color,marker=marker,s=30,label=kind)
            for x in pos:axs[0].axvline(x,color=color,alpha=.5,lw=.6)
        axs[1].plot(raw,color='#aaa',lw=.5,alpha=.6,label='Média 18–71')
        axs[1].plot(smooth,color='#333',lw=.8,label='Suavização gaussiana σ=16 px')
        axs[1].legend(ncols=4,fontsize=8);axs[1].set_ylabel('Intensidade média')
        # Exemplo determinado pela proximidade ao centro, sem usar wavelets.
        pos=candidates['minimo'];c=int(pos[np.argmin(abs(pos-len(raw)/2))]) if len(pos) else len(raw)//2
        axs[2].imshow(pano[:,c-1024:c+1024],cmap='gray',vmin=0,vmax=255,aspect='auto',extent=[c-1024,c+1024,80,0])
        axs[2].axvline(c,color='#137f8b',lw=2);axs[2].set_title(f'Mínimo mais próximo do centro: coluna {c}; contexto de 2.048 px')
        for ax in axs:ax.set_xlabel('Coluna do panorama (px)')
        fig.suptitle(f'{run} | {inv.set_index("run").loc[run,"pattern"]} | Piloto de viabilidade óptica',fontsize=14)
        fig.savefig(ROOT/'figures/extensao_churn_anular/piloto'/f'{run}.png',dpi=125);plt.close(fig)
    df=pd.DataFrame(rows);df.to_csv(out/'piloto_optico.csv',index=False)
    print(df.groupby(['aquisicao','tipo']).agg(candidatos=('ancora','size'),prom_mediana=('prominencia','median'),largura_mediana=('largura_meia_prominencia','median')).to_string())

if __name__=='__main__':main()
