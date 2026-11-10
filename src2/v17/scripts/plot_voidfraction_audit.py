#!/usr/bin/env python3
"""Figuras de diagnóstico da auditoria v14.1; somente saídas em v17."""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'src2/v17/resultados/voidfraction'
os.environ['MPLCONFIGDIR'] = str(OUT/'cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

a=pd.read_csv(OUT/'arithmetic_and_seam_checks.csv')
n=pd.read_csv(OUT/'raw_reconstruction_checks.csv')
order=['dispersed bubbles','plug','slug','churn','wavy annular']
labels=['Bolhas\nn=1','Plug\nn=3','Slug\nn=3','Churn\nn=10','Anular\nn=8']
colors=['#d69a00','#0072b2','#009e73','#d55e00','#aa66a4']
fig,axes=plt.subplots(2,2,figsize=(12.5,9),layout='constrained')
for k,(p,c) in enumerate(zip(order,colors)):
    sub=a[(a.pattern==p)&a.eligible]
    raw=n[(n.pattern==p)&n.eligible]
    axes[0,0].scatter(sub.alpha,1000*sub.flip_rate,c=c,label=p,s=46,edgecolor='white',linewidth=.4)
    jitter=[(i-(len(sub)-1)/2)*.045 for i in range(len(sub))]
    axes[0,1].scatter([k+j for j in jitter],sub.flip_over_rowwise_shuffle,c=c,s=38)
    axes[1,0].scatter([k+j for j in jitter],raw.common_1ms_switches_per_pixel_s,c=c,s=38)
    axes[1,1].scatter(100*sub.seam_comparison_share,100*sub.seam_flip_share,c=c,s=42)
    for ax,key,frame in [(axes[0,1],'flip_over_rowwise_shuffle',sub),
                          (axes[1,0],'common_1ms_switches_per_pixel_s',raw)]:
        ax.plot([k-.18,k+.18],[frame[key].median()]*2,color='black',lw=1.5)
axes[0,0].set(xlabel='Área circular condicional α',ylabel='1000 F [transições / comparação]',
              title='A. Ocupação e flips da mesma máscara')
axes[0,0].legend(fontsize=8,loc='upper left')
axes[0,1].set(ylabel='F / E(F | ocupação por linha)',title='B. Flips relativos à ocupação')
axes[1,0].set(ylabel='Mudanças por pixel por segundo\nobservadas em passos de 1 ms',
              title='C. Quadros nativos: janela comum de 70 ms')
for ax in [axes[0,1],axes[1,0]]:
    ax.set_xticks(range(5),labels)
axes[1,1].plot([0,25],[0,25],color='gray',ls='--',label='Participação proporcional')
axes[1,1].set(xlabel='Comparações situadas em junções [%]',ylabel='Flips situados em junções [%]',
              title='D. Concentração de transições nas junções')
axes[1,1].legend(fontsize=8)
for ax in axes.ravel():
    ax.grid(alpha=.2)
fig.suptitle('Diagnósticos exploratórios sem aprendizado de máquina — 25 aquisições elegíveis',fontsize=13)
fig.savefig(OUT/'diagnosticos_voidfraction.png',dpi=180)
fig.savefig(OUT/'diagnosticos_voidfraction.pdf')
