"""Render the saved numerical maps; update labels without changing estimates."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_analysis import ROOT,FIGURES,LAM,load,plot_pair,small_maps
from wavelet_core import pair,prepare
from plot_labels import pattern_label, method_label

def main():
    pairs=pd.read_csv(ROOT/'results/pair_methods.csv')
    pairs=pairs[pairs.method=='v16_biased']
    for i,row in enumerate(pairs.itertuples()):
        name=f'{row.run_x}__{row.run_y}'
        with np.load(ROOT/'intermediate'/f'{name}_maps.npz') as data:
            q={key:data[key] for key in data.files if key!='wavelengths'}
        plot_pair(q,row.run_x,row.run_y,row.pattern,FIGURES/'pairs'/f'{name}.png')
        if i%10==0: print('Figuras renderizadas:',i+1,'de',len(pairs),flush=True)
    for pattern,g in pairs.groupby('pattern'):
        row=g.iloc[0];a,b=load(row.run_x),load(row.run_y);n=min(len(a),len(b))
        p=pair(prepare(a[:n]),prepare(b[:n]),LAM,estimator='unbiased')
        plot_pair(small_maps(p),row.run_x,row.run_y,pattern,FIGURES/'diagnostics'/f'{pattern}_unbiased.png',unbiased=True)
    render_run_profiles()

def render_run_profiles():
    """Refaz somente os 31 perfis, a partir das estatísticas integrais salvas."""
    profiles=pd.read_csv(ROOT/'results/profiles.csv');inv=pd.read_csv(ROOT/'results/input_audit.csv').set_index('run')
    for run,g in profiles.groupby('run'):
        fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
        for method,gm in g.groupby('method'):
            for ax,field in zip(axes,['zero_lag_wps','lag_mean_all']):
                ax.plot(gm.wavelength_px,gm[field]/gm[field].max(),label=method_label(method))
        for ax,title in zip(axes,['alternative WPS na defasagem zero','Média da alternative WPS em todas as defasagens']):
            ax.set(xscale='log',xlabel='Comprimento de onda (px)',ylabel='Módulo normalizado pelo máximo',title=title)
            ax.legend(fontsize=8);ax.grid(alpha=.2)
        row=inv.loc[run]
        fig.suptitle(f'{run} | {pattern_label(row.pattern)} | N={row.n}'+('\nFORA DO ESCOPO DE 26 REGISTROS' if not row.included else ''))
        fig.savefig(FIGURES/'runs'/f'{run}.png',dpi=130);plt.close(fig)

if __name__=='__main__':main()
