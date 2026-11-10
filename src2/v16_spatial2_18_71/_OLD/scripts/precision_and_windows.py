"""Isolate numerical CWT precision and finite-record window effects."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from wavelet_core import *
from plot_labels import pattern_label, WINDOW_PT

ROOT=Path(__file__).resolve().parents[1]
LAM=np.geomspace(50,3000,200)

def main():
    inv=pd.read_csv(ROOT/'results/input_audit.csv'); selected=inv[inv.included]
    rows=[]
    for name in ['m04q40_001','m07q01_001','m03q05_001','m07q05_001']:
        raw=np.load(ROOT/'intermediate'/f'{name}_input.npz')['raw'];x=prepare(raw);n=len(x);r=correlation(x,x)
        for precision in (10,12,16):
            w,scales=transform_correlation(r,LAM,'legacy',precision=precision)
            co,_=safe_coherence(w,w,w,n-1)
            idx=np.unravel_index(np.nanargmax(co),co.shape)
            rows.append(dict(run=name,precision=precision,maximum_self_wco=float(co[idx]),
                wavelength_at_maximum=float(LAM[idx[0]]),lag_at_maximum=int(idx[1]-n+1),
                zero_wps_relative_to_peak_at_maximum=float(abs(w[idx[0],n-1])/abs(w[:,n-1]).max())))
    pd.DataFrame(rows).to_csv(ROOT/'results/precision_controls.csv',index=False)
    windows=[]
    # Same 3171-pixel length for every included run; preserve actual pixel spacing.
    common=int(selected.n.min())
    for item in selected.itertuples():
        raw=np.load(ROOT/'intermediate'/f'{item.run}_input.npz')['raw']; full=prepare(raw);n=len(raw)
        starts={'prefix':0,'center':(n-common)//2,'suffix':n-common}
        for window,start in starts.items():
            x=prepare(raw[start:start+common]);legacy=full[start:start+common]
            w=morlet_cwt(correlation(x,x),LAM/FOURIER_FACTOR)[:,common-1]
            for j,lam in enumerate(LAM):
                windows.append(dict(run=item.run,pattern=item.pattern,window=window,n=common,
                    wavelength_px=lam,zero_lag_wps=float(abs(w[j])),
                    old_preprocessed_window_mean=float(legacy.mean()),old_preprocessed_window_std=float(legacy.std())))
    df=pd.DataFrame(windows);df.to_csv(ROOT/'results/equal_length_windows.csv',index=False)
    fig,axes=plt.subplots(1,5,figsize=(18,4),layout='constrained')
    for ax,(pattern,g) in zip(axes,df.groupby('pattern')):
        for window,gw in g.groupby('window'):
            pivot=gw.pivot(index='run',columns='wavelength_px',values='zero_lag_wps')
            normalized=pivot.div(pivot.max(axis=1),axis=0)
            ax.plot(pivot.columns,normalized.mean(axis=0),label=WINDOW_PT[window])
        ax.set(xscale='log',xlabel='Comprimento de onda (px)',ylabel='Média da alternative WPS normalizada\nna defasagem zero',title=pattern_label(pattern),xlim=(50,3000))
        ax.legend(fontsize=7);ax.grid(alpha=.2)
    fig.suptitle(f'Sensibilidade à posição do segmento, comprimento comum N={common} px\nPerfis descritivos de alternative WPS')
    fig.savefig(ROOT/'figures/diagnostics/window_sensitivity.png',dpi=130);plt.close(fig)
    print('precision controls:',len(rows),'equal-length window profiles:',len(df))

if __name__=='__main__': main()
