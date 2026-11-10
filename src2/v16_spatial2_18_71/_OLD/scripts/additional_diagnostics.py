"""Finite-record, panorama seam, and ensemble-definition diagnostics."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from wavelet_core import *

ROOT=Path(__file__).resolve().parents[1]
LAM=np.geomspace(50,3000,200)

def main():
    inv=pd.read_csv(ROOT/'results/input_audit.csv')
    rows=[]
    for item in inv.itertuples():
        z=np.load(ROOT/'intermediate'/f'{item.run}_input.npz')
        raw=z['raw']; n=len(raw); x=prepare(raw)
        scales=LAM/FOURIER_FACTOR
        base=morlet_cwt(correlation(x,x),scales)[:,n-1]
        # A sensitivity experiment, not an automatic correction to the input.
        good=np.ones(n,dtype=bool)
        for seam in z['seams']:
            good[max(0,int(seam)-2):min(n,int(seam)+2)]=False
        fraction=float((~good).mean())
        error=float('nan')
        if good.sum()>2 and fraction<.5:
            alternative=np.interp(np.arange(n),np.flatnonzero(good),raw[good])
            w=morlet_cwt(correlation(prepare(alternative),prepare(alternative)),scales)[:,n-1]
            error=float(np.linalg.norm(np.abs(w)-np.abs(base))/np.linalg.norm(np.abs(base)))
        rows.append(dict(run=item.run,included=item.included,pattern=item.pattern,seam_neighbour_fraction=fraction,
            seam_interpolation_zero_lag_relative_l2=error,interpolation_test_performed=np.isfinite(error)))
    pd.DataFrame(rows).to_csv(ROOT/'results/seam_sensitivity.csv',index=False)

    # Boundary/estimator diagnostic: first pair from the confirmed churn scope.
    name='m03q05_001'
    raw=np.load(ROOT/'intermediate'/f'{name}_input.npz')['raw']; x=prepare(raw); n=len(x)
    lines={}
    fig,axes=plt.subplots(1,3,figsize=(15,4),layout='constrained')
    for engine,est,label in [('legacy','biased','v5'),('morlet','biased','v16 biased'),('morlet','unbiased','Schuck unbiased')]:
        r=correlation(x,x,est)
        w,s=transform_correlation(r,LAM,engine)
        co,_=safe_coherence(w,w,w,n-1)
        maxima=np.nanmax(co,axis=1)
        label_pt={'v5':'v5','v16 biased':'v16 enviesado','Schuck unbiased':'Schuck não enviesado'}[label]
        axes[0].plot(LAM,maxima,label=label_pt)
        axes[1].plot(np.arange(-n+1,n),r,label=label_pt,alpha=.7)
        axes[2].plot(LAM,np.abs(w[:,n-1])/np.abs(w[:,n-1]).max(),label=label_pt)
        lines[label]=dict(max_self_wco=float(np.nanmax(co)))
    axes[0].axhline(1,color='k',ls=':');axes[0].set(xscale='log',yscale='log',xlabel='Comprimento de onda (px)',ylabel='Máximo da alternative WCO\n(sinal consigo mesmo)',title='Máximo em todas as defasagens')
    axes[1].set(xlabel='Defasagem (px)',ylabel='Autocorrelação',title='Estimadores enviesado e não enviesado')
    axes[2].set(xscale='log',xlabel='Comprimento de onda (px)',ylabel='alternative WPS normalizada\nna defasagem zero',title='Forma do perfil por escala')
    for ax in axes: ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle(name);fig.savefig(ROOT/'figures/diagnostics/autocorrelation_estimators.png',dpi=140);plt.close(fig)

    # K>1 ambiguity in the printed chapter. Deliberately controlled simulation;
    # not claimed as a pixel-for-pixel reproduction of undocumented MATLAB code.
    rng=np.random.default_rng(32017); n=2048; k=3; delay=500
    t=np.arange(n); waves=np.geomspace(8,160,100);scales=waves/FOURIER_FACTOR
    xx=[];yy=[];rrx=[];rry=[];rrxy=[];pcs=[];pxs=[];pys=[]
    for i,center in enumerate([350,650,950]):
        u=(t-center)/25
        signal=np.exp(-u*u)*np.cos(6*u)
        shifted=np.r_[np.zeros(delay),signal[:-delay]]
        # 5 dB defined here using whole-record mean square, explicitly.
        sigma=np.sqrt(np.mean(signal**2)/10**(.5))
        a=signal+sigma*rng.normal(size=n); b=shifted+sigma*rng.normal(size=n)
        a-=a.mean();b-=b.mean()
        xx.append(a); yy.append(b)
        rx=correlation(a,a);ry=correlation(b,b);rxy=correlation(a,b)
        rrx.append(rx);rry.append(ry);rrxy.append(rxy)
        pxs.append(morlet_cwt(rx,scales));pys.append(morlet_cwt(ry,scales));pcs.append(morlet_cwt(rxy,scales))
    mean_complex=np.mean(pcs,axis=0)
    mean_abs=np.mean(np.abs(pcs),axis=0)
    denominator=np.mean(np.abs(pxs),axis=0)[:,[n-1]]*np.mean(np.abs(pys),axis=0)[:,[n-1]]
    coherent=np.abs(mean_complex)**2/denominator
    magnitude_averaged=mean_abs**2/denominator
    commutation=float(np.max(np.abs(mean_complex-morlet_cwt(np.mean(rrxy,axis=0),scales))))
    magnitude_noncommutation=float(np.linalg.norm(np.mean(np.abs(pxs),axis=0)-np.abs(morlet_cwt(np.mean(rrx,axis=0),scales)))/np.linalg.norm(np.mean(np.abs(pxs),axis=0)))
    lags=np.arange(-n+1,n)
    peak=int(lags[np.unravel_index(np.argmax(np.abs(mean_complex)),mean_complex.shape)[1]])
    assert abs(peak+delay)<=2
    assert commutation<1e-12
    fig,axes=plt.subplots(2,2,figsize=(12,7),layout='constrained')
    for ax,data,title in zip(axes.flat,[np.abs(mean_complex),mean_abs,coherent,magnitude_averaged],
        ['Eq. 32.16: módulo da média complexa','Eq. 32.20: média dos módulos',
         'alternative WCO usando a Eq. 32.16','alternative WCO usando a Eq. 32.20 impressa']):
        im=ax.pcolormesh(lags[::3],waves,data[:,::3],shading='auto',cmap='viridis')
        ax.axvline(-500,color='white',ls=':');ax.set(xlabel='Defasagem na convenção do artigo (amostras)',ylabel='Comprimento de onda (amostras)',yscale='log',title=title)
        fig.colorbar(im,ax=ax,label='Módulo' if ax in axes[0] else 'alternative WCO')
    fig.suptitle('Controle inspirado em Schuck: 3 realizações, atraso de 500 amostras\nRelação sinal–ruído de 5 dB calculada sobre o registro inteiro')
    fig.savefig(ROOT/'figures/diagnostics/ensemble_definition_control.png',dpi=130);plt.close(fig)
    result=dict(self_example=lines,ensemble_complex_mean_commutation_error=commutation,
        ensemble_wps_magnitude_noncommutation_relative_l2=magnitude_noncommutation,ensemble_peak_lag=peak,
        note='Controlled simulation, not an exact reproduction of the chapter figures.')
    (ROOT/'results/additional_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
