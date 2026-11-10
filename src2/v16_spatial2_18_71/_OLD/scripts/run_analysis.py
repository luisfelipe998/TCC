"""Reproducible Schuck-focused analysis of the user-confirmed 26-run scope."""
import argparse
from itertools import combinations
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from wavelet_core import *
from plot_labels import pattern_label, method_label

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT/'results'
FIGURES = ROOT/'figures'
LAM = np.geomspace(50, 3000, 200)

def load(run, field='raw'):
    return np.load(ROOT/'intermediate'/f'{run}_input.npz')[field]

def small_maps(p):
    # Full-resolution calculation; only stored/rendered maps are decimated.
    step = max(1,int(np.ceil(len(p['lags'])/1800)))
    idx = np.unique(np.r_[np.arange(0,len(p['lags']),step),len(p['lags'])//2])
    q = {k:p[k][:,idx] for k in ('wxx','wyy','wxy','wco')}
    q['lags'] = p['lags'][idx]
    q['mean_x_all'] = np.abs(p['wxx']).mean(axis=1)
    q['mean_y_all'] = np.abs(p['wyy']).mean(axis=1)
    return q

def plot_pair(q, name1, name2, pattern, target, unbiased=False):
    fig = plt.figure(figsize=(17,11),layout='constrained')
    gs = fig.add_gridspec(3,4,width_ratios=[4,1.2,4,1.2])
    axes = [fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,2]),
            fig.add_subplot(gs[1,:]),fig.add_subplot(gs[2,:])]
    common_max = max(np.abs(q['wxx']).max(), np.abs(q['wyy']).max())
    data = [np.abs(q['wxx']),np.abs(q['wyy']),np.abs(q['wxy']),q['wco']]
    titles = [f'alternative WPS x — {name1}',f'alternative WPS y — {name2}',
              'Módulo da alternative WCS — CWT da correlação cruzada', 'alternative WCO — Schuck, Eq. 32.17']
    for i,(ax,z,title) in enumerate(zip(axes,data,titles)):
        vmax = common_max if i<2 else (None if i==2 or unbiased else 1)
        im = ax.pcolormesh(q['lags'],LAM,z,shading='auto',cmap='viridis' if i<3 else 'magma',vmin=0,vmax=vmax,rasterized=True)
        ax.axvline(0,color='white',lw=.5)
        xlabel = 'Defasagem (px)' if i<2 else 'Defasagem (px); positiva: x deslocado à direita de y'
        ax.set(title=title,xlabel=xlabel,ylabel='Comprimento de onda (px)',yscale='log',ylim=(3000,50))
        ax.set_yticks([50,100,200,500,1000,2000,3000],labels=['50','100','200','500','1000','2000','3000'])
        fig.colorbar(im,ax=ax,label='Módulo' if i<3 else 'alternative WCO (adimensional)')
    lim = max(np.nanmax(q['mean_x_all']),np.nanmax(q['mean_y_all']))*1.05
    for col,key in [(1,'x'),(3,'y')]:
        ax=fig.add_subplot(gs[0,col])
        ax.plot(q[f'mean_{key}_all'],LAM,color='tab:blue')
        ax.set(title='Média da\nalternative WPS',yscale='log',ylim=(3000,50),xlim=(0,lim))
        ax.set_xlabel('Módulo médio\n(todas as defasagens)',fontsize=9)
        ax.set_yticks([])
    mode = 'NÃO ENVIESADO — diagnóstico; pode exceder 1' if unbiased else 'ENVIESADO — estimador principal'
    fig.suptitle(f'{pattern_label(pattern)} | {mode}\nMorlet de energia unitária, ω₀ = 6; plano completo de defasagens e escalas',fontsize=12)
    fig.savefig(target,dpi=110)
    plt.close(fig)

def analyze_runs(inventory):
    summary,profiles,sensitivity = [],[],[]
    for row in inventory.itertuples():
        start=time.perf_counter()
        raw=load(row.run); x=prepare(raw); n=len(x)
        spectra={}
        for engine,est,label in [('legacy','biased','v5'),('morlet','biased','v16_biased'),('morlet','unbiased','schuck_unbiased')]:
            r=correlation(x,x,est)
            w,s=transform_correlation(r,LAM,engine)
            mag=np.abs(w); zero=mag[:,n-1]
            co,valid=safe_coherence(w,w,w,n-1)
            allmean=mag.mean(axis=1)
            spectra[label]=(zero,allmean)
            summary.append(dict(run=row.run,pattern=row.pattern,included=row.included,method=label,n=n,
                self_wco_zero_error=float(np.nanmax(np.abs(co[:,n-1]-1))),
                self_wco_max=float(np.nanmax(co)),
                wco_fraction_above_1=float(np.mean(co>1+1e-8)),
                peak_zero_lag_px=float(LAM[np.argmax(zero)]),
                peak_lag_mean_px=float(LAM[np.argmax(allmean)])))
            for j,lam in enumerate(LAM):
                profiles.append(dict(run=row.run,pattern=row.pattern,included=row.included,method=label,wavelength_px=lam,
                    zero_lag_wps=zero[j],lag_mean_all=allmean[j]))
            del w,mag,co
        # Isolate the discretization effect at identical physical scales.
        r=correlation(x,x)
        scales=LAM/FOURIER_FACTOR
        exact=morlet_cwt(r,scales)[:,n-1]
        for precision in (12,16):
            numerical=pywt.cwt(r,scales,LEGACY_WAVELET,method='fft',precision=precision)[0][:,n-1]/LEGACY_AMPLITUDE
            # The carrier 0.9549 differs slightly from exactly 6/(2*pi).
            relative=np.linalg.norm(np.abs(numerical)-np.abs(exact))/np.linalg.norm(np.abs(exact))
            sensitivity.append(dict(run=row.run,comparison=f'pywt_precision_{precision}_same_scales',relative_l2=float(relative)))
        for linear in (False,):
            alt=prepare(raw,linear=linear)
            altw=morlet_cwt(correlation(alt,alt),scales)[:,n-1]
            sensitivity.append(dict(run=row.run,comparison='demean_vs_linear_detrend_zero_lag',
                relative_l2=float(np.linalg.norm(np.abs(altw)-np.abs(exact))/np.linalg.norm(np.abs(exact)))))
        for field in ('center','upper','lower'):
            alt=prepare(load(row.run,field))
            altw=morlet_cwt(correlation(alt,alt),scales)[:,n-1]
            sensitivity.append(dict(run=row.run,comparison=f'{field}_vs_core_mean_zero_lag',
                relative_l2=float(np.linalg.norm(np.abs(altw)-np.abs(exact))/np.linalg.norm(np.abs(exact)))))
        fig,ax=plt.subplots(1,2,figsize=(11,4),layout='constrained')
        for label,(zero,mean) in spectra.items():
            ax[0].plot(LAM,zero/zero.max(),label=method_label(label))
            ax[1].plot(LAM,mean/mean.max(),label=method_label(label))
        for a,title in zip(ax,['alternative WPS na defasagem zero','Média da alternative WPS em todas as defasagens']):
            a.set(xscale='log',xlabel='Comprimento de onda (px)',ylabel='Módulo normalizado pelo máximo',title=title)
            a.legend(fontsize=8);a.grid(alpha=.2)
        fig.suptitle(f'{row.run} | {pattern_label(row.pattern)} | N={n}'+('\nFORA DO ESCOPO DE 26 REGISTROS' if not row.included else ''))
        fig.savefig(FIGURES/'runs'/f'{row.run}.png',dpi=130);plt.close(fig)
        print('run',row.run,f'{time.perf_counter()-start:.2f}s',flush=True)
    pd.DataFrame(summary).to_csv(RESULTS/'per_run_methods.csv',index=False)
    pd.DataFrame(profiles).to_csv(RESULTS/'profiles.csv',index=False)
    pd.DataFrame(sensitivity).to_csv(RESULTS/'sensitivity.csv',index=False)

def analyze_pairs(inventory,limit=None):
    selected=inventory[inventory.included]
    tasks=[(pattern,a,b) for pattern,g in selected.groupby('pattern') for a,b in combinations(sorted(g.run),2)]
    if limit is not None: tasks=tasks[:limit]
    summary=[]
    for k,(pattern,a,b) in enumerate(tasks):
        start=time.perf_counter()
        ra,rb=load(a),load(b);n=min(len(ra),len(rb))
        # Shared prefix is explicit; detrend and normalize AFTER selecting it.
        x,y=prepare(ra[:n]),prepare(rb[:n])
        baseline={}
        for engine,est,label in [('legacy','biased','v5_recomputed'),('morlet','biased','v16_biased'),('morlet','unbiased','schuck_unbiased')]:
            if engine=='legacy':
                # Literal v5 preprocessing order: preprocess full run, then truncate.
                xx,yy=prepare(ra)[:n],prepare(rb)[:n]
            else: xx,yy=x,y
            p=pair(xx,yy,LAM,estimator=est,engine=engine)
            co=p['wco']; j0=n-1
            metrics=dict(pattern=pattern,run_x=a,run_y=b,method=label,n=n,discarded_x=len(ra)-n,discarded_y=len(rb)-n,
                max_wco=float(np.nanmax(co)),
                wco_fraction_above_1=float(np.mean(co>1+1e-8)),
                max_zero_lag_wco=float(np.nanmax(co[:,j0])),median_zero_lag_wco=float(np.nanmedian(co[:,j0])))
            summary.append(metrics)
            if label=='v16_biased':
                assert metrics['max_wco']<=1+1e-8, metrics
                q=small_maps(p)
                plot_pair(q,a,b,pattern,FIGURES/'pairs'/f'{a}__{b}.png')
                # Compact rendering data retain complex coefficients and all-lag means.
                np.savez_compressed(ROOT/'intermediate'/f'{a}__{b}_maps.npz',wavelengths=LAM,**q)
            elif label=='schuck_unbiased' and not any(t[0]==pattern for t in tasks[:k]):
                plot_pair(small_maps(p),a,b,pattern,FIGURES/'diagnostics'/f'{pattern}_unbiased.png',unbiased=True)
            del p,co
        print('pair',k+1,'/',len(tasks),a,b,f'{time.perf_counter()-start:.2f}s',flush=True)
    pd.DataFrame(summary).to_csv(RESULTS/'pair_methods.csv',index=False)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--stage',choices=['runs','pairs','all'],default='all')
    ap.add_argument('--limit',type=int)
    args=ap.parse_args()
    for folder in ('runs','pairs','diagnostics'): (FIGURES/folder).mkdir(parents=True,exist_ok=True)
    inventory=pd.read_csv(RESULTS/'input_audit.csv')
    start=time.perf_counter()
    if args.stage in ('runs','all'): analyze_runs(inventory.head(args.limit) if args.limit else inventory)
    if args.stage in ('pairs','all'): analyze_pairs(inventory,args.limit)
    print('elapsed seconds',time.perf_counter()-start,flush=True)

if __name__=='__main__':
    main()
