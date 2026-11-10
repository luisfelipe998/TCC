"""Compare fixed preprocessing choices, keeping Schuck alternative formulas intact.

Same-pair power weights are fixed from the baseline. All 325 experimental
pairs include 87 same-pattern and 238 different-pattern controls. No alignment,
length warping, favorable-segment search or coherence clipping is performed.
"""
import argparse
from itertools import combinations
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
from PIL import Image
from scipy.ndimage import gaussian_filter1d, median_filter
from scipy.stats import rankdata
from wavelet_core import *

ROOT=Path(__file__).resolve().parents[1]
VARIANTS=['baseline','demean_only','center_row','central_band20','core_row_median',
          'core_trimmed_mean20','spatial_gaussian2','spatial_median5','rank_uniform','winsor_1pct']

def extract(run):
    image=np.asarray(Image.open(ROOT.parent/'data_processed'/run/'__panorama_v11.png').convert('L'),dtype=float)
    core=image[20:60];raw=core.mean(axis=0)
    sorted_rows=np.sort(core,axis=0)
    return dict(baseline=raw,demean_only=raw,center_row=image[40],central_band20=image[30:50].mean(axis=0),
                core_row_median=np.median(core,axis=0),core_trimmed_mean20=sorted_rows[8:32].mean(axis=0),
                spatial_gaussian2=raw,spatial_median5=raw,rank_uniform=raw,winsor_1pct=raw)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--scales',type=int,default=64)
    ap.add_argument('--variants',nargs='+',default=VARIANTS,choices=VARIANTS)
    ap.add_argument('--limit',type=int)
    ap.add_argument('--tag',default='screen')
    args=ap.parse_args()
    assert 'baseline' in args.variants
    lam=np.geomspace(50,3000,args.scales);scales=lam/FOURIER_FACTOR
    inv=pd.read_csv(ROOT/'results/input_audit.csv');inv=inv[inv.included].sort_values('run')
    inputs={r:extract(r) for r in inv.run};labels=inv.set_index('run').pattern.to_dict()
    tasks=list(combinations(inv.run,2))
    if args.limit:tasks=tasks[:args.limit]
    cache={};rows=[];start=time.perf_counter()
    def get(run,n,variant):
        key=(run,n,variant)
        if key not in cache:
            raw=inputs[run][variant][:n]
            # Filter only the observed common segment; no access past its end.
            if variant=='spatial_gaussian2':raw=gaussian_filter1d(raw,2,mode='reflect')
            if variant=='spatial_median5':raw=median_filter(raw,size=5,mode='reflect')
            if variant=='rank_uniform':raw=rankdata(raw,method='average')/(len(raw)+1)
            if variant=='winsor_1pct':raw=np.clip(raw,*np.quantile(raw,[.01,.99]))
            x=prepare(raw,linear=variant!='demean_only')
            zero=np.abs(morlet_cwt(correlation(x,x),scales)[:,n-1])
            cache[key]=(x,zero)
        return cache[key]
    for k,(a,b) in enumerate(tasks):
        n=min(len(inputs[a]['baseline']),len(inputs[b]['baseline']))
        _,bx=get(a,n,'baseline');_,by=get(b,n,'baseline')
        weight=np.sqrt(bx*by)
        # Baseline power thresholds are fixed across candidates in each pair.
        powered={cut:(bx>=cut*bx.max())&(by>=cut*by.max()) for cut in [.01,.05,.1]}
        for variant in args.variants:
            x,zx=get(a,n,variant);y,zy=get(b,n,variant)
            wxy=morlet_cwt(correlation(x,y),scales)
            valid=(zx>1e-10*zx.max())&(zy>1e-10*zy.max())
            co=np.divide(np.abs(wxy)**2,(zx*zy)[:,None],out=np.full(wxy.shape,np.nan),where=valid[:,None])
            assert np.nanmax(co)<=1+1e-8
            peak=np.full(len(scales),np.nan)
            peak[valid]=np.max(co[valid],axis=1)
            rec=dict(run_x=a,run_y=b,pattern_x=labels[a],pattern_y=labels[b],same_pattern=labels[a]==labels[b],
                     pair_group=labels[a] if labels[a]==labels[b] else ' / '.join(sorted([labels[a],labels[b]])),
                     variant=variant,n=n,max_wco_all=float(np.nanmax(peak)),
                     weighted_zero_lag_wco=float(np.nansum(weight*co[:,n-1])/weight[valid].sum()))
            for cut,active in powered.items():
                w=np.where(active&valid,weight,0)
                key=f'weighted_scale_peak_power{int(cut*100):02d}'
                rec[key]=float(np.nansum(w*peak)/w.sum()) if w.sum()>0 else float('nan')
            rows.append(rec)
        if k%25==0:print(k+1,'/',len(tasks),f'{time.perf_counter()-start:.1f}s',flush=True)
    data=pd.DataFrame(rows)
    data.to_csv(ROOT/'results'/f'preprocessing_{args.tag}_pairs.csv',index=False)
    meta=dict(variants=args.variants,scale_count=args.scales,elapsed_seconds=time.perf_counter()-start,
              total_pairs=len(tasks),same_pairs=sum(labels[a]==labels[b] for a,b in tasks),
              wavelengths_px=[50,3000],estimator='biased',wavelet='Schuck Morlet omega0=6',
              invariant='Alternative WPS/WCS/WCO formulas unchanged in every candidate',
              domain='Todos os lags, sem COI e sem máscara conservadora',
              main_score='Média ponderada dos máximos em todos os lags da alternative WCO por escala; pesos pela média geométrica das alternative WPS do baseline no lag zero; ambas >=5% de seus próprios máximos',
              score_sensitivity='Identical calculation with 1% and 10% baseline power floors',
              limitation='Scale-wise maximization is a descriptive similarity score, not evidence of common events or causality. Correlated pairs are not independent samples.')
    (ROOT/'results'/f'preprocessing_{args.tag}_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(json.dumps(meta,indent=2),flush=True)

if __name__=='__main__':main()
