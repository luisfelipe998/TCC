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
from extraction_config import (ROW_START,ROW_END,ROW_COUNT,ROW_SLICE,FIXED_ROW,
    FIXED_BAND_START,FIXED_BAND_END,FIXED_BAND_SLICE,TRIM_FRACTION,TRIM_COUNT,TRIM_SLICE)

ROOT=Path(__file__).resolve().parents[1]
VARIANTS=['baseline','demean_only','center_row','central_band20','core_row_median',
          'core_trimmed_mean20','spatial_gaussian2','spatial_median5','rank_uniform','winsor_1pct']

def extract(run):
    image=np.asarray(Image.open(ROOT.parent/'data_processed'/run/'__panorama_v11.png').convert('L'),dtype=float)
    assert image.ndim==2 and image.shape[0]>ROW_END
    core=image[ROW_SLICE];raw=core.mean(axis=0)
    sorted_rows=np.sort(core,axis=0)
    return dict(baseline=raw,demean_only=raw,center_row=image[FIXED_ROW],central_band20=image[FIXED_BAND_SLICE].mean(axis=0),
                core_row_median=np.median(core,axis=0),core_trimmed_mean20=sorted_rows[TRIM_SLICE].mean(axis=0),
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
    reference_inputs={}
    for run in inv.run:
        with np.load(ROOT/'intermediate'/f'{run}_input.npz') as old:
            reference_inputs[run]=old['raw_20_59'].copy()
    cache={};reference_cache={};rows=[];start=time.perf_counter()
    def get_reference_zero(run,n):
        key=(run,n)
        if key not in reference_cache:
            x=prepare(reference_inputs[run][:n])
            reference_cache[key]=np.abs(morlet_cwt(correlation(x,x),scales)[:,n-1])
        return reference_cache[key]
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
        rx,ry=get_reference_zero(a,n),get_reference_zero(b,n)
        reference_weight=np.sqrt(rx*ry)
        # Baseline power thresholds are fixed across candidates in each pair.
        powered={cut:(bx>=cut*bx.max())&(by>=cut*by.max()) for cut in [.01,.05,.1]}
        reference_powered={cut:(rx>=cut*rx.max())&(ry>=cut*ry.max()) for cut in [.01,.05,.1]}
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
                     denominator_valid_scales=int(valid.sum()),
                     eligible_scales_power05=int((powered[.05]&valid).sum()),
                     fixed20_59_eligible_scales_power05=int((reference_powered[.05]&valid).sum()),
                     weighted_zero_lag_wco=float(np.nansum(weight*co[:,n-1])/weight[valid].sum()),
                     fixed20_59_weighted_zero_lag_wco=float(np.nansum(reference_weight*co[:,n-1])/reference_weight[valid].sum()))
            for cut,active in powered.items():
                w=np.where(active&valid,weight,0)
                key=f'weighted_scale_peak_power{int(cut*100):02d}'
                rec[key]=float(np.nansum(w*peak)/w.sum()) if w.sum()>0 else float('nan')
                rw=np.where(reference_powered[cut]&valid,reference_weight,0)
                rec['fixed20_59_'+key]=float(np.nansum(rw*peak)/rw.sum()) if rw.sum()>0 else float('nan')
            rows.append(rec)
        if k%25==0:print(k+1,'/',len(tasks),f'{time.perf_counter()-start:.1f}s',flush=True)
    data=pd.DataFrame(rows)
    data.to_csv(ROOT/'results'/f'preprocessing_{args.tag}_pairs.csv',index=False)
    meta=dict(variants=args.variants,scale_count=args.scales,elapsed_seconds=time.perf_counter()-start,
              total_pairs=len(tasks),same_pairs=sum(labels[a]==labels[b] for a,b in tasks),
              wavelengths_px=[50,3000],estimator='biased',wavelet='Schuck Morlet omega0=6',
              extraction_rows_inclusive=[ROW_START,ROW_END],extraction_row_count=ROW_COUNT,
              fixed_row_control=FIXED_ROW,fixed_band_control_inclusive=[FIXED_BAND_START,FIXED_BAND_END],
              trimmed_fraction_requested_each_end=TRIM_FRACTION,trimmed_rows_each_end=TRIM_COUNT,
              trimmed_fraction_actual_each_end=TRIM_COUNT/ROW_COUNT,trimmed_retained_rows=ROW_COUNT-2*TRIM_COUNT,
              invariant='Alternative WPS/WCS/WCO formulas unchanged in every candidate',
              domain='Todos os lags, sem COI e sem máscara conservadora',
              main_score='Média ponderada dos máximos em todos os lags da alternative WCO por escala; pesos pela média geométrica das alternative WPS do baseline no lag zero; ambas >=5% de seus próprios máximos',
              score_sensitivity='Identical calculation with 1% and 10% baseline power floors',
              fixed20_59_scores='Mesmas alternative WCO desta faixa; pesos e pisos fixos calculados com a referência corrigida 20–59 do cache anterior, cortada antes do preparo. Não misturar com escores nativos 18–71.',
              denominator_valid_scales='Contagem de escalas que passaram o piso numérico de 1e-10 para ambos os denominadores do par e da variante.',
              eligible_scales_power05='Contagem após combinar denominadores válidos com ambos os pisos nativos de 5%; a coluna fixed20_59 usa os pisos antigos.',
              limitation='Scale-wise maximization is a descriptive similarity score, not evidence of common events or causality. Correlated pairs are not independent samples.')
    (ROOT/'results'/f'preprocessing_{args.tag}_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(json.dumps(meta,indent=2),flush=True)

if __name__=='__main__':main()
