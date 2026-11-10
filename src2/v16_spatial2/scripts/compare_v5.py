"""Recompute v5 on all 325 pairs using exactly the preprocessing-study summaries."""
import argparse
from itertools import combinations
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
from wavelet_core import *

ROOT=Path(__file__).resolve().parents[1]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--limit',type=int);ap.add_argument('--tag',default='v5');args=ap.parse_args()
    lam=np.geomspace(50,3000,200);scales=lam/FOURIER_FACTOR
    inv=pd.read_csv(ROOT/'results/input_audit.csv');inv=inv[inv.included].sort_values('run')
    raw={r:np.load(ROOT/'intermediate'/f'{r}_input.npz')['raw'] for r in inv.run}
    full={r:prepare(x) for r,x in raw.items()};labels=inv.set_index('run').pattern.to_dict()
    cache={};rows=[];start=time.perf_counter()
    def get(run,n):
        key=(run,n)
        if key not in cache:
            x=full[run][:n]
            zero=abs(transform_correlation(correlation(x,x),lam,'legacy')[0][:,n-1])
            reference=prepare(raw[run][:n])
            basezero=abs(morlet_cwt(correlation(reference,reference),scales)[:,n-1])
            cache[key]=(x,zero,basezero)
        return cache[key]
    tasks=list(combinations(inv.run,2))
    if args.limit:tasks=tasks[:args.limit]
    for k,(a,b) in enumerate(tasks):
        n=min(len(raw[a]),len(raw[b]));x,zx,bx=get(a,n);y,zy,by=get(b,n)
        weight=np.sqrt(bx*by)
        wxy=transform_correlation(correlation(y,x),lam,'legacy')[0]
        valid=(zx>1e-10*zx.max())&(zy>1e-10*zy.max())
        co=np.divide(abs(wxy)**2,(zx*zy)[:,None],out=np.full(wxy.shape,np.nan),where=valid[:,None])
        peak=np.full(len(scales),np.nan)
        peak[valid]=np.max(co[valid],axis=1)
        rec=dict(run_x=a,run_y=b,pattern_x=labels[a],pattern_y=labels[b],same_pattern=labels[a]==labels[b],
            pair_group=labels[a] if labels[a]==labels[b] else ' / '.join(sorted([labels[a],labels[b]])),
            variant='v5',n=n,max_wco_all=float(np.nanmax(peak)),
            weighted_zero_lag_wco=float(np.nansum(weight*co[:,n-1])/weight[valid].sum()))
        for cut in [.01,.05,.1]:
            active=(bx>=cut*bx.max())&(by>=cut*by.max())
            w=np.where(active&valid,weight,0)
            rec[f'weighted_scale_peak_power{int(cut*100):02d}']=float(np.nansum(w*peak)/w.sum()) if w.sum()>0 else float('nan')
        rows.append(rec)
        if k%25==0:print(k+1,'/',len(tasks),f'{time.perf_counter()-start:.1f}s',flush=True)
    pd.DataFrame(rows).to_csv(ROOT/'results'/f'preprocessing_{args.tag}_pairs.csv',index=False)
    (ROOT/'results'/f'preprocessing_{args.tag}_metadata.json').write_text(json.dumps(dict(
        pairs=len(tasks),scales=200,elapsed_seconds=time.perf_counter()-start,
        baseline='Literal v5 estimator, legacy cmor scale conversion and preprocessing-before-crop; no coherence clipping',
        domain='Todos os lags, sem COI e sem máscara conservadora',
        comparison='Mesmos pesos de alternative WPS do baseline Morlet e domínio integral de preprocessing_confirm_pairs.csv',
        coordinate='Original v5 lag sign; all summaries are at lag zero or symmetric lag searches'),indent=2)+'\n')
    print('elapsed',time.perf_counter()-start,flush=True)

if __name__=='__main__':main()
