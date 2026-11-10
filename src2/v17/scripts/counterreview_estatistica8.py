#!/usr/bin/env python3
"""Contraprova cruzada independente; escreve somente em v17."""
from pathlib import Path
import hashlib
import importlib
import json
import os
import sys

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'src2/v17/resultados/contraprova_estatistica8'
OUT.mkdir(parents=True,exist_ok=True)
sys.dont_write_bytecode=True
os.environ['MPLCONFIGDIR']=str(OUT/'cache')
sys.path.insert(0,str(ROOT/'src2/v16.3_ancoras/scripts'))
import numpy as np
import pandas as pd
from PIL import Image
from scipy.optimize import brentq

selection=importlib.import_module('select_anchors')
inv=pd.read_csv(ROOT/'src2/estatistica8/dados/inventario.csv')
rows=[]
for r in inv.itertuples():
    if r.padrao not in ['churn','wavy annular']:continue
    with Image.open(ROOT/'src2/data_processed'/r.ensaio/'__panorama_v11.png') as im:
        raw=np.asarray(im.convert('L'))[18:72].mean(axis=0)
    chosen=[]
    for label in ['churn','wavy annular']:
        candidates,best=selection.texture_candidates(r.ensaio,label,raw)
        assert best is not None
        group={'pattern':label,'rows':candidates,'best':best}
        window=selection.common_window([group])
        best=group['best']
        chosen.append((best['ancora'],window,best['anchor_id'],len(candidates)))
    rows.append(dict(run=r.ensaio,original_pattern=r.padrao,
                     anchor_if_churn=chosen[0][0],anchor_if_annular=chosen[1][0],
                     native_window_if_churn=chosen[0][1],native_window_if_annular=chosen[1][1],
                     identical_anchor_and_window=chosen[0]==chosen[1],
                     matches_saved_center=chosen[0][0]==r.ancora_px))
pd.DataFrame(rows).to_csv(OUT/'label_swap_anchor_checks.csv',index=False)

# Independently recompute group contrasts from the new per-window widths.
w=pd.read_csv(ROOT/'src2/v17/resultados/estatistica8/larguras_janelas_distribuidas.csv')
w=w[w.pattern.isin(['churn','wavy annular'])]

def effects(x,y):
    dist=np.abs(x[:,None]-x[None,:])
    s=[]
    for i in range(len(x)):
        own=y==y[i]; own[i]=False
        aa=dist[i,own].mean();bb=dist[i,y!=y[i]].mean()
        s.append((bb-aa)/max(aa,bb) if max(aa,bb) else 0)
    s=np.array(s)
    groups=[s[y==g].mean() for g in np.unique(y)]
    c=x[y=='churn'];a=x[y=='wavy annular']
    delta=((c[:,None]>a).sum()-(c[:,None]<a).sum())/(len(c)*len(a))
    return float(np.mean(groups)),float(min(groups)),float(delta)

checks=[]
for position,g in w.groupby('position',sort=True):
    values=effects(g.width.to_numpy(),g.pattern.to_numpy())
    checks.append(dict(condition='janela_uniforme_'+str(position),macro_silhouette=values[0],
                       min_class_silhouette=values[1],cliff_delta=values[2]))
med=w.groupby(['run','pattern'],sort=True).width.median().reset_index()
values=effects(med.width.to_numpy(),med.pattern.to_numpy())
checks.append(dict(condition='mediana_5_janelas_uniformes',macro_silhouette=values[0],
                   min_class_silhouette=values[1],cliff_delta=values[2]))
pd.DataFrame(checks).to_csv(OUT/'window_statistics_recheck.csv',index=False)
reported=pd.read_csv(ROOT/'src2/v17/resultados/estatistica8/sensibilidades_descritivas.csv').set_index('condition')
maxerr=max(abs(row[k]-reported.loc[row['condition'],k])
           for row in checks for k in ['macro_silhouette','min_class_silhouette','cliff_delta'])

# Independent inversion by bracketed root of the piecewise-linear integral.
with np.load(ROOT/'src2/v17/resultados/estatistica8/wps_soma_direta.npz') as z:
    p=z['p_200'];t=np.log(z['lambda_200']);runs=z['runs']
widths=[]
for profile in p:
    dx=np.diff(t);mass=dx*(profile[:-1]+profile[1:])/2
    cumulative=np.r_[0,np.cumsum(mass)]
    qq=[]
    for q in [.25,.75]:
        target=q*cumulative[-1]
        k=int(np.searchsorted(cumulative,target,side='right')-1)
        slope=(profile[k+1]-profile[k])/dx[k]
        def objective(u):return cumulative[k]+profile[k]*u+slope*u*u/2-target
        qq.append(t[k]+brentq(objective,0,dx[k],xtol=1e-14))
    widths.append(qq[1]-qq[0])
original=pd.read_csv(ROOT/'src2/v17/resultados/estatistica8/larguras_por_ensaio.csv').set_index('run')
quantile_error=float(np.max(np.abs(np.array(widths)-original.loc[runs,'width_continuous_200'].to_numpy())))
result=dict(anchor_label_swaps_checked=len(rows),
            all_swapped_anchors_and_native_windows_identical=all(r['identical_anchor_and_window'] for r in rows),
            all_centers_match_archived_anchors=all(r['matches_saved_center'] for r in rows),
            max_independent_window_statistics_error=maxerr,
            max_independent_continuous_quantile_width_error=quantile_error,
            local_64px_overlap_fraction=(2048-64)/2048,
            interpretation='Invariância churn/anular não garante permutabilidade observacional nem valida extração entre famílias.',
            reviewed_report_sha256=hashlib.sha256((ROOT/'src2/v17/findings/auditoria_estatistica8.md').read_bytes()).hexdigest())
(OUT/'summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(result,ensure_ascii=False,indent=2))
