"""Resumo explícito dos mapas; não altera as equações do núcleo."""
from pathlib import Path
import json
import numpy as np
from PIL import Image
from wavelet_core import pair, prepare, FOURIER_FACTOR

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT.parent/'v16_spatial2_18_71'
DATA=ROOT.parent/'data_processed'
CFG=json.loads((ROOT/'config.json').read_text())
LAM=np.geomspace(*CFG['comprimentos_onda_px'],CFG['numero_escalas'])

def crop_prepare(raw,start,n=2048,linear=True):
    assert isinstance(start,(int,np.integer)) and start>=0 and start+n<=len(raw)
    return prepare(raw[start:start+n],linear=linear)

def evaluate(x,y,engine='morlet',wavelengths=LAM):
    m=pair(x,y,wavelengths,engine=engine)
    if engine=='legacy':
        m['wxy']=m['wxy'][:,::-1].copy();m['wco']=m['wco'][:,::-1].copy()
    return m,summarize(m,len(x))

def summarize(m,n):
    co=m['wco'];lags=m['lags'];a=m['scales'];valid=m['denominator_valid']
    domain=np.abs(lags)<=CFG['lag_resumo_px']
    geom=(n-1-np.abs(lags[domain])[None,:])>=np.sqrt(2)*a[:,None]
    # A mesma faixa nominal de lambdas é usada no v5; sua escala difere.
    assert geom.all(), 'O domínio principal não tem o suporte geométrico declarado.'
    values=co[:,domain];center=n-1
    assert lags[center]==0
    scale_means=np.mean(values[valid],axis=1)
    cons=(np.abs(lags[domain])[None,:]+3*a[:,None]<=n/2)&valid[:,None]
    cnt=cons.sum(axis=1);ok=cnt>0
    cons_scale=np.zeros(len(a));np.divide(np.where(cons,values,0).sum(axis=1),cnt,out=cons_scale,where=ok)
    den=np.abs(m['wxx'][:,center])*np.abs(m['wyy'][:,center])
    recomputed=np.abs(m['wxy'][valid])**2/den[valid,None]
    error=float(np.max(np.abs(co[valid]-recomputed)))
    return dict(escore=float(scale_means.mean()) if len(scale_means) else np.nan,
        media_lag_zero=float(np.mean(co[valid,center])) if valid.any() else np.nan,
        maximo_integral=float(np.nanmax(co)),minimo_integral=float(np.nanmin(co)),
        escore_conservador=float(np.mean(cons_scale[ok])) if ok.any() else np.nan,
        escalas_validas=int(valid.sum()),escalas_invalidas=int((~valid).sum()),
        escalas_conservadoras=int(ok.sum()),lags_principais=int(domain.sum()),
        suporte_geometrico_fracao=float(geom.mean()),
        fracao_integral_acima_1=float(np.mean(co[valid]>1+1e-8)),
        denominador_minimo=float(den[valid].min()) if valid.any() else np.nan,
        erro_equacao=error)

def load_signals():
    import pandas as pd
    inv=pd.read_csv(BASE/'results/input_audit.csv');inv=inv[inv.included].sort_values('run')
    signals={}
    for run in inv.run:
        pano=np.asarray(Image.open(DATA/run/'__panorama_v11.png').convert('L'))
        signals[run]={'principal':pano[18:72].mean(axis=0),'faixa30_49':pano[30:50].mean(axis=0)}
    return inv,signals

def save_map(path,m,x,y,**meta):
    np.savez_compressed(path,**m,x=x,y=y,wavelengths=LAM,
                        meta_json=np.array(json.dumps(meta,ensure_ascii=False)))
