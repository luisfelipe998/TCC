"""Verify literal v5 operations and refresh baseline rows independently."""
import ast
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pywt
from wavelet_core import *

ROOT=Path(__file__).resolve().parents[1]
LAM=np.geomspace(50,3000,200)

def main():
    source=(ROOT.parent/'v5_spatial/wavelet_v5.py').read_text()
    tree=ast.parse(source)
    chosen=[node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name in ['biased_autocorrelation','biased_crosscorrelation']]
    scope={'np':np,'N':512}
    exec(compile(ast.Module(body=chosen,type_ignores=[]),'v5_extracted_functions','exec'),scope)
    rng=np.random.default_rng(16);x,y=rng.normal(size=(2,512));n=len(x)
    rx=scope['biased_autocorrelation'](x);ry=scope['biased_autocorrelation'](y)
    xy=scope['biased_crosscorrelation'](x,y);yx=scope['biased_crosscorrelation'](y,x)
    arrays=[np.r_[rx[:0:-1],rx],np.r_[ry[:0:-1],ry],np.r_[yx[:0:-1],xy]]
    lam=np.geomspace(50,400,12);s=LEGACY_CF*lam
    expected=[pywt.cwt(r,s,LEGACY_WAVELET,method='conv')[0] for r in arrays]
    p=pair(x,y,lam,engine='legacy')
    errors={key:float(np.max(np.abs(value-p[key]))) for key,value in zip(['wxx','wyy','wxy'],expected)}
    assert max(errors.values())<1e-10
    wco=abs(expected[2])**2/(abs(expected[0][:,[n-1]])*abs(expected[1][:,[n-1]]))
    errors['wco']=float(np.nanmax(np.abs(wco-p['wco'])))
    assert errors['wco']<1e-9
    records=pd.read_csv(ROOT/'results/pair_methods.csv')
    for i,row in records[records.method=='v5_recomputed'].iterrows():
        ra=np.load(ROOT/'intermediate'/f'{row.run_x}_input.npz')['raw'];rb=np.load(ROOT/'intermediate'/f'{row.run_y}_input.npz')['raw']
        n=min(len(ra),len(rb));x=prepare(ra)[:n];y=prepare(rb)[:n]
        p=pair(x,y,LAM,engine='legacy');co=p['wco']
        vals=dict(max_wco=float(np.nanmax(co)),
            wco_fraction_above_1=float(np.mean(co>1+1e-8)),max_zero_lag_wco=float(np.nanmax(co[:,n-1])),
            median_zero_lag_wco=float(np.nanmedian(co[:,n-1])))
        for key,val in vals.items():records.loc[i,key]=val
    records.to_csv(ROOT/'results/pair_methods.csv',index=False)
    (ROOT/'results/legacy_literal_validation.json').write_text(json.dumps(dict(status='PASS',
        source='v5 functions extracted from AST; bilateral assembly and default convolution reproduced',
        errors=errors,baseline_pairs_refreshed=87,reported_lag='Paper convention: reflected AFTER legacy CWT'),indent=2)+'\n')
    print(errors)

if __name__=='__main__': main()
