#!/usr/bin/env python3
"""Recompute selected saved real-data NCC peaks by direct masked Pearson.

This checks numerical implementation and data alignment, not physical accuracy.
"""
import hashlib
import numpy as np
import pandas as pd
from PIL import Image
import residual_motion as study

def run():
    saved=pd.read_csv(study.BASE/'results/residual_lags.csv')
    results=[]
    for run,radius,block in [('m05q10_001',1,0),('m06q10_001',3,2),('m04q80_001',3,4)]:
        files=sorted((study.SRC2/'data_processed'/run).glob('ImgA*.jpg'))
        images=np.stack([np.array(Image.open(f).convert('L'),dtype=np.float32) for f in files])
        with np.load(study.MASKS/(run+'.npz'),allow_pickle=False) as z:
            mm=[np.unpackbits(z[k],axis=-1,count=1016,bitorder='big').astype(bool) for k in ('masks_060','masks_070','masks_080')]
        gas=np.logical_and.reduce(mm);liquid=~np.logical_or.reduce(mm)
        liquid[:,:18]=False;liquid[:,71:]=False
        transformed,_=study.crossfit_predictions(images[:,18:71],study.distance_bins(mm[1])[:,18:71])
        ids=study.block_indices(501)[block]
        for phase,mask in zip(study.PHASES,(gas,liquid)):
            mask=study.erode(mask,radius)[:,18:71]
            for method in study.METHODS:
                data=transformed[method]
                for lag in (1,2,4):
                    row=saved[(saved.run==run)&(saved.erosion_radius==radius)&(saved.block==block)&
                        (saved.phase==phase)&(saved.method==method)&(saved.lag==lag)].iloc[0]
                    if not np.isfinite(row.peak):continue
                    d=int(np.rint(row['shift']))
                    xs=slice(max(0,-d),min(1016,1016-d));xt=slice(max(0,d),min(1016,1016+d))
                    use=mask[ids[:-lag],:,xs]&mask[ids[lag:],:,xt]
                    a=data[ids[:-lag],:,xs][use].astype(np.float64)
                    b=data[ids[lag:],:,xt][use].astype(np.float64)
                    direct=float(np.corrcoef(a,b)[0,1]);error=abs(direct-row.peak)
                    std=float(min(a.std(),b.std()));std_error=abs(std-row['std'])
                    relative_overlap_error=abs(len(a)-row.overlap)/max(len(a),1)
                    assert error<3e-5,(run,phase,method,lag,error)
                    assert std_error<.003,(run,phase,method,lag,std_error)
                    assert relative_overlap_error<1e-6,(run,phase,method,lag,relative_overlap_error)
                    results.append(dict(run=run,erosion_radius=radius,block=block,phase=phase,method=method,lag=lag,
                        integer_shift=d,matched_pixel_pairs=len(a),direct_correlation=direct,saved_fft_correlation=row.peak,
                        correlation_absolute_error=error,std_absolute_error=std_error,overlap_relative_error=relative_overlap_error))
        print(run,flush=True)
    study.write_csv('residual_real_numerical_check.csv',results)
    study.dump('residual_real_numerical_check.json',dict(status='passed',
        records=3,checks=len(results),maximum_pearson_absolute_error=max(r['correlation_absolute_error'] for r in results),
        maximum_texture_std_absolute_error=max(r['std_absolute_error'] for r in results),
        maximum_overlap_relative_error=max(r['overlap_relative_error'] for r in results),
        interpretation='Direct calculation from native images and masks verifies selected saved numerical peaks; no independent experimental velocity truth.',
        script_sha256=hashlib.sha256(study.Path(__file__).read_bytes()).hexdigest()))
    print('passed',len(results))

if __name__=='__main__':run()
