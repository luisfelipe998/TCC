#!/usr/bin/env python3
"""Targeted topology sensitivity for the one dispersed-bubble recording.

The case was selected because the inherited open-envelope mask creates visible
rectangular regions in clear image areas. All primary optical gates are kept.
"""
import hashlib
import numpy as np
from PIL import Image
import residual_motion as study

def run():
    run='m07q01_001';fps=7000.
    files=sorted((study.SRC2/'data_processed'/run).glob('ImgA*.jpg'))
    images=np.stack([np.array(Image.open(f).convert('L'),dtype=np.float32) for f in files])
    file=study.SRC2/'v15_gas_liquid_velocity/phase_conditioned_motion/results/closed_contour_m07q01_masks.npz'
    with np.load(file,allow_pickle=False) as z:
        assert np.array_equal(z['frame_ids'],np.arange(501))
        mm=[np.unpackbits(z[k],axis=-1,count=1016,bitorder='big').astype(bool) for k in ('masks_060','masks_070','masks_080')]
    gas=np.logical_and.reduce(mm);liquid=~np.logical_or.reduce(mm)
    liquid[:,:18]=False;liquid[:,71:]=False
    original=images[:,18:71]
    prepared,predictions=study.crossfit_predictions(original,study.distance_bins(mm[1])[:,18:71])
    summaries=[];blocks=[];lags=[];nulls=[]
    for radius in (3,1):
        for phase,mask in zip(study.PHASES,(gas,liquid)):
            mask=study.erode(mask,radius)[:,18:71]
            for method in study.METHODS:
                out,bs,ls,ns=study.analyze(prepared[method],mask,original)
                keys=dict(run=run,pattern='dispersed bubbles',fps=fps,erosion_radius=radius,phase=phase,method=method)
                summaries.append(dict(**keys,signed_px_s=out['signed_px_frame']*fps,**out))
                blocks.extend(dict(**keys,**x) for x in bs)
                lags.extend(dict(**keys,**x) for x in ls)
                nulls.extend(dict(**keys,**x) for x in ns)
    pairs=study.pair_rows([dict(run=run,pattern='dispersed bubbles',fps=fps)],summaries,blocks)
    for name,data in [('summary',summaries),('blocks',blocks),('lags',lags),('nulls',nulls),('pairs',pairs)]:
        study.write_csv('residual_closed_m07q01_'+name+'.csv',data)
    study.dump('residual_closed_m07q01_metadata.json',dict(run=run,native_frames=501,
        targeted_scope='One dispersed-bubble run, selected by prior known topology artifact; not pooled with primary 26-run outcomes.',
        mask_source=str(file.relative_to(study.SRC2)),mask_source_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(study.Path(__file__).read_bytes()).hexdigest(),
        full_script_sha256=hashlib.sha256(study.Path(study.__file__).read_bytes()).hexdigest(),
        predictor=predictions,pairs=pairs,material_pair_identified=False))
    print([(s['erosion_radius'],s['method'],s['phase'],s['accepted_blocks']) for s in summaries])

if __name__=='__main__':run()
