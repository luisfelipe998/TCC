#!/usr/bin/env python3
"""Harder moving-support and deformed-appearance controls for phase-local NCC."""
from pathlib import Path
import importlib.util
import json
import numpy as np
from scipy.ndimage import gaussian_filter1d, binary_erosion

BASE=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('phase_ncc',BASE/'scripts/phase_local_experiment.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)

def run():
    rng=np.random.default_rng(9025)
    n,h,w=101,20,1016
    xx=np.arange(w);tt=np.arange(n)
    masks=np.stack([((xx+12*t)%510)<330 for t in tt])[:,None,:]
    masks=np.repeat(masks,h,axis=1)
    velocity=(-21.3,-5.7)
    textures=[gaussian_filter1d(rng.normal(size=(h,12000)),1.2,axis=1)*35 for _ in range(2)]
    content=[]
    for p,v in enumerate(velocity):
        content.append(np.stack([np.stack([np.interp(xx-v*t+4500,np.arange(12000),textures[p][row]) for row in range(h)]) for t in tt]))
    image=np.where(masks,140+content[0],80+content[1]).astype(np.float32)
    results=[]
    for p,m in enumerate((masks,~masks)):
        eroded=binary_erosion(m,structure=np.ones((1,7,7)))
        o,_,_,_=mod.analyze(image,eroded)
        assert o['optical_supported'] and abs(o['signed_px_frame']-velocity[p])<.15
        results.append(dict(test='moving_mask_independent_phase_textures',phase=mod.PHASES[p],interface_speed=-12.,
            true_texture_speed=velocity[p],estimate=o['signed_px_frame'],error=abs(o['signed_px_frame']-velocity[p]),
            passed=o['optical_supported'],reasons=o['reasons']))
    # Static radial shading prevents the raw phase-variance guard alone from
    # rejecting this moving binary interface; temporal nulls must still reject.
    shading=np.linspace(0,30,h)[None,:,None]
    bare=(np.where(masks,170.,45.)+shading).astype(np.float32)
    for p,m in enumerate((masks,~masks)):
        eroded=binary_erosion(m,structure=np.ones((1,7,7)))
        o,_,_,_=mod.analyze(bare,eroded)
        assert not o['optical_supported']
        results.append(dict(test='untextured_phases_with_radial_shading',phase=mod.PHASES[p],interface_speed=-12.,
            true_texture_speed=None,estimate=o['signed_px_frame'],error=None,passed=o['optical_supported'],reasons=o['reasons']))
    # An independently randomized intensity at each frame preserves changing
    # masks and marginal brightness but removes material texture correspondence.
    random=(np.where(masks,140.,80.)+rng.normal(size=(n,h,w))*15).astype(np.float32)
    for p,m in enumerate((masks,~masks)):
        eroded=binary_erosion(m,structure=np.ones((1,7,7)))
        o,_,_,_=mod.analyze(random,eroded)
        assert not o['optical_supported']
        results.append(dict(test='frame_independent_noise_with_moving_mask',phase=mod.PHASES[p],interface_speed=-12.,
            true_texture_speed=None,estimate=o['signed_px_frame'],error=None,passed=o['optical_supported'],reasons=o['reasons']))
    out=dict(results=results,limitation='These synthetic controls test observable texture recovery under moving support. Real silhouettes need not identify material phases.')
    (BASE/'results/phase_local_challenge.json').write_text(json.dumps(mod.clean(out),indent=2)+'\n')
    print(json.dumps(mod.clean(out),indent=2))

if __name__=='__main__':run()
