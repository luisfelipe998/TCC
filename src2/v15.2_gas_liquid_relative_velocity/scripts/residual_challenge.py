#!/usr/bin/env python3
"""Independent numerical and adversarial controls for residual motion."""
import hashlib
import json
import numpy as np
from scipy.ndimage import gaussian_filter1d
import residual_motion as study

def run():
    rng=np.random.default_rng(771924)
    n,h,w=181,24,1016;xx=np.arange(w);t=np.arange(n)
    gas=np.stack([((xx-12*i)%500)<250 for i in t])[:,None,:]
    gas=np.repeat(gas,h,axis=1)
    bins=study.distance_bins(gas)
    halo=(110+2*(bins-12)+20*np.exp(-((bins-12)/3.)**2)).astype(np.float32)
    fields=[]
    speeds=(-17.3,-7.4)
    for v in speeds:
        tex=gaussian_filter1d(rng.normal(size=(h,14000)),1.2,axis=-1)*45
        arr=np.empty_like(halo)
        for row in range(h):
            arr[:,row]=np.stack([np.interp(xx-v*i+4000,np.arange(14000),tex[row]) for i in t])
        fields.append(arr)
    im=halo+np.where(gas,fields[0],fields[1])
    prepared,_=study.crossfit_predictions(im,bins)
    results=[]
    for method in study.METHODS:
        for pi,phase in enumerate(study.PHASES):
            mask=study.erode(gas if pi==0 else ~gas,3)
            out,_,_,_=study.analyze(prepared[method],mask,im)
            error=abs(out['signed_px_frame']-speeds[pi]) if np.isfinite(out['signed_px_frame']) else np.nan
            results.append(dict(case='interface_12_and_independent_interior_textures',method=method,phase=phase,
                prescribed_texture_px_frame=speeds[pi],interface_px_frame=12,estimate_px_frame=out['signed_px_frame'],
                error_px_frame=error,accepted_blocks=out['accepted_blocks'],optical_supported=out['optical_supported']))
            if method=='interface_distance':
                assert out['optical_supported'] and error<.15,(method,phase,out,error)
    # Reverse the actual input sequence: the recovered sign must reverse.
    reverse,_=study.crossfit_predictions(im[::-1].copy(),bins[::-1].copy())
    for pi,phase in enumerate(study.PHASES):
        mask=study.erode((gas if pi==0 else ~gas)[::-1].copy(),3)
        out,_,_,_=study.analyze(reverse['interface_distance'],mask,im[::-1])
        error=abs(out['signed_px_frame']+speeds[pi])
        assert out['optical_supported'] and error<.15,(phase,out,error)
        results.append(dict(case='time_reversal',method='interface_distance',phase=phase,
            prescribed_texture_px_frame=-speeds[pi],interface_px_frame=-12,estimate_px_frame=out['signed_px_frame'],
            error_px_frame=error,accepted_blocks=out['accepted_blocks'],optical_supported=out['optical_supported']))
    # Deliberately alter evaluation intensities while keeping training images
    # and geometry unchanged: no predicted evaluation value may depend on it.
    mutated=im.copy();mutated[:n//2]+=rng.normal(0,30,size=mutated[:n//2].shape).astype(np.float32)
    altered,_=study.crossfit_predictions(mutated,bins)
    leakage_errors={}
    for method in ('fixed_background','interface_distance'):
        before=im[:n//2]-prepared[method][:n//2]
        after=mutated[:n//2]-altered[method][:n//2]
        error=float(np.max(abs(before-after)))
        leakage_errors[method]=error
        assert error<5e-5,(method,error)
    # Direct Pearson calculation is independent of the FFT accumulation.
    a=rng.uniform(40,210,size=(3,8,1016)).astype(np.float32)
    masks=rng.random(a.shape)>.3
    ch=study.math.channels(a,masks)
    sums=study.math.sufficient(ch,[0],[1])
    curve,_,_,_=study.math.curve_from_sums(sums)
    direct_errors=[]
    for d in (-199,-17,0,23,211):
        xs=slice(max(0,-d),min(1016,1016-d));xt=slice(max(0,d),min(1016,1016+d))
        use=masks[0,:,xs]&masks[1,:,xt]
        av=a[0,:,xs][use].astype(float);bv=a[1,:,xt][use].astype(float)
        direct=float(np.corrcoef(av,bv)[0,1]);estimated=float(curve[d+240])
        err=abs(direct-estimated);direct_errors.append(err)
        assert err<1e-5,(d,direct,estimated)
    extra=[]
    static=np.repeat(im[:1],n,axis=0)
    independent_frames=(120+rng.normal(0,15,size=im.shape)).astype(np.float32)
    for case,control_images in [('stationary_image_moving_mask',static),('independent_frames_moving_mask',independent_frames)]:
        transformed,_=study.crossfit_predictions(control_images,bins)
        for method in study.METHODS:
            for pi,phase in enumerate(study.PHASES):
                mask=study.erode(gas if pi==0 else ~gas,3)
                out,_,_,_=study.analyze(transformed[method],mask,control_images)
                extra.append(dict(case=case,method=method,phase=phase,accepted_blocks=out['accepted_blocks'],
                    optical_supported=out['optical_supported'],estimated_px_frame=out['signed_px_frame']))
                assert not out['optical_supported'],(case,method,phase,out)
    study.write_csv('residual_challenge.csv',results)
    study.write_csv('residual_challenge_controls.csv',extra)
    study.dump('residual_challenge.json',dict(status='passed',
        numerical_scope='Prescribed synthetic textures and held-out prediction invariance; no experimental velocity truth.',
        maximum_dynamic_interface_texture_error_px_frame=max(r['error_px_frame'] for r in results if r['method']=='interface_distance'),
        prediction_invariance_error=leakage_errors,direct_pearson_max_error=max(direct_errors),
        static_and_independent_frame_combinations_rejected=sum(not e['optical_supported'] for e in extra),
        script_sha256=hashlib.sha256(__import__('pathlib').Path(__file__).read_bytes()).hexdigest(),
        study_script_sha256=hashlib.sha256(study.Path(study.__file__).read_bytes()).hexdigest()))
    print(json.dumps(study.math.clean(results),indent=2))

if __name__=='__main__':run()
