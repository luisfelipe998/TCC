"""Valida o motor equivalente e fixa o custo do modelo nulo antes da triagem."""
import time,platform,importlib.metadata
from datetime import datetime,timezone
from core import *
from scipy.signal import fftconvolve


def main():
    for name in ['results','validation','figures']:(ROOT/name).mkdir(exist_ok=True)
    runs,labels,inv,raw=load_inputs()
    source_paths=[SRC/'estatistica5/todos_contra_todos/results/descritores_por_ensaio.npz',SRC/'estatistica4/results/inventory.csv',
                  V16/'scripts/wavelet_core.py',V16/'results/analysis_config.json',V16/'results/input_audit.csv',
                  V16/'results/preprocessing_confirm_pairs.csv',V16/'results/profiles.csv',V16/'EQUATION_AUDIT.md']
    source_paths += [V16/'intermediate'/f'{run}_input.npz' for run in runs]
    previous=json.loads((SRC/'estatistica6/results/manifest_started.json').read_text())
    for key,want in previous['sources'].items():assert sha(SRC/key)==want,key
    sources={str(p.relative_to(SRC)):sha(p) for p in source_paths}
    json_write(ROOT/'results/manifest_started.json',dict(started_utc=datetime.now(timezone.utc).isoformat(),
               config_sha256=sha(ROOT/'configs/wco_v1.json'),sources=sources,
               python=platform.python_version(),packages={k:importlib.metadata.version(k) for k in ['numpy','scipy','pandas','matplotlib','PyWavelets','Pillow']}))
    errors=[];symmetry=[];times=[]
    for n in [3171,5071,7001]:
        x=prepare(raw[runs[0]][:n]);y=prepare(raw[runs[1]][:n])
        wave=Wavelet(n)
        r=correlation(x,y)
        fast=wave.transform(r)
        original=morlet_cwt(r,LAM/FOURIER_FACTOR)
        error=float(np.max(abs(fast-original)));assert error<1e-10
        errors.append(dict(n=n,max_abs_error=error))
        px,cx=wave.auto(x);py,cy=wave.auto(y)
        co,valid=wave.cross(x,y,px,py)
        reverse,_=wave.cross(y,x,py,px)
        sy=float(np.nanmax(abs(co-reverse[:,::-1])));assert sy<1e-10;symmetry.append(sy)
        np.testing.assert_allclose(cx[:,n-1],1,atol=1e-12)
        # Soma/convolução temporal independente em três escalas extremas.
        for k in [0,99,199]:
            lag=np.arange(-len(r)+1,len(r));a=LAM[k]/FOURIER_FACTOR;t=lag/a
            kernel=np.pi**(-.25)/np.sqrt(a)*np.exp(-.5*t*t)*np.exp(6j*t)
            expected=fftconvolve(r,kernel,mode='same')
            np.testing.assert_allclose(fast[k],expected,atol=1e-11)
        if n==3171:
            start=time.perf_counter()
            for _ in range(5):
                cc,_=wave.cross(x,y,px,py);map_summary(cc,px,py,n,quantiles=False)
            pair_seconds=(time.perf_counter()-start)/5
            start=time.perf_counter()
            for _ in range(5):wave.auto(x)
            auto_seconds=(time.perf_counter()-start)/5
    predicted499=499*(186*pair_seconds+25*auto_seconds)
    B=499 if predicted499<=900 else 199
    json_write(ROOT/'configs/runtime.json',dict(null_replicates=B,pair_seconds=pair_seconds,auto_seconds=auto_seconds,
               projected_499_seconds=predicted499,decision='somente custo; nenhum escore experimental usado',configured_utc=datetime.now(timezone.utc).isoformat()))
    # Controles de agrupamento e limites; não são resultados experimentais.
    d=squareform(pdist(np.array([[0],[.1],[5],[5.1]])))
    assert metric_summary(d,np.array(['a','a','b','b']))['macro_silhouette']>.95
    assert metric_summary(d,np.array(['a','b','a','b']))['macro_silhouette']<0
    np.testing.assert_allclose(ari(np.array(['a','a','b','b']),np.array([0,1,0,1])),-.5)
    np.testing.assert_allclose(holm([.01,.04,.03]),[.03,.06,.06])
    json_write(ROOT/'validation/engine.json',dict(status='PASS',original_fft_errors=errors,
               maximum_pair_reversal_error=max(symmetry),temporal_convolution_checked=True,
               zero_lag_self_identity_checked=True,synthetic_silhouette_ari_holm_checked=True))
    print(json.dumps(dict(engine_errors=errors,pair_seconds=pair_seconds,auto_seconds=auto_seconds,null_replicates=B,projected_null_seconds=B*(186*pair_seconds+25*auto_seconds)),indent=2),flush=True)


if __name__=='__main__':main()
