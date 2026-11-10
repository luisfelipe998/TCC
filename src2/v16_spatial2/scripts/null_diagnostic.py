"""Conditional random-phase null diagnostic for Schuck WCO.

This is not the Torrence/Grinsted local-coherence significance test. It tests
the maximum over the entire lag-scale plane under a stationary random-phase
model. Such a model is not asserted to describe the stitched boiling images.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
from scipy.signal import detrend
from wavelet_core import *

ROOT=Path(__file__).resolve().parents[1]
# Same 200-scale grid as the primary analysis.
LAM=np.geomspace(50,3000,200)

def scramble(x,rng):
    z=np.fft.rfft(x)
    phases=rng.uniform(0,2*np.pi,len(z))
    phases[0]=0
    if len(x)%2==0: phases[-1]=0
    y=np.fft.irfft(z*np.exp(1j*phases),n=len(x))
    # Same removal of fitted line and normalization as the observed records.
    y=detrend(y);return y/y.std()

def statistic(x,y):
    p=pair(x,y,LAM)
    use=np.isfinite(p['wco'])
    return float(p['wco'][use].max())

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--surrogates',type=int,default=199);args=ap.parse_args()
    inv=pd.read_csv(ROOT/'results/input_audit.csv');inv=inv[inv.included]
    tasks=[]
    for pattern,g in inv.groupby('pattern'):
        names=sorted(g.run)
        if len(names)>=2: tasks.append((pattern,names[0],names[1]))
    tasks.append(('cross-pattern control','m04q01_001','m05q01_001'))
    rng=np.random.default_rng(32172026); rows=[]; samples=[]; start=time.perf_counter()
    for pattern,a,b in tasks:
        ra=np.load(ROOT/'intermediate'/f'{a}_input.npz')['raw'];rb=np.load(ROOT/'intermediate'/f'{b}_input.npz')['raw']
        n=min(len(ra),len(rb));x=prepare(ra[:n]);y=prepare(rb[:n]);observed=statistic(x,y)
        null=[]
        for i in range(args.surrogates):
            value=statistic(scramble(x,rng),scramble(y,rng));null.append(value)
            samples.append(dict(run_x=a,run_y=b,replicate=i,maximum_wco=value))
        count=int(np.count_nonzero(np.array(null)>=observed))
        rows.append(dict(pattern=pattern,run_x=a,run_y=b,n=n,surrogates=args.surrogates,
            observed_max=observed,null_median=float(np.median(null)),null_q95=float(np.quantile(null,.95)),
            exceedances=count,conditional_monte_carlo_p=(1+count)/(1+args.surrogates)))
        print(a,b,rows[-1],f'elapsed {time.perf_counter()-start:.1f}s',flush=True)
    pd.DataFrame(rows).to_csv(ROOT/'results/null_diagnostic.csv',index=False)
    pd.DataFrame(samples).to_csv(ROOT/'results/null_replicates.csv',index=False)
    (ROOT/'results/null_metadata.json').write_text(json.dumps(dict(seed=32172026,
        scales=200,wavelength_min=50,wavelength_max=3000,
        domain='Todos os lags e as 200 escalas, sem COI e sem máscara conservadora; excluídos apenas denominadores numericamente degenerados',
        statistic='maximum WCO across all included lags and scales',
        model='Independent random phases preserving each record circular Fourier amplitudes before linear detrending',
        limitation='Conditional diagnostic only. Does not reproduce non-Gaussian morphology, nonstationarity or stitching. No physical linkage or causality claim.',
        between_pair_multiple_testing='No rejection decision; p values are unadjusted across the five prespecified examples.',
        elapsed_seconds=time.perf_counter()-start),indent=2)+'\n')

if __name__=='__main__':main()
