"""Numerical controls independent of the experimental flow labels."""
import json
import time
from pathlib import Path
import numpy as np
import pywt
from scipy.signal import fftconvolve
from wavelet_core import *

ROOT = Path(__file__).resolve().parents[1]

def main():
    start = time.perf_counter()
    rng = np.random.default_rng(20261023)
    checks = {}
    x, y = rng.normal(size=(2, 128))
    for estimator in ('biased', 'unbiased'):
        r = correlation(x,y,estimator)
        expected = []
        for m in range(-127,128):
            terms = [x[k+m]*y[k] for k in range(128) if 0<=k+m<128]
            expected.append(sum(terms)/(128 if estimator=='biased' else len(terms)))
        error = np.max(np.abs(r-expected))
        assert error<1e-12
        checks[f'{estimator}_literal_sum_max_error'] = float(error)
    x = rng.normal(size=2048)
    scales = np.array([20.,50.,100.])
    wc = morlet_cwt(x,scales)
    errs = []
    for i,a in enumerate(scales):
        t = np.arange(-int(np.ceil(8*a)), int(np.ceil(8*a))+1)
        kernel = np.pi**(-.25)/np.sqrt(a)*np.exp(-.5*(t/a)**2)*np.exp(1j*6*t/a)
        expected = fftconvolve(x,kernel,mode='same')
        errs.append(float(np.max(np.abs(wc[i]-expected))))
    assert max(errs)<1e-10
    checks['morlet_independent_time_convolution_max_error'] = max(errs)
    padded = morlet_cwt(x,scales,padding_sigma=12)
    checks['padding_8_vs_12_sigma_max_error'] = float(np.max(np.abs(wc-padded)))
    assert checks['padding_8_vs_12_sigma_max_error']<1e-10
    wavelengths = np.geomspace(50,500,40)
    p = pair(x,x,wavelengths)
    checks['self_zero_lag_max_error'] = float(np.nanmax(np.abs(p['wco'][:,len(x)-1]-1)))
    checks['self_global_max'] = float(np.nanmax(p['wco']))
    assert checks['self_zero_lag_max_error']<1e-12
    assert checks['self_global_max']<=1+1e-10
    q = pair(x,rng.normal(size=len(x)),wavelengths)
    checks['independent_signals_max_wco'] = float(np.nanmax(q['wco']))
    assert checks['independent_signals_max_wco']<1+1e-10
    swapped = pair(rng.normal(size=1024),rng.normal(size=1024),wavelengths)
    # A single modulated Gaussian transient with a known 100-pixel delay.
    t = np.arange(2048)
    x = np.exp(-.5*((t-800)/80)**2)*np.cos(2*np.pi*(t-800)/100)
    y = np.r_[np.zeros(100),x[:-100]]
    p = pair(x,y,wavelengths)
    j = np.argmax(np.abs(p['wxy']).max(axis=1))
    lag = int(p['lags'][np.argmax(np.abs(p['wxy'][j]))])
    checks['known_delay_100_peak_paper_lag'] = lag
    assert lag==-100
    reverse = pair(y,x,wavelengths)
    err = float(np.max(np.abs(np.abs(reverse['wxy'][:,::-1])-np.abs(p['wxy']))))
    checks['swap_reflect_magnitude_max_error'] = err
    assert err<1e-10
    pscaled = pair(2*x,-3*y,wavelengths)
    err = float(np.nanmax(np.abs(p['wco']-pscaled['wco'])))
    checks['coherence_amplitude_sign_invariance_error'] = err
    assert err<1e-9
    # Reproduce the v5 CWT engine without invoking its plotting side effects.
    r = correlation(x,x)
    s = LEGACY_CF*np.array([50.,250.,500.])
    a = pywt.cwt(r,s,LEGACY_WAVELET,method='conv',precision=12)[0]
    b = pywt.cwt(r,s,LEGACY_WAVELET,method='fft',precision=12)[0]
    checks['legacy_conv_vs_fft_max_error'] = float(np.max(np.abs(a-b)))
    assert checks['legacy_conv_vs_fft_max_error']<1e-10
    # Harmonic ridge must follow Eq.32.8; compare a dense scale grid.
    t = np.arange(8192)
    lam = np.geomspace(80,125,160)
    w = morlet_cwt(np.cos(2*np.pi*t/100),lam/FOURIER_FACTOR)
    observed = float(lam[np.argmax(np.abs(w[:,4096])**2)])
    checks['100_px_sinusoid_ridge_px'] = observed
    assert abs(observed-100)<.5
    checks['constant_signal_rejected'] = False
    try:
        prepare(np.ones(128))
    except ValueError:
        checks['constant_signal_rejected'] = True
    assert checks['constant_signal_rejected']
    checks['elapsed_seconds'] = time.perf_counter()-start
    checks['status'] = 'PASS'
    (ROOT/'results/numerical_validation.json').write_text(json.dumps(checks,indent=2)+'\n')
    print(json.dumps(checks,indent=2))

if __name__=='__main__':
    main()
