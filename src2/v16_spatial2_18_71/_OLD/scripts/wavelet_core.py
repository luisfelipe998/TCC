"""Schuck–Bodmann single-realization wavelet analysis (lag, scale).

The recommended estimator is explicitly a biased adaptation of Eqs. 32.21–23.
The unbiased option reproduces those equations; neither is silently substituted.
"""
import numpy as np
from scipy.fft import fft, ifft, fftfreq, next_fast_len
from scipy.signal import correlate, detrend
import pywt

OMEGA0 = 6.0
FOURIER_FACTOR = 4*np.pi/(OMEGA0+np.sqrt(OMEGA0**2+2))
LEGACY_WAVELET = 'cmor2.0-0.9549'
LEGACY_CF = pywt.central_frequency(LEGACY_WAVELET)
# Wavefunction amplitude ratio: cmor B=2 versus unit-energy Morlet.
LEGACY_AMPLITUDE = 1/(np.sqrt(2)*np.pi**0.25)

def prepare(raw, linear=True):
    raw = np.asarray(raw, dtype=float)
    if raw.ndim != 1 or len(raw)<8 or not np.isfinite(raw).all():
        raise ValueError('A finite, uniformly sampled one-dimensional record is required.')
    x = 255-raw
    x = detrend(x, type='linear') if linear else x-x.mean()
    x -= x.mean()
    std = x.std()
    if std <= 1e-12*max(1., np.max(np.abs(raw))):
        raise ValueError('Constant or numerically degenerate signal.')
    return x/std

def correlation(x, y, estimator='biased'):
    """Paper convention Rxy[m]=sum x[n+m] conj(y[n])/divisor.

    For y delayed by d>0 relative to x, the peak is at m=-d.
    """
    x, y = np.asarray(x), np.asarray(y)
    if x.ndim != 1 or x.shape != y.shape:
        raise ValueError('Equal-length records required; crop explicitly before preprocessing.')
    n = len(x)
    lags = np.arange(-n+1, n)
    divisor = n if estimator == 'biased' else n-np.abs(lags)
    if estimator not in ('biased', 'unbiased'):
        raise ValueError(estimator)
    return correlate(x, y, mode='full', method='fft')/divisor

def morlet_cwt(x, scales, padding_sigma=8):
    """Unit-energy Morlet CWT via a real nonnegative Fourier multiplier.

    H_a(w)=sqrt(2*pi*a)*pi^(-1/4)*exp(-(a*w-6)^2/2).
    This is the sampled Fourier form of Schuck Eq.32.6, with the usual
    non-admissible Morlet approximation (DC ~exp(-18)); no hard frequency cut.
    Padding separates periodic images by at least 16 maximum-scale sigmas.
    The scale and pixel sampling are both expressed in samples (dx=1 px).
    """
    x, scales = np.asarray(x), np.asarray(scales, dtype=float)
    if x.ndim != 1 or not np.isfinite(x).all() or np.any(scales<=0):
        raise ValueError('Invalid CWT input.')
    nfft = next_fast_len(len(x)+2*int(np.ceil(padding_sigma*scales.max())))
    omega = 2*np.pi*fftfreq(nfft)
    xf = fft(x, nfft)
    out = np.empty((len(scales), len(x)), dtype=complex)
    for j, a in enumerate(scales):
        h = np.sqrt(2*np.pi*a)*np.pi**(-0.25)*np.exp(-0.5*(a*omega-OMEGA0)**2)
        out[j] = ifft(xf*h)[:len(x)]
    return out

def transform_correlation(r, wavelengths, engine='morlet', precision=12):
    wavelengths = np.asarray(wavelengths)
    if engine == 'morlet':
        scales = wavelengths/FOURIER_FACTOR
        return morlet_cwt(r, scales), scales
    if engine == 'legacy':
        scales = LEGACY_CF*wavelengths
        return pywt.cwt(r, scales, LEGACY_WAVELET, method='fft', precision=precision)[0], scales
    raise ValueError(engine)

def safe_coherence(wxy, wxx, wyy, center, relative_floor=1e-10):
    dx, dy = np.abs(wxx[:, center]), np.abs(wyy[:, center])
    valid = (dx > relative_floor*dx.max()) & (dy > relative_floor*dy.max())
    den = dx*dy
    out = np.full(wxy.shape, np.nan)
    np.divide(np.abs(wxy)**2, den[:, None], out=out, where=valid[:, None])
    # Never clip: unbiased estimates and discrete kernels may exceed unity.
    return out, valid

def pair(x, y, wavelengths, estimator='biased', engine='morlet', precision=12):
    n = len(x)
    rxx, ryy, rxy = [correlation(a,b,estimator) for a,b in [(x,x),(y,y),(x,y)]]
    if engine == 'legacy':
        # Reproduce v5's actual Rxy convention before applying its discrete CWT.
        # Reflect the reported lag afterwards when comparing with the paper.
        rxy = correlation(y,x,estimator)
    wxx, scales = transform_correlation(rxx, wavelengths, engine, precision)
    wyy, _ = transform_correlation(ryy, wavelengths, engine, precision)
    wxy, _ = transform_correlation(rxy, wavelengths, engine, precision)
    wco, denominator_valid = safe_coherence(wxy,wxx,wyy,n-1)
    return dict(wxx=wxx,wyy=wyy,wxy=wxy,wco=wco,scales=scales,
                lags=np.arange(-n+1,n),denominator_valid=denominator_valid)
