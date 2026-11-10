"""Operações novas para perfis; as equações de Schuck permanecem na base."""
from pathlib import Path
import sys
import numpy as np
from scipy.fft import fft, fftfreq, next_fast_len

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT.parent / 'v16_spatial2_18_71'
sys.dont_write_bytecode = True
sys.path.insert(0, str(BASE / 'scripts'))
import wavelet_core as wc

LAM = np.geomspace(50., 3000., 200)


def quadrature(lam):
    lam = np.asarray(lam, float)
    if lam.ndim != 1 or len(lam) < 2 or not np.isfinite(lam).all() or np.any(lam <= 0) or np.any(np.diff(lam) <= 0):
        raise ValueError('Eixo positivo, finito e crescente obrigatório.')
    d = np.diff(np.log(lam))
    return np.r_[d[0]/2, (d[:-1] + d[1:])/2, d[-1]/2]


def distribution(profile, lam):
    p = np.asarray(profile, float)
    w = quadrature(lam)
    if p.shape != w.shape or not np.isfinite(p).all() or np.any(p < 0) or p.max() <= 0:
        raise ValueError('Perfil inválido ou área degenerada.')
    # Reescalonamento intermediário evita overflow sem alterar a distribuição.
    mass = (p / p.max()) * w
    area = mass.sum()
    if not np.isfinite(area) or area <= 0:
        raise ValueError('Área degenerada.')
    return mass / area


def hellinger(p, q):
    return float(np.sqrt(np.sum((np.sqrt(p)-np.sqrt(q))**2)/2))


def cosine_distance(p, q, lam):
    w = quadrature(lam)
    p, q = np.asarray(p, float), np.asarray(q, float)
    distribution(p, lam)
    distribution(q, lam)
    p, q = p/p.max(), q/q.max()
    return float(1 - np.dot(w*p, q)/np.sqrt(np.dot(w*p, p)*np.dot(w*q, q)))


def zero_profile(x, lam=LAM):
    """Avalia apenas a coluna central da mesma FFT/CWT do núcleo herdado.

    O preenchimento e os multiplicadores são idênticos a morlet_cwt.
    A soma inversa de Fourier substitui apenas o cálculo dos demais lags.
    """
    x = np.asarray(x)
    r = wc.correlation(x, x)
    scales = np.asarray(lam)/wc.FOURIER_FACTOR
    nfft = next_fast_len(len(r) + 2*int(np.ceil(8*scales.max())))
    omega = 2*np.pi*fftfreq(nfft)
    kernel = fft(r, nfft) * np.exp(1j*omega*(len(x)-1)) / nfft
    out = []
    for a in scales:
        h = np.sqrt(2*np.pi*a)*np.pi**(-.25)*np.exp(-.5*(a*omega-6)**2)
        out.append(abs(np.sum(kernel*h)))
    return np.asarray(out)


def legacy_zero(x, lam=LAM, same_scales=False):
    axis = np.asarray(lam)
    if same_scales:
        axis = axis/(wc.FOURIER_FACTOR*wc.LEGACY_CF)
    w, _ = wc.transform_correlation(wc.correlation(x, x), axis, 'legacy')
    return np.abs(w[:, len(x)-1])


def support_counts(n, lam=LAM):
    b = np.arange(-n+1, n)
    a = np.asarray(lam)/wc.FOURIER_FACTOR
    coi = (n-1-np.abs(b)[None, :]) >= np.sqrt(2)*a[:, None]
    conservative = np.abs(b)[None, :] + 3*a[:, None] <= n/2
    return coi, conservative
