"""Estatísticas de ensaio e extração independente do coeficiente wavelet em zero."""
from itertools import combinations, product

import numpy as np
from scipy.signal import correlate
from scipy.stats import rankdata

FOURIER_FACTOR = 4 * np.pi / (6 + np.sqrt(38))
WAVELENGTHS = np.geomspace(50, 3000, 200)


def zero_lag_profile(x, wavelengths=WAVELENGTHS):
    """Soma temporal da Eq. 32.4 na autocorrelação enviesada, no lag zero.

    O kernel completo preserva os lags existentes; não há máscara/COI. Esta
    função calcula apenas o coeficiente necessário, sem criar mapas descartados.
    """
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or not np.isfinite(x).all() or len(x) < 8:
        raise ValueError('Sinal inválido.')
    n = len(x)
    lag = np.arange(-n + 1, n, dtype=float)
    r = correlate(x, x, mode='full', method='fft') / n
    result = []
    for a in np.asarray(wavelengths) / FOURIER_FACTOR:
        t = lag / a
        kernel = np.pi ** (-0.25) / np.sqrt(a) * np.exp(-0.5 * t*t)
        # Para autocorrelação real e par, a parte imaginária se cancela.
        result.append(abs(np.sum(r * kernel * np.cos(6*t))))
    return np.asarray(result)


def residual_descriptor(profiles, wavelengths=WAVELENGTHS):
    p = np.asarray(profiles, dtype=float)
    if p.ndim != 2 or not np.isfinite(p).all() or np.any(p < 0):
        raise ValueError('Perfis inválidos.')
    logp = np.log(p + 1e-12)
    design = np.column_stack([np.ones(len(wavelengths)), np.log(wavelengths)])
    coefficients = np.linalg.lstsq(design, logp.T, rcond=None)[0]
    residuals = logp - (design @ coefficients).T
    return residuals, coefficients.T, logp


def unit_rows(z):
    z = np.asarray(z, dtype=float)
    if z.ndim != 2 or z.shape[1] < 2 or not np.isfinite(z).all():
        raise ValueError('Descritores inválidos.')
    centered = z - z.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(centered, axis=1, keepdims=True)
    if np.any(norms <= 1e-12 * np.maximum(1., np.linalg.norm(z, axis=1, keepdims=True))):
        raise ValueError('Pearson indefinido para descritor constante ou degenerado.')
    return centered / norms


def distances(z, spearman=False):
    if spearman:
        z = rankdata(z, axis=1, method='average')
    u = unit_rows(z)
    d = 1 - u @ u.T
    # Corrigir somente arredondamento de produto de vetores unitários.
    if d.min() < -1e-12 or d.max() > 2 + 1e-12:
        raise ArithmeticError('Dissimilaridade fora do domínio.')
    d = np.clip(d, 0, 2)
    np.fill_diagonal(d, 0.)
    return d


def anosim(d, groups):
    groups = np.asarray(groups)
    iu = np.triu_indices(len(groups), 1)
    same = groups[iu[0]] == groups[iu[1]]
    if not same.any() or same.all():
        raise ValueError('ANOSIM precisa de pares intra e inter.')
    ranks = rankdata(d[iu], method='average')
    return float((ranks[~same].mean() - ranks[same].mean()) / (len(ranks)/2))


def binary_arrangements(groups, strata):
    groups = np.asarray(groups, dtype=bool)
    strata = np.asarray(strata)
    options = []
    for value in np.unique(strata):
        members = np.flatnonzero(strata == value)
        options.append(list(combinations(members, int(groups[members].sum()))))
    arrangements = []
    for choices in product(*options):
        assignment = np.zeros(len(groups), dtype=bool)
        for choice in choices:
            assignment[list(choice)] = True
        arrangements.append(assignment)
    return np.asarray(arrangements)


def exact_upper_p(statistics, observed):
    return float(np.mean(np.asarray(statistics) >= observed - 1e-12))


def permanova(u, groups):
    u, groups = np.asarray(u), np.asarray(groups)
    n, k = len(groups), len(np.unique(groups))
    total = float(np.sum((u-u.mean(axis=0))**2))
    within = sum(float(np.sum((u[groups == g]-u[groups == g].mean(axis=0))**2))
                 for g in np.unique(groups))
    between = total - within
    if total <= 0 or within <= 0 or not 1 < k < n:
        raise ValueError('PERMANOVA degenerada.')
    return float((between/(k-1)) / (within/(n-k))), float(between/total)


def holm(pvalues):
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    adjusted = np.minimum(1., np.maximum.accumulate(p[order] * np.arange(len(p), 0, -1)))
    out = np.empty_like(p)
    out[order] = adjusted
    return out
