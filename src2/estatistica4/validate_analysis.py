"""Controles de fórmulas; dados sintéticos nunca entram nos resultados experimentais."""
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import numpy as np
from scipy.signal import fftconvolve
from scipy.stats import f_oneway, pearsonr

from statistics_core import (FOURIER_FACTOR, anosim, binary_arrangements, distances,
                             exact_upper_p, holm, permanova, residual_descriptor,
                             unit_rows, zero_lag_profile)

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'v16_spatial2_18_71/scripts'))
import wavelet_core as wavelet


def validate():
    rng = np.random.default_rng(20261024)
    checks = {}
    z = rng.normal(size=(7, 15))
    d = distances(z)
    expected = np.array([[1-pearsonr(a,b).statistic for b in z] for a in z])
    checks['pearson_vs_scipy_max_abs_error'] = float(np.max(abs(d-expected)))
    assert checks['pearson_vs_scipy_max_abs_error'] < 1e-12
    u = unit_rows(z)
    chord_squared = np.sum((u[:,None,:]-u[None,:,:])**2, axis=2)
    assert np.allclose(chord_squared, 2*d, atol=1e-12)
    assert np.allclose(distances(7*z+4), d, atol=1e-12)
    try:
        distances(np.ones((3,10)))
    except ValueError:
        checks['constant_descriptor_rejected'] = True
    else:
        raise AssertionError('Descritor constante aceito.')
    groups = np.array([0,0,0,1,1,1,1])
    iu = np.triu_indices(len(groups),1)
    same = groups[iu[0]] == groups[iu[1]]
    intra, inter = d[iu][same], d[iu][~same]
    auc = np.mean((inter[:,None] > intra[None,:]) + .5*(inter[:,None] == intra[None,:]))
    assert abs(anosim(d,groups)-(2*auc-1)) < 1e-12
    # Todos os pares empatados devem produzir R=0, independentemente da ordem.
    tied = np.ones((7,7))-np.eye(7)
    assert anosim(tied,groups) == 0.
    perfect = np.where(groups[:,None] == groups[None,:], .1, .9)
    np.fill_diagonal(perfect,0.)
    assert abs(anosim(perfect,groups)-1.) < 1e-12
    checks['anosim_auc_ties_and_perfect_separation'] = 'PASS'
    y = rng.normal(size=7)
    stat, r2 = permanova(y[:,None],groups)
    assert np.isclose(stat, f_oneway(y[:3],y[3:]).statistic, rtol=1e-12)
    # Identidade independente de somas de quadrados calculadas das distâncias.
    total = chord_squared[iu].sum()/len(groups)
    within = 0.
    for group in np.unique(groups):
        idx = np.flatnonzero(groups == group)
        sub = chord_squared[np.ix_(idx,idx)]
        within += sub[np.triu_indices(len(idx),1)].sum()/len(idx)
    stat, r2 = permanova(u,groups)
    assert np.isclose(stat, (total-within)/(within/5), rtol=1e-12)
    assert np.isclose(r2,(total-within)/total,rtol=1e-12)
    checks['permanova_vs_univariate_anova_and_distance_identity'] = 'PASS'
    arrangements = binary_arrangements([1,0,1,0], [0,0,1,1])
    assert len(arrangements) == 4 and len(np.unique(arrangements,axis=0)) == 4
    assert np.all(arrangements[:,:2].sum(axis=1) == 1)
    assert exact_upper_p([.1,.2,.3,.4],.4) == .25
    assert np.allclose(holm([.04,.01,.03]),[.06,.03,.06])
    checks['enumeration_exact_tail_and_holm'] = 'PASS'
    lam = np.geomspace(50,3000,200)
    profiles = np.exp(rng.normal(size=(7,200)))
    residual,_,_ = residual_descriptor(profiles,lam)
    assert np.max(abs(residual.sum(axis=1))) < 1e-10
    assert np.max(abs(residual@np.log(lam))) < 1e-10
    residual2,_,_ = residual_descriptor(profiles,lam/FOURIER_FACTOR)
    assert np.allclose(residual,residual2,atol=1e-12)
    checks['log_linear_residual_orthogonality_and_scale_conversion'] = 'PASS'
    x,y = rng.normal(size=(2,128))
    for estimator in ['biased','unbiased']:
        direct = []
        for lag in range(-127,128):
            terms = [x[i+lag]*y[i] for i in range(128) if 0 <= i+lag < 128]
            direct.append(sum(terms)/(128 if estimator == 'biased' else len(terms)))
        error = float(np.max(abs(wavelet.correlation(x,y,estimator)-direct)))
        assert error < 1e-12
        checks[f'correlation_{estimator}_vs_literal_sum_max_error'] = error
    x = rng.normal(size=1024)
    r = wavelet.correlation(x,x)
    lam = np.array([50.,200.,1000.,3000.])
    w,_ = wavelet.transform_correlation(r,lam)
    reference = abs(w[:,len(x)-1])
    error = float(np.max(abs(reference-zero_lag_profile(x,lam))))
    assert error < 1e-11
    checks['zero_lag_direct_sum_vs_v16_fft_max_error'] = error
    a = 100.
    t = np.arange(-800,801)
    kernel = np.pi**(-.25)/np.sqrt(a)*np.exp(-.5*(t/a)**2)*np.exp(1j*6*t/a)
    direct = fftconvolve(x,kernel,mode='same')
    fft = wavelet.morlet_cwt(x,np.array([a]))[0]
    error = float(np.max(abs(fft-direct)))
    assert error < 1e-11
    checks['morlet_time_convolution_vs_fft_max_error'] = error
    signal = np.cos(2*np.pi*np.arange(8192)/100)
    lam = np.geomspace(90,110,200)
    w,_ = wavelet.transform_correlation(signal,lam)
    peak = float(lam[np.argmax(abs(w[:,4096]))])
    assert abs(peak-100) < .2
    checks['known_100_px_sinusoid_ridge_px'] = peak
    checks['status'] = 'PASS'
    checks['scope'] = 'Controles sintéticos de fórmulas, sem resultados experimentais.'
    return checks


if __name__ == '__main__':
    result = validate()
    (HERE/'results').mkdir(exist_ok=True)
    (HERE/'results/formula_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
