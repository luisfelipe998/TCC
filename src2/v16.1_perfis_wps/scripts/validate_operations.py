"""Controles independentes das novas operações e da interpretação dos perfis."""
import json
import numpy as np
from scipy.integrate import trapezoid
from profile_core import (ROOT, LAM, wc, quadrature, distribution, hellinger,
                          cosine_distance, zero_profile, support_counts)
from evaluation import concordance


def main():
    checks = {}
    lam = np.array([5., 11., 30., 100.])
    profile = 2+3*np.log(lam)
    w = quadrature(lam)
    analytic = 2*np.log(100/5)+1.5*(np.log(100)**2-np.log(5)**2)
    checks['quadratura_erro_funcao_linear_log'] = abs(w@profile-analytic)
    assert np.isclose(w@profile, analytic, atol=1e-12)
    assert np.isclose(w@profile, trapezoid(profile, x=np.log(lam)))
    p = distribution(profile, lam)
    checks['soma_massas_erro'] = abs(p.sum()-1)
    assert hellinger(p,p) == 0
    assert np.isclose(hellinger(np.array([1.,0.]), np.array([0.,1.])), 1)
    checks['amplitude_perfil_erro'] = hellinger(p, distribution(profile*1e200, lam))
    assert checks['amplitude_perfil_erro'] < 1e-14
    rejected = 0
    for bad in [np.zeros(4), np.full(4, np.nan), np.array([1,2,-1,4]), np.full(4,np.inf)]:
        try:
            distribution(bad,lam)
        except ValueError:
            rejected += 1
    assert rejected == 4
    checks['perfis_degenerados_rejeitados'] = rejected
    assert np.isclose(cosine_distance(profile,profile*3,lam), 0, atol=1e-14)
    # Enumeração manual: (1,0), (1,2), (2,0), (2,2): 1+0+1+0,5 / 4.
    assert concordance([1,2], [0,2]) == .625
    assert concordance([1,2], [0,2], [1,3], [2,1]) == 9.5/12
    # Índices das máscaras comparados a contagens por enumeração independente.
    coi, cons = support_counts(101, [10.,40.])
    for j, la in enumerate([10.,40.]):
        a = la/wc.FOURIER_FACTOR
        assert coi[j].sum() == sum(100-abs(b) >= np.sqrt(2)*a for b in range(-100,101))
        assert cons[j].sum() == sum(abs(b)+3*a <= 101/2 for b in range(-100,101))
    rng = np.random.default_rng(20261024)
    x = rng.normal(size=1024)
    shifted = np.r_[np.zeros(100),x[:-100]]
    peak = int(np.argmax(wc.correlation(x,shifted))-len(x)+1)
    assert peak == -100
    checks['lag_atraso_100'] = peak
    # Componentes bem separadas na faixa; sem seleção de parâmetros por rótulos.
    t = np.arange(8192)
    xa, xb = np.sin(2*np.pi*t/100), np.sin(2*np.pi*t/500)
    pa, pb = zero_profile(wc.prepare(xa)), zero_profile(wc.prepare(xb))
    checks['H_componentes_100_500'] = hellinger(distribution(pa,LAM), distribution(pb,LAM))
    assert checks['H_componentes_100_500'] > .8
    checks['H_amplitude_sinal'] = hellinger(distribution(pa,LAM), distribution(zero_profile(wc.prepare(7*xa)),LAM))
    assert checks['H_amplitude_sinal'] < 1e-10
    # Espectros circulares iguais e fases distintas; correlações lineares não são idênticas.
    f = np.fft.rfftfreq(len(t))
    amplitude = np.zeros_like(f)
    band = (f >= 1/800) & (f <= 1/60)
    amplitude[band] = 1/np.sqrt(f[band])
    xx = np.fft.irfft(amplitude*np.exp(1j*rng.uniform(-np.pi,np.pi,len(f))),n=len(t))
    yy = np.fft.irfft(amplitude*np.exp(1j*rng.uniform(-np.pi,np.pi,len(f))),n=len(t))
    q = wc.pair(wc.prepare(xx),wc.prepare(yy),LAM)
    zx,zy = np.abs(q['wxx'][:,len(t)-1]),np.abs(q['wyy'][:,len(t)-1])
    px,py = distribution(zx,LAM),distribution(zy,LAM)
    weight = np.sqrt(zx*zy)
    active = (zx >= .05*zx.max()) & (zy >= .05*zy.max()) & q['denominator_valid']
    weight = np.where(active,weight,0)
    phase_score = float(np.average(np.max(q['wco'][active],axis=1),weights=weight[active]))
    checks['fases_distintas'] = dict(H=hellinger(px,py), escore_WCO=phase_score,
        WCO_zero_ponderada=float(np.average(q['wco'][active,len(t)-1],weights=weight[active])),
        erro_espectros_circulares=float(np.max(np.abs(np.abs(np.fft.rfft(xx))-np.abs(np.fft.rfft(yy))))),
        nota='Espectros circulares iguais antes do preparo; WPS linear com bordas e detrend pode diferir.')
    assert checks['fases_distintas']['H'] < .15 and phase_score < .5
    np.savez_compressed(ROOT/'results/synthetic_controls.npz', wavelengths=LAM,
                        seno100=pa,seno500=pb,phase_x=zx,phase_y=zy,
                        signal_x=xx,signal_y=yy,
                        phase_wco_scale_peak=np.max(q['wco'],axis=1))
    checks['status'] = 'PASS'
    (ROOT/'results/operations_validation.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(checks,ensure_ascii=False,indent=2))


if __name__ == '__main__':
    main()
