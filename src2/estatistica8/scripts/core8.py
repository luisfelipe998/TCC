"""Matemática e leitura somente de fontes; todas as saídas ficam em estatistica8."""
from pathlib import Path
from itertools import combinations
import hashlib
import json
import os
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT.parent
SOURCE = SRC / 'v16.3_ancoras'
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / 'cache_matplotlib'))
for name in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS']:
    os.environ.setdefault(name, '1')
import numpy as np
import pandas as pd
from scipy.fft import fft, ifft, next_fast_len
from scipy.signal import correlate, detrend
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, cut_tree

CFG = json.loads((ROOT / 'protocolo.json').read_text())
PATTERNS = ['churn', 'wavy annular', 'plug', 'slug']
SCOPES = {'quatro_padroes': PATTERNS,
          'churn_annular': PATTERNS[:2], 'plug_slug': PATTERNS[2:],
          'churn_plug': ['churn', 'plug'], 'churn_slug': ['churn', 'slug'],
          'annular_plug': ['wavy annular', 'plug'], 'annular_slug': ['wavy annular', 'slug']}
CROSS = ['global_max', 'zero_weighted', 'scale_peak_weighted', 'scale_peak_power05',
         'common_lag_peak', 'map_mean', 'map_mean_weighted', 'lag_q90_weighted',
         'area05_weighted', 'mean_short', 'mean_middle', 'mean_long',
         'common_peak_short', 'common_peak_middle', 'common_peak_long', 'central_mean']
N = 2048
FACTOR = 4*np.pi / (6 + np.sqrt(38))


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def json_write(path, obj):
    def native(x):
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.generic): return x.item()
        if isinstance(x, Path): return str(x)
        raise TypeError(type(x).__name__)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2,
                                   default=native, allow_nan=False)+'\n')


def save_csv(name, rows):
    f = rows if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows)
    (ROOT/'dados').mkdir(exist_ok=True)
    f.to_csv(ROOT/'dados'/name, index=False, float_format='%.17g')
    return f


def prepare(raw):
    x = detrend(255-np.asarray(raw, float), type='linear')
    x -= x.mean()
    assert x.std() > 1e-12
    return x/x.std()


def weights(lam):
    t = np.log(lam)
    return np.r_[np.diff(t)[0]/2, (t[2:]-t[:-2])/2, np.diff(t)[-1]/2]


def bands(lam):
    return {'short': lam < 200, 'middle': (lam >= 200) & (lam < 800), 'long': lam >= 800}


class Wavelet:
    """Convolução linear exata no suporte requerido, conferida contra mapas da fonte."""
    def __init__(self, lam, n=N):
        self.lam = lam
        self.n = n
        self.m = 2*n-1
        self.nfft = next_fast_len(2*self.m-1)
        lag = np.arange(-self.m+1, self.m)
        scale = lam / FACTOR
        tau = lag[None, :]/scale[:, None]
        kernel = np.pi**(-.25)/np.sqrt(scale[:, None])*np.exp(-.5*tau*tau+6j*tau)
        wrapped = np.zeros((len(lam), self.nfft), complex)
        wrapped[:, lag % self.nfft] = kernel
        self.kernel = fft(wrapped, axis=1, workers=1)
        lags = np.arange(-n+1, n)
        self.support = (n-1-np.abs(lags)[None, :]) >= np.sqrt(2)*scale[:, None]

    def transform(self, x, y):
        r = correlate(x, y, mode='full', method='fft')/self.n
        return ifft(self.kernel*fft(r, self.nfft)[None, :], axis=1, workers=1)[:, :self.m]

    def auto(self, x):
        w = self.transform(x, x)
        p = abs(w[:, self.n-1])
        assert np.all(p > 1e-10*p.max())
        co = abs(w)**2/p[:, None]**2
        mp = np.sum(abs(w)*self.support, axis=1)/self.support.sum(axis=1)
        return p, co, mp

    def cross(self, x, y, px, py):
        return abs(self.transform(x, y))**2/(px*py)[:, None]


def summary(co, px, py, lam):
    assert np.isfinite(co).all() and co.min() >= 0 and co.max() < 1+1e-8
    h = weights(lam)
    u = h/h.sum()
    w = h*np.sqrt(px*py)
    w /= w.sum()
    center = co.shape[1]//2
    mean, peak, zero = co.mean(axis=1), co.max(axis=1), co[:, center]
    top = np.unravel_index(np.argmax(co), co.shape)
    power_weight = np.sqrt(px*py)*(px >= .05*px.max())*(py >= .05*py.max())
    assert power_weight.sum() > 0
    result = dict(global_max=co[top], zero_weighted=w@zero,
                  scale_peak_weighted=w@peak, scale_peak_power05=power_weight@peak/power_weight.sum(),
                  common_lag_peak=(w@co).max(), map_mean=u@mean,
                  map_mean_weighted=w@mean, lag_q90_weighted=w@np.quantile(co, .9, axis=1),
                  area05_weighted=w@(co >= .5).mean(axis=1),
                  central_mean=co[:, center-32:center+33].mean(),
                  peak_lambda_px=lam[top[0]], peak_lag_px=top[1]-center)
    for name, mask in bands(lam).items():
        result['mean_'+name] = (h[mask]@mean[mask])/h[mask].sum()
        wb = h[mask]*np.sqrt(px[mask]*py[mask])
        result['common_peak_'+name] = np.max(wb@co[mask])/wb.sum()
    return result, mean, peak, zero


def self_features(co, lam):
    positive = co[:, co.shape[1]//2:]
    cross = positive < .5
    assert np.all(cross.any(axis=1)), 'Meia largura censurada'
    return co.mean(axis=1), np.argmax(cross, axis=1)/lam


def spectral(profiles, lam):
    h = weights(lam)
    t = np.log(lam)
    p = profiles*h
    p /= p.sum(axis=1, keepdims=True)
    cdf = p.cumsum(axis=1)
    q25, q75 = [(cdf >= v).argmax(axis=1) for v in [.25, .75]]
    design = np.column_stack([np.ones(len(t)), t])
    coef = np.linalg.lstsq(design, np.log(profiles+1e-12).T, rcond=None)[0].T
    residual = np.log(profiles+1e-12)-coef@design.T
    unit = residual/np.linalg.norm(residual, axis=1, keepdims=True)
    peak = profiles.argmax(axis=1)
    features = dict(log_peak=t[peak], peak_px=lam[peak], slope=coef[:, 1],
                    width=t[q75]-t[q25], centroid=p@t,
                    entropy=-np.sum(p*np.log(np.maximum(p, 1e-300)), axis=1)/np.log(len(t)),
                    peak_at_edge=(peak == 0) | (peak == len(lam)-1))
    d = {'pearson': squareform(pdist(unit)),
         'hellinger': squareform(pdist(np.sqrt(p)))/np.sqrt(2),
         'wasserstein': np.sum(abs(cdf[:, None, :-1]-cdf[None, :, :-1])*np.diff(t), axis=2)}
    for name in ['log_peak', 'slope', 'width']:
        d[name] = squareform(pdist(features[name][:, None]))
    trio = np.column_stack([features[name] for name in ['log_peak', 'slope', 'width']])
    return d, features, trio


def scaled_distance(x, train=None):
    if train is None: train = np.arange(len(x))
    scale = np.quantile(x[train], .75, axis=0)-np.quantile(x[train], .25, axis=0)
    # Fixed fallback for degenerate training IQR; no dependence on labels/test values.
    scale = np.where(scale > 1e-12, scale, np.std(x[train], axis=0))
    scale = np.where(scale > 1e-12, scale, 1.)
    return squareform(pdist(x/scale))


def matrices(data):
    lam, p = data['lam'], data['p']
    n = len(p)
    out, fitted, geometry = {}, {}, []
    for k, name in enumerate(CROSS):
        d = np.zeros((n, n))
        i, j = data['pairs'].T
        d[i, j] = d[j, i] = 1-data['scores'][:, k]
        out['wco_'+name] = d
    u = weights(lam); u /= u.sum()
    for name in ['area', 'width']:
        out['self_'+name] = squareform(pdist(np.log(data['self_'+name])*np.sqrt(u)))
    bf = np.column_stack([np.log(data['self_area'][:, b])@(u[b]/u[b].sum()) for b in bands(lam).values()])
    out['self_bands'] = scaled_distance(bf)
    fitted['self_bands'] = bf
    for name, band in {'full': np.ones(len(lam), bool), **bands(lam)}.items():
        w = u[band]/u[band].sum()
        k = np.diag(data['self_area'][:, band]@w)
        i, j = data['pairs'].T
        k[i, j] = k[j, i] = data['mean'][:, band]@w
        norm = k/np.sqrt(np.outer(k.diagonal(), k.diagonal()))
        ev = np.linalg.eigvalsh(k).min()
        valid = bool(ev >= -1e-10 and norm.max() <= 1+1e-10)
        key = 'kernel_area_'+name
        geometry.append(dict(metric=key, min_eigenvalue=ev, max_normalized=norm.max(), valid=valid))
        if valid:
            d = np.sqrt(np.maximum(0, 2*(1-norm)))
            np.fill_diagonal(d, 0)
            out[key] = d
    k = np.eye(n)
    i, j = data['pairs'].T
    k[i, j] = k[j, i] = data['zero']@u
    ev = np.linalg.eigvalsh(k).min()
    valid = bool(ev >= -1e-10 and k.max() <= 1+1e-10)
    geometry.append(dict(metric='kernel_zero', min_eigenvalue=ev, max_normalized=k.max(), valid=valid))
    if valid: out['kernel_zero'] = np.sqrt(np.maximum(0, 2*(1-k)))
    feature_rows = []
    for family, power in [('wps', p), ('mean_power', data['mp'])]:
        dist, features, trio = spectral(power, lam)
        for name, d in dist.items(): out[family+'_'+name] = d
        out[family+'_trio'] = scaled_distance(trio)
        fitted[family+'_trio'] = trio
        for i in range(n):
            feature_rows.append(dict(family=family, run=str(data['runs'][i]), pattern=str(data['labels'][i]),
                                     **{k: v[i] for k, v in features.items()}))
    for d in out.values():
        assert np.isfinite(d).all() and np.allclose(d, d.T) and np.allclose(d.diagonal(), 0)
        assert d.min() >= -1e-10
        d[d < 0] = 0  # Only floating point roundoff, never alter WCO maps.
    return out, fitted, geometry, feature_rows


def silhouette(d, labels):
    groups = np.unique(labels)
    a = np.zeros(len(labels)); b = np.full(len(labels), np.inf)
    for group in groups:
        ix = labels == group
        a[ix] = d[np.ix_(ix, ix)].sum(axis=1)/(ix.sum()-1)
        b[~ix] = np.minimum(b[~ix], d[np.ix_(~ix, ix)].mean(axis=1))
    return np.divide(b-a, np.maximum(a, b), out=np.zeros_like(a), where=np.maximum(a, b)>0)


def neighbors(d, features=None, blocked=None, base=None):
    """LOO neighbor weights. Scalers and combination medians use training only."""
    n = len(d)
    w = np.zeros_like(d)
    for i in range(n):
        train = np.arange(n) != i
        if blocked is not None: train &= blocked != blocked[i]
        assert train.any()
        use = scaled_distance(features, np.flatnonzero(train)) if features is not None else d
        row = use[i].copy()
        if base is not None:
            iu = np.triu_indices(train.sum(), 1)
            scale = np.median(use[np.ix_(train, train)][iu])
            bscale = np.median(base[np.ix_(train, train)][iu])
            assert scale > 1e-12 and bscale > 1e-12
            row = .5*row/scale + .5*base[i]/bscale
        row[~train] = np.inf
        tied = np.isclose(row, row.min(), rtol=1e-12, atol=1e-12)
        w[i] = tied/tied.sum()
    return w


def balanced_credit(w, labels):
    each = np.sum(w*(labels[:, None] == labels[None, :]), axis=1)
    return np.mean([each[labels == g].mean() for g in np.unique(labels)]), each


def ari(labels, clusters):
    i, j = np.triu_indices(len(labels), 1)
    a, b = labels[i] == labels[j], clusters[i] == clusters[j]
    expected = a.sum()*b.sum()/len(a)
    den = (a.sum()+b.sum())/2-expected
    return float(((a & b).sum()-expected)/den) if den else 1.


def metrics(d, labels, w=None):
    s = silhouette(d, labels)
    cs = [s[labels == g].mean() for g in np.unique(labels)]
    w = neighbors(d) if w is None else w
    nn, credits = balanced_credit(w, labels)
    chance = np.mean([(sum(labels == g)-1)/(len(labels)-1) for g in np.unique(labels)])
    cl = cut_tree(linkage(squareform(d), method='average'), n_clusters=len(cs)).ravel()
    macro, minimum = np.mean(cs), min(cs)
    return dict(macro_silhouette=macro, min_class_silhouette=minimum, balanced_neighbor=nn,
                chance_neighbor=chance, neighbor_lift=nn-chance, ari=ari(labels, cl),
                screen_pass=bool(macro >= .1-1e-12 and minimum > 0 and nn-chance >= .1-1e-12)), s, credits


def holm(p):
    p = np.asarray(p)
    order = np.argsort(p)
    corrected = np.minimum(1, np.maximum.accumulate(p[order]*np.arange(len(p), 0, -1)))
    out = np.empty_like(p); out[order] = corrected
    return out
