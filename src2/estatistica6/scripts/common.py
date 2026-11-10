"""Funções numéricas e leitura da estatistica6, sem escrita na origem."""
from pathlib import Path
from itertools import combinations
import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, cut_tree

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT.parent
V16 = SRC / 'v16_spatial2_18_71'
S5 = SRC / 'estatistica5/todos_contra_todos/results'
S4 = SRC / 'estatistica4/results'
sys.dont_write_bytecode = True
sys.path.insert(0, str(SRC / 'estatistica4'))
sys.path.insert(0, str(V16 / 'scripts'))
from statistics_core import residual_descriptor, unit_rows, zero_lag_profile
from wavelet_core import prepare, correlation, morlet_cwt, FOURIER_FACTOR

CONFIG = json.loads((ROOT / 'configs/triagem_v1.json').read_text())
PATTERNS = ['churn', 'wavy annular', 'plug', 'slug', 'dispersed bubbles']
COLORS = dict(zip(PATTERNS, ['#2878a0', '#e4a138', '#7460a7', '#31977b', '#69737c']))


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, obj):
    def native(x):
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.generic): return x.item()
        if isinstance(x, Path): return str(x)
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2, default=native, allow_nan=False) + '\n')


def save_csv(name, rows):
    frame = rows if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows)
    frame.to_csv(ROOT / 'results' / name, index=False, float_format='%.17g')
    return frame


def read_csv(path):
    return pd.read_csv(path, float_precision='round_trip')


def bools(values):
    out = pd.Series(values).astype(str).str.lower().map({'true': True, 'false': False})
    if out.isna().any(): raise ValueError('Booleano inválido')
    return out.to_numpy(dtype=bool)


def load_inputs():
    with np.load(S5 / 'descritores_por_ensaio.npz', allow_pickle=False) as z:
        source = {k: z[k].copy() for k in z.files}
    runs = source['runs'].astype(str)
    labels = source['patterns'].astype(str)
    lam = source['wavelengths_px']
    inventory = read_csv(S4 / 'inventory.csv').set_index('run').loc[runs].reset_index()
    assert np.array_equal(inventory.pattern, labels)
    assert bools(inventory.included).all()
    windows = read_csv(V16 / 'results/equal_length_windows.csv')
    conditions = {}
    for name in ['prefix', 'center', 'suffix']:
        table = windows[windows.window == name].pivot(index='run', columns='wavelength_px', values='zero_lag_wps').loc[runs]
        np.testing.assert_allclose(table.columns, lam, rtol=1e-14)
        conditions[name] = dict(runs=runs, labels=labels, profiles=table.to_numpy(), lam=lam)
    conditions['full'] = dict(runs=runs, labels=labels, profiles=source['profile_full'], lam=lam)
    with np.load(S4 / 'descriptors_principal_7001_inicio.npz', allow_pickle=False) as z:
        order = np.array([list(z['runs']).index(r) for r in runs if r in z['runs']])
        conditions['7001_prefix'] = dict(runs=z['runs'][order].astype(str), labels=z['labels'][order].astype(str), profiles=z['profiles'][order], lam=lam)
    return source, inventory, conditions


def log_weights(lam):
    t = np.log(np.asarray(lam, dtype=float))
    assert np.all(np.diff(t) > 0)
    h = np.empty_like(t)
    h[0] = (t[1] - t[0]) / 2
    h[-1] = (t[-1] - t[-2]) / 2
    h[1:-1] = (t[2:] - t[:-2]) / 2
    return t, h


def probabilities(profiles, lam):
    p = np.asarray(profiles, dtype=float)
    if p.ndim != 2 or not np.isfinite(p).all() or np.any(p < 0):
        raise ValueError('Perfis inválidos')
    t, h = log_weights(lam)
    mass = p * h
    total = mass.sum(axis=1, keepdims=True)
    if np.any(total <= 0): raise ValueError('Massa nula')
    return mass / total, t


def entropy(p):
    terms = np.zeros_like(p)
    np.log(p, out=terms, where=p > 0)
    return -(p * terms).sum(axis=-1) / np.log(p.shape[-1])


def features(profiles, lam):
    p, t = probabilities(profiles, lam)
    residual, coefficient, _ = residual_descriptor(profiles, lam)
    cdf = p.cumsum(axis=1)
    q25 = np.array([np.searchsorted(c, .25) for c in cdf])
    q75 = np.array([np.searchsorted(c, .75) for c in cdf])
    peak = np.argmax(profiles, axis=1)
    return pd.DataFrame(dict(
        peak_px=lam[peak], log_peak=t[peak], slope=coefficient[:, 1],
        width=t[q75] - t[q25], centroid=p @ t, entropy=entropy(p),
        band_short=p[:, lam < 200].sum(axis=1),
        band_middle=p[:, (lam >= 200) & (lam < 800)].sum(axis=1),
        band_long=p[:, lam >= 800].sum(axis=1),
        peak_ties=np.sum(profiles == np.max(profiles, axis=1, keepdims=True), axis=1),
        peak_at_edge=(peak == 0) | (peak == len(lam) - 1),
        local_maxima=np.sum((profiles[:, 1:-1] > profiles[:, :-2]) & (profiles[:, 1:-1] >= profiles[:, 2:]), axis=1)
    ))


def fit_scaler(frame):
    x = frame[['log_peak', 'slope', 'width']].to_numpy()
    center = np.median(x, axis=0)
    scale = np.quantile(x, .75, axis=0) - np.quantile(x, .25, axis=0)
    if np.any(scale <= 1e-12): raise ValueError('IQR degenerado no trio')
    return center, scale


def distances(profiles, lam, scaler=None, metrics=None):
    f = features(profiles, lam)
    p, t = probabilities(profiles, lam)
    residual, _, _ = residual_descriptor(profiles, lam)
    u = unit_rows(residual)
    cdf = p.cumsum(axis=1)
    values = {
        'pearson_residual': squareform(pdist(u)),
        'hellinger': squareform(pdist(np.sqrt(p))) / np.sqrt(2),
        'wasserstein': np.sum(np.abs(cdf[:, None, :-1] - cdf[None, :, :-1]) * np.diff(t), axis=2)
    }
    if scaler is None: scaler = fit_scaler(f)
    x = f[['log_peak', 'slope', 'width']].to_numpy()
    values['trio'] = squareform(pdist((x - scaler[0]) / scaler[1]))
    for name in ['log_peak', 'slope', 'width']:
        values[name] = squareform(pdist(f[[name]]))
    return {k: v for k, v in values.items() if metrics is None or k in metrics}, f


def silhouette_by_run(d, labels):
    labels = np.asarray(labels)
    out = np.full(len(labels), np.nan)
    for i, label in enumerate(labels):
        same = (labels == label)
        same[i] = False
        if not same.any(): continue
        a = d[i, same].mean()
        b = min(d[i, labels == other].mean() for other in np.unique(labels) if other != label)
        out[i] = (b - a) / max(a, b) if max(a, b) else 0
    return out


def nearest_agreement(d, labels):
    labels = np.asarray(labels)
    work = d.copy()
    np.fill_diagonal(work, np.inf)
    tied = np.isclose(work, work.min(axis=1, keepdims=True), atol=1e-12, rtol=1e-12)
    return (tied * (labels[:, None] == labels[None, :])).sum(axis=1) / tied.sum(axis=1)


def metric_summary(d, labels):
    labels = np.asarray(labels)
    iu = np.triu_indices(len(d), 1)
    same = labels[iu[0]] == labels[iu[1]]
    intra, inter = d[iu][same], d[iu][~same]
    s = silhouette_by_run(d, labels)
    nn = nearest_agreement(d, labels)
    groups = np.unique(labels)
    by_class = [np.mean(s[labels == g]) for g in groups if np.sum(labels == g) > 1]
    balanced_nn = np.mean([nn[labels == g].mean() for g in groups])
    chance = np.mean([(np.sum(labels == g) - 1) / (len(labels) - 1) for g in groups])
    return dict(n_runs=len(d), n_intra=len(intra), n_inter=len(inter),
                median_intra=np.median(intra), median_inter=np.median(inter),
                distance_gap=np.median(inter) - np.median(intra),
                macro_silhouette=np.nanmean(by_class), min_class_silhouette=np.nanmin(by_class),
                positive_silhouette_fraction=np.nanmean(s > 0),
                balanced_neighbor=balanced_nn, neighbor_chance=chance,
                neighbor_lift=balanced_nn - chance)


def passes(row):
    gate = CONFIG['prefix_gate']
    return bool(row['macro_silhouette'] + 1e-12 >= gate['min_macro_silhouette'] and
                row['min_class_silhouette'] > gate['min_class_silhouette_strict'] and
                row['neighbor_lift'] + 1e-12 >= gate['min_nearest_neighbor_lift'])


def scope_mask(labels, scope):
    return np.isin(labels, CONFIG['scopes'][scope])


def adjusted_rand(a, b):
    a, b = np.asarray(a), np.asarray(b)
    _, ia = np.unique(a, return_inverse=True)
    _, ib = np.unique(b, return_inverse=True)
    counts = np.zeros((ia.max()+1, ib.max()+1), dtype=int)
    np.add.at(counts, (ia, ib), 1)
    choose = lambda x: np.sum(x * (x - 1) / 2)
    total = len(a) * (len(a)-1) / 2
    observed = choose(counts)
    row, col = choose(counts.sum(axis=1)), choose(counts.sum(axis=0))
    expected = row * col / total
    denom = (row + col) / 2 - expected
    return float((observed - expected) / denom) if denom else 1.0


def clustering(d, labels):
    tree = linkage(squareform(d, checks=True), method='average', optimal_ordering=True)
    assigned = cut_tree(tree, n_clusters=len(np.unique(labels))).ravel()
    return assigned, adjusted_rand(labels, assigned), tree


def cliff(x, y):
    return float(np.sign(np.asarray(x)[:, None] - np.asarray(y)[None, :]).mean())


def weighted_median(values, weights):
    values, weights = np.asarray(values), np.asarray(weights)
    if not np.isfinite(values).all() or not np.isfinite(weights).all() or np.any(weights <= 0):
        raise ValueError('Valores/pesos inválidos')
    order = np.argsort(values, kind='stable')
    x, w = values[order], weights[order]
    c = w.cumsum() / w.sum()
    k = int(np.searchsorted(c, .5 - 1e-12))
    if abs(c[k] - .5) <= 1e-12 and k+1 < len(x): return float((x[k] + x[k+1]) / 2)
    return float(x[k])
