"""Triagem congelada antes dos cálculos: descritores, fontes e métricas."""
import importlib.metadata
from datetime import datetime, timezone
import platform
import time
from itertools import combinations

from common import *
from PIL import Image
from scipy.stats import wasserstein_distance, pearsonr


def source_paths():
    paths = [S5 / n for n in ['descritores_por_ensaio.npz', 'indice_ensaios.csv', 'matriz_pearson.csv', 'pearson_325_pares.csv', 'summary.json', 'manifest.json']]
    paths += [S4 / n for n in ['inventory.csv', 'descriptors_principal_7001_inicio.npz']]
    paths += [V16 / 'results' / n for n in ['input_audit.csv', 'analysis_config.json', 'profiles.csv', 'equal_length_windows.csv', 'preprocessing_confirm_pairs.csv', 'preprocessing_confirm_metadata.json', 'null_diagnostic.csv']]
    paths += [V16 / 'scripts/wavelet_core.py', SRC / 'estatistica4/statistics_core.py', ROOT / 'PLANO.md']
    previous = json.loads((S5 / 'manifest.json').read_text())
    paths += [SRC / k for k in previous['source_sha256']]
    return sorted(set(paths))


def begin_manifest():
    for name in ['results', 'figures', 'validation']:
        (ROOT / name).mkdir(exist_ok=True)
    paths = source_paths()
    manifest = dict(started_utc=datetime.now(timezone.utc).isoformat(), python=platform.python_version(),
                    packages={k: importlib.metadata.version(k) for k in ['numpy', 'scipy', 'pandas', 'matplotlib', 'Pillow', 'PyWavelets']},
                    config_sha256=sha(ROOT / 'configs/triagem_v1.json'),
                    sources={str(p.relative_to(SRC)): sha(p) for p in paths})
    write_json(ROOT / 'results/manifest_started.json', manifest)
    previous = json.loads((S5 / 'manifest.json').read_text())
    for k, v in previous['source_sha256'].items():
        if sha(SRC / k) != v: raise AssertionError(f'Fonte alterada desde estatistica5: {k}')
    for k in ['todos_contra_todos/results/descritores_por_ensaio.npz', 'todos_contra_todos/results/summary.json', 'todos_contra_todos/results/matriz_pearson.csv']:
        assert sha(SRC / 'estatistica5' / k) == previous['output_sha256'][k]
    return manifest


def validate(source, inventory, conditions):
    rng = np.random.default_rng(np.random.SeedSequence([CONFIG['seed'], 100]))
    rows = []
    runs, labels, lam = source['runs'], source['patterns'], source['wavelengths_px']
    assert len(runs) == len(set(runs)) == 26
    assert dict(pd.Series(labels).value_counts()) == {'churn': 10, 'wavy annular': 9, 'plug': 3, 'slug': 3, 'dispersed bubbles': 1}
    audit = read_csv(V16 / 'results/input_audit.csv')
    assert set(audit.loc[~bools(audit.included), 'run']) == set(CONFIG['excluded_runs'])
    assert set(audit.loc[bools(audit.included), 'run']) == set(runs)
    np.testing.assert_allclose(lam, np.geomspace(50, 3000, 200), rtol=1e-14)
    assert int(source['n_common_px']) == 3171
    np.testing.assert_allclose(source['profile'], conditions['prefix']['profiles'], atol=1e-14)
    raw = {}
    for run, n in zip(runs, source['n_original']):
        with np.load(V16 / 'intermediate' / f'{run}_input.npz', allow_pickle=False) as z:
            raw[run] = z['raw'].copy()
        image = np.asarray(Image.open(SRC / 'data_processed' / run / '__panorama_v11.png').convert('L'), dtype=float)
        np.testing.assert_array_equal(image[18:72].mean(axis=0), raw[run])
        assert len(raw[run]) == int(n)
    for name, condition in conditions.items():
        for i, run in enumerate(condition['runs']):
            n = 7001 if name == '7001_prefix' else (len(raw[run]) if name == 'full' else 3171)
            start = (len(raw[run])-n)//2 if name == 'center' else (len(raw[run])-n if name == 'suffix' else 0)
            actual = zero_lag_profile(prepare(raw[run][start:start+n]), lam)
            cached = condition['profiles'][i]
            rel = np.linalg.norm(actual-cached) / np.linalg.norm(cached)
            assert rel < 1e-10
            rows.append(dict(condition=name, run=run, n=n, start=start, end_exclusive=start+n, relative_profile_error=rel))
    pd.DataFrame(rows).to_csv(ROOT / 'validation/cache_profiles.csv', index=False, float_format='%.17g')
    residual, coeff, _ = residual_descriptor(source['profile'], lam)
    np.testing.assert_allclose(residual, source['residual'], atol=1e-12)
    u = unit_rows(residual)
    r = u @ u.T
    saved = read_csv(S5 / 'matriz_pearson.csv').set_index('run').loc[runs, runs].to_numpy()
    np.testing.assert_allclose(r, saved, atol=2e-14)
    pearson_error = max(abs(r[i, j]-pearsonr(residual[i], residual[j]).statistic) for i, j in combinations(range(26), 2))
    assert pearson_error < 1e-12
    iu = np.triu_indices(26, 1)
    same = labels[iu[0]] == labels[iu[1]]
    delta = np.median(r[iu][same]) - np.median(r[iu][~same])
    np.testing.assert_allclose(delta, .07555100435301754, atol=1e-13)
    assert same.sum() == 87 and (~same).sum() == 238
    p, t = probabilities(source['profile'], lam)
    np.testing.assert_allclose(p.sum(axis=1), 1)
    np.testing.assert_allclose(probabilities(13*source['profile'], lam)[0], p, atol=1e-15)
    ds, f = distances(source['profile'], lam)
    w_error = max(abs(ds['wasserstein'][i, j] - wasserstein_distance(t, t, p[i], p[j])) for i, j in combinations(range(26), 2))
    assert w_error < 1e-12
    independent_h = np.sqrt(np.maximum(0, 1-np.sqrt(p) @ np.sqrt(p).T))
    np.testing.assert_allclose(ds['hellinger'][iu], independent_h[iu], atol=1e-12)
    for d in ds.values():
        np.testing.assert_allclose(d, d.T, atol=1e-14)
        np.testing.assert_allclose(np.diag(d), 0, atol=1e-14)
        assert d.min() >= 0 and np.isfinite(d).all()
    assert ds['hellinger'].max() <= 1 and f.width.min() >= 0
    np.testing.assert_allclose(f[['band_short','band_middle','band_long']].sum(axis=1), 1)
    point = np.zeros((1,200)); point[0,10]=1
    np.testing.assert_allclose(entropy(point), 0)
    np.testing.assert_allclose(entropy(np.ones((1,200))/200), 1)
    point2 = np.zeros((2,200)); point2[0,10]=1; point2[1,100]=1
    wpoint = np.sum(abs(np.cumsum(point2[0])[:-1]-np.cumsum(point2[1])[:-1])*np.diff(t))
    np.testing.assert_allclose(wpoint, abs(t[10]-t[100]))
    np.testing.assert_allclose(np.linalg.norm(np.sqrt(point2[0])-np.sqrt(point2[1]))/np.sqrt(2), 1)
    for n in range(1, 14):
        x = rng.normal(size=n)
        np.testing.assert_allclose(weighted_median(x, np.ones(n)), np.median(x))
    try:
        unit_rows(np.ones((2,20)))
        raise AssertionError('Descritor degenerado aceito')
    except ValueError:
        pass
    # Centros separados e mistura deliberada detectam inversões da silhueta/ARI.
    d = squareform(pdist(np.array([[0.],[.1],[5.],[5.1]])))
    assert metric_summary(d, np.array(['a','a','b','b']))['macro_silhouette'] > .95
    assert metric_summary(d, np.array(['a','b','a','b']))['macro_silhouette'] < 0
    assert adjusted_rand([0,0,1,1], [7,7,2,2]) == 1
    np.testing.assert_allclose(adjusted_rand([0,0,1,1], [0,1,0,1]), -.5)
    assert passes(dict(macro_silhouette=.1, min_class_silhouette=.05, neighbor_lift=.5-.4))
    permutation = rng.permutation(26)
    before = metric_summary(ds['hellinger'], labels)
    after = metric_summary(ds['hellinger'][np.ix_(permutation, permutation)], labels[permutation])
    np.testing.assert_allclose(list(before.values()), list(after.values()))
    # Independência temporal/FFT no coeficiente central, incluindo a ordem do preparo.
    x = prepare(raw[runs[0]][:3171])
    w = morlet_cwt(correlation(x,x), lam/FOURIER_FACTOR)
    fft_error = np.max(abs(abs(w[:,3170])-source['profile'][0]))
    assert fft_error < 1e-10
    wrong = prepare(raw[runs[0]])[:3171]
    assert not np.allclose(x, wrong)
    z = dict(status='PASS', n_cache_profiles=len(rows), max_cache_relative_error=max(v['relative_profile_error'] for v in rows),
             n_extractions=26, pearson_max_error=pearson_error, wasserstein_max_error=w_error,
             fft_zero_max_error=fft_error, reference_delta=delta,
             controls=['extração PNG', 'elegibilidade', 'hashes anteriores', 'ordenação', 'pares distintos', 'corte antes do preparo',
                       'Pearson SciPy', 'Wasserstein SciPy', 'Hellinger independente', 'entropia/massas/bandas',
                       'mediana ponderada', 'vetor degenerado', 'silhueta sintética', 'ARI sintético'])
    write_json(ROOT / 'validation/numerical.json', z)
    return raw


def mean_profiles(raw, runs, lam, condition):
    out = []
    started = time.perf_counter()
    for i, run in enumerate(runs):
        length = 7001 if condition == '7001_prefix' else 3171
        start = (len(raw[run])-length)//2 if condition == 'center' else (len(raw[run])-length if condition == 'suffix' else 0)
        x = prepare(raw[run][start:start+length])
        w = morlet_cwt(correlation(x,x), lam/FOURIER_FACTOR)
        out.append(np.abs(w).mean(axis=1))
        if (i+1) % 10 == 0: print(f'mean power {condition}: {i+1}/{len(runs)}; {time.perf_counter()-started:.1f}s', flush=True)
    return np.asarray(out)


def wcs_zero(raw, source):
    runs, lam = source['runs'], source['wavelengths_px']
    x = [prepare(raw[r][:3171]) for r in runs]
    pairs = list(combinations(range(len(runs)), 2))
    rxy = np.array([correlation(x[i], x[j]) for i,j in pairs])
    lag = np.arange(-3170,3171)
    a = lam / FOURIER_FACTOR
    tau = lag[None,:] / a[:,None]
    kernel = np.pi**(-.25)/np.sqrt(a[:,None])*np.exp(-.5*tau*tau)*np.exp(-6j*tau)
    coefficients = rxy @ kernel.T
    magnitude = abs(coefficients)
    _, h = log_weights(lam)
    similarity = np.eye(len(runs))
    for k, (i,j) in enumerate(pairs):
        p, q = source['profile'][i], source['profile'][j]
        valid = (p>1e-10*p.max()) & (q>1e-10*q.max())
        similarity[i,j] = similarity[j,i] = np.sum(magnitude[k,valid]*h[valid])/np.sum(np.sqrt(p[valid]*q[valid])*h[valid])
    assert similarity.min() >= 0 and similarity.max() <= 1+1e-10
    # Verificação de dois pares contra a CWT integral original, sem aproximar a grade.
    errors = []
    for k in [0, len(pairs)-1]:
        full = morlet_cwt(rxy[k], a)[:,3170]
        errors.append(float(np.max(abs(coefficients[k]-full))))
    assert max(errors) < 1e-10
    np.savez_compressed(ROOT/'results/wcs_zero_3171.npz', runs=runs, pairs=np.array(pairs), wavelengths_px=lam, coefficients=coefficients, similarity=similarity)
    write_json(ROOT/'validation/wcs_zero.json', dict(status='PASS', independent_fft_max_error=max(errors), pairs=len(pairs), n=3171))
    return 1-similarity


def main():
    manifest = begin_manifest()
    source, inventory, conditions = load_inputs()
    raw = validate(source, inventory, conditions)
    print('Fontes e controles numéricos: PASS', flush=True)
    scaler = fit_scaler(features(source['profile'], source['wavelengths_px']))
    write_json(ROOT/'results/scaling.json', dict(center=scaler[0], iqr=scaler[1]))
    matrices, summaries, feature_rows, individual = {}, [], [], []

    def add(family, condition, metric, d, runs, labels):
        key = f'{family}__{condition}__{metric}'
        matrices[key] = d
        for scope in CONFIG['scopes']:
            mask = scope_mask(labels, scope)
            if np.sum(mask) < 4 or set(labels[mask]) != set(CONFIG['scopes'][scope]): continue
            sub = d[np.ix_(mask,mask)]
            row = dict(family=family, condition=condition, metric=metric, scope=scope, **metric_summary(sub, labels[mask]))
            row['prefix_gate'] = passes(row)
            summaries.append(row)
            for run, label, s, nn in zip(runs[mask], labels[mask], silhouette_by_run(sub,labels[mask]), nearest_agreement(sub,labels[mask])):
                individual.append(dict(family=family,condition=condition,metric=metric,scope=scope,run=run,pattern=label,silhouette=s,neighbor_agreement=nn))

    for condition, data in conditions.items():
        values, f = distances(data['profiles'], data['lam'], scaler)
        for metric, d in values.items(): add('wps', condition, metric, d, data['runs'], data['labels'])
        f.insert(0, 'run', data['runs']);f.insert(1, 'pattern', data['labels']);f.insert(2, 'condition', condition)
        feature_rows.append(f)
    # Mean power comparável: mesmo comprimento para todos os ensaios.
    mp = mean_profiles(raw, source['runs'], source['wavelengths_px'], 'prefix')
    np.savez_compressed(ROOT/'results/mean_power_prefix.npz', runs=source['runs'], labels=source['patterns'], wavelengths_px=source['wavelengths_px'], profiles=mp)
    for metric, d in distances(mp, source['wavelengths_px'], metrics=CONFIG['mean_power_candidates'])[0].items():
        add('mean_power','prefix',metric,d,source['runs'],source['patterns'])
    add('wcs_zero','prefix','normalized_magnitude',wcs_zero(raw,source),source['runs'],source['patterns'])
    # Escores arquivados WCO mantêm N variável por par; apenas triagem histórica.
    cached = read_csv(V16/'results/preprocessing_confirm_pairs.csv')
    cached = cached[cached.variant=='baseline'].copy()
    meta = json.loads((V16/'results/preprocessing_confirm_metadata.json').read_text())
    assert meta['scale_count']==200 and len(cached)==325 and bools(cached.same_pattern).sum()==87
    index = {r:i for i,r in enumerate(source['runs'])}
    assert len({tuple(sorted([r.run_x,r.run_y])) for r in cached.itertuples()})==325
    for metric in CONFIG['wco_cached_candidates']:
        similarity = np.eye(26)
        for r in cached.itertuples():
            i,j=index[r.run_x],index[r.run_y]
            assert r.n==min(len(raw[r.run_x]),len(raw[r.run_y]))
            assert r.pattern_x==source['patterns'][i] and r.pattern_y==source['patterns'][j]
            similarity[i,j]=similarity[j,i]=getattr(r,metric)
        assert np.isfinite(similarity).all() and similarity.max()<=1+1e-10
        add('wco_cached','pair_min_length',metric,1-similarity,source['runs'],source['patterns'])
    frame = save_csv('triagem.csv', summaries)
    save_csv('triagem_por_ensaio.csv', individual)
    save_csv('a3_features.csv', pd.concat(feature_rows,ignore_index=True))
    np.savez_compressed(ROOT/'results/distances_screen.npz', **matrices)
    save_csv('inventory.csv', inventory)
    decisions=[]
    for r in frame[frame.condition.isin(['prefix','pair_min_length'])].to_dict('records'):
        decisions.append(dict(family=r['family'], metric=r['metric'], scope=r['scope'],
                              action='avaliar_estabilidade' if r['prefix_gate'] else 'encerrar_na_triagem',
                              macro_silhouette=r['macro_silhouette'], min_class_silhouette=r['min_class_silhouette'], neighbor_lift=r['neighbor_lift']))
    save_csv('decisoes_preliminares.csv',decisions)
    assert all(sha(SRC/k)==v for k,v in manifest['sources'].items())
    write_json(ROOT/'results/screen_completed.json', dict(completed_utc=datetime.now(timezone.utc).isoformat(), source_files_unchanged=True, config_sha256=manifest['config_sha256'],
                                                         candidates=len(decisions), pass_count=sum(r['action']=='avaliar_estabilidade' for r in decisions)))
    print(frame[frame.condition.isin(['prefix','pair_min_length'])][['family','metric','scope','macro_silhouette','min_class_silhouette','neighbor_lift','prefix_gate']].to_string(index=False),flush=True)


if __name__=='__main__': main()
