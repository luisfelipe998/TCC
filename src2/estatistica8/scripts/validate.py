"""Auditoria independente dos cálculos e integridade das fontes utilizadas."""
from core8 import *
from analyze import load
import ast
import platform
import scipy
import time


def main():
    start = time.monotonic()
    data = load('ancora')
    dist, fitted, geometry, feat = matrices(data)
    metrics_frame = pd.read_csv(ROOT/'dados/metricas.csv').query('condition == "ancora"')
    # Reuse only pure reference expressions through AST, never execute old scripts.
    reference = SRC/'estatistica7/scripts/core.py'
    wanted = ['silhouette', 'neighbor', 'ari', 'holm']
    tree = ast.parse(reference.read_text())
    pure = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in wanted]
    assert len(pure) == len(wanted)
    ns = {'np': np, 'combinations': combinations}
    exec(compile(ast.Module(body=pure, type_ignores=[]), str(reference), 'exec'), ns)
    max_sil_error = 0.
    for row in metrics_frame.itertuples():
        mk = np.isin(data['labels'], SCOPES[row.scope])
        d, y = dist[row.metric][np.ix_(mk, mk)], data['labels'][mk]
        ref = ns['silhouette'](d, y)
        macro = np.mean([ref[y == g].mean() for g in np.unique(y)])
        max_sil_error = max(max_sil_error, abs(macro-row.macro_silhouette))
        np.testing.assert_allclose(macro, row.macro_silhouette, atol=1e-12)
        np.testing.assert_allclose(min(ref[y == g].mean() for g in np.unique(y)), row.min_class_silhouette, atol=1e-12)
        if row.metric not in fitted:
            credit = ns['neighbor'](d, y)
            value = np.mean([credit[y == g].mean() for g in np.unique(y)])
            np.testing.assert_allclose(value, row.balanced_neighbor, atol=1e-12)
        clusters = cut_tree(linkage(squareform(d), method='average'), n_clusters=len(np.unique(y))).ravel()
        np.testing.assert_allclose(ns['ari'](y, clusters), row.ari, atol=1e-12)
    # Independent weighted-quantile loop for each run; compare against saved widths.
    f = pd.DataFrame(feat).query('family == "wps"').set_index('run').loc[data['runs']]
    lam = data['lam']; t = np.log(lam)
    h = weights(lam)
    for i, power in enumerate(data['p']):
        mass = power*h/np.dot(power, h)
        cumulative = 0.; q25 = q75 = None
        for k, value in enumerate(mass):
            cumulative += value
            if q25 is None and cumulative >= .25: q25 = k
            if q75 is None and cumulative >= .75: q75 = k
        np.testing.assert_allclose(t[q75]-t[q25], f.width.iloc[i], atol=1e-14)
    # Spectral distance symmetry and exact identity with the archived Pearson r.
    with np.load(SOURCE/'intermediario/estatistica/descritores.npz') as z:
        np.testing.assert_allclose(data['p'], z['p'], rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(dist['wps_pearson']**2, 2*(1-z['r']), atol=1e-12)
    for name, d in dist.items():
        assert np.isfinite(d).all() and np.min(d) >= 0
        np.testing.assert_allclose(d, d.T, atol=1e-14)
        np.testing.assert_allclose(d.diagonal(), 0, atol=1e-14)
    # Check all stored permutation tails directly against saved null distributions.
    tests = pd.read_csv(ROOT/'dados/permutacoes.csv')
    gains = pd.read_csv(ROOT/'dados/ganho_pearson.csv').set_index(['scope', 'metric'])
    obs = metrics_frame.set_index(['scope', 'metric'])
    with np.load(ROOT/'dados/distribuicoes_nulas.npz') as z:
        for row in tests.itertuples():
            key = (row.scope, row.metric)
            values = z[row.scope+'__'+row.metric]
            assert values.shape == (row.permutations, 4)
            g = gains.loc[key]; o = obs.loc[key]
            thresholds = [o.macro_silhouette, o.balanced_neighbor, g.gain_direct, g.gain_combination]
            expected = []
            for column, value in enumerate(thresholds):
                count = np.count_nonzero(values[:, column] >= value-1e-12)
                expected.append(count/len(values) if row.exact else (count+1)/(len(values)+1))
            np.testing.assert_allclose(expected, [row.p_silhouette, row.p_neighbor, row.p_gain_direct, row.p_gain_combination], atol=1e-12)
        assert z['churn_annular__wps_width'].shape[0] == 92378
        assert z['plug_slug__wps_width'].shape[0] == 20
    for raw, corrected in [(['p_silhouette', 'p_neighbor'], ['p_silhouette_holm', 'p_neighbor_holm']),
                           (['p_gain_direct', 'p_gain_combination'], ['p_gain_direct_holm', 'p_gain_combination_holm'])]:
        ref = ns['holm'](tests[raw].to_numpy().ravel()).reshape(-1, 2)
        np.testing.assert_allclose(ref, tests[corrected].to_numpy(), atol=1e-12)
    # Independently reconstruct the highlighted combined-1NN result, with medians
    # calculated only in the 18 training acquisitions of each held-out run.
    mk = np.isin(data['labels'], SCOPES['churn_annular'])
    base = dist['wps_pearson'][np.ix_(mk, mk)]
    width = f.width.to_numpy()[mk]; y = data['labels'][mk]
    credits = []
    for i in range(len(y)):
        train = [j for j in range(len(y)) if i != j]
        ps = np.median([base[j, k] for j, k in combinations(train, 2)])
        ws = np.median([abs(width[j]-width[k]) for j, k in combinations(train, 2)])
        distances = np.array([.5*base[i, j]/ps+.5*abs(width[i]-width[j])/ws for j in train])
        best = np.isclose(distances, distances.min(), atol=1e-12, rtol=1e-12)
        credits.append(np.mean([y[j] == y[i] for j, tied in zip(train, best) if tied]))
    independent_combo = np.mean([np.asarray(credits)[y == group].mean() for group in np.unique(y)])
    np.testing.assert_allclose(independent_combo, gains.loc[('churn_annular', 'wps_width')].combination_neighbor, atol=1e-12)
    # Verify actual files again, not only names/mtimes. No source is written.
    sources = json.loads((ROOT/'dados/fontes_sha256.json').read_text())
    for i, (relative, expected) in enumerate(sources.items()):
        assert sha(SRC.parent/relative) == expected, relative
        if (i+1) % 150 == 0: print(f'Fontes preservadas: {i+1}/{len(sources)}', flush=True)
    import matplotlib
    json_write(ROOT/'dados/validacao_final.json', dict(status='PASS', sources_unchanged=len(sources),
        source_reference=str(reference.relative_to(SRC.parent)), source_reference_sha256=sha(reference),
        configurations_checked=252, max_silhouette_error=max_sil_error,
        combined_neighbor_independent=independent_combo,
        probability_tails_checked=1008, geometry_rejected=[g['metric'] for g in geometry if not g['valid']],
        python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__,
        matplotlib=matplotlib.__version__, elapsed_s=time.monotonic()-start))
    print('Validação independente e fontes: PASS', flush=True)


if __name__ == '__main__': main()
