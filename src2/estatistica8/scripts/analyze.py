"""Separação, influência, controles operacionais e permutações por aquisição."""
from core8 import *
from math import comb
from scipy.stats import spearmanr
import time


def permutations(y, seed):
    groups, codes = np.unique(y, return_inverse=True)
    n = len(y)
    if len(groups) == 2 and comb(n, sum(codes == 0)) <= 100000:
        assignments = []
        for chosen in combinations(range(n), sum(codes == 0)):
            row = np.ones(n, int); row[list(chosen)] = 0
            assignments.append(row)
        return np.asarray(assignments), True, codes
    rng = np.random.default_rng(seed)
    return np.asarray([rng.permutation(codes) for _ in range(9999)]), False, codes


def permuted_statistics(d, w, yp):
    """All labels travel together across matrix rows/columns, never permute pairs."""
    k = int(yp.max())+1
    oh = yp[:, :, None] == np.arange(k)
    counts = oh.sum(axis=1)
    assert np.all(counts >= 2)
    means = np.einsum('ij,bjk->bik', d, oh, optimize=True)/counts[:, None, :]
    own = np.take_along_axis(means, yp[:, :, None], axis=2)[:, :, 0]
    own_n = np.take_along_axis(counts, yp, axis=1)
    a = own*own_n/(own_n-1)
    means[oh] = np.inf
    b = means.min(axis=2)
    s = np.divide(b-a, np.maximum(a, b), out=np.zeros_like(a), where=np.maximum(a, b)>0)
    sil = (s/own_n).sum(axis=1)/k
    credit = np.einsum('ij,bij->bi', w, yp[:, :, None] == yp[:, None, :], optimize=True)
    nn = (credit/own_n).sum(axis=1)/k
    return sil, nn


def tail(null, observed, exact):
    greater = np.count_nonzero(null >= observed-1e-12)
    return greater/len(null) if exact else (greater+1)/(len(null)+1)


def load(name):
    with np.load(ROOT/f'dados/perfis_{name}.npz') as z:
        return {k: z[k].copy() for k in z.files}


def main():
    start = time.monotonic()
    all_rows, classes, influence, geometries, features, individual, controls, combos = [], [], [], [], [], [], [], []
    data_primary = load('ancora')
    inv = pd.read_csv(ROOT/'dados/inventario.csv')
    m = inv.ensaio.str.extract(r'm(\d+)')[0].astype(int).to_numpy()
    q = inv.ensaio.str.extract(r'q(\d+)')[0].astype(int).to_numpy()
    np.savez_compressed(ROOT/'dados/covariaveis.npz', m=m, q=q)
    primary_matrices = primary_fitted = None
    for condition in ['ancora', 'desloc_m64', 'desloc_p64']:
        data = load(condition)
        dist, fitted, geometry, f = matrices(data)
        if condition == 'ancora': primary_matrices, primary_fitted = dist, fitted
        geometries += [dict(condition=condition, **g) for g in geometry]
        features += [dict(condition=condition, **row) for row in f]
        np.savez_compressed(ROOT/f'dados/distancias_{condition}.npz', **dist)
        for scope, patterns in SCOPES.items():
            mask = np.isin(data['labels'], patterns)
            labels, runs = data['labels'][mask], data['runs'][mask]
            base = dist['wps_pearson'][np.ix_(mask, mask)]
            for name, full in dist.items():
                d = full[np.ix_(mask, mask)]
                x = fitted[name][mask] if name in fitted else None
                w = neighbors(d, features=x)
                stat, s, credit = metrics(d, labels, w)
                row = dict(condition=condition, scope=scope, metric=name, n=len(labels), **stat)
                all_rows.append(row)
                for g in patterns:
                    classes.append(dict(condition=condition, scope=scope, metric=name, pattern=g,
                                        silhouette=float(s[labels == g].mean()), neighbor=float(credit[labels == g].mean())))
                for r, g, si, ci in zip(runs, labels, s, credit):
                    individual.append(dict(condition=condition, scope=scope, metric=name,
                                           run=r, pattern=g, silhouette=si, neighbor_credit=ci))
                if condition != 'ancora': continue
                wc = neighbors(d, features=x, base=base)
                bc, _ = balanced_credit(wc, labels)
                bb, _ = balanced_credit(neighbors(base), labels)
                combos.append(dict(scope=scope, metric=name, candidate_neighbor=stat['balanced_neighbor'],
                                   pearson_neighbor=bb, combination_neighbor=bc,
                                   gain_direct=stat['balanced_neighbor']-bb, gain_combination=bc-bb))
                for omitted in range(len(labels)):
                    keep = np.arange(len(labels)) != omitted
                    sd = scaled_distance(x[keep]) if x is not None else d[np.ix_(keep, keep)]
                    yy = labels[keep]
                    si = silhouette(sd, yy)
                    cc = [si[yy == g].mean() for g in np.unique(yy)]
                    influence.append(dict(scope=scope, metric=name, omitted_run=runs[omitted],
                                          macro_silhouette=np.mean(cc), min_class_silhouette=min(cc)))
                # Exclude all acquisitions sharing m with each query; no training reuse.
                wm = neighbors(d, features=x, blocked=m[mask])
                bm, _ = balanced_credit(wm, labels)
                controls.append(dict(scope=scope, metric=name, control='leave_m_out', value=bm, n=len(labels)))
                # Descriptive same-q control only when every class has >=2 runs.
                for value in np.unique(q[mask]):
                    sameq = q[mask] == value
                    if min(np.sum(labels[sameq] == g) for g in patterns) < 2: continue
                    qq = scaled_distance(x[sameq]) if x is not None else d[np.ix_(sameq, sameq)]
                    sx = x[sameq] if x is not None else None
                    ss, _, _ = metrics(qq, labels[sameq], neighbors(qq, features=sx))
                    controls.append(dict(scope=scope, metric=name, control='same_q_'+str(value),
                                         value=ss['balanced_neighbor'], n=int(sameq.sum()),
                                         macro_silhouette=ss['macro_silhouette'], min_class_silhouette=ss['min_class_silhouette']))
        print(f'{condition}: {len(dist)} representações, {time.monotonic()-start:.1f}s', flush=True)
    met = save_csv('metricas.csv', all_rows)
    save_csv('por_classe.csv', classes); save_csv('por_ensaio.csv', individual)
    inf = save_csv('influencia.csv', influence)
    save_csv('geometria_kernels.csv', geometries); feat = save_csv('caracteristicas.csv', features)
    save_csv('controles_operacionais.csv', controls)
    gain = pd.DataFrame(combos)
    # Test all observed candidates, not only the favorable rows. Exact binary tests
    # include all assignments (20 for plug/slug); Monte Carlo uses the +1 formula.
    tests = []
    null_cache = {}
    for sid, (scope, patterns) in enumerate(SCOPES.items()):
        mask = np.isin(data_primary['labels'], patterns)
        y = data_primary['labels'][mask]
        yp, exact, codes = permutations(y, CFG['seed']+sid)
        base = primary_matrices['wps_pearson'][np.ix_(mask, mask)]
        base_w = neighbors(base)
        base_s, base_nn = permuted_statistics(base, base_w, yp)
        observed_base = balanced_credit(base_w, codes)[0]
        for name, full in primary_matrices.items():
            d = full[np.ix_(mask, mask)]
            x = primary_fitted[name][mask] if name in primary_fitted else None
            w = neighbors(d, features=x)
            sil, nn = permuted_statistics(d, w, yp)
            observed, _, _ = metrics(d, codes, w)
            # Check the vectorized formula against independent scalar computation.
            scheck, ncheck = permuted_statistics(d, w, codes[None, :])
            np.testing.assert_allclose(scheck[0], observed['macro_silhouette'], atol=1e-12)
            np.testing.assert_allclose(ncheck[0], observed['balanced_neighbor'], atol=1e-12)
            wc = neighbors(d, features=x, base=base)
            _, nc = permuted_statistics(d, wc, yp)
            obc = balanced_credit(wc, codes)[0]
            tests.append(dict(scope=scope, metric=name, permutations=len(yp), exact=exact,
                              p_silhouette=tail(sil, observed['macro_silhouette'], exact),
                              p_neighbor=tail(nn, observed['balanced_neighbor'], exact),
                              p_gain_direct=tail(nn-base_nn, observed['balanced_neighbor']-observed_base, exact),
                              p_gain_combination=tail(nc-base_nn, obc-observed_base, exact)))
            null_cache[scope+'__'+name] = np.column_stack([sil, nn, nn-base_nn, nc-base_nn])
        print(f'Permutações {scope}: {len(yp)} alocações, exato={exact}; {time.monotonic()-start:.1f}s', flush=True)
    tf = pd.DataFrame(tests)
    tf[['p_silhouette_holm', 'p_neighbor_holm']] = holm(tf[['p_silhouette', 'p_neighbor']].to_numpy().ravel()).reshape(-1, 2)
    tf[['p_gain_direct_holm', 'p_gain_combination_holm']] = holm(tf[['p_gain_direct', 'p_gain_combination']].to_numpy().ravel()).reshape(-1, 2)
    save_csv('permutacoes.csv', tf)
    save_csv('ganho_pearson.csv', gain.merge(tf, on=['scope', 'metric'], validate='one_to_one'))
    np.savez_compressed(ROOT/'dados/distribuicoes_nulas.npz', **null_cache)
    # Stable screening; both shifts required, including geometry admissibility.
    decisions = []
    for row in met[met.condition == 'ancora'].itertuples():
        sub = met[(met.scope == row.scope) & (met.metric == row.metric)]
        ii = inf[(inf.scope == row.scope) & (inf.metric == row.metric)]
        stable = len(sub) == 3 and (sub.min_class_silhouette > 0).all() and (sub.neighbor_lift > 0).all()
        decisions.append(dict(scope=row.scope, metric=row.metric, screen_pass=row.screen_pass,
                              stable_local=stable, min_silhouette_shifts=sub.min_class_silhouette.min(),
                              min_macro_loo=ii.macro_silhouette.min(), min_class_loo=ii.min_class_silhouette.min(),
                              robust_pass=bool(row.screen_pass and stable and ii.macro_silhouette.min() > 0)))
    dec = save_csv('decisoes.csv', decisions)
    # Full panorama only for WPS; changed content and length, not a pure length test.
    with np.load(ROOT/'dados/controle_completo.npz') as z:
        dfull, ffull, xfull = spectral(z['p'], z['lam'])
    dfull['trio'] = scaled_distance(xfull)
    control_full = []
    for scope, patterns in SCOPES.items():
        mask = np.isin(data_primary['labels'], patterns); yy = data_primary['labels'][mask]
        for name, full in dfull.items():
            xx = xfull[mask] if name == 'trio' else None
            dd = full[np.ix_(mask, mask)]
            st, _, _ = metrics(dd, yy, neighbors(dd, features=xx))
            control_full.append(dict(scope=scope, metric='wps_'+name, **st))
    save_csv('sensibilidade_panorama_completo.csv', control_full)
    # Report individual-feature magnitudes without pseudo-replication.
    effects, associations = [], []
    for (condition, family), frame in feat.groupby(['condition', 'family']):
        frame = frame.set_index('run').loc[data_primary['runs']].reset_index()
        for name in ['log_peak', 'slope', 'width', 'centroid', 'entropy']:
            for scope, patterns in SCOPES.items():
                if len(patterns) != 2: continue
                a = frame.loc[frame.pattern == patterns[0], name].to_numpy()
                b = frame.loc[frame.pattern == patterns[1], name].to_numpy()
                cliff = np.sign(a[:, None]-b[None, :]).mean()
                effects.append(dict(condition=condition, family=family, feature=name, scope=scope,
                                    cliff_delta=cliff, rank_auc=(cliff+1)/2, n_a=len(a), n_b=len(b),
                                    median_a=np.median(a), median_b=np.median(b),
                                    min_a=a.min(), max_a=a.max(), min_b=b.min(), max_b=b.max()))
            if condition == 'ancora':
                for cname, cv in [('m', m), ('q', q), ('comprimento_original', inv.panorama_px.to_numpy())]:
                    associations.append(dict(family=family, feature=name, covariate=cname,
                                             spearman=spearmanr(frame[name], cv).statistic))
    save_csv('efeitos_individuais.csv', effects); save_csv('associacoes_operacionais.csv', associations)
    print('TRIAGEM\n'+dec[dec.screen_pass].to_string(index=False), flush=True)
    print('ROBUSTOS\n'+dec[dec.robust_pass].to_string(index=False), flush=True)
    json_write(ROOT/'dados/analise.json', dict(status='PASS', candidates=len(primary_matrices),
               tests_separation=2*len(tf), tests_gain=2*len(tf), elapsed_s=time.monotonic()-start,
               screen_pass=int(dec.screen_pass.sum()), robust_pass=int(dec.robust_pass.sum()),
               min_corrected_separation=float(tf[['p_silhouette_holm', 'p_neighbor_holm']].min().min()),
               min_corrected_gain=float(tf[['p_gain_direct_holm', 'p_gain_combination_holm']].min().min())))


if __name__ == '__main__': main()
