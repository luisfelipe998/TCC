"""Independent numerical audit and descriptive sensitivities; writes only v17.

Run from repository root: PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1
    OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 src2/v17/scripts/audit_estatistica8.py
No audited module is imported. 1NN is reconstructed only to check the pre-existing
95% claim; new analyses contain no prediction or classifier.
"""
from pathlib import Path
from itertools import combinations
import hashlib
import json
import platform
import numpy as np
import pandas as pd
import scipy
from scipy.spatial.distance import pdist, squareform
from PIL import Image

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'src2/estatistica8'
OUT = ROOT / 'src2/v17/resultados/estatistica8'
OUT.mkdir(parents=True, exist_ok=True)
FACTOR = 4*np.pi/(6+np.sqrt(38))
N = 2048


def save_json(name, obj):
    (OUT/name).write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def load_npz(path):
    with np.load(path) as z:
        return {k: z[k].copy() for k in z.files}


def silhouette(d, y):
    values = []
    for i in range(len(y)):
        same = np.flatnonzero(y == y[i])
        same = same[same != i]
        a = d[i, same].mean()
        b = min(d[i, y == group].mean() for group in np.unique(y) if group != y[i])
        values.append((b-a)/max(a, b) if max(a, b) else 0.)
    values = np.asarray(values)
    means = [values[y == g].mean() for g in np.unique(y)]
    return float(np.mean(means)), float(min(means))


def trapezoid_weights(t):
    edges = np.r_[t[0], (t[1:]+t[:-1])/2, t[-1]]
    return np.diff(edges)


def widths(power, lam, continuous=False):
    t = np.log(lam)
    if not continuous:
        cdf = np.cumsum(power*trapezoid_weights(t), axis=1)
        cdf /= cdf[:, -1:]
        return t[(cdf >= .75).argmax(axis=1)]-t[(cdf >= .25).argmax(axis=1)]
    result = []
    # Exact inversion of the integral of a piecewise-linear density in log lambda.
    for p in power:
        dx = np.diff(t)
        mass = (p[1:]+p[:-1])*dx/2
        cumulative = np.r_[0., np.cumsum(mass)]
        quantiles = []
        for q in [.25, .75]:
            target = q*cumulative[-1]
            k = min(int(np.searchsorted(cumulative, target, side='right')-1), len(dx)-1)
            needed = target-cumulative[k]
            slope = (p[k+1]-p[k])/dx[k]
            u = 2*needed/(p[k]+np.sqrt(max(0., p[k]**2+2*slope*needed)))
            quantiles.append(t[k]+u)
        result.append(quantiles[1]-quantiles[0])
    return np.asarray(result)


def prepare(raw):
    t = np.arange(len(raw), dtype=float)
    a = np.column_stack([np.ones(len(raw)), t-t.mean()])
    x = 255-np.asarray(raw, dtype=float)
    x -= a@np.linalg.lstsq(a, x, rcond=None)[0]
    x -= x.mean()
    return x/np.sqrt(np.mean(x*x))


def zero_kernel(lam):
    a = lam/FACTOR
    tau = np.arange(-N+1, N)[None, :]/a[:, None]
    # Autocorrelation is even, so the imaginary part cancels exactly.
    return np.pi**(-.25)/np.sqrt(a[:, None])*np.exp(-tau*tau/2)*np.cos(6*tau)


def direct_power(raw, kernel):
    x = prepare(raw)
    return abs(kernel@(np.correlate(x, x, mode='full')/len(x)))


def effect(values, labels, condition, details=None):
    mk = np.isin(labels, ['churn', 'wavy annular'])
    v, y = values[mk], labels[mk]
    a, b = v[y == 'churn'], v[y == 'wavy annular']
    macro, minimum = silhouette(abs(v[:, None]-v), y)
    return dict(condition=condition, n=len(y), macro_silhouette=macro,
                min_class_silhouette=minimum, cliff_delta=float(np.sign(a[:, None]-b).mean()),
                median_churn=float(np.median(a)), median_annular=float(np.median(b)),
                median_difference=float(np.median(a)-np.median(b)), **(details or {}))


def main():
    data = load_npz(BASE/'dados/perfis_ancora.npz')
    lam, y, runs = data['lam'], data['labels'], data['runs']
    inv = pd.read_csv(BASE/'dados/inventario.csv')
    assert np.array_equal(inv.ensaio.to_numpy(), runs)
    checks = {}
    expected_hashes = json.loads((BASE/'dados/fontes_sha256.json').read_text())
    mismatches = []
    for relative, expected in expected_hashes.items():
        p = ROOT/relative
        if not p.exists() or hashlib.file_digest(p.open('rb'), 'sha256').hexdigest() != expected:
            mismatches.append(relative)
    checks['archived_source_hashes_checked'] = len(expected_hashes)
    checks['archived_source_hash_mismatches'] = mismatches

    scopes = {'quatro_padroes': ['churn', 'wavy annular', 'plug', 'slug'],
              'churn_annular': ['churn', 'wavy annular'], 'plug_slug': ['plug', 'slug'],
              'churn_plug': ['churn', 'plug'], 'churn_slug': ['churn', 'slug'],
              'annular_plug': ['wavy annular', 'plug'], 'annular_slug': ['wavy annular', 'slug']}
    metrics = pd.read_csv(BASE/'dados/metricas.csv')
    differences = []
    for condition in ['ancora', 'desloc_m64', 'desloc_p64']:
        distances = load_npz(BASE/f'dados/distancias_{condition}.npz')
        for row in metrics[metrics.condition == condition].itertuples():
            mk = np.isin(y, scopes[row.scope])
            macro, minimum = silhouette(distances[row.metric][np.ix_(mk, mk)], y[mk])
            differences.append(dict(condition=condition, scope=row.scope, metric=row.metric,
                                    macro_error=abs(macro-row.macro_silhouette),
                                    min_class_error=abs(minimum-row.min_class_silhouette)))
    dif = pd.DataFrame(differences)
    dif.to_csv(OUT/'reproducao_silhuetas.csv', index=False)
    checks['silhouette_configurations'] = len(dif)
    checks['max_silhouette_error'] = float(dif[['macro_error','min_class_error']].to_numpy().max())
    assert checks['max_silhouette_error'] < 1e-12

    features = pd.read_csv(BASE/'dados/caracteristicas.csv')
    f = features.query('condition == "ancora" and family == "wps"').set_index('run').loc[runs]
    original_width = widths(data['p'], lam)
    checks['max_width_error'] = float(abs(original_width-f.width.to_numpy()).max())
    assert checks['max_width_error'] < 1e-12
    dist = load_npz(BASE/'dados/distancias_ancora.npz')
    design = np.column_stack([np.ones(len(lam)), np.log(lam)])
    lp = np.log(data['p']+1e-12)
    residual = lp-(design@np.linalg.lstsq(design, lp.T, rcond=None)[0]).T
    r = np.corrcoef(residual)
    checks['max_pearson_identity_error'] = float(abs(dist['wps_pearson']**2-2*(1-r)).max())

    mk = np.isin(y, scopes['churn_annular'])
    yy = y[mk]
    dv = abs(original_width[mk, None]-original_width[None, mk])
    dp = dist['wps_pearson'][np.ix_(mk, mk)]
    credits = []
    for i in range(len(yy)):
        train = np.array([j for j in range(len(yy)) if j != i])
        pairs = list(combinations(train, 2))
        a = np.median([dv[j,k] for j,k in pairs])
        b = np.median([dp[j,k] for j,k in pairs])
        score = .5*dv[i, train]/a+.5*dp[i, train]/b
        close = np.isclose(score, min(score), rtol=1e-12, atol=1e-12)
        credits.append(float(np.mean(yy[train[close]] == yy[i])))
    checks['existing_combo_balanced_neighbor'] = float(np.mean([np.mean(np.array(credits)[yy == g]) for g in np.unique(yy)]))
    checks['existing_combo_total_credit'] = float(np.sum(credits))
    assert abs(checks['existing_combo_balanced_neighbor']-.95) < 1e-12

    # Independent enumeration of the binary allocations (not saved null tails).
    masks = np.zeros((92378, len(yy)), dtype=bool)
    for j, chosen in enumerate(combinations(range(len(yy)), sum(yy == 'churn'))):
        masks[j, list(chosen)] = True
    null = []
    for z in np.array_split(masks, 100):
        na, nb = z.sum(axis=1)[:, None], (~z).sum(axis=1)[:, None]
        suma, sumb = z@dv.T, (~z)@dv.T
        within = np.where(z, suma/(na-1), sumb/(nb-1))
        between = np.where(z, sumb/nb, suma/na)
        s = (between-within)/np.maximum(between, within)
        null.extend(((s*z).sum(axis=1)/na[:, 0]+(s*(~z)).sum(axis=1)/nb[:, 0])/2)
    null = np.asarray(null)
    obs = silhouette(dv, yy)[0]
    count = int(np.count_nonzero(null >= obs-1e-12))
    checks['exact_width_tail_count'] = count
    checks['exact_width_allocations'] = len(null)
    checks['exact_width_p'] = count/len(null)
    np.save(OUT/'nulo_exato_largura_independente.npy', null)
    tests = pd.read_csv(BASE/'dados/permutacoes.csv')
    row = tests.query('scope == "churn_annular" and metric == "wps_width"').iloc[0]
    assert abs(checks['exact_width_p']-row.p_silhouette) < 1e-12
    for raw, adjusted in [(['p_silhouette', 'p_neighbor'], ['p_silhouette_holm', 'p_neighbor_holm']),
                          (['p_gain_direct', 'p_gain_combination'], ['p_gain_direct_holm', 'p_gain_combination_holm'])]:
        p = tests[raw].to_numpy().ravel()
        order = np.argsort(p)
        adjusted_values = np.empty_like(p)
        running = 0.
        for rank, idx in enumerate(order):
            running = max(running, (len(p)-rank)*p[idx])
            adjusted_values[idx] = min(1., running)
        error = float(abs(adjusted_values-tests[adjusted].to_numpy().ravel()).max())
        checks['max_holm_error_'+raw[0]] = error
        assert error < 1e-12

    kernel = zero_kernel(lam)
    dense_lam = np.geomspace(50, 1024, 800)
    dense_kernel = zero_kernel(dense_lam)
    p_direct, p_dense, window_rows = [], [], []
    src_hashes = {}
    for i, row in enumerate(inv.itertuples()):
        path = ROOT/'src2/data_processed'/row.ensaio/'__panorama_v11.png'
        src_hashes[str(path.relative_to(ROOT))] = hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()
        with Image.open(path) as image:
            raw = np.asarray(image.convert('L'))[18:72].mean(axis=0)
        anchor = raw[row.inicio_px:row.fim_exclusivo_px]
        p_direct.append(direct_power(anchor, kernel))
        p_dense.append(direct_power(anchor, dense_kernel))
        starts = np.rint(np.linspace(0, len(raw)-N, 5)).astype(int)
        for position, start in enumerate(starts):
            power = direct_power(raw[start:start+N], kernel)[None, :]
            window_rows.append(dict(run=row.ensaio, pattern=row.padrao, position=position,
                                    start=int(start), stop=int(start+N), width=float(widths(power, lam, True)[0])))
    p_direct, p_dense = np.asarray(p_direct), np.asarray(p_dense)
    checks['max_direct_power_error'] = float(abs(p_direct-data['p']).max())
    checks['max_direct_power_relative_error'] = float((abs(p_direct-data['p'])/data['p']).max())
    assert np.allclose(p_direct, data['p'], rtol=1e-10, atol=1e-10)
    save_json('fontes_sensibilidades_sha256.json', src_hashes)
    np.savez_compressed(OUT/'wps_soma_direta.npz', p_200=p_direct, lambda_200=lam,
                        p_800=p_dense, lambda_800=dense_lam, runs=runs, labels=y)
    sensitivity = [effect(original_width, y, 'ancora_discreta_200'),
                   effect(widths(p_direct, lam, True), y, 'ancora_continua_200'),
                   effect(widths(p_dense, dense_lam, True), y, 'ancora_continua_800')]
    per_run = pd.DataFrame(dict(run=runs, pattern=y, width_original=original_width,
                              width_continuous_200=widths(p_direct, lam, True),
                              width_continuous_800=widths(p_dense, dense_lam, True)))
    for upper in [512, 768]:
        use = lam <= upper
        width = widths(p_direct[:, use], lam[use], True)
        sensitivity.append(effect(width, y, 'ancora_continua_teto_'+str(upper),
                                   dict(lambda_max_effective=float(lam[use][-1]))))
        per_run['width_upper_'+str(upper)] = width
    windows = pd.DataFrame(window_rows)
    windows.to_csv(OUT/'larguras_janelas_distribuidas.csv', index=False)
    for position in range(5):
        v = windows[windows.position == position].set_index('run').loc[runs, 'width'].to_numpy()
        sensitivity.append(effect(v, y, 'janela_uniforme_'+str(position)))
    med = windows.groupby('run').width.median().loc[runs].to_numpy()
    sensitivity.append(effect(med, y, 'mediana_5_janelas_uniformes'))
    spatial = windows.groupby('run').width.agg(['min', 'max', 'median', 'std']).loc[runs]
    per_run['spatial_min'] = spatial['min'].to_numpy()
    per_run['spatial_max'] = spatial['max'].to_numpy()
    per_run['spatial_median'] = spatial['median'].to_numpy()
    per_run['spatial_std'] = spatial['std'].to_numpy()
    per_run.to_csv(OUT/'larguras_por_ensaio.csv', index=False)
    pd.DataFrame(sensitivity).to_csv(OUT/'sensibilidades_descritivas.csv', index=False)

    # Operational support: independent evidence of the scope of exchangeability.
    support = inv[['ensaio','padrao']].copy()
    support['m'] = support.ensaio.str.extract(r'm(\d+)')[0].astype(int)
    support['q'] = support.ensaio.str.extract(r'q(\d+)')[0].astype(int)
    support.to_csv(OUT/'suporte_operacional.csv', index=False)
    checks['unique_joint_mq_blocks'] = int(support.groupby(['m','q']).ngroups)
    checks['n_runs'] = len(runs)
    checks['quantile_log_grid_step'] = float(np.diff(np.log(lam)).mean())
    checks['overlap_fraction_anchor_shift64'] = (N-64)/N
    checks['python'] = platform.python_version()
    checks['numpy'] = np.__version__
    checks['pandas'] = pd.__version__
    checks['scipy'] = scipy.__version__
    checks['script_sha256'] = hashlib.file_digest(Path(__file__).open('rb'),'sha256').hexdigest()
    save_json('verificacao_independente.json', checks)
    print(json.dumps(checks, ensure_ascii=False, indent=2))
    print(pd.DataFrame(sensitivity).to_string(index=False))


if __name__ == '__main__':
    main()
