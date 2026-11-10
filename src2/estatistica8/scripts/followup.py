"""Controles previstos e aprofundamento exploratório após a triagem inicial.

Não muda a família inferencial principal nem redefine seus valores-p.
"""
from core8 import *
from analyze import permutations, permuted_statistics, tail, load
from itertools import product
from scipy.signal import lfilter
from scipy.stats import spearmanr
import time


def conditional_labels(y, blocks):
    _, codes = np.unique(y, return_inverse=True)
    slots, options = [], []
    for block in np.unique(blocks):
        ix = np.flatnonzero(blocks == block)
        slots.append(ix)
        assignments = []
        for chosen in combinations(range(len(ix)), sum(codes[ix] == 0)):
            a = np.ones(len(ix), int); a[list(chosen)] = 0
            assignments.append(a)
        options.append(assignments)
    result = []
    for choice in product(*options):
        row = codes.copy()
        for ix, a in zip(slots, choice): row[ix] = a
        result.append(row)
    return np.asarray(result), codes


def main():
    start = time.monotonic()
    data = load('ancora'); dist, fitted, _, _ = matrices(data)
    inv = pd.read_csv(ROOT/'dados/inventario.csv')
    m = inv.ensaio.str.extract(r'm(\d+)')[0].astype(int).to_numpy()
    q = inv.ensaio.str.extract(r'q(\d+)')[0].astype(int).to_numpy()
    mask = np.isin(data['labels'], SCOPES['churn_annular'])
    labels = data['labels'][mask]
    conditional, covariates = [], []
    for block_name, blocks in [('m', m[mask]), ('q', q[mask])]:
        yp, codes = conditional_labels(labels, blocks)
        base = dist['wps_pearson'][np.ix_(mask, mask)]
        _, base_null = permuted_statistics(base, neighbors(base), yp)
        observed_base = balanced_credit(neighbors(base), labels)[0]
        for name, full in dist.items():
            d = full[np.ix_(mask, mask)]
            x = fitted[name][mask] if name in fitted else None
            w = neighbors(d, features=x)
            st, _, _ = metrics(d, labels, w)
            sil, nn = permuted_statistics(d, w, yp)
            wc = neighbors(d, features=x, base=base)
            _, cn = permuted_statistics(d, wc, yp)
            obs_combo = balanced_credit(wc, labels)[0]
            conditional.append(dict(block=block_name, scope='churn_annular', metric=name,
                permutations=len(yp), exact=True, p_silhouette=tail(sil, st['macro_silhouette'], True),
                p_neighbor=tail(nn, st['balanced_neighbor'], True),
                p_gain_combination=tail(cn-base_null, obs_combo-observed_base, True)))
        print(f'Controle condicional {block_name}: {len(yp)} alocações exatas', flush=True)
    cf = pd.DataFrame(conditional)
    cf[['p_silhouette_holm', 'p_neighbor_holm', 'p_gain_combination_holm']] = holm(
        cf[['p_silhouette', 'p_neighbor', 'p_gain_combination']].to_numpy().ravel()).reshape(-1, 3)
    save_csv('permutacoes_condicionais_exploratorias.csv', cf)
    # Evaluate every candidate combination in both shifted windows and operational
    # holdouts, so the followup does not silently discard unfavorable candidates.
    combo_rows, prediction_rows, covariates = [], [], []
    for condition in ['ancora', 'desloc_m64', 'desloc_p64']:
        dd = load(condition); distances, fit, _, frows = matrices(dd)
        for scope, patterns in SCOPES.items():
            mk = np.isin(dd['labels'], patterns)
            yy, rr = dd['labels'][mk], dd['runs'][mk]
            base = distances['wps_pearson'][np.ix_(mk, mk)]
            for name, full in distances.items():
                d = full[np.ix_(mk, mk)]
                x = fit[name][mk] if name in fit else None
                for block_name, block in [('leave_run_out', None), ('leave_m_out', m[mk]), ('leave_q_out', q[mk])]:
                    # These controls are evaluated only if training retains all classes.
                    if block is not None and any(set(yy[block != block[i]]) != set(yy) for i in range(len(yy))):
                        continue
                    for combination in [False, True]:
                        w = neighbors(d, features=x, blocked=block, base=base if combination else None)
                        value, credit = balanced_credit(w, yy)
                        combo_rows.append(dict(condition=condition, scope=scope, metric=name,
                            combination=combination, control=block_name, balanced_neighbor=value))
                        if scope == 'churn_annular' and name in ['wps_width', 'wps_pearson']:
                            for r, y, c in zip(rr, yy, credit):
                                prediction_rows.append(dict(condition=condition, metric=name, combination=combination,
                                    control=block_name, run=r, pattern=y, credit=c))
        if condition == 'ancora':
            frame = pd.DataFrame(frows).query('family == "wps"').set_index('run').loc[dd['runs']]
            for scope, patterns in SCOPES.items():
                mk = np.isin(dd['labels'], patterns)
                for variable, cv in [('m', m), ('q', q), ('comprimento_original', inv.panorama_px.to_numpy())]:
                    covariates.append(dict(scope=scope, feature='wps_width', covariate=variable,
                                          spearman=spearmanr(frame.width.to_numpy()[mk], cv[mk]).statistic))
    save_csv('ganho_sensibilidades.csv', combo_rows)
    save_csv('previsoes_largura_pearson.csv', prediction_rows)
    save_csv('largura_covariaveis_por_escopo.csv', covariates)
    # Known-independent pairs on exactly the current grid and length.
    wave = Wavelet(data['lam'])
    noise_rows = []
    for mid, model in enumerate(['branco', 'ar1_095']):
        for rep in range(99):
            rng = np.random.default_rng(np.random.SeedSequence([CFG['seed'], 50, mid, rep]))
            xy = rng.normal(size=(2, N+1000))
            if mid: xy = lfilter([1], [1, -.95], xy, axis=1)
            x, y = [prepare(row[-N:]) for row in xy]
            px, _, _ = wave.auto(x); py, _, _ = wave.auto(y)
            co = wave.cross(x, y, px, py)
            s, _, _, _ = summary(co, px, py, data['lam'])
            noise_rows.append(dict(model=model, replicate=rep, **s))
        print(f'Ruído {model}: 99 pares concluídos', flush=True)
    nf = save_csv('ruidos_independentes.csv', noise_rows)
    summary_rows = []
    for model, f in nf.groupby('model'):
        summary_rows.append(dict(model=model, pairs=len(f), median_max=f.global_max.median(),
            q05_max=f.global_max.quantile(.05), q95_max=f.global_max.quantile(.95),
            fraction_above09=(f.global_max >= .9).mean(),
            fraction_upper_edge=(f.peak_lambda_px == 1024).mean()))
    save_csv('resumo_ruidos.csv', summary_rows)
    dec = pd.read_csv(ROOT/'dados/decisoes.csv')
    wco = dec[dec.metric.str.startswith(('wco_', 'kernel_')) & dec.scope.isin(CFG['escopos_prioritarios'])]
    assert not wco.screen_pass.any(), 'Ativou-se condição para executar nulo de fases.'
    json_write(ROOT/'dados/controles.json', dict(status='PASS', tempo_s=time.monotonic()-start,
        fases_executadas=False,
        motivo='Nenhum resumo WCO passou a triagem prioritária. Conforme protocolo, o nulo de fases não foi executado: não pode transformar ausência de separação em identificação útil.',
        aprofundamento='Permutações condicionais e combinações nas sensibilidades são exploração posterior à triagem; não substituem família principal nem validação externa.'))


if __name__ == '__main__': main()
