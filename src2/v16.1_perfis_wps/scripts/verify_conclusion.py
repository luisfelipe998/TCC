"""Auditoria adicional da conclusão; não altera sinais, métodos ou resultados originais."""
from contextlib import redirect_stdout
from datetime import datetime, timezone
from hashlib import sha256
import importlib.util
from itertools import combinations
import io
import json
from pathlib import Path
import sys
import tempfile
import time

import numpy as np
import pandas as pd
from PIL import Image
from scipy.fft import next_fast_len
from scipy.stats import mannwhitneyu, spearmanr

sys.dont_write_bytecode = True
from profile_core import ROOT, BASE, LAM, wc, zero_profile


def digest(path):
    h = sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def independent_prepare(raw):
    """Retira reta por projeção explícita, sem scipy.signal.detrend."""
    y = 255.0 - np.asarray(raw, float)
    t = np.arange(len(y), dtype=float)
    t -= t.mean()
    y -= y.mean()
    y -= t * np.dot(t, y) / np.dot(t, t)
    y -= y.mean()
    return y / np.sqrt(np.mean(y*y))


def independent_zero(x):
    """WPS(0,a) por soma temporal da Morlet, sem a CWT de produção.

    A autocorrelação linear é construída a partir de |FFT(x)|², com
    preenchimento suficiente; o divisor N é a adaptação declarada.
    A parte imaginária cancela porque a autocorrelação real é par.
    """
    n = len(x)
    size = next_fast_len(2*n-1)
    spectrum = np.fft.rfft(x, n=size)
    positive = np.fft.irfft(np.abs(spectrum)**2, n=size)[:n] / n
    correlation = np.r_[positive[:0:-1], positive]
    lag = np.arange(-n+1, n, dtype=float)
    result = []
    for scale in LAM / (4*np.pi/(6+np.sqrt(38))):
        kernel = np.pi**(-.25)/np.sqrt(scale)
        kernel = kernel*np.exp(-.5*(lag/scale)**2)*np.cos(6*lag/scale)
        result.append(abs(np.sum(correlation*kernel)))
    return np.asarray(result)


def density(values, axis):
    return values / np.trapezoid(values, x=np.log(axis))


def distance(px, py, axis):
    return float(np.sqrt(np.trapezoid(
        (np.sqrt(density(px, axis))-np.sqrt(density(py, axis)))**2,
        x=np.log(axis))/2))


def rank_values(frame):
    """Ordenação por U de Mann–Whitney; nenhum valor-p é utilizado."""
    same = frame[frame.mesmo_padrao]
    cross = frame[~frame.mesmo_padrao]
    def rank(s, c):
        return float(mannwhitneyu(s, c).statistic/(len(s)*len(c)))
    pair = rank(same.similaridade, cross.similaridade)
    balanced = np.mean([
        rank(s.similaridade, c.similaridade)
        for _, s in same.groupby('grupo') for _, c in cross.groupby('grupo')
    ])
    return np.asarray([pair, balanced])


def inherited_controls(path, result_name):
    """Executa controles existentes em pasta temporária, preservando os relatórios antigos."""
    spec = importlib.util.spec_from_file_location('conclusion_control', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    with tempfile.TemporaryDirectory(prefix='tcc-conclusion-check-') as tmp:
        module.ROOT = Path(tmp)
        (module.ROOT/'results').mkdir()
        with redirect_stdout(io.StringIO()):
            module.main()
        return json.loads((module.ROOT/'results'/result_name).read_text())


def main():
    start = time.perf_counter()
    out = ROOT/'results'
    checks = {}
    cfg = json.loads((ROOT/'config.json').read_text())
    assert digest(ROOT/'PROTOCOL.md') == cfg['protocol_sha256']
    checks['core_controls'] = inherited_controls(
        BASE/'scripts/validate_wavelets.py', 'numerical_validation.json')
    checks['profile_controls'] = inherited_controls(
        ROOT/'scripts/validate_operations.py', 'operations_validation.json')
    print('Controles numéricos das fórmulas: PASS', flush=True)

    inv = pd.read_csv(BASE/'results/input_audit.csv')
    selected = inv[inv.included].sort_values('run')
    assert len(inv) == 31 and len(selected) == 26
    assert selected.run.tolist() == cfg['expected_included_runs']
    assert sorted(inv.loc[~inv.included, 'run']) == sorted(cfg['excluded_runs'])
    raw = {}
    for row in inv.itertuples():
        path = ROOT.parent/'data_processed'/row.run/'__panorama_v11.png'
        assert digest(path) == row.panorama_sha256
        with Image.open(path) as image:
            arr = np.asarray(image.convert('L'))
        signal = arr[18:72].sum(axis=0, dtype=np.int64)/54
        with np.load(BASE/'intermediate'/f'{row.run}_input.npz') as cache:
            np.testing.assert_array_equal(signal, cache['raw'])
            assert len(signal) == row.n
        label_files = [p.stem.lstrip('_') for p in path.parent.glob('*.txt')
                       if p.stem.lstrip('_') != 'dados']
        assert label_files == [row.pattern]
        if row.included:
            raw[row.run] = signal
    checks['panoramas_e_extracoes_verificados'] = 31
    nmin = min(map(len, raw.values()))
    primary_mask = (nmin-1 >= np.sqrt(2)*LAM/wc.FOURIER_FACTOR)
    conservative_mask = 3*LAM/wc.FOURIER_FACTOR <= nmin/2
    np.testing.assert_array_equal(primary_mask, cfg['primary_mask'])
    np.testing.assert_array_equal(conservative_mask, cfg['conservative_mask'])
    assert nmin == cfg['common_length'] == 3171
    checks['dominio_principal'] = dict(pontos=int(primary_mask.sum()),
        minimo_px=float(LAM[primary_mask].min()), maximo_px=float(LAM[primary_mask].max()))

    stored_profiles = pd.read_csv(out/'profiles_raw.csv')
    primary = stored_profiles[stored_profiles.variante == 'principal'].pivot(
        index='aquisicao', columns='comprimento_onda_px', values='perfil')
    errors = []
    preparation_errors = []
    independent = {}
    for run, signal in raw.items():
        x = independent_prepare(signal)
        preparation_errors.append(float(np.max(abs(x-wc.prepare(signal)))))
        independent[run] = independent_zero(x)
        expected = primary.loc[run].to_numpy()
        errors.append(float(np.max(abs(independent[run]-expected))))
        np.testing.assert_allclose(independent[run], expected, rtol=1e-10, atol=1e-10)
    checks['perfis_principais_por_soma_temporal'] = 26
    checks['erro_maximo_perfil_soma_temporal'] = max(errors)
    checks['erro_maximo_preparo_independente'] = max(preparation_errors)
    print('26 perfis conferidos por soma temporal independente: PASS', flush=True)

    # Auditar normalização sem usar distribution()/quadrature() da produção.
    norm = pd.read_csv(out/'profiles_normalized.csv')
    norm_error = 0.
    for _, group in norm.groupby(['aquisicao', 'variante', 'dominio']):
        group = group.sort_values('comprimento_onda_px')
        axis = group.comprimento_onda_px.to_numpy()
        log_axis = np.log(axis)
        p = group.perfil.to_numpy()
        area = np.trapezoid(p, x=log_axis)
        intervals = np.diff(log_axis)
        weights = np.r_[intervals[0], intervals[:-1]+intervals[1:], intervals[-1]]/2
        mass = p*weights/area
        np.testing.assert_allclose(mass, group.massa, rtol=1e-10, atol=1e-12)
        np.testing.assert_allclose(p/area, group.densidade_log, rtol=1e-10, atol=1e-12)
        assert abs(mass.sum()-1) < 1e-12
        norm_error = max(norm_error, float(np.max(abs(mass-group.massa))))
    checks['erro_maximo_massas_normalizadas'] = norm_error

    pairs = pd.read_csv(out/'all_pairs.csv')
    mainpairs = pairs[pairs.metodo == 'wps_principal']
    distance_errors = []
    for row in mainpairs.itertuples():
        d = distance(independent[row.aquisicao_x][primary_mask],
                     independent[row.aquisicao_y][primary_mask], LAM[primary_mask])
        distance_errors.append(abs(d-row.valor))
        assert abs(d-row.valor) < 1e-10
    checks['distancias_principais_recalculadas_do_sinal'] = len(distance_errors)
    checks['erro_maximo_distancias_do_sinal'] = max(distance_errors)

    # Recalcular os 350 prefixos e depois os 325 escores WCO originais.
    prefix = pd.read_csv(out/'pair_prefix_profiles.csv')
    cache = {}
    prefix_errors = []
    for (run, n), group in prefix.groupby(['aquisicao', 'n']):
        x = wc.prepare(raw[run][:n])
        p = zero_profile(x)
        expected = group.sort_values('comprimento_onda_px').perfil_v16.to_numpy()
        np.testing.assert_allclose(p, expected, rtol=1e-10, atol=1e-11)
        prefix_errors.append(float(np.max(abs(p-expected))))
        cache[run, int(n)] = (x, p)
    assert len(cache) == 350
    checks['prefixos_recalculados'] = len(cache)
    checks['erro_maximo_perfis_prefixo'] = max(prefix_errors)
    original = pd.read_csv(BASE/'results/preprocessing_comparison_pairs.csv')
    original = original[original.variant == 'baseline'].set_index(['run_x', 'run_y'])
    wco_rows = pairs[pairs.metodo == 'wco_v16_18_71']
    prefix_rows = pairs[pairs.metodo == 'wps_prefixo_par'].set_index(['aquisicao_x', 'aquisicao_y'])
    wco_errors = []
    max_coherence = 0.
    for i, row in enumerate(wco_rows.itertuples(), start=1):
        a, b = row.aquisicao_x, row.aquisicao_y
        n = min(len(raw[a]), len(raw[b]))
        x, px = cache[a, n]
        y, py = cache[b, n]
        wxy = wc.morlet_cwt(wc.correlation(x, y), LAM/wc.FOURIER_FACTOR)
        valid = (px > 1e-10*px.max()) & (py > 1e-10*py.max())
        coherence = abs(wxy[valid])**2/(px[valid]*py[valid])[:, None]
        max_coherence = max(max_coherence, float(coherence.max()))
        assert coherence.min() >= 0 and coherence.max() <= 1+1e-8
        peaks = np.full(len(LAM), np.nan)
        peaks[valid] = coherence.max(axis=1)
        active = valid & (px >= .05*px.max()) & (py >= .05*py.max())
        score = float(np.average(peaks[active], weights=np.sqrt(px[active]*py[active])))
        wco_errors.append(abs(score-row.valor))
        assert abs(score-row.valor) < 1e-10
        assert abs(score-original.loc[(a, b), 'weighted_scale_peak_power05']) < 1e-10
        assert abs(distance(px, py, LAM)-prefix_rows.loc[(a, b), 'valor']) < 1e-10
        del wxy, coherence
        if i % 25 == 0 or i == len(wco_rows):
            print(f'WCO recalculada: {i}/325; {time.perf_counter()-start:.1f} s', flush=True)
    checks['escores_WCO_recalculados_do_sinal'] = len(wco_errors)
    checks['erro_maximo_escores_WCO'] = max(wco_errors)
    checks['maior_WCO_no_plano_integral'] = max_coherence

    summary = pd.read_csv(out/'ranking_summary.csv').set_index('metodo')
    deletion = pd.read_csv(out/'acquisition_deletion.csv')
    expected_pairs = set(combinations(sorted(raw), 2))
    rank_error = deletion_error = 0.
    for method, group in pairs.groupby('metodo'):
        assert len(group) == 325 and int(group.mesmo_padrao.sum()) == 87
        assert set(zip(group.aquisicao_x, group.aquisicao_y)) == expected_pairs
        assert group.mesmo_padrao.eq(group.padrao_x == group.padrao_y).all()
        sign = 1 if method.startswith('wco_') else -1
        np.testing.assert_array_equal(group.similaridade, sign*group.valor)
        expected = summary.loc[method, ['A_pares', 'A_balanceada']].to_numpy(float)
        actual = rank_values(group)
        rank_error = max(rank_error, float(np.max(abs(actual-expected))))
        np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-13)
        for row in deletion[deletion.metodo == method].itertuples():
            sub = group[(group.aquisicao_x != row.aquisicao_retirada) &
                        (group.aquisicao_y != row.aquisicao_retirada)]
            assert len(sub) == 300
            actual = rank_values(sub)
            error = float(np.max(abs(actual-[row.A_pares, row.A_balanceada])))
            deletion_error = max(deletion_error, error)
            assert error < 1e-13
    checks['ordenacoes_recalculadas'] = len(summary)
    checks['erro_maximo_ordenacao'] = rank_error
    checks['retiradas_verificadas_nas_duas_ponderacoes'] = len(deletion)
    checks['erro_maximo_retiradas'] = deletion_error

    pattern_table = pd.read_csv(out/'ranking_by_pattern.csv')
    for row in pattern_table.itertuples():
        group = pairs[pairs.metodo == row.metodo]
        same = group[group.mesmo_padrao & (group.padrao_x == row.padrao)]
        cross = group[~group.mesmo_padrao & ((group.padrao_x == row.padrao) |
                                            (group.padrao_y == row.padrao))]
        assert len(same) == row.pares_intra and len(cross) == row.pares_entre
        if len(same):
            actual = rank_values(pd.concat([same, cross]))
            np.testing.assert_allclose(actual, [row.A_pares, row.A_balanceada], rtol=0, atol=1e-13)
        else:
            assert np.isnan(row.A_pares) and np.isnan(row.A_balanceada)
    checks['resultados_por_padrao_verificados'] = len(pattern_table)
    chosen = pd.read_csv(out/'selected_cases.csv')
    for title, row in [
        ('Controle cruzado mais semelhante', mainpairs[~mainpairs.mesmo_padrao].sort_values('valor').iloc[0]),
        ('Par intrapadrão menos semelhante', mainpairs[mainpairs.mesmo_padrao].sort_values('valor').iloc[-1]),
    ]:
        expected = chosen[chosen.tipo == title].iloc[0]
        assert row.aquisicao_x == expected.aquisicao_x and row.aquisicao_y == expected.aquisicao_y
        assert abs(row.valor-expected.H) < 1e-14
    ratio = np.maximum(mainpairs.n_x, mainpairs.n_y)/np.minimum(mainpairs.n_x, mainpairs.n_y)
    checks['correlacao_descritiva_distancia_comprimento'] = float(
        spearmanr(mainpairs.valor, np.log(ratio)).statistic)
    checks['nota_direcao_WCO'] = (
        'A < 0,5 indica ordenação contrária à direção de semelhança definida. '
        'Não significa ausência de informação discriminativa. Inverter a direção '
        'após observar os rótulos muda a regra; não foi aplicado ao protocolo.')
    checks['nota_fonte_primaria'] = (
        'Páginas impressas 411, 413 e 414 de Schuck Jr. e Bodmann conferidas '
        'visualmente em 25/10/2026: CWT com 1/sqrt(a), Morlet omega0=6, '
        'WPS=|CWT(Rxx)|; WCS complexa; WCO=|WCS|²/[WPSx(0,a)WPSy(0,a)]. '
        'Divisor N é adaptação declarada; o artigo imprime N-|m|.')
    checks['limites'] = (
        'Auditoria computacional. Não reconstrói panoramas, não valida rótulos '
        'fisicamente e não constitui teste externo nem independência dos pares. '
        'A faixa geométrica é escolha de suporte, não garantia de ausência de bordas.')
    files = [ROOT/'config.json', ROOT/'PROTOCOL.md', Path(__file__),
             ROOT/'scripts/profile_core.py', ROOT/'scripts/evaluation.py',
             BASE/'scripts/wavelet_core.py', BASE/'scripts/validate_wavelets.py',
             ROOT/'scripts/validate_operations.py']
    files += [out/name for name in ['profiles_raw.csv', 'profiles_normalized.csv',
              'pair_prefix_profiles.csv', 'all_pairs.csv', 'ranking_summary.csv',
              'ranking_by_pattern.csv', 'acquisition_deletion.csv', 'selected_cases.csv']]
    checks['sha256_fontes'] = {str(p.relative_to(ROOT.parent)): digest(p) for p in files}
    checks['status'] = 'PASS'
    checks['completed_utc'] = datetime.now(timezone.utc).isoformat()
    checks['elapsed_seconds'] = time.perf_counter()-start
    (out/'conclusion_audit.json').write_text(json.dumps(checks, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({k: v for k, v in checks.items() if k not in
                     ['core_controls', 'profile_controls', 'sha256_fontes']}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
