"""Extrai os mapas congelados e recalcula somente as duas sensibilidades locais."""
from core8 import *
import time
from PIL import Image


def main():
    start = time.monotonic()
    (ROOT/'dados').mkdir(exist_ok=True)
    inputs = {}
    def source(path, expected=None):
        value = sha(path)
        if expected is not None: assert value == expected, path
        inputs[str(path.relative_to(SRC.parent))] = value
        return path
    index = pd.read_csv(source(SOURCE/'04_estatistica/dados/indice_ensaios.csv'))
    with np.load(source(SOURCE/'intermediario/estatistica/descritores.npz')) as z:
        frozen = {k: z[k].copy() for k in z.files}
    assert np.array_equal(index.ensaio, frozen['runs'])
    runs, labels, lam = frozen['runs'], index.padrao.to_numpy(str), frozen['wavelengths']
    assert len(runs) == 25 and len(lam) == 200 and lam[0] == 50 and lam[-1] == 1024
    wave = Wavelet(lam)
    raw_full, autos, signals, power, means, areas, widths, checks = [], [], [], [], [], [], [], []
    for i, row in enumerate(index.itertuples()):
        folder = SOURCE/'intermediario/precalculo_wavelets/individuais'/row.cache_individual
        meta = json.loads(source(folder/'metadata.json').read_text())
        entry = {k: np.load(source(folder/(k+'.npy'), h), allow_pickle=False)
                 for k, h in meta['arrays_sha256'].items()}
        image = source(SRC/'data_processed'/row.ensaio/'__panorama_v11.png', row.panorama_sha256)
        with Image.open(image) as f:
            full = np.asarray(f.convert('L'))[18:72].mean(axis=0)
        assert len(full) == row.panorama_px
        np.testing.assert_array_equal(full[row.inicio_px:row.fim_exclusivo_px], entry['raw'])
        x = prepare(entry['raw'])
        np.testing.assert_allclose(x, entry['x'], atol=1e-12)
        p, co, mp = wave.auto(x)
        error_p = abs(p-frozen['p'][i]).max()
        error_auto = abs(co-abs(entry['wxx'])**2/frozen['p'][i, :, None]**2).max()
        error_mp = abs(mp-entry['power_mean']).max()
        assert max(error_p, error_auto, error_mp) < 1e-10, (row.ensaio, error_p, error_auto, error_mp)
        area, width = self_features(co, lam)
        raw_full.append(full); signals.append(x); power.append(p); means.append(mp)
        areas.append(area); widths.append(width)
        checks.append(dict(run=row.ensaio, power_error=error_p, auto_error=error_auto, mean_power_error=error_mp))
    signals, power = np.asarray(signals), np.asarray(power)
    pair_records = []
    for part in ['lote_intra', 'lote_inter']:
        pair_records += json.loads(source(SOURCE/'intermediario'/part/'metricas_wavelets.json').read_text())
    assert len(pair_records) == 300
    loc = {r: i for i, r in enumerate(runs)}
    pair_records.sort(key=lambda r: tuple(sorted([loc[r['run_x']], loc[r['run_y']]])))
    pairs, scores, prof_mean, prof_peak, prof_zero, detail = [], [], [], [], [], []
    pair_checks = []
    for count, row in enumerate(pair_records):
        path = source(SOURCE/row['arquivo'], row['mapa_sha256'])
        i, j = loc[row['run_x']], loc[row['run_y']]
        with np.load(path) as z:
            co = z['wco']
            meta = json.loads(str(z['meta_json']))
            assert z['denominator_valid'].all()
        assert co.shape == (200, 4095)
        assert meta['comparacao']['run_x'] == runs[i] and meta['comparacao']['run_y'] == runs[j]
        assert meta['cache_individual_x'] == index.iloc[i].cache_individual
        assert meta['cache_individual_y'] == index.iloc[j].cache_individual
        assert row['comprimento'] == N and row['perfil'] == 'comparavel'
        # Six fixed positions in the complete pair list, independent of outcomes.
        if count in [0, 19, 77, 149, 224, 299]:
            error = abs(wave.cross(signals[i], signals[j], power[i], power[j])-co).max()
            assert error < 1e-10
            pair_checks.append(dict(pair=count, run_x=runs[i], run_y=runs[j], max_error=error))
        score, mean, peak, zero = summary(co, power[i], power[j], lam)
        assert abs(score['central_mean']-row['wco_media']) < 1e-12
        pairs.append([i, j]); scores.append([score[k] for k in CROSS])
        prof_mean.append(mean); prof_peak.append(peak); prof_zero.append(zero)
        detail.append(dict(condition='ancora', run_x=runs[i], run_y=runs[j],
                           pattern_x=labels[i], pattern_y=labels[j], **score))
        if (count+1) % 75 == 0:
            print(f'Originais: {count+1}/300 mapas, {time.monotonic()-start:.1f}s', flush=True)
    pairs = np.asarray(pairs)
    assert len({tuple(sorted(p)) for p in pairs}) == 300
    data = dict(runs=runs, labels=labels, lam=lam, pairs=pairs, scores=np.asarray(scores),
                p=power, mp=np.asarray(means), self_area=np.asarray(areas), self_width=np.asarray(widths),
                mean=np.asarray(prof_mean), peak=np.asarray(prof_peak), zero=np.asarray(prof_zero))
    np.savez_compressed(ROOT/'dados/perfis_ancora.npz', **data)
    np.testing.assert_allclose(spectral(power, lam)[0]['pearson']**2, 2*(1-frozen['r']), atol=1e-12)
    save_csv('inventario.csv', index)
    np.savez_compressed(ROOT/'dados/controle_completo.npz', runs=runs, labels=labels, lam=lam,
                        p=frozen['p_full'], r=frozen['r_full'])
    for shift in [-64, 64]:
        condition = 'desloc_m64' if shift < 0 else 'desloc_p64'
        xs, ps, mps, areas, widths = [], [], [], [], []
        for row, full in zip(index.itertuples(), raw_full):
            lo, hi = row.inicio_px+shift, row.fim_exclusivo_px+shift
            assert 0 <= lo < hi <= len(full)
            x = prepare(full[lo:hi])
            p, co, mp = wave.auto(x)
            area, width = self_features(co, lam)
            xs.append(x); ps.append(p); mps.append(mp); areas.append(area); widths.append(width)
        scores, prof_mean, prof_peak, prof_zero = [], [], [], []
        for count, (i, j) in enumerate(pairs):
            co = wave.cross(xs[i], xs[j], ps[i], ps[j])
            score, mean, peak, zero = summary(co, ps[i], ps[j], lam)
            scores.append([score[k] for k in CROSS]); prof_mean.append(mean); prof_peak.append(peak); prof_zero.append(zero)
            detail.append(dict(condition=condition, run_x=runs[i], run_y=runs[j],
                               pattern_x=labels[i], pattern_y=labels[j], **score))
            if (count+1) % 100 == 0:
                print(f'{condition}: {count+1}/300 pares, {time.monotonic()-start:.1f}s', flush=True)
        np.savez_compressed(ROOT/f'dados/perfis_{condition}.npz', runs=runs, labels=labels, lam=lam,
                            pairs=pairs, scores=scores, p=ps, mp=mps, self_area=areas, self_width=widths,
                            mean=prof_mean, peak=prof_peak, zero=prof_zero)
    save_csv('resumos_wco.csv', detail)
    save_csv('verificacao_individuais.csv', checks)
    save_csv('verificacao_mapas.csv', pair_checks)
    json_write(ROOT/'dados/fontes_sha256.json', inputs)
    json_write(ROOT/'dados/extracao.json', dict(status='PASS', ensaios=25, pares=300,
               condicoes=3, snapshots_sha256=len(inputs), tempo_s=time.monotonic()-start,
               max_auto_error=max(v['auto_error'] for v in checks),
               max_cross_error=max(v['max_error'] for v in pair_checks)))
    print('Extração e sensibilidades: PASS', flush=True)


if __name__ == '__main__': main()
