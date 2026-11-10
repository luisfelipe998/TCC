"""Audita entradas, congelamento, tabelas, mapas e entrega sem refazer o lote."""
from pathlib import Path
from datetime import datetime, timezone
import contextlib
import hashlib
import importlib
import importlib.metadata
import io
import itertools
import json
import re
import shutil
import subprocess
import sys
import tempfile
import time
import numpy as np
import pandas as pd
from PIL import Image
from analysis_common import ROOT, BASE, DATA, CFG, LAM, crop_prepare, evaluate
from run_analysis import VARIANTS


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def close(a, b, tol=1e-11):
    np.testing.assert_allclose(a, b, atol=tol, rtol=tol, equal_nan=True)


def main():
    start = time.perf_counter()
    out = ROOT / 'results'
    checks = {}
    extension = out / 'extensao_churn_anular/execucao.json'
    extension_ready = extension.exists() and json.loads(extension.read_text()).get('status') == 'PASS'
    if extension_ready:
        from verify_texture_extension import main as verify_extension
        verify_extension()
        checks['extensao_churn_anular'] = json.loads((out / 'extensao_churn_anular/verificacao_final.json').read_text())
    read = lambda name: pd.read_csv(out / name)
    freeze = json.loads((out / 'congelamento.json').read_text())
    for filename, expected in freeze['sha256'].items():
        assert digest(ROOT / filename) == expected, filename
    checks['arquivos_congelados_inalterados'] = len(freeze['sha256'])
    assert digest(ROOT / 'scripts/wavelet_core.py') == digest(BASE / 'scripts/wavelet_core.py')
    checks['nucleo_identico_a_base'] = True
    original = (ROOT / 'PLAN_ORIGINAL.md').read_bytes()
    assert (ROOT / 'PLAN.md').read_bytes().endswith(original)
    checks['plano_original_sha256'] = digest(ROOT / 'PLAN_ORIGINAL.md')

    sources = json.loads((out / 'fontes_sha256.json').read_text())
    anchors = json.loads((out / 'quadros_catalogo_sha256.json').read_text())
    for filename, expected in {**sources, **anchors}.items():
        assert digest(ROOT.parent / filename) == expected, filename
    checks['fontes_e_quadros_inalterados'] = len(sources.keys() | anchors.keys())
    base_manifest = json.loads((BASE / 'results/artifact_manifest.json').read_text())
    for filename, expected in base_manifest['output_sha256'].items():
        assert digest(BASE / filename) == expected, f'Base: {filename}'
    for filename, expected in base_manifest['source_sha256'].items():
        assert digest(filename) == expected, f'Fonte anterior: {filename}'
    checks['artefatos_anteriores_preservados'] = len(base_manifest['output_sha256'])
    checks['fontes_anteriores_preservadas'] = len(base_manifest['source_sha256'])

    inv = pd.read_csv(BASE / 'results/input_audit.csv')
    audit = read('auditoria_entradas.csv')
    assert len(inv) == len(audit) == 31 and inv.included.sum() == 26
    excluded = {'m03q01_001', 'm03q80_001', 'm04q05_001', 'm05q05_002', 'm07q10_001'}
    assert set(inv.loc[~inv.included, 'run']) == excluded
    close(audit.set_index('aquisicao').loc[inv.run, 'incluida'], inv.included)
    images, signals, caches = {}, {}, {}
    for r in inv.itertuples():
        im = np.asarray(Image.open(DATA / r.run / '__panorama_v11.png').convert('L'))
        with np.load(BASE / 'intermediate' / f'{r.run}_input.npz') as z:
            cache = {k: z[k] for k in ['raw', 'shifts', 'seams']}
        assert np.array_equal(im[18:72].mean(axis=0), cache['raw'])
        assert np.array_equal(np.cumsum(cache['shifts'])[:-1], cache['seams'])
        assert int(cache['shifts'].sum()) == im.shape[1] == r.n
        images[r.run], caches[r.run] = im, cache
        signals[r.run] = {'principal': im[18:72].mean(axis=0), 'faixa30_49': im[30:50].mean(axis=0)}
    checks['extracoes_exatas_18_71_e_costuras'] = len(images)

    cat = read('catalogo_eventos.csv').set_index('event_id')
    accepted = cat[cat.incluida]
    assert len(cat) == 163 and len(accepted) == 26 and accepted.aquisicao.nunique() == 6
    assert accepted.qualidade.eq('revisada').all()
    assert not set(cat.aquisicao) & excluded
    for eid, r in accepted.iterrows():
        assert r.extensao == r.fim - r.inicio and 96 <= r.extensao <= 1536
        assert r.ancora == (r.inicio + r.fim) // 2
        assert r.segmento_fim - r.segmento_inicio == 2048
        assert r.segmento_inicio >= 16 and r.segmento_fim <= images[r.aquisicao].shape[1] - 16
        assert min(r.contexto_esquerdo, r.contexto_direito) >= 256
        assert r.contexto_esquerdo == r.inicio - r.segmento_inicio
        assert r.contexto_direito == r.segmento_fim - r.fim
        g = accepted[accepted.aquisicao == r.aquisicao]
        assert r.outros_eventos_no_contexto == ((g.inicio < r.segmento_fim) & (g.fim > r.segmento_inicio)).sum() - 1
    # Refaz o vínculo da âncora com a faixa da montagem e verifica seus pixels.
    pixel_checks = 0
    for eid, r in cat.iterrows():
        shifts = caches[r.aquisicao]['shifts']
        bounds = np.r_[0, np.cumsum(shifts)]
        i = int(np.searchsorted(bounds, r.ancora, side='right') - 1)
        assert r.quadro_ancora == f'ImgA{i:06d}.jpg'
        frame = np.asarray(Image.open(DATA / r.aquisicao / r.quadro_ancora).convert('L'))
        col = frame.shape[1] // 2 - int(shifts[i]) // 2 + int(r.ancora - bounds[i])
        assert r.coluna_quadro_ancora == col and 0 <= col < frame.shape[1]
        if not r.ancora_em_zona_mistura:
            assert np.array_equal(images[r.aquisicao][:, r.ancora], frame[:, col]), eid
            pixel_checks += 1
    checks['ancoras_rastreaveis'] = len(cat)
    checks['ancoras_fora_de_mistura_pixels_identicos'] = pixel_checks

    candidates = read('pares_candidatos.csv')
    selected = read('pares_eventos.csv')
    assert len(candidates) == 264 and len(selected) == 165
    assert not candidates.pair_id.duplicated().any()
    chosen = set()
    for r in candidates.itertuples():
        x, y = accepted.loc[r.event_x], accepted.loc[r.event_y]
        assert x.aquisicao == r.run_x and y.aquisicao == r.run_y and r.run_x < r.run_y
        ratio = max(x.extensao, y.extensao) / min(x.extensao, y.extensao)
        close(r.razao_extensoes, ratio)
        relation = 'equivalentes' if ratio <= 1.25 else 'distintos' if ratio >= 1.5 else 'intermediarios'
        scope = 'intra' if x.padrao == y.padrao else 'entre'
        assert r.grupo == f'{scope}_{relation}'
        assert r.prioridade == hashlib.sha256(f'{CFG["semente"]}|{r.pair_id}|{r.grupo}'.encode()).hexdigest()
    for (_, _, group), g in candidates.groupby(['run_x', 'run_y', 'grupo']):
        if not group.endswith('intermediarios'):
            chosen.update(g.sort_values('prioridade').head(12).pair_id)
    assert chosen == set(selected.pair_id) == set(candidates.loc[candidates.selecionado, 'pair_id'])
    for eid, r in accepted.iterrows():
        eq = candidates[candidates.grupo.str.endswith('equivalentes')]
        assert r.parceiros_equivalentes_outra_aquisicao == ((eq.event_x == eid) | (eq.event_y == eid)).sum()
    checks['pareamento_sem_wavelets_recalculado'] = {'candidatos': len(candidates), 'selecionados': len(chosen)}

    # Reprodução do detector/catalogador em diretório isolado; revisão é insumo.
    import catalog_events
    original_root = catalog_events.ROOT
    with tempfile.TemporaryDirectory(prefix='.verificacao_eventos_', dir=ROOT.parent) as tmp:
        temp = Path(tmp)
        (temp / 'results').mkdir()
        shutil.copyfile(out / 'revisao_visual.csv', temp / 'results/revisao_visual.csv')
        catalog_events.ROOT = temp
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                catalog_events.main()
            reproducible = ['catalogo_eventos.csv', 'pares_candidatos.csv', 'pares_eventos.csv', 'segmentos_fixos.csv', 'observabilidade.csv']
            for filename in reproducible:
                assert digest(temp / 'results' / filename) == digest(out / filename), filename
        finally:
            catalog_events.ROOT = original_root
    checks['catalogo_reproduzido_byte_a_byte'] = reproducible

    events, fixed = read('resultados_eventos.csv'), read('resultados_fixos.csv')
    assert len(events) == 1320 and len(fixed) == 1110
    assert set(events.variante) == set(VARIANTS)
    assert not events.duplicated(['pair_id', 'variante']).any()
    for _, g in events.groupby('variante'):
        assert set(g.pair_id) == chosen
    assert fixed.groupby('variante').size().to_dict() == {'principal': 975, 'v5_mesmos_segmentos': 45, 'faixa30_49': 45, 'somente_media': 45}
    primary_fixed = fixed[fixed.variante == 'principal']
    acq_pairs = primary_fixed[['run_x', 'run_y', 'grupo']].drop_duplicates()
    assert acq_pairs.grupo.value_counts().to_dict() == {'entre_fixo': 238, 'intra_fixo': 87}
    assert set(zip(acq_pairs.run_x, acq_pairs.run_y)) == set(itertools.combinations(sorted(inv.loc[inv.included, 'run']), 2))
    assert primary_fixed.groupby(['run_x', 'run_y']).size().eq(3).all()
    for frame in [events, fixed]:
        assert frame.escalas_validas.eq(200).all() and frame.escalas_invalidas.eq(0).all()
        assert frame.lags_principais.eq(65).all() and frame.suporte_geometrico_fracao.eq(1).all()
        assert np.isfinite(frame.escore).all() and frame.erro_equacao.max() < 1e-12
        assert frame.minimo_integral.min() >= 0 and frame.maximo_integral.max() <= 1 + 1e-8
    diag = read('controles_mesma_aquisicao.csv')
    assert diag.tipo.value_counts().to_dict() == {'self': 26, 'mesma_aquisicao': 5}
    assert diag.loc[diag.tipo == 'self', 'erro_self_zero'].max() < 1e-12
    checks['comparacoes_e_escopo'] = {'eventos': len(events), 'fixos': len(fixed), 'self_e_mesma_aquisicao': len(diag), 'pares_intra': 87, 'pares_entre': 238}

    # Recalcula uma comparação em cada variante a partir do PNG, inclusive v5.
    sample_id = sorted(chosen)[0]
    sample = selected.set_index('pair_id').loc[sample_id]
    max_error = 0.
    for variant, (band, linear, engine, dx, dy) in VARIANTS.items():
        x = crop_prepare(signals[sample.run_x][band], int(accepted.loc[sample.event_x].segmento_inicio) + dx, linear=linear)
        y = crop_prepare(signals[sample.run_y][band], int(accepted.loc[sample.event_y].segmento_inicio) + dy, linear=linear)
        _, stats = evaluate(x, y, engine=engine)
        row = events[(events.pair_id == sample_id) & (events.variante == variant)].iloc[0]
        for key, val in stats.items():
            close(val, row[key])
        max_error = max(max_error, abs(stats['escore'] - row.escore))
    checks['recalculo_oito_variantes_erro_maximo'] = float(max_error)

    # Usa os coeficientes integrais, sem chamar o resumidor para conferir mapas.
    map_errors = []
    for path in sorted((ROOT / 'maps').glob('*.npz')):
        with np.load(path) as z:
            assert z['wco'].shape == (200, 4095)
            assert np.array_equal(z['lags'], np.arange(-2047, 2048))
            close(z['wavelengths'], LAM)
            close(z['scales'], LAM * (6 + np.sqrt(38)) / (4 * np.pi))
            denx, deny = abs(z['wxx'][:, 2047]), abs(z['wyy'][:, 2047])
            valid = (denx > denx.max() * 1e-10) & (deny > deny.max() * 1e-10)
            assert np.array_equal(valid, z['denominator_valid'])
            wco = abs(z['wxy'][valid]) ** 2 / (denx[valid] * deny[valid])[:, None]
            close(wco, z['wco'][valid])
            assert np.isnan(z['wco'][~valid]).all()
            row = events[(events.pair_id == path.stem) & (events.variante == 'principal')].iloc[0]
            score = wco[:, np.abs(z['lags']) <= 32].mean(axis=1).mean()
            close(score, row.escore)
            close(wco.max(), row.maximo_integral)
            for key, eid in [('x', row.event_x), ('y', row.event_y)]:
                r = accepted.loc[eid]
                raw = signals[r.aquisicao]['principal'][r.segmento_inicio:r.segmento_fim]
                independent = 255 - raw
                t = np.arange(len(raw))
                independent -= np.polyval(np.polyfit(t, independent, 1), t)
                independent -= independent.mean()
                independent /= independent.std()
                close(z[key], independent)
            map_errors.append(abs(score - row.escore))
    checks['mapas_integrais_recalculados'] = len(map_errors)
    checks['erro_maximo_escore_mapas'] = float(max(map_errors))
    with np.load(out / 'perfis_eventos.npz') as z:
        assert z['mean_lag'].shape == (165, 200)
        close(z['mean_lag'].mean(axis=1), events[events.variante == 'principal'].set_index('pair_id').loc[z['pair_ids'], 'escore'])
    checks['perfis_dos_165_pares_conferidos'] = True

    # Confere pesos por aquisição e os controles de cada contraste exportado.
    for raw, name in [(events, 'agregados_eventos.csv'), (fixed, 'agregados_fixos.csv')]:
        keyed = raw.groupby(['variante', 'grupo', 'run_x', 'run_y'])
        for r in read(name).itertuples():
            g = keyed.get_group((r.variante, r.grupo, r.run_x, r.run_y))
            close(r.escore, g.escore.mean())
            assert r.n_segmentos == len(g)
    agg = read('agregados_eventos.csv').set_index(['variante', 'grupo', 'run_x', 'run_y'])
    agf = read('agregados_fixos.csv').set_index(['variante', 'grupo', 'run_x', 'run_y'])
    contrasts = read('contrastes_pareados.csv')
    for r in contrasts.itertuples():
        eq = agg.loc[(r.variante, r.escopo + '_equivalentes', r.run_x, r.run_y)]
        if r.comparador == 'distintos':
            ref = agg.loc[(r.variante, r.escopo + '_distintos', r.run_x, r.run_y)]
        else:
            var = 'principal' if r.variante.startswith(('bordas_', 'ancoras_')) else r.variante
            ref = agf.loc[(var, r.escopo + '_fixo', r.run_x, r.run_y)]
        close([r.equivalentes, r.controle, r.diferenca, r.diferenca_conservador],
              [eq.escore, ref.escore, eq.escore - ref.escore, eq.conservador - ref.conservador])
    by_pattern = read('contrastes_por_padrao.csv')
    for r in by_pattern.itertuples():
        g = contrasts[(contrasts.variante == r.variante) & (contrasts.escopo == r.escopo) & (contrasts.comparador == r.comparador) & (contrasts.padrao == r.padrao)]
        close(r.diferenca, g.diferenca.mean())
        assert r.pares_aquisicoes == len(g)
    for r in read('contrastes_balanceados.csv').itertuples():
        g = by_pattern[(by_pattern.variante == r.variante) & (by_pattern.escopo == r.escopo) & (by_pattern.comparador == r.comparador)]
        close(r.diferenca, g.diferenca.mean())
    for r in read('retirada_aquisicao.csv').itertuples():
        g = contrasts[(contrasts.variante == 'principal') & (contrasts.escopo == r.escopo) & (contrasts.comparador == r.comparador) & (contrasts.padrao == r.padrao) & (contrasts.run_x != r.retirada) & (contrasts.run_y != r.retirada)]
        assert r.pares_restantes == len(g)
        close(r.diferenca, g.diferenca.mean())
    checks['agregacoes_pesos_contrastes_e_retiradas_conferidos'] = True

    links = 0
    relocations_path = out / 'caminhos_historicos.json'
    relocations = json.loads(relocations_path.read_text()) if relocations_path.exists() else {}
    relocated_links = []
    for path in ROOT.glob('*.md'):
        for target in re.findall(r'\[[^\]]*\]\(([^)]+)\)', path.read_text()):
            if '://' in target or target.startswith('#'):
                continue
            local = target.split('#')[0]
            if local in ['results/verificacao_final.json', 'results/manifesto.json']:
                continue  # Gerados ao final desta mesma verificação.
            if not (path.parent / local).exists() and path.name in {'PLAN.md', 'PLAN_ORIGINAL.md'} and local in relocations:
                replacement = relocations[local]
                assert (ROOT / replacement).is_file(), f'Destino histórico ausente: {replacement}'
                relocated_links.append({'documento': path.name, 'caminho_original': local, 'destino_atual': replacement})
                continue  # Texto histórico preservado; destino atual explicitado no cabeçalho.
            assert (path.parent / local).exists(), f'Link quebrado: {path.name}: {target}'
            links += 1
    checks['links_locais_existentes'] = links
    checks['links_historicos_com_destino_atual_documentado'] = relocated_links
    pngs = list((ROOT / 'figures').rglob('*.png'))
    for path in pngs:
        with Image.open(path) as im:
            assert min(im.size) >= 400
            im.verify()
    checks['imagens_png_integras'] = len(pngs)
    pdf = ROOT / ('RESUMO_EXECUTIVO_ETAPA1.pdf' if extension_ready else 'RESUMO_EXECUTIVO.pdf')
    blob = pdf.read_bytes()
    assert blob.startswith(b'%PDF') and blob.rstrip().endswith(b'%%EOF')
    checks['paginas_pdf'] = len(re.findall(rb'/Type /Page\b', blob))
    assert checks['paginas_pdf'] == 7
    checks['arquivo_pdf_etapa1_conferido'] = pdf.name
    if shutil.which('pdftotext'):
        text = subprocess.run(['pdftotext', str(pdf), '-'], check=True, capture_output=True).stdout.decode()
        assert 'Não conseguimos melhorar de forma convincente' in text
        assert '-0,0427' in text and 'Limitações decisivas' in text
        assert 'Qual padrão é qual' in text and 'plug × plug' in text and 'plug × slug' in text
        checks['texto_e_conclusao_pdf_conferidos'] = True

    checks['status'] = 'PASS'
    checks['data_utc'] = datetime.now(timezone.utc).isoformat()
    checks['segundos'] = time.perf_counter() - start
    checks['limite_da_verificacao'] = 'Testes de integridade e consistência computacional; não validam identidade física, anotação independente ou generalização experimental.'
    (out / 'verificacao_final.json').write_text(json.dumps(checks, indent=2, ensure_ascii=False) + '\n')
    versions, modules, module_hashes = {}, {}, {}
    for dist, module in [('numpy', 'numpy'), ('scipy', 'scipy'), ('pandas', 'pandas'), ('PyWavelets', 'pywt'), ('matplotlib', 'matplotlib'), ('Pillow', 'PIL'), ('openpyxl', 'openpyxl')]:
        versions[dist] = importlib.metadata.version(dist)
        mod = importlib.import_module(module)
        modules[module] = mod.__version__
        module_hashes[module] = {'arquivo': str(Path(mod.__file__)), 'sha256': digest(mod.__file__)}
    cwt = importlib.import_module('pywt._cwt')
    module_hashes['pywt._cwt'] = {'arquivo': cwt.__file__, 'sha256': digest(cwt.__file__)}
    output_hashes = {str(p.relative_to(ROOT)): digest(p) for p in sorted(ROOT.rglob('*'))
                     if p.is_file() and '__pycache__' not in p.parts and p != out / 'manifesto.json' and not p.name.startswith('.')}
    manifest = dict(data_utc=datetime.now(timezone.utc).isoformat(), pasta=ROOT.name,
                    base='../v16_spatial2_18_71', python=sys.version, executavel=sys.executable,
                    distribuicoes=versions, modulos=modules, codigo_modulos=module_hashes,
                    discrepancia_pywavelets=versions['PyWavelets'] != modules['pywt'],
                    fontes='fontes_sha256.json e quadros_catalogo_sha256.json; caminhos relativos a src2',
                    parametros=CFG, saidas_sha256=output_hashes)
    (out / 'manifesto.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + '\n')
    for filename, expected in output_hashes.items():
        assert digest(ROOT / filename) == expected, filename
    print(json.dumps(checks, indent=2, ensure_ascii=False))
    print(f'Manifesto: {len(output_hashes)} arquivos conferidos após gravação.')


if __name__ == '__main__':
    main()
