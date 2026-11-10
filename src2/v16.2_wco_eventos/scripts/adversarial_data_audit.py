"""Auditoria independente das entradas; não altera a base nem os panoramas."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import time
import numpy as np
import pandas as pd
from PIL import Image
from scipy.signal import correlate

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT.parent / 'v16_spatial2_18_71'
DATA = ROOT.parent / 'data_processed'
OUT = ROOT / 'results/auditoria_adversarial'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def shifts_from_profiles(profiles):
    """Pearson em todos os lags, por momentos locais; nenhum import do v11."""
    n = profiles.shape[1]
    lag = np.arange(-min(200, n // 4), min(200, n // 4) + 1)
    count = n - abs(lag)
    ia, ib = np.maximum(lag, 0), np.maximum(-lag, 0)
    shifts, scores = [], []
    for a, b in zip(profiles[:-1], profiles[1:]):
        ac, bc = np.r_[0., np.cumsum(a)], np.r_[0., np.cumsum(b)]
        aq, bq = np.r_[0., np.cumsum(a*a)], np.r_[0., np.cumsum(b*b)]
        sa, sb = ac[ia+count]-ac[ia], bc[ib+count]-bc[ib]
        va = aq[ia+count]-aq[ia]-sa*sa/count
        vb = bq[ib+count]-bq[ib]-sb*sb/count
        num = correlate(a, b, mode='full', method='fft')[lag+n-1]-sa*sb/count
        den = np.sqrt(np.maximum(va, 0)*np.maximum(vb, 0))
        corr = np.divide(num, den, out=np.full(len(lag), -np.inf), where=den > 0)
        k = int(np.argmax(corr)); shifts.append(int(lag[k])); scores.append(float(corr[k]))
    return np.array(shifts), np.array(scores)


def main():
    t0 = time.perf_counter()
    inv = pd.read_csv(BASE / 'results/input_audit.csv').sort_values('run')
    rows, hashes, images_by_hash = [], {}, {}
    signals = {}
    for r in inv.itertuples():
        folder = DATA / r.run
        labels = sorted(p.stem.lstrip('_') for p in folder.glob('*.txt') if p.stem.lstrip('_') != 'dados')
        assert labels == [r.pattern], (r.run, labels, r.pattern)
        names = sorted(folder.glob('ImgA*.jpg'))
        assert [p.name for p in names] == [f'ImgA{i:06d}.jpg' for i in range(501)]
        cache = np.load(BASE / 'intermediate' / f'{r.run}_input.npz')
        widths = cache['shifts'].astype(int)
        pieces, profiles = [], []
        for i, path in enumerate(names):
            im = np.asarray(Image.open(path).convert('L'))
            assert im.shape[0] == 80
            center = im.shape[1] // 2
            pieces.append(im[:, center-widths[i]//2:center-widths[i]//2+widths[i]])
            profiles.append(im[20:60].sum(axis=0, dtype=np.uint64)/40)
            h = digest(path); hashes[str(path.relative_to(ROOT.parent))] = h
            images_by_hash.setdefault(h, []).append(str(path.relative_to(DATA)))
        measured, scores = shifts_from_profiles(np.asarray(profiles))
        assert np.array_equal(measured, cache['raw_shifts']), r.run
        assert np.array_equal(widths, np.maximum(np.r_[measured, measured[-1]], 1))
        # Reprodução literal do efeito observado: para largura de feather 2,
        # os pesos são [0,1]; somente a última coluna à esquerda é substituída.
        original = np.concatenate(pieces, axis=1)
        rebuilt = original.copy()
        seams = np.cumsum(widths)[:-1]
        modified = []
        for i, boundary in enumerate(seams):
            feather = min(2, min(widths[i], widths[i+1]) // 2)
            if feather == 2:
                rebuilt[:, boundary-1] = pieces[i+1][:, 1]
                modified.append(int(boundary-1))
        pano_path = folder / '__panorama_v11.png'
        pano = np.asarray(Image.open(pano_path).convert('L'))
        assert np.array_equal(pano, rebuilt), r.run
        raw = pano[18:72].sum(axis=0, dtype=np.uint64)/54
        assert np.array_equal(raw, cache['raw'])
        book = folder / '__dados_cropped_v11.xlsx'
        table = pd.read_excel(book, sheet_name='Sheet')
        assert np.array_equal(table.position_px, np.arange(pano.shape[1]))
        workbook_error = float(np.max(abs(table.Voltage_2.to_numpy()-pano[20:60].sum(axis=0, dtype=np.uint64)/40)))
        assert workbook_error < 1e-10
        for path in [pano_path, book, *folder.glob('*.txt')]:
            hashes[str(path.relative_to(ROOT.parent))] = digest(path)
        no_feather = original[18:72].sum(axis=0, dtype=np.uint64)/54
        signals[r.run+'_sem_feather'] = no_feather
        signals[r.run+'_original'] = raw
        signals[r.run+'_colunas_alteradas'] = np.asarray(modified)
        rows.append(dict(aquisicao=r.run, padrao=r.pattern, incluida=r.included, quadros=len(names),
            largura=pano.shape[1], deslocamentos_conferidos=len(measured), reconstruido_exatamente=True,
            erro_planilha=workbook_error, rotulo_confere=True, nao_positivos=int((measured<=0).sum()),
            correlacao_minima=float(scores.min()), correlacao_mediana=float(np.median(scores)),
            correlacoes_abaixo_08=int((scores<.8).sum()), feather_colunas_substituidas=len(modified),
            fracao_colunas_substituidas=len(modified)/pano.shape[1],
            diferenca_rms_perfil=float(np.sqrt(np.mean((raw-no_feather)**2))),
            diferenca_maxima_perfil=float(abs(raw-no_feather).max())))
        print(f'{r.run}: 501 quadros, 500 deslocamentos e panorama conferidos.', flush=True)
    pd.DataFrame(rows).to_csv(OUT/'dados_por_aquisicao.csv', index=False)
    np.savez_compressed(OUT/'perfis_sem_feather.npz', **signals)
    (OUT/'fontes_sha256.json').write_text(json.dumps(hashes, indent=2, ensure_ascii=False)+'\n')
    repeated = {h:p for h,p in images_by_hash.items() if len(p)>1}
    (OUT/'quadros_duplicados.json').write_text(json.dumps(repeated, indent=2, ensure_ascii=False)+'\n')
    report = dict(status='PASS', data_utc=datetime.now(timezone.utc).isoformat(), aquisicoes=len(rows),
        quadros=sum(r['quadros'] for r in rows), deslocamentos=sum(r['deslocamentos_conferidos'] for r in rows),
        panoramas_reproduzidos=len(rows), planilhas_conferidas=len(rows), rotulos_conferidos=len(rows),
        grupos_quadros_duplicados=len(repeated), fontes_com_hash=len(hashes), segundos=time.perf_counter()-t0,
        limite='Concordância com a montagem v11 e seus rótulos, não validação física de deslocamento, interfaces ou classificação.')
    (OUT/'dados.json').write_text(json.dumps(report, indent=2, ensure_ascii=False)+'\n')
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
