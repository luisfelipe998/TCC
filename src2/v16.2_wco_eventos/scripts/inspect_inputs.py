"""Auditoria e atlas de imagens, sem importar ou calcular wavelets."""
from pathlib import Path
import hashlib
import json
import os
import time
os.environ.setdefault('MPLCONFIGDIR', '/tmp/wco-eventos-mpl')
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT.parent / 'v16_spatial2_18_71'
DATA = ROOT.parent / 'data_processed'
LABELS = {'dispersed bubbles': 'Bolhas dispersas', 'plug': 'Pistonado (plug)',
          'slug': 'Golfadas (slug)', 'churn': 'Agitado (churn)',
          'wavy annular': 'Anular ondulado'}
PILOT = ['m07q01_001', 'm05q01_001', 'm04q01_001', 'm03q05_001', 'm03q10_001']

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    t0 = time.perf_counter()
    out = ROOT / 'results'; out.mkdir(exist_ok=True)
    figs = ROOT / 'figures' / 'observabilidade'; figs.mkdir(parents=True, exist_ok=True)
    inventory = pd.read_csv(BASE / 'results/input_audit.csv')
    assert len(inventory) == 31 and inventory.included.sum() == 26
    records, source_hashes = [], {}
    # Freeze reference inputs without executing the base analysis.
    for path in BASE.rglob('*'):
        if path.is_file() and path.suffix in {'.py', '.md', '.csv', '.json'}:
            source_hashes[str(path.relative_to(ROOT.parent))] = digest(path)
    for r in inventory.itertuples(index=False):
        folder = DATA / r.run
        for filename, expected in [('__panorama_v11.png', r.panorama_sha256),
                                   ('__dados_cropped_v11.xlsx', r.xlsx_sha256)]:
            f = folder / filename
            assert digest(f) == expected, f'Fonte modificada: {f}'
            source_hashes[str(f.relative_to(ROOT.parent))] = expected
        pano = np.asarray(Image.open(folder / '__panorama_v11.png').convert('L'))
        cachefile = BASE / 'intermediate' / f'{r.run}_input.npz'
        cache = np.load(cachefile)
        source_hashes[str(cachefile.relative_to(ROOT.parent))] = digest(cachefile)
        raw = pano[18:72].mean(axis=0)
        assert np.array_equal(raw, cache['raw'])
        assert np.array_equal(cache['seams'], np.cumsum(cache['shifts'])[:-1])
        assert cache['shifts'].sum() == pano.shape[1] == r.n
        frames = sorted(folder.glob('ImgA*.jpg'))
        assert len(frames) == r.frames == 501
        record = dict(aquisicao=r.run, padrao=r.pattern, incluida=bool(r.included),
                      n=r.n, linhas='18–71', erro_extracao=0.,
                      deslocamentos_nao_positivos=r.nonpositive_shifts,
                      quadros=len(frames), piloto=r.run in PILOT)
        records.append(record)
        if not r.included:
            continue
        fig, axs = plt.subplots(7, 1, figsize=(15, 10.5), constrained_layout=True)
        fig.suptitle(f'{r.run} | {LABELS[r.pattern]} | Inspeção anterior à WCO', fontsize=14)
        axs[0].imshow(pano, cmap='gray', vmin=0, vmax=255, aspect='auto', extent=[0,r.n,80,0])
        axs[0].set_title('Panorama completo; compressão horizontal somente nesta visualização', fontsize=10)
        for ax, start, name in zip(axs[1:4], [0, (r.n-1024)//2, r.n-1024], ['Início','Centro','Fim']):
            ax.imshow(pano[:, start:start+1024], cmap='gray', vmin=0,vmax=255,aspect='auto',extent=[start,start+1024,80,0])
            ax.set_title(f'{name}: recorte fixo de 1.024 px',fontsize=10)
        for ax, idx in zip(axs[4:], [0,250,500]):
            f = frames[idx]; im=np.asarray(Image.open(f).convert('L'))
            source_hashes[str(f.relative_to(ROOT.parent))] = digest(f)
            ax.imshow(im,cmap='gray',vmin=0,vmax=255,aspect='auto',extent=[0,im.shape[1],80,0])
            ax.axvline(im.shape[1]//2, color='#e49717', lw=1)
            ax.set_title(f'Quadro original {f.name}; traço vertical = centro usado na montagem',fontsize=10)
        for ax in axs:
            ax.axhline(18, color='#168b9a', lw=.8); ax.axhline(72, color='#168b9a', lw=.8)
            ax.set_ylabel('Linha');ax.set_xlabel('Coluna (px)')
        fig.savefig(figs / f'{r.run}.png', dpi=120);plt.close(fig)
    pd.DataFrame(records).to_csv(out/'auditoria_entradas.csv',index=False)
    (out/'fontes_sha256.json').write_text(json.dumps(source_hashes,indent=2,ensure_ascii=False)+'\n')
    (out/'inspecao_metadados.json').write_text(json.dumps(dict(
        status='PASS',auditadas=31,elegiveis=26,piloto=PILOT,
        panoramas_planilhas_verificados=62,quadros_fixos_por_aquisicao=3,
        faixa=[18,71],reconstrucao_reexecutada=False,
        segundos=time.perf_counter()-t0),indent=2,ensure_ascii=False)+'\n')
    print(f'Auditoria PASS; atlas de 26 aquisições; {time.perf_counter()-t0:.1f} s.')

if __name__=='__main__':main()
