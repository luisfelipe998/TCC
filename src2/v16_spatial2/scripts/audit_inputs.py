"""Read-only reconstruction audit; all generated artifacts stay in v16."""
from pathlib import Path
import hashlib
import importlib.util
import json
import numpy as np
import pandas as pd
from PIL import Image
from scipy.signal import correlate, correlation_lags

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT.parent / 'data_processed'
EXCLUDED = {'m03q01_001', 'm03q80_001', 'm04q05_001', 'm05q05_002', 'm07q10_001'}

def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def shifts_fast(a, b):
    """Same overlapping Pearson objective as v11, via FFT and prefix sums."""
    n = len(a)
    lag = np.arange(-min(200, n // 4), min(200, n // 4) + 1)
    lo_a, lo_b = np.maximum(lag, 0), np.maximum(-lag, 0)
    count = n - np.abs(lag)
    def moments(x, lo):
        c = np.r_[0., np.cumsum(x)]
        c2 = np.r_[0., np.cumsum(x*x)]
        return c[lo+count]-c[lo], c2[lo+count]-c2[lo]
    sa, qa = moments(a, lo_a)
    sb, qb = moments(b, lo_b)
    numer = correlate(a, b, mode='full', method='fft')[lag+n-1] - sa*sb/count
    denom = np.sqrt(np.maximum(qa-sa*sa/count, 0)*np.maximum(qb-sb*sb/count, 0))
    score = np.divide(numer, denom, out=np.full(len(lag), -np.inf), where=denom>0)
    if not np.isfinite(score).any():
        return 0, float('nan')
    i = np.argmax(score)
    return int(lag[i]), float(score[i])

def main():
    out = ROOT / 'results'
    cache = ROOT / 'intermediate'
    cache.mkdir(exist_ok=True)
    spec = importlib.util.spec_from_file_location('pipeline_v11', ROOT.parent/'pipeline_v11.py')
    pipeline = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pipeline)
    rows = []
    for path in sorted(DATA.glob('m*/__panorama_v11.png')):
        folder = path.parent
        files = sorted(folder.glob('ImgA*.jpg'))
        pano = np.asarray(Image.open(path).convert('L'))
        profiles = np.array([np.asarray(Image.open(f).convert('L'), dtype=float)[20:60].mean(axis=0) for f in files])
        measured = [shifts_fast(a,b) for a,b in zip(profiles[:-1], profiles[1:])]
        raw_shifts = np.array([s for s,_ in measured])
        scores = np.array([s for _,s in measured])
        # Check every run against the literal source implementation at three pairs.
        check_indices = sorted(set([0, len(measured)//2, len(measured)-1]))
        literal = [pipeline._best_shift(profiles[i], profiles[i+1], min(200, profiles.shape[1]//4)) for i in check_indices]
        if not np.array_equal(raw_shifts[check_indices], literal):
            raise AssertionError(f'Pearson shift implementation mismatch: {folder.name}')
        shifts = np.maximum(np.r_[raw_shifts, raw_shifts[-1]], 1)
        rebuilt, _ = pipeline._build_canvas([f.name for f in files], shifts, str(folder))
        exact = rebuilt.shape == pano.shape and np.array_equal(rebuilt, pano)
        book = folder/'__dados_cropped_v11.xlsx'
        frame = pd.read_excel(book, sheet_name='Sheet')
        raw = frame['Voltage_2'].to_numpy(float)
        expected = pano[20:60].mean(axis=0)
        if len(raw) != pano.shape[1] or not np.isfinite(raw).all():
            raise AssertionError(f'Invalid signal length/values: {folder.name}')
        err = float(np.max(np.abs(raw-expected)))
        if err > 1e-10 or not np.array_equal(frame.position_px, np.arange(len(raw))):
            raise AssertionError(f'Panorama/xlsx mismatch: {folder.name}')
        labels = [p.stem.lstrip('_') for p in folder.glob('*.txt') if p.stem.lstrip('_') != 'dados']
        if len(labels) != 1:
            raise AssertionError(f'Ambiguous label: {folder.name}')
        meta = dict(pd.read_excel(book, sheet_name='metadata').values)
        seams = np.cumsum(shifts)[:-1]
        np.savez_compressed(cache/f'{folder.name}_input.npz', raw=raw, center=pano[40].astype(float),
                            upper=pano[20:40].mean(axis=0), lower=pano[40:60].mean(axis=0),
                            raw_shifts=raw_shifts, shifts=shifts, scores=scores, seams=seams)
        rows.append(dict(run=folder.name, pattern=labels[0], included=folder.name not in EXCLUDED,
                         n=len(raw), height=pano.shape[0], frames=len(files), panorama_sha256=sha256(path),
                         xlsx_sha256=sha256(book), xlsx_max_error=err, reconstruction_exact=exact,
                         rebuilt_width=rebuilt.shape[1], nonpositive_shifts=int((raw_shifts<=0).sum()),
                         shift_median=float(np.median(shifts)), shift_min=int(shifts.min()), shift_max=int(shifts.max()),
                         score_min=float(np.nanmin(scores)), score_median=float(np.nanmedian(scores)),
                         metadata_shift_mean=float(meta['shift_mean_px']), recomputed_shift_mean=float(shifts.mean()),
                         raw_std=float(raw.std()), saturation_fraction=float((pano[20:60]==255).mean())))
        print(folder.name, 'exact reconstruction:', exact, 'nonpositive shifts:', rows[-1]['nonpositive_shifts'], flush=True)
    pd.DataFrame(rows).to_csv(out/'input_audit.csv', index=False)
    (out/'input_audit_summary.json').write_text(json.dumps(dict(runs=len(rows), included=sum(r['included'] for r in rows),
        all_xlsx_match=all(r['xlsx_max_error']<1e-10 for r in rows), all_reconstruction_exact=all(r['reconstruction_exact'] for r in rows),
        literal_shift_checks=3*len(rows)), indent=2)+'\n')

if __name__ == '__main__':
    main()
