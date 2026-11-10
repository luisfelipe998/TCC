"""Preserve v11 strips and shifts; omit its column-reordering seam overwrite."""
from datetime import datetime, timezone
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.signal import correlate

from v14_paths import ROOT, DATA, BASELINE, PANORAMAS, CACHE, artifact_path
from v14_selection import selected_folders


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def measured_shift(a, b):
    """Original overlap Pearson objective, vectorized without altering lag order."""
    n = len(a)
    limit = min(200, n // 4)
    lag = np.arange(-limit, limit + 1)
    la, lb = np.maximum(lag, 0), np.maximum(-lag, 0)
    count = n - np.abs(lag)
    def moments(x, start):
        total = np.r_[0., np.cumsum(x)]
        square = np.r_[0., np.cumsum(x*x)]
        return total[start+count]-total[start], square[start+count]-square[start]
    sa, qa = moments(a, la)
    sb, qb = moments(b, lb)
    covariance = correlate(a, b, mode='full', method='fft')[lag+n-1] - sa*sb/count
    denominator = np.sqrt(np.maximum(qa-sa*sa/count, 0) * np.maximum(qb-sb*sb/count, 0))
    score = np.divide(covariance, denominator, out=np.full(len(lag), -np.inf), where=denominator>0)
    return int(lag[np.argmax(score)]) if np.isfinite(score).any() else 0


def center_strips(frames, shifts):
    if len(frames) != len(shifts) or not len(frames):
        raise ValueError('One positive strip width per frame is required')
    width = frames[0].shape[1]
    result = []
    for frame, shift in zip(frames, shifts):
        if frame.shape != frames[0].shape or not 0 < shift <= width:
            raise ValueError('Invalid frame shape or strip width')
        start = width//2 - int(shift)//2
        result.append(frame[:, start:start+int(shift)])
    return result


def stitch_strips(strips):
    """Every output pixel retains its original within-strip order and intensity."""
    return np.concatenate(strips, axis=1).copy()


def legacy_stitch(strips):
    """Audit-only reproduction of the defect; never used as the corrected input."""
    out = stitch_strips(strips).astype(np.float32)
    boundaries = np.cumsum([part.shape[1] for part in strips])[:-1]
    for i, x in enumerate(boundaries, 1):
        size = min(2, min(strips[i-1].shape[1], strips[i].shape[1])//2)
        if size:
            weight = np.linspace(0., 1., size)[None, :]
            out[:, x-size:x] = strips[i-1][:, -size:]*(1-weight) + strips[i][:, :size]*weight
    return np.clip(out, 0, 255).astype(np.uint8)


def main():
    os.environ.setdefault('MPLCONFIGDIR', str(CACHE/'matplotlib'))
    original = ROOT.parent / 'pipeline_v11.py'
    spec = importlib.util.spec_from_file_location('original_v11', original)
    source = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source)
    source_hashes = {str(original.relative_to(ROOT.parent)): sha(original)}
    with (BASELINE/'results/validation/data_inventory.csv').open() as stream:
        baseline = {r['run']: r for r in csv.DictReader(stream)}
    runs, shifts_table = [], []
    for folder in selected_folders():
        files = sorted(folder.glob('ImgA*.jpg'))
        frames = []
        for path in files:
            source_hashes[str(path.relative_to(ROOT.parent))] = sha(path)
            with Image.open(path) as image:
                frames.append(np.asarray(image.convert('L')))
        profiles = np.asarray([im[20:60].mean(axis=0) for im in frames])
        raw = np.array([measured_shift(a, b) for a, b in zip(profiles[:-1], profiles[1:])])
        for i in sorted({0, len(raw)//2, len(raw)-1}):
            assert raw[i] == source._best_shift(profiles[i], profiles[i+1], min(200, len(profiles[i])//4))
        shifts = np.maximum(np.r_[raw, raw[-1]], 1)
        strips = center_strips(frames, shifts)
        old_path = folder/'__panorama_v11.png'
        assert sha(old_path) == baseline[folder.name]['sha256'], 'Original v14 input changed'
        with Image.open(old_path) as image:
            old = np.asarray(image.convert('L'))
        assert np.array_equal(legacy_stitch(strips), old), f'Cannot reproduce original input: {folder.name}'
        new = stitch_strips(strips)
        assert new.shape == old.shape and new.dtype == np.uint8
        target = PANORAMAS/folder.name/'__panorama_v11.png'
        target.parent.mkdir(parents=True, exist_ok=True)
        Image.fromarray(new).save(target)
        for path in [old_path, folder/'_config.xsv', *folder.glob('__*.txt')]:
            source_hashes[str(path.relative_to(ROOT.parent))] = sha(path)
        for i, shift in enumerate(shifts):
            shifts_table.append(dict(run=folder.name,frame=files[i].name,strip_width_px=int(shift),
                                     raw_shift_px=int(raw[min(i,len(raw)-1)]),
                                     last_shift_reused=bool(i==len(shifts)-1)))
        runs.append(dict(run=folder.name,n_frames=len(frames),width_px=new.shape[1],height_px=new.shape[0],
                         nonpositive_measured_shifts=int((raw<=0).sum()),
                         legacy_reconstruction_exact=True,changed_pixels=int(np.count_nonzero(new!=old)),
                         original_panorama_sha256=sha(old_path),corrected_panorama_sha256=sha(target)))
        print(f'{folder.name}: original reproduced; corrected input saved', flush=True)
    with artifact_path('reconstruction_shifts.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(shifts_table[0])); writer.writeheader(); writer.writerows(shifts_table)
    artifact_path('source_input_hashes.json').write_text(json.dumps(source_hashes, indent=2)+'\n')
    result = dict(passed=True,created_utc=datetime.now(timezone.utc).isoformat(),runs=runs,
                  source_files=len(source_hashes),source_files_unchanged=all(sha(ROOT.parent/p)==h for p,h in source_hashes.items()),
                  correction='Concatenate unchanged center strips without the v11 seam overwrite; no new interpolation or smoothing.',
                  preserved='Frame order, Pearson shift objective, clamp to >=1 px, last-shift reuse, strip widths, lumen and all downstream parameters.',
                  remaining_limit='Frozen-advection reconstruction and measured shifts remain assumptions; no phase or physical void calibration.')
    assert result['source_files_unchanged']
    artifact_path('input_reconstruction.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
