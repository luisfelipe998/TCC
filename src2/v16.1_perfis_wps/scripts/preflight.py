"""Auditoria somente de leitura das fontes e medição antes do lote."""
from pathlib import Path
import hashlib
import json
import time
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from PIL import Image
from profile_core import ROOT, BASE, LAM, wc, zero_profile, legacy_zero


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    out = ROOT/'results'
    out.mkdir(exist_ok=True)
    inv = pd.read_csv(BASE/'results/input_audit.csv')
    assert len(inv) == 31 and inv.included.sum() == 26
    assert set(inv.loc[~inv.included, 'run']) == {'m03q01_001','m03q80_001','m04q05_001','m05q05_002','m07q10_001'}
    manifest = json.loads((BASE/'results/artifact_manifest.json').read_text())
    snapshot = {str(f.relative_to(BASE)): sha(f) for f in BASE.rglob('*') if f.is_file()}
    failures = []
    checked = 0
    for group in ['output_sha256', 'script_sha256']:
        for rel, expected in manifest[group].items():
            if sha(BASE/rel) != expected:
                failures.append(rel)
            checked += 1
    assert not failures, f'Manifesto da base divergiu: {failures}'
    rows = []
    inputs = {}
    for r in inv.itertuples():
        folder = ROOT.parent/'data_processed'/r.run
        for name, expected in [('__panorama_v11.png', r.panorama_sha256), ('__dados_cropped_v11.xlsx', r.xlsx_sha256)]:
            path = folder/name
            assert sha(path) == expected, path
            inputs[str(path.relative_to(ROOT.parent))] = expected
        with Image.open(folder/'__panorama_v11.png') as im:
            pano = np.asarray(im.convert('L'))
        with np.load(BASE/'intermediate'/f'{r.run}_input.npz') as z:
            error = float(np.max(np.abs(z['raw']-pano[18:72].mean(axis=0))))
            assert error == 0 and np.array_equal(z['raw_20_59'], pano[20:60].mean(axis=0))
            assert hashlib.sha256(z['raw'].astype('<f8').tobytes()).hexdigest() == r.raw_array_sha256
            assert np.count_nonzero(z['raw_shifts'] <= 0) == r.nonpositive_shifts
        labels = [p.stem.lstrip('_') for p in folder.glob('*.txt') if p.stem.lstrip('_') != 'dados']
        assert labels == [r.pattern], (r.run, labels, r.pattern)
        rows.append(dict(aquisicao=r.run, incluida=bool(r.included), n=r.n,
                         erro_extracao=error, deslocamentos_nao_positivos=r.nonpositive_shifts))
    (out/'source_snapshot.json').write_text(json.dumps(snapshot, indent=2)+'\n')
    (out/'input_hashes.json').write_text(json.dumps(inputs, indent=2)+'\n')
    pd.DataFrame(rows).to_csv(out/'input_verification.csv', index=False)
    trials = []
    eligible = inv[inv.included].sort_values('n')
    for r in [eligible.iloc[0], eligible.iloc[len(eligible)//2], eligible.iloc[-1]]:
        raw = np.load(BASE/'intermediate'/f'{r.run}_input.npz')['raw']
        x = wc.prepare(raw)
        t = time.perf_counter()
        full = wc.morlet_cwt(wc.correlation(x,x), LAM/wc.FOURIER_FACTOR)
        elapsed = time.perf_counter()-t
        t = time.perf_counter()
        fast = zero_profile(x)
        fast_elapsed = time.perf_counter()-t
        err = float(np.max(np.abs(fast-np.abs(full[:,len(x)-1]))))
        assert np.allclose(fast, np.abs(full[:,len(x)-1]), rtol=1e-10, atol=1e-11)
        size = full.nbytes
        del full
        t = time.perf_counter()
        legacy_zero(x)
        trials.append(dict(n=len(x), run=r.run, full_seconds=elapsed, zero_seconds=fast_elapsed,
                           legacy_seconds=time.perf_counter()-t, zero_max_abs_error=err,
                           complex_map_bytes=size))
    report = dict(status='PASS', created_utc=datetime.now(timezone.utc).isoformat(),
                  source_files=len(snapshot), manifest_entries_checked=checked,
                  audited_inputs=31, eligible_runs=26, extraction_rows=[18,71],
                  original_reconstruction_reexecuted=False, benchmark=trials)
    (out/'preflight.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
