"""Independent saved-output checks. Does not import the motion implementation."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT.parent
DATA = BASE.parent / 'data_processed'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def direct_ncc(a, b):
    # Explicit candidate windows, without the FFT/prefix-sum implementation.
    t = a[200:-200].copy(); t -= t.mean()
    windows = np.lib.stride_tricks.sliding_window_view(b, len(t))
    centered = windows - windows.mean(axis=1, keepdims=True)
    denominator = np.sqrt(np.sum(t * t) * np.sum(centered * centered, axis=1))
    values = np.full(len(windows), np.nan)
    np.divide(centered @ t, denominator, out=values, where=denominator > 1e-9)
    if not np.isfinite(values).any(): return np.nan, np.nan
    k = int(np.nanargmax(values)); sub = 0.
    if 0 < k < len(values) - 1:
        curvature = values[k - 1] - 2 * values[k] + values[k + 1]
        if np.isfinite(curvature) and curvature < -1e-12:
            sub = float(np.clip(.5 * (values[k - 1] - values[k + 1]) / curvature, -.5, .5))
    return k - 200 + sub, values[k]


def main():
    inventory = pd.read_csv(BASE / 'results/input_audit.csv').set_index('run')
    pairs = pd.read_csv(ROOT / 'results/displacements.csv')
    temporal = pd.read_csv(ROOT / 'results/temporal_checks.csv')
    summary = pd.read_csv(ROOT / 'results/run_summary.csv').set_index('run')
    assert set(summary.index) == set(inventory.index) and len(summary) == 31
    assert len(pairs) == 31 * 500 and len(temporal) == 31 * 499
    assert not pairs.duplicated(['run', 'frame_i']).any()
    assert not temporal.duplicated(['run', 'frame_i']).any()
    finite_new = np.isfinite(pairs.new_dx)
    assert np.array_equal(finite_new, pairs.supported)
    assert np.allclose(pairs.loc[finite_new, 'new_dx'], pairs.loc[finite_new, 'candidate_dx'])
    assert np.allclose(pairs.agreement_tolerance, np.maximum(1., .15 * abs(pairs.candidate_dx)), equal_nan=True)
    quality = ((pairs.fit_peak >= .6) & (pairs.peak_margin >= .03)
               & (pairs.texture_std >= 1) & ~pairs.search_boundary)
    assert np.array_equal(quality, pairs.fit_pass)
    base = (quality & pairs.reverse_pass & pairs.check_rows_pass
            & (pairs.row_difference <= pairs.agreement_tolerance)
            & (pairs.reverse_sum <= pairs.agreement_tolerance))
    assert np.array_equal(base, pairs.base_supported)
    expected_support = base & (pairs.temporal_checks > 0) & (pairs.temporal_failures == 0)
    assert np.array_equal(expected_support, pairs.supported)
    largest_ncc_error = 0.; largest_peak_error = 0.; sampled_pairs = 0; panoramas = 0
    for run, group in pairs.groupby('run'):
        group = group.sort_values('frame_i').reset_index(drop=True)
        assert np.array_equal(group.frame_i, np.arange(500))
        assert np.array_equal(group.frame_j, np.arange(1, 501))
        assert group.pattern.eq(inventory.loc[run, 'pattern']).all()
        assert group.included.eq(inventory.loc[run, 'included']).all()
        old = np.load(BASE / 'intermediate' / f'{run}_input.npz')
        assert np.array_equal(group.old_raw_s, old['raw_shifts'])
        assert np.array_equal(group.old_dx, -old['shifts'][:500])
        assert np.array_equal(group.old_clamped, old['raw_shifts'] <= 0)
        checks = np.zeros(500, dtype=int); failures = np.zeros(500, dtype=int)
        for row in temporal[temporal.run == run].itertuples():
            i = row.frame_i
            total = group.candidate_dx.iloc[i:i + 2].sum(skipna=False)
            valid = bool(group.base_supported.iloc[i:i + 2].all() and row.direct_pass)
            assert np.isclose(total, row.adjacent_sum_dx, equal_nan=True)
            assert np.isclose(abs(row.direct_dx - total), row.absolute_difference, equal_nan=True)
            assert bool(row.comparable) == valid
            assert bool(row.agrees) == bool(valid and row.absolute_difference <= row.tolerance)
            if valid:
                checks[i:i + 2] += 1
                if not row.agrees: failures[i:i + 2] += 1
        assert np.array_equal(checks, group.temporal_checks)
        assert np.array_equal(failures, group.temporal_failures)
        supported = group[group.supported]
        r = summary.loc[run]
        assert len(supported) == r.supported_count
        assert np.isclose(len(supported) / 500, r.supported_fraction)
        assert int(group.old_clamped.sum()) == r.old_clamped_count
        assert int((group.old_clamped & group.supported).sum()) == r.supported_old_clamps
        assert np.isclose(supported.old_dx.median(), r.old_median_on_supported, equal_nan=True)
        assert np.isclose(supported.new_dx.median(), r.new_supported_median, equal_nan=True)
        assert np.isclose(supported.dynamic_correlation_change.mean(), r.mean_check_correlation_change, equal_nan=True)
        delta = supported.dynamic_correlation_change
        assert int((delta > 1e-6).sum()) == r.alignment_improved
        assert int((delta < -1e-6).sum()) == r.alignment_worsened
        files = sorted((DATA / run).glob('ImgA*.jpg'))
        frames = np.stack([np.asarray(Image.open(f).convert('L')) for f in files])
        profiles = frames[:, np.arange(18, 72, 2)].mean(axis=1)
        profiles -= profiles.mean(axis=0, keepdims=True)
        indices = {0, 250, 499}
        clipped = np.flatnonzero(old['raw_shifts'] <= 0)
        if len(clipped): indices.add(int(clipped[0]))
        for i in sorted(indices):
            dx, peak = direct_ncc(profiles[i], profiles[i + 1])
            saved = group.iloc[i]
            assert np.isclose(dx, saved.candidate_dx, rtol=0, atol=1e-7, equal_nan=True), (run, i, dx, saved.candidate_dx)
            assert np.isclose(peak, saved.fit_peak, rtol=0, atol=1e-10, equal_nan=True)
            if np.isfinite(dx): largest_ncc_error = max(largest_ncc_error, abs(dx - saved.candidate_dx))
            if np.isfinite(peak): largest_peak_error = max(largest_peak_error, abs(peak - saved.fit_peak))
            sampled_pairs += 1
        if r.segment_intervals:
            first, last = int(r.segment_first), int(r.segment_last)
            assert last - first == r.segment_intervals and group.supported.iloc[first:last].all()
            with np.load(ROOT / 'intermediate' / f'{run}_segment_coordinates.npz') as z:
                dx = group.new_dx.iloc[first:last].to_numpy()
                positions = np.r_[0., np.cumsum(-dx)]
                assert (np.diff(positions) > 0).all() or (np.diff(positions) < 0).all()
                order = np.argsort(positions); pos = positions[order] - positions.min()
                grid = np.arange(0., pos[-1], 1.)
                owners_sorted = np.searchsorted((pos[:-1] + pos[1:]) / 2, grid, side='right')
                owners = order[owners_sorted]
                columns = 507.5 + grid - pos[owners_sorted]
                assert np.array_equal(z['source_frame'], owners)
                assert np.array_equal(z['position_px'], grid)
                assert np.allclose(z['source_column'], columns, rtol=0, atol=1e-8)
                lo = np.floor(z['source_column']).astype(int); hi = np.minimum(lo + 1, 1015)
                part = z['source_column'] - lo
                rebuilt = np.rint(np.clip(frames[first + owners, :, lo].T * (1 - part)
                                         + frames[first + owners, :, hi].T * part, 0, 255)).astype('uint8')
            before = np.asarray(Image.open(ROOT / 'panoramas' / f'{run}_before.png'))
            after = np.asarray(Image.open(ROOT / 'panoramas' / f'{run}_after.png'))
            offsets = np.r_[0, np.cumsum(old['shifts'])]
            original = np.asarray(Image.open(DATA / run / '__panorama_v11.png'))
            assert np.array_equal(before, original[:, offsets[first]:offsets[last + 1]])
            assert np.array_equal(after, rebuilt)
            assert before.shape[1] == r.old_segment_width_px and after.shape[1] == r.new_segment_width_px
            for label, image in [('before', before), ('after', after)]:
                signal = pd.read_csv(ROOT / 'signals' / f'{run}_{label}.csv')
                assert np.array_equal(signal.position_px, np.arange(image.shape[1]))
                assert np.allclose(signal.gray_mean_rows_18_71, image[18:72].sum(axis=0, dtype=np.uint64) / 54,
                                   rtol=0, atol=1e-12)
            panoramas += 1
        assert (ROOT / 'figures' / f'{run}_comparison.png').exists()
    for path in (ROOT / 'figures').glob('*.png'):
        with Image.open(path) as im: im.verify()
    frozen = json.loads((ROOT / 'results/parent_snapshot.json').read_text())
    for name, digest in frozen.items():
        assert sha(BASE / name) == digest, f'Parent artifact modified: {name}'
    current = {str(p.relative_to(BASE)) for p in BASE.rglob('*')
               if p.is_file() and ROOT not in p.parents and '__pycache__' not in p.parts}
    assert current == set(frozen), 'Unexpected artifact change outside displacement/.'
    numerical = json.loads((ROOT / 'results/numerical_validation.json').read_text())
    assert numerical['status'] == 'PASS'
    result = dict(status='PASS', runs=31, intervals=15500,
                  temporal_checks=15469, independent_direct_ncc_pairs=sampled_pairs,
                  maximum_independent_shift_difference=largest_ncc_error,
                  maximum_independent_peak_difference=largest_peak_error,
                  panorama_pairs_verified=panoramas, exact_old_panorama_crops=True,
                  exact_new_panorama_reconstructions=True, exact_extraction_rows=[18, 71],
                  rejected_intervals_are_nan=True, parent_artifacts_unchanged=len(frozen),
                  same_labels_and_inclusion_as_parent=True,
                  scope='Implementation, saved arithmetic, construction and provenance; not physical ground-truth accuracy.')
    (ROOT / 'results/verification.json').write_text(json.dumps(result, indent=2) + '\n')
    manifest = {str(p.relative_to(ROOT)): sha(p) for p in sorted(ROOT.rglob('*'))
                if p.is_file() and p.name != 'artifact_manifest.json' and '__pycache__' not in p.parts}
    (ROOT / 'results/artifact_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__': main()
