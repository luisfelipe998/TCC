"""Known-motion controls for displacement sign, interpolation and rejection."""
from pathlib import Path
import json
import numpy as np
from scipy.ndimage import gaussian_filter1d, shift
from motion_core import (ncc, estimate, longest_segment, reconstruct_segment,
                         SEARCH)

ROOT = Path(__file__).resolve().parents[1]


def run():
    rng = np.random.default_rng(20261024)
    a = gaussian_filter1d(rng.normal(size=1016), 2)
    a = 25 * a / a.std()
    translations = []
    for truth in [-120., -52.3, -12.25, -6., 0., .7, 7.4, 53.2, 120.]:
        b = shift(a, truth, order=3, mode='constant')
        r = ncc(a, b)
        error = abs(r['dx'] - truth)
        assert r['passed'] and error < .08, (truth, r, error)
        integer = int(round(r['dx']))
        direct = np.corrcoef(a[SEARCH:-SEARCH],
                             b[SEARCH + integer:len(a) - SEARCH + integer])[0, 1]
        assert abs(direct - r['peak']) < 1e-10
        translations.append(dict(true_dx=truth, estimated_dx=r['dx'], error_px=error))
    assert not ncc(np.ones(1016), np.ones(1016))['passed']
    periodic = 25 * np.sin(2 * np.pi * np.arange(1016) / 40)
    assert not ncc(periodic, shift(periodic, -8., mode='wrap'))['passed']

    # Known advected texture under a stronger, stationary optical background.
    speed = -12.25; count = 41
    large = gaussian_filter1d(rng.normal(size=3000), 3)
    large = 15 * large / large.std()
    stationary = gaussian_filter1d(rng.normal(size=1016), 5)
    stationary = 50 * stationary / stationary.std()
    cube = np.stack([np.tile(128 + stationary + np.interp(
        np.arange(1016) - speed * i + 300, np.arange(3000), large), (80, 1))
        + rng.normal(0, .2, (80, 1016)) for i in range(count)])
    records, _ = estimate(cube)
    good = np.array([r['new_dx'] for r in records if r['supported']])
    assert len(good) >= 30 and np.max(abs(good - speed)) < .35
    static = np.repeat(cube[:1], 12, axis=0)
    stationary_records, _ = estimate(static)
    assert not any(r['supported'] for r in stationary_records)

    # Linear intensities allow an exact geometric interpolation check in both
    # directions, independently of the NCC estimator.
    reconstruction = []
    for dx in [np.array([-3.2, -5.7, -2.4, -6.1]),
               np.array([3.2, 5.7, 2.4, 6.1])]:
        pos = np.r_[0., np.cumsum(-dx)]
        frames = np.stack([np.tile(25 + .1 * (np.arange(1016) + s), (80, 1)) for s in pos])
        image, coordinates = reconstruct_segment(frames, dx)
        expected = np.rint(25 + .1 * (507.5 + coordinates['position_px'] + pos.min()))
        assert np.array_equal(image, np.tile(expected, (80, 1)))
        assert coordinates['position_px'][-1] < abs(pos[-1])
        assert np.isclose(coordinates['total_travel_px'], abs(dx.sum()))
        reconstruction.append(dict(direction='left' if dx[0] < 0 else 'right',
                                   pixel_error=0, total_travel_px=abs(dx.sum())))
    synthetic_records = [dict(supported=ok, new_dx=-2. if ok else np.nan)
                         for ok in [True, True, False, True, True, True, False]]
    assert longest_segment(synthetic_records) == (3, 6)
    result = dict(status='PASS', translations=translations,
                  maximum_translation_error_px=max(r['error_px'] for r in translations),
                  constant_and_ambiguous_periodic_profiles_rejected=True,
                  stationary_sequence_rejected=True,
                  static_background_sequence=dict(true_dx=speed, supported=len(good),
                                                  intervals=count - 1,
                                                  max_error_px=float(max(abs(good - speed)))),
                  reconstruction=reconstruction, unsupported_gaps_not_bridged=True,
                  limitation='Known translations verify implementation, not real-flow measurement accuracy.')
    (ROOT / 'results/numerical_validation.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


if __name__ == '__main__':
    print(json.dumps(run(), indent=2))
