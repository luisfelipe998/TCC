#!/usr/bin/env python3
"""Conditional conservation bounds for a rigidly translating isolated bubble cell.

This module does not infer experimental velocities. Its controls are mathematical
examples in dimensionless units. In particular, their parameters are NOT R600a
properties or measurements from the recordings.
"""
from pathlib import Path
import hashlib
import json
import numpy as np

BASE = Path(__file__).resolve().parents[1]


def slug_ratio_bounds(alpha_body_lower):
    """Bounds on U_s/c if 0 <= U_f <= c at a bubble-body section.

alpha_body_lower must be a valid lower bound on the true gas AREA fraction
at that section, not a projected mask fraction. No unit cell is needed for
this local liquid-slug relation; an all-liquid section and rigid traveling
shape connected to that body section ARE needed.
    """
    a = float(alpha_body_lower)
    if not np.isfinite(a) or not 0 <= a < 1:
        raise ValueError('alpha_body_lower must be finite and in [0, 1)')
    return a, 1.0


def cell_ratio_bounds(alpha_body_lower, alpha_cell_upper):
    """Conservative bounds on the intrinsic cell mean liquid velocity / c.

If exact maximum gas area is known, it is the strongest eligible body
section. A lower body fraction also yields a weaker valid bound, PROVIDED
the same co-current film hypothesis holds at that section. The mean gas
fraction must refer to a COMPLETE translating unit cell, including caps
and its associated full-liquid slug. No phi*body approximation is made.
    """
    a, _ = slug_ratio_bounds(alpha_body_lower)
    b = float(alpha_cell_upper)
    if not np.isfinite(b) or not 0 <= b < 1:
        raise ValueError('alpha_cell_upper must be finite and in [0, 1)')
    return max(0.0, (a - b) / (1.0 - b)), 1.0


def gas_liquid_slip_bounds(liquid_ratio_lower):
    """Bounds on c/U_l for the SAME eligible unit cell.

Infinity expresses an unbounded conditional upper end, never a result
assigned to an experimental recording.
    """
    r = float(liquid_ratio_lower)
    if not np.isfinite(r) or not 0 <= r <= 1:
        raise ValueError('liquid_ratio_lower must be finite and in [0, 1]')
    return 1.0, 1.0 / r if r > 0 else float('inf')


def inter_run_liquid_ratio_bounds(c_ratio_lower, c_ratio_upper, r_i_lower, r_j_lower):
    """Bounds on U_liquid_i/U_liquid_j given c_i/c_j and cell bounds.

This requires physically comparable gas velocities and calibrated relative
axial scales (or an explicitly declared equal-magnification assumption).
    """
    lo, hi, ri, rj = map(float, (c_ratio_lower, c_ratio_upper, r_i_lower, r_j_lower))
    if not all(np.isfinite(v) for v in (lo, hi, ri, rj)):
        raise ValueError('all inputs must be finite')
    if not (0 < lo <= hi and 0 <= ri <= 1 and 0 <= rj <= 1):
        raise ValueError('invalid positive c ratio or liquid-ratio bounds')
    return lo * ri, hi / rj if rj > 0 else float('inf')


def controls():
    rng = np.random.default_rng(152)
    max_conservation_residual = 0.0
    max_mean_residual = 0.0
    max_film_residual = 0.0
    n_cell = 0
    x = (np.arange(8192) + .5) / 8192
    # Smooth caps, finite plateau and a truly gas-free liquid slug.
    profile = np.zeros_like(x)
    rise = (x >= .15) & (x < .25)
    plateau = (x >= .25) & (x <= .65)
    fall = (x > .65) & (x <= .75)
    profile[rise] = np.sin(np.pi * (x[rise] - .15) / .2) ** 2
    profile[plateau] = 1.0
    profile[fall] = np.cos(np.pi * (x[fall] - .65) / .2) ** 2
    for a_body in (.1, .35, .6, .85, .98):
        alpha = a_body * profile
        alpha_mean = float(alpha.mean())
        c = 2.3  # Dimensionless test value only.
        for f in (0., .2, .5, .8, 1.):
            u_f = f * c
            u_s = a_body * c + (1 - a_body) * u_f
            u_l = (u_s - c * alpha) / (1 - alpha)
            j_g = alpha * c
            j_l = (1 - alpha) * u_l
            weighted_l = float(j_l.mean() / (1 - alpha_mean))
            exact_l = (u_s - c * alpha_mean) / (1 - alpha_mean)
            lo, hi = cell_ratio_bounds(a_body, alpha_mean)
            assert lo - 1e-13 <= weighted_l / c <= hi + 1e-13
            assert np.min(u_l) >= -1e-13
            assert np.max(u_l) <= c + 1e-13
            assert np.allclose(u_l[plateau], u_f, atol=1e-13)
            # For a rigid traveling field, the liquid flux relative to c is constant.
            residual = float(np.max(np.abs((1 - alpha) * (u_l - c) - (u_s - c))))
            max_conservation_residual = max(max_conservation_residual, residual,
                                            float(np.max(np.abs(j_g + j_l - u_s))))
            max_mean_residual = max(max_mean_residual, abs(weighted_l - exact_l))
            n_cell += 1
    # Test interval monotonicity without asserting experimental uncertainty.
    for _ in range(2000):
        a = rng.uniform(.01, .999)
        b = rng.uniform(0, a)
        f = rng.uniform(0, 1)
        ratio = (a + (1 - a) * f - b) / (1 - b)
        a_lo = max(0., a - rng.uniform(0, .2))
        b_hi = min(.99999, b + rng.uniform(0, .2))
        lo, hi = cell_ratio_bounds(a_lo, b_hi)
        assert lo - 1e-13 <= ratio <= hi + 1e-13
    # Independent inter-run interval arithmetic.
    for _ in range(2000):
        c = rng.uniform(.05, 20)
        r_i, r_j = rng.uniform(.001, 1, 2)
        r_i_lo, r_j_lo = r_i * rng.random(), r_j * rng.random()
        lo, hi = inter_run_liquid_ratio_bounds(.9*c, 1.1*c, r_i_lo, r_j_lo)
        assert lo <= c * r_i / r_j <= hi
    # The often suggested c/2--c/3 film interval needs nonnegative,
    # thickness-independent shear and pressure drives in the planar model.
    for shear_fraction in np.linspace(0, 1, 101):
        # Contributions to c sum to one: tau*h/mu + K*h^2/mu = 1.
        shear_term = shear_fraction
        pressure_term = 1 - shear_fraction
        u_over_c = shear_term / 2 + pressure_term / 3
        exact = 1/3 + shear_fraction / 6
        max_film_residual = max(max_film_residual, abs(u_over_c - exact))
        assert 1/3 - 1e-14 <= u_over_c <= 1/2 + 1e-14
    # Counterexamples are deliberate violations of the branch assumption,
    # proving that the interval is conditional. They are NOT real run data.
    a = .8
    counterexamples = []
    for f in (-.5, 1.5):
        ratio = a + (1-a) * f
        assert ratio < a if f < 0 else ratio > 1
        counterexamples.append(dict(film_over_bubble=f, liquid_slug_over_bubble=ratio,
                                     violates='co-current film' if f < 0 else 'film no faster than bubble'))
    return dict(
        scope='Mathematical controls only; no measured phase velocities or R600a properties.',
        traveling_unit_cell_controls=n_cell,
        interval_monotonicity_controls=2000,
        inter_run_interval_controls=2000,
        planar_film_controls=101,
        max_conservation_residual=max_conservation_residual,
        max_cell_mean_residual=max_mean_residual,
        max_planar_film_residual=max_film_residual,
        conditional_branch_counterexamples=counterexamples,
        required_assumptions=[
            'Rigidly translating isolated gas bubble and adjoining full-liquid slug.',
            'Constant phase densities, negligible phase change, no coalescence or breakup.',
            'Valid gas cross-sectional AREA fractions, not projected image fractions.',
            'At the nominated bubble-body section: 0 <= area-mean film velocity <= bubble celerity.',
            'Complete unit cell and its own mean gas fraction for the intrinsic liquid mean.',
            'Inter-run physical speed ratios require a relative axial optical scale.'
        ],
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


if __name__ == '__main__':
    result = controls()
    out = BASE / 'results' / 'physical_bounds_controls.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(result, indent=2, allow_nan=False))
