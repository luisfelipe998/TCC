#!/usr/bin/env python3
"""Shared-motion information in complementary masks, not experimental slip."""
from pathlib import Path
import json
import numpy as np

BASE=Path(__file__).resolve().parents[1]


def direct_correlations(a,b):
    template=a[50:250]
    centered=template-template.mean()
    result=[]
    for displacement in range(-40,41):
        candidate=b[50+displacement:250+displacement]
        candidate=candidate-candidate.mean()
        result.append(float(np.dot(centered,candidate)/np.sqrt(np.dot(centered,centered)*np.dot(candidate,candidate))))
    return np.array(result)


def run():
    gas_a=np.zeros(300)
    for lo,hi in [(40,90),(140,167),(210,225)]:
        gas_a[lo:hi]=1
    truth=17
    gas_b=np.zeros_like(gas_a)
    gas_b[truth:]=gas_a[:-truth]
    liquid_a=1-gas_a
    liquid_b=1-gas_b
    cgas=direct_correlations(gas_a,gas_b)
    cliquid=direct_correlations(liquid_a,liquid_b)
    shifts=np.arange(-40,41)
    gas_shift=int(shifts[np.argmax(cgas)])
    liquid_shift=int(shifts[np.argmax(cliquid)])
    discrepancy=float(np.max(np.abs(cgas-cliquid)))
    assert gas_shift==liquid_shift==truth
    assert discrepancy<1e-12
    result=dict(control='Complementary binary masks of the same translating interfaces',
        true_interface_shift_px=truth,gas_mask_estimated_shift_px=gas_shift,
        liquid_complement_estimated_shift_px=liquid_shift,
        inferred_mask_motion_ratio=gas_shift/liquid_shift,
        maximum_correlation_curve_difference=discrepancy,passed=True,
        interpretation='Both labels contain the same interface-motion information. Ratio=1 is imposed by complementary masks and does not identify either material phase velocity.',
        limitation='This control does not predict experimental conditional-grayscale ratios; those can differ because they sample different intervals/features.')
    (BASE/'results').mkdir(parents=True,exist_ok=True)
    (BASE/'results'/'mask_complement_control.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    run()
