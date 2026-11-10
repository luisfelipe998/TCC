#!/usr/bin/env python3
"""Challenge image-based liquid inference; none of these rows is a measurement."""
from pathlib import Path
import hashlib
import json
import math
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
SRC2 = BASE.parent
OLD = SRC2/'v15_gas_liquid_velocity'
PHASE = OLD/'phase_conditioned_motion'
DATA = SRC2/'data_processed'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def geometry(mask):
    """For each frame/x, sum squared lengths of disjoint radial gas runs."""
    counts = np.zeros((mask.shape[0], mask.shape[2]), dtype=np.int16)
    square_sum = np.zeros_like(counts, dtype=np.int32)
    for y in range(18, 71):
        active = mask[:, y, :]
        square_sum += np.where(active, 0, counts.astype(np.int32)**2)
        counts = np.where(active, counts+1, 0)
    square_sum += counts.astype(np.int32)**2
    return square_sum / 53**2


def run():
    inventory = pd.read_csv(PHASE/'results/native_phase_inventory.csv').set_index('run')
    optical = pd.read_csv(OLD/'frame_motion_evaluation/results/temporal_per_run.csv').set_index('run')
    flow = pd.read_csv(OLD/'results/massflow_per_run.csv').set_index('run')
    model = json.loads((OLD/'results/validation/model_metadata.json').read_text())
    props = model['properties'] if 'properties' in model else model
    # Existing property values are explicit assumptions at 27 C, not new measurements.
    if 'rho_l' not in props:
        props = next(value for value in model.values() if isinstance(value,dict) and 'rho_l' in value)
    rho_l, rho_g = float(props['rho_l']), float(props['rho_g'])
    area = math.pi*.001**2/4
    geometries, trials, sources = [], [], []
    for name, inv in inventory.iterrows():
        digest = hashlib.sha256()
        files = sorted((DATA/name).glob('ImgA*.jpg'))
        assert len(files) == 501
        for path in files:
            digest.update(path.name.encode())
            digest.update(hashlib.sha256(path.read_bytes()).digest())
        assert digest.hexdigest() == inv.raw_frame_sha256, name
        path = PHASE/inv.packed_file
        assert sha(path) == inv.packed_sha256, name
        saved = np.load(path, allow_pickle=False)
        assert tuple(saved['mask_shape']) == (501,80,1016)
        assert tuple(saved['lumen_rows']) == (18,71)
        proxy_px_s = float(optical.loc[name,'median_signed_px_s'])
        for suffix in ('060','070','080'):
            mask = np.unpackbits(saved['masks_'+suffix],axis=2,count=1016,bitorder='big').astype(bool)
            alpha_x = geometry(mask)
            alpha = float(alpha_x.mean())
            fraction = float(mask[:,18:71].mean())
            assert np.isclose(fraction, inv['masks_'+suffix+'_lumen_fraction'], atol=1e-14,rtol=0)
            assert 0 <= alpha <= 1
            geom = dict(run=name, pattern=inv.pattern, threshold=float(suffix)/100,
                        native_circle_alpha=alpha, native_projected_fraction=fraction,
                        full_height_column_fraction=float((alpha_x==1).mean()),
                        alpha_frame_min=float(alpha_x.mean(axis=1).min()),
                        alpha_frame_max=float(alpha_x.mean(axis=1).max()),
                        high_saturation=bool(optical.loc[name,'high_saturation']))
            geometries.append(geom)
            for flow_key in ('mean','p05','p95'):
                G = float(flow.loc[name,'massflow_'+flow_key+'_g_s'])*1e-3/area
                for scale_multiplier in (.9,1.,1.1):
                    scale = .001/53*scale_multiplier
                    assumed_ug = abs(proxy_px_s)*scale
                    denominator = rho_l*(1-alpha)
                    liquid = (G-rho_g*alpha*assumed_ug)/denominator if denominator>0 else math.nan
                    jg = alpha*assumed_ug
                    jl = (1-alpha)*liquid
                    residual = rho_g*jg+rho_l*jl-G
                    reasons = ['unverified_material_identity','conditional_native_geometry',
                               'condition_only_massflow','product_of_nonmatched_means',
                               'assumed_axial_scale_and_properties']
                    if not bool(optical.loc[name,'base_coverage_and_saturation_pass']):
                        reasons.append('low_optical_support_or_saturation')
                    if alpha >= .98:
                        reasons.append('small_liquid_area_denominator')
                    if alpha == 1:
                        reasons.append('singular_alpha')
                    if math.isfinite(liquid) and liquid<0:
                        reasons.append('negative_liquid_inference')
                    if math.isfinite(residual):
                        assert abs(residual) < 1e-9
                    trials.append(dict(**geom, massflow_statistic=flow_key, G_kg_m2_s=G,
                                       optical_proxy_signed_px_s=proxy_px_s, axial_scale_multiplier=scale_multiplier,
                                       assumed_axial_scale_m_px=scale, substituted_ug_m_s=assumed_ug,
                                       inferred_ul_m_s=liquid, implied_Jg_m_s=jg, implied_Jl_m_s=jl,
                                       conservation_residual_kg_m2_s=residual, rho_l_kg_m3=rho_l,rho_g_kg_m3=rho_g,
                                       accepted_as_phase_velocity=False, limitation_flags=';'.join(reasons)))
        saved.close()
        sources.append(dict(run=name, raw_frames=len(files),raw_frames_sha256=digest.hexdigest(),
                            mask_file=str(path.relative_to(SRC2)),mask_sha256=inv.packed_sha256))
        print(name, 'native masks and source hashes checked', flush=True)
    geometry_frame = pd.DataFrame(geometries)
    frame = pd.DataFrame(trials)
    geometry_frame.to_csv(BASE/'results/native_geometry_stress.csv',index=False)
    frame.to_csv(BASE/'results/mass_balance_stress.csv',index=False)
    baseline = frame[(frame.threshold==.7)&(frame.massflow_statistic=='mean')&(frame.axial_scale_multiplier==1)]
    threshold_only = frame[(frame.massflow_statistic=='mean')&(frame.axial_scale_multiplier==1)]
    summary = []
    for name, group in threshold_only.groupby('run',sort=False):
        values = group.inferred_ul_m_s.dropna()
        summary.append(dict(run=name,pattern=group.pattern.iloc[0],
                            alpha_min=float(group.native_circle_alpha.min()),
                            alpha_max=float(group.native_circle_alpha.max()),
                            candidate_liquid_min_m_s=float(values.min()),
                            candidate_liquid_max_m_s=float(values.max()),
                            threshold_max_min_factor=float(values.max()/values.min()) if len(values) and values.min()>0 else math.nan,
                            any_singular_alpha=bool((group.native_circle_alpha==1).any()),
                            accepted_as_phase_velocity=False))
    pd.DataFrame(summary).to_csv(BASE/'results/mass_balance_stress_per_run.csv',index=False)
    metadata = dict(status='arithmetic_checked_only',native_frames_hash_checked=26*501,
                    native_masks_reused=26*501*3,geometry_rows=len(geometries),stress_rows=len(trials),
                    accepted_physical_phase_velocity_pairs=0,
                    baseline_candidate_liquid_min_m_s=float(baseline.inferred_ul_m_s.min()),
                    baseline_candidate_liquid_max_m_s=float(baseline.inferred_ul_m_s.max()),
                    baseline_alpha_min=float(baseline.native_circle_alpha.min()),
                    baseline_alpha_max=float(baseline.native_circle_alpha.max()),
                    max_threshold_liquid_factor=max(r['threshold_max_min_factor'] for r in summary),
                    singular_rows=int(frame.inferred_ul_m_s.isna().sum()),
                    negative_liquid_rows=int((frame.inferred_ul_m_s<0).sum()),
                    formula='ul=(G-rho_g*alpha*substituted_ug)/(rho_l*(1-alpha))',
                    optical_proxy='Absolute value of median signed whole-profile image motion; unverified substitution for intrinsic gas speed.',
                    averaging_limitation='Native circle alpha is a space-time mean; selected optical median and condition-level G are not a matched phase flux measurement.',
                    scale_variants='0.9/1.0/1.1 times D/53 are hypothetical sensitivity, not calibration uncertainty.',
                    flow_variants='Mean/p05/p95 are sample distribution statistics from a condition; they are not synchronized measurements or confidence bounds.',
                    sources=sources,script_sha256=sha(Path(__file__)))
    (BASE/'results/mass_balance_stress_checks.json').write_text(json.dumps(metadata,indent=2)+chr(10))
    print(json.dumps({k:v for k,v in metadata.items() if k not in ('sources',)},indent=2))


if __name__=='__main__':
    run()
