#!/usr/bin/env python3
"""Independent identifiability controls; no experimental phase speed is inferred.

Reads the closed dataset and prior saved results. Writes only adversarial_*
artifacts beside this script. Synthetic values are explicitly distinguished.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.fft import irfft, rfft
from scipy.ndimage import gaussian_filter1d

BASE = Path(__file__).resolve().parents[1]
V151 = BASE.parent
SRC2 = V151.parent
V15 = SRC2 / 'v15_gas_liquid_velocity'
PHASE = V15 / 'phase_conditioned_motion'
OUT = BASE / 'results'
SEED = 20261024
SOURCES = {}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    SOURCES[str(path.relative_to(SRC2))] = sha(path)
    return pd.read_csv(path)


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean(v) for v in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    return value


def save_csv(stem, records):
    frame = pd.DataFrame(records)
    frame.to_csv(OUT / ('adversarial_' + stem + '.csv'), index=False)
    return frame


def ncc_curve(a, b, start=380, width=256, maxshift=96):
    """Independent direct centered NCC over identical pixel supports."""
    t = np.asarray(a[start:start + width], dtype=float)
    t -= t.mean()
    windows = np.lib.stride_tricks.sliding_window_view(
        b[start-maxshift:start+width+maxshift], width).copy()
    windows -= windows.mean(axis=1, keepdims=True)
    denominator = np.sqrt(np.sum(t*t) * np.sum(windows*windows, axis=1))
    return np.divide(windows @ t, denominator,
                     out=np.full(2*maxshift+1, np.nan), where=denominator>1e-12)


def complement_control(inventory):
    records = []
    indices = np.linspace(0, 499, 10, dtype=int)
    for inv in inventory.itertuples():
        path = PHASE / 'intermediate/masks' / (inv.run + '.npz')
        SOURCES[str(path.relative_to(SRC2))] = sha(path)
        with np.load(path, allow_pickle=False) as stored:
            assert tuple(stored['mask_shape']) == (501,80,1016)
            for suffix in ('060','070','080'):
                selected = np.unique(np.r_[indices,indices+1])
                masks = np.unpackbits(stored['masks_'+suffix][selected,18:71],
                                      axis=-1,count=1016,bitorder='big').astype(float)
                profiles = masks.mean(axis=1)
                for first in indices:
                    a = profiles[np.flatnonzero(selected==first)[0]]
                    b = profiles[np.flatnonzero(selected==first+1)[0]]
                    cg, cl = ncc_curve(a.copy(),b.copy()), ncc_curve(1-a,1-b)
                    finite = np.isfinite(cg)&np.isfinite(cl)
                    valid = finite.any()
                    kg = int(np.nanargmax(cg)) if valid else None
                    kl = int(np.nanargmax(cl)) if valid else None
                    difference = float(np.max(abs(cg[finite]-cl[finite]))) if valid else np.nan
                    # Different indices are possible only for numerical ties.
                    assert not valid or difference < 1e-10, (inv.run,suffix,first,difference)
                    records.append(dict(run=inv.run,pattern=inv.pattern,threshold=int(suffix)/100,
                        first_frame=int(first),second_frame=int(first+1),
                        has_nonconstant_match=bool(valid),finite_lag_cells=int(finite.sum()),
                        gas_candidate_shift_px=kg-96 if valid else np.nan,
                        complement_candidate_shift_px=kl-96 if valid else np.nan,
                        identical_argmax=kg==kl if valid else False,
                        gas_loss_at_complement_max=(float(cg[kg]-cg[kl]) if valid else np.nan),
                        complement_loss_at_gas_max=(float(cl[kl]-cl[kg]) if valid else np.nan),
                        gas_lags_within_1e_10_of_max=(int(np.sum(abs(cg-cg[kg])<1e-10)) if valid else 0),
                        max_curve_difference=difference,
                        observation='complementary_mask_geometry_not_independent_phase_motion'))
        print('complement',inv.run,flush=True)
    frame = save_csv('real_mask_complement',records)
    valid = frame[frame.has_nonconstant_match]
    different = valid[~valid.identical_argmax]
    if len(different):
        assert max(different.gas_loss_at_complement_max.max(),different.complement_loss_at_gas_max.max())<1e-10
    return dict(tested_pairs=len(frame),nonconstant_pairs=len(valid),
                constant_or_degenerate_pairs=len(frame)-len(valid),
                identical_argmax_pairs=int(valid.identical_argmax.sum()),
                different_argmax_only_numerical_ties=len(different),
                unique_peak_pairs_at_1e_10=int((valid.gas_lags_within_1e_10_of_max==1).sum()),
                max_curve_difference=valid.max_curve_difference.max(),
                all_curves_agree_within_1e_10=True,
                limitation='Identity applies to normalized correlation of complementary masks over identical support, not to independent textured grayscale pixels within each phase.')


def weighted_common_motion():
    temporal = read_csv(V15/'frame_motion_evaluation/results/temporal_pairs.csv')
    fractions = read_csv(PHASE/'results/native_phase_frame_fractions.csv')
    rng = np.random.default_rng(SEED)
    records = []
    for name, group in temporal[temporal.lag_frames==1].groupby('run'):
        group = group.sort_values('first_frame_id')
        assert len(group)==500
        fr = fractions[fractions.run==name].set_index('frame_id')
        for suffix in ('060','070','080'):
            key = 'masks_'+suffix+'_lumen_fraction'
            weight = (fr.loc[group.first_frame_id,key].to_numpy()+
                      fr.loc[group.second_frame_id,key].to_numpy())/2
            velocity = group.structure_celerity_px_per_s.to_numpy()
            for gate in ('all_finite_candidates','legacy_heuristic_pass'):
                keep = np.isfinite(velocity)
                if gate=='legacy_heuristic_pass':
                    keep &= group.heuristic_pass.to_numpy(dtype=bool)
                w,v = weight[keep],velocity[keep]
                assert np.all((w>=0)&(w<=1))
                valid = len(v)>0 and w.sum()>0 and (1-w).sum()>0
                if not valid:
                    records.append(dict(run=name,pattern=group.pattern.iloc[0],threshold=int(suffix)/100,
                        gate=gate,selected_pairs=len(v),has_two_weighted_means=False))
                    continue
                gas,liquid = np.average(v,weights=w),np.average(v,weights=1-w)
                covariance = np.mean((w-w.mean())*(v-v.mean()))
                algebraic_difference = covariance/(w.mean()*(1-w.mean()))
                residual = gas-liquid-algebraic_difference
                assert abs(residual)<1e-6,(name,suffix,gate,residual)
                ratios = []
                for _ in range(199):
                    shuffled = rng.permutation(w)
                    denominator = np.average(v,weights=1-shuffled)
                    ratios.append(np.average(v,weights=shuffled)/denominator if denominator else np.nan)
                # Equal masks are not inferred experimentally: this is a control.
                constant_v = np.full(len(v),-12345.)
                constant_ratio = np.average(constant_v,weights=w)/np.average(constant_v,weights=1-w)
                assert abs(constant_ratio-1)<1e-12
                records.append(dict(run=name,pattern=group.pattern.iloc[0],threshold=int(suffix)/100,
                    gate=gate,selected_pairs=len(v),has_two_weighted_means=True,
                    mean_mask_fraction=w.mean(),gas_weight_effective_pairs=w.sum()**2/(w*w).sum(),
                    complement_weight_effective_pairs=(1-w).sum()**2/((1-w)**2).sum(),
                    weighted_gas_candidate_px_s=gas,weighted_complement_candidate_px_s=liquid,
                    signed_ratio=gas/liquid if liquid else np.nan,common_velocity_mean_px_s=v.mean(),
                    difference_px_s=gas-liquid,covariance_identity_difference_px_s=algebraic_difference,
                    covariance_identity_residual_px_s=residual,constant_motion_control_ratio=constant_ratio,
                    shuffled_ratio_median=np.nanmedian(ratios),shuffled_ratio_p025=np.nanquantile(ratios,.025),
                    shuffled_ratio_p975=np.nanquantile(ratios,.975),
                    interpretation='Two occupancy-weighted averages of one global motion signal; not independently identified phase velocities'))
    frame = save_csv('weighted_common_motion',records)
    baseline = frame[(frame.threshold==.7)&(frame.gate=='all_finite_candidates')]
    valid = frame[frame.has_two_weighted_means]
    return dict(trials=len(frame),baseline_finite_pairs=int(baseline.has_two_weighted_means.sum()),
                baseline_ratio_min=baseline.signed_ratio.min(),baseline_ratio_max=baseline.signed_ratio.max(),
                max_covariance_identity_residual_px_s=abs(valid.covariance_identity_residual_px_s).max(),
                max_shuffled_ratio_median_deviation_from_one=abs(valid.shuffled_ratio_median-1).max(),
                identity='weighted_gas-weighted_complement = population_covariance(mask_fraction,global_motion)/(mean_fraction*(1-mean_fraction))',
                shuffle_role='Break appearance/motion association; percentiles are permutation diagnostics, not experimental confidence intervals.')


def batch_ncc(sequence, lag):
    """Independently implement the original 536-px template/481-lag geometry."""
    a,b = sequence[:-lag],sequence[lag:]
    t = a[:,240:776]-a[:,240:776].mean(axis=1,keepdims=True)
    numerator = irfft(rfft(b,n=2048,axis=1)*rfft(t,n=2048,axis=1).conj(),
                       n=2048,axis=1)[:,:481]
    sums = np.pad(np.cumsum(b,axis=1),((0,0),(1,0)))
    sums2 = np.pad(np.cumsum(b*b,axis=1),((0,0),(1,0)))
    total = sums[:,536:]-sums[:,:-536]
    var = np.maximum(sums2[:,536:]-sums2[:,:-536]-total*total/536,0)
    denominator = np.sqrt(np.sum(t*t,axis=1)[:,None]*var)
    return np.divide(numerator,denominator,out=np.full_like(numerator,np.nan),where=denominator>1e-12)


def peak_shift(curve):
    k = int(np.nanargmax(curve))
    offset = 0.
    if 0<k<len(curve)-1:
        denominator = curve[k-1]-2*curve[k]+curve[k+1]
        if denominator<0:
            offset=float(np.clip(.5*(curve[k-1]-curve[k+1])/denominator,-.5,.5))
    return k-240+offset,float(curve[k])


def nonstationary_control():
    rng = np.random.default_rng(SEED+1)
    source = gaussian_filter1d(rng.normal(size=24000),2)*60
    cases = dict(constant=np.full(500,-16.),
                 two_real_speeds=np.r_[np.full(300,-8.),np.full(200,-24.)],
                 continuous_acceleration=np.linspace(-4.,-28.,500))
    rows, pair_rows = [],[]
    for name,truth in cases.items():
        position = np.r_[0,np.cumsum(truth)]
        images = np.array([np.interp(np.arange(1016)-s+1000,np.arange(len(source)),source) for s in position])
        dynamic = images-images.mean(axis=0,keepdims=True)
        lag_speeds = {}
        for lag in (1,2,4):
            corr = batch_ncc(dynamic,lag)
            lag_speeds[lag],peak = peak_shift(np.mean(corr,axis=0))
            lag_speeds[lag] /= lag
            if lag==1:
                individual=np.array([peak_shift(c)[0] for c in corr])
                for i,(actual,estimate) in enumerate(zip(truth,individual)):
                    pair_rows.append(dict(synthetic_case=name,first_frame=i,
                                          truth_shift_px=actual,estimated_shift_px=estimate))
        speed = lag_speeds[1]
        tol = max(2.,.2*abs(speed))
        blocks=[]
        for ids in np.array_split(np.arange(501),5):
            corr=batch_ncc(dynamic[ids],1)
            blocks.append(peak_shift(np.mean(corr,axis=0))[0])
        agreement=sum(abs(v-speed)<=tol for v in blocks)
        lag_consistent=all(abs(lag_speeds[k]-speed)<=tol for k in (2,4))
        rows.append(dict(synthetic_case=name,frames=501,true_time_mean_px_frame=truth.mean(),
            true_min_px_frame=truth.min(),true_max_px_frame=truth.max(),
            mean_individual_displacements_px_frame=individual.mean(),
            individual_mean_absolute_error_px_frame=np.mean(abs(individual-truth)),
            individual_max_absolute_error_px_frame=np.max(abs(individual-truth)),
            ensemble_peak_px_frame=speed,ensemble_peak_error_vs_time_mean=abs(speed-truth.mean()),
            lag2_ensemble_px_frame=lag_speeds[2],lag4_ensemble_px_frame=lag_speeds[4],
            temporal_tolerance_px_frame=tol,block_agree_count=agreement,
            block_consistency_gate=agreement>=4,lag_consistency_gate=lag_consistent,
            rejected_by_either_stationarity_gate=agreement<4 or not lag_consistent,
            block_speeds=';'.join(f'{x:.12g}' for x in blocks),
            status='synthetic_truth_known; not an experimental correction'))
    frame=save_csv('nonstationary_synthetic',rows)
    save_csv('nonstationary_synthetic_pairs',pair_rows)
    # The control must succeed before drawing the narrower methodological claim.
    assert frame.individual_mean_absolute_error_px_frame.max()<.1
    assert bool(frame.loc[frame.synthetic_case=='two_real_speeds','rejected_by_either_stationarity_gate'].iloc[0])
    return dict(cases=rows,
        conclusion='A single dominant correlation peak and stationarity gates can reject real varying optical motion; mean pair displacements recover the synthetic temporal mean. This does not establish phase identity or validate any rejected real-data peak.')


def physical_nonuniqueness():
    # Kinematic interface condition only. It supplies no momentum closure.
    x=np.linspace(0,1000,1001)
    c,h0,amplitude,wavelength=12.,20.,2.,200.
    derivative=amplitude*2*np.pi/wavelength*np.cos(2*np.pi*x/wavelength)
    kinematic=[]
    for name,ug,ul in [('kinematic_A',18.,6.),('kinematic_B',30.,4.)]:
        vg,vl=(ug-c)*derivative,(ul-c)*derivative
        # F=y-h(x-ct): F_t+u F_x+v F_y = c h_x-u h_x+v.
        residual_g=c*derivative-ug*derivative+vg
        residual_l=c*derivative-ul*derivative+vl
        kinematic.append(dict(case=name,identical_observed_interface_c_px_frame=c,
            assumed_axial_gas_px_frame=ug,assumed_axial_liquid_px_frame=ul,assumed_axial_ratio=ug/ul,
            max_normal_constraint_residual=max(abs(residual_g).max(),abs(residual_l).max()),
            interpretation='Two kinematically admissible interfacial motions; not full hydrodynamic solutions or measured phase-average speeds'))
    save_csv('interface_kinematic_nonuniqueness',kinematic)
    # Locally planar, no-slip laminar film, constant shear/pressure coefficients.
    # q(h)=tau*h^2/(2*mu)+K*h^3/(3*mu). c=dq/dh differs from q/h.
    h0,mu,c=1e-4,1e-3,.6
    films=[]
    for fraction in (0.,.25,.5,.75,1.):
        tau=fraction*mu*c/h0
        K=(1-fraction)*mu*c/h0**2
        q=tau*h0*h0/(2*mu)+K*h0**3/(3*mu)
        speed=tau*h0/mu+K*h0*h0/mu
        assert abs(speed-c)<1e-12
        films.append(dict(synthetic_shear_fraction=fraction,film_thickness_m=h0,viscosity_Pa_s=mu,
            interfacial_shear_Pa=tau,pressure_gravity_drive_Pa_per_m=K,
            identical_linear_wave_speed_m_s=speed,material_film_mean_m_s=q/h0,
            wave_to_material_factor=speed/(q/h0),
            status='Illustrative planar laminar model, not fitted to R600a recordings or a valid annular bound'))
    save_csv('film_wave_nonuniqueness',films)
    # Even a perfectly measured traveling volumetric alpha gives flux gradients
    # but leaves two boundary-flux constants. Mean total mass flux removes only one.
    x=np.linspace(0,2*np.pi,1001,endpoint=False)
    alpha=.5+.1*np.sin(x)
    alpha_x=.1*np.cos(x)
    speed,rhog,rhol,mean_G=2.,1.,10.,20.
    continuity=[]
    for Cg in (1.,3.):
        Cl=(mean_G-speed*(rhog*alpha.mean()+rhol*(1-alpha.mean()))-rhog*Cg)/rhol
        jg=speed*alpha+Cg
        jl=speed*(1-alpha)+Cl
        gas_residual=-speed*alpha_x+speed*alpha_x
        liquid_residual=speed*alpha_x-speed*alpha_x
        mean_mass=np.mean(rhog*jg+rhol*jl)
        assert abs(mean_mass-mean_G)<1e-12
        continuity.append(dict(synthetic_gas_flux_constant=Cg,synthetic_liquid_flux_constant=Cl,
            identical_alpha_min=alpha.min(),identical_alpha_max=alpha.max(),
            identical_void_structure_speed=speed,identical_total_mean_mass_flux=mean_mass,
            gas_flux_weighted_intrinsic_mean=jg.mean()/alpha.mean(),
            liquid_flux_weighted_intrinsic_mean=jl.mean()/(1-alpha.mean()),
            max_phase_continuity_residual=max(abs(gas_residual).max(),abs(liquid_residual).max()),
            status='Dimensionless synthetic 1D continuity solutions; momentum and constitutive closure not supplied'))
    save_csv('void_continuity_nonuniqueness',continuity)
    # Even granting the existing assumed alpha/density/G, one equation admits a line.
    data=read_csv(V151/'results/mass_balance_stress.csv')
    reference=data[(data.threshold==.7)&(data.massflow_statistic=='mean')&(data.axial_scale_multiplier==1)]
    assert len(reference)==26
    candidates=[]
    for r in reference.itertuples():
        for quality in (.1,.3,.5,.7,.9):
            ug=r.G_kg_m2_s*quality/(r.rho_g_kg_m3*r.native_circle_alpha)
            ul=r.G_kg_m2_s*(1-quality)/(r.rho_l_kg_m3*(1-r.native_circle_alpha))
            residual=r.rho_g_kg_m3*r.native_circle_alpha*ug+r.rho_l_kg_m3*(1-r.native_circle_alpha)*ul-r.G_kg_m2_s
            assert abs(residual)<1e-9
            candidates.append(dict(run=r.run,pattern=r.pattern,
                hypothetical_gas_mass_fraction=quality,assumed_alpha=r.native_circle_alpha,
                candidate_gas_m_s=ug,candidate_liquid_m_s=ul,ratio=ug/ul,
                conservation_residual_kg_m2_s=residual,
                gas_nonnegative_flux_upper_bound_m_s=r.G_kg_m2_s/(r.rho_g_kg_m3*r.native_circle_alpha),
                liquid_nonnegative_flux_upper_bound_m_s=r.G_kg_m2_s/(r.rho_l_kg_m3*(1-r.native_circle_alpha)),
                status='Algebraic family under unverified inputs; not experimental bounds or accepted estimates'))
    frame=save_csv('conservation_nonuniqueness',candidates)
    return dict(kinematic=kinematic,film_examples=films,void_continuity_examples=continuity,
                conservation_family_rows=len(frame),conservation_runs=frame.run.nunique(),
                max_mass_residual_kg_m2_s=abs(frame.conservation_residual_kg_m2_s).max(),
                interpretation='Conservation alone constrains a locus, not a unique point. Positive flux endpoints require extra assumptions and are not experimental confidence bounds.')


def real_gate_attrition():
    regions=read_csv(V151/'results/annular_region_summary.csv')
    rows=[]
    for r in regions.itertuples():
        stationarity_reasons=('lag_inconsistency','block_inconsistency')
        reasons=[] if pd.isna(r.withheld_reasons) else r.withheld_reasons.split(';')
        nonstationarity_only=bool(reasons) and all(x in stationarity_reasons for x in reasons)
        remaining_other_checks=bool(r.time_reversal_consistent and r.axial_template_consistent and r.all_lag_peaks_supported)
        rows.append(dict(run=r.run,pattern=r.pattern,region=r.region,
            robust_ensemble_supported=bool(r.robust_ensemble_supported),
            original_reasons=';'.join(reasons),baseline_rejection_only_stationarity=nonstationarity_only,
            other_extended_checks_pass=remaining_other_checks,
            would_only_relax_stationarity=bool(nonstationarity_only and remaining_other_checks),
            interpretation='Candidate for nonstationary follow-up; not automatically accepted'))
    frame=save_csv('legacy_gate_attrition',rows)
    return dict(region_rows=len(frame),
                baseline_rejected_only_stationarity=int(frame.baseline_rejection_only_stationarity.sum()),
                only_stationarity_with_other_extended_checks_pass=int(frame.would_only_relax_stationarity.sum()),
                warning='These counts flag investigable cases; loosening gates is not validation.')


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    inventory=read_csv(PHASE/'results/native_phase_inventory.csv')
    assert len(inventory)==26
    summary=dict(seed=SEED,scope='26 selected recordings; original exclusions preserved',
                 status='adversarial_controls_executed_not_phase_velocity_validation')
    summary['real_complement_masks']=complement_control(inventory)
    summary['weighted_common_motion']=weighted_common_motion()
    summary['synthetic_nonstationarity']=nonstationary_control()
    summary['physical_nonuniqueness']=physical_nonuniqueness()
    summary['legacy_gate_attrition']=real_gate_attrition()
    summary['sources_sha256']=SOURCES
    summary['script_sha256']=sha(Path(__file__))
    (OUT/'adversarial_identifiability.json').write_text(json.dumps(clean(summary),indent=2,allow_nan=False)+'\n')
    print(json.dumps(clean({k:v for k,v in summary.items() if k not in ('sources_sha256','physical_nonuniqueness')}),indent=2))


if __name__=='__main__':
    main()
