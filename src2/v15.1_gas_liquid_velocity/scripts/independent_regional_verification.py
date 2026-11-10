#!/usr/bin/env python3
"""Audit regional saved evidence with direct Pearson calculations and gate reconstruction."""
from pathlib import Path
from collections import defaultdict
import hashlib
import json
import math
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / 'results'
COUNTS = defaultdict(int)
DIFF = defaultdict(float)
INPUTS = {}
ROIS = {'core':(35,54), 'upper_inner':(24,35), 'lower_inner':(54,65),
        'upper_outer':(18,29), 'lower_outer':(60,71)}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def remember(path):
    INPUTS[str(path.relative_to(BASE))] = sha(path)
    return path


def load(name):
    return pd.read_csv(remember(OUT/name))


def same(actual, expected, name, atol=1e-9, rtol=2e-12):
    a,b = np.asarray(actual,dtype=float),np.asarray(expected,dtype=float)
    assert a.shape == b.shape
    assert np.array_equal(np.isnan(a),np.isnan(b)), name
    assert np.allclose(a,b,atol=atol,rtol=rtol,equal_nan=True),(name,a,b)
    finite=np.isfinite(a)&np.isfinite(b)
    if finite.any(): DIFF[name]=max(DIFF[name],float(np.max(np.abs(a[finite]-b[finite]))))
    COUNTS[name] += a.size


def point_correlations(source, target, destination_x, source_x=240):
    a=source[:,source_x:source_x+536].astype(float)
    b=target[:,destination_x:destination_x+536].astype(float)
    a-=a.mean(axis=1)[:,None]
    b-=b.mean(axis=1)[:,None]
    numerator=np.einsum('ij,ij->i',a,b)
    denominator=np.sqrt(np.einsum('ij,ij->i',a,a)*np.einsum('ij,ij->i',b,b))
    return np.divide(numerator,denominator,out=np.full(len(a),np.nan),where=denominator>1e-9)


def peak_of_saved_curve(curve):
    if not np.isfinite(curve).any(): return math.nan,math.nan,math.nan,True
    k=int(np.nanargmax(curve));offset=0.
    if 0<k<480:
        left,center,right=map(float,curve[k-1:k+2])
        curvature=left-2*center+right
        if math.isfinite(curvature) and curvature < -1e-12:
            offset=max(-.5,min(.5,(left-right)/(2*curvature)))
    competitors=curve[np.abs(np.arange(481)-k)>5]
    return k-240+offset,float(curve[k]),float(curve[k]-np.nanmax(competitors)),bool(k<=2 or k>=478)


def main():
    regions=load('annular_region_summary.csv').set_index(['run','region'])
    lags=load('annular_lag_summary.csv').set_index(['run','region','lag'])
    blocks=load('annular_block_summary.csv').set_index(['run','region','block_index'])
    nulls=load('annular_null_trials.csv')
    additional=load('annular_additional_checks.csv').set_index(['run','region','check'])
    spatial=load('annular_spatial_pairs.csv').set_index('run')
    inventory=load('scenario_inventory.csv').set_index('run')
    metadata=json.loads(remember(OUT/'annular_metadata.json').read_text())
    assert regions.index.is_unique and lags.index.is_unique and blocks.index.is_unique
    assert len(regions)==148 and len(lags)==444 and len(blocks)==740 and len(spatial)==26
    assert len(nulls)==148*3*19 and not nulls.duplicated(['run','region','lag','shuffle_index']).any()
    assert set(spatial.index)==set(inventory[inventory.included].index)
    profile_path=remember(BASE/'intermediate/annular_profiles.npz')
    arrays_path=remember(OUT/'annular_pair_results.npz')
    direct_null_regions=0
    with np.load(profile_path,allow_pickle=False) as profiles, np.load(arrays_path,allow_pickle=False) as arrays:
        for n,((run,region),row) in enumerate(regions.iterrows()):
            prefix=run+'__'+region
            profile=profiles[prefix].astype(float)
            assert profile.shape==(501,1016)
            if region in ROIS:
                # Original mean profiles are integer row sums / band height.
                # This exactly reverses their lossless-at-integer-sum float32 archive.
                height=ROIS[region][1]-ROIS[region][0]
                profile=np.rint(profile*height)/height
            else:
                assert np.array_equal(profile*2,np.rint(profile*2))
            dynamic=profile-profile.mean(axis=0)[None,:]
            threshold=.25 if 'edge' in region else 1.
            checks_by_lag=[]
            for lag in (1,2,4):
                l=lags.loc[(run,region,lag)]
                pair=arrays[prefix+'__lag'+str(lag)].astype(float)
                curve=arrays[prefix+'__curve'+str(lag)].astype(float)
                assert pair.shape==(501-lag,5) and curve.shape==(481,)
                texture=np.minimum(dynamic[:-lag].std(axis=1),dynamic[lag:].std(axis=1))
                valid=texture>=threshold
                same(pair[:,3],texture,'pair_texture_float32',atol=1e-5,rtol=1e-7)
                k=np.rint(pair[:,0]).astype(int)+240
                accepted=valid&(pair[:,1]>=.6)&(pair[:,2]>=.03)&(k>2)&(k<478)
                assert np.array_equal(accepted,pair[:,4].astype(bool)),(run,region,lag)
                same([l.comparisons,l.valid_texture_pairs,l.individual_accepted_pairs,l.individual_coverage],
                     [501-lag,valid.sum(),accepted.sum(),accepted.mean()],'lag_counts')
                median=np.median(pair[accepted,0])/lag if accepted.any() else math.nan
                same(l.individual_median_px_frame,median,'lag_median_float32',atol=2e-5,rtol=1e-7)
                shift,p,m,boundary=peak_of_saved_curve(curve)
                same([l.ensemble_shift_px,l.ensemble_velocity_px_frame,l.ensemble_peak,l.ensemble_margin],
                     [shift,shift/lag,p,m],'ensemble_peak_float32',atol=.003,rtol=1e-7)
                assert bool(l.boundary)==boundary
                # Independent direct Pearson averages validate the peak neighborhood
                # and a separated target location across every original pair.
                best=int(np.nanargmax(curve))
                destinations=sorted(set([max(0,best-1),best,min(480,best+1),(37*n+53*lag)%481]))
                for destination in destinations:
                    corr=point_correlations(dynamic[:-lag],dynamic[lag:],destination)
                    value=float(np.nanmean(corr[valid])) if valid.any() else math.nan
                    same(curve[destination],value,'direct_ensemble_points',atol=1e-7,rtol=2e-7)
                # One original pair per region/lag: directly check its stored peak.
                sample=(17*n+41*lag)%(501-lag)
                point=int(k[sample])
                observed=point_correlations(dynamic[sample:sample+1],dynamic[sample+lag:sample+lag+1],point)[0]
                same(pair[sample,1],observed,'direct_individual_peak',atol=1e-7,rtol=2e-7)
                null=nulls[(nulls.run==run)&(nulls.region==region)&(nulls.lag==lag)].sort_values('shuffle_index')
                assert list(null.shuffle_index)==list(range(19))
                maximum=float(null.peak.max())
                same(l.shuffle_max_peak,maximum,'null_maxima')
                support=bool(np.isfinite(maximum) and l.ensemble_peak>=.10 and
                             l.ensemble_peak>=maximum+.03 and l.valid_texture_pairs>=250 and not l.boundary)
                assert bool(l.peak_and_shuffle_supported)==support
                checks_by_lag.append(support)
                # Recreate one shuffled trial for all core regions at every lag.
                if region=='core':
                    permutation=np.random.default_rng(20261023).permutation(501)
                    shuffled=dynamic[permutation]
                    destination=round(float(null.iloc[0].signed_px_frame)*lag)+240
                    corr=point_correlations(shuffled[:-lag],shuffled[lag:],destination)
                    valid_null=np.minimum(shuffled[:-lag].std(axis=1),shuffled[lag:].std(axis=1))>=threshold
                    value=float(np.nanmean(corr[valid_null]))
                    same(null.iloc[0].peak,value,'direct_null_peak',atol=1e-8,rtol=2e-10)
                    direct_null_regions+=1
            first=lags.loc[(run,region,1)]
            speed=float(first.ensemble_velocity_px_frame)
            tolerance=max(2.,.2*abs(speed))
            b=blocks.loc[(run,region)].sort_index()
            consistent=all(abs(float(lags.loc[(run,region,lag)].ensemble_velocity_px_frame)-speed)<=tolerance
                           and not bool(lags.loc[(run,region,lag)].boundary) for lag in (2,4))
            count=sum(math.isfinite(v) and abs(v-speed)<=tolerance for v in b.signed_px_frame)
            same([row.signed_velocity_px_frame,row.signed_velocity_px_s,row.temporal_tolerance_px_frame,
                  row.block_agree_count,row.min_block_px_frame,row.max_block_px_frame],
                 [speed,speed*float(inventory.loc[run].fps),tolerance,count,b.signed_px_frame.min(),b.signed_px_frame.max()],
                 'region_motion_summary')
            assert bool(row.lag_consistent)==consistent
            reasons=[]
            if first.valid_texture_pairs<250: reasons.append('low_texture_coverage')
            if not math.isfinite(speed) or abs(speed)<1: reasons.append('missing_or_near_zero_motion')
            if first.boundary: reasons.append('search_boundary')
            if not math.isfinite(first.ensemble_peak) or first.ensemble_peak<.10: reasons.append('weak_ensemble_peak')
            if not math.isfinite(first.shuffle_max_peak) or first.ensemble_peak<first.shuffle_max_peak+.03: reasons.append('weak_shuffle_separation')
            if not consistent: reasons.append('lag_inconsistency')
            if count<4: reasons.append('block_inconsistency')
            if row.saturation_fraction>.2: reasons.append('region_saturation')
            if 'edge' in region and (row.weak_edge_fraction>.2 or row.edge_search_boundary_fraction>.2): reasons.append('ambiguous_edge_location')
            assert bool(row.ensemble_supported)==(not reasons),(run,region,reasons)
            actual_reason=row.withheld_reasons if pd.notna(row.withheld_reasons) else ''
            assert set(filter(None,actual_reason.split(';')))==set(reasons)
            checks=additional.loc[(run,region)]
            for name,s in checks.iterrows():
                difference=abs(s.signed_px_frame+speed) if name=='time_reversed' else abs(s.signed_px_frame-speed)
                same(s.difference_px_frame,difference,'additional_differences')
                assert bool(s.within_tolerance)==bool(math.isfinite(difference) and difference<=tolerance and not s.boundary)
                assert s.min_shift_px==-s.source_x and s.max_shift_px==480-s.source_x and s.source_width==536
            reverse=bool(checks.loc['time_reversed'].within_tolerance)
            offsets=all(bool(checks.loc[name].within_tolerance) for name in ('source_x160','source_x320'))
            robust=not reasons and reverse and offsets and all(checks_by_lag)
            assert bool(row.time_reversal_consistent)==reverse
            assert bool(row.axial_template_consistent)==offsets
            assert bool(row.all_lag_peaks_supported)==all(checks_by_lag)
            assert bool(row.robust_ensemble_supported)==robust
            COUNTS['region_gate_checks']+=1
    for run,p in spatial.iterrows():
        regs=regions.loc[run]
        c,u,l=(regs.loc[name] for name in ('core','upper_inner','lower_inner'))
        wall=(float(u.signed_velocity_px_frame)+float(l.signed_velocity_px_frame))/2
        sign=np.sign(u.signed_velocity_px_frame)==np.sign(l.signed_velocity_px_frame)
        agree=bool(sign and abs(u.signed_velocity_px_frame-l.signed_velocity_px_frame)<=max(2.,.25*abs(wall)))
        basic=bool(c.ensemble_supported and u.ensemble_supported and l.ensemble_supported and agree and abs(wall)>=1)
        robust=bool(basic and c.robust_ensemble_supported and u.robust_ensemble_supported and l.robust_ensemble_supported)
        same([p.core_signed_px_frame,p.nearwall_signed_px_frame,p.upper_signed_px_frame,p.lower_signed_px_frame,
              p.core_signed_px_s,p.nearwall_signed_px_s],
             [c.signed_velocity_px_frame,wall,u.signed_velocity_px_frame,l.signed_velocity_px_frame,
              c.signed_velocity_px_s,wall*p.fps],'spatial_pair')
        assert bool(p.upper_lower_agreement)==agree
        assert bool(p.eligible_spatial_descriptor)==basic and bool(p.robust_spatial_descriptor)==robust
        if robust: assert p.core_signed_px_s<0 and p.nearwall_signed_px_s<0
    for path,digest in INPUTS.items():
        assert sha(BASE/path)==digest,('Input changed during verification',path)
    result=dict(status='passed',counts=dict(COUNTS),max_absolute_differences=dict(DIFF),
                regions=len(regions),baseline_supported_regions=int(regions.ensemble_supported.sum()),
                final_supported_regions=int(regions.robust_ensemble_supported.sum()),
                baseline_spatial_pairs=int(spatial.eligible_spatial_descriptor.sum()),
                final_spatial_pairs=int(spatial.robust_spatial_descriptor.sum()),
                direct_null_region_lags=direct_null_regions,input_sha256=INPUTS,verifier_sha256=sha(Path(__file__)),
                independence='No generator imported. Direct centered Pearson products verify ensemble peak neighborhoods for every region/lag and sampled individual/null peaks; gates and summaries independently reconstructed.',
                limitation='Validates saved arithmetic and specified support gates. It does not validate material identity, calibrated m/s, the entire correlation surface independently, or physical measurement accuracy.')
    (OUT/'independent_regional_verification.json').write_text(json.dumps(result,indent=2)+chr(10))
    print(json.dumps({k:v for k,v in result.items() if k!='input_sha256'},indent=2))


if __name__=='__main__': main()
