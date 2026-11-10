#!/usr/bin/env python3
"""Fixed-rule local-in-time optical propagation, not material velocimetry."""
from pathlib import Path
import sys, json, hashlib, time
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from PIL import Image

BASE = Path(__file__).resolve().parents[1]
OLD = BASE.parent
sys.path.insert(0, str(OLD/'scripts'))
from annular_motion import correlations, peak, mean_curve
REGIONS = ('core','upper_outer','lower_outer')
ROIS = {'core':(35,54),'upper_outer':(18,29),'lower_outer':(60,71)}
SEED = 20261024


def curve_result(p, lag=1, source=240):
    c,t = correlations(p[:-lag],p[lag:],source_x=source)
    valid=t>=1.
    r=peak(mean_curve(c,valid),np.arange(481)-source)
    r['coverage']=float(valid.mean())
    r['velocity']=r['shift']/lag
    return r


def local(p, seed, nulls=9):
    one=curve_result(p)
    tol=max(2.,.2*abs(one['velocity']))
    checks={'lag2':curve_result(p,2), 'reverse':curve_result(p[::-1]),
            'x160':curve_result(p,source=160),'x320':curve_result(p,source=320)}
    rng=np.random.default_rng(seed)
    null=[curve_result(p[rng.permutation(len(p))])['peak'] for _ in range(nulls)]
    maximum=float(np.max(null))
    reasons=[]
    if one['coverage']<.6: reasons.append('texture')
    if not np.isfinite(one['peak']) or one['peak']<.3: reasons.append('weak_peak')
    if one['boundary']: reasons.append('search_boundary')
    if not np.isfinite(maximum) or one['peak']<maximum+.05: reasons.append('shuffle')
    row={'velocity_px_frame':one['velocity'],'peak':one['peak'],'texture_coverage':one['coverage'],
         'null_max_peak':maximum,'tolerance_px_frame':tol,'boundary':one['boundary']}
    for name,r in checks.items():
        expected=-one['velocity'] if name=='reverse' else one['velocity']
        ok=bool(np.isfinite(r['velocity']) and abs(r['velocity']-expected)<=tol and
                not r['boundary'] and r['peak']>=.15 and r['coverage']>=.6)
        row[name+'_velocity_px_frame']=r['velocity']
        row[name+'_peak']=r['peak']
        row[name+'_boundary']=r['boundary']
        row[name+'_texture_coverage']=r['coverage']
        row[name+'_ok']=ok
        if not ok: reasons.append(name)
    row.update(supported=not reasons,reasons=';'.join(reasons))
    return row


def analyze(profile, nblocks, seed):
    result=[]
    for b,indices in enumerate(np.array_split(np.arange(len(profile)),nblocks)):
        r=local(profile[indices],seed+b)
        r.update(block=b,first_frame=int(indices[0]),last_frame=int(indices[-1]),
                 intervals=len(indices)-1)
        result.append(r)
    return result


def summarize(table):
    records=[]
    for (run,variant),g in table.groupby(['run','variant'],sort=False):
        d={region:g[g.region==region].sort_values('block').reset_index(drop=True) for region in REGIONS}
        n=len(d['core'])
        ok=np.logical_and.reduce([d[r].supported.to_numpy(bool) for r in REGIONS])
        weight=d['core'].intervals.to_numpy()
        midpoint=(d['core'].first_frame.to_numpy()+d['core'].last_frame.to_numpy())/2
        fifths=len(set(np.minimum((midpoint[ok]/501*5).astype(int),4)))
        vals=np.array([d[r].velocity_px_frame.to_numpy(float) for r in REGIONS])
        accepted=np.average(vals[:,ok],axis=1,weights=weight[ok]) if ok.any() else np.full(3,np.nan)
        all_values=np.average(vals,axis=1,weights=weight)
        fraction=float(ok.mean())
        supported=bool(fraction>=.6 and fifths>=4 and np.isfinite(accepted).all())
        fps=float(g.fps.iloc[0])
        row={'run':run,'pattern':g.pattern.iloc[0],'fps':fps,'variant':variant,'blocks':n,
             'joint_accepted_blocks':int(ok.sum()),'joint_block_fraction':fraction,'temporal_fifths':fifths,
             'joint_accepted_intervals':int(weight[ok].sum()),'total_internal_intervals':int(weight.sum()),
             'joint_interval_fraction_of_500':float(weight[ok].sum()/500),'supported':supported,
             'core_signed_px_s':accepted[0]*fps,'wall_signed_px_s':float(np.mean(accepted[1:]))*fps,
             'upper_signed_px_s':accepted[1]*fps,'lower_signed_px_s':accepted[2]*fps,
             'all_blocks_core_candidate_px_s':all_values[0]*fps,
             'all_blocks_wall_candidate_px_s':float(np.mean(all_values[1:]))*fps,
             'material_phase_pair':False}
        records.append(row)
    return pd.DataFrame(records)


def controls():
    rng=np.random.default_rng(SEED)
    texture=gaussian_filter1d(rng.normal(size=30000),2.)
    texture=texture/texture.std()*40
    x=np.arange(1016)
    rows=[]
    for speed in (-35.,-13.4,22.75):
        frames=np.stack([128+np.interp(x-speed*t+12000,np.arange(len(texture)),texture) for t in range(201)])
        for clip_fraction in (0.,.2,.5,.8):
            seq=frames.copy()
            if clip_fraction:
                cap=np.quantile(seq,1-clip_fraction)
                seq=np.minimum(seq,cap)
            seq-=seq.mean(axis=0)
            block=analyze(seq,5,SEED)
            est=np.mean([b['velocity_px_frame'] for b in block])
            error=abs(est-speed)
            rows.append({'control':'known_shift_with_clipping','truth_px_frame':speed,
                'upper_clipped_fraction':clip_fraction,'estimate_px_frame':est,'absolute_error_px_frame':error,
                'supported_blocks':sum(b['supported'] for b in block),'blocks':5,
                'within_declared_error':bool(error<=max(.5,.05*abs(speed)))})
    # Mean subtraction should remove a perfectly stationary scene.
    static=np.tile(rng.normal(size=1016),(51,1))
    static-=static.mean(axis=0)
    r=local(static,SEED)
    rows.append({'control':'static','supported':r['supported'],'reasons':r['reasons']})
    constant=np.ones((51,1016))
    r=local(constant,SEED)
    rows.append({'control':'constant','supported':r['supported'],'reasons':r['reasons']})
    return rows


def main():
    start=time.monotonic()
    scope=pd.read_csv(OLD/'results/scenario_inventory.csv')
    scope=scope[scope.included]
    archive=np.load(OLD/'intermediate/annular_profiles.npz')
    source_hash=hashlib.sha256((OLD/'intermediate/annular_profiles.npz').read_bytes()).hexdigest()
    rows=[];provenance=[]
    for index,entry in enumerate(scope.itertuples()):
        # Verify cached profiles against three original native images per run.
        files=sorted((OLD.parent/'data_processed'/entry.run).glob('ImgA*.jpg'))
        assert len(files)==501
        max_difference=0.
        for t in (0,250,500):
            frame=np.asarray(Image.open(files[t]).convert('L'))
            assert frame.shape==(80,1016)
            for region,(lo,hi) in ROIS.items():
                expected=frame[lo:hi].mean(axis=0)
                diff=np.max(abs(expected-archive[f'{entry.run}__{region}'][t]))
                max_difference=max(max_difference,float(diff))
                assert diff<2e-5
        provenance.append({'run':entry.run,'native_frames':501,'checked_frames':[0,250,500],
                           'max_float32_profile_difference':max_difference})
        for region in REGIONS:
            raw=archive[f'{entry.run}__{region}'].astype(float)
            p=raw-raw.mean(axis=0)
            alternatives=[('blocks10',10,p),('blocks5',5,p),('blocks20',20,p),
                          ('raw_blocks10',10,raw),
                          ('highpass8_blocks10',10,p-gaussian_filter1d(p,8.,axis=1)),
                          ('highpass32_blocks10',10,p-gaussian_filter1d(p,32.,axis=1))]
            for variant,nb,sequence in alternatives:
                for row in analyze(sequence,nb,SEED+index*100):
                    rows.append(dict(run=entry.run,pattern=entry.pattern,fps=entry.fps,region=region,variant=variant,**row))
        pd.DataFrame(rows).to_csv(BASE/'results/block_ncc_details.csv',index=False)
        print(entry.run,'done',round(time.monotonic()-start,1),'s',flush=True)
    data=pd.DataFrame(rows)
    summary=summarize(data)
    summary.to_csv(BASE/'results/block_ncc_per_run.csv',index=False)
    synth=controls()
    pd.DataFrame(synth).to_csv(BASE/'results/block_ncc_controls.csv',index=False)
    counts=summary.groupby('variant').supported.sum().astype(int).to_dict()
    metadata={'status':'complete','runtime_s':time.monotonic()-start,'counts_supported':counts,
        'source_profiles':str(OLD/'intermediate/annular_profiles.npz'),'source_profiles_sha256':source_hash,
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'sampled_native_profile_verification':provenance,'seed':SEED,'null_repeats_per_block':9,
        'controls':synth,'scope':26,'phase_material_pairs':0,
        'interpretation':'Interval-weighted mean optical propagation over jointly supported blocks; walls averaged equally; no phase-material identification',
        'note':'Block/null/support checks are empirical consistency gates, not confidence intervals or experimental accuracy.'}
    (BASE/'results/block_ncc_metadata.json').write_text(json.dumps(metadata,indent=2,allow_nan=False))
    print(counts,flush=True)


if __name__=='__main__':
    main()
