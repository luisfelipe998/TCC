#!/usr/bin/env python3
"""Native-frame structure motion and temporal consistency, not phase velocity."""
from __future__ import annotations
import csv
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import re
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter1d, shift

HERE = Path(__file__).resolve().parents[1]
V15 = HERE.parent
DATA = V15.parent / 'data_processed'
for sub in ('results','figures','findings','intermediate'):
    (HERE/sub).mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR',str(HERE/'intermediate'/'mplcache'))
spec=importlib.util.spec_from_file_location('v15_ncc',V15/'scripts'/'optical_investigation.py')
NCC=importlib.util.module_from_spec(spec)
spec.loader.exec_module(NCC)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LAGS=(1,2,4)
# Fixed diagnostic gates; these do not describe measurement uncertainty.
ABS_TOL=1.0
REL_TOL=0.15
MIN_BASE_COVERAGE=0.50
MIN_JOINT_COVERAGE=0.25
MIN_AGREEMENT=0.80


def save_csv(name,rows):
    with (HERE/'results'/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def med(a):
    return float(np.median(a)) if len(a) else None


def quantile(a,q):
    return float(np.quantile(a,q)) if len(a) else None


def synthetic_check():
    rng=np.random.default_rng(20261023)
    base=gaussian_filter1d(rng.normal(size=1016),2)*30
    results=[]
    for speed in (-16.25,-6.3,0,7.4,20.1):
        seq=[shift(base,speed*k,order=3,mode='constant') for k in range(5)]
        for lag in LAGS:
            r=NCC.ncc(seq[0],seq[lag])
            err=abs(r['shift']/lag-speed)
            assert err<0.1,(speed,lag,err)
            results.append(dict(speed_px_frame=speed,lag=lag,error_px_frame=err))
    r=NCC.ncc(np.ones(1016),np.ones(1016))
    assert not NCC.quality(r,np.ones(1016),np.ones(1016),'ncc_dynamic')[0]
    return dict(trials=results,featureless_rejected=True,
                max_error_px_frame=max(r['error_px_frame'] for r in results))


def run():
    synthetic=synthetic_check()
    old={r['run']:r for r in csv.DictReader(
        (V15/'results'/'optical_preferred_candidate_per_run.csv').open())}
    pairs,inventory,summaries,excluded=[],[],[],[]
    examples={}
    for folder in sorted(DATA.glob('m??q??_???')):
        markers=list(folder.glob('__*.txt'))
        assert len(markers)==1,folder
        pattern=markers[0].stem[2:].replace('_',' ')
        if folder.name=='m05q05_002' or pattern not in NCC.PATTERNS:
            excluded.append(dict(run=folder.name,pattern=pattern))
            continue
        xsv=(folder/'_config.xsv').read_text()
        mii=(folder/'_config.mii').read_text()
        def field(name):
            return float(re.search(r'^'+re.escape(name)+r'=([^\r\n]+)',xsv,re.M)[1])
        fps=field('Rate')
        exposure_s=field('ExposureNs')*1e-9
        assert fps==float(re.search(r'Film speed 1\s*:([0-9.]+)',mii)[1])
        files=sorted(folder.glob('ImgA*.jpg'))
        ids=[int(re.search(r'(\d+)\.jpg$',p.name)[1]) for p in files]
        assert len(files)==field('Images')==501
        assert np.all(np.diff(ids)==1)
        profiles,saturation=[],[]
        digest=hashlib.sha256((folder/'_config.xsv').read_bytes())
        for file in files:
            raw=file.read_bytes()
            digest.update(file.name.encode())
            digest.update(hashlib.sha256(raw).digest())
            a=np.asarray(Image.open(io.BytesIO(raw)).convert('L'))
            assert a.shape==(80,1016)
            profiles.append(a[20:60].mean(axis=0))
            saturation.append(float(np.mean(a[20:60]>=250)))
        profiles=np.asarray(profiles)
        dynamic=profiles-profiles.mean(axis=0,keepdims=True)
        bylag={lag:NCC.quantify(dynamic,'ncc_dynamic',lag,range(501-lag),fps)
               for lag in LAGS}
        # Distinct central templates: reversal sensitivity is not a same-patch
        # round trip. Spatially different motion can cause disagreement.
        reverse=NCC.quantify(dynamic[::-1],'ncc_dynamic',1,range(500),fps)[::-1]
        unit=bylag[1]
        shifts=np.array([r['shift_px'] for r in unit])
        ok=np.array([r['heuristic_pass'] for r in unit])
        accepted=shifts[ok]
        high_sat=float(np.mean(saturation))>.20
        rsum=dict(run=folder.name,pattern=pattern,frames=501,
            first_frame_id=ids[0],last_frame_id=ids[-1],fps=fps,
            first_to_last_s=500/fps,exposure_us=exposure_s*1e6,
            core_saturation_fraction=float(np.mean(saturation)),high_saturation=high_sat,
            adjacent_pairs=500,adjacent_pass_pairs=int(ok.sum()),
            adjacent_pass_fraction=float(ok.mean()),median_signed_px_frame=med(accepted),
            q25_signed_px_frame=quantile(accepted,.25),q75_signed_px_frame=quantile(accepted,.75),
            median_signed_px_s=med(accepted)*fps,q25_signed_px_s=quantile(accepted,.25)*fps,
            q75_signed_px_s=quantile(accepted,.75)*fps,
            estimated_motion_during_exposure_px=abs(med(accepted))*fps*exposure_s,
            base_coverage_and_saturation_pass=bool(ok.mean()>=MIN_BASE_COVERAGE and not high_sat))
        assert int(old[folder.name]['heuristic_pass_pairs'])==int(ok.sum())
        assert abs(float(old[folder.name]['median_pass_px_per_frame'])-med(accepted))<1e-10
        reverse_errors=[]
        for lag,rows in bylag.items():
            diffs,relative,agree,longspeeds=[],[],[],[]
            for i,r in enumerate(rows):
                short_pass=bool(np.all(ok[i:i+lag]))
                short_speed=float(np.mean(shifts[i:i+lag]))
                speed=r['displacement_px_per_frame']
                delta=abs(speed-short_speed)
                rel=delta/max(abs(short_speed),1.)
                tol=max(ABS_TOL,REL_TOL*abs(short_speed))
                joint=bool(r['heuristic_pass'] and short_pass)
                rev=reverse[i] if lag==1 else None
                rev_joint=bool(joint and rev['heuristic_pass']) if rev else None
                rev_error=abs(r['shift_px']+rev['shift_px']) if rev else None
                if rev_joint: reverse_errors.append(rev_error)
                if joint:
                    diffs.append(delta)
                    relative.append(rel)
                    agree.append(delta<=tol)
                    longspeeds.append(speed)
                pairs.append(dict(run=folder.name,pattern=pattern,fps=fps,
                    first_frame_id=ids[i],second_frame_id=ids[i+lag],**r,
                    adjacent_path_pass=short_pass,adjacent_path_mean_px_frame=short_speed,
                    temporal_joint_pass=joint,temporal_abs_difference_px_frame=delta,
                    temporal_relative_difference=rel,temporal_tolerance_px_frame=tol,
                    temporal_within_tolerance=bool(joint and delta<=tol),
                    reversal_shift_px=rev['shift_px'] if rev else None,
                    reversal_pass=rev['heuristic_pass'] if rev else None,
                    reversal_joint_pass=rev_joint,reversal_abs_sum_px=rev_error))
            if lag>1:
                rsum.update({f'lag{lag}_pairs':len(rows),
                    f'lag{lag}_pass_pairs':sum(r['heuristic_pass'] for r in rows),
                    f'lag{lag}_joint_pairs':len(diffs),f'lag{lag}_joint_fraction':len(diffs)/len(rows),
                    f'lag{lag}_median_abs_difference_px_frame':med(diffs),
                    f'lag{lag}_median_relative_difference':med(relative),
                    f'lag{lag}_within_tolerance_fraction':float(np.mean(agree)) if agree else None,
                    f'lag{lag}_joint_median_px_frame':med(longspeeds)})
        rsum['reversal_joint_pairs']=len(reverse_errors)
        rsum['reversal_median_abs_sum_px']=med(reverse_errors)
        rsum['stable_under_lag_diagnostics']=bool(rsum['base_coverage_and_saturation_pass'] and all(
            rsum[f'lag{k}_joint_fraction']>=MIN_JOINT_COVERAGE and
            rsum[f'lag{k}_within_tolerance_fraction'] is not None and
            rsum[f'lag{k}_within_tolerance_fraction']>=MIN_AGREEMENT for k in (2,4)))
        summaries.append(rsum)
        inventory.append(dict(run=folder.name,pattern=pattern,fps=fps,frames=501,
            first_frame_id=ids[0],last_frame_id=ids[-1],shape=[80,1016],
            source_sha256=digest.hexdigest(),exposure_us=exposure_s*1e6))
        if folder.name in ('m07q01_001','m05q01_001','m05q05_001','m05q20_001','m06q80_001'):
            examples[folder.name]=(fps,shifts,ok,dynamic)
        print(folder.name,'accepted',int(ok.sum()),'/500','lag stable',rsum['stable_under_lag_diagnostics'],flush=True)
    assert len(summaries)==26 and len(excluded)==5
    save_csv('temporal_pairs.csv',pairs)
    save_csv('temporal_per_run.csv',summaries)
    metadata=dict(definition='Axial image-structure displacement; positive is rightward.',
        frame_lags=list(LAGS),source_roi_rows=[20,60],
        preprocessing='Subtract per-pixel mean across the complete 501-frame sequence.',
        ncc_source_sha256=hashlib.sha256((V15/'scripts'/'optical_investigation.py').read_bytes()).hexdigest(),
        gates=dict(ncc_peak=.60,peak_competitor_margin=.03,min_profile_std=1.,
            max_abs_search_px=240,excluded_boundary_px=2,base_min_coverage=MIN_BASE_COVERAGE,
            max_core_saturation=.20,temporal_abs_tolerance_px_frame=ABS_TOL,
            temporal_relative_tolerance=REL_TOL,
            temporal_tolerance_rule='max(absolute, relative * abs(adjacent-path mean))',
            lag_joint_min_fraction=MIN_JOINT_COVERAGE,lag_min_agreement_fraction=MIN_AGREEMENT),
        limitations=[
            'Temporal and fixed-window reversal discrepancies also reflect spatially varying motion and deformation.',
            'Diagnostic gates are heuristic, not measurement error bars.',
            'Adjacent comparisons overlap and are not independent experiment replicates.',
            'No independent phase-velocity truth or axial optical calibration is available.',
            'A 501-frame record spans only 0.0714 to 0.125 seconds.',
            'Exposure displacement is estimated from apparent motion, not independently measured blur.'
        ],inventory=inventory,exclusions=excluded,synthetic=synthetic,
        totals=dict(runs=26,adjacent_pairs=13000,adjacent_pass=sum(r['adjacent_pass_pairs'] for r in summaries),
            base_supported=sum(r['base_coverage_and_saturation_pass'] for r in summaries),
            lag_stable=sum(r['stable_under_lag_diagnostics'] for r in summaries),temporal_pairs=len(pairs)),
        old_preferred_reproduced=True)
    (HERE/'results'/'temporal_metadata.json').write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n')
    plot(summaries,examples)
    print(json.dumps(metadata['totals']),flush=True)


def plot(rows,examples):
    fig,axes=plt.subplots(1,3,figsize=(13,10),sharey=True,gridspec_kw={'width_ratios':[2.2,1,1.5]})
    for i,r in enumerate(rows):
        color=NCC.COLORS[r['pattern']]
        if r['base_coverage_and_saturation_pass']:
            x=r['median_signed_px_s']/1000
            axes[0].errorbar(x,i,xerr=[[x-r['q25_signed_px_s']/1000],[r['q75_signed_px_s']/1000-x]],
                fmt='o',color=color,ms=5,capsize=2)
        else:
            axes[0].text(.98,i,'unsupported',transform=axes[0].get_yaxis_transform(),ha='right',va='center',fontsize=8,color='0.45')
        axes[1].barh(i,r['adjacent_pass_fraction']*100,color=color,alpha=.8)
        for k,marker,dy in [(2,'o',-.12),(4,'s',.12)]:
            a=r[f'lag{k}_joint_fraction']
            if a is not None:
                axes[2].plot(a*100,i+dy,marker,color=color,ms=4,label=f'{k}-frame lag' if i==0 else None)
    axes[0].set_yticks(range(len(rows)),[r['run'] for r in rows],fontsize=8)
    axes[0].invert_yaxis()
    axes[0].set_xlabel('Signed structure motion (thousand px/s)\nBars: accepted-pair interquartile range')
    axes[1].set_xlabel('Accepted adjacent pairs (%)')
    axes[1].axvline(50,color='0.4',ls='--',lw=.8)
    axes[1].set_xlim(0,105)
    axes[2].set_xlabel('Joint comparison coverage (%)\nLonger lag and adjacent path both accepted')
    axes[2].set_xlim(0,105)
    axes[2].axvline(25,color='0.4',ls='--',lw=.8)
    axes[2].legend(fontsize=8,loc='lower left')
    for ax in axes: ax.grid(axis='x',alpha=.2)
    fig.suptitle('Native-frame structure motion: coverage across frame intervals\nPositive = rightward; no conversion to material phase velocity',fontsize=12)
    from matplotlib.lines import Line2D
    handles=[Line2D([0],[0],color=NCC.COLORS[p],marker='o',linestyle='',label=p) for p in NCC.PATTERNS]
    fig.legend(handles=handles,ncol=5,loc='upper center',bbox_to_anchor=(.55,.925),fontsize=8)
    fig.tight_layout(rect=(0,0,1,.90))
    for ext in ('png','svg'): fig.savefig(HERE/'figures'/f'temporal_motion_overview.{ext}',dpi=180)
    plt.close(fig)
    if not examples:
        return
    fig,axes=plt.subplots(len(examples),2,figsize=(12,10))
    for (run,(fps,d,ok,profiles)),axs in zip(examples.items(),axes):
        t=np.arange(len(d))/fps*1000
        axs[0].plot(t,d*fps/1000,color='0.75',lw=.7)
        axs[0].scatter(t[ok],d[ok]*fps/1000,s=3,color='#2457a6')
        axs[0].set_ylabel(run+'\n1000 px/s',fontsize=9)
        axs[0].grid(alpha=.2)
        lim=np.quantile(np.abs(profiles),.99)
        axs[1].imshow(profiles,aspect='auto',origin='upper',cmap='gray',vmin=-lim,vmax=lim,
            extent=[0,1016,500/fps*1000,0])
        axs[1].set_ylabel('Time (ms)')
    axes[-1,0].set_xlabel('Frame time (ms)')
    axes[-1,1].set_xlabel('Axial position (px)')
    axes[0,0].set_title('Adjacent estimates (blue = accepted)')
    axes[0,1].set_title('Space–time image profiles')
    fig.suptitle('Examples of structure correspondence and rejected motion estimates')
    fig.tight_layout(rect=(0,0,1,.97))
    fig.savefig(HERE/'figures'/'temporal_examples.png',dpi=180)
    plt.close(fig)


if __name__=='__main__':
    run()
