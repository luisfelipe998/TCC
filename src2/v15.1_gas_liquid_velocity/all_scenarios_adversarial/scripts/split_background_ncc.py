#!/usr/bin/env python3
"""Remove assessment-half leakage from optical-background estimation."""
from pathlib import Path
import json,hashlib,time
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from block_ncc import analyze,REGIONS,SEED

BASE=Path(__file__).resolve().parents[1]


def main():
    start=time.monotonic()
    old=BASE.parent
    archive=np.load(old/'intermediate/annular_profiles.npz')
    scope=pd.read_csv(old/'results/scenario_inventory.csv')
    scope=scope[scope.included]
    details=[];pairs=[]
    for index,entry in enumerate(scope.itertuples()):
        for direction,train,test in [('forward',list(range(5)),list(range(5,10))),
                                     ('reverse',list(range(5,10)),list(range(5)))]:
            train_frames=np.arange(251) if direction=='forward' else np.arange(251,501)
            for filt,sigma in [('training_background',None),('training_background_hp32',32.),('training_background_hp8',8.)]:
                regions={}
                for region in REGIONS:
                    raw=archive[f'{entry.run}__{region}'].astype(float)
                    p=raw-raw[train_frames].mean(axis=0)
                    if sigma is not None:p=p-gaussian_filter1d(p,sigma,axis=1)
                    rows=analyze(p,10,SEED+index*100)
                    frame=pd.DataFrame(rows).set_index('block')
                    regions[region]=frame
                    for r in rows:details.append(dict(run=entry.run,pattern=entry.pattern,fps=entry.fps,
                        region=region,variant=filt,selection_direction=direction,**r))
                core=regions['core']
                def rank(region):
                    w=regions[region].loc[train]
                    mean_peak=float(w.peak.mean())
                    return int(w.supported.sum()),mean_peak if np.isfinite(mean_peak) else -np.inf
                selected=max(('upper_outer','lower_outer'),key=rank)
                wall=regions[selected]
                train_ok=core.loc[train].supported&wall.loc[train].supported
                test_ok=core.loc[test].supported&wall.loc[test].supported
                ids=list(test_ok[test_ok].index)
                def mean(frame,which):
                    return float(np.average(frame.loc[which].velocity_px_frame,weights=core.loc[which].intervals)*entry.fps) if len(which) else np.nan
                pairs.append(dict(run=entry.run,pattern=entry.pattern,fps=entry.fps,variant=filt,selection_direction=direction,
                    background_first_frame=int(train_frames[0]),background_last_frame=int(train_frames[-1]),
                    selected_wall=selected,training_joint_blocks=int(train_ok.sum()),test_joint_blocks=int(test_ok.sum()),
                    supported=bool(train_ok.sum()>=3 and test_ok.sum()>=3),accepted_block_ids=','.join(map(str,ids)),
                    test_intervals_accepted=int(core.loc[ids].intervals.sum()),
                    test_intervals_total=int(core.loc[test].intervals.sum()),
                    core_signed_px_s=mean(core,ids),wall_signed_px_s=mean(wall,ids),
                    all_test_core_candidate_px_s=mean(core,test),all_test_wall_candidate_px_s=mean(wall,test),
                    material_phase_pair=False))
        print(entry.run,round(time.monotonic()-start,1),'s',flush=True)
    data=pd.DataFrame(details);out=pd.DataFrame(pairs)
    data.to_csv(BASE/'results/split_background_details.csv',index=False)
    out.to_csv(BASE/'results/split_background_pairs.csv',index=False)
    counts=out.groupby(['variant','selection_direction']).supported.sum()
    meta={'runtime_s':time.monotonic()-start,'counts':{v+'__'+d:int(n) for (v,d),n in counts.items()},
        'source_profiles_sha256':hashlib.sha256((old/'intermediate/annular_profiles.npz').read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'background':'Each direction estimates the temporal mean only from its training half, frozen for all assessment frames.',
        'selection':'Wall ranked using training-only supported-block count then mean peak; speed magnitude and pattern are unused.',
        'limitations':'Internally separated processing/selection, but the experiment design explored these recordings; not new-data or independent physical validation. Coordinate means include only accepted assessment blocks.'}
    (BASE/'results/split_background_metadata.json').write_text(json.dumps(meta,indent=2))
    print(counts.to_string(),flush=True)


if __name__=='__main__':main()
