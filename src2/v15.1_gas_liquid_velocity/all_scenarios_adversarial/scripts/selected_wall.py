#!/usr/bin/env python3
"""Exploratory temporal selection/assessment of one observable wall."""
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd

BASE=Path(__file__).resolve().parents[1]


def main():
    source=BASE/'results/block_ncc_details.csv'
    details=pd.read_csv(source)
    rows=[]
    for variant in ('blocks10','raw_blocks10','highpass32_blocks10','highpass8_blocks10'):
        for run,g in details[details.variant==variant].groupby('run',sort=False):
            core=g[g.region=='core'].set_index('block')
            walls={r:g[g.region==r].set_index('block') for r in ('upper_outer','lower_outer')}
            for direction,train,test in [('forward',list(range(5)),list(range(5,10))),
                                         ('reverse',list(range(5,10)),list(range(5)))]:
                def rank(region):
                    w=walls[region].loc[train]
                    mean_peak=float(w.peak.mean())
                    return int(w.supported.sum()),mean_peak if np.isfinite(mean_peak) else -np.inf
                selected=max(walls,key=rank)
                wall=walls[selected]
                train_ok=core.loc[train].supported&wall.loc[train].supported
                test_ok=core.loc[test].supported&wall.loc[test].supported
                accepted=list(test_ok[test_ok].index)
                weights=core.loc[accepted].intervals
                fps=float(core.fps.iloc[0])
                def mean(frame,indices):
                    if not len(indices):return np.nan
                    return float(np.average(frame.loc[indices].velocity_px_frame,weights=core.loc[indices].intervals)*fps)
                other=walls['lower_outer' if selected=='upper_outer' else 'upper_outer']
                supported=bool(train_ok.sum()>=3 and test_ok.sum()>=3)
                rows.append({'run':run,'pattern':core.pattern.iloc[0],'variant':variant,'selection_direction':direction,
                    'fps':fps,'selected_wall':selected,'training_joint_blocks':int(train_ok.sum()),
                    'test_joint_blocks':int(test_ok.sum()),'test_blocks':5,
                    'test_intervals_accepted':int(weights.sum()),'test_intervals_total':int(core.loc[test].intervals.sum()),
                    'accepted_block_ids':','.join(map(str,accepted)),
                    'supported':supported,'core_signed_px_s':mean(core,accepted),'wall_signed_px_s':mean(wall,accepted),
                    'all_test_core_candidate_px_s':mean(core,test),'all_test_wall_candidate_px_s':mean(wall,test),
                    'other_wall_same_epochs_candidate_px_s':mean(other,accepted),
                    'material_phase_pair':False,
                    'observable':'Optical propagation in core and one training-selected wall, averaged only over jointly accepted assessment blocks; exploratory reused data'})
    output=pd.DataFrame(rows)
    output.to_csv(BASE/'results/selected_wall_per_run.csv',index=False)
    count=output.groupby(['variant','selection_direction']).supported.sum()
    metadata={'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'counts':{v+'__'+d:int(n) for (v,d),n in count.items()},
        'selection':'Training wall supported-block count, then mean NCC peak; no velocity magnitude or pattern used',
        'support':'At least 3/5 simultaneous core/wall blocks in training AND assessment halves',
        'limitation':'Exploratory analysis of reused recordings, not external validation. A selected wall is not the average film and optical motion is not gas/liquid material velocity.'}
    (BASE/'results/selected_wall_metadata.json').write_text(json.dumps(metadata,indent=2))
    print(count.to_string())


if __name__=='__main__':main()
