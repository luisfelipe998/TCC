#!/usr/bin/env python3
"""Test probe position and liquid-rich gate relaxation; retain all outcomes."""
from pathlib import Path
import csv
from collections import defaultdict
import importlib.util
import json
import numpy as np
BASE=Path(__file__).resolve().parents[1]
s=importlib.util.spec_from_file_location('conditioned',BASE/'scripts'/'conditioned_motion.py')
C=importlib.util.module_from_spec(s);s.loader.exec_module(C)
VARIANTS=[('center_strict',496,520,.05),('left_strict',208,232,.05),
          ('right_strict',624,648,.05),('center_liquidrich_020',496,520,.20),
          ('center_liquidrich_040',496,520,.40)]


def run():
    base={}
    for r in C.read_csv(C.V15/'results'/'optical_preferred_candidate_pairs.csv'):
        base[(r['run'],int(r['frame_index']))]=r
    temporal={}
    for r in C.read_csv(C.FRAME/'results'/'temporal_pairs.csv'):
        if int(r['lag_frames'])==2:temporal[(r['run'],int(r['frame_index']))]=r
    ref={(r['run'],r['motion_gate']):r for r in C.read_csv(C.OUT/'passage_ratios.csv') if r['config']=='consensus'}
    results=[];episodes=[]
    for inv in C.read_csv(C.OUT/'native_phase_inventory.csv'):
        name=inv['run'];sat=float(inv['core_saturation_fraction'])>.20
        with np.load(BASE/'intermediate'/'masks'/(name+'.npz'),allow_pickle=False) as data:
            for variant,lo,hi,liqmax in VARIANTS:
                occupancy=np.array([np.unpackbits(data[k][:,18:71,lo//8:hi//8],axis=-1,count=hi-lo,bitorder='big').mean(axis=(1,2)) for k in C.KEYS])
                states=np.full(501,-1,dtype=int);states[occupancy.min(axis=0)>=.60]=1
                states[occupancy.max(axis=0)<=liqmax]=0
                states,epids,eps=C.episodes(states)
                for gate in C.GATES:
                    values={p:[] for p in [0,1]};byep={p:defaultdict(list) for p in [0,1]}
                    for i in range(500):
                        phase=int(states[i]) if states[i]>=0 and states[i]==states[i+1] else -1
                        b=base[(name,i)];lg=temporal.get((name,i))
                        ok=phase>=0 and C.truth(b['heuristic_pass'])
                        if gate=='lag2_confirmed':
                            ok=bool(ok and lg and C.truth(lg['temporal_within_tolerance']) and i<499 and states[i+2]==phase and epids[i+2]==epids[i])
                        if ok:
                            v=float(b['displacement_px_per_frame']);values[phase].append(v);byep[phase][int(epids[i])].append(v)
                    stats={p:C.phase_stats(values[p],byep[p]) for p in [0,1]}
                    rr=C.ratio_row(stats[1],stats[0],sat)
                    if variant=='center_strict':
                        prior=ref[(name,gate)]
                        assert rr['eligible_descriptive_ratio']==C.truth(prior['eligible_descriptive_ratio'])
                        for k in ['gas_mean_px_frame','liquid_mean_px_frame','unfiltered_ratio_of_means']:
                            assert np.isclose(rr[k],float(prior[k]),equal_nan=True,rtol=1e-12,atol=1e-12),(name,gate,k)
                    results.append(dict(run=name,pattern=inv['pattern'],variant=variant,probe_start=lo,probe_stop=hi,
                        liquid_max_occupancy=liqmax,gas_min_occupancy=.60,motion_gate=gate,**rr,
                        gas_pairs=stats[1]['n_accepted_pairs'],liquid_pairs=stats[0]['n_accepted_pairs'],
                        gas_episodes=stats[1]['n_accepted_episodes'],liquid_episodes=stats[0]['n_accepted_episodes']))
                    for e in eps:
                        phase=1 if e['phase']=='gas_candidate' else 0
                        vals=byep[phase].get(e['episode_id'],[])
                        episodes.append(dict(run=name,variant=variant,motion_gate=gate,**e,
                            accepted_pairs=len(vals),mean_px_frame=C.mean(vals)))
    C.save_csv('passage_sensitivity_ratios.csv',results)
    C.save_csv('passage_sensitivity_episodes.csv',episodes)
    metadata=dict(variants=[dict(name=v,probe=[lo,hi],liquid_max_fraction=l) for v,lo,hi,l in VARIANTS],
        gas_min_fraction=.60,baseline_reproduced=True,ratio_rows=len(results),episode_rows=len(episodes),
        interpretation='Relaxed gates label liquid-rich projected appearance, not pure liquid or material velocity.')
    (C.OUT/'passage_sensitivity_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    for variant,*_ in VARIANTS:
        good=[r for r in results if r['variant']==variant and r['motion_gate']=='adjacent' and r['eligible_descriptive_ratio']]
        print(variant,[(r['run'],round(r['descriptive_ratio'],3)) for r in good],flush=True)


if __name__=='__main__':run()
