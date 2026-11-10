#!/usr/bin/env python3
"""Sensitivity only: never select a threshold to fill the plot."""
from pathlib import Path
import pandas as pd
import numpy as np
BASE=Path(__file__).resolve().parents[1]


def main():
    d=pd.read_csv(BASE/'results/split_background_details.csv')
    d=d[(d.variant=='training_background')&(d.selection_direction=='forward')]
    rows=[]
    for required in (1,2,3,4,5):
        for tolerance_fraction in (.1,.2,.3):
            q=d.copy()
            tolerance=np.maximum(2.,abs(q.velocity_px_frame)*tolerance_fraction)
            ok=(q.texture_coverage>=.6)&(q.peak>=.3)&(~q.boundary)&(q.peak>=q.null_max_peak+.05)
            for check in ('lag2','reverse','x160','x320'):
                expected=-q.velocity_px_frame if check=='reverse' else q.velocity_px_frame
                ok &= ((q[check+'_velocity_px_frame']-expected).abs()<=tolerance)&(q[check+'_peak']>=.15)&(~q[check+'_boundary'])&(q[check+'_texture_coverage']>=.6)
            q['local_supported']=ok
            for run,g in q.groupby('run',sort=False):
                c=g[g.region=='core'].set_index('block')
                ws={r:g[g.region==r].set_index('block') for r in ('upper_outer','lower_outer')}
                def rank(region):
                    w=ws[region].loc[:4]
                    mean=float(w.peak.mean())
                    return int(w.local_supported.sum()),mean if np.isfinite(mean) else -np.inf
                chosen=max(ws,key=rank);w=ws[chosen]
                a=int((c.loc[:4].local_supported&w.loc[:4].local_supported).sum())
                b=int((c.loc[5:].local_supported&w.loc[5:].local_supported).sum())
                rows.append(dict(run=run,required_blocks_per_half=required,tolerance_fraction=tolerance_fraction,
                    selected_wall=chosen,training_blocks=a,assessment_blocks=b,supported=a>=required and b>=required,
                    primary_setting=required==3 and tolerance_fraction==.2))
    out=pd.DataFrame(rows)
    out.to_csv(BASE/'results/reporting_gate_sensitivity.csv',index=False)
    counts=out.groupby(['required_blocks_per_half','tolerance_fraction']).supported.sum().reset_index(name='supported_pairs')
    counts.to_csv(BASE/'results/reporting_gate_sensitivity_counts.csv',index=False)
    print(counts.to_string(index=False))


if __name__=='__main__':main()
