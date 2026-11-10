"""Robust endpoint regression; exploratory response to noisy meniscus localization.

Velocity here remains gas-envelope celerity. Fit residuals are diagnostics,
not physical accuracy estimates. No phase-liquid value is inferred.
"""
import gzip,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import theilslopes
import object_motion as obj
BASE=obj.BASE

def slope(t,x):
    v=float(theilslopes(x,t)[0]);inter=float(np.median(x-v*t));res=x-inter-v*t
    return v,float(np.median(abs(res))),float(np.quantile(abs(res),.9))

def assess(nodes):
    t=np.array([o['frame'] for o in nodes]);left=np.array([o['xmin'] for o in nodes],float);right=np.array([o['xmax'] for o in nodes],float)
    vl,el,ql=slope(t,left);vr,er,qr=slope(t,right);v=(vl+vr)/2;tol=max(1.5,.2*abs(v))
    half=[]
    for idx in np.array_split(np.arange(len(t)),2):
        if len(idx)<3:half.append(np.nan)
        else:half.append((slope(t[idx],left[idx])[0]+slope(t[idx],right[idx])[0])/2)
    reasons=[]
    if t[-1]-t[0]<20:reasons.append('fewer_than_20_intervals')
    if abs(v)<=2 or abs(v)*(t[-1]-t[0])<50:reasons.append('insufficient_displacement')
    if abs(vl-vr)>tol:reasons.append('front_rear_inconsistent')
    if max(el,er)>3 or max(ql,qr)>12:reasons.append('endpoint_fit_residual')
    if not all(np.isfinite(h) and abs(h-v)<=tol for h in half):reasons.append('half_velocity_inconsistent')
    return dict(robust_px_frame=v,robust_left_px_frame=vl,robust_right_px_frame=vr,
        endpoint_median_error_px=max(el,er),endpoint_q90_error_px=max(ql,qr),
        half1_px_frame=half[0],half2_px_frame=half[1],robust_accepted=not reasons,robust_rejection=';'.join(reasons))

def controls():
    rng=np.random.default_rng(152);rows=[]
    for case in ('rigid','jittered','growth','stationary','accelerated'):
        nodes=[]
        for i in range(61):
            shift=-8*i if case!='stationary' else 0
            if case=='accelerated':shift=-2*i-.15*i*i
            noise=rng.normal(0,2,2) if case=='jittered' else [0.,0.]
            nodes.append(dict(frame=i,xmin=800+shift+noise[0],xmax=900+shift+(3*i if case=='growth' else 0)+noise[1]))
        out=assess(nodes);rows.append(dict(case=case,**out))
        assert out['robust_accepted']==(case in ('rigid','jittered')),(case,out)
        if case in ('rigid','jittered'):assert abs(out['robust_px_frame']+8)<.1
    return rows

def run():
    tests=controls();tracks=pd.read_csv(BASE/'results/meniscus_tracks.csv');out=[];summary=[]
    for (name,th),group in tracks.groupby(['run','threshold'],sort=False):
        with gzip.open(BASE/'intermediate'/f'meniscus_{name}_{th:.1f}.json.gz','rt') as f:full=json.load(f)
        for row in group.to_dict('records'):
            nodes=full[int(row['track_id'])]['nodes'];r=assess(nodes);out.append(dict(**row,**r,robust_px_s=r['robust_px_frame']*row['fps']))
    df=pd.DataFrame(out)
    inv=pd.read_csv(obj.SRC2/'v15.1_gas_liquid_velocity/results/scenario_inventory.csv');inv=inv[inv.included]
    for r in inv.to_dict('records'):
        for th in obj.THRESHOLDS:
            sel=df[(df.run==r['run'])&(df.threshold==th)&df.robust_accepted]
            intervals=set(i for t in sel.to_dict('records') for i in range(t['first_frame'],t['last_frame']))
            blocks={min(i//100,4) for i in intervals}
            summary.append(dict(run=r['run'],pattern=r['pattern'],fps=r['fps'],threshold=th,accepted_tracks=len(sel),
                unique_intervals=len(intervals),observed_blocks=len(blocks),median_px_frame=sel.robust_px_frame.median(),
                median_px_s=sel.robust_px_s.median(),q25_px_s=sel.robust_px_s.quantile(.25),q75_px_s=sel.robust_px_s.quantile(.75),
                supported=len(sel)>=3 and len(intervals)>=50 and len(blocks)>=3))
    obj.write_csv('robust_meniscus_tracks.csv',out);obj.write_csv('robust_meniscus_summary.csv',summary)
    obj.write_csv('robust_meniscus_controls.csv',tests)
    (BASE/'results/robust_meniscus_metadata.json').write_text(json.dumps(dict(exploratory=True,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),interpretation='Gas-envelope endpoint motion; not liquid material speed.'),indent=2)+'\n')
    print(pd.DataFrame(summary).query('threshold==0.7')[['run','accepted_tracks','median_px_s','supported']].to_string(index=False))

if __name__=='__main__':run()
