"""Descriptive object comparisons and a falsifiable zero-size extrapolation probe."""
import gzip,json
import numpy as np
import pandas as pd
import object_motion as obj
BASE=obj.BASE

def eligible(r):return bool(r['accepted'] and abs(r['signed_px_frame'])>2 and abs(r['signed_px_frame'])*r['intervals']>=20)

def stability(rows):
    a=rows.set_index('threshold');v=a.loc[.7,'median_px_frame'];tol=max(1.5,.2*abs(v))
    return bool(a.supported.all() and np.isfinite(a.median_px_frame).all() and
        (np.sign(a.median_px_frame)==np.sign(v)).all() and (abs(a.median_px_frame-v)<=tol).all())

def run():
    tracks=pd.read_csv(BASE/'results/object_tracks.csv');camera=pd.read_csv(obj.SRC2/'v15.1_gas_liquid_velocity/results/calibration_inventory.csv').set_index('run')
    inv=pd.read_csv(obj.SRC2/'v15.1_gas_liquid_velocity/results/scenario_inventory.csv');inv=inv[inv.included]
    summaries=[];resolved=[]
    for r in inv.to_dict('records'):
        for th in obj.THRESHOLDS:
            with gzip.open(BASE/'intermediate'/f"tracks_{r['run']}_{th:.1f}.json.gz",'rt') as f:ts=json.load(f)
            for t in ts:
                t['accepted']=eligible(t)
                if th==.7 and t['accepted'] and t['object_class']=='small':
                    nodes=t['nodes'];y=float(np.median([n['y'] for n in nodes]));h=t['median_height_px']
                    resolved.append(dict(run=r['run'],pattern=r['pattern'],track_id=t['track_id'],first_frame=t['first_frame'],last_frame=t['last_frame'],
                        signed_px_s=t['signed_px_frame']*r['fps'],height_px=h,width_px=t['median_width_px'],y_px=y,
                        block=min(int((t['first_frame']+t['last_frame'])/2)//100,4),isolated_small=t['isolated_small'],
                        blur_travel_over_height=abs(t['signed_px_frame'])*r['fps']*camera.loc[r['run'],'exposure_ns']*1e-9/h))
            for group in ('small','isolated_small'):
                summaries.append(dict(run=r['run'],pattern=r['pattern'],fps=r['fps'],threshold=th,**obj.summarize(ts,group,r['fps'])))
    ss=pd.DataFrame(summaries);mm=pd.read_csv(BASE/'results/robust_meniscus_summary.csv');output=[]
    for r in inv.to_dict('records'):
        name=r['run'];s=ss[(ss.run==name)&(ss.group=='small')];i=ss[(ss.run==name)&(ss.group=='isolated_small')];m=mm[mm.run==name]
        a=s[s.threshold==.7].iloc[0];b=m[m.threshold==.7].iloc[0];st=stability(s);mt=stability(m)
        output.append(dict(run=name,pattern=r['pattern'],small_supported=bool(a.supported),small_threshold_stable=st,
            small_px_s=a.median_px_s if a.supported else np.nan,small_tracks=int(a.accepted_tracks),small_intervals=int(a.unique_intervals),
            isolated_small_stable=stability(i),meniscus_supported=bool(b.supported),meniscus_threshold_stable=mt,
            meniscus_px_s=b.median_px_s if b.supported else np.nan,meniscus_tracks=int(b.accepted_tracks),meniscus_intervals=int(b.unique_intervals),
            small_over_meniscus_ratio=a.median_px_s/b.median_px_s if st and mt else np.nan,
            ratio_interpretation='different object populations and selected intervals; not phase slip',material_phase_pair=False))
    obj.write_csv('motion_resolved_object_summary.csv',summaries);obj.write_csv('motion_resolved_small_tracks.csv',resolved);obj.write_csv('object_relative_per_run.csv',output)
    # Competing low-order extrapolations. The intercept at size zero is NOT liquid velocity.
    df=pd.DataFrame(resolved);fits=[]
    for name,group in df.groupby('run'):
        if len(group)<12 or group.block.nunique()<3:continue
        y=group.signed_px_s.to_numpy();h=group.height_px.to_numpy()/20;rad=(group.y_px.to_numpy()-44.5)/26.5
        for model in ('constant','height','height_squared','height_squared_radial'):
            columns=[np.ones(len(group))]
            if model=='height':columns.append(h)
            if 'squared' in model:columns.append(h*h)
            if 'radial' in model:columns.append(rad*rad)
            x=np.column_stack(columns);coef=np.linalg.lstsq(x,y,rcond=None)[0];pred=np.full(len(y),np.nan)
            intercepts=[]
            for block in sorted(group.block.unique()):
                test=group.block.to_numpy()==block;train=~test
                # Purge training tracks sharing any native interval with held-out tracks.
                test_rows=group[test]
                for t in test_rows.to_dict('records'):
                    overlap=(group.first_frame.to_numpy()<t['last_frame'])&(group.last_frame.to_numpy()>t['first_frame'])
                    train &= ~overlap
                if train.sum()<x.shape[1] or np.linalg.matrix_rank(x[train])<x.shape[1]:continue
                b=np.linalg.lstsq(x[train],y[train],rcond=None)[0];pred[test]=x[test]@b;intercepts.append(b[0])
            fits.append(dict(run=name,pattern=group.pattern.iloc[0],model=model,tracks=len(group),minimum_height_px=group.height_px.min(),maximum_height_px=group.height_px.max(),
                zero_size_intercept_px_s=coef[0],heldout_mae_px_s=float(np.nanmean(abs(y-pred))),
                heldout_predicted_tracks=int(np.isfinite(pred).sum()),
                leave_block_intercept_min_px_s=min(intercepts) if intercepts else np.nan,leave_block_intercept_max_px_s=max(intercepts) if intercepts else np.nan,
                design_condition_number=float(np.linalg.cond(x)),physical_liquid_identity_validated=False))
    obj.write_csv('zero_size_extrapolation.csv',fits)
    print(pd.DataFrame(output).to_string(index=False))

if __name__=='__main__':run()
