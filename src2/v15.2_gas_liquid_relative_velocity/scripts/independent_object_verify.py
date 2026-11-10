"""Independent arithmetic/correspondence audit, not physical validation."""
from pathlib import Path
import ast, csv, gzip, hashlib, importlib.util, json, re, sys, time
import numpy as np
import pandas as pd
from PIL import Image
from scipy.ndimage import label, find_objects, distance_transform_edt

BASE=Path(__file__).resolve().parents[1]; SRC2=BASE.parent
sys.path.insert(0,str(SRC2/'v14_new_voidfraction/scripts'))
import panorama_void_alternatives as detector
THRESHOLDS=(.6,.7,.8)
COUNTS={}; ERRORS={}

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def bump(name,n=1):COUNTS[name]=COUNTS.get(name,0)+n
def check(a,b,name,tol=1e-8):
    if a is None or isinstance(a,float) and not np.isfinite(a):
        assert b is None or not np.isfinite(b),(name,a,b);return
    if b is None or isinstance(b,float) and not np.isfinite(b):
        assert not np.isfinite(a),(name,a,b);return
    if isinstance(a,(bool,np.bool_,str)):
        assert a==b,(name,a,b)
    else:
        error=float(abs(a-b));ERRORS[name]=max(ERRORS.get(name,0.),error)
        assert error<=tol,(name,a,b,error)
    bump('checked_scalar_values')

def square_spans(mask):
    values=[]
    for col in mask.T:
        changes=np.diff(np.r_[0,col.astype(int),0]);starts=np.flatnonzero(changes==1);stops=np.flatnonzero(changes==-1)
        values.append(np.sum((stops-starts)**2))
    return np.asarray(values,float)

def independent_objects(mask):
    labs,_=label(mask);out=[];host_ids=[];sizes=np.bincount(labs.ravel())
    for idx,box in enumerate(find_objects(labs),1):
        if box is None:continue
        yy,xx=box;w=xx.stop-xx.start;h=yy.stop-yy.start;area=int(sizes[idx])
        if area>1000 or w>60 or h>35:host_ids.append(idx)
        small=15<=area<=1000 and 3<=w<=60 and 3<=h<=35 and w/h<=4
        large=area>1000 and 60<=w<=900 and h>=20
        if not (small or large) or yy.start<20 or yy.stop>69 or xx.start<10 or xx.stop>1006:continue
        component=labs[box]==idx;y,x=np.where(component);s=square_spans(component)
        out.append(dict(x=float(np.mean(x)+xx.start),y=float(np.mean(y)+yy.start),area=area,
            xmin=xx.start,xmax=xx.stop,ymin=yy.start,ymax=yy.stop,width=w,height=h,
            object_class='small' if small else 'large',label_id=idx,host_distance_px=0.,
            volume_centroid_x=float(np.average(np.arange(xx.start,xx.stop),weights=s)),
            volume_proxy_px3=float(np.sum(s)*np.pi/4),
            central_alpha_conditional=float(np.median(s[w//4:w-w//4])/53**2)))
    host=np.isin(labs,host_ids)
    distance=distance_transform_edt(~host) if host.any() else np.full(mask.shape,np.inf)
    for o in out:
        if o['object_class']=='small':
            o['host_distance_px']=float(min(distance[labs==o['label_id']].min(),o['ymin']-18,71-o['ymax']))
    return out

def independent_envelopes(mask):
    active=mask[35:54].sum(axis=0)>=10
    diff=np.diff(np.r_[False,active,False].astype(int));starts=np.where(diff==1)[0];ends=np.where(diff==-1)[0];out=[]
    for left,right in zip(starts,ends):
        w=right-left
        if not (60<=w<=900 and left>=10 and right<=1006):continue
        m=mask[18:71,left:right];s=square_spans(m)
        out.append(dict(x=float((left+right-1)/2),y=44.,area=int(m.sum()),xmin=int(left),xmax=int(right),
            ymin=18,ymax=71,width=int(w),height=53,object_class='large',label_id=len(out),host_distance_px=0.,
            volume_centroid_x=float(np.average(np.arange(left,right),weights=s)),
            volume_proxy_px3=float(np.sum(s)*np.pi/4),central_alpha_conditional=float(np.median(s[w//4:w-w//4])/53**2),
            radial_contact_fraction=float(np.count_nonzero(m[0]|m[-1])/w)))
    return out

def independent_pairs(a,b,lag=1):
    if not a or not b:return set()
    costs=np.full((len(a),len(b)),np.inf)
    for i,aa in enumerate(a):
        for j,bb in enumerate(b):
            ratio=bb['area']/aa['area'];dy=bb['y']-aa['y'];distance=float(np.hypot(bb['x']-aa['x'],dy))
            if aa['object_class']==bb['object_class'] and distance<=60*lag and abs(dy)<=5*lag and .67**lag<=ratio<=1.5**lag:
                costs[i,j]=distance
    result=set()
    for i in range(len(a)):
        j=int(np.argmin(costs[i]));best=costs[i,j]
        if not np.isfinite(best) or int(np.argmin(costs[:,j]))!=i:continue
        alternative=list(costs[i,:j])+list(costs[i,j+1:])+list(costs[:i,j])+list(costs[i+1:,j])
        if not alternative or min(alternative)-best>=3:result.add((i,j))
    return result

def ols(nodes,key):
    t=np.asarray([n['frame'] for n in nodes],float);x=np.asarray([n[key] for n in nodes],float)
    design=np.column_stack([np.ones(len(t)),t-t.mean()]);coef=np.linalg.lstsq(design,x,rcond=None)[0]
    return float(coef[1]),float(np.linalg.norm(x-design@coef)/np.sqrt(len(t)))

def audit_tracks(sequence,tracks,method):
    links=[independent_pairs(a,b) for a,b in zip(sequence[:-1],sequence[1:])]
    reverse=[independent_pairs(b,a) for a,b in zip(sequence[:-1],sequence[1:])]
    lag2=[independent_pairs(a,b,2) for a,b in zip(sequence[:-2],sequence[2:])]
    expected={(f,i,j) for f,ps in enumerate(links) for i,j in ps};seen=set();seen_nodes=set();computed=[]
    for t in tracks:
        nodes=t['nodes'];frames=np.array([n['frame'] for n in nodes]);assert (np.diff(frames)==1).all()
        assert t['first_frame']==frames[0] and t['last_frame']==frames[-1] and t['intervals']==len(nodes)-1
        for node in nodes:
            key=(node['frame'],node['index']);assert key not in seen_nodes,(method,'duplicate node',key);seen_nodes.add(key)
            source=sequence[node['frame']][node['index']]
            for field,value in source.items():check(node[field],value,'node_'+field)
        for a,b in zip(nodes[:-1],nodes[1:]):
            link=(a['frame'],a['index'],b['index']);assert link not in seen;seen.add(link)
        v,rmse=ols(nodes,'x');vl,_=ols(nodes,'xmin');vr,_=ols(nodes,'xmax');vv,vrmse=ols(nodes,'volume_centroid_x')
        check(t['signed_px_frame'],v,'OLS_velocity');check(t['centroid_rmse_px'],rmse,'OLS_rmse')
        check(t['left_px_frame'],vl,'OLS_left');check(t['right_px_frame'],vr,'OLS_right')
        check(t['volume_centroid_px_frame'],vv,'OLS_volume');check(t['volume_centroid_rmse_px'],vrmse,'OLS_volume_rmse')
        dx=np.diff([n['x'] for n in nodes]);dy=np.diff([n['y'] for n in nodes]);acc=np.inf
        if len(dx)>1:acc=max(np.percentile(abs(np.diff(dx)),95),np.percentile(abs(np.diff(dy)),95))
        check(t['acceleration95_px'],acc,'acceleration95')
        lags=[(a['index'],b['index']) in lag2[a['frame']] for a,b in zip(nodes[:-2],nodes[2:])]
        reverse_ok=all((b['index'],a['index']) in reverse[a['frame']] for a,b in zip(nodes[:-1],nodes[1:]))
        lag_fraction=float(np.mean(lags)) if lags else 0.
        check(t['lag2_correspondence_fraction'],lag_fraction,'lag2_fraction');check(t['reverse_all'],reverse_ok,'reverse_all')
        reasons=[]
        if len(dx)<5:reasons.append('short_track')
        if acc>3:reasons.append('centroid_acceleration')
        if rmse>2:reasons.append('nonlinear_centroid')
        if abs(vl-vr)>max(2,.25*abs(v)):reasons.append('front_rear_deformation')
        if lag_fraction<.8:reasons.append('lag2_correspondence')
        if not reverse_ok:reasons.append('reverse_correspondence')
        original_reasons=t['rejection'].split(';') if t['rejection'] else []
        changed=set(reasons)^set(original_reasons)
        if changed:
            # Distinguish machine-rounding at an exact cutoff from an audit failure.
            assert changed=={'front_rear_deformation'} and abs(abs(vl-vr)-max(2,.25*abs(v)))<1e-9,(changed,t)
            bump('roundoff_boundary_rejection_differences')
            if bool(reasons)!=bool(original_reasons):bump('roundoff_boundary_acceptance_differences')
            reasons=original_reasons
        check(t['rejection'],';'.join(reasons),'rejection');check(t['accepted'],not reasons,'accepted')
        strict=bool(not reasons and abs(vl-vr)<=max(1,.1*abs(v)))
        if strict!=t['strict_deformation']:
            assert not reasons and abs(abs(vl-vr)-max(1,.1*abs(v)))<1e-9
            bump('roundoff_boundary_strict_differences');strict=t['strict_deformation']
        check(t['strict_deformation'],strict,'strict_deformation')
        isolated=t['object_class']=='small' and min(n['host_distance_px'] for n in nodes)>=5
        check(t['isolated_small'],isolated,'isolated_small')
        for field,key in [('area_max_min_ratio','area'),('volume_max_min_ratio','volume_proxy_px3')]:
            values=[n[key] for n in nodes];check(t[field],max(values)/min(values),field)
        check(t['volume_area_speed_difference_px_frame'],abs(vv-v),'volume_area_difference')
        check(t['central_alpha_conditional'],np.median([n['central_alpha_conditional'] for n in nodes]),'central_alpha')
        if method=='meniscus':
            check(v,(vl+vr)/2,'midpoint_is_mean_of_endpoint_velocities')
            check(t['radial_contact_fraction'],np.mean([n['radial_contact_fraction'] for n in nodes]),'radial_contact')
        computed.append(dict(t,accepted=not reasons,strict_deformation=strict,signed_px_frame=v))
    assert seen==expected,(method,'missing/extra matched links',len(seen),len(expected))
    bump(method+'_tracks',len(tracks));bump(method+'_links',len(seen));bump(method+'_unique_nodes',len(seen_nodes))
    return computed

def summary(tracks,group,fps):
    selected=[t for t in tracks if t['accepted'] and (t['object_class']==group or group=='isolated_small' and t['isolated_small'])]
    indices={f for t in selected for f in range(t['first_frame'],t['last_frame'])};blocks={f//100 for f in indices}
    v=np.array([t['signed_px_frame'] for t in selected]);strict=np.array([t['signed_px_frame'] for t in selected if t['strict_deformation']])
    weighted=np.repeat(v,[t['intervals'] for t in selected]);median=float(np.median(v)) if len(v) else np.nan
    return dict(accepted_tracks=len(selected),unique_intervals=len(indices),observed_blocks=len(blocks),median_px_frame=median,
        median_px_s=median*fps,q25_px_s=np.quantile(v,.25)*fps if len(v) else np.nan,q75_px_s=np.quantile(v,.75)*fps if len(v) else np.nan,
        interval_weighted_median_px_s=np.median(weighted)*fps if len(v) else np.nan,strict_tracks=len(strict),
        strict_median_px_s=np.median(strict)*fps if len(strict) else np.nan,supported=len(selected)>=3 and len(indices)>=50 and len(blocks)>=3)

def run():
    begun=time.monotonic();inv=pd.read_csv(SRC2/'v15.1_gas_liquid_velocity/results/scenario_inventory.csv');inv=inv[inv.included]
    assert len(inv)==26 and set(inv['frames'])=={501}
    result_tables={m:pd.read_csv(BASE/'results'/f'{prefix}_summary.csv') for m,prefix in [('object','object'),('meniscus','meniscus')]}
    track_tables={m:pd.read_csv(BASE/'results'/f'{prefix}_tracks.csv') for m,prefix in [('object','object'),('meniscus','meniscus')]}
    for m,t in result_tables.items():assert set(t.run)==set(inv.run) and set(t.threshold)==set(THRESHOLDS)
    refpath=SRC2/'v15_gas_liquid_velocity/phase_conditioned_motion/results/native_phase_reference.json'
    reference=np.asarray(json.loads(refpath.read_text())['radial_reference'])
    validation_files=['object_validation.json','meniscus_validation.json']
    for f in validation_files:json.loads((BASE/'results'/f).read_text())
    robust_rows=[];static_rows=[];cache_missing_objects_hash=[];manifest={};scope_rows=[]
    source=(BASE/'scripts/object_motion.py').read_text();lines=source.splitlines(keepends=True);tree=ast.parse(source)
    sources={node.name:''.join(lines[node.lineno-1:node.end_lineno]) for node in tree.body if isinstance(node,ast.FunctionDef)}
    expected_algorithm_hash=hashlib.sha256((sources['objects']+sources['weights']).encode()).hexdigest()
    for record in inv.to_dict('records'):
        name=record['run'];fps=record['fps'];folder=SRC2/'data_processed'/name
        xsv=(folder/'_config.xsv').read_text();mii=(folder/'_config.mii').read_text()
        native=float(re.search(r'^Rate=([\d.]+)$',xsv,re.M).group(1))
        second=float(re.search(r'^Film speed 1\s*:([\d.]+)$',mii,re.M).group(1));assert native==second==fps
        files=sorted(folder.glob('ImgA*.jpg'));assert len(files)==501
        images_sha=hashlib.sha256()
        for i,path in enumerate(files):
            assert path.name==f'ImgA{i:06d}.jpg';images_sha.update(path.name.encode());images_sha.update(hashlib.sha256(path.read_bytes()).digest())
        cache_path=BASE/'intermediate'/f'objects_{name}.json.gz'
        with gzip.open(cache_path,'rt') as f:cache=json.load(f)
        assert cache['images_sha256']==images_sha.hexdigest();assert cache['detector_sha256']==digest(Path(detector.__file__))
        assert cache['reference_sha256']==digest(refpath)
        if 'algorithm_sha256' not in cache:cache_missing_objects_hash.append(name)
        else:assert cache['algorithm_sha256']==expected_algorithm_hash
        manifest[str(cache_path.relative_to(BASE))]=digest(cache_path)
        for frame in (0,125,250,375,500):
            im=np.asarray(Image.open(files[frame]).convert('L'),float);assert im.shape==(80,1016)
            for th in THRESHOLDS:
                mask,_,_=detector.silhouette(im,reference,threshold=th,min_open_interval=1000000)
                actual=independent_objects(mask);stored=cache['sequences'][str(th)][frame]
                assert len(actual)==len(stored)
                for a,b in zip(actual,stored):
                    assert set(a)==set(b)
                    for key in a:check(a[key],b[key],'redetected_'+key)
                bump('native_frames_redetected_by_threshold');bump('native_objects_redetected',len(actual))
        maskpath=SRC2/'v15_gas_liquid_velocity/phase_conditioned_motion/intermediate/masks'/f'{name}.npz'
        all_summaries=[]
        with np.load(maskpath) as masks:
            for th,key in zip(THRESHOLDS,('masks_060','masks_070','masks_080')):
                for method in ('object','meniscus'):
                    if method=='object':sequence=cache['sequences'][str(th)];stem='tracks'
                    else:
                        unpacked=np.unpackbits(masks[key],axis=-1,count=1016,bitorder='big').astype(bool)
                        sequence=[independent_envelopes(m) for m in unpacked];stem='meniscus'
                        bump('meniscus_masks_reprocessed',len(sequence))
                    path=BASE/'intermediate'/f'{stem}_{name}_{th:.1f}.json.gz'
                    with gzip.open(path,'rt') as f:tracks=json.load(f)
                    manifest[str(path.relative_to(BASE))]=digest(path)
                    recalculated=audit_tracks(sequence,tracks,method)
                    tcsv=track_tables[method];tcsv=tcsv[(tcsv.run==name)&(tcsv.threshold==th)].set_index('track_id')
                    assert len(tcsv)==len(tracks)
                    for t in recalculated:
                        row=tcsv.loc[t['track_id']]
                        check(row.signed_px_s,t['signed_px_frame']*fps,'native_fps_conversion',tol=1e-6)
                        check(row.accepted,t['accepted'],'csv_acceptance')
                    groups=('small','isolated_small','large') if method=='object' else ('large',)
                    for group in groups:
                        values=summary(recalculated,group,fps)
                        table=result_tables[method];row=table[(table.run==name)&(table.threshold==th)&(table.group==group)]
                        assert len(row)==1;row=row.iloc[0]
                        for field,value in values.items():check(row[field],value,'summary_'+field,tol=1e-6)
                        all_summaries.append(dict(method=method,group=group,threshold=th,**values))
                        if method=='object' and group in ('small','isolated_small'):
                            selected=[t for t in recalculated if t['accepted'] and (t['object_class']==group or group=='isolated_small' and t['isolated_small'])]
                            moving=[t for t in selected if abs(t['signed_px_frame'])>2]
                            static_rows.append(dict(run=name,pattern=record['pattern'],threshold=th,group=group,
                                all_tracks=len(selected),near_zero_tracks=sum(abs(t['signed_px_frame'])<=2 for t in selected),
                                baseline_median_px_s=values['median_px_s'],
                                motion_above_2_median_px_s=float(np.median([t['signed_px_frame'] for t in moving])*fps) if moving else np.nan,
                                motion_above_2_tracks=len(moving),
                                status='Diagnostic post-hoc sensitivity, not physical validation or revised baseline'))
        frame=pd.DataFrame(all_summaries)
        for (method,group),grouped in frame.groupby(['method','group']):
            base=grouped[grouped.threshold==.7].iloc[0];v=float(base.median_px_frame);tol=max(1.5,.2*abs(v))
            ok=bool(grouped.supported.all() and np.isfinite(v) and np.sign(grouped.median_px_frame).eq(np.sign(v)).all()
                    and ((grouped.median_px_frame-v).abs()<=tol).all())
            robust_rows.append(dict(run=name,pattern=record['pattern'],method=method,group=group,
                baseline_supported=bool(base.supported),all_thresholds_supported=bool(grouped.supported.all()),
                threshold_robust=ok,reference_median_px_s=float(base.median_px_s),tolerance_px_frame=tol))
        scope_rows.append(dict(run=name,pattern=record['pattern'],fps=fps,native_frames=len(files),images_sha256=images_sha.hexdigest()))
        print(name,'verified',flush=True)
    for method,table in track_tables.items():assert len(table)==COUNTS[method+'_tracks']
    robust=pd.DataFrame(robust_rows);static=pd.DataFrame(static_rows)
    robust.to_csv(BASE/'results/object_independent_threshold_robustness.csv',index=False)
    static.to_csv(BASE/'results/object_stationary_sensitivity.csv',index=False)
    base=robust.groupby(['method','group']).agg(baseline_supported=('baseline_supported','sum'),threshold_robust=('threshold_robust','sum')).reset_index()
    accepted_counts={m:track_tables[m].query('threshold == 0.7 and accepted').groupby('object_class').size().to_dict() for m in track_tables}
    for path in [Path(__file__),BASE/'scripts/object_motion.py',BASE/'scripts/meniscus_motion.py',BASE/'findings/PROTOCOL.md']:
        manifest[str(path.relative_to(BASE))]=digest(path)
    for path in (BASE/'results').glob('object_*.csv'):manifest[str(path.relative_to(BASE))]=digest(path)
    for path in (BASE/'results').glob('meniscus_*.csv'):manifest[str(path.relative_to(BASE))]=digest(path)
    payload=dict(status='passed_numerical_audit_with_documented_interpretation_limitations',elapsed_s=time.monotonic()-begun,
        counts=COUNTS,max_errors=ERRORS,scope=scope_rows,support=base.to_dict('records'),
        baseline_accepted_trajectory_counts=accepted_counts,cache_without_objects_function_hash=cache_missing_objects_hash,
        algorithm_function_hash=expected_algorithm_hash,
        controls='All saved nodes, links, fits, gates, summary cells and FPS conversions re-evaluated independently.',
        physical_phase_velocity_validated=False,manifest=manifest)
    (BASE/'results/object_independent_verification.json').write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps(dict(counts=COUNTS,support=base.to_dict('records'),accepted_counts=accepted_counts),indent=2))

if __name__=='__main__':run()
