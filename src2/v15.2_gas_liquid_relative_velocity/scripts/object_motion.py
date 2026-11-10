"""Closed-object translation and shape change; no automatic liquid assignment."""
from pathlib import Path
import csv, gzip, hashlib, inspect, json, sys, time
import numpy as np
from PIL import Image
from scipy.ndimage import label, find_objects, distance_transform_edt
from scipy.spatial.distance import cdist
BASE=Path(__file__).resolve().parents[1]; SRC2=BASE.parent
sys.path.insert(0,str(SRC2/'v14_new_voidfraction/scripts'))
import panorama_void_alternatives as detector
THRESHOLDS=(.6,.7,.8)

def clean(v):
    if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [clean(x) for x in v]
    if isinstance(v,np.ndarray):return clean(v.tolist())
    if isinstance(v,np.bool_):return bool(v)
    if isinstance(v,np.integer):return int(v)
    if isinstance(v,(float,np.floating)):return float(v) if np.isfinite(v) else None
    return v

def write_csv(name,rows):
    with (BASE/'results'/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def weights(m):
    counts=np.zeros(m.shape[1]);s=np.zeros_like(counts)
    for row in m:
        s+=np.where(~row,counts**2,0);counts=np.where(row,counts+1,0)
    return s+counts**2

def objects(mask):
    labs,n=label(mask);boxes=find_objects(labs);sizes=np.bincount(labs.ravel())
    host=np.zeros_like(mask);found=[]
    for i,box in enumerate(boxes,1):
        if box is None:continue
        yy,xx=box;h=yy.stop-yy.start;w=xx.stop-xx.start;area=int(sizes[i])
        small=15<=area<=1000 and 3<=w<=60 and 3<=h<=35 and w/h<=4
        large=area>1000 and 60<=w<=900 and h>=20
        if area>1000 or w>60 or h>35:host[box]|=labs[box]==i
        if not(small or large):continue
        if yy.start<20 or yy.stop>69 or xx.start<10 or xx.stop>1006:continue
        m=labs[box]==i;y,x=np.nonzero(m);a=weights(m)
        found.append(dict(x=float(x.mean()+xx.start),y=float(y.mean()+yy.start),area=area,
            xmin=xx.start,xmax=xx.stop,ymin=yy.start,ymax=yy.stop,width=w,height=h,
            object_class='small' if small else 'large',label_id=i,host_distance_px=0.,
            volume_centroid_x=float(np.dot(np.arange(xx.start,xx.stop),a)/a.sum()),
            volume_proxy_px3=float(np.pi/4*a.sum()),
            central_alpha_conditional=float(np.median(a[w//4:w-w//4]/53**2))))
    distance=distance_transform_edt(~host) if host.any() else np.full(mask.shape,np.inf)
    for o in found:
        if o['object_class']=='small':
            box=(slice(o['ymin'],o['ymax']),slice(o['xmin'],o['xmax']))
            o['host_distance_px']=float(min(o['ymin']-18,71-o['ymax'],distance[box][labs[box]==o['label_id']].min()))
    return found

def matches(a,b,lag=1):
    if not a or not b:return []
    aa=np.array([[o['x'],o['y']] for o in a]);bb=np.array([[o['x'],o['y']] for o in b])
    dist=cdist(aa,bb);dy=abs(aa[:,None,1]-bb[None,:,1])
    ar=np.array([o['area'] for o in b])[None,:]/np.array([o['area'] for o in a])[:,None]
    same=np.array([o['object_class'] for o in a])[:,None]==np.array([o['object_class'] for o in b])[None,:]
    ok=(dist<=60*lag)&(dy<=5*lag)&(ar>=.67**lag)&(ar<=1.5**lag)&same
    cost=np.where(ok,dist,np.inf);fw=cost.argmin(1);rv=cost.argmin(0);out=[]
    for i,j in enumerate(fw):
        if rv[j]!=i or not np.isfinite(cost[i,j]):continue
        alt=np.concatenate([np.delete(cost[i],j),np.delete(cost[:,j],i)])
        margin=float(alt.min()-cost[i,j]) if len(alt) else np.inf
        if margin>=3:out.append((i,int(j)))
    return out

def fit(nodes,key):
    t=np.array([n['frame'] for n in nodes],float);t-=t.mean();x=np.array([n[key] for n in nodes])
    v=float(np.dot(t,x)/np.dot(t,t));return v,float(np.sqrt(np.mean((x-x.mean()-v*t)**2)))

def track(sequence):
    chains=[];incoming={}
    pairs=[matches(a,b) for a,b in zip(sequence[:-1],sequence[1:])]
    longer=[set(matches(a,b,2)) for a,b in zip(sequence[:-2],sequence[2:])]
    backward=[set(matches(b,a)) for a,b in zip(sequence[:-1],sequence[1:])]
    for t,ps in enumerate(pairs):
        for i,j in ps:
            if (t,i) not in incoming:
                tid=len(chains);chains.append([dict(frame=t,index=i,**sequence[t][i])])
            else:tid=incoming[t,i]
            incoming[t+1,j]=tid;chains[tid].append(dict(frame=t+1,index=j,**sequence[t+1][j]))
    output=[]
    for tid,nodes in enumerate(chains):
        v,rmse=fit(nodes,'x');vl,_=fit(nodes,'xmin');vr,_=fit(nodes,'xmax');vv,vrmse=fit(nodes,'volume_centroid_x')
        dx=np.diff([o['x'] for o in nodes]);dy=np.diff([o['y'] for o in nodes])
        acc=float(max(np.quantile(abs(np.diff(dx)),.95),np.quantile(abs(np.diff(dy)),.95))) if len(dx)>1 else np.inf
        reasons=[]
        if len(dx)<5:reasons.append('short_track')
        if acc>3:reasons.append('centroid_acceleration')
        if rmse>2:reasons.append('nonlinear_centroid')
        if abs(vl-vr)>max(2,.25*abs(v)):reasons.append('front_rear_deformation')
        lg=[(a['index'],b['index']) in longer[a['frame']] for a,b in zip(nodes[:-2],nodes[2:])]
        rev=[(b['index'],a['index']) in backward[a['frame']] for a,b in zip(nodes[:-1],nodes[1:])]
        lf=float(np.mean(lg)) if lg else 0.
        if lf<.8:reasons.append('lag2_correspondence')
        if not all(rev):reasons.append('reverse_correspondence')
        area=[o['area'] for o in nodes];volume=[o['volume_proxy_px3'] for o in nodes]
        output.append(dict(track_id=tid,object_class=nodes[0]['object_class'],first_frame=nodes[0]['frame'],
            last_frame=nodes[-1]['frame'],intervals=len(dx),signed_px_frame=v,
            volume_centroid_px_frame=vv,volume_centroid_rmse_px=vrmse,
            left_px_frame=vl,right_px_frame=vr,centroid_rmse_px=rmse,acceleration95_px=acc,
            lag2_correspondence_fraction=lf,reverse_all=all(rev),accepted=not reasons,rejection=';'.join(reasons),
            strict_deformation=bool(not reasons and abs(vl-vr)<=max(1,.1*abs(v))),
            isolated_small=bool(nodes[0]['object_class']=='small' and min(o['host_distance_px'] for o in nodes)>=5),
            host_distance_min_px=min(o['host_distance_px'] for o in nodes),
            median_height_px=float(np.median([o['height'] for o in nodes])),
            median_width_px=float(np.median([o['width'] for o in nodes])),
            area_max_min_ratio=max(area)/min(area),volume_max_min_ratio=max(volume)/min(volume),
            central_alpha_conditional=float(np.median([o['central_alpha_conditional'] for o in nodes])),
            volume_area_speed_difference_px_frame=abs(vv-v),nodes=nodes))
    return output

def summarize(tracks,group,fps):
    selected=[t for t in tracks if t['accepted'] and (t['object_class']==group or group=='isolated_small' and t['isolated_small'])]
    intervals=set(i for t in selected for i in range(t['first_frame'],t['last_frame']))
    blocks={min(i//100,4) for i in intervals};v=[t['signed_px_frame'] for t in selected]
    value=float(np.median(v)) if v else np.nan
    weighted=[t['signed_px_frame'] for t in selected for _ in range(t['intervals'])]
    strict=[t['signed_px_frame'] for t in selected if t['strict_deformation']]
    return dict(group=group,accepted_tracks=len(selected),unique_intervals=len(intervals),observed_blocks=len(blocks),
        median_px_frame=value,median_px_s=value*fps,
        q25_px_s=float(np.quantile(v,.25)*fps) if v else np.nan,q75_px_s=float(np.quantile(v,.75)*fps) if v else np.nan,
        interval_weighted_median_px_s=float(np.median(weighted)*fps) if v else np.nan,
        strict_tracks=len(strict),strict_median_px_s=float(np.median(strict)*fps) if strict else np.nan,
        supported=len(selected)>=3 and len(intervals)>=50 and len(blocks)>=3)

def synthetic():
    seq=[]
    for t in range(20):
        m=np.zeros((80,1016),bool);m[30:40,600-7*t:610-7*t]=1;seq.append(objects(m))
    ts=track(seq);assert len(ts)==1 and ts[0]['accepted']
    assert abs(ts[0]['signed_px_frame']+7)<1e-12 and abs(ts[0]['volume_centroid_px_frame']+7)<1e-12
    rev=track(seq[::-1]);assert abs(rev[0]['signed_px_frame']-7)<1e-12
    grow=[]
    for t in range(12):
        m=np.zeros((80,1016),bool);m[30:48,400-5*t:420-t]=1;grow.append(objects(m))
    assert not any(t['accepted'] for t in track(grow))
    a=[dict(x=100.,y=40.,area=100,object_class='small')]
    b=[dict(x=x,y=40.,area=100,object_class='small') for x in (96.,104.)]
    assert not matches(a,b)
    assert not track([objects(np.zeros((80,1016),bool)) for _ in range(12)])
    return dict(known_translation_error=abs(ts[0]['signed_px_frame']+7),reversal_error=abs(rev[0]['signed_px_frame']-7),
        growth_rejected=True,equal_competitors_rejected=True,empty_image_has_no_tracks=True,
        interpretation='Arithmetic and screening controls; not experimental validation.')

def run():
    start=time.monotonic();validation=synthetic()
    inv=[r for r in csv.DictReader((SRC2/'v15.1_gas_liquid_velocity/results/scenario_inventory.csv').open()) if r['included']=='True']
    refpath=SRC2/'v15_gas_liquid_velocity/phase_conditioned_motion/results/native_phase_reference.json'
    reference=np.array(json.loads(refpath.read_text())['radial_reference']);rows=[];flat=[];provenance=[]
    detector_hash=hashlib.sha256(Path(detector.__file__).read_bytes()).hexdigest()
    reference_hash=hashlib.sha256(refpath.read_bytes()).hexdigest()
    algorithm_hash=hashlib.sha256((inspect.getsource(objects)+inspect.getsource(weights)).encode()).hexdigest()
    for r in inv:
        begun=time.monotonic();name=r['run'];fps=float(r['fps']);path=BASE/'intermediate'/f'objects_{name}.json.gz'
        cache=None
        if path.exists():
            with gzip.open(path,'rt') as f:cache=json.load(f)
        files=sorted((SRC2/'data_processed'/name).glob('ImgA*.jpg'));assert len(files)==501
        check=hashlib.sha256()
        for i,f in enumerate(files):
            assert f.name==f'ImgA{i:06d}.jpg'
            check.update(f.name.encode());check.update(hashlib.sha256(f.read_bytes()).digest())
        current_digest=check.hexdigest()
        if (cache and cache.get('detector_sha256')==detector_hash and cache.get('reference_sha256')==reference_hash
                and cache.get('algorithm_sha256')==algorithm_hash and cache.get('images_sha256')==current_digest):
            sequences=cache['sequences'];digest=cache['images_sha256']
        else:
            files=sorted((SRC2/'data_processed'/name).glob('ImgA*.jpg'));assert len(files)==501
            sequences={str(th):[] for th in THRESHOLDS};sha=hashlib.sha256()
            for i,f in enumerate(files):
                assert f.name==f'ImgA{i:06d}.jpg';raw=f.read_bytes();sha.update(f.name.encode());sha.update(hashlib.sha256(raw).digest())
                im=np.array(Image.open(f).convert('L'),float);assert im.shape==(80,1016)
                for th in THRESHOLDS:
                    mask,_,_=detector.silhouette(im,reference,threshold=th,min_open_interval=1000000)
                    sequences[str(th)].append(objects(mask))
            digest=sha.hexdigest();assert digest==current_digest
            cache=dict(run=name,images_sha256=digest,detector_sha256=detector_hash,algorithm_sha256=algorithm_hash,
                reference_sha256=reference_hash,sequences=sequences)
            with gzip.open(path,'wt') as f:json.dump(clean(cache),f)
        for th in THRESHOLDS:
            tracks=track(sequences[str(th)])
            for group in ('small','isolated_small','large'):
                rows.append(dict(run=name,pattern=r['pattern'],fps=fps,threshold=th,
                    detected_objects=sum(len(os) for os in sequences[str(th)]),**summarize(tracks,group,fps)))
            for t in tracks:
                flat.append(dict(run=name,pattern=r['pattern'],fps=fps,threshold=th,
                    signed_px_s=t['signed_px_frame']*fps,**{k:v for k,v in t.items() if k!='nodes'}))
            with gzip.open(BASE/'intermediate'/f'tracks_{name}_{th:.1f}.json.gz','wt') as f:json.dump(clean(tracks),f)
        provenance.append(dict(run=name,frames=501,images_sha256=digest,detections_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        write_csv('object_summary.csv',rows);write_csv('object_tracks.csv',flat)
        print(name,round(time.monotonic()-begun,1),'seconds',[(x['group'],x['accepted_tracks'],x['supported']) for x in rows[-6:-3]],flush=True)
    write_csv('object_inputs.csv',provenance)
    (BASE/'results/object_validation.json').write_text(json.dumps(clean(dict(controls=validation,
        elapsed_s=time.monotonic()-start,recordings=len(inv),frames_per_run=501,thresholds=THRESHOLDS,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        detector_sha256=detector_hash,reference_sha256=reference_hash)),indent=2)+'\n')

if __name__=='__main__':run()
