"""Axially closed core-envelope endpoints, including radially clipped masks.

This exploratory alternative addresses the closed-contour detector's failure
on wall-connected shadows. It measures interface translation, not liquid speed.
"""
import gzip,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.ndimage import label,find_objects
import object_motion as obj
BASE=obj.BASE;SRC2=obj.SRC2

def envelopes(mask):
    core=mask[35:54].mean(axis=0)>.5
    labs,n=label(core);out=[]
    for box in find_objects(labs):
        xx=box[0];w=xx.stop-xx.start
        if not (60<=w<=900 and xx.start>=10 and xx.stop<=1006):continue
        m=mask[18:71,xx];a=obj.weights(m);y,x=np.nonzero(m)
        out.append(dict(x=(xx.start+xx.stop-1)/2,y=44.,area=int(m.sum()),xmin=xx.start,xmax=xx.stop,
            ymin=18,ymax=71,width=w,height=53,object_class='large',label_id=len(out),host_distance_px=0.,
            volume_centroid_x=float(np.dot(np.arange(xx.start,xx.stop),a)/a.sum()),
            volume_proxy_px3=float(np.pi/4*a.sum()),central_alpha_conditional=float(np.median(a[w//4:w-w//4]/53**2)),
            radial_contact_fraction=float(np.mean(m[0]|m[-1]))))
    return out

def run():
    inv=pd.read_csv(SRC2/'v15.1_gas_liquid_velocity/results/scenario_inventory.csv');inv=inv[inv.included]
    rows=[];flat=[];inputs=[]
    for r in inv.to_dict('records'):
        name=r['run'];fps=r['fps'];path=SRC2/'v15_gas_liquid_velocity/phase_conditioned_motion/intermediate/masks'/f'{name}.npz'
        with np.load(path) as z:
            for th,key in zip(obj.THRESHOLDS,('masks_060','masks_070','masks_080')):
                m=np.unpackbits(z[key],axis=-1,count=1016,bitorder='big').astype(bool)
                sequence=[envelopes(mask) for mask in m];tracks=obj.track(sequence)
                for t in tracks:
                    t['radial_contact_fraction']=float(np.mean([n['radial_contact_fraction'] for n in t['nodes']]))
                    flat.append(dict(run=name,pattern=r['pattern'],fps=fps,threshold=th,signed_px_s=t['signed_px_frame']*fps,
                        **{k:v for k,v in t.items() if k!='nodes'}))
                rows.append(dict(run=name,pattern=r['pattern'],fps=fps,threshold=th,**obj.summarize(tracks,'large',fps)))
                with gzip.open(BASE/'intermediate'/f'meniscus_{name}_{th:.1f}.json.gz','wt') as f:json.dump(obj.clean(tracks),f)
        inputs.append(dict(run=name,mask_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    obj.write_csv('meniscus_summary.csv',rows);obj.write_csv('meniscus_tracks.csv',flat);obj.write_csv('meniscus_inputs.csv',inputs)
    # Independent endpoint arithmetic on a rigid translating synthetic core.
    seq=[]
    for t in range(20):
        m=np.zeros((80,1016),bool);m[24:66,600-9*t:800-9*t]=1;seq.append(envelopes(m))
    control=obj.track(seq)
    assert len(control)==1 and control[0]['accepted'] and abs(control[0]['signed_px_frame']+9)<1e-12
    (BASE/'results/meniscus_validation.json').write_text(json.dumps(dict(known_translation_error=abs(control[0]['signed_px_frame']+9),
        recordings=len(inv),frames_per_recording=501,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        interpretation='Interface endpoint translation; radial mask contact is retained and flagged.'),indent=2)+'\n')
    print(pd.DataFrame(rows).query('threshold==0.7')[['run','accepted_tracks','supported']].to_string(index=False))

if __name__=='__main__':run()
