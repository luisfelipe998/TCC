#!/usr/bin/env python3
"""Exploratory isolated-object tracking after closed-envelope restriction.

Targets the one manually classified dispersed-bubble recording only. Component
filters and reciprocal geometric matches identify candidate bubble trajectories,
not phase ground truth or liquid velocity. Baseline inputs are unchanged.
"""
from pathlib import Path
import csv
import json
import numpy as np
from scipy.ndimage import label, find_objects
from scipy.spatial.distance import cdist

BASE=Path(__file__).resolve().parents[1]
OUT=BASE/'results'
RUN='m07q01_001'
FPS=7000.


def write_csv(name,rows):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def objects(mask):
    labels,n=label(mask)
    found=[]
    for component,box in enumerate(find_objects(labels),start=1):
        if box is None:continue
        yy,xx=box
        h,w=yy.stop-yy.start,xx.stop-xx.start
        area=int((labels[box]==component).sum())
        if not (15<=area<=1000 and 3<=w<=60 and 3<=h<=35 and w/h<=4):continue
        if yy.start<20 or yy.stop>69 or xx.start<10 or xx.stop>1006:continue
        y,x=np.nonzero(labels[box]==component)
        found.append(dict(x=float(x.mean()+xx.start),y=float(y.mean()+yy.start),area=area,
                          xmin=xx.start,xmax=xx.stop,ymin=yy.start,ymax=yy.stop))
    return found


def matches(a,b):
    if not a or not b:return []
    xy_a=np.array([[o['x'],o['y']] for o in a]);xy_b=np.array([[o['x'],o['y']] for o in b])
    distance=cdist(xy_a,xy_b)
    dy=np.abs(xy_a[:,None,1]-xy_b[None,:,1])
    areas=np.array([o['area'] for o in b])[None,:]/np.array([o['area'] for o in a])[:,None]
    allowed=(distance<=30)&(dy<=5)&(areas>=.5)&(areas<=2.)
    costs=np.where(allowed,distance,np.inf)
    forward=np.argmin(costs,axis=1);reverse=np.argmin(costs,axis=0)
    out=[]
    for i,j in enumerate(forward):
        if reverse[j]!=i or not np.isfinite(costs[i,j]):continue
        alt_a=np.delete(costs[i],j);alt_b=np.delete(costs[:,j],i)
        margin_a=float(alt_a.min()-costs[i,j]) if len(alt_a) else np.inf
        margin_b=float(alt_b.min()-costs[i,j]) if len(alt_b) else np.inf
        if margin_a<2 or margin_b<2:continue
        out.append((i,int(j),float(costs[i,j]),margin_a,margin_b))
    return out


def run():
    # A basic known translation verifies centroid/matching direction arithmetic.
    m=np.zeros((80,1016),bool);m[30:40,200:210]=1;m[35:44,500:509]=1
    n=np.zeros_like(m);n[:,:-7]=m[:,7:]
    a,b=objects(m),objects(n);pairs=matches(a,b)
    assert len(pairs)==2 and all(abs((b[j]['x']-a[i]['x'])+7)<1e-12 for i,j,*_ in pairs)
    rows=[];summary=[];object_rows=[]
    with np.load(OUT/'closed_contour_m07q01_masks.npz',allow_pickle=False) as z:
        for key in ('masks_060','masks_070','masks_080'):
            masks=np.unpackbits(z[key],axis=-1,count=1016,bitorder='big').astype(bool)
            sequence=[objects(mask) for mask in masks]
            chains={};next_id=0
            incoming={};links=[]
            for t in range(500):
                for i,j,d,ma,mb in matches(sequence[t],sequence[t+1]):
                    node=(t,i)
                    if node not in incoming:
                        chain=next_id;next_id+=1;chains[chain]=[]
                    else:chain=incoming[node]
                    incoming[(t+1,j)]=chain
                    a=sequence[t][i];b=sequence[t+1][j]
                    row=dict(run=RUN,threshold=key[-3:],frame_id=t,source_object=i,target_object=j,chain_id=chain,
                        x=a['x'],y=a['y'],next_x=b['x'],next_y=b['y'],source_area=a['area'],target_area=b['area'],
                        dx_px=b['x']-a['x'],dy_px=b['y']-a['y'],distance_px=d,
                        forward_margin_px=ma,reverse_margin_px=mb)
                    links.append(row);chains[chain].append(row)
            stable={}
            for chain,ls in chains.items():
                dx=np.array([r['dx_px'] for r in ls]);dy=np.array([r['dy_px'] for r in ls])
                # Five matched frame intervals and no large adjacent centroid jump.
                stable[chain]=bool(len(ls)>=5 and np.max(np.abs(np.diff(dx)))<=2 and np.max(np.abs(np.diff(dy)))<=2)
            accepted=[]
            for r in links:
                r['chain_intervals']=len(chains[r['chain_id']]);r['stable_chain']=stable[r['chain_id']]
                r['velocity_px_s']=r['dx_px']*FPS;rows.append(r)
                if r['stable_chain']:accepted.append(r['dx_px'])
            for t,os in enumerate(sequence):
                for i,o in enumerate(os):object_rows.append(dict(run=RUN,threshold=key[-3:],frame_id=t,object_id=i,**o))
            summary.append(dict(threshold=key[-3:],candidate_objects=sum(map(len,sequence)),
                matched_intervals=len(links),total_chains=len(chains),stable_chains=sum(stable.values()),
                stable_intervals=len(accepted),mean_px_frame=float(np.mean(accepted)) if accepted else np.nan,
                median_px_frame=float(np.median(accepted)) if accepted else np.nan,
                median_px_s=float(np.median(accepted))*FPS if accepted else np.nan,
                q25_px_frame=float(np.quantile(accepted,.25)) if accepted else np.nan,
                q75_px_frame=float(np.quantile(accepted,.75)) if accepted else np.nan))
    write_csv('closed_contour_bubble_objects.csv',object_rows)
    write_csv('closed_contour_bubble_links.csv',rows)
    write_csv('closed_contour_bubble_summary.csv',summary)
    metadata=dict(run=RUN,fps=FPS,area_px=[15,1000],width_px=[3,60],height_px=[3,35],maximum_aspect_ratio=4,
        valid_rows=[20,69],valid_columns=[10,1006],maximum_step_px=30,maximum_vertical_step_px=5,
        area_ratio=[.5,2.],reciprocal_nearest=True,minimum_competitor_distance_margin_px=2,
        minimum_chain_intervals=5,maximum_consecutive_dx_or_dy_change_px=2,
        synthetic_centroid_translation_passed=True,
        interpretation='Exploratory motion of selected small closed gas-silhouette candidates; not verified gas-phase mass velocity.',
        limitations=['Filters omit small, touching, coalescing, cropped and deformed objects.',
            'Entire-chain stability gate selects persistent easy objects.',
            'No independent bubble labels or velocity ground truth; comparison with global NCC is not accuracy validation.',
            'No liquid parcel motion is measured.'])
    (OUT/'closed_contour_bubble_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':run()
