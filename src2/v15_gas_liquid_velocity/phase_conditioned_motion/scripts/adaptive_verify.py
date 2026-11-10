#!/usr/bin/env python3
"""Independently check saved adaptive phase assignments, gates and summaries."""
from pathlib import Path
from collections import defaultdict
import csv
import json
import math
import numpy as np
from scipy.ndimage import minimum_filter
BASE=Path(__file__).resolve().parents[1]
OUT=BASE/'results'


def read(name):
    with (OUT/name).open() as f:return list(csv.DictReader(f))


def true(v):
    assert v in ('True','False');return v=='True'


def run():
    tracks=read('adaptive_tracks.csv');opps=read('adaptive_opportunities.csv')
    summary=read('adaptive_phase_summary.csv');ratios=read('adaptive_ratios.csv')
    inv={r['run']:r for r in read('native_phase_inventory.csv')}
    groups=defaultdict(list)
    for r in tracks:groups[(r['run'],int(r['width']),r['phase'],int(r['frame_index']))].append(r)
    numeric=0
    def same(a,b):
        nonlocal numeric
        numeric+=1;a=float(a);b=float(b)
        assert np.isclose(a,b,equal_nan=True,rtol=1e-11,atol=1e-10),(a,b)
    for run in sorted(inv):
        current=[r for r in tracks if r['run']==run]
        tids=sorted({int(r['frame_index']) for r in current}|{int(r['next_frame_index']) for r in current})
        index={t:i for i,t in enumerate(tids)}
        with np.load(BASE/'intermediate'/'masks'/(run+'.npz'),allow_pickle=False) as z:
            masks=[np.unpackbits(z[k][tids],axis=-1,count=1016,bitorder='big').astype(bool) for k in ('masks_060','masks_070','masks_080')]
        lumen=np.zeros((80,1016),bool);lumen[18:71]=True
        gas=np.logical_and.reduce(masks);liquid=(~np.logical_or.reduce(masks))&lumen
        interiors={p:minimum_filter(m.astype(np.uint8),size=(1,7,7),mode='constant',cval=0).astype(bool)
                   for p,m in [('gas_candidate',gas),('liquid_candidate',liquid)]}
        for r in current:
            t=int(r['frame_index']);u=int(r['next_frame_index']);x=int(r['source_x']);y=int(r['source_y'])
            bx=int(r['matched_x']);by=int(r['matched_y']);w=int(r['width']);h=int(r['height'])
            masks=interiors[r['phase']]
            assert masks[index[t],y:y+h,x:x+w].all()
            inside=bx>=0 and bx+w<=1016 and by>=18 and by+h<=71
            phase=bool(inside and masks[index[u],by:by+h,bx:bx+w].all())
            assert true(r['destination_phase_pass'])==phase
            forward=float(r['peak'])>=.6 and float(r['margin'])>=.03 and float(r['texture'])>=1 and not true(r['boundary'])
            assert forward==true(r['forward_pass'])
            rb=bx-64>=0 and bx+w+64<=1016 and by-4>=0 and by+h+4<=80
            assert rb==true(r['reverse_search_in_bounds'])
            if forward and rb:
                rt=math.hypot((bx-x)+float(r['reverse_dx']),(by-y)+float(r['reverse_dy']))
                same(r['roundtrip_px'],rt)
            motion=forward and rb and true(r['reverse_pass']) and float(r['roundtrip_px'])<=1
            assert motion==true(r['motion_pass'])
            accepted=motion and phase and float(r['target_saturation'])<=.20
            assert accepted==true(r['accepted'])
            same(r['velocity_px_s'],float(r['dx_px'])*float(inv[run]['fps']))
            assert float(r['source_std'])>=1 and float(r['source_saturation'])<=.2
    for o in opps:
        g=groups[(o['run'],int(o['width']),o['phase'],int(o['frame_index']))]
        assert len(g)==int(o['selected_patches'])<=2
        if len(g)==2:
            a,b=g;x,y=int(a['source_x']),int(a['source_y']);xx,yy=int(b['source_x']),int(b['source_y'])
            w,h=int(a['width']),int(a['height'])
            assert not(x<xx+w and x+w>xx and y<yy+h and y+h>yy)
        same(o['accepted_patches'],sum(true(r['accepted']) for r in g))
        same(o['motion_failures'],sum(r['reason']=='motion_failure' for r in g))
        same(o['destination_phase_failures'],sum(r['reason']=='destination_phase_failure' for r in g))
        same(o['target_saturation_failures'],sum(r['reason']=='target_saturation' for r in g))
    derived={}
    for r in summary:
        good=[t for t in tracks if t['run']==r['run'] and t['width']==r['width'] and t['phase']==r['phase'] and true(t['accepted'])]
        byframe=defaultdict(list)
        for t in good:byframe[int(t['frame_index'])].append(float(t['dx_px']))
        v=[np.mean(x) for x in byframe.values()]
        m=float(np.mean(v)) if v else np.nan
        same(r['accepted_tracks'],len(good));same(r['accepted_frames'],len(byframe))
        same(r['mean_px_frame'],m);same(r['median_px_frame'],np.median(v) if v else np.nan)
        same(r['mean_px_s'],m*float(inv[r['run']]['fps']))
        same(r['patch_weighted_mean_px_frame'],np.mean([float(t['dx_px']) for t in good]) if good else np.nan)
        derived[(r['run'],r['width'],r['phase'])]=(m,len(good),len(byframe))
    for r in ratios:
        gm,gt,gf=derived[(r['run'],r['width'],'gas_candidate')]
        lm,lt,lf=derived[(r['run'],r['width'],'liquid_candidate')]
        finite=np.isfinite(gm) and np.isfinite(lm)
        okay=gt>=20 and lt>=20 and gf>=10 and lf>=10 and finite and abs(lm)>=1 and gm*lm>0 and float(inv[r['run']]['core_saturation_fraction'])<=.2
        assert okay==true(r['eligible_descriptive_ratio'])
        same(r['unfiltered_ratio'],gm/lm if finite and lm!=0 else np.nan)
        same(r['descriptive_ratio'],gm/lm if okay else np.nan)
    result=dict(passed=True,numeric_checks=numeric,tracks=len(tracks),opportunities=len(opps),
        phase_summaries=len(summary),ratios=len(ratios),phase_masks_checked_with_independent_minimum_filter=True,
        material_phase_identity_validated=False,source_selection_note='Saved source geometry, nonoverlap and texture gates checked; source texture ranking is specified by the analysis implementation.')
    (OUT/'adaptive_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':run()
