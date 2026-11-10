#!/usr/bin/env python3
"""Independent saved-result and mask-membership checks for adaptive patches."""
from pathlib import Path
import csv
import json
import math
from collections import defaultdict
import numpy as np
from scipy.ndimage import minimum_filter
BASE=Path(__file__).resolve().parents[1]
OUT=BASE/'results'


def read(name):
    with (OUT/name).open() as f:return list(csv.DictReader(f))


def yes(x):assert x in ('True','False');return x=='True'


def run():
    tracks=read('adaptive_tracks.csv');summary=read('adaptive_phase_summary.csv')
    ratios=read('adaptive_ratios.csv');opp=read('adaptive_opportunities.csv')
    rates={r['run']:float(r['fps']) for r in read('native_phase_inventory.csv')}
    sats={r['run']:float(r['core_saturation_fraction'])>.20 for r in read('native_phase_inventory.csv')}
    checks=0;maximum=0.;groups=defaultdict(list);byrun=defaultdict(list)
    def same(a,b):
        nonlocal checks,maximum
        a=float(a);b=float(b);checks+=1
        if math.isnan(a) or math.isnan(b):assert math.isnan(a) and math.isnan(b);return
        maximum=max(maximum,abs(a-b));assert math.isclose(a,b,abs_tol=1e-10,rel_tol=1e-11),(a,b)
    for r in tracks:
        groups[(r['run'],r['width'],r['phase'])].append(r);byrun[r['run']].append(r)
    for name,rs in byrun.items():
        ids=sorted({int(r['frame_index']) for r in rs}|{int(r['next_frame_index']) for r in rs})
        lookup={v:i for i,v in enumerate(ids)}
        with np.load(BASE/'intermediate'/'masks'/(name+'.npz'),allow_pickle=False) as z:
            ms=[np.unpackbits(z[k][ids],axis=-1,count=1016,bitorder='big').astype(bool) for k in ('masks_060','masks_070','masks_080')]
        lumen=np.zeros((80,1016),bool);lumen[18:71]=True
        candidates={'gas_candidate':np.logical_and.reduce(ms),'liquid_candidate':~np.logical_or.reduce(ms)&lumen}
        masks={k:minimum_filter(v.astype(np.uint8),size=(1,7,7),mode='constant',cval=0).astype(bool) for k,v in candidates.items()}
        for r in rs:
            i=int(r['frame_index']);j=int(r['next_frame_index']);assert j==i+1
            x=int(r['source_x']);y=int(r['source_y']);bx=int(r['matched_x']);by=int(r['matched_y'])
            w=int(r['width']);h=int(r['height']);m=masks[r['phase']]
            assert m[lookup[i],y:y+h,x:x+w].shape==(h,w) and m[lookup[i],y:y+h,x:x+w].all()
            dst=m[lookup[j],by:by+h,bx:bx+w] if 0<=by<=80-h and 0<=bx<=1016-w else np.zeros((0,0),bool)
            destination=dst.shape==(h,w) and dst.all()
            assert yes(r['destination_phase_pass'])==destination
            dx=float(r['dx_px']);dy=float(r['dy_px'])
            assert bx==x+round(dx) and by==y+round(dy)
            same(r['velocity_px_s'],dx*rates[name])
            forward=float(r['peak'])>=.6 and float(r['margin'])>=.03 and float(r['texture'])>=1 and not yes(r['boundary'])
            assert yes(r['forward_pass'])==forward
            if yes(r['reverse_pass']):
                assert yes(r['reverse_search_in_bounds']) and float(r['reverse_peak'])>=.6 and float(r['reverse_margin'])>=.03
            if forward and yes(r['reverse_search_in_bounds']):
                same(r['roundtrip_px'],math.hypot(round(dx)+float(r['reverse_dx']),round(dy)+float(r['reverse_dy'])))
            motion=forward and yes(r['reverse_pass']) and float(r['roundtrip_px'])<=1
            assert yes(r['motion_pass'])==motion
            accepted=motion and destination and float(r['source_saturation'])<=.2 and float(r['target_saturation'])<=.2
            assert yes(r['accepted'])==accepted
    summaries={}
    for r in summary:
        key=(r['run'],r['width'],r['phase']);rs=groups[key];good=[x for x in rs if yes(x['accepted'])]
        frame=defaultdict(list)
        for t in good:frame[int(t['frame_index'])].append(float(t['dx_px']))
        values=[float(np.mean(v)) for v in frame.values()]
        avg=float(np.mean(values)) if values else float('nan')
        same(r['selected_tracks'],len(rs));same(r['accepted_tracks'],len(good));same(r['accepted_frames'],len(frame))
        same(r['mean_px_frame'],avg);same(r['mean_px_s'],avg*rates[r['run']])
        same(r['median_px_frame'],np.median(values) if values else float('nan'))
        same(r['patch_weighted_mean_px_frame'],np.mean([float(t['dx_px']) for t in good]) if good else float('nan'))
        assert yes(r['high_saturation'])==sats[r['run']]
        summaries[key]=r
    for r in ratios:
        g=summaries[(r['run'],r['width'],'gas_candidate')];l=summaries[(r['run'],r['width'],'liquid_candidate')]
        gm=float(g['mean_px_frame']);lm=float(l['mean_px_frame'])
        finite=math.isfinite(gm) and math.isfinite(lm)
        expected=(int(g['accepted_tracks'])>=20 and int(l['accepted_tracks'])>=20 and
            int(g['accepted_frames'])>=10 and int(l['accepted_frames'])>=10 and finite and
            abs(lm)>=1 and gm*lm>0 and not sats[r['run']])
        assert yes(r['eligible_descriptive_ratio'])==expected
        same(r['gas_mean_px_frame'],gm);same(r['liquid_mean_px_frame'],lm)
        same(r['unfiltered_ratio'],gm/lm if finite and lm!=0 else float('nan'))
        same(r['descriptive_ratio'],gm/lm if expected else float('nan'))
    opgroups=defaultdict(list)
    for r in tracks:opgroups[(r['run'],r['width'],r['phase'],r['frame_index'])].append(r)
    for r in opp:
        rs=opgroups[(r['run'],r['width'],r['phase'],r['frame_index'])]
        same(r['selected_patches'],len(rs));same(r['accepted_patches'],sum(yes(x['accepted']) for x in rs))
        assert len(rs)<=2
        if len(rs)==2:
            a,b=rs;ax=int(a['source_x']);ay=int(a['source_y']);bx=int(b['source_x']);by=int(b['source_y'])
            w=int(a['width']);h=int(a['height']);assert ax+w<=bx or bx+w<=ax or ay+h<=by or by+h<=ay
    totals={}
    for width in ('64','32'):
        totals[width]={phase:sum(int(r['accepted_tracks']) for r in summary if r['width']==width and r['phase']==phase) for phase in ('gas_candidate','liquid_candidate')}
    result=dict(passed=True,tracks=len(tracks),opportunities=len(opp),summary_rows=len(summary),ratio_rows=len(ratios),
        numeric_checks=checks,maximum_difference=maximum,all_source_and_target_phase_memberships_verified=True,
        accepted_track_totals=totals,eligible_ratios=sum(yes(r['eligible_descriptive_ratio']) for r in ratios),
        scope='Reconstructs eroded mask membership and stored acceptance/summary arithmetic; does not validate phase identity or independently recompute all correlation surfaces.')
    (OUT/'adaptive_independent_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':run()
