"""Compare small objects and large envelopes on exactly common native intervals."""
import gzip,json
import numpy as np
import pandas as pd
import object_motion as obj
from object_comparisons import eligible
BASE=obj.BASE

def run():
    robust=pd.read_csv(BASE/'results/robust_meniscus_tracks.csv')
    inv=pd.read_csv(obj.SRC2/'v15.1_gas_liquid_velocity/results/scenario_inventory.csv');inv=inv[inv.included]
    output=[];intervalrows=[]
    for r in inv.to_dict('records'):
        for th in obj.THRESHOLDS:
            with gzip.open(BASE/'intermediate'/f"tracks_{r['run']}_{th:.1f}.json.gz",'rt') as f:small=json.load(f)
            small=[t for t in small if t['object_class']=='small' and eligible(t)]
            large=robust[(robust.run==r['run'])&(robust.threshold==th)&robust.robust_accepted].to_dict('records')
            sv=[];gv=[];ids=[];sids=set();gids=set()
            for i in range(500):
                st=[t for t in small if t['first_frame']<=i<t['last_frame']]
                gt=[t for t in large if t['first_frame']<=i<t['last_frame']]
                if not st or not gt:continue
                a=float(np.mean([t['signed_px_frame'] for t in st])*r['fps'])
                b=float(np.mean([t['robust_px_frame'] for t in gt])*r['fps'])
                sv.append(a);gv.append(b);ids.append(i);sids.update(t['track_id'] for t in st);gids.update(t['track_id'] for t in gt)
                intervalrows.append(dict(run=r['run'],threshold=th,interval=i,small_px_s=a,meniscus_px_s=b,n_small=len(st),n_menisci=len(gt)))
            blocks={min(i//100,4) for i in ids};supported=len(ids)>=50 and len(blocks)>=3 and len(sids)>=3 and len(gids)>=3
            a=float(np.mean(sv)) if sv else np.nan;b=float(np.mean(gv)) if gv else np.nan
            ratio=a/b if supported and abs(b)>2*r['fps'] else np.nan
            output.append(dict(run=r['run'],pattern=r['pattern'],fps=r['fps'],threshold=th,common_intervals=len(ids),common_blocks=len(blocks),
                small_tracks=len(sids),meniscus_tracks=len(gids),supported=supported,small_common_px_s=a,meniscus_common_px_s=b,ratio=ratio))
    df=pd.DataFrame(output);consensus=[]
    for name,d in df.groupby('run',sort=False):
        a=d[d.threshold==.7].iloc[0];stable=bool(d.supported.all())
        for col in ('small_common_px_s','meniscus_common_px_s'):
            stable=stable and bool(np.isfinite(d[col]).all() and (np.sign(d[col])==np.sign(a[col])).all() and (abs(d[col]-a[col])<=max(1.5*a.fps,.2*abs(a[col]))).all())
        consensus.append(dict(run=name,pattern=a.pattern,baseline_supported=bool(a.supported),threshold_stable=stable,common_intervals=int(a.common_intervals),
            small_px_s=a.small_common_px_s if a.supported else np.nan,meniscus_px_s=a.meniscus_common_px_s if a.supported else np.nan,
            ratio=a.ratio if stable else np.nan,ratio_sensitivity_min=d.ratio.min() if stable else np.nan,ratio_sensitivity_max=d.ratio.max() if stable else np.nan,
            material_phase_pair=False))
    obj.write_csv('paired_object_variants.csv',output);obj.write_csv('paired_object_intervals.csv',intervalrows);obj.write_csv('paired_object_consensus.csv',consensus)
    # Inter-recording comparisons retain one common observable, with no inference about pattern effects.
    s=pd.read_csv(BASE/'results/motion_resolved_object_summary.csv');per=pd.read_csv(BASE/'results/object_relative_per_run.csv')
    names=per.loc[per.small_threshold_stable,'run'].tolist();pairs=[]
    for k,i in enumerate(names):
        for j in names[:k]:
            a=s[(s.run==i)&(s.group=='small')];b=s[(s.run==j)&(s.group=='small')]
            av=abs(a[a.threshold==.7].median_px_s.iloc[0]);bv=abs(b[b.threshold==.7].median_px_s.iloc[0])
            pairs.append(dict(numerator_run=i,denominator_run=j,numerator_pattern=a.pattern.iloc[0],denominator_pattern=b.pattern.iloc[0],
                same_pattern=a.pattern.iloc[0]==b.pattern.iloc[0],image_speed_ratio=av/bv,
                segmentation_envelope_min=abs(a.median_px_s).min()/abs(b.median_px_s).max(),
                segmentation_envelope_max=abs(a.median_px_s).max()/abs(b.median_px_s).min(),
                physical_ratio_requires_equal_axial_magnification=True,phase_liquid_ratio=False))
    obj.write_csv('inter_recording_object_ratios.csv',pairs)
    print(pd.DataFrame(consensus).query('baseline_supported').to_string(index=False))

if __name__=='__main__':run()
