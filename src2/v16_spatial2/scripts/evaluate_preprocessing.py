"""Balance pattern groups and select preprocessing without the held-out run."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
MAIN='weighted_scale_peak_power05'
METRICS=[MAIN,'weighted_scale_peak_power01','weighted_scale_peak_power10',
         'weighted_zero_lag_wco','max_wco_all']

def summarize(data,metric):
    out=[]
    for variant,g in data.groupby('variant'):
        # Equal weight for each same-pattern class and each different-class pair.
        grouped=g.groupby(['same_pattern','pair_group'])[metric].mean()
        same=float(grouped.loc[True].mean()) if True in grouped.index.get_level_values(0) else float('nan')
        different=float(grouped.loc[False].mean()) if False in grouped.index.get_level_values(0) else float('nan')
        out.append(dict(variant=variant,same=same,different=different,gap=same-different,
                        valid_pairs=int(g[metric].notna().sum())))
    frame=pd.DataFrame(out).set_index('variant')
    base=frame.loc['baseline']
    for field in ['same','different','gap']:frame['delta_'+field]=frame[field]-base[field]
    return frame.sort_values('gap',ascending=False)

def choose(summary):
    eligible=summary[(summary.delta_same>0)&(summary.delta_gap>0)]
    return str(eligible.index[0]) if len(eligible) else 'baseline'

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--tag',default='screen');args=ap.parse_args()
    data=pd.read_csv(ROOT/'results'/f'preprocessing_{args.tag}_pairs.csv')
    expected=len(data.variant.unique())*325
    assert len(data)==expected
    assert data.groupby('variant').same_pattern.sum().eq(87).all()
    rankings=[];choices={}
    for metric in METRICS:
        s=summarize(data,metric);choices[metric]=choose(s)
        rankings.append(s.reset_index().assign(metric=metric))
    ranks=pd.concat(rankings,ignore_index=True)
    ranks.to_csv(ROOT/'results'/f'preprocessing_{args.tag}_rankings.csv',index=False)
    held=[]
    runs=sorted(set(data.run_x)|set(data.run_y))
    for run in runs:
        is_test=(data.run_x==run)|(data.run_y==run)
        train=data[~is_test];test=data[is_test]
        choice=choose(summarize(train,MAIN))
        score=summarize(test,MAIN)
        b=score.loc['baseline'];c=score.loc[choice]
        pattern=str(test.loc[test.run_x==run,'pattern_x'].iloc[0]) if (test.run_x==run).any() else str(test.loc[test.run_y==run,'pattern_y'].iloc[0])
        held.append(dict(run=run,pattern=pattern,selected=choice,train_pairs=len(train[train.variant=='baseline']),
                         test_pairs=len(test[test.variant=='baseline']),same_base=b.same,same_selected=c.same,
                         different_base=b.different,different_selected=c.different,delta_same=c.same-b.same,
                         delta_different=c.different-b.different,delta_gap=(c.same-c.different)-(b.same-b.different)))
    held=pd.DataFrame(held)
    held.to_csv(ROOT/'results'/f'preprocessing_{args.tag}_heldout.csv',index=False)
    balanced=held[held.delta_same.notna()].groupby('pattern')[['delta_same','delta_different','delta_gap']].mean().mean()
    out=dict(main_metric=MAIN,choice_on_full_dataset=choices[MAIN],choices_by_metric=choices,
        selection_rule='Improve same-pattern mean AND improve same-minus-different gap; choose highest training gap; otherwise baseline',
        validation='Leave one acquisition out of preprocessing selection; all 25 pairs involving it are withheld; no classifier',
        heldout_pattern_balanced_deltas={k:float(v) for k,v in balanced.items()},
        heldout_selection_counts={k:int(v) for k,v in held.selected.value_counts().items()},
        heldout_runs_with_same_pattern_comparator=int(held.delta_same.notna().sum()),
        caveat='Exploratory reuse of 26 acquisitions. Held-out pairs share training partners and folds overlap; no independent test-set or p-value claim.')
    (ROOT/'results'/f'preprocessing_{args.tag}_evaluation.json').write_text(json.dumps(out,indent=2)+'\n')
    print(ranks[ranks.metric==MAIN].to_string(index=False));print(json.dumps(out,indent=2))

if __name__=='__main__':main()
