"""Compara faixas com os mesmos pesos de referência 20–59, sem mudar WCO."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from evaluate_preprocessing import MAIN,METRICS,summarize

ROOT=Path(__file__).resolve().parents[1]
PREVIOUS=ROOT.parent/'v16_spatial2'

def read_band(root,band):
    data=pd.concat([pd.read_csv(root/'results/preprocessing_confirm_pairs.csv'),
                    pd.read_csv(root/'results/preprocessing_v5_pairs.csv')],ignore_index=True)
    out=data[['run_x','run_y','pattern_x','pattern_y','same_pattern','pair_group','variant','n']].copy()
    for metric in METRICS:
        field=('fixed20_59_'+metric) if band=='18_71' and metric!='max_wco_all' else metric
        out[metric]=data[field]
    out['band']=band
    return out

def main():
    old=read_band(PREVIOUS,'20_59');new=read_band(ROOT,'18_71')
    both=pd.concat([old,new],ignore_index=True)
    assert len(old)==len(new)==4*325
    both.to_csv(ROOT/'results/band_comparison_pairs.csv',index=False)
    summary=[]
    for band,data in both.groupby('band'):
        for metric in METRICS:
            rows=summarize(data,metric).reset_index()
            for row in rows.itertuples():
                g=data[data.variant==row.variant]
                summary.append(dict(band=band,variant=row.variant,metric=metric,
                    pooled_same=float(g[g.same_pattern][metric].mean()),
                    balanced_same=float(row.same),balanced_different=float(row.different),
                    balanced_gap=float(row.gap)))
    pd.DataFrame(summary).to_csv(ROOT/'results/band_comparison_summary.csv',index=False)
    changes=[];all_changes=[]
    for variant in sorted(new.variant.unique()):
        x=new[new.variant==variant].merge(old[old.variant==variant],on=['run_x','run_y'],
                                          suffixes=('_new','_old'),validate='one_to_one')
        assert len(x)==325 and x.same_pattern_new.equals(x.same_pattern_old)
        for metric in METRICS:
            delta=x[metric+'_new']-x[metric+'_old'];same=x.same_pattern_new
            groups=x.assign(delta=delta).groupby(['same_pattern_new','pair_group_new']).delta.mean()
            r=dict(variant=variant,metric=metric,old_pooled_same=float(x.loc[same,metric+'_old'].mean()),
                new_pooled_same=float(x.loc[same,metric+'_new'].mean()),
                delta_pooled_same=float(delta[same].mean()),
                same_pairs_improved=int((delta[same]>1e-12).sum()),same_pairs_worsened=int((delta[same]<-1e-12).sum()),
                same_pairs_unchanged=int((abs(delta[same])<=1e-12).sum()),
                delta_balanced_same=float(groups.loc[True].mean()),
                delta_balanced_different=float(groups.loc[False].mean()),
                delta_balanced_gap=float(groups.loc[True].mean()-groups.loc[False].mean()))
            changes.append(r)
        part=x[['run_x','run_y','same_pattern_new','pair_group_new']].copy();part['variant']=variant
        part['old_score']=x[MAIN+'_old'];part['new_score']=x[MAIN+'_new'];part['delta']=part.new_score-part.old_score
        all_changes.append(part)
    pd.DataFrame(changes).to_csv(ROOT/'results/band_change_summary.csv',index=False)
    pd.concat(all_changes,ignore_index=True).to_csv(ROOT/'results/band_change_pairs.csv',index=False)
    control=new[new.variant=='central_band20'].merge(old[old.variant=='central_band20'],
        on=['run_x','run_y'],suffixes=('_new','_old'),validate='one_to_one')
    error=float(np.max(abs(control[MAIN+'_new']-control[MAIN+'_old'])))
    assert error<1e-10, 'O controle fixo 30–49 deve reproduzir a análise anterior com os mesmos pesos.'
    principal={r['variant']:r for r in changes if r['metric']==MAIN}
    old_v5=old[(old.variant=='v5')&old.same_pattern][MAIN].mean()
    against={v:float(new[(new.variant==v)&new.same_pattern][MAIN].mean()-old_v5) for v in new.variant.unique()}
    result=dict(status='PASS',metric=MAIN,weights='Referência corrigida 20–59, pesos e pisos fixos entre faixas',
        baseline_change=principal['baseline'],variant_changes=principal,
        original_v5_pooled_same=float(old_v5),current_delta_vs_original_v5=against,
        fixed_control_30_49_score_max_error=error,
        caveat='Escores nativos 18–71 usam outros pesos; para atribuir diferenças à extração, consultar esta comparação fixa. Pares compartilham aquisições; avaliação interna, sem teste externo.')
    (ROOT/'results/band_comparison.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
