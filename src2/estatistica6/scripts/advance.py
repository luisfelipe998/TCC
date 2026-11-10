"""Aprofunda somente candidatos da triagem e os diagnósticos necessários."""
from itertools import combinations
from datetime import datetime, timezone
from math import comb
import time

from common import *
from screen import mean_profiles
from scipy.stats import spearmanr


def pearson_contrast(r, labels):
    iu = np.triu_indices(len(labels), 1)
    same = labels[iu[0]] == labels[iu[1]]
    return dict(n_runs=len(labels), n_intra=int(same.sum()), n_inter=int((~same).sum()),
                median_intra=float(np.median(r[iu][same])), median_inter=float(np.median(r[iu][~same])),
                delta=float(np.median(r[iu][same])-np.median(r[iu][~same])))


def a1(source, conditions):
    runs, labels = source['runs'], source['patterns']
    r = source['unit_descriptor'] @ source['unit_descriptor'].T
    out, influence = [], []
    masks = {'todos26': np.ones(26,dtype=bool), **{s:scope_mask(labels,s) for s in CONFIG['scopes']}}
    for scope, mask in masks.items():
        sub = r[np.ix_(mask,mask)]; rr=runs[mask]; ll=labels[mask]
        base = pearson_contrast(sub,ll)
        out.append(dict(scope=scope,condition='prefix',**base))
        for i, run in enumerate(rr):
            keep = np.arange(len(rr))!=i
            result = pearson_contrast(sub[np.ix_(keep,keep)],ll[keep])
            influence.append(dict(scope=scope,omitted_run=run,omitted_pattern=ll[i],delta_change=result['delta']-base['delta'],**result))
    for condition,data in conditions.items():
        if condition=='prefix':continue
        residual,_,_=residual_descriptor(data['profiles'],data['lam'])
        u=unit_rows(residual);mat=u@u.T
        for scope in CONFIG['scopes']:
            mask=scope_mask(data['labels'],scope)
            if set(data['labels'][mask])!=set(CONFIG['scopes'][scope]):continue
            out.append(dict(scope=scope,condition=condition,**pearson_contrast(mat[np.ix_(mask,mask)],data['labels'][mask])))
    combinations_rows=[]
    for i,g in enumerate(PATTERNS):
        for h in PATTERNS[i:]:
            vals=[r[a,b] for a,b in combinations(range(26),2) if sorted([labels[a],labels[b]])==sorted([g,h])]
            if not vals:continue
            combinations_rows.append(dict(pattern_x=g,pattern_y=h,n_pairs=len(vals),median_r=np.median(vals),min_r=np.min(vals),max_r=np.max(vals)))
    groups=pd.DataFrame(combinations_rows)
    repeated=groups[~groups.pattern_x.eq('dispersed bubbles') & ~groups.pattern_y.eq('dispersed bubbles')]
    balanced=float(repeated[repeated.pattern_x==repeated.pattern_y].median_r.mean()-repeated[repeated.pattern_x!=repeated.pattern_y].median_r.mean())
    save_csv('a1_contrastes.csv',out);save_csv('a1_leave_one_out.csv',influence);save_csv('a1_combinacoes.csv',groups)
    write_json(ROOT/'results/a1_balanced.json',dict(delta_balanced_four_classes=balanced,n_intra_class_medians=4,n_inter_class_medians=6))
    perturb=[];draws=[]
    for stream,scope in enumerate(['churn_annular','plug_slug']):
        mask=masks[scope];rr=r[np.ix_(mask,mask)];ll=labels[mask]
        iu=np.triu_indices(len(ll),1);same=ll[iu[0]]==ll[iu[1]]
        rng=np.random.default_rng(np.random.SeedSequence([CONFIG['seed'],1,stream]))
        deltas=[]
        for k in range(CONFIG['a1_weight_sensitivity']['replicates']):
            weight=rng.exponential(size=len(ll))
            for g in np.unique(ll):weight[ll==g]/=weight[ll==g].mean()
            pairweight=weight[iu[0]]*weight[iu[1]]
            delta=weighted_median(rr[iu][same],pairweight[same])-weighted_median(rr[iu][~same],pairweight[~same])
            deltas.append(delta);draws.append(dict(scope=scope,replicate=k,delta=delta))
        perturb.append(dict(scope=scope,replicates=len(deltas),q025=np.quantile(deltas,.025),median=np.median(deltas),q975=np.quantile(deltas,.975),positive_fraction=np.mean(np.array(deltas)>0),interval_type='sensibilidade à ponderação de ensaios, não IC'))
    save_csv('a1_weight_sensitivity.csv',perturb);save_csv('a1_weight_draws.csv',draws)


def a2(source, conditions):
    coordinates=[];rows=[];lam=source['wavelengths_px']
    for i,(run,n) in enumerate(zip(source['runs'],source['n_original'])):
        starts={'prefix':0,'center':(int(n)-3171)//2,'suffix':int(n)-3171}
        unique={}
        for position,start in starts.items():unique.setdefault(start,[]).append(position)
        for start,aliases in unique.items():
            coordinates.append(dict(run=run,pattern=source['patterns'][i],windows='|'.join(aliases),start=start,end_exclusive=start+3171,n_original=n))
        for (start1,aliases1),(start2,aliases2) in combinations(unique.items(),2):
            p1=conditions[aliases1[0]]['profiles'][i];p2=conditions[aliases2[0]]['profiles'][i]
            ds,f=distances(np.array([p1,p2]),lam,scaler=(np.zeros(3),np.ones(3)))
            overlap=max(0,min(start1+3171,start2+3171)-max(start1,start2))/3171
            rows.append(dict(run=run,pattern=source['patterns'][i],window_a=aliases1[0],window_b=aliases2[0],overlap_fraction=overlap,
                             pearson_r=1-ds['pearson_residual'][0,1]**2/2,hellinger=ds['hellinger'][0,1],wasserstein=ds['wasserstein'][0,1],
                             absolute_log_peak_change=abs(f.log_peak.iloc[0]-f.log_peak.iloc[1]),absolute_width_change=abs(f.width.iloc[0]-f.width.iloc[1])))
    assert len(coordinates)==76 and len(rows)==75
    save_csv('a2_window_coordinates.csv',coordinates);save_csv('a2_within_run_windows.csv',rows)


def a3_a6(source, inventory, feature_frame):
    effects=[];within_q=[];association=[]
    names=['log_peak','slope','width','centroid','entropy','band_short','band_middle','band_long']
    joined=feature_frame.merge(inventory[['run','m_code','q_code','n','recomputed_shift_mean']],on='run',validate='many_to_one')
    for condition, frame in joined.groupby('condition',sort=False):
        for scope in ['churn_annular','plug_slug']:
            a,b=CONFIG['scopes'][scope]; sub=frame[frame.pattern.isin([a,b])]
            if set(sub.pattern)!={a,b}:continue
            for feature in names:
                x=sub.loc[sub.pattern==a,feature].to_numpy();y=sub.loc[sub.pattern==b,feature].to_numpy()
                effects.append(dict(scope=scope,condition=condition,feature=feature,pattern_a=a,pattern_b=b,n_a=len(x),n_b=len(y),median_a=np.median(x),median_b=np.median(y),median_difference_a_minus_b=np.median(x)-np.median(y),cliff_delta=cliff(x,y)))
                for q, group in sub.groupby('q_code'):
                    xx=group.loc[group.pattern==a,feature].to_numpy();yy=group.loc[group.pattern==b,feature].to_numpy()
                    if len(xx) and len(yy):within_q.append(dict(scope=scope,condition=condition,feature=feature,q_code=q,n_a=len(xx),n_b=len(yy),median_difference_a_minus_b=np.median(xx)-np.median(yy),cliff_delta=cliff(xx,yy)))
                for covariate in ['m_code','q_code','n','recomputed_shift_mean']:
                    rho=spearmanr(sub[feature],sub[covariate]).statistic
                    association.append(dict(scope=scope,condition=condition,feature=feature,covariate=covariate,spearman=rho,n_runs=len(sub)))
    save_csv('a3_effects.csv',effects);save_csv('a6_within_q_effects.csv',within_q);save_csv('a6_associations.csv',association)
    save_csv('a6_coverage.csv',inventory.groupby(['m_code','q_code','pattern']).size().rename('n_runs').reset_index())
    coverage=[]
    for scope in ['churn_annular','plug_slug']:
        frame=inventory[inventory.pattern.isin(CONFIG['scopes'][scope])]
        a,b=CONFIG['scopes'][scope]
        counts=frame.groupby(['q_code','pattern']).size().unstack(fill_value=0).reindex(columns=[a,b],fill_value=0)
        q_arrangements=int(np.prod([comb(int(row.sum()),int(row[a])) for _,row in counts.iterrows()]))
        mq_arrangements=1
        for _,group in frame.groupby(['m_code','q_code']):mq_arrangements*=comb(len(group),int((group.pattern==a).sum()))
        coverage.append(dict(scope=scope,n_runs=len(frame),q_arrangements=q_arrangements,mq_arrangements=mq_arrangements,
                             length_shift_pearson=float(np.corrcoef(frame.n,frame.recomputed_shift_mean)[0,1]),
                             shared_q=[int(q) for q,row in counts.iterrows() if (row>0).all()],
                             tests_performed=False))
    write_json(ROOT/'results/a6_design.json',coverage)


def a7(source):
    u=source['unit_descriptor'];labels=source['patterns'];rows=[]
    for g in PATTERNS:
        mask=labels==g;n=int(mask.sum())
        if n<2:continue
        centroid=u[mask].mean(axis=0)
        radius=np.linalg.norm(u[mask]-centroid,axis=1)
        for run,a in zip(source['runs'][mask],radius):rows.append(dict(run=run,pattern=g,n_class=n,radius=a,adjusted_radius=np.sqrt(n/(n-1))*a))
    frame=save_csv('a7_radii.csv',rows)
    save_csv('a7_summary.csv',frame.groupby('pattern').agg(n_runs=('run','size'),mean_radius=('radius','mean'),mean_adjusted_radius=('adjusted_radius','mean'),median_adjusted_radius=('adjusted_radius','median')).reset_index())


def main():
    source,inventory,conditions=load_inputs()
    screen=read_csv(ROOT/'results/triagem.csv')
    with np.load(ROOT/'results/distances_screen.npz',allow_pickle=False) as z:matrices={k:z[k].copy() for k in z.files}
    advanced=[];individual=[]
    candidates=screen[(screen.condition=='prefix') & bools(screen.prefix_gate)].copy()
    assert not screen[(screen.family=='wco_cached') & bools(screen.prefix_gate)].shape[0], 'WCO passou: exige implementação fixa antes da conclusão.'
    # Somente os participantes do contraste que justificou a extensão de mean power.
    mean_candidates=candidates[candidates.family=='mean_power']
    if len(mean_candidates):
        needed=sorted(set(g for scope in mean_candidates.scope for g in CONFIG['scopes'][scope]),key=PATTERNS.index)
        mask=np.isin(source['patterns'],needed);runs=source['runs'][mask];labels=source['patterns'][mask]
        raw={}
        for run in runs:
            with np.load(V16/'intermediate'/f'{run}_input.npz',allow_pickle=False) as z:raw[run]=z['raw'].copy()
        for condition in ['center','suffix']:
            profiles=mean_profiles(raw,runs,source['wavelengths_px'],condition)
            np.savez_compressed(ROOT/f'results/mean_power_{condition}.npz',runs=runs,labels=labels,wavelengths_px=source['wavelengths_px'],profiles=profiles)
            for metric,d in distances(profiles,source['wavelengths_px'],metrics=CONFIG['mean_power_candidates'])[0].items():
                matrices[f'mean_power__{condition}__{metric}']=d
                for scope in mean_candidates.scope.unique():
                    sm=scope_mask(labels,scope);sub=d[np.ix_(sm,sm)]
                    advanced.append(dict(family='mean_power',condition=condition,metric=metric,scope=scope,**metric_summary(sub,labels[sm])))
                    for run,g,s,nn in zip(runs[sm],labels[sm],silhouette_by_run(sub,labels[sm]),nearest_agreement(sub,labels[sm])):
                        individual.append(dict(family='mean_power',condition=condition,metric=metric,scope=scope,run=run,pattern=g,silhouette=s,neighbor_agreement=nn))
    save_csv('extensao_mean_power.csv',advanced);save_csv('extensao_mean_power_por_ensaio.csv',individual)
    combined=pd.concat([screen,pd.DataFrame(advanced)],ignore_index=True)
    decisions=[];influence=[];clusters=[];cluster_summary=[]
    for rec in candidates.itertuples():
        prefix=matrices[f'{rec.family}__prefix__{rec.metric}']
        mask=scope_mask(source['patterns'],rec.scope)
        d=prefix[np.ix_(mask,mask)];runs=source['runs'][mask];labels=source['patterns'][mask]
        loo=[]
        for i,run in enumerate(runs):
            keep=np.arange(len(runs))!=i
            stats=metric_summary(d[np.ix_(keep,keep)],labels[keep])
            loo.append(stats['macro_silhouette'])
            influence.append(dict(family=rec.family,metric=rec.metric,scope=rec.scope,omitted_run=run,**stats))
        rows=combined[(combined.family==rec.family)&(combined.metric==rec.metric)&(combined.scope==rec.scope)]
        reasons=[]
        for condition in ['center','suffix']:
            row=rows[rows.condition==condition].iloc[0]
            if row.macro_silhouette<=0:reasons.append(f'silhueta não positiva em {condition}')
            if row.neighbor_lift<=1e-12:reasons.append(f'vizinhança sem vantagem em {condition}')
        if min(loo)<=0:reasons.append('silhueta não positiva ao retirar ensaio')
        if rec.scope=='churn_annular':
            control=rows[rows.condition=='7001_prefix']
            if len(control)!=1 or control.iloc[0].macro_silhouette<=0:reasons.append('controle 7001 sem sinal positivo')
        decisions.append(dict(family=rec.family,metric=rec.metric,scope=rec.scope,retained=not reasons,loo_silhouette_min=min(loo),loo_silhouette_max=max(loo),
                              reasons='; '.join(reasons) if reasons else 'passou os critérios exploratórios de estabilidade'))
        # Agrupamento fixo, sem rótulos no ajuste; compara com rótulos somente ao resumir.
        for condition in ['prefix','center','suffix']:
            full=matrices[f'{rec.family}__{condition}__{rec.metric}']
            if len(full)==26:
                sub=full[np.ix_(mask,mask)]
            else:
                assert len(full)==len(runs)
                sub=full
            assigned,ari,tree=clustering(sub,labels)
            cluster_summary.append(dict(family=rec.family,metric=rec.metric,scope=rec.scope,condition=condition,k=len(np.unique(labels)),ari=ari))
            for run,g,cl in zip(runs,labels,assigned):clusters.append(dict(family=rec.family,metric=rec.metric,scope=rec.scope,condition=condition,run=run,pattern=g,cluster=int(cl)))
    decision_frame=save_csv('decisoes_finais.csv',decisions)
    save_csv('candidatos_leave_one_out.csv',influence);save_csv('clusters.csv',clusters);save_csv('cluster_summary.csv',cluster_summary)
    np.savez_compressed(ROOT/'results/distances_final.npz',**matrices)
    save_csv('metricas_finais.csv',combined)
    a1(source,conditions);a2(source,conditions)
    a3_a6(source,inventory,read_csv(ROOT/'results/a3_features.csv'))
    a7(source)
    # Comparação funcional descritiva: todas as escalas, sem selecionar máximo favorável.
    functional=[]
    for scope in ['churn_annular','plug_slug']:
        a,b=CONFIG['scopes'][scope]
        for condition,data in conditions.items():
            if a not in data['labels'] or b not in data['labels']:continue
            p=data['profiles'];logp=np.log(p+1e-12)
            aa=logp[data['labels']==a];bb=logp[data['labels']==b]
            for k,lam in enumerate(data['lam']):
                functional.append(dict(scope=scope,condition=condition,wavelength_px=lam,mean_log_a=aa[:,k].mean(),mean_log_b=bb[:,k].mean(),mean_log_difference=aa[:,k].mean()-bb[:,k].mean(),sd_a=aa[:,k].std(ddof=1),sd_b=bb[:,k].std(ddof=1)))
    save_csv('a4_functional_descriptive.csv',functional)
    # Resumos A5 mantêm todos os padrões, inclusive o singleton entre classes.
    pairs=[]
    for metric in ['hellinger','wasserstein']:
        d=matrices[f'wps__prefix__{metric}']
        for i,j in combinations(range(26),2):pairs.append(dict(metric=metric,run_x=source['runs'][i],run_y=source['runs'][j],pattern_x=source['patterns'][i],pattern_y=source['patterns'][j],same_pattern=source['patterns'][i]==source['patterns'][j],distance=d[i,j]))
    save_csv('a5_pairs.csv',pairs)
    manifest=json.loads((ROOT/'results/manifest_started.json').read_text())
    assert all(sha(SRC/k)==v for k,v in manifest['sources'].items())
    write_json(ROOT/'results/analysis_completed.json',dict(completed_utc=datetime.now(timezone.utc).isoformat(),source_files_unchanged=True,
                                                        config_sha256=sha(ROOT/'configs/triagem_v1.json'),retained_candidates=int(decision_frame.retained.sum()),
                                                        new_null_replicates=0,formal_tests=0,spatial_bootstraps=0))
    print(decision_frame.to_string(index=False),flush=True)


if __name__=='__main__':main()
