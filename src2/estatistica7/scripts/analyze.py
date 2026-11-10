"""Compara padrões sem tratar pares ou escalas como aquisições independentes."""
from core import *
from scipy.stats import spearmanr


def band_features(area):
    return np.array([np.log(area)[:,band]@(H[band]/H[band].sum()) for band in BANDS.values()]).T


def get_distances(scores,z,scaler):
    runs=z['runs'];pairs=z['pairs'];out={};index={r:i for i,r in enumerate(runs)}
    assert len(scores)==len(pairs)
    for metric in CFG['cross_candidates']:
        d=np.zeros((len(runs),len(runs)))
        for row in scores.itertuples():
            i,j=index[row.run_x],index[row.run_y]
            d[i,j]=d[j,i]=1-getattr(row,metric)
        out[metric]=d
    for name in ['area','width']:
        x=np.log(z['self_'+name])*np.sqrt(H/H.sum())
        out[f'self_{name}_log_l2']=squareform(pdist(x))
    bands=band_features(z['self_area'])
    out['self_area_bands']=squareform(pdist((bands-scaler[0])/scaler[1]))
    return out,bands


def mean_contrast_profile(values,pairs,labels):
    groups=np.unique(labels[np.unique(pairs)]);within=[];between=[]
    for k,g in enumerate(groups):
        mask=(labels[pairs[:,0]]==g)&(labels[pairs[:,1]]==g)
        assert mask.any();within.append(values[mask].mean(axis=0))
        for h in groups[k+1:]:
            mask=((labels[pairs[:,0]]==g)&(labels[pairs[:,1]]==h))|((labels[pairs[:,0]]==h)&(labels[pairs[:,1]]==g))
            between.append(values[mask].mean(axis=0))
    return np.mean(within,axis=0)-np.mean(between,axis=0)


def main():
    inv=csv_read(ROOT/'results/inventory.csv').set_index('run')
    with np.load(ROOT/'results/profiles_3171_prefix.npz',allow_pickle=False) as z:
        f=band_features(z['self_area'][z['labels']!='dispersed bubbles'])
    scaler=(np.median(f,axis=0),np.quantile(f,.75,axis=0)-np.quantile(f,.25,axis=0))
    assert np.all(scaler[1]>1e-12)
    json_write(ROOT/'results/self_scaler.json',dict(center=scaler[0],iqr=scaler[1]))
    metrics=[];class_rows=[];individual=[];clusters=[];influence=[];covariates=[];within_q=[];functional=[];self_rows=[]
    all_matrices={};cache={};pair_relations=[]
    for condition in CFG['conditions']:
        with np.load(ROOT/f'results/profiles_{condition}.npz',allow_pickle=False) as zz:z={k:zz[k].copy() for k in zz.files}
        scores=csv_read(ROOT/f'results/scores_{condition}.csv')
        rr,ll,pairs=z['runs'],z['labels'],z['pairs'];cache[condition]=z
        ds,bands=get_distances(scores,z,scaler)
        for i,run in enumerate(rr):
            self_rows.append(dict(condition=condition,run=run,pattern=ll[i],
                log_area_short=bands[i,0],log_area_middle=bands[i,1],log_area_long=bands[i,2],
                mean_area=float((H/H.sum())@z['self_area'][i]),
                mean_width_lambda=float((H/H.sum())@z['self_width'][i]),
                censored_scales=int(z['self_width_censored'][i].sum())))
        for metric,d in ds.items():
            all_matrices[f'{condition}__{metric}']=d
            family='autocoerencia' if metric.startswith('self_') else 'WCO_cruzada'
            for scope,patterns in CFG['scopes'].items():
                mask=np.isin(ll,patterns)
                if set(ll[mask])!=set(patterns):continue
                sub=d[np.ix_(mask,mask)];labels=ll[mask];runs=rr[mask]
                stats=metric_summary(sub,labels)
                statsrow=dict(condition=condition,metric=metric,family=family,scope=scope,**stats)
                statsrow['screen_pass']=screen_pass(statsrow)
                metrics.append(statsrow)
                s=silhouette(sub,labels);nn=neighbor(sub,labels)
                for run,label,si,ne in zip(runs,labels,s,nn):
                    individual.append(dict(condition=condition,metric=metric,scope=scope,run=run,pattern=label,silhouette=si,neighbor_agreement=ne))
                for g in patterns:class_rows.append(dict(condition=condition,metric=metric,scope=scope,pattern=g,n_runs=int(sum(labels==g)),mean_silhouette=float(s[labels==g].mean()),mean_neighbor=float(nn[labels==g].mean())))
                assigned,score=cluster(sub,labels)
                clusters.append(dict(condition=condition,metric=metric,scope=scope,ari=score,k=len(patterns)))
                if condition=='3171_prefix':
                    for omitted,run in enumerate(runs):
                        keep=np.arange(len(runs))!=omitted
                        result=metric_summary(sub[np.ix_(keep,keep)],labels[keep])
                        influence.append(dict(metric=metric,scope=scope,omitted_run=run,macro_silhouette=result['macro_silhouette'],min_class_silhouette=result['min_class_silhouette']))
                # Controle simétrico de q só é estimável com >=2 ensaios de cada classe.
                if len(patterns)==2:
                    q=inv.loc[runs,'q_code'].to_numpy()
                    for value in np.unique(q):
                        qm=q==value;counts={g:int(sum(labels[qm]==g)) for g in patterns}
                        estimable=min(counts.values())>=2
                        rec=dict(condition=condition,metric=metric,scope=scope,q_code=int(value),n_a=counts[patterns[0]],n_b=counts[patterns[1]],estimable=estimable)
                        if estimable:rec.update(metric_summary(sub[np.ix_(qm,qm)],labels[qm]))
                        within_q.append(rec)
        # Redundância descritiva das proximidades com WPS; sem p para linhas de pares.
        p=z['zero_wps']*H;p=p/p.sum(axis=1,keepdims=True)
        hdist=squareform(pdist(np.sqrt(p)))/np.sqrt(2)
        raw_inv=inv.loc[rr]
        for scope in CFG['primary_scopes']:
            pm=np.isin(ll,CFG['scopes'][scope])
            if set(ll[pm])!=set(CFG['scopes'][scope]):continue
            i,j=np.triu_indices(pm.sum(),1)
            hh=hdist[np.ix_(pm,pm)][i,j]
            for metric,d in ds.items():
                sub=d[np.ix_(pm,pm)][i,j]
                if np.std(sub)==0:continue
                pair_relations.append(dict(condition=condition,scope=scope,metric=metric,relation='distance_WPS_hellinger',spearman=float(spearmanr(sub,hh).statistic),n_pairs=len(sub)))
                for variable in ['m_code','q_code','n','recomputed_shift_mean']:
                    cov=raw_inv.loc[rr[pm],variable].to_numpy()
                    dif=abs(cov[i]-cov[j])
                    if np.std(dif)>0:pair_relations.append(dict(condition=condition,scope=scope,metric=metric,relation='absolute_difference_'+variable,spearman=float(spearmanr(sub,dif).statistic),n_pairs=len(sub)))
        # Perfis funcionais e envelope de influência por ensaio, não IC nem banda de significância.
        for scope in CFG['primary_scopes']:
            mask=np.isin(ll,CFG['scopes'][scope])
            if set(ll[mask])!=set(CFG['scopes'][scope]):continue
            selected=mask[pairs[:,0]]&mask[pairs[:,1]]
            selected_pairs=pairs[selected]
            for observable in ['mean','peak','zero']:
                vals=z[observable][selected]
                obs=mean_contrast_profile(vals,selected_pairs,ll)
                loo=[]
                for omitted in np.flatnonzero(mask):
                    keep=(selected_pairs[:,0]!=omitted)&(selected_pairs[:,1]!=omitted)
                    # labels só das classes do escopo; a função recebe o vetor global e infere grupos dos pares.
                    loo.append(mean_contrast_profile(vals[keep],selected_pairs[keep],ll))
                lo=np.min(loo,axis=0);hi=np.max(loo,axis=0)
                for k,lam in enumerate(LAM):functional.append(dict(condition=condition,scope=scope,observable=observable,wavelength_px=lam,balanced_contrast=obs[k],loo_min=lo[k],loo_max=hi[k]))
    metric_frame=csv_save('metrics.csv',metrics);csv_save('metrics_by_class.csv',class_rows);csv_save('metrics_by_run.csv',individual)
    csv_save('clustering.csv',clusters);inf=csv_save('influence.csv',influence)
    csv_save('within_q.csv',within_q);csv_save('pair_associations.csv',pair_relations);csv_save('functional.csv',functional)
    sf=csv_save('self_features.csv',self_rows)
    for (condition,scope) in metric_frame[['condition','scope']].drop_duplicates().itertuples(index=False,name=None):
        sub=sf[(sf.condition==condition)&sf.pattern.isin(CFG['scopes'][scope])].merge(inv.reset_index()[['run','m_code','q_code','n','recomputed_shift_mean']],on='run',validate='one_to_one')
        for feat in ['mean_area','mean_width_lambda','log_area_short','log_area_middle','log_area_long']:
            for cov in ['m_code','q_code','n','recomputed_shift_mean']:
                covariates.append(dict(condition=condition,scope=scope,feature=feat,covariate=cov,spearman=float(spearmanr(sub[feat],sub[cov]).statistic),n_runs=len(sub)))
    csv_save('self_covariates.csv',covariates)
    np.savez_compressed(ROOT/'results/distance_matrices.npz',**all_matrices)
    decisions=[]
    for row in metric_frame[metric_frame.condition=='3171_prefix'].itertuples():
        reasons=[]
        if not row.screen_pass:reasons.append('não passou triagem inicial')
        relevant=[c for c in CFG['conditions'] if not c.startswith('7001') or row.scope=='churn_annular']
        subset=metric_frame[(metric_frame.metric==row.metric)&(metric_frame.scope==row.scope)&metric_frame.condition.isin(relevant)]
        for rr in subset.itertuples():
            if rr.min_class_silhouette<=0:reasons.append('classe sem compactação em '+rr.condition)
            if rr.neighbor_lift<=1e-12:reasons.append('vizinhança sem vantagem em '+rr.condition)
        loo=inf[(inf.metric==row.metric)&(inf.scope==row.scope)]
        if loo.macro_silhouette.min()<=0:reasons.append('silhueta macro não positiva após retirar ensaio')
        decisions.append(dict(metric=row.metric,family=row.family,scope=row.scope,screen_pass=row.screen_pass,
                              robust_pass=not reasons,min_class_silhouette_all_conditions=float(subset.min_class_silhouette.min()),
                              min_macro_silhouette_loo=float(loo.macro_silhouette.min()),reasons='; '.join(reasons) if reasons else 'passou os critérios exploratórios pré-fixados'))
    dec=csv_save('decisions.csv',decisions)
    print('TRIAGEM INICIAL\n'+metric_frame[(metric_frame.condition=='3171_prefix')&metric_frame.screen_pass][['metric','scope','macro_silhouette','min_class_silhouette','balanced_neighbor']].to_string(index=False),flush=True)
    print('CRITÉRIO COMPLETO\n'+dec[dec.robust_pass].to_string(index=False),flush=True)


if __name__=='__main__':main()
