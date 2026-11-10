"""Controles de normalização, suporte e ruído, identificados como exploração adicional."""
from scipy.signal import lfilter
from core import *
import time


def normalized_kernels():
    frames=[];classes=[];clusters=[];influence=[];checks=[];matrices={}
    for condition in CFG['conditions']:
        with np.load(ROOT/f'results/profiles_{condition}.npz',allow_pickle=False) as z:
            runs=z['runs'];labels=z['labels'];pairs=z['pairs'];area=z['self_area'];means=z['mean'];zeros=z['zero']
        candidates={}
        for name,band in {'full':np.ones(len(LAM),bool),**BANDS}.items():
            weight=H[band]/H[band].sum()
            k=np.diag(area[:,band]@weight)
            pairvalues=means[:,band]@weight
            for (i,j),v in zip(pairs,pairvalues):k[i,j]=k[j,i]=v
            candidates['area_normalized_'+name]=k
        k=np.eye(len(runs))
        for (i,j),value in zip(pairs,zeros@(H/H.sum())):k[i,j]=k[j,i]=value
        candidates['zero_uniform_kernel']=k
        for metric,k in candidates.items():
            minimum=float(np.linalg.eigvalsh(k).min())
            normalized=k/np.sqrt(np.outer(np.diag(k),np.diag(k)))
            excess=float(normalized.max()-1)
            valid=minimum>=-1e-10 and excess<=1e-10
            checks.append(dict(condition=condition,metric=metric,minimum_eigenvalue=minimum,maximum_above_one=excess,geometry_valid=valid))
            if not valid:continue
            # Somente arredondamento de diagonal; sem reparo de autovalores ou seleção de eixos.
            d=np.sqrt(np.maximum(0,2*(1-normalized)));np.fill_diagonal(d,0)
            matrices[f'{condition}__{metric}']=d
            for scope,patterns in CFG['scopes'].items():
                mask=np.isin(labels,patterns)
                if set(labels[mask])!=set(patterns):continue
                sub=d[np.ix_(mask,mask)];ll=labels[mask];rr=runs[mask]
                stats=metric_summary(sub,ll)
                frames.append(dict(condition=condition,metric=metric,family='WCO_normalizada_adendo',scope=scope,**stats,screen_pass=screen_pass(stats)))
                s=silhouette(sub,ll)
                for g in patterns:classes.append(dict(condition=condition,metric=metric,scope=scope,pattern=g,mean_silhouette=float(s[ll==g].mean())))
                _,ar=cluster(sub,ll);clusters.append(dict(condition=condition,metric=metric,scope=scope,k=len(patterns),ari=ar))
                if condition=='3171_prefix':
                    for i,run in enumerate(rr):
                        keep=np.arange(len(rr))!=i
                        loo=metric_summary(sub[np.ix_(keep,keep)],ll[keep])
                        influence.append(dict(metric=metric,scope=scope,omitted_run=run,macro_silhouette=loo['macro_silhouette']))
    pd.DataFrame(checks).to_csv(ROOT/'validation/kernel_geometry.csv',index=False,float_format='%.17g')
    f=csv_save('metrics_addendum.csv',frames);csv_save('classes_addendum.csv',classes);csv_save('clustering_addendum.csv',clusters)
    inf=csv_save('influence_addendum.csv',influence)
    np.savez_compressed(ROOT/'results/distance_addendum.npz',**matrices)
    decisions=[]
    for row in f[f.condition=='3171_prefix'].itertuples():
        sub=f[(f.metric==row.metric)&(f.scope==row.scope)]
        loo=inf[(inf.metric==row.metric)&(inf.scope==row.scope)]
        passes=bool(row.screen_pass and (sub.min_class_silhouette>0).all() and (sub.neighbor_lift>1e-12).all() and loo.macro_silhouette.min()>0)
        decisions.append(dict(metric=row.metric,scope=row.scope,screen_pass=row.screen_pass,robust_pass=passes,
                              min_class_silhouette_all_conditions=float(sub.min_class_silhouette.min()),min_macro_silhouette_loo=float(loo.macro_silhouette.min())))
    csv_save('decisions_addendum.csv',decisions)
    print('NORMALIZAÇÃO: triagem\n'+f[(f.condition=='3171_prefix')&f.screen_pass].to_string(index=False),flush=True)
    print('NORMALIZAÇÃO: estabilidade completa\n'+pd.DataFrame(decisions)[lambda z:z.robust_pass].to_string(index=False),flush=True)


def noise_controls():
    ad=json.loads((ROOT/'configs/diagnostics_addendum.json').read_text())['noise_diagnostic']
    rows=[];started=time.perf_counter()
    for n in ad['lengths']:
        wave=Wavelet(n)
        for model_id,model in enumerate(ad['models']):
            for replicate in range(ad['independent_pairs_each']):
                rng=np.random.default_rng(np.random.SeedSequence([CFG['seed'],9,n,model_id,replicate]))
                data=rng.normal(size=(2,n+ad['burnin_ar1']))
                if model=='gaussian_ar1_095':data=lfilter([1],[1,-.95],data,axis=1)
                x,y=[prepare(v[-n:]) for v in data]
                px,_=wave.auto(x);py,_=wave.auto(y);co,_=wave.cross(x,y,px,py)
                peak=np.unravel_index(np.argmax(co),co.shape)
                weight=normalized_weights(px,py,np.ones(len(LAM),bool))
                rows.append(dict(n=n,model=model,replicate=replicate,global_max=float(co[peak]),peak_lambda_px=float(LAM[peak[0]]),
                                 common_lag_peak=float(np.max(weight@co)),map_mean=float((H/H.sum())@co.mean(axis=1)),
                                 maximum_up_to_800px=float(co[LAM<=800].max())))
            print(f'Ruído independente: {model}, N={n}, {ad["independent_pairs_each"]} pares; {time.perf_counter()-started:.1f}s',flush=True)
    csv_save('independent_noise_controls.csv',rows)
    summary=[]
    for (n,model),f in pd.DataFrame(rows).groupby(['n','model']):
        summary.append(dict(n=int(n),model=model,pairs=len(f),median_max=float(f.global_max.median()),q05_max=float(f.global_max.quantile(.05)),q95_max=float(f.global_max.quantile(.95)),
                            fraction_max_above09=float((f.global_max>=.9).mean()),fraction_peak_upper_boundary=float(np.isclose(f.peak_lambda_px,3000).mean()),
                            median_max_800=float(f.maximum_up_to_800px.median())))
    csv_save('independent_noise_summary.csv',summary)
    print(pd.DataFrame(summary).to_string(index=False),flush=True)


if __name__=='__main__':
    normalized_kernels()
    noise_controls()
