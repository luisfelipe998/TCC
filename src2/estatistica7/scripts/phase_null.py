"""Referência de fases para contraste WCO, preservando dependência entre pares."""
import time
from datetime import datetime,timezone
from scipy.stats import beta
from core import *


def phase_surrogate(x,rng):
    f=np.fft.rfft(x);phase=np.ones(len(f),dtype=complex)
    stop=len(f)-1 if len(x)%2==0 else len(f)
    phase[1:stop]=np.exp(1j*rng.uniform(0,2*np.pi,stop-1))
    before=np.fft.irfft(f*phase,n=len(x))
    # Mesmo detrending e padronização do diagnóstico arquivado; sem reinverter sinal.
    after=detrend(before,type='linear');after-=after.mean();after/=after.std()
    e0=abs(f)**2;e1=abs(np.fft.rfft(after))**2
    distortion=np.linalg.norm(e1-e0)/np.linalg.norm(e0)
    exact_error=np.linalg.norm(abs(np.fft.rfft(before))-abs(f))/np.linalg.norm(abs(f))
    return after,float(distortion),float(exact_error)


def primary_scores(co,px,py):
    valid=np.isfinite(co).all(axis=1)
    work=np.where(valid[:,None],co,0)
    weight=normalized_weights(px,py,valid,True)
    uniform=normalized_weights(px,py,valid,False)
    return np.array([np.nanmax(co),np.max(weight@work),uniform@work.mean(axis=1)])


def main():
    started=time.perf_counter()
    runs,labels,inv,raw=load_inputs()
    rr,ll,x,_=condition_inputs('3171_prefix',runs,labels,raw)
    select=ll!='dispersed bubbles';rr,ll,x=rr[select],ll[select],x[select]
    n=x.shape[1];wave=Wavelet(n)
    runtime=json.loads((ROOT/'configs/runtime.json').read_text());B=int(runtime['null_replicates'])
    metrics=CFG['null']['primary_metrics'];assert metrics==['global_max','common_lag_peak','map_mean']
    scopes=['churn_annular','plug_slug'];allpairs=np.array(list(combinations(range(len(rr)),2)))
    allowed=np.zeros(len(allpairs),bool);scope_pair_masks={}
    for scope in scopes:
        member=np.isin(ll,CFG['scopes'][scope]);mask=member[allpairs[:,0]]&member[allpairs[:,1]]
        allowed|=mask
    pairs=allpairs[allowed];assert len(pairs)==186
    for scope in scopes:
        member=np.isin(ll,CFG['scopes'][scope]);scope_pair_masks[scope]=member[pairs[:,0]]&member[pairs[:,1]]
    lag=np.arange(-n+1,n)
    t=lag[None,:]/(LAM[:,None]/FOURIER_FACTOR)
    zero_kernel=np.pi**(-.25)/np.sqrt(LAM[:,None]/FOURIER_FACTOR)*np.exp(-.5*t*t)*np.cos(6*t)
    # Denominadores centrais por soma temporal, economizando mapas de autocoerência descartados.
    zero_many=lambda xx: np.abs(np.array([correlation(v,v) for v in xx])@zero_kernel.T)
    with np.load(ROOT/'results/profiles_3171_prefix.npz',allow_pickle=False) as z:
        np.testing.assert_allclose(zero_many(x),z['zero_wps'][z['labels']!='dispersed bubbles'],rtol=1e-11,atol=1e-12)
    observed=csv_read(ROOT/'results/scores_3171_prefix.csv').set_index(['run_x','run_y'])
    observed_scores=np.array([[observed.loc[(rr[i],rr[j]),metric] for metric in metrics] for i,j in pairs])
    observed_stats={}
    for scope in scopes:
        mask=scope_pair_masks[scope];pp=pairs[mask]
        # A função recebe apenas rótulos das duas classes pelo índice dos pares.
        groups=CFG['scopes'][scope]
        for k,metric in enumerate(metrics):
            intra=[observed_scores[mask,k][(ll[pp[:,0]]==g)&(ll[pp[:,1]]==g)].mean() for g in groups]
            inter=observed_scores[mask,k][ll[pp[:,0]]!=ll[pp[:,1]]].mean()
            observed_stats[(scope,metric)]=float(np.mean(intra)-inter)
    rows=[];distortions=[];largest_spectral_error=0.
    for b in range(B):
        sx=[]
        for i,signal in enumerate(x):
            rng=np.random.default_rng(np.random.SeedSequence([CFG['seed'],8,b,i]))
            s,distortion,error=phase_surrogate(signal,rng);sx.append(s)
            distortions.append(dict(replicate=b,run=rr[i],relative_power_spectrum_change=distortion))
            largest_spectral_error=max(largest_spectral_error,error)
        sx=np.array(sx);p=zero_many(sx);values=[]
        for i,j in pairs:
            co,_=wave.cross(sx[i],sx[j],p[i],p[j]);values.append(primary_scores(co,p[i],p[j]))
        values=np.array(values)
        for scope in scopes:
            mask=scope_pair_masks[scope];pp=pairs[mask]
            groups=CFG['scopes'][scope]
            for k,metric in enumerate(metrics):
                intra=[values[mask,k][(ll[pp[:,0]]==g)&(ll[pp[:,1]]==g)].mean() for g in groups]
                inter=values[mask,k][ll[pp[:,0]]!=ll[pp[:,1]]].mean()
                rows.append(dict(replicate=b,scope=scope,metric=metric,balanced_contrast=float(np.mean(intra)-inter),
                                 mean_intra_a=float(intra[0]),mean_intra_b=float(intra[1]),mean_inter=float(inter)))
        if (b+1)%5==0 or b+1==B:
            elapsed=time.perf_counter()-started
            print(f'Nulo de fases: {b+1}/{B} réplicas de coorte; {elapsed:.1f}s; previsão restante {elapsed/(b+1)*(B-b-1):.1f}s',flush=True)
        if (b+1)%25==0:csv_save('phase_null_replicates.csv',rows)
    frame=csv_save('phase_null_replicates.csv',rows);csv_save('phase_null_spectral_distortion.csv',distortions)
    results=[]
    for scope in scopes:
        for metric in metrics:
            values=frame[(frame.scope==scope)&(frame.metric==metric)].balanced_contrast.to_numpy()
            obs=observed_stats[(scope,metric)];exceed=int(np.sum(values>=obs-1e-12));pval=(1+exceed)/(B+1)
            lower=0 if exceed==0 else beta.ppf(.025,exceed,B-exceed+1)
            upper=1 if exceed==B else beta.ppf(.975,exceed+1,B-exceed)
            results.append(dict(scope=scope,metric=metric,observed_contrast=obs,null_mean=float(values.mean()),
                                null_q025=float(np.quantile(values,.025)),null_q975=float(np.quantile(values,.975)),
                                exceedances=exceed,replicates=B,p_model=pval,mc_tail_ci_low=float(lower),mc_tail_ci_high=float(upper)))
    for row,padj in zip(results,holm([r['p_model'] for r in results])):row['p_model_holm']=float(padj)
    csv_save('phase_null_summary.csv',results)
    assert largest_spectral_error<1e-12
    dd=pd.DataFrame(distortions).relative_power_spectrum_change
    meta=dict(status='PASS',completed_utc=datetime.now(timezone.utc).isoformat(),replicates_requested=B,replicates_completed=B,
              cohort_pairs_per_replicate=len(pairs),tests=6,elapsed_seconds=time.perf_counter()-started,
              spectrum_error_before_redetrending=largest_spectral_error,
              relative_spectrum_change_after_redetrending_median=float(dd.median()),
              relative_spectrum_change_after_redetrending_q95=float(dd.quantile(.95)),
              relative_spectrum_change_after_redetrending_max=float(dd.max()),
              interpretation='Modelo de fases aproximado com espectros estimados nos dados e novo detrending; não é teste de rótulos, não prova ausência de informação, generalização ou causalidade. IC descreve erro Monte Carlo da cauda, não incerteza populacional.')
    json_write(ROOT/'results/phase_null_metadata.json',meta)
    print(pd.DataFrame(results).to_string(index=False),flush=True)


if __name__=='__main__':main()
