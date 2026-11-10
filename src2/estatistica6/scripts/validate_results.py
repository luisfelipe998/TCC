"""Conferência independente dos resultados exportados e da integridade das fontes."""
from common import *
from scipy.signal import fftconvolve


def main():
    source,inventory,conditions=load_inputs()
    results=ROOT/'results'
    numeric=json.loads((ROOT/'validation/numerical.json').read_text())
    assert numeric['status']=='PASS'
    r=source['unit_descriptor']@source['unit_descriptor'].T
    influence=read_csv(results/'a1_leave_one_out.csv')
    max_error=0.
    for row in influence.itertuples():
        members=[i for i,run in enumerate(source['runs']) if run!=row.omitted_run and (row.scope=='todos26' or source['patterns'][i] in CONFIG['scopes'][row.scope])]
        intra=[];inter=[]
        for i,j in combinations(members,2):
            (intra if source['patterns'][i]==source['patterns'][j] else inter).append(r[i,j])
        delta=np.median(intra)-np.median(inter)
        max_error=max(max_error,abs(delta-row.delta))
        assert len(intra)==row.n_intra and len(inter)==row.n_inter
    assert max_error<1e-13
    mp_errors=[]
    # Convolução temporal independente confere mean power em ambos os extremos da banda.
    for condition in ['prefix','center','suffix']:
        with np.load(results/f'mean_power_{condition}.npz',allow_pickle=False) as z:
            run=z['runs'][0];profile=z['profiles'][0];lam=z['wavelengths_px']
        with np.load(V16/'intermediate'/f'{run}_input.npz',allow_pickle=False) as z:raw=z['raw']
        start=0 if condition=='prefix' else ((len(raw)-3171)//2 if condition=='center' else len(raw)-3171)
        corr=correlation(prepare(raw[start:start+3171]),prepare(raw[start:start+3171]))
        lag=np.arange(-len(corr)+1,len(corr))
        for k in [0,99,199]:
            a=lam[k]/FOURIER_FACTOR;t=lag/a
            kernel=np.pi**(-.25)/np.sqrt(a)*np.exp(-.5*t*t)*np.exp(6j*t)
            temporal=fftconvolve(corr,kernel,mode='same')
            actual=abs(temporal).mean()
            error=abs(actual-profile[k])
            assert error<1e-10
            mp_errors.append(error)
    with np.load(results/'mean_power_prefix.npz',allow_pickle=False) as z:
        i=list(z['runs']).index('m07q01_001');saved=z['profiles'][i]
    archived=read_csv(V16/'results/profiles.csv')
    archived=archived[(archived.run=='m07q01_001')&(archived.method=='v16_biased')].sort_values('wavelength_px').lag_mean_all.to_numpy()
    np.testing.assert_allclose(saved,archived,rtol=1e-11,atol=1e-12)
    assert len(read_csv(results/'a2_window_coordinates.csv'))==76
    windows=read_csv(results/'a2_within_run_windows.csv')
    assert len(windows)==75 and not (windows.run=='m07q01_001').any()
    assert windows.overlap_fraction.between(0,1).all()
    design=json.loads((results/'a6_design.json').read_text())
    assert [x['q_arrangements'] for x in design]==[200,9]
    assert all(x['mq_arrangements']==1 for x in design)
    draws=read_csv(results/'a1_weight_draws.csv')
    assert draws.groupby('scope').size().to_dict()=={'churn_annular':5000,'plug_slug':5000}
    for row in read_csv(results/'a1_weight_sensitivity.csv').itertuples():
        data=draws[draws.scope==row.scope].delta
        np.testing.assert_allclose([row.q025,row.median,row.q975],np.quantile(data,[.025,.5,.975]))
    # ARI por contagem literal de pares concordantes, sem reutilizar o agrupador.
    clusters=read_csv(results/'clusters.csv');summary=read_csv(results/'cluster_summary.csv')
    ari_error=0.
    for row in summary.itertuples():
        data=clusters[(clusters.family==row.family)&(clusters.metric==row.metric)&(clusters.scope==row.scope)&(clusters.condition==row.condition)]
        pairs=list(combinations(range(len(data)),2))
        same_a=np.array([data.pattern.iloc[i]==data.pattern.iloc[j] for i,j in pairs])
        same_b=np.array([data.cluster.iloc[i]==data.cluster.iloc[j] for i,j in pairs])
        expected=same_a.sum()*same_b.sum()/len(pairs)
        ari=((same_a&same_b).sum()-expected)/((same_a.sum()+same_b.sum())/2-expected)
        ari_error=max(ari_error,abs(ari-row.ari))
    assert ari_error<1e-12
    allmetrics=read_csv(results/'metricas_finais.csv')
    assert not ((allmetrics.scope=='quatro_padroes')&(allmetrics.n_runs!=25)).any()
    assert not ((allmetrics.scope=='plug_slug')&(allmetrics.n_runs!=6)).any()
    assert not ((allmetrics.scope=='churn_annular')&(allmetrics.n_runs!=19)).any()
    manifest=json.loads((results/'manifest_started.json').read_text())
    assert sha(ROOT/'configs/triagem_v1.json')==manifest['config_sha256']
    assert all(sha(SRC/k)==v for k,v in manifest['sources'].items())
    out=dict(status='PASS',sources_checked=len(manifest['sources']),sources_unchanged=True,
             a1_independent_leave_one_out_max_error=max_error,mean_power_temporal_max_error=max(mp_errors),
             ari_independent_max_error=ari_error,valid_window_pairs=len(windows),
             weighting_replicates=10000,condition_memberships_verified=True,
             q_arrangements=[200,9],mq_arrangements=[1,1])
    write_json(ROOT/'validation/final_numeric.json',out)
    print(json.dumps(out,indent=2,ensure_ascii=False))


if __name__=='__main__':main()
