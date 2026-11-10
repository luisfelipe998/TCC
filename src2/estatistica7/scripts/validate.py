"""Auditoria independente do cálculo WCO, dos resumos, da unidade de análise e do nulo."""
from core import *


def main():
    runs,labels,inv,raw=load_inputs();errors=[]
    cases=[('3171_prefix','m03q05_001','m03q10_001'),('3171_center','m05q01_001','m04q01_001'),
           ('5071_suffix','m06q01_001','m06q05_001'),('7001_center','m04q20_001','m04q80_001')]
    for condition,a,b in cases:
        rr,ll,x,_=condition_inputs(condition,runs,labels,raw);n=x.shape[1]
        i,j=list(rr).index(a),list(rr).index(b)
        xx=morlet_cwt(correlation(x[i],x[i]),LAM/FOURIER_FACTOR)
        yy=morlet_cwt(correlation(x[j],x[j]),LAM/FOURIER_FACTOR)
        xy=morlet_cwt(correlation(x[i],x[j]),LAM/FOURIER_FACTOR)
        px,py=abs(xx[:,n-1]),abs(yy[:,n-1])
        independent=abs(xy)**2/(px*py)[:,None]
        fast,valid=Wavelet(n).cross(x[i],x[j],px,py)
        error=float(np.max(abs(independent-fast)));assert error<1e-10
        scores=csv_read(ROOT/f'results/scores_{condition}.csv')
        row=scores[((scores.run_x==a)&(scores.run_y==b))|((scores.run_x==b)&(scores.run_y==a))].iloc[0]
        expected_max=float(independent.max())
        expected_mean=float(np.average(independent.mean(axis=1),weights=H))
        w=H*np.sqrt(px*py);w/=w.sum()
        expected_common=float(np.max(np.sum(independent*w[:,None],axis=0)))
        np.testing.assert_allclose([row.global_max,row.map_mean,row.common_lag_peak],[expected_max,expected_mean,expected_common],atol=1e-11)
        errors.append(dict(condition=condition,run_a=a,run_b=b,max_wco_map_error=error))
    # Um descritor por ensaio em cada condição; todos os pares distintos.
    count=0
    for condition in CFG['conditions']:
        scores=csv_read(ROOT/f'results/scores_{condition}.csv')
        with np.load(ROOT/f'results/profiles_{condition}.npz',allow_pickle=False) as z:
            n_runs=len(z['runs']);pairs=z['pairs']
            assert len(pairs)==n_runs*(n_runs-1)//2==len(scores)
            assert len({tuple(p) for p in pairs})==len(pairs)
            assert (pairs[:,0]<pairs[:,1]).all()
            assert np.isfinite(z['mean']).all() and np.isfinite(z['self_area']).all()
            assert (z['self_area']>0).all() and (z['self_width']>0).all()
            assert np.all(z['peak']>=z['zero']-1e-12)
        assert scores.n.nunique()==1
        np.testing.assert_array_less(scores.common_lag_peak.to_numpy()-scores.scale_peak_weighted.to_numpy(),np.full(len(scores),1e-11))
        assert scores.global_max.between(0,1+1e-10).all()
        count+=len(scores)
    assert count==2388
    # Comparação de silhuetas por expressão independente e contagem literal de vizinhos.
    combined=pd.concat([csv_read(ROOT/'results/metrics.csv'),csv_read(ROOT/'results/metrics_addendum.csv')],ignore_index=True)
    with np.load(ROOT/'results/distance_matrices.npz',allow_pickle=False) as z:
        d=z['3171_prefix__global_max']
    mask=np.isin(labels,['churn','wavy annular']);d=d[np.ix_(mask,mask)];ll=labels[mask]
    class_values=[]
    for g in ['churn','wavy annular']:
        values=[]
        for i in np.flatnonzero(ll==g):
            others=np.flatnonzero((ll==g)&(np.arange(len(ll))!=i))
            aa=sum(d[i,j] for j in others)/len(others)
            otherclass=np.flatnonzero(ll!=g);bb=sum(d[i,j] for j in otherclass)/len(otherclass)
            values.append((bb-aa)/max(aa,bb))
        class_values.append(np.mean(values))
    row=combined[(combined.condition=='3171_prefix')&(combined.metric=='global_max')&(combined.scope=='churn_annular')].iloc[0]
    np.testing.assert_allclose(row.macro_silhouette,np.mean(class_values),atol=1e-13)
    np.testing.assert_allclose(row.min_class_silhouette,min(class_values),atol=1e-13)
    neighbor_values=[]
    for i in range(len(ll)):
        choices=[j for j in range(len(ll)) if j!=i]
        nearest=min(d[i,j] for j in choices)
        ties=[j for j in choices if np.isclose(d[i,j],nearest,rtol=1e-12,atol=1e-12)]
        neighbor_values.append(sum(ll[j]==ll[i] for j in ties)/len(ties))
    np.testing.assert_allclose(row.balanced_neighbor,np.mean([np.mean(np.array(neighbor_values)[ll==g]) for g in np.unique(ll)]),atol=1e-13)
    # A inexistência de separação estável é conferida sem depender do limite 0.10.
    stability=combined[~combined.condition.str.startswith('7001')].groupby(['metric','scope']).agg(
        min_class=('min_class_silhouette','min'),min_neighbor=('neighbor_lift','min'),conditions=('condition','nunique'))
    positive_classes=int(((stability.min_class>0)&(stability.conditions==6)).sum())
    stable_count=int(((stability.min_class>0)&(stability.min_neighbor>0)&(stability.conditions==6)).sum())
    # Replicações de coorte, não 186 testes por réplica.
    null=csv_read(ROOT/'results/phase_null_replicates.csv');summary=csv_read(ROOT/'results/phase_null_summary.csv')
    B=json.loads((ROOT/'configs/runtime.json').read_text())['null_replicates']
    assert len(null)==6*B and len(summary)==6
    values=[]
    for row in summary.itertuples():
        f=null[(null.scope==row.scope)&(null.metric==row.metric)]
        assert f.replicate.nunique()==B
        tail=int((f.balanced_contrast>=row.observed_contrast-1e-12).sum())
        expected=(tail+1)/(B+1)
        np.testing.assert_allclose(row.p_model,expected,atol=1e-14);values.append(expected)
    # Holm por definição de máximo sucessivo, independente da função usada.
    order=np.argsort(values);expected=np.zeros(6)
    for rank,index in enumerate(order):expected[index]=min(1,max((6-k)*values[order[k]] for k in range(rank+1)))
    np.testing.assert_allclose(summary.p_model_holm,expected,atol=1e-14)
    noise=csv_read(ROOT/'results/independent_noise_controls.csv')
    assert len(noise)==594 and (noise.groupby(['n','model']).size()==99).all()
    geometry=csv_read(ROOT/'validation/kernel_geometry.csv')
    assert len(geometry)==45
    # A banda longa pode ser indefinida no domínio finito. Rejeitá-la é parte
    # da auditoria: não exigir PSD de matrizes que o adendo mandou descartar.
    expected_valid=(geometry.minimum_eigenvalue>=-1e-10)&(geometry.maximum_above_one<=1e-10)
    assert np.array_equal(expected_valid,geometry.geometry_valid)
    additional=csv_read(ROOT/'results/metrics_addendum.csv')
    with np.load(ROOT/'results/distance_addendum.npz',allow_pickle=False) as matrices:
        for row in geometry.itertuples():
            key=f'{row.condition}__{row.metric}'
            included=((additional.condition==row.condition)&(additional.metric==row.metric)).any()
            assert (key in matrices.files)==row.geometry_valid==included
            if row.geometry_valid:
                d=matrices[key]
                np.testing.assert_allclose(d,d.T,atol=1e-13)
                np.testing.assert_allclose(np.diag(d),0,atol=1e-13)
                assert np.isfinite(d).all() and (d>=0).all()
    manifest=json.loads((ROOT/'results/manifest_started.json').read_text())
    assert sha(ROOT/'configs/wco_v1.json')==manifest['config_sha256']
    assert all(sha(SRC/k)==v for k,v in manifest['sources'].items())
    output=dict(status='PASS',original_engine_map_checks=errors,original_engine_max_error=max(r['max_wco_map_error'] for r in errors),
                pair_maps=count,independent_noise_pairs=len(noise),null_cohort_replicates=B,null_hypotheses=6,
                empirical_kernel_matrices_checked=len(geometry),accepted_kernel_matrices=int(geometry.geometry_valid.sum()),
                rejected_kernel_matrices=geometry[~geometry.geometry_valid][['condition','metric','minimum_eigenvalue']].to_dict('records'),
                all_accepted_kernel_matrices_psd_within_tolerance=True,
                stable_positive_classes_without_010_threshold=positive_classes,
                stable_methods_without_010_threshold=stable_count,source_files_unchanged=len(manifest['sources']),
                scope='Verificação numérica; não valida representatividade, hipóteses físicas ou generalização.')
    json_write(ROOT/'validation/final.json',output)
    print(json.dumps(output,indent=2,ensure_ascii=False),flush=True)


if __name__=='__main__':main()
