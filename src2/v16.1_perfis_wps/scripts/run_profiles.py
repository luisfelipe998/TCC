"""Executa o protocolo fixado sem alterar as análises de origem."""
from itertools import combinations
from datetime import datetime, timezone
import hashlib
import json
import time
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from PIL import Image
from profile_core import (ROOT, BASE, LAM, wc, quadrature, distribution,
                          hellinger, cosine_distance, zero_profile, legacy_zero, support_counts)
from evaluation import summarize, summarize_pattern

PT = {'dispersed bubbles':'Bolhas dispersas', 'plug':'Bolhas alongadas (plug)',
      'slug':'Golfadas (slug)', 'churn':'Agitado (churn)', 'wavy annular':'Anular ondulado'}


def main():
    start = time.perf_counter()
    out = ROOT/'results'
    cfg = json.loads((ROOT/'config.json').read_text())
    assert hashlib.sha256((ROOT/'PROTOCOL.md').read_bytes()).hexdigest() == cfg['protocol_sha256']
    assert json.loads((out/'operations_validation.json').read_text())['status'] == 'PASS'
    inv = pd.read_csv(BASE/'results/input_audit.csv')
    inv = inv[inv.included].sort_values('run').copy()
    assert sorted(inv.run) == cfg['expected_included_runs']
    runs = inv.run.tolist()
    labels = inv.set_index('run').pattern.map(PT).to_dict()
    lengths = inv.set_index('run').n.to_dict()
    common = cfg['common_length']
    domains = {'principal':np.array(cfg['primary_mask']),
               'integral':np.ones(len(LAM),bool),
               'conservador':np.array(cfg['conservative_mask'])}
    source_profiles = pd.read_csv(BASE/'results/profiles.csv')
    source_windows = pd.read_csv(BASE/'results/equal_length_windows.csv')
    raw = {r:np.load(BASE/'intermediate'/f'{r}_input.npz')['raw'] for r in runs}
    prepared = {r:wc.prepare(raw[r]) for r in runs}
    profiles, rows, lagrows = {}, [], []
    errors = []

    def put(run, variant, values, n, origin):
        values = np.asarray(values, float)
        assert values.shape == LAM.shape
        profiles[run,variant] = values
        for la,val in zip(LAM,values):
            rows.append(dict(aquisicao=run,padrao=labels[run],variante=variant,n=n,
                             comprimento_onda_px=la,perfil=val,origem=origin))

    for i, run in enumerate(runs):
        n = lengths[run]
        for method, variant in [('v16_biased','principal'),('v5','v5_nominal')]:
            g = source_profiles[(source_profiles.run==run)&(source_profiles.method==method)].sort_values('wavelength_px')
            assert len(g)==200 and np.allclose(g.wavelength_px,LAM,rtol=1e-14)
            put(run,variant,g.zero_lag_wps.to_numpy(),n,'Perfis herdados, hashes verificados')
            if method == 'v16_biased':
                put(run,'media_lags_todos',g.lag_mean_all.to_numpy(),n,'Média integral herdada')
        w = wc.morlet_cwt(wc.correlation(prepared[run],prepared[run]),LAM/wc.FOURIER_FACTOR)
        mag = np.abs(w)
        err = float(np.max(np.abs(mag[:,n-1]-profiles[run,'principal'])))
        meanerr = float(np.max(np.abs(mag.mean(axis=1)-profiles[run,'media_lags_todos'])))
        assert np.allclose(mag[:,n-1],profiles[run,'principal'],atol=1e-12,rtol=1e-10)
        assert np.allclose(mag.mean(axis=1),profiles[run,'media_lags_todos'],atol=1e-12,rtol=1e-10)
        errors.append(dict(aquisicao=run,erro_zero=err,erro_media_integral=meanerr))
        coi, cons = support_counts(n)
        for mask, name in [(coi,'media_lags_coi'),(cons,'media_lags_conservador')]:
            count = mask.sum(axis=1)
            values = np.divide(np.sum(np.where(mask,mag,0),axis=1),count,
                               out=np.full(len(LAM),np.nan),where=count>0)
            put(run,name,values,n,'Média recalculada integralmente com suporte explícito')
        for j,la in enumerate(LAM):
            lagrows.append(dict(aquisicao=run,padrao=labels[run],comprimento_onda_px=la,
                lags_todos=2*n-1,lags_coi=int(coi[j].sum()),lags_conservador=int(cons[j].sum()),
                zero_dentro_coi=bool(coi[j,n-1]),zero_dentro_conservador=bool(cons[j,n-1])))
        del w,mag,coi,cons
        with Image.open(ROOT.parent/'data_processed'/run/'__panorama_v11.png') as im:
            pano=np.asarray(im.convert('L'))
        for var,signal,linear in [('faixa30_49',pano[30:50].mean(axis=0),True),
                                  ('faixa20_59',pano[20:60].mean(axis=0),True),
                                  ('somente_media',raw[run],False)]:
            put(run,var,zero_profile(wc.prepare(signal,linear=linear)),n,'Perfil novo; extração e preparo declarados')
        put(run,'v5_mesmas_escalas',legacy_zero(prepared[run],same_scales=True),n,
            'Mesmo sinal preparado e mesmos valores de a do v16; discretização PyWavelets')
        for window, variant in [('prefix','inicio3171'),('center','centro3171'),('suffix','fim3171')]:
            g=source_windows[(source_windows.run==run)&(source_windows.window==window)].sort_values('wavelength_px')
            assert len(g)==200 and set(g.n)=={common}
            put(run,variant,g.zero_lag_wps.to_numpy(),common,'Segmentos herdados, corte antes do preparo')
        # Auditoria independente dos três segmentos para menor/mediano/maior comprimento.
        if n in sorted(lengths.values())[::12] or n==max(lengths.values()):
            for variant,offset in [('inicio3171',0),('centro3171',(n-common)//2),('fim3171',n-common)]:
                actual=zero_profile(wc.prepare(raw[run][offset:offset+common]))
                assert np.allclose(actual,profiles[run,variant],rtol=1e-10,atol=1e-12)
        print(f'Perfis por aquisição {i+1}/26: {run}; {time.perf_counter()-start:.1f} s',flush=True)

    pd.DataFrame(rows).to_csv(out/'profiles_raw.csv',index=False)
    pd.DataFrame(lagrows).to_csv(out/'lag_support_counts.csv',index=False)
    pd.DataFrame(errors).to_csv(out/'reused_profiles_verification.csv',index=False)
    normalized=[]
    for (run,var),p in profiles.items():
        for domain,mask in domains.items():
            if not np.isfinite(p[mask]).all():
                continue
            axis=LAM[mask];w=quadrature(axis);mass=distribution(p[mask],axis)
            for la,v,ww,mm in zip(axis,p[mask],w,mass):
                normalized.append(dict(aquisicao=run,padrao=labels[run],variante=var,dominio=domain,
                    comprimento_onda_px=la,perfil=v,peso_quadratura=ww,massa=mm,
                    area=float(w@p[mask]),densidade_log=v/(w@p[mask])))
    pd.DataFrame(normalized).to_csv(out/'profiles_normalized.csv',index=False)

    specs = {
        'wps_principal':('principal','principal','hellinger','WPS, registro inteiro'),
        'wps_cosseno':('principal','principal','cosine','WPS, distância cosseno'),
        'wps_integral':('principal','integral','hellinger','WPS, faixa integral'),
        'wps_conservador':('principal','conservador','hellinger','WPS, faixa conservadora'),
        'media_lags_todos':('media_lags_todos','principal','hellinger','Média em todos os lags'),
        'media_lags_coi':('media_lags_coi','principal','hellinger','Média em lags sob COI'),
        'media_lags_conservador':('media_lags_conservador','conservador','hellinger','Média sob suporte conservador'),
        'faixa30_49':('faixa30_49','principal','hellinger','WPS, linhas 30–49'),
        'faixa20_59':('faixa20_59','principal','hellinger','WPS, linhas 20–59'),
        'somente_media':('somente_media','principal','hellinger','WPS, somente remoção da média'),
        'v5_nominal':('v5_nominal','principal','hellinger','Perfil v5, eixo nominal'),
        'v5_mesmas_escalas':('v5_mesmas_escalas','principal','hellinger','Perfil v5, mesmas escalas a'),
        'inicio3171':('inicio3171','principal','hellinger','Segmento inicial, 3.171 px'),
        'centro3171':('centro3171','principal','hellinger','Segmento central, 3.171 px'),
        'fim3171':('fim3171','principal','hellinger','Segmento final, 3.171 px')}
    tasks=list(combinations(runs,2))
    distances=[]

    def record(a,b,method,value,direction=-1):
        distances.append(dict(aquisicao_x=a,aquisicao_y=b,padrao_x=labels[a],padrao_y=labels[b],
            mesmo_padrao=labels[a]==labels[b],grupo=labels[a] if labels[a]==labels[b] else ' / '.join(sorted([labels[a],labels[b]])),
            metodo=method,valor=float(value),similaridade=float(direction*value),n_x=lengths[a],n_y=lengths[b]))

    for method,(variant,domain,metric,label) in specs.items():
        mask=domains[domain];axis=LAM[mask]
        masses={r:distribution(profiles[r,variant][mask],axis) for r in runs}
        for a,b in tasks:
            d=hellinger(masses[a],masses[b]) if metric=='hellinger' else cosine_distance(profiles[a,variant][mask],profiles[b,variant][mask],axis)
            record(a,b,method,d)
    # Prefixos específicos de cada par: representação e WCO com o mesmo segmento/domínio nominal.
    cache={};prefix_rows=[]
    def get(run,n):
        if (run,n) not in cache:
            if n==lengths[run]:
                corr=profiles[run,'principal'];old=profiles[run,'v5_nominal']
            else:
                corr=zero_profile(wc.prepare(raw[run][:n]))
                old=legacy_zero(prepared[run][:n])
            cache[run,n]=(distribution(corr,LAM),distribution(old,LAM))
            for j,la in enumerate(LAM):
                prefix_rows.append(dict(aquisicao=run,n=n,comprimento_onda_px=la,
                    perfil_v16=corr[j],perfil_v5_literal=old[j]))
        return cache[run,n]
    for i,(a,b) in enumerate(tasks):
        n=min(lengths[a],lengths[b]);px,ox=get(a,n);py,oy=get(b,n)
        record(a,b,'wps_prefixo_par',hellinger(px,py))
        record(a,b,'v5_prefixo_par',hellinger(ox,oy))
        if i%25==0:
            print(f'Prefixos por par {i+1}/325; {len(cache)} perfis; {time.perf_counter()-start:.1f} s',flush=True)
    pd.DataFrame(prefix_rows).to_csv(out/'pair_prefix_profiles.csv',index=False)

    score='weighted_scale_peak_power05'
    source_scores=pd.read_csv(BASE/'results/preprocessing_comparison_pairs.csv')
    old_scores=pd.read_csv(ROOT.parent/'v16_spatial2/results/preprocessing_comparison_pairs.csv')
    for source,variant,method in [(source_scores,'baseline','wco_v16_18_71'),
            (source_scores,'v5','wco_v5_18_71'),(old_scores,'v5','wco_v5_20_59'),
            (old_scores,'baseline','wco_v16_20_59')]:
        g=source[source.variant==variant]
        assert len(g)==325
        for r in g.itertuples():
            assert PT[r.pattern_x]==labels[r.run_x] and PT[r.pattern_y]==labels[r.run_y]
            record(r.run_x,r.run_y,method,getattr(r,score),1)
    data=pd.DataFrame(distances)
    assert data.groupby('metodo').size().eq(325).all()
    assert data.groupby('metodo').mesmo_padrao.sum().eq(87).all()
    assert np.isfinite(data[['valor','similaridade']]).all().all()
    data.to_csv(out/'all_pairs.csv',index=False)
    summaries=[];patterns=[];influence=[];groups=[]
    for method,g in data.groupby('metodo',sort=False):
        summaries.append(dict(metodo=method,**summarize(g)))
        for pattern in PT.values():
            patterns.append(dict(metodo=method,padrao=pattern,**summarize_pattern(g,pattern)))
        for run in runs:
            sub=g[(g.aquisicao_x!=run)&(g.aquisicao_y!=run)]
            influence.append(dict(metodo=method,aquisicao_retirada=run,**summarize(sub)))
        for group,sub in g.groupby('grupo'):
            groups.append(dict(metodo=method,grupo=group,mesmo_padrao=bool(sub.mesmo_padrao.iloc[0]),
                pares=len(sub),media=float(sub.valor.mean()),mediana=float(sub.valor.median()),
                minimo=float(sub.valor.min()),maximo=float(sub.valor.max())))
    pd.DataFrame(summaries).to_csv(out/'ranking_summary.csv',index=False)
    pd.DataFrame(patterns).to_csv(out/'ranking_by_pattern.csv',index=False)
    pd.DataFrame(influence).to_csv(out/'acquisition_deletion.csv',index=False)
    pd.DataFrame(groups).to_csv(out/'pair_group_summary.csv',index=False)
    sensitivity=[]
    mask=domains['principal'];axis=LAM[mask]
    for r in runs:
        ref=distribution(profiles[r,'principal'][mask],axis)
        for variant in ['inicio3171','centro3171','fim3171','faixa30_49','faixa20_59','somente_media','v5_mesmas_escalas']:
            sensitivity.append(dict(aquisicao=r,padrao=labels[r],variante=variant,n=lengths[r],
                H_mesma_aquisicao=hellinger(ref,distribution(profiles[r,variant][mask],axis))))
    pd.DataFrame(sensitivity).to_csv(out/'within_acquisition_sensitivity.csv',index=False)
    mainpairs=data[data.metodo=='wps_principal'].copy()
    mainpairs['razao_comprimentos']=np.maximum(mainpairs.n_x,mainpairs.n_y)/np.minimum(mainpairs.n_x,mainpairs.n_y)
    rho=float(spearmanr(mainpairs.valor, np.log(mainpairs.razao_comprimentos)).statistic)
    metadata=dict(status='PASS',completed_utc=datetime.now(timezone.utc).isoformat(),
        elapsed_seconds=time.perf_counter()-start,acquisitions=26,pairs_per_method=325,
        methods=data.metodo.nunique(),unique_pair_prefixes=len(cache),
        protocol_sha256=cfg['protocol_sha256'],configuration_sha256=hashlib.sha256((ROOT/'config.json').read_bytes()).hexdigest(),
        max_reused_zero_error=max(r['erro_zero'] for r in errors),
        max_reused_mean_error=max(r['erro_media_integral'] for r in errors),
        H_vs_log_length_ratio_spearman_descriptive=rho,
        method_specs={k:dict(profile=v[0],domain=v[1],metric=v[2],label=v[3]) for k,v in specs.items()},
        inference='Sem valor-p; pares dependentes, retirada de aquisições é análise de influência.')
    (out/'run_metadata.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2)+'\n')
    print(pd.DataFrame(summaries)[['metodo','A_pares','A_balanceada']].to_string(index=False))
    print(f'Concluído em {metadata["elapsed_seconds"]:.1f} s',flush=True)


if __name__ == '__main__':
    main()
