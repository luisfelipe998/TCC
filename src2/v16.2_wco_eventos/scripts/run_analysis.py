"""Executa somente a análise nova com catálogo congelado e K=1."""
import argparse
import hashlib
import itertools
import json
import resource
import time
from datetime import datetime,timezone
import numpy as np
import pandas as pd
from analysis_common import ROOT,CFG,LAM,load_signals,crop_prepare,evaluate,save_map

VARIANTS={
 'principal':('principal',True,'morlet',0,0),
 'v5_mesmos_segmentos':('principal',True,'legacy',0,0),
 'faixa30_49':('faixa30_49',True,'morlet',0,0),
 'somente_media':('principal',False,'morlet',0,0),
 'bordas_menos16':('principal',True,'morlet',-16,-16),
 'bordas_mais16':('principal',True,'morlet',16,16),
 'ancoras_opostas_a':('principal',True,'morlet',-16,16),
 'ancoras_opostas_b':('principal',True,'morlet',16,-16)}

def verify_freeze():
    f=json.loads((ROOT/'results/congelamento.json').read_text())
    for path,digest in f['sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,path
    assert json.loads((ROOT/'results/validacao_numerica.json').read_text())['status']=='PASS'

def checked_eval(x,y,engine):
    m,s=evaluate(x,y,engine=engine)
    assert s['erro_equacao']<1e-10
    if engine=='morlet':assert s['maximo_integral']<=1+1e-8
    return m,s

def events(signals):
    t0=time.perf_counter();out=ROOT/'results';maps=ROOT/'maps';maps.mkdir(exist_ok=True)
    cat=pd.read_csv(out/'catalogo_eventos.csv').set_index('event_id')
    pairs=pd.read_csv(out/'pares_eventos.csv').sort_values('pair_id')
    examples=set(pairs.groupby('grupo',sort=True).first().pair_id)
    results=[];profiles=[]
    for variant,(band,linear,engine,shiftx,shifty) in VARIANTS.items():
        start=time.perf_counter()
        for j,r in enumerate(pairs.itertuples(index=False)):
            cx=cat.loc[r.event_x];cy=cat.loc[r.event_y]
            x=crop_prepare(signals[r.run_x][band],int(cx.segmento_inicio)+shiftx,linear=linear)
            y=crop_prepare(signals[r.run_y][band],int(cy.segmento_inicio)+shifty,linear=linear)
            m,s=checked_eval(x,y,engine)
            row={k:getattr(r,k) for k in ['pair_id','event_x','event_y','run_x','run_y','pattern_x','pattern_y','grupo','razao_extensoes']}
            row.update(variante=variant,**s);results.append(row)
            if variant=='principal':
                domain=np.abs(m['lags'])<=CFG['lag_resumo_px']
                profiles.append(m['wco'][:,domain].mean(axis=1))
                if r.pair_id in examples:
                    save_map(maps/f'{r.pair_id}.npz',m,x,y,**row)
            del m
            if (j+1)%50==0:print(f'Eventos {variant}: {j+1}/{len(pairs)}',flush=True)
        print(f'Eventos {variant}: concluído em {time.perf_counter()-start:.1f} s',flush=True)
        pd.DataFrame(results).to_csv(out/'resultados_eventos.csv',index=False)
    np.savez_compressed(out/'perfis_eventos.npz',pair_ids=pairs.pair_id.to_numpy(str),wavelengths=LAM,mean_lag=np.stack(profiles))
    return dict(segundos=time.perf_counter()-t0,comparacoes=len(results),mapas_predeterminados=len(examples))

def fixed(inv,signals):
    t0=time.perf_counter();out=ROOT/'results';segments=pd.read_csv(out/'segmentos_fixos.csv').set_index('segment_id')
    obs=set(pd.read_csv(out/'catalogo_eventos.csv').query('incluida').aquisicao)
    rows=[];runs=list(inv.run);patterns=inv.set_index('run').pattern.to_dict()
    variants=['principal','v5_mesmos_segmentos','faixa30_49','somente_media']
    for variant in variants:
        band,linear,engine,_,_=VARIANTS[variant];count=0
        for rx,ry in itertools.combinations(runs,2):
            if variant!='principal' and not {rx,ry}<=obs:continue
            for k in range(3):
                sx=segments.loc[f'{rx}_F{k}'];sy=segments.loc[f'{ry}_F{k}']
                x=crop_prepare(signals[rx][band],int(sx.inicio),linear=linear)
                y=crop_prepare(signals[ry][band],int(sy.inicio),linear=linear)
                m,s=checked_eval(x,y,engine);del m
                rows.append(dict(pair_id=f'{rx}_F{k}__{ry}_F{k}',run_x=rx,run_y=ry,
                    pattern_x=patterns[rx],pattern_y=patterns[ry],grupo='intra_fixo' if patterns[rx]==patterns[ry] else 'entre_fixo',
                    posicao=sx.posicao,variante=variant,**s));count+=1
                if count%150==0:print(f'Referências {variant}: {count}',flush=True)
        print(f'Referências {variant}: {count} concluídas',flush=True)
        pd.DataFrame(rows).to_csv(out/'resultados_fixos.csv',index=False)
    return dict(segundos=time.perf_counter()-t0,comparacoes=len(rows))

def diagnostics(signals):
    t0=time.perf_counter();out=ROOT/'results';cat=pd.read_csv(out/'catalogo_eventos.csv').query('incluida')
    rows=[]
    # Todos os selfs; primeiro par de eventos em ordem de coordenadas por aquisição.
    for run,g in cat.groupby('aquisicao'):
        g=g.sort_values('inicio')
        for r in g.itertuples(index=False):
            x=crop_prepare(signals[run]['principal'],int(r.segmento_inicio));m,s=checked_eval(x,x,'morlet')
            err=float(np.nanmax(np.abs(m['wco'][:,2047]-1)));assert err<1e-12
            rows.append(dict(aquisicao=run,event_x=r.event_id,event_y=r.event_id,tipo='self',erro_self_zero=err,**s))
            del m
        if len(g)>=2:
            a,b=g.iloc[0],g.iloc[1]
            x=crop_prepare(signals[run]['principal'],int(a.segmento_inicio));y=crop_prepare(signals[run]['principal'],int(b.segmento_inicio))
            m,s=checked_eval(x,y,'morlet');del m
            rows.append(dict(aquisicao=run,event_x=a.event_id,event_y=b.event_id,tipo='mesma_aquisicao',erro_self_zero=np.nan,**s))
    pd.DataFrame(rows).to_csv(out/'controles_mesma_aquisicao.csv',index=False)
    return dict(segundos=time.perf_counter()-t0,comparacoes=len(rows))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['all','events','fixed','diagnostics'],default='all');args=parser.parse_args()
    verify_freeze();inv,signals=load_signals();t0=time.perf_counter()
    meta=dict(inicio_utc=datetime.now(timezone.utc).isoformat(),etapas={})
    for name,fun in [('events',lambda:events(signals)),('fixed',lambda:fixed(inv,signals)),('diagnostics',lambda:diagnostics(signals))]:
        if args.stage in ['all',name]:
            meta['etapas'][name]=fun()
            (ROOT/'results/execucao.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2)+'\n')
    verify_freeze();meta.update(fim_utc=datetime.now(timezone.utc).isoformat(),segundos_total=time.perf_counter()-t0,
        pico_memoria_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,status='PASS')
    (ROOT/'results/execucao.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(meta,ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':main()
