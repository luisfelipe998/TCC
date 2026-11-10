"""Extensão churn/anular: catálogo congelado, oito variantes e controles."""
import argparse
import itertools
import json
import resource
import time
from datetime import datetime,timezone
import numpy as np
import pandas as pd
from analysis_common import ROOT,load_signals,save_map,LAM
from texture_common import OUT,MAPS,CFG,Evaluator,spec,verify_freeze,sha

def anchors(signals):
    t0=time.perf_counter();ev=Evaluator(signals)
    cat=pd.read_csv(OUT/'catalogo_ancoras.csv').query('selecionada').set_index('anchor_id')
    pairs=pd.read_csv(OUT/'pares_ancoras.csv').sort_values(['run_x','anchor_x','run_y','anchor_y'])
    examples=set(pairs.sort_values('pair_id').groupby(['pattern_x','pattern_y','grupo']).first().pair_id)
    rows=[];profiles=[];profile_ids=[]
    for variant in CFG['variantes']:
        start=time.perf_counter();ev.cache.clear()
        for j,r in enumerate(pairs.itertuples(index=False)):
            sx=spec(cat.loc[r.anchor_x],variant,'x');sy=spec(cat.loc[r.anchor_y],variant,'y')
            row=r._asdict();row.update(variante=variant,status='calculado',motivo='')
            if sx is None or sy is None:
                row.update(status='nao_estimavel',motivo='sem_correspondencia_optica_em_64px')
            else:
                m,s,x,y=ev.pair(sx,sy)
                row.update(centro_x=sx[1]+1024,centro_y=sy[1]+1024,**s)
                if variant=='principal':
                    profiles.append(m['wco'][:,abs(m['lags'])<=32].mean(axis=1));profile_ids.append(r.pair_id)
                    if r.pair_id in examples:save_map(MAPS/f'{r.pair_id}.npz',m,x,y,**row,escolha='primeiro_ID_antes_dos_resultados')
                del m
            rows.append(row)
            if (j+1)%100==0:print(f'Âncoras {variant}: {j+1}/{len(pairs)}',flush=True)
        pd.DataFrame(rows).to_csv(OUT/'resultados_ancoras.csv',index=False)
        print(f'Âncoras {variant}: {time.perf_counter()-start:.1f} s',flush=True)
    np.savez_compressed(OUT/'perfis_ancoras.npz',pair_ids=np.array(profile_ids,str),wavelengths=LAM,mean_lag=np.stack(profiles))
    return dict(segundos=time.perf_counter()-t0,linhas=len(rows),calculados=sum(r['status']=='calculado' for r in rows),
                cache_acertos=ev.hits,cache_faltas=ev.misses,mapas_predeterminados=len(examples))

def fixed(inv,signals):
    t0=time.perf_counter();ev=Evaluator(signals)
    runs=sorted(inv[inv.pattern.isin(CFG['padroes'])].run);patterns=inv.set_index('run').pattern.to_dict()
    original=pd.read_csv(ROOT/'results/resultados_fixos.csv')
    reused=original[(original.variante=='principal')&original.run_x.isin(runs)&original.run_y.isin(runs)].copy()
    assert len(reused)==513
    reused['origem']='etapa1_reutilizada_sem_alteracao';reused['status']='calculado'
    rows=reused.to_dict('records')
    segments=pd.read_csv(ROOT/'results/segmentos_fixos.csv').set_index('segment_id')
    for variant in ['faixa30_49','somente_media','v5_mesmos_segmentos']:
        count=0;ev.cache.clear();start=time.perf_counter()
        for rx,ry in itertools.combinations(runs,2):
            for k in range(3):
                a=segments.loc[f'{rx}_F{k}'];b=segments.loc[f'{ry}_F{k}']
                band='faixa30_49' if variant=='faixa30_49' else 'principal'
                linear=variant!='somente_media';engine='legacy' if variant=='v5_mesmos_segmentos' else 'morlet'
                m,s,_,_=ev.pair((rx,int(a.inicio),band,linear,engine),(ry,int(b.inicio),band,linear,engine));del m
                rows.append(dict(pair_id=f'{rx}_F{k}__{ry}_F{k}',run_x=rx,run_y=ry,pattern_x=patterns[rx],pattern_y=patterns[ry],
                    grupo='intra_fixo' if patterns[rx]==patterns[ry] else 'entre_fixo',posicao=a.posicao,variante=variant,
                    status='calculado',origem='extensao_calculada',**s));count+=1
                if count%150==0:print(f'Fixos {variant}: {count}/513',flush=True)
        pd.DataFrame(rows).to_csv(OUT/'resultados_fixos.csv',index=False)
        print(f'Fixos {variant}: {time.perf_counter()-start:.1f} s',flush=True)
    return dict(segundos=time.perf_counter()-t0,comparacoes_novas=len(rows)-len(reused),reutilizadas=len(reused))

def diagnostics(signals):
    t0=time.perf_counter();ev=Evaluator(signals);rows=[]
    cat=pd.read_csv(OUT/'catalogo_ancoras.csv').query('selecionada')
    for run,g in cat.groupby('aquisicao'):
        for _,r in g.iterrows():
            k=spec(r,'principal');m,s,_,_=ev.pair(k,k)
            err=float(np.max(abs(m['wco'][m['denominator_valid'],2047]-1)))
            assert err<1e-12;del m
            rows.append(dict(aquisicao=run,padrao=r.padrao,anchor_x=r.anchor_id,anchor_y=r.anchor_id,tipo='self',erro_self_zero=err,**s))
        h=g[g.tipo=='minimo'].sort_values('ancora');a,b=h.iloc[0],h.iloc[1]
        m,s,_,_=ev.pair(spec(a,'principal'),spec(b,'principal'));del m
        rows.append(dict(aquisicao=run,padrao=a.padrao,anchor_x=a.anchor_id,anchor_y=b.anchor_id,tipo='mesma_aquisicao',erro_self_zero=np.nan,**s))
    pd.DataFrame(rows).to_csv(OUT/'controles_mesma_aquisicao.csv',index=False)
    return dict(segundos=time.perf_counter()-t0,comparacoes=len(rows))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['all','anchors','fixed','diagnostics'],default='all');args=parser.parse_args()
    verify_freeze();assert json.loads((OUT/'validacao_numerica.json').read_text())['status']=='PASS'
    inv,signals=load_signals();MAPS.mkdir(exist_ok=True)
    t0=time.perf_counter();meta=dict(inicio_utc=datetime.now(timezone.utc).isoformat(),etapas={},status='EM_EXECUCAO',
        scripts_sha256={p:sha(ROOT/'scripts'/p) for p in ['texture_common.py','run_texture_analysis.py','wavelet_core.py']})
    for name,fun in [('anchors',lambda:anchors(signals)),('fixed',lambda:fixed(inv,signals)),('diagnostics',lambda:diagnostics(signals))]:
        if args.stage in ['all',name]:
            meta['etapas'][name]=fun();(OUT/'execucao.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False)+'\n')
    verify_freeze();meta.update(status='PASS',fim_utc=datetime.now(timezone.utc).isoformat(),segundos_total=time.perf_counter()-t0,
        pico_memoria_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2)
    (OUT/'execucao.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False)+'\n');print(json.dumps(meta,indent=2),flush=True)

if __name__=='__main__':main()
