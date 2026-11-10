"""Validação independente do recorte/detector e equivalência do cache; benchmark."""
import json
import resource
import time
from datetime import datetime,timezone
import numpy as np
import pandas as pd
from analysis_common import ROOT,LAM,load_signals,evaluate,crop_prepare
from wavelet_core import correlation,prepare
from catalog_texture_anchors import detect,select
from texture_common import OUT,CFG,Evaluator,spec,verify_freeze

def main():
    verify_freeze();checks={};times={};inv,signals=load_signals()
    cat=pd.read_csv(OUT/'catalogo_ancoras.csv').query('selecionada').set_index('anchor_id')
    pairs=pd.read_csv(OUT/'pares_ancoras.csv').sort_values(['run_x','anchor_x','run_y','anchor_y'])
    ev=Evaluator(signals)
    c,_=detect(np.ones(4096)*120);assert c.empty
    t=np.arange(4096);raw=120-50*np.exp(-.5*((t-1500)/45)**2)
    c,_=detect(raw);a=c[(c.tipo=='minimo')&c.elegivel_contexto];assert list(a.ancora)==[1500]
    shifted=120-50*np.exp(-.5*((t-1520)/45)**2)
    c2,_=detect(shifted);assert list(c2[c2.tipo=='minimo'].ancora)==[1520]
    invc,_=detect(240-raw);assert list(invc[invc.tipo=='maximo'].ancora)==[1500]
    assert select(c)==select(c.copy());checks['detector_sinal_constante_centro_translacao_polaridade']='PASS'
    maxerr=0.;scoreerr=0.
    for variant in CFG['variantes']:
        r=pairs.iloc[0];kx=spec(cat.loc[r.anchor_x],variant,'x');ky=spec(cat.loc[r.anchor_y],variant,'y')
        assert kx is not None and ky is not None
        m,s,x,y=ev.pair(kx,ky);ref,rs=evaluate(x,y,engine=kx[-1])
        for key in ['wxx','wyy','wxy','wco']:
            err=float(np.nanmax(abs(m[key]-ref[key])));maxerr=max(maxerr,err);assert err<1e-10
        scoreerr=max(scoreerr,abs(s['escore']-rs['escore']))
        del m,ref
    checks['cache_vs_nucleo_erro_maximo']=maxerr;checks['cache_vs_nucleo_erro_escore']=scoreerr
    r=pairs.iloc[0];kx=spec(cat.loc[r.anchor_x],'principal');ky=spec(cat.loc[r.anchor_y],'principal')
    x=ev.raw_prepared(kx);y=ev.raw_prepared(ky)
    direct=np.correlate(x,y,'full')/len(x)
    err=float(np.max(abs(direct-correlation(x,y))));assert err<1e-12
    checks['correlacao_direta_erro_maximo']=err
    # Pré-processamento independente: regressão linear por mínimos quadrados.
    a=255-signals[kx[0]][kx[2]][kx[1]:kx[1]+2048]
    design=np.column_stack([np.arange(2048),np.ones(2048)])
    z=a-design@np.linalg.lstsq(design,a,rcond=None)[0];z-=z.mean();z/=z.std()
    err=float(np.max(abs(x-z)));assert err<1e-10;checks['recorte_antes_regressao_erro_maximo']=err
    m,s=evaluate(x,y);neg,_=evaluate(x,-y)
    err=float(np.nanmax(abs(m['wco']-neg['wco'])));assert err<1e-12
    checks['invariancia_WCO_inversao_sinal_erro_maximo']=err
    checks['fase_WCS_inversao_sinal_erro_maximo']=float(np.max(abs(m['wxy']+neg['wxy'])))
    del m,neg
    impulse=np.zeros(2048);impulse[800]=1;delayed=np.zeros(2048);delayed[900]=1
    lag=int(np.argmax(correlation(impulse,delayed))-2047);assert lag==-100
    checks['atraso_y_100px_lag_pico']=lag
    selfs,_=evaluate(x,x);err=float(np.max(abs(selfs['wco'][selfs['denominator_valid'],2047]-1)))
    assert err<1e-12;checks['self_zero_erro_maximo']=err;del selfs
    # Rechecagem de uma referência fixa reutilizada: mesmo domínio, núcleo e entrada.
    old=pd.read_csv(ROOT/'results/resultados_fixos.csv');runs=set(cat.aquisicao)
    ref=old[(old.variante=='principal')&old.run_x.isin(runs)&old.run_y.isin(runs)].iloc[0]
    seg=pd.read_csv(ROOT/'results/segmentos_fixos.csv').set_index('segment_id')
    ix,iy=ref.pair_id.split('__')
    m,s=evaluate(crop_prepare(signals[ref.run_x]['principal'],int(seg.loc[ix].inicio)),crop_prepare(signals[ref.run_y]['principal'],int(seg.loc[iy].inicio)))
    err=abs(s['escore']-ref.escore);assert err<1e-12;checks['referencia_reutilizada_erro_escore']=err;del m
    q=np.arange(2048)-1024
    valley=lambda width:50*np.exp(-.5*(q/width)**2)
    background=25*np.sin(2*np.pi*q/180)+20*np.sin(2*np.pi*q/420)
    synth=[]
    for name,rx,ry in [
        ('minimos_compativeis',120-valley(40),120-valley(48)),
        ('minimos_larguras_distintas',120-valley(40),120-valley(120)),
        ('fundo_comum_larguras_distintas',120+background-valley(40),120+background-valley(120)),
        ('polaridade_invertida',120-valley(40),120+valley(40))]:
        m,s=evaluate(prepare(rx),prepare(ry));synth.append(dict(caso=name,**s));del m
    pd.DataFrame(synth).to_csv(OUT/'controles_sinteticos.csv',index=False)
    # Custos representativos, cache inicialmente vazio; oito pares por variante de motor.
    for variant in ['principal','v5_mesmos_segmentos']:
        ev.cache.clear();t0=time.perf_counter()
        for _,r in pairs.head(12).iterrows():
            m,_,_,_=ev.pair(spec(cat.loc[r.anchor_x],variant,'x'),spec(cat.loc[r.anchor_y],variant,'y'));del m
        times[variant]=(time.perf_counter()-t0)/12
    planned=len(pairs)*len(CFG['variantes'])+1539+len(cat)+19
    naive=(len(pairs)*7+1026+len(cat)+19)*times['principal']+(len(pairs)+513)*times['v5_mesmos_segmentos']
    bench=dict(data_utc=datetime.now(timezone.utc).isoformat(),segundos_por_par=times,comparacoes_planejadas_novas=planned,
        referencias_reutilizadas=513,estimativa_segundos=naive,intervalo_planejado_minutos=[naive/60*.8,naive/60*1.8],
        pico_memoria_benchmark_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,
        nota='Estimativa de tempo, não garantia. Ausências de correspondência óptica reduzem o número executável.')
    (OUT/'benchmark.json').write_text(json.dumps(bench,indent=2,ensure_ascii=False)+'\n')
    report=dict(status='PASS',data_utc=datetime.now(timezone.utc).isoformat(),checks=checks)
    (OUT/'validacao_numerica.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(report,indent=2,ensure_ascii=False));print(json.dumps(bench,indent=2,ensure_ascii=False))

if __name__=='__main__':main()
