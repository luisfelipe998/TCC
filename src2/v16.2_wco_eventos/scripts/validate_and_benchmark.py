"""Controles independentes e orçamento antes do lote experimental."""
import ast
import hashlib
import json
import time
from datetime import datetime,timezone
import numpy as np
import pandas as pd
import pywt
from scipy.signal import fftconvolve
from wavelet_core import correlation,morlet_cwt,prepare,pair,LEGACY_CF,LEGACY_WAVELET,safe_coherence
from analysis_common import ROOT,BASE,CFG,LAM,evaluate,crop_prepare

def main():
    t0=time.perf_counter();out=ROOT/'results';rng=np.random.default_rng(CFG['semente']);checks={}
    # Congelamento material de regras, revisão e parceiros anterior à WCO.
    files=['config.json','PROTOCOL.md','results/catalogo_eventos.csv','results/revisao_visual.csv',
           'results/pares_candidatos.csv','results/pares_eventos.csv','results/segmentos_fixos.csv',
           'scripts/catalog_events.py','scripts/wavelet_core.py']
    frozen={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files}
    freeze=out/'congelamento.json'
    if freeze.exists():
        old=json.loads(freeze.read_text());assert frozen==old['sha256'],'Regra ou catálogo mudou após congelamento.'
    else:freeze.write_text(json.dumps(dict(data_utc=datetime.now(timezone.utc).isoformat(),sha256=frozen),indent=2)+'\n')
    assert (ROOT/'scripts/wavelet_core.py').read_bytes()==(BASE/'scripts/wavelet_core.py').read_bytes()
    x,y=rng.normal(size=(2,128))
    expected=[sum(x[k+m]*y[k] for k in range(128) if 0<=k+m<128)/128 for m in range(-127,128)]
    checks['erro_correlacao_soma_direta']=float(np.max(np.abs(correlation(x,y)-expected)))
    assert checks['erro_correlacao_soma_direta']<1e-12
    err=[];sc=np.array([20.,50.,100.]);x=rng.normal(size=2048);w=morlet_cwt(x,sc)
    for j,a in enumerate(sc):
        t=np.arange(-int(np.ceil(8*a)),int(np.ceil(8*a))+1)
        k=np.pi**(-.25)/np.sqrt(a)*np.exp(-.5*(t/a)**2)*np.exp(1j*6*t/a)
        err.append(np.max(np.abs(w[j]-fftconvolve(x,k,mode='same'))))
    checks['erro_convolucao_independente']=float(max(err));assert max(err)<1e-10
    checks['erro_padding_8_12']=float(np.max(np.abs(w-morlet_cwt(x,sc,12))))
    assert checks['erro_padding_8_12']<1e-10
    m=pair(x,x,LAM);checks['erro_self_zero']=float(np.max(np.abs(m['wco'][:,2047]-1)))
    assert checks['erro_self_zero']<1e-12 and np.nanmax(m['wco'])<=1+1e-10
    t=np.arange(2048);event=np.exp(-.5*((t-800)/80)**2)*np.cos(2*np.pi*(t-800)/100)
    delayed=np.r_[np.zeros(100),event[:-100]];m=pair(event,delayed,LAM)
    j=int(np.argmax(np.abs(m['wxy']).max(axis=1)))
    lag=int(m['lags'][np.argmax(np.abs(m['wxy'][j]))]);checks['lag_atraso_100']=lag;assert lag==-100
    for bad in [np.ones(128),np.full(128,np.nan),np.arange(128.)]:
        try:prepare(bad)
        except ValueError:pass
        else:raise AssertionError('Sinal degenerado não rejeitado')
    checks['constante_nao_finito_rampa_rejeitados']=True
    co,valid=safe_coherence(np.ones((3,5),complex),np.zeros((3,5),complex),np.ones((3,5),complex),2)
    assert not valid.any() and np.isnan(co).all();checks['denominador_zero_nao_virou_zero']=True
    try:correlation(np.ones(8),np.ones(9))
    except ValueError:pass
    else:raise AssertionError('Comprimentos distintos não rejeitados')
    raw=120+.03*np.arange(4096)+rng.normal(size=4096);seg=raw[300:2348]
    xx=255-seg;xx=xx-np.polyval(np.polyfit(np.arange(len(xx)),xx,1),np.arange(len(xx)));xx-=xx.mean();xx/=xx.std()
    checks['erro_recorte_antes_preparo']=float(np.max(np.abs(crop_prepare(raw,300)-xx)))
    assert checks['erro_recorte_antes_preparo']<1e-11
    checks['diferenca_preparo_antes_corte']=float(np.max(np.abs(prepare(raw)[300:2348]-xx)))
    # Prova de igualdade com as operações v5 extraídas da fonte sem executá-la.
    tree=ast.parse((ROOT.parent/'v5_spatial/wavelet_v5.py').read_text())
    chosen=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['biased_autocorrelation','biased_crosscorrelation']]
    scope={'np':np,'N':512};exec(compile(ast.Module(body=chosen,type_ignores=[]),'v5_literal','exec'),scope)
    x,y=rng.normal(size=(2,512));rx=scope['biased_autocorrelation'](x);ry=scope['biased_autocorrelation'](y)
    xy=scope['biased_crosscorrelation'](x,y);yx=scope['biased_crosscorrelation'](y,x)
    corr=[np.r_[rx[:0:-1],rx],np.r_[ry[:0:-1],ry],np.r_[yx[:0:-1],xy]]
    wl=np.geomspace(50,400,12);expected=[pywt.cwt(r,LEGACY_CF*wl,LEGACY_WAVELET,method='conv')[0] for r in corr]
    m=pair(x,y,wl,engine='legacy')
    checks['erro_v5_literal']=max(float(np.max(np.abs(v-m[k]))) for k,v in zip(['wxx','wyy','wxy'],expected))
    assert checks['erro_v5_literal']<1e-10
    # Não há limiar de sucesso para os controles negativos: uma WCO alta é reportada.
    smooth=np.exp(-.5*((t-1024)/160)**2);different=np.exp(-.5*((t-1024)/65)**2)
    background=np.sin(2*np.pi*t/500)+.5*np.cos(2*np.pi*t/170)
    cases={'evento_identico':(event,event),'atraso_100':(event,delayed),
           'evento_com_ruido':(event,event+.2*rng.normal(size=len(t))),
           'eventos_suaves_distintos':(smooth,different),
           'fundo_comum_eventos_distintos':(background+.25*smooth,background+.25*np.roll(different,350)),
           'ruidos_independentes':tuple(rng.normal(size=(2,len(t)))),
           'evento_cortado_na_borda':(smooth,np.exp(-.5*((t-20)/160)**2))}
    results=[]
    arrays={}
    for name,(u,v) in cases.items():
        a=prepare(120+u);b=prepare(120+v);m,s=evaluate(a,b)
        assert s['erro_equacao']<1e-12 and s['maximo_integral']<=1+1e-9
        results.append(dict(controle=name,**s));arrays[name+'_x']=a;arrays[name+'_y']=b
    pd.DataFrame(results).to_csv(out/'controles_sinteticos.csv',index=False)
    np.savez_compressed(out/'controles_sinteticos.npz',**arrays)
    checks['status']='PASS';checks['segundos']=time.perf_counter()-t0
    (out/'validacao_numerica.json').write_text(json.dumps(checks,indent=2,ensure_ascii=False)+'\n')
    timings=[]
    for n in [512,1024,2048]:
        x,y=rng.normal(size=(2,n));wl=np.geomspace(50,n/2,200)
        for engine in ['morlet','legacy']:
            samples=[]
            for _ in range(3):
                t0=time.perf_counter();pair(x,y,wl,engine=engine);samples.append(time.perf_counter()-t0)
            timings.append(dict(n=n,nucleo=engine,mediana_s=float(np.median(samples)),maximo_s=max(samples),
                                matrizes_principais_mib=200*(2*n-1)*56/1024**2))
    df=pd.DataFrame(timings);df.to_csv(out/'benchmark.csv',index=False)
    npairs=len(pd.read_csv(out/'pares_eventos.csv'))
    # 8 versões de eventos; 975 referências principais e 45 em outras 3 versões.
    nmorlet=npairs*7+975+90;nlegacy=npairs+45
    sec=nmorlet*df.query('n==2048 and nucleo=="morlet"').mediana_s.iloc[0]+nlegacy*df.query('n==2048 and nucleo=="legacy"').mediana_s.iloc[0]
    estimate=dict(pares_eventos=npairs,escalas=200,comparacoes_morlet=int(nmorlet),comparacoes_v5=int(nlegacy),
        segundos_estimados_nucleos=float(sec),minutos_com_margem_2x=float(sec*2/60),
        memoria_planejada_mib=350,nota='Estimativa inclui 8 versões dos eventos; controles self/intra e renderização separados.')
    (out/'estimativa.json').write_text(json.dumps(estimate,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(checks,indent=2));print(df.to_string(index=False));print(json.dumps(estimate,indent=2))

if __name__=='__main__':main()
