"""Conferência das fórmulas por somas explícitas, sem reutilizar a CWT de referência."""
from datetime import datetime, timezone
import hashlib
import json
import numpy as np
import pandas as pd
from analysis_common import ROOT, LAM, load_signals, evaluate
from wavelet_core import correlation, morlet_cwt, transform_correlation, pair

OUT=ROOT/'results/auditoria_adversarial'


def independent_prepare(raw):
    z=255-np.asarray(raw,float)
    design=np.column_stack((np.ones(len(z)),np.arange(len(z))))
    z=z-design@np.linalg.lstsq(design,z,rcond=None)[0]
    z-=z.mean()
    return z/np.sqrt(np.dot(z,z)/len(z))


def direct_cwt(corr, lags, at, scale):
    u=(lags-at)/scale
    kernel=np.pi**(-.25)/np.sqrt(scale)*np.exp(-.5*u*u)*np.exp(-6j*u)
    return np.dot(corr,kernel)


def main():
    _,signals=load_signals()
    chosen=[]
    for stage,folder,file,id_x,id_y,catfile,catid in [
        ('etapa1',ROOT/'results','resultados_eventos.csv','event_x','event_y','catalogo_eventos.csv','event_id'),
        ('extensao',ROOT/'results/extensao_churn_anular','resultados_ancoras.csv','anchor_x','anchor_y','catalogo_ancoras.csv','anchor_id')]:
        table=pd.read_csv(folder/file).query("variante=='principal' and pattern_x==pattern_y")
        catalog=pd.read_csv(folder/catfile).set_index(catid)
        for pattern,g in table.groupby('pattern_x'):
            for _,r in g.sort_values('pair_id').iloc[[0,-1]].iterrows():
                a,b=catalog.loc[r[id_x]],catalog.loc[r[id_y]]
                x=independent_prepare(signals[a.aquisicao]['principal'][int(a.segmento_inicio):int(a.segmento_fim)])
                y=independent_prepare(signals[b.aquisicao]['principal'][int(b.segmento_inicio):int(b.segmento_fim)])
                chosen.append((stage,pattern,r,x,y))
    errors=[];lams=LAM[[0,41,98,150,199]]
    coords=np.array([-2047,-1024,-32,0,32,1024,2047]);lags=np.arange(-2047,2048)
    for stage,pattern,r,x,y in chosen:
        maps,score=evaluate(x,y,wavelengths=lams)
        records=[np.correlate(a,b,'full')/2048 for a,b in [(x,x),(y,y),(x,y)]]
        maxerr=0.
        for key,rr in zip(['wxx','wyy','wxy'],records):
            direct=np.array([[direct_cwt(rr,lags,b,a) for b in coords] for a in maps['scales']])
            maxerr=max(maxerr,float(np.max(abs(direct-maps[key][:,coords+2047]))))
        assert maxerr<1e-10,(r.pair_id,maxerr)
        center=[np.array([direct_cwt(rr,lags,0,a) for a in maps['scales']]) for rr in records[:2]]
        cross=np.array([[direct_cwt(records[2],lags,b,a) for b in range(-32,33)] for a in maps['scales']])
        expected=abs(cross)**2/(abs(center[0])*abs(center[1]))[:,None]
        wcoerr=float(abs(expected-maps['wco'][:,2015:2080]).max())
        assert wcoerr<1e-10
        full,summary=evaluate(x,y)
        scoreerr=abs(summary['escore']-r.escore);assert scoreerr<1e-10
        swapped,_=evaluate(y,x,wavelengths=lams)
        symmetry=float(abs(maps['wxy']-np.conj(swapped['wxy'][:,::-1])).max());assert symmetry<1e-10
        pad=float(abs(morlet_cwt(records[2],maps['scales'],8)-morlet_cwt(records[2],maps['scales'],12)).max())
        assert pad<1e-10
        errors.append(dict(etapa=stage,padrao=pattern,pair_id=r.pair_id,erro_CWT_soma_direta=maxerr,
            erro_WCO_soma_direta=wcoerr,erro_escore_salvo=scoreerr,erro_troca_xy=symmetry,erro_padding=pad))
        del full,maps,swapped
    # O comentário direcional do v11 é contrário ao sinal que ele calcula.
    a=np.zeros(128);a[50]=1;b=np.zeros(128);b[57]=1
    pearson={s:np.corrcoef(a[s:] if s>=0 else a[:128+s],b[:128-s] if s>=0 else b[-s:])[0,1] for s in range(-20,21)}
    right_shift=max(pearson,key=pearson.get);assert right_shift==-7
    # Trocar a orientação do lag sem conjugar mantém a magnitude, mas inverte
    # a fase. Verificação tanto ideal quanto no kernel discreto legado.
    x,y=chosen[0][-2:]
    legacy,_=evaluate(x,y,engine='legacy',wavelengths=lams)
    canonical,_=transform_correlation(np.correlate(x,y,'full')/2048,lams,engine='legacy')
    phase=dict(erro_relativo_sem_conjugar=float(np.linalg.norm(legacy['wxy']-canonical)/np.linalg.norm(canonical)),
        erro_relativo_conjugando=float(np.linalg.norm(np.conj(legacy['wxy'])-canonical)/np.linalg.norm(canonical)),
        nota='A discretização PyWavelets também introduz assimetrias; a fase legada refletida não representa a convenção complexa principal.')
    # Relação exata entre os estimadores; diferença é metodológica declarada.
    rb=np.correlate(x,y,'full')/len(x);ru=np.correlate(x,y,'full')/(len(x)-abs(lags))
    rel=float(abs(rb-ru*(1-abs(lags)/len(x))).max());assert rel<1e-12
    # Ensemble de duas realizações de fases opostas: exemplo da divergência
    # entre E(CWT) e E(|CWT|) impressa nas Eqs.32.16/32.20.
    v=1+2j;ensemble_complex=abs((v-v)/2);ensemble_magnitude=(abs(v)+abs(-v))/2
    report=dict(status='PASS',data_utc=datetime.now(timezone.utc).isoformat(),
        comparacoes_independentes=errors,escalas_por_comparacao=5,lags_CWT_por_comparacao=7,lags_WCO_por_comparacao=65,
        v11_deslocamento_real_direita_7_retorna=right_shift,fase_legado=phase,relacao_enviesado_nao_enviesado_erro=rel,
        ensemble_fases_opostas=dict(modulo_media_complexa=ensemble_complex,media_modulos=ensemble_magnitude),
        fonte='../../assets/schuck-bardo-wavelet-coherence.pdf; páginas impressas 411, 413 e 414, conferidas visualmente',
        fonte_sha256=hashlib.sha256((ROOT.parent.parent/'assets/schuck-bardo-wavelet-coherence.pdf').read_bytes()).hexdigest(),
        limite='Conferência numérica e formal do estimador finito declarado. Não valida estacionariedade, validade física do panorama ou inferência populacional.')
    pd.DataFrame(errors).to_csv(OUT/'formulas_comparacoes.csv',index=False)
    (OUT/'formulas.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(report,indent=2,ensure_ascii=False))


if __name__=='__main__':main()
