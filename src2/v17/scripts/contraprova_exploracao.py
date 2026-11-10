"""Independent review of the coordinator's Welch result; no source imports.

Checks numerical identity using manual Welch periodograms, then challenges
quantile and mixture conventions and segment length/detrending. These are
post-result adversarial sensitivities, not a search for a replacement winner.
Writes exclusively v17/resultados/contraprova_exploracao.
"""
from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT/'src2/v17/resultados/contraprova_exploracao'
OUT.mkdir(parents=True, exist_ok=True)
SOURCE = ROOT/'src2/v17/resultados/exploracao'


def save(name, obj):
    (OUT/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n')


def prepare(raw):
    raw = 255-np.asarray(raw,dtype=float)
    t = np.linspace(-1.,1.,len(raw))
    z = raw-np.polynomial.polynomial.polyval(t,np.polynomial.polynomial.polyfit(t,raw,1))
    z -= z.mean()
    return z/np.sqrt(np.mean(z*z))


def manual_welch(x, length=1024, detrend_segment=False):
    w = .5-.5*np.cos(2*np.pi*np.arange(length)/length)
    spectra = []
    powers = []
    for start in range(0,len(x)-length+1,length//2):
        v = x[start:start+length].copy()
        if detrend_segment: v -= v.mean()
        v = v*w
        spectra.append(abs(np.fft.rfft(v,n=8192))**2 / np.dot(w,w))
        powers.append(float(np.dot(v,v)/np.dot(w,w)))
    psd = np.mean(spectra,axis=0)
    psd[1:-1] *= 2
    assert abs(psd.sum()/8192-np.mean(powers)) < 1e-12
    return np.fft.rfftfreq(8192),psd


def weights(t):
    return np.diff(np.r_[t[0],(t[:-1]+t[1:])/2,t[-1]])


def widths_cdf(profile,lam):
    t=np.log(lam)
    masses=profile*weights(t)
    masses/=masses.sum(axis=-1,keepdims=True)
    result=[]
    for mass in masses:
        c=np.r_[0,np.cumsum(mass)]
        q=np.interp([.25,.75],c,np.r_[t[0],t])
        result.append(float(q[1]-q[0]))
    return np.asarray(result)


def widths_continuous(profile,lam):
    t=np.log(lam);dx=np.diff(t);result=[]
    for p in profile:
        c=np.r_[0,np.cumsum(dx*(p[:-1]+p[1:])/2)]
        q=[]
        for target in [.25*c[-1],.75*c[-1]]:
            k=int(np.searchsorted(c,target,side='right')-1)
            k=min(k,len(dx)-1)
            need=target-c[k];s=(p[k+1]-p[k])/dx[k]
            u=2*need/(p[k]+np.sqrt(max(0,p[k]**2+2*s*need)))
            q.append(t[k]+u)
        result.append(q[1]-q[0])
    return np.asarray(result)


def effect(width,y,name):
    sel=np.isin(y,['churn','wavy annular']);v=width[sel];y=y[sel]
    a=v[y=='churn'];b=v[y=='wavy annular']
    s=[]
    for i in range(len(v)):
        own=v[(y==y[i])&(np.arange(len(v))!=i)]
        other=v[y!=y[i]]
        ai=np.mean(abs(v[i]-own));bi=np.mean(abs(v[i]-other))
        s.append((bi-ai)/max(ai,bi) if max(ai,bi)>0 else 0)
    sil=[np.mean(np.asarray(s)[y==g]) for g in np.unique(y)]
    return dict(scenario=name,macro_silhouette=float(np.mean(sil)),min_class_silhouette=float(min(sil)),
                churn_min=float(a.min()),churn_max=float(a.max()),annular_min=float(b.min()),annular_max=float(b.max()),
                empirical_probability_churn_gt_annular=float(((a[:,None]>b)+.5*(a[:,None]==b)).mean()),
                cliff_delta=float(np.sign(a[:,None]-b).mean()))


def main():
    save('protocolo_contraprova.json',dict(status='Pós-resultado positivo, antes destes cálculos de contraprova',
        checks=['Welch manual FFT versus espectros arquivados','Parseval por segmento','mistura de espectros normalizados versus pooling antes de normalizar',
                'quantil interpolado da CDF versus integral do interpolante linear e grade 2000',
                'segmentos 512/1024/2048 com overlap 50% e Hann periódica','remoção da média por segmento como sensibilidade'],
        limit='Somente descrições, sem novos p-valores ou classificadores; não escolher a maior separação como novo resultado principal'))
    index=pd.read_csv(ROOT/'src2/estatistica8/dados/inventario.csv')
    labels=index.padrao.to_numpy();runs=index.ensaio.to_numpy()
    stored=np.load(SOURCE/'espectros_janelas.npz')
    records=[];output=[];width_rows=[];errors=[]
    configs=[(1024,False),(512,False),(2048,False),(1024,True)]
    allprofiles={}
    for nseg,detr in configs:
        for cap in [512,1024]:
            lam=np.geomspace(50,cap,200);dense=np.geomspace(50,cap,2000)
            profiles=[];dense_profiles=[]
            for i,row in enumerate(index.itertuples()):
                p=ROOT/'src2/data_processed'/row.ensaio/'__panorama_v11.png'
                with Image.open(p) as im: raw=np.asarray(im.convert('L'))[18:72].mean(axis=0)
                windows=[];dwindows=[]
                for name,start in [('prefixo',0),('centro',(len(raw)-2048)//2),('sufixo',len(raw)-2048)]:
                    f,psd=manual_welch(prepare(raw[start:start+2048]),nseg,detr)
                    profile=np.interp(1/lam,f,psd)/lam
                    dense_profile=np.interp(1/dense,f,psd)/dense
                    windows.append(profile);dwindows.append(dense_profile)
                    if nseg==1024 and not detr:
                        errors.append(float(abs(profile-stored[f'{cap}_{name}_welch_logenergia'][i]).max()))
                profiles.append(windows);dense_profiles.append(dwindows)
            profiles=np.asarray(profiles);dense_profiles=np.asarray(dense_profiles)
            hw=weights(np.log(lam));hd=weights(np.log(dense))
            norm=profiles/np.sum(profiles*hw,axis=2,keepdims=True)
            mixture=norm.mean(axis=1)
            pooled=profiles.mean(axis=1)
            dense_norm=dense_profiles/np.sum(dense_profiles*hd,axis=2,keepdims=True)
            dense_mixture=dense_norm.mean(axis=1)
            tag=f'cap{cap}_nseg{nseg}_detrend{detr}'
            w=widths_cdf(mixture,lam)
            output.append(effect(w,labels,tag+'_original_mixture_cdf'))
            if nseg==1024 and not detr:
                for method,vals in [('continuous200',widths_continuous(mixture,lam)),
                                    ('continuous2000',widths_continuous(dense_mixture,dense)),
                                    ('pool_then_normalize',widths_cdf(pooled,lam))]:
                    output.append(effect(vals,labels,tag+'_'+method))
                expected=pd.read_csv(SOURCE/'largura_todas_representacoes.csv')
                expected=expected[(expected.teto_lambda==cap)&(expected.janela=='media_tres_janelas')&(expected.metodo=='welch_logenergia')]
                ev=expected.set_index('ensaio').loc[runs,'largura'].to_numpy()
                assert np.max(abs(w-ev)) < 1e-12
                # Mixture must give exactly 1/3 normalized probability to each window.
                masses=norm*hw
                assert np.allclose(masses.sum(axis=2),1)
                assert np.allclose((mixture*hw).sum(axis=1),1)
                for run,label,value in zip(runs,labels,w):
                    width_rows.append(dict(run=run,pattern=label,cap=cap,width=value))
    pd.DataFrame(output).to_csv(OUT/'contraprovas_welch.csv',index=False)
    pd.DataFrame(width_rows).to_csv(OUT/'reproducao_largura_welch.csv',index=False)
    # Exact algebra: a lambda^-1/2 reweighting removes only a fitted linear term
    # in log-lambda for individual-window Pearson residuals (apart from epsilon).
    duplicate_errors=[]
    for cap in [512,1024]:
        lam=np.geomspace(50,cap,200);a=np.column_stack([np.ones(200),np.log(lam)])
        for w in ['ancora','prefixo','centro','sufixo']:
            residuals=[]
            for method in ['wps','wps_ponderada']:
                p=stored[f'{cap}_{w}_{method}'];lp=np.log(p+1e-12)
                z=lp-(a@np.linalg.lstsq(a,lp.T,rcond=None)[0]).T
                residuals.append(np.corrcoef(z))
            duplicate_errors.append(dict(cap=cap,window=w,max_abs_pearson_difference=float(abs(residuals[0]-residuals[1]).max())))
    pd.DataFrame(duplicate_errors).to_csv(OUT/'redundancia_pearson_reponderacao.csv',index=False)
    save('verificacao.json',dict(max_abs_manual_welch_error=max(errors),profiles_checked=len(errors),
        manual_welch_matches=max(errors)<1e-10,parseval_checks='All passed',n_new_classifiers=0,
        script_sha256=hashlib.file_digest(Path(__file__).open('rb'),'sha256').hexdigest()))
    print(pd.DataFrame(output).to_string(index=False))
    print('max error',max(errors))


if __name__=='__main__': main()
