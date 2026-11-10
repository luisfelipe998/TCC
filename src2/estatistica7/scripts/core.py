"""WCO integral e resumos por pares, sem efeitos de escrita nas fontes."""
from pathlib import Path
from itertools import combinations
import json
import hashlib
import sys

import numpy as np
import pandas as pd
from scipy.fft import fft, ifft, next_fast_len
from scipy.signal import detrend
from scipy.spatial.distance import pdist,squareform
from scipy.cluster.hierarchy import linkage,cut_tree

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT.parent
V16=SRC/'v16_spatial2_18_71'
sys.path.insert(0,str(V16/'scripts'))
from wavelet_core import prepare,correlation,morlet_cwt,FOURIER_FACTOR
CFG=json.loads((ROOT/'configs/wco_v1.json').read_text())
LAM=np.geomspace(50,3000,200)
PATTERNS=['churn','wavy annular','plug','slug','dispersed bubbles']
COLORS=dict(zip(PATTERNS,['#277da1','#dfaa3c','#8064a8','#31977b','#667580']))
T=np.log(LAM)
H=np.r_[(T[1]-T[0])/2,(T[2:]-T[:-2])/2,(T[-1]-T[-2])/2]
BANDS={'short':LAM<200,'middle':(LAM>=200)&(LAM<800),'long':LAM>=800}


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):digest.update(chunk)
    return digest.hexdigest()


def json_write(path,obj):
    def default(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,Path):return str(x)
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=default,allow_nan=False)+'\n')


def csv_read(path):return pd.read_csv(path,float_precision='round_trip')


def csv_save(name,rows):
    frame=rows if isinstance(rows,pd.DataFrame) else pd.DataFrame(rows)
    frame.to_csv(ROOT/'results'/name,index=False,float_format='%.17g')
    return frame


def load_inputs():
    with np.load(SRC/'estatistica5/todos_contra_todos/results/descritores_por_ensaio.npz',allow_pickle=False) as z:
        runs=z['runs'].astype(str);labels=z['patterns'].astype(str)
    inv=csv_read(SRC/'estatistica4/results/inventory.csv').set_index('run').loc[runs].reset_index()
    assert np.array_equal(inv.pattern,labels)
    raw={}
    for run in runs:
        with np.load(V16/'intermediate'/f'{run}_input.npz',allow_pickle=False) as z:raw[run]=z['raw'].copy()
    return runs,labels,inv,raw


def condition_inputs(condition,runs,labels,raw):
    length=int(condition.split('_')[0]);position=condition.split('_')[1]
    if length==3171:mask=np.ones(len(runs),dtype=bool)
    elif length==5071:mask=labels!='dispersed bubbles'
    elif length==7001:mask=np.isin(labels,['churn','wavy annular'])
    else:raise ValueError(condition)
    rr,ll=runs[mask],labels[mask]
    signals=[];coords=[]
    for run,label in zip(rr,ll):
        n=len(raw[run]);assert n>=length
        start=0 if position=='prefix' else ((n-length)//2 if position=='center' else n-length)
        signals.append(prepare(raw[run][start:start+length]))
        coords.append(dict(run=run,pattern=label,condition=condition,start=start,end_exclusive=start+length,n_original=n))
    return rr,ll,np.asarray(signals),coords


class Wavelet:
    """Convolução exata nas posições observadas, sem calcular regiões descartadas.

    Para correlação com M=2N-1 elementos, só são usados deslocamentos de kernel
    em [-(M-1),M-1]. Uma FFT de tamanho >=2M-1 representa esses deslocamentos
    sem colisão circular; os primeiros M valores são a convolução linear exata.
    Não é uma aproximação que encurta o suporte relevante para esses valores.
    """
    def __init__(self,n,wavelengths=LAM):
        self.n=n;self.m=2*n-1;self.nfft=next_fast_len(2*self.m-1)
        self.wavelengths=np.asarray(wavelengths);self.scales=self.wavelengths/FOURIER_FACTOR
        lag=np.arange(-self.m+1,self.m)
        tau=lag[None,:]/self.scales[:,None]
        kernel=np.pi**(-.25)/np.sqrt(self.scales[:,None])*np.exp(-.5*tau*tau)*np.exp(6j*tau)
        wrapped=np.zeros((len(self.scales),self.nfft),dtype=complex)
        wrapped[:,lag % self.nfft]=kernel
        self.kernel_fft=fft(wrapped,axis=1,workers=4)
        self.center=n-1
    def transform(self,corr):
        assert len(corr)==self.m
        return ifft(self.kernel_fft*fft(corr,self.nfft)[None,:],axis=1,workers=4)[:,:self.m]
    def auto(self,x):
        w=self.transform(correlation(x,x));p=np.abs(w[:,self.center])
        valid=p>1e-10*p.max()
        co=np.full(w.shape,np.nan)
        np.divide(np.abs(w)**2,p[:,None]**2,out=co,where=valid[:,None])
        return p,co
    def cross(self,x,y,px,py):
        w=self.transform(correlation(x,y))
        valid=(px>1e-10*px.max())&(py>1e-10*py.max())
        co=np.full(w.shape,np.nan)
        np.divide(np.abs(w)**2,(px*py)[:,None],out=co,where=valid[:,None])
        return co,valid


def normalized_weights(px,py,valid,powered=True,band=None):
    w=H.copy()
    if powered:w*=np.sqrt(px*py)
    w[~valid]=0
    if band is not None:w[~band]=0
    if w.sum()<=0:raise ValueError('Nenhuma escala válida no resumo')
    return w/w.sum()


def map_summary(co,px,py,n,quantiles=True):
    valid=np.isfinite(co).all(axis=1)
    if not valid.any():raise ValueError('Mapa sem escalas válidas')
    # Pesos zero nas linhas inválidas; não transformar NaN em significância.
    work=np.where(valid[:,None],co,0.)
    uniform=normalized_weights(px,py,valid,False)
    weight=normalized_weights(px,py,valid,True)
    peak=work.max(axis=1);mean=work.mean(axis=1);zero=work[:,n-1]
    curve=weight@work
    top=np.unravel_index(np.nanargmax(co),co.shape)
    out=dict(global_max=float(co[top]),zero_weighted=float(weight@zero),
             scale_peak_weighted=float(weight@peak),common_lag_peak=float(curve.max()),
             map_mean=float(uniform@mean),map_mean_weighted=float(weight@mean),
             area05_weighted=float(weight@(work>=.5).mean(axis=1)),
             peak_lambda_px=float(LAM[top[0]]),peak_lag_px=int(top[1]-n+1),
             common_peak_lag_px=int(np.argmax(curve)-n+1),
             fraction_valid_scales=float(valid.mean()),
             global_peak_power_x_fraction=float(px[top[0]]/px.max()),
             global_peak_power_y_fraction=float(py[top[0]]/py.max()),
             common_to_scale_peak=float(curve.max()/(weight@peak)))
    for threshold in [.01,.05,.10]:
        w=np.sqrt(px*py)*valid*(px>=threshold*px.max())*(py>=threshold*py.max())
        if w.sum()<=0:raise ValueError('Faixa de potência vazia')
        out[f'scale_peak_power{int(threshold*100):02d}']=float(w@peak/w.sum())
    if quantiles:out['lag_q90_weighted']=float(weight@np.quantile(work,.9,axis=1))
    for name,band in BANDS.items():
        u=normalized_weights(px,py,valid,False,band)
        w=normalized_weights(px,py,valid,True,band)
        out[f'mean_{name}']=float(u@mean)
        out[f'common_peak_{name}']=float(np.max(w@work))
    for key in CFG['cross_candidates']:
        if key in out and (out[key]<-1e-10 or out[key]>1+1e-8):raise ArithmeticError(f'{key} fora do limite teórico do estimador enviesado')
    profiles=dict(mean=mean,peak=peak,zero=zero,peak_lag=np.argmax(work,axis=1)-n+1)
    return out,profiles


def self_descriptors(co,n):
    assert np.isfinite(co).all()
    area=co.mean(axis=1)
    positive=co[:,n-1:]
    crossed=positive<.5
    censored=~crossed.any(axis=1)
    half=np.where(censored,n-1,np.argmax(crossed,axis=1)).astype(float)
    return area,half/LAM,censored


def silhouette(d,labels):
    labels=np.asarray(labels);s=[]
    for i,g in enumerate(labels):
        own=labels==g;own[i]=False
        if not own.any():s.append(np.nan);continue
        a=float(d[i,own].mean())
        b=min(float(d[i,labels==other].mean()) for other in np.unique(labels) if other!=g)
        s.append((b-a)/max(a,b) if max(a,b)>0 else 0.)
    return np.array(s)


def neighbor(d,labels):
    dd=d.copy();np.fill_diagonal(dd,np.inf)
    ties=np.isclose(dd,dd.min(axis=1,keepdims=True),rtol=1e-12,atol=1e-12)
    return (ties*(labels[:,None]==labels[None,:])).sum(axis=1)/ties.sum(axis=1)


def metric_summary(d,labels):
    groups=np.unique(labels);s=silhouette(d,labels);nn=neighbor(d,labels)
    class_s=[s[labels==g].mean() for g in groups]
    iu=np.triu_indices(len(labels),1);same=labels[iu[0]]==labels[iu[1]]
    chance=np.mean([(sum(labels==g)-1)/(len(labels)-1) for g in groups])
    balanced_nn=np.mean([nn[labels==g].mean() for g in groups])
    by_run=[]
    for i,g in enumerate(labels):
        own=labels==g;own[i]=False
        # Distância menor dentro de cada classe é o sinal favorável.
        by_run.append(min(d[i,labels==h].mean() for h in groups if h!=g)-d[i,own].mean())
    return dict(n_runs=len(labels),n_intra=int(same.sum()),n_inter=int((~same).sum()),
                macro_silhouette=float(np.mean(class_s)),min_class_silhouette=float(min(class_s)),
                balanced_neighbor=float(balanced_nn),neighbor_chance=float(chance),neighbor_lift=float(balanced_nn-chance),
                median_intra_distance=float(np.median(d[iu][same])),median_inter_distance=float(np.median(d[iu][~same])),
                balanced_mean_gap=float(np.mean([np.mean(np.array(by_run)[labels==g]) for g in groups])),
                positive_run_fraction=float(np.mean(np.array(by_run)>0)))


def screen_pass(row):
    gate=CFG['screen_gate']
    return bool(row['macro_silhouette']+1e-12>=gate['macro_silhouette_at_least'] and
                row['min_class_silhouette']>gate['every_class_silhouette_above'] and
                row['neighbor_lift']+1e-12>=gate['balanced_neighbor_lift_at_least'])


def ari(labels,clusters):
    a=np.array([labels[i]==labels[j] for i,j in combinations(range(len(labels)),2)])
    b=np.array([clusters[i]==clusters[j] for i,j in combinations(range(len(labels)),2)])
    expected=a.sum()*b.sum()/len(a);den=(a.sum()+b.sum())/2-expected
    return float(((a&b).sum()-expected)/den) if den else 1.


def cluster(d,labels):
    cl=cut_tree(linkage(squareform(d,checks=True),method='average'),n_clusters=len(np.unique(labels))).ravel()
    return cl,ari(labels,cl)


def holm(p):
    p=np.asarray(p);order=np.argsort(p)
    ranked=np.minimum(1,np.maximum.accumulate(p[order]*np.arange(len(p),0,-1)))
    out=np.empty_like(p);out[order]=ranked
    return out


def balanced_similarity_contrast(scores,pairs,labels):
    groups=np.unique(labels);within=[]
    for g in groups:
        mask=(labels[pairs[:,0]]==g)&(labels[pairs[:,1]]==g)
        within.append(scores[mask].mean())
    inter=labels[pairs[:,0]]!=labels[pairs[:,1]]
    return float(np.mean(within)-scores[inter].mean())
