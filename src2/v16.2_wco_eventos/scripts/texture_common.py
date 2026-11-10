"""Recortes da extensão e cache limitado; equações herdadas sem alteração."""
from collections import OrderedDict
import hashlib
import json
import numpy as np
import pandas as pd
from analysis_common import ROOT, LAM, crop_prepare, evaluate, summarize
from wavelet_core import correlation, transform_correlation, safe_coherence

OUT=ROOT/'results/extensao_churn_anular'
MAPS=ROOT/'maps/extensao_churn_anular'
CFG=json.loads((ROOT/'config_ancoras_textura.json').read_text())
LABELS={'churn':'Agitado (churn)','wavy annular':'Anular ondulado (wavy annular)',
        'entre':'Churn × anular ondulado'}
GROUPS={'minimos_compativeis':'Mínimos de larguras compatíveis',
        'minimos_distintos':'Mínimos de larguras distintas',
        'polaridades_opostas':'Mínimo × máximo'}
VARIANT_LABELS={'principal':'Principal: linhas 18–71', 'faixa30_49':'Linhas 30–49',
    'somente_media':'Remoção apenas da média', 'centros_menos_mais16':'Centros x −16 / y +16 px',
    'centros_mais_menos16':'Centros x +16 / y −16 px', 'detector_sigma12':'Detector σ=12 px',
    'detector_sigma20':'Detector σ=20 px', 'v5_mesmos_segmentos':'CWT da v5; mesmos recortes'}

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def verify_freeze():
    for freeze in [ROOT/'results/congelamento.json',OUT/'congelamento.json']:
        for p,h in json.loads(freeze.read_text())['sha256'].items():
            assert sha(ROOT/p)==h,p

def spec(row,variant,side='x'):
    band='faixa30_49' if variant=='faixa30_49' else 'principal'
    linear=variant!='somente_media'
    engine='legacy' if variant=='v5_mesmos_segmentos' else 'morlet'
    center=row.ancora
    if variant.startswith('detector_sigma'):
        center=row['ancora_'+variant.replace('detector_','')]
    if pd.isna(center): return None
    if variant=='centros_menos_mais16':center+=-16 if side=='x' else 16
    if variant=='centros_mais_menos16':center+=16 if side=='x' else -16
    return (row.aquisicao,int(center)-1024,band,linear,engine)

class Evaluator:
    """LRU de até 16 autocorrelações transformadas (~200 MiB).

    Não guarda pares nem escores entre chamadas: apenas reutiliza Wx,x do
    mesmo recorte e tratamento. A implementação v5 é avaliada diretamente.
    """
    def __init__(self,signals,capacity=16):
        self.signals=signals;self.capacity=capacity;self.cache=OrderedDict()
        self.hits=0;self.misses=0

    def raw_prepared(self,key):
        run,start,band,linear,_=key
        return crop_prepare(self.signals[run][band],start,linear=linear)

    def auto(self,key):
        if key in self.cache:
            self.hits+=1;self.cache.move_to_end(key);return self.cache[key]
        self.misses+=1
        x=self.raw_prepared(key)
        w,a=transform_correlation(correlation(x,x),LAM,engine='morlet')
        for arr in [x,w,a]:arr.setflags(write=False)
        self.cache[key]=(x,w,a)
        if len(self.cache)>self.capacity:self.cache.popitem(last=False)
        return x,w,a

    def pair(self,kx,ky):
        assert kx[-1]==ky[-1]
        if kx[-1]=='legacy':
            x=self.raw_prepared(kx);y=self.raw_prepared(ky);m,s=evaluate(x,y,engine='legacy')
        else:
            x,wxx,a=self.auto(kx);y,wyy,b=self.auto(ky)
            assert np.array_equal(a,b)
            wxy,_=transform_correlation(correlation(x,y),LAM)
            co,valid=safe_coherence(wxy,wxx,wyy,len(x)-1)
            m=dict(wxx=wxx,wyy=wyy,wxy=wxy,wco=co,scales=a,
                   lags=np.arange(-len(x)+1,len(x)),denominator_valid=valid)
            s=summarize(m,len(x))
        assert s['erro_equacao']<1e-10
        assert s['escalas_validas']>0
        if kx[-1]=='morlet':assert s['maximo_integral']<=1+1e-8
        return m,s,x,y

def class_name(row):
    return row.pattern_x if row.pattern_x==row.pattern_y else 'entre'
