"""Diagnósticos posteriores aos resultados, sem substituir o experimento congelado."""
from datetime import datetime, timezone
import argparse
import hashlib
import itertools
import json
import time
import numpy as np
import pandas as pd
from scipy.signal import lfilter
from analysis_common import ROOT, load_signals
from texture_common import Evaluator
from catalog_texture_anchors import detect, select

OUT=ROOT/'results/auditoria_adversarial'
EXT=ROOT/'results/extensao_churn_anular'
NSEEDS=100
NNULL=40
NULL_SEED=20261026


def datasets():
    _,signals=load_signals()
    with np.load(OUT/'perfis_sem_feather.npz') as z:
        for run in signals:signals[run]['sem_feather']=z[run+'_sem_feather'].copy()
    stages=[]
    for stage,folder,file,catfile,catid,xid,yid,seed,cap,compatible,distinct in [
        ('etapa1',ROOT/'results','resultados_eventos.csv','catalogo_eventos.csv','event_id','event_x','event_y',20261024,12,'intra_equivalentes','intra_distintos'),
        ('extensao',EXT,'resultados_ancoras.csv','catalogo_ancoras.csv','anchor_id','anchor_x','anchor_y',20261025,2,'minimos_compativeis','minimos_distintos')]:
        raw=pd.read_csv(folder/file).query("variante=='principal' and pattern_x==pattern_y").copy()
        candidates=pd.read_csv(folder/'pares_candidatos.csv')
        candidates=candidates[(candidates.pattern_x==candidates.pattern_y)&candidates.grupo.isin([compatible,distinct])].copy()
        stages.append(dict(stage=stage,raw=raw,candidates=candidates,cat=pd.read_csv(folder/catfile).set_index(catid),
            xid=xid,yid=yid,seed=seed,cap=cap,compatible=compatible,distinct=distinct))
    fixed=pd.read_csv(ROOT/'results/resultados_fixos.csv').query("variante=='principal' and pattern_x==pattern_y").copy()
    return signals,stages,fixed


def key(row,stage,side,band):
    anchor=stage['cat'].loc[row[stage[side+'id']]]
    return (anchor.aquisicao,int(anchor.segmento_inicio),band,True,'morlet')


def mean_by_run_pair(df):
    return df.groupby(['pattern_x','grupo','run_x','run_y']).escore.mean()


def contrast_rows(raw,fixed,stage,kind):
    a=raw.groupby(['pattern_x','grupo','run_x','run_y'],as_index=False).escore.mean()
    f=fixed.groupby(['pattern_x','run_x','run_y'],as_index=False).escore.mean()
    rows=[]
    for pattern,g in a.groupby('pattern_x'):
        target=g[g.grupo==stage['compatible']]
        for comparator in ['fixos',stage['distinct'],'polaridades_opostas']:
            control=f[f.pattern_x==pattern] if comparator=='fixos' else g[g.grupo==comparator]
            paired=target.merge(control,on=['run_x','run_y'],suffixes=('_a','_c'),validate='one_to_one')
            if len(paired):rows.append(dict(etapa=stage['stage'],padrao=pattern,comparador=comparator,variante=kind,
                diferenca=float((paired.escore_a-paired.escore_c).mean()),pares=len(paired)))
    return rows


def admissible(signals,stages,fixed):
    ev=Evaluator(signals);rows=[];contrasts=[]
    for stage in stages:
        stored=stage['raw'].set_index('pair_id')
        for i,(_,r) in enumerate(stage['candidates'].sort_values('pair_id').iterrows()):
            if r.pair_id in stored.index:
                score=float(stored.loc[r.pair_id].escore);origin='resultado_congelado'
            else:
                m,s,_,_=ev.pair(key(r,stage,'x','principal'),key(r,stage,'y','principal'));del m
                score=s['escore'];origin='diagnostico_posterior'
            rows.append(dict(etapa=stage['stage'],pair_id=r.pair_id,pattern_x=r.pattern_x,
                run_x=r.run_x,run_y=r.run_y,grupo=r.grupo,escore=score,origem=origin))
            if (i+1)%100==0:print(f'Todos os admissíveis {stage["stage"]}: {i+1}',flush=True)
        all_scores=pd.DataFrame([r for r in rows if r['etapa']==stage['stage']])
        contrasts+=contrast_rows(all_scores,fixed,stage,'todos_admissiveis')
        groups=list(stage['candidates'].groupby(['run_x','run_y','grupo']))
        for seed in range(stage['seed'],stage['seed']+NSEEDS+1):
            selected=[]
            for _,g in groups:
                def priority(r):
                    msg=f'{seed}|{r.pair_id}|{r.grupo}' if stage['stage']=='etapa1' else f'{seed}|{r.grupo}|{r.pair_id}'
                    return hashlib.sha256(msg.encode()).hexdigest()
                h=g.copy();h['hash_auditoria']=[priority(r) for r in h.itertuples()]
                selected.extend(h.sort_values('hash_auditoria').head(stage['cap']).pair_id)
            if seed==stage['seed']:
                expected=set(stage['raw'][stage['raw'].grupo.isin([stage['compatible'],stage['distinct']])].pair_id)
                assert set(selected)==expected
            chosen=all_scores[all_scores.pair_id.isin(selected)]
            new=contrast_rows(chosen,fixed,stage,'semente_original' if seed==stage['seed'] else 'semente_alternativa')
            for row in new:row['semente']=seed
            contrasts+=new
        pd.DataFrame(rows).to_csv(OUT/'todos_pares_admissiveis.csv',index=False)
        pd.DataFrame(contrasts).to_csv(OUT/'sensibilidade_semente.csv',index=False)
    df=pd.DataFrame(contrasts);g=df[df.variante=='semente_alternativa']
    summary=g.groupby(['etapa','padrao','comparador'],as_index=False).agg(
        sementes=('semente','size'),media=('diferenca','mean'),minimo=('diferenca','min'),maximo=('diferenca','max'),
        p025=('diferenca',lambda x:x.quantile(.025)),p975=('diferenca',lambda x:x.quantile(.975)),
        fracao_positiva=('diferenca',lambda x:float((x>0).mean())))
    summary.to_csv(OUT/'sensibilidade_semente_resumo.csv',index=False)
    return dict(pares_admissiveis=len(rows),novos=sum(r['origem']=='diagnostico_posterior' for r in rows),
        sementes_alternativas=NSEEDS,nota='Distribuição de decisões de amostragem de parceiros; não intervalo de confiança nem replicações experimentais.')


def seams(signals,stages,fixed):
    ev=Evaluator(signals);rows=[];contrasts=[];fixed_rows=[]
    segments=pd.read_csv(ROOT/'results/segmentos_fixos.csv').set_index('segment_id')
    for _,r in fixed.sort_values('pair_id').iterrows():
        x,y=r.pair_id.split('__');a,b=segments.loc[x],segments.loc[y]
        m,s,_,_=ev.pair((a.aquisicao,int(a.inicio),'sem_feather',True,'morlet'),(b.aquisicao,int(b.inicio),'sem_feather',True,'morlet'));del m
        row=r.to_dict();row['escore_original']=r.escore;row['escore']=s['escore'];fixed_rows.append(row)
    new_fixed=pd.DataFrame(fixed_rows);new_fixed.to_csv(OUT/'fixos_sem_feather.csv',index=False)
    for stage in stages:
        contrasts+=contrast_rows(stage['raw'],fixed,stage,'original')
        for i,(_,r) in enumerate(stage['raw'].sort_values('pair_id').iterrows()):
            m,s,_,_=ev.pair(key(r,stage,'x','sem_feather'),key(r,stage,'y','sem_feather'));del m
            row=r.to_dict();row.update(etapa=stage['stage'],escore_original=r.escore,escore=s['escore']);rows.append(row)
            if (i+1)%100==0:print(f'Sem feather {stage["stage"]}: {i+1}',flush=True)
        new=pd.DataFrame([r for r in rows if r['etapa']==stage['stage']])
        contrasts+=contrast_rows(new,new_fixed,stage,'sem_feather')
    pd.DataFrame(rows).to_csv(OUT/'ancoras_sem_feather.csv',index=False)
    pd.DataFrame(contrasts).to_csv(OUT/'contrastes_sem_feather.csv',index=False)
    return dict(pares_ancorados=len(rows),pares_fixos=len(fixed_rows),
        nota='Contrafactual da operação de feather, com deslocamentos, âncoras e parceiros fixos; não reconstrói o escoamento real nem refaz a seleção óptica.')


def synthetic():
    rng=np.random.default_rng(NULL_SEED);records=[]
    for trial in range(NNULL):
        local={};anchors={}
        for name in ['A','B']:
            x=lfilter([1.],[1.,-.95],rng.normal(size=8192+1024))[1024:]
            raw=128+20*(x-x.mean())/x.std()
            local[name]={'principal':raw}
            candidates,_=detect(raw,16)
            selected=candidates.loc[sorted(select(candidates))]
            anchors[name]=selected[selected.tipo=='minimo']
        pairs=[]
        for a,b in itertools.product(anchors['A'].itertuples(),anchors['B'].itertuples()):
            if max(a.largura,b.largura)/min(a.largura,b.largura)<=1.5:
                pid=f'sim{trial:03d}A_MIN{a.ancora:05d}__sim{trial:03d}B_MIN{b.ancora:05d}'
                pairs.append((hashlib.sha256(f'20261025|minimos_compativeis|{pid}'.encode()).hexdigest(),a.ancora,b.ancora))
        pairs=sorted(pairs)[:2];ev=Evaluator(local)
        anchored=[];controls=[]
        for _,a,b in pairs:
            m,s,_,_=ev.pair(('A',a-1024,'principal',True,'morlet'),('B',b-1024,'principal',True,'morlet'));del m
            anchored.append(s['escore'])
        for start in [0,3072,6144]:
            m,s,_,_=ev.pair(('A',start,'principal',True,'morlet'),('B',start,'principal',True,'morlet'));del m
            controls.append(s['escore'])
        aa=float(np.mean(anchored)) if anchored else np.nan;ff=float(np.mean(controls))
        records.append(dict(repeticao=trial,minimos_A=len(anchors['A']),minimos_B=len(anchors['B']),
            pares_ancorados=len(pairs),ancoras=aa,fixos=ff,diferenca=aa-ff))
        if (trial+1)%10==0:print(f'Nulo com detecção completa: {trial+1}/{NNULL}',flush=True)
    df=pd.DataFrame(records);df.to_csv(OUT/'nulo_selecao_completa.csv',index=False)
    valid=df.dropna(subset=['diferenca'])
    return dict(repeticoes=NNULL,calculaveis=len(valid),semente=NULL_SEED,phi=.95,N_panorama=8192,
        ancoras_media=float(valid.ancoras.mean()),fixos_media=float(valid.fixos.mean()),
        diferenca_media=float(valid.diferenca.mean()),diferenca_minima=float(valid.diferenca.min()),diferenca_maxima=float(valid.diferenca.max()),
        fracao_positiva=float((valid.diferenca>0).mean()),
        nota='Dois AR(1) gaussianos independentes por repetição, sem evento físico compartilhado. Detector e seleção completos antes da WCO. Diagnóstico posterior; não teste nulo calibrado para os dados experimentais.')


def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=['all','admissible','seams','synthetic'],default='all');args=p.parse_args()
    signals,stages,fixed=datasets()
    meta_path=OUT/'diagnosticos.json';meta=json.loads(meta_path.read_text()) if meta_path.exists() else {'etapas':{}}
    for name,fun in [('admissible',lambda:admissible(signals,stages,fixed)),('seams',lambda:seams(signals,stages,fixed)),('synthetic',synthetic)]:
        if args.stage in ['all',name]:
            t=time.perf_counter();result=fun();result['segundos']=time.perf_counter()-t
            meta['etapas'][name]=result;meta['status']='PASS' if set(meta['etapas'])=={'admissible','seams','synthetic'} else 'PARCIAL'
            meta['data_utc']=datetime.now(timezone.utc).isoformat()
            meta_path.write_text(json.dumps(meta,indent=2,ensure_ascii=False)+'\n')
            print(name,json.dumps(result,ensure_ascii=False),flush=True)
    print(json.dumps(meta,indent=2,ensure_ascii=False))


if __name__=='__main__':main()
