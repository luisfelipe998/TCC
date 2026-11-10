"""Confere preservação, evidências e agregações da revisão adversarial."""
from pathlib import Path
from datetime import datetime,timezone
from collections import defaultdict
import hashlib
import json
import re
import subprocess
import numpy as np
import pandas as pd
from PIL import Image
from analysis_common import ROOT,load_signals,evaluate,crop_prepare

OUT=ROOT/'results/auditoria_adversarial'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(n):return json.loads((OUT/n).read_text())
def close(a,b):np.testing.assert_allclose(a,b,rtol=1e-11,atol=1e-11,equal_nan=True)


def main():
    checks={}
    prior=read('estado_anterior/results/manifesto.json')['saidas_sha256']
    protected={p:h for p,h in prior.items() if p.startswith(('results/','maps/')) and Path(p).suffix in ['.csv','.npz','.npy']}
    for p,h in protected.items():assert sha(ROOT/p)==h,p
    checks['tabelas_e_matrizes_anteriores_preservadas']=len(protected)
    for p,h in read('fontes_sha256.json').items():assert sha(ROOT.parent/p)==h,p
    checks['fontes_atuais_conferidas']=len(read('fontes_sha256.json'))
    for n in ['dados.json','ancoras.json','formulas.json','diagnosticos.json']:assert read(n)['status']=='PASS',n
    data=read('dados.json');assert (data['quadros'],data['deslocamentos'],data['panoramas_reproduzidos'])==(15531,15500,31)
    assert data['grupos_quadros_duplicados']==0
    plan=read('plano_diagnosticos.json');export=read('revisao_exportacao.json')
    for name,h in plan['scripts_sha256'].items():
        path=ROOT/export['codigo_executado'] if name=='adversarial_diagnostics.py' else ROOT/'scripts'/name
        assert sha(path)==h,name
    assert sha(ROOT/'scripts/adversarial_diagnostics.py')==export['codigo_entregue_sha256']
    for item in export['exportacoes']:
        assert item['escores_inalterados'] and sha(OUT/item['arquivo'])==item['sha256_depois']
    math=read('formulas.json')
    assert sha(ROOT.parent.parent/'assets/schuck-bardo-wavelet-coherence.pdf')==math['fonte_sha256']
    for row in math['comparacoes_independentes']:
        for k,v in row.items():
            if k.startswith('erro_'):assert v<1e-10,(k,v)
    assert len(math['comparacoes_independentes'])==8
    checks['evidencia_independente_das_formulas_conferida']=True
    selection=pd.read_csv(OUT/'repetibilidade_selecao_detalhes.csv')
    summary=pd.read_csv(OUT/'repetibilidade_selecao.csv')
    for r in summary.itertuples():
        g=selection[(selection.padrao==r.padrao)&(selection.tipo==r.tipo)&(selection.alternativa==r.alternativa)]
        assert len(g)==r.total
        assert g.candidato_em_64.sum()==r.candidatos_em_64
        assert g.selecionado_em_64.sum()==r.selecionados_em_64
        assert g.selecionado_em_16.sum()==r.selecionados_em_16
        assert (g.selecionado_em_64==(g.distancia_selecionado<=64)).all()
    assert len(selection)==152*3
    counts=summary[summary.tipo=='minimo'].groupby('alternativa').selecionados_em_64.sum().to_dict()
    assert counts=={'faixa30_49':55,'sigma12':66,'sigma20':68}
    checks['correspondencia_e_selecao_distintas_conferidas']=counts
    cp=pd.read_csv(ROOT/'results/extensao_churn_anular/contrastes_pareados.csv').query("variante=='principal'")
    for r in pd.read_csv(OUT/'cobertura_denominadores.csv').itertuples():
        g=cp[(cp.classe==r.padrao)&(cp.comparador==r.controle)]
        assert r.aquisicoes==len(set(g.run_x)|set(g.run_y)) and r.pares==len(g)
        close(r.fracao_pares,r.pares/r.pares_possiveis);close(r.fracao_aquisicoes,r.aquisicoes/r.aquisicoes_possiveis)
    scores=pd.read_csv(OUT/'todos_pares_admissiveis.csv');assert len(scores)==1095 and not scores.duplicated(['etapa','pair_id']).any()
    assert scores.escore.between(0,1).all()
    fixed=pd.read_csv(ROOT/'results/resultados_fixos.csv').query("variante=='principal' and pattern_x==pattern_y")
    fixed_values=defaultdict(list)
    for r in fixed.itertuples():fixed_values[(r.pattern_x,r.run_x,r.run_y)].append(r.escore)
    fixed_means={k:sum(v)/len(v) for k,v in fixed_values.items()}
    seeds=pd.read_csv(OUT/'sensibilidade_semente.csv')
    # Reconstrução independente: seleção por hash e médias em dicionários.
    tested=0
    for stage,folder,file,compatible,distinct,cap in [
        ('etapa1',ROOT/'results','resultados_eventos.csv','intra_equivalentes','intra_distintos',12),
        ('extensao',ROOT/'results/extensao_churn_anular','resultados_ancoras.csv','minimos_compativeis','minimos_distintos',2)]:
        candidates=pd.read_csv(folder/'pares_candidatos.csv')
        candidates=candidates[(candidates.pattern_x==candidates.pattern_y)&candidates.grupo.isin([compatible,distinct])]
        ss=scores[scores.etapa==stage].set_index('pair_id');assert set(candidates.pair_id)==set(ss.index)
        old=pd.read_csv(folder/file).query("variante=='principal'").set_index('pair_id')
        for pid in ss[ss.origem=='resultado_congelado'].index:close(ss.loc[pid,'escore'],old.loc[pid,'escore'])
        basegroups=defaultdict(list)
        for r in candidates.itertuples():basegroups[(r.run_x,r.run_y,r.grupo)].append(r)
        needed=seeds[seeds.etapa==stage]
        for (variant,seed),block in needed.groupby(['variante','semente'],dropna=False):
            selected=[]
            for g in basegroups.values():
                if variant=='todos_admissiveis':selected+=g
                else:
                    def priority(r):
                        fields=[str(int(seed)),r.pair_id,r.grupo] if stage=='etapa1' else [str(int(seed)),r.grupo,r.pair_id]
                        return hashlib.sha256('|'.join(fields).encode()).hexdigest()
                    selected+=sorted(g,key=priority)[:cap]
            values=defaultdict(list)
            for r in selected:values[(r.pattern_x,r.grupo,r.run_x,r.run_y)].append(ss.loc[r.pair_id,'escore'])
            means={k:sum(v)/len(v) for k,v in values.items()}
            for row in block.itertuples():
                differences=[]
                for (pat,group,rx,ry),value in means.items():
                    if pat!=row.padrao or group!=compatible:continue
                    control=fixed_means.get((pat,rx,ry)) if row.comparador=='fixos' else means.get((pat,distinct,rx,ry))
                    if control is not None:differences.append(value-control)
                assert len(differences)==row.pares
                close(np.mean(differences),row.diferenca);tested+=1
    for r in pd.read_csv(OUT/'sensibilidade_semente_resumo.csv').itertuples():
        x=seeds[(seeds.variante=='semente_alternativa')&(seeds.etapa==r.etapa)&(seeds.padrao==r.padrao)&(seeds.comparador==r.comparador)].diferenca.to_numpy()
        assert len(x)==r.sementes==100
        close([x.mean(),x.min(),x.max(),np.quantile(x,.025),np.quantile(x,.975),(x>0).mean()],
            [r.media,r.minimo,r.maximo,r.p025,r.p975,r.fracao_positiva])
    checks['contrastes_sementes_e_todos_pares_recalculados']=tested
    changed=pd.read_csv(OUT/'ancoras_sem_feather.csv');fx=pd.read_csv(OUT/'fixos_sem_feather.csv')
    assert len(changed)==458 and len(fx)==261
    assert set(fx.columns)=={'pair_id','run_x','run_y','pattern_x','pattern_y','grupo','representacao','escore_original','escore'}
    assert set(changed.columns)==set(fx.columns)|{'etapa'}
    seam_contrasts=pd.read_csv(OUT/'contrastes_sem_feather.csv')
    for r in seam_contrasts.itertuples():
        metric='escore_original' if r.variante=='original' else 'escore'
        g=changed[(changed.etapa==r.etapa)&(changed.pattern_x==r.padrao)]
        compatible='intra_equivalentes' if r.etapa=='etapa1' else 'minimos_compativeis'
        a=g[g.grupo==compatible].groupby(['run_x','run_y'])[metric].mean()
        b=fx[fx.pattern_x==r.padrao].groupby(['run_x','run_y'])[metric].mean() if r.comparador=='fixos' else g[g.grupo==r.comparador].groupby(['run_x','run_y'])[metric].mean()
        ids=a.index.intersection(b.index);assert len(ids)==r.pares;close((a.loc[ids]-b.loc[ids]).mean(),r.diferenca)
    null=pd.read_csv(OUT/'nulo_selecao_completa.csv');n=read('diagnosticos.json')['etapas']['synthetic']
    assert len(null)==null.repeticao.nunique()==n['repeticoes']==40
    assert null.pares_ancorados.between(1,2).all() and null[['ancoras','fixos']].notna().all().all()
    close(null.ancoras-null.fixos,null.diferenca)
    close([null.ancoras.mean(),null.fixos.mean(),null.diferenca.mean(),(null.diferenca>0).mean()],
          [n['ancoras_media'],n['fixos_media'],n['diferenca_media'],n['fracao_positiva']])
    checks['costura_e_nulo_conferidos']=True
    for p in (ROOT/'figures/auditoria_adversarial').glob('*.png'):
        with Image.open(p) as im:assert min(im.size)>=400;im.verify()
    for name in ['EXECUTIVE_SUMMARY.md','RESULTS.md','PASSO_A_PASSO.md','ANCORAS_E_AUTOCORRELACAO.md','README.md','RETOMADA.md','DECISIONS.md']:
        text=(ROOT/name).read_text();assert 'AUDITORIA_ADVERSARIAL.md' in text,name
    text=(ROOT/'EXECUTIVE_SUMMARY.md').read_text()
    for phrase in ['66/76','68/76','55/76','10/10','8/9','0,0169']:assert phrase in text,phrase
    assert 'abaixo da referência de cobertura de 2/3' not in text
    pdf=ROOT/'RESUMO_EXECUTIVO.pdf';blob=pdf.read_bytes()
    pages=len(re.findall(rb'/Type\s*/Page\b',blob));assert pages==25,pages
    extracted=subprocess.run(['pdftotext',str(pdf),'-'],check=True,capture_output=True).stdout.decode()
    for phrase in ['Revisão adversarial','66/76','68/76','0,0169','20/100','15.531']:assert phrase in extracted,phrase
    checks['documentos_principais_e_pdf_atualizados']=True;checks['paginas_pdf']=pages
    review=read('revisao_apresentacao.json');assert review['estado']=='REVISADO'
    for p,h in review['sha256'].items():assert sha(ROOT/p)==h,p
    checks.update(status='PASS',data_utc=datetime.now(timezone.utc).isoformat(),
        conclusao='Análise exploratória com ressalvas. Integridade e consistência verificadas não equivalem a validade física ou ausência garantida de falhas.')
    (OUT/'verificacao_final.json').write_text(json.dumps(checks,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(checks,indent=2,ensure_ascii=False))


if __name__=='__main__':main()
