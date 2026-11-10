"""Auditoria final da extensão, sem refazer o lote completo."""
from pathlib import Path
from datetime import datetime,timezone
import contextlib
import io
import json
import re
import subprocess
import tempfile
import numpy as np
import pandas as pd
from PIL import Image
from analysis_common import ROOT,DATA,BASE,load_signals,evaluate,crop_prepare
from texture_common import OUT,MAPS,CFG,sha,verify_freeze,spec,class_name
from summarize_texture_analysis import fixed_variant

def close(x,y):np.testing.assert_allclose(x,y,atol=1e-11,rtol=1e-11,equal_nan=True)

def main():
    verify_freeze();checks={}
    previous=json.loads((OUT/'estado_anterior.json').read_text())['manifesto_antes_extensao']['saidas_sha256']
    protected={p:h for p,h in previous.items() if (p.startswith(('results/','maps/')) and Path(p).suffix in ['.csv','.npz','.npy']) or p.startswith('figures/')}
    for p,h in protected.items():assert sha(ROOT/p)==h,p
    checks['tabelas_mapas_e_figuras_etapa1_inalterados']=len(protected)
    for p,h in json.loads((OUT/'fontes_sha256.json').read_text()).items():assert sha(ROOT.parent/p)==h,p
    checks['fontes_da_extensao_inalteradas']=True
    execution=json.loads((OUT/'execucao.json').read_text());assert execution['status']=='PASS'
    checkpoint=json.loads((OUT/'checkpoint_pausa.json').read_text())
    checkpoint_checked=json.loads((OUT/'verificacao_checkpoint_retomada.json').read_text())
    assert checkpoint_checked['status']=='PASS' and checkpoint_checked['arquivos_conferidos']==46
    preserved=0
    for p,h in checkpoint['sha256'].items():
        if (p.startswith('results/extensao_churn_anular/') and p!='results/extensao_churn_anular/execucao.json') or p.startswith('maps/extensao_churn_anular/'):
            assert sha(ROOT/p)==h,p;preserved+=1
    assert sha(ROOT/'RETOMADA_CHECKPOINT.md')==checkpoint['sha256']['RETOMADA.md']
    assert sha(OUT/'execucao_checkpoint_pausa.json')==checkpoint['sha256']['results/extensao_churn_anular/execucao.json']
    assert sha(ROOT/'scripts/resume_texture_analysis.py')==execution['retomada_script_sha256']
    close(execution['segundos_total'],sum(v['segundos'] for v in execution['etapas'].values()))
    checks['checkpoint_retomada']={'conferidos_antes_retomada':46,'resultados_e_mapas_preservados':preserved,'historico_pausa_preservado':True}
    for p,h in execution['scripts_sha256'].items():assert sha(ROOT/'scripts'/p)==h,p
    assert json.loads((OUT/'validacao_numerica.json').read_text())['status']=='PASS'
    checks['codigo_executado_e_testes_conferidos']=True
    cat=pd.read_csv(OUT/'catalogo_ancoras.csv');chosen=cat[cat.selecionada].copy()
    assert len(chosen)==152 and chosen.aquisicao.nunique()==19
    assert chosen.groupby(['aquisicao','tipo']).size().eq(4).all()
    inv,signals=load_signals();lengths=inv.set_index('run').n
    pixels=0
    for r in chosen.itertuples():
        assert r.segmento_inicio==r.ancora-1024 and r.segmento_fim==r.ancora+1024
        assert r.segmento_inicio>=16 and r.segmento_fim<=lengths[r.aquisicao]-16
        assert 16<=r.largura<=512 and r.prominencia>=r.limiar_prominencia
        if not r.ancora_em_zona_mistura:
            p=np.asarray(Image.open(DATA/r.aquisicao/'__panorama_v11.png').convert('L'))
            f=np.asarray(Image.open(DATA/r.aquisicao/r.quadro_ancora).convert('L'))
            assert np.array_equal(p[:,r.ancora],f[:,r.coluna_quadro_ancora]);pixels+=1
    checks['ancoras_contexto_e_procedencia']={'total':len(chosen),'pixels_identicos_fora_mistura':pixels,'em_mistura':int(chosen.ancora_em_zona_mistura.sum())}
    import catalog_texture_anchors as catalog
    original_out=catalog.OUT
    with tempfile.TemporaryDirectory(prefix='wco_textura_verificacao_') as tmp:
        catalog.OUT=Path(tmp)
        try:
            with contextlib.redirect_stdout(io.StringIO()):catalog.main(render=False)
            for name in ['catalogo_ancoras.csv','pares_candidatos.csv','pares_ancoras.csv','cobertura_ancoras.csv','fontes_sha256.json']:
                assert sha(Path(tmp)/name)==sha(OUT/name),name
        finally:catalog.OUT=original_out
    checks['catalogo_e_pares_reproduzidos_byte_a_byte']=True
    pairs=pd.read_csv(OUT/'pares_ancoras.csv');assert len(pairs)==830 and not pairs.pair_id.duplicated().any()
    chosen=chosen.set_index('anchor_id')
    for r in pairs.itertuples():
        a,b=chosen.loc[r.anchor_x],chosen.loc[r.anchor_y]
        assert a.aquisicao==r.run_x and b.aquisicao==r.run_y and r.run_x<r.run_y
        ratio=max(a.largura,b.largura)/min(a.largura,b.largura);close(r.razao_larguras,ratio)
        if r.grupo=='minimos_compativeis':assert a.tipo==b.tipo=='minimo' and ratio<=1.5
        elif r.grupo=='minimos_distintos':assert a.tipo==b.tipo=='minimo' and ratio>=2
        else:assert r.grupo=='polaridades_opostas' and a.tipo!=b.tipo
    ev=pd.read_csv(OUT/'resultados_ancoras.csv');fixed=pd.read_csv(OUT/'resultados_fixos.csv')
    assert len(ev)==6640 and len(fixed)==2052
    assert not ev.duplicated(['pair_id','variante']).any()
    for v,g in ev.groupby('variante'):
        assert set(g.pair_id)==set(pairs.pair_id)
        for r in g.itertuples():
            a=spec(chosen.loc[r.anchor_x],v,'x');b=spec(chosen.loc[r.anchor_y],v,'y')
            absent=a is None or b is None
            assert (r.status=='nao_estimavel')==absent
            if absent:assert pd.isna(r.escore)
            else:assert r.centro_x==a[1]+1024 and r.centro_y==b[1]+1024
    ok=ev[ev.status=='calculado'];assert ok.escore.notna().all()
    for d in [ok,fixed]:
        assert d.escalas_validas.eq(200).all() and d.escalas_invalidas.eq(0).all()
        assert d.lags_principais.eq(65).all() and d.suporte_geometrico_fracao.eq(1).all()
        assert d.erro_equacao.max()<1e-11 and d.minimo_integral.min()>=0
        assert d[d.variante!='v5_mesmos_segmentos'].maximo_integral.max()<=1+1e-8
    legacy=json.loads((OUT/'auditoria_ultrapassagens_v5.json').read_text())
    assert legacy['status']=='PASS' and legacy['erro_maximo']<1e-9
    unexpected=[]
    for frame,kind in [(ok,'ancoras'),(fixed,'fixos')]:
        unexpected.extend((kind,r.pair_id) for r in frame[frame.maximo_integral>1+1e-8].itertuples())
    assert set(unexpected)=={(r['origem'],r['pair_id']) for r in legacy['casos_detalhados']}
    assert all(r['celulas_acima_1_lags_principais']==0 for r in legacy['casos_detalhados'])
    checks['ultrapassagens_legado_auditadas']=legacy
    original=pd.read_csv(ROOT/'results/resultados_fixos.csv').set_index(['pair_id','variante'])
    for r in fixed[fixed.origem=='etapa1_reutilizada_sem_alteracao'].itertuples():
        close(r.escore,original.loc[(r.pair_id,'principal')].escore)
    diag=pd.read_csv(OUT/'controles_mesma_aquisicao.csv')
    assert diag.tipo.value_counts().to_dict()=={'self':152,'mesma_aquisicao':19}
    assert diag.loc[diag.tipo=='self','erro_self_zero'].max()<1e-12
    checks['comparacoes']={'ancoras_calculadas':len(ok),'ancoras_ausentes':len(ev)-len(ok),'fixos':len(fixed),'diagnosticos':len(diag)}
    # Uma comparação experimental recalculada de forma direta em cada variante.
    errors=[]
    for variant in CFG['variantes']:
        r=ok[ok.variante==variant].sort_values('pair_id').iloc[0]
        sx=spec(chosen.loc[r.anchor_x],variant,'x');sy=spec(chosen.loc[r.anchor_y],variant,'y')
        x=crop_prepare(signals[sx[0]][sx[2]],sx[1],linear=sx[3]);y=crop_prepare(signals[sy[0]][sy[2]],sy[1],linear=sy[3])
        m,s=evaluate(x,y,engine=sx[-1]);del m
        for k,v in s.items():close(v,r[k])
        errors.append(abs(s['escore']-r.escore))
    checks['recalculo_oito_variantes_erro_maximo']=max(errors)
    with np.load(OUT/'perfis_ancoras.npz') as z:
        assert z['mean_lag'].shape==(830,200)
        close(z['mean_lag'].mean(axis=1),ok[ok.variante=='principal'].set_index('pair_id').loc[z['pair_ids'],'escore'])
    nmap=0
    for p in MAPS.glob('*.npz'):
        with np.load(p) as z:
            assert z['wco'].shape==(200,4095) and z['x'].shape==z['y'].shape==(2048,)
            valid=z['denominator_valid'];den=abs(z['wxx'][:,2047])*abs(z['wyy'][:,2047])
            close(z['wco'][valid],abs(z['wxy'][valid])**2/den[valid,None])
            meta=json.loads(str(z['meta_json']));score=np.mean(z['wco'][valid][:,abs(z['lags'])<=32])
            for axis,key in [('x','anchor_x'),('y','anchor_y')]:
                anchor=chosen.loc[meta[key]]
                expected=crop_prepare(signals[anchor.aquisicao]['principal'],int(anchor.segmento_inicio))
                close(z[axis],expected)
            close(score,meta['escore']);close(score,ok[(ok.pair_id==meta['pair_id'])&(ok.variante=='principal')].escore.iloc[0]);nmap+=1
    checks['mapas_integrais_conferidos']=nmap
    from audit_example_scores import main as verify_examples
    example_review=verify_examples()
    checks['extremos_dos_exemplos_conferidos']={
        'exemplos':len(example_review['exemplos']),
        'universos':example_review['universos'],
        'erro_maximo_recalculo_mapa':max(e['erro_maximo_mapa'] for e in example_review['exemplos'])}
    for raw,name in [(ok,'agregados_ancoras.csv'),(fixed,'agregados_fixos.csv')]:
        indexed=raw.groupby(['variante','grupo','run_x','run_y'])
        for r in pd.read_csv(OUT/name).itertuples():
            g=indexed.get_group((r.variante,r.grupo,r.run_x,r.run_y));close(r.escore,g.escore.mean());assert r.n_segmentos==len(g)
    ag=pd.read_csv(OUT/'agregados_ancoras.csv').set_index(['variante','grupo','run_x','run_y'])
    af=pd.read_csv(OUT/'agregados_fixos.csv').set_index(['variante','run_x','run_y'])
    cp=pd.read_csv(OUT/'contrastes_pareados.csv')
    for r in cp.itertuples():
        a=ag.loc[(r.variante,'minimos_compativeis',r.run_x,r.run_y)]
        b=af.loc[(fixed_variant(r.variante),r.run_x,r.run_y)] if r.comparador=='fixos' else ag.loc[(r.variante,r.comparador,r.run_x,r.run_y)]
        close([r.ancoras,r.controle,r.diferenca,r.diferenca_conservador],[a.escore,b.escore,a.escore-b.escore,a.conservador-b.conservador])
    cs=pd.read_csv(OUT/'contrastes_por_padrao.csv')
    for r in cs.itertuples():
        g=cp[(cp.variante==r.variante)&(cp.classe==r.classe)&(cp.comparador==r.comparador)]
        close(r.diferenca,g.diferenca.mean());assert r.pares_aquisicoes==len(g)
    for r in pd.read_csv(OUT/'retirada_aquisicao.csv').itertuples():
        g=cp[(cp.variante=='principal')&(cp.classe==r.classe)&(cp.comparador==r.comparador)&(cp.run_x!=r.retirada)&(cp.run_y!=r.retirada)]
        close(r.diferenca,g.diferenca.mean());assert r.pares_restantes==len(g)
    # Recalcula as sensibilidades diretamente dos pares locais comuns,
    # incluindo perdas de correspondência em grupos de controle.
    def direct_contrast(data,variant,cl,co,ids):
        d=data[(data.variante==variant)&data.pair_id.isin(ids)].copy()
        same=d.pattern_x.eq(d.pattern_y)
        d=d[~same] if cl=='entre' else d[same & d.pattern_x.eq(cl)]
        eq=d[d.grupo=='minimos_compativeis'].groupby(['run_x','run_y']).escore.mean()
        control=fixed[fixed.variante==fixed_variant(variant)] if co=='fixos' else d[d.grupo==co]
        ref=control.groupby(['run_x','run_y']).escore.mean()
        index=eq.index.intersection(ref.index)
        return eq.loc[index]-ref.loc[index]
    base_ids=set(ok[ok.variante=='principal'].pair_id)
    for r in pd.read_csv(OUT/'sensibilidade_conjunto_comum.csv').itertuples():
        ids=base_ids & set(ok[ok.variante==r.variante].pair_id)
        b=direct_contrast(ok,'principal',r.classe,r.comparador,ids)
        a=direct_contrast(ok,r.variante,r.classe,r.comparador,ids)
        common=b.index.intersection(a.index);assert len(common)==r.pares_aquisicoes
        close([r.principal_comum,r.variante_comum,r.mudanca],
              [b.loc[common].mean(),a.loc[common].mean(),(a.loc[common]-b.loc[common]).mean()])
    checks['sensibilidades_mesmos_ids_conferidas']=True
    global_ids=set.intersection(*(set(ok[ok.variante==v].pair_id) for v in CFG['variantes']))
    assert len(global_ids)==815
    for r in pd.read_csv(OUT/'sensibilidade_conjunto_global.csv').itertuples():
        direct=direct_contrast(ok,r.variante,r.classe,r.comparador,global_ids)
        close(r.diferenca,direct.mean());assert r.pares_aquisicoes==len(direct)
        assert r.pares_locais_comuns_todas_classes==len(global_ids)
    checks['painel_sensibilidades_mesmos_815_pares_locais']=True
    checks['agregacoes_contrastes_e_retiradas_conferidos']=True
    for p in (ROOT/'figures/extensao_churn_anular').rglob('*.png'):
        with Image.open(p) as im:assert min(im.size)>=400;im.verify()
    pdf=ROOT/'RESUMO_EXECUTIVO.pdf';blob=pdf.read_bytes();assert blob.startswith(b'%PDF') and blob.rstrip().endswith(b'%%EOF')
    pages=len(re.findall(rb'/Type\s*/Page\b',blob));assert pages>=15
    text=subprocess.run(['pdftotext',str(pdf),'-'],check=True,capture_output=True).stdout.decode()
    for phrase in ['152','Churn','Anular','Limitações','2.048']:assert phrase in text,phrase
    checks['paginas_pdf_integrado']=pages
    checks.update(status='PASS',data_utc=datetime.now(timezone.utc).isoformat(),
        limite='Consistência e rastreabilidade computacionais; não valida identidade de ondas/bolhas, rótulos locais, independência estatística ou generalização.')
    (OUT/'verificacao_final.json').write_text(json.dumps(checks,indent=2,ensure_ascii=False)+'\n');print(json.dumps(checks,indent=2,ensure_ascii=False))

if __name__=='__main__':main()
