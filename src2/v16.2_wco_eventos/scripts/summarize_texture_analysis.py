"""Agregação por par de aquisições, contrastes pareados e perdas de cobertura."""
import json
import numpy as np
import pandas as pd
from texture_common import OUT,CFG,class_name

KEY=['variante','classe','grupo','run_x','run_y']

def aggregate(df):
    return df.groupby(KEY,as_index=False).agg(escore=('escore','mean'),conservador=('escore_conservador','mean'),n_segmentos=('pair_id','size'))

def fixed_variant(v):
    return v if v in ['principal','faixa30_49','somente_media','v5_mesmos_segmentos'] else 'principal'

def contrasts(events,fixed):
    rows=[]
    for (v,cl),g in events.groupby(['variante','classe']):
        eq=g[g.grupo=='minimos_compativeis']
        for comparator in ['fixos','minimos_distintos','polaridades_opostas']:
            control=(fixed[(fixed.variante==fixed_variant(v))&(fixed.classe==cl)] if comparator=='fixos'
                     else g[g.grupo==comparator])
            common=eq.merge(control,on=['run_x','run_y'],suffixes=('_a','_c'),validate='one_to_one')
            for r in common.itertuples():
                rows.append(dict(variante=v,classe=cl,comparador=comparator,run_x=r.run_x,run_y=r.run_y,
                    ancoras=r.escore_a,controle=r.escore_c,diferenca=r.escore_a-r.escore_c,
                    diferenca_conservador=r.conservador_a-r.conservador_c,n_ancoras=r.n_segmentos_a,n_controle=r.n_segmentos_c))
    return pd.DataFrame(rows)

def pattern_summary(df):
    return df.groupby(['variante','classe','comparador'],as_index=False).agg(
        ancoras=('ancoras','mean'),controle=('controle','mean'),diferenca=('diferenca','mean'),
        diferenca_conservador=('diferenca_conservador','mean'),pares_aquisicoes=('diferenca','size'),
        fracao_positiva=('diferenca',lambda x:float((x>0).mean())),minimo=('diferenca','min'),maximo=('diferenca','max'))

def main():
    execution=json.loads((OUT/'execucao.json').read_text());assert execution['status']=='PASS'
    all_ev=pd.read_csv(OUT/'resultados_ancoras.csv');ev=all_ev[all_ev.status=='calculado'].copy()
    fixed=pd.read_csv(OUT/'resultados_fixos.csv')
    for df in [ev,fixed]:df['classe']=df.apply(class_name,axis=1)
    ag=aggregate(ev);af=aggregate(fixed)
    ag.to_csv(OUT/'agregados_ancoras.csv',index=False);af.to_csv(OUT/'agregados_fixos.csv',index=False)
    cp=contrasts(ag,af);cs=pattern_summary(cp)
    cp.to_csv(OUT/'contrastes_pareados.csv',index=False);cs.to_csv(OUT/'contrastes_por_padrao.csv',index=False)
    bal=cs[cs.classe.isin(CFG['padroes'])].groupby(['variante','comparador'],as_index=False).agg(
        diferenca=('diferenca','mean'),padroes=('classe','nunique'))
    bal.to_csv(OUT/'contrastes_balanceados.csv',index=False)
    # Comparação de variantes usa os mesmos IDs de pares de âncoras e aquisições.
    common_rows=[];base=ev[ev.variante=='principal']
    for variant in CFG['variantes'][1:]:
        alt=ev[ev.variante==variant];ids=set(alt.pair_id)&set(base.pair_id)
        b=base[base.pair_id.isin(ids)].copy();a=alt[alt.pair_id.isin(ids)].copy()
        cb=contrasts(aggregate(b),af);ca=contrasts(aggregate(a),af)
        merged=cb.merge(ca,on=['classe','comparador','run_x','run_y'],suffixes=('_base','_alt'),validate='one_to_one')
        for (cl,co),g in merged.groupby(['classe','comparador']):
            common_rows.append(dict(variante=variant,classe=cl,comparador=co,pares_aquisicoes=len(g),
                pares_ancoras_disponiveis_classe=int((a.classe==cl).sum()),
                principal_comum=float(g.diferenca_base.mean()),variante_comum=float(g.diferenca_alt.mean()),
                mudanca=float((g.diferenca_alt-g.diferenca_base).mean())))
    pd.DataFrame(common_rows).to_csv(OUT/'sensibilidade_conjunto_comum.csv',index=False)
    global_ids=set.intersection(*(set(ev[ev.variante==v].pair_id) for v in CFG['variantes']))
    global_summary=pattern_summary(contrasts(aggregate(ev[ev.pair_id.isin(global_ids)]),af))
    global_summary['pares_locais_comuns_todas_classes']=len(global_ids)
    global_summary.to_csv(OUT/'sensibilidade_conjunto_global.csv',index=False)
    leave=[]
    for (cl,co),g in cp[cp.variante=='principal'].groupby(['classe','comparador']):
        for run in sorted(set(g.run_x)|set(g.run_y)):
            rest=g[(g.run_x!=run)&(g.run_y!=run)]
            leave.append(dict(classe=cl,comparador=co,retirada=run,pares_restantes=len(rest),
                diferenca=float(rest.diferenca.mean()),diferenca_completa=float(g.diferenca.mean()),
                fracao_positiva=float((rest.diferenca>0).mean()) if len(rest) else np.nan))
    pd.DataFrame(leave).to_csv(OUT/'retirada_aquisicao.csv',index=False)
    cat=pd.read_csv(OUT/'catalogo_ancoras.csv').query('selecionada');stability=[]
    for (pattern,kind),g in cat.groupby(['padrao','tipo']):
        for alt in ['sigma12','sigma20','faixa30_49']:
            d=g['deslocamento_'+alt].abs().dropna()
            stability.append(dict(padrao=pattern,tipo=kind,alternativa=alt,total=len(g),disponiveis=len(d),ausentes=len(g)-len(d),
                mediana_abs_px=float(d.median()),p90_abs_px=float(d.quantile(.9)),max_abs_px=float(d.max()),
                fracao_ate16_total=float((d<=16).sum()/len(g))))
    pd.DataFrame(stability).to_csv(OUT/'estabilidade_ancoras.csv',index=False)
    selected=pd.read_csv(OUT/'pares_ancoras.csv');selected['classe']=selected.apply(class_name,axis=1)
    coverage=[]
    for (cl,group),g in selected.groupby(['classe','grupo']):
        expected={'churn':45,'wavy annular':36,'entre':90}[cl]
        coverage.append(dict(classe=cl,grupo=group,pares_ancoras=len(g),pares_aquisicoes=len(g[['run_x','run_y']].drop_duplicates()),
            pares_aquisicoes_possiveis=expected,aquisicoes=len(set(g.run_x)|set(g.run_y))))
    pd.DataFrame(coverage).to_csv(OUT/'cobertura_grupos.csv',index=False)
    examples=[]
    for cl in CFG['padroes']:
        group=base[(base.classe==cl)&(base.grupo=='minimos_compativeis')]
        for name,row in [('predeterminado',group.sort_values('pair_id').iloc[0]),
                         ('maior_escore',group.sort_values(['escore','pair_id']).iloc[-1]),
                         ('menor_escore',group.sort_values(['escore','pair_id']).iloc[0])]:
            examples.append(dict(classe=cl,tipo=name,pair_id=row.pair_id,escore=row.escore,
                escolha='primeiro_ID_sem_WCO' if name=='predeterminado' else 'ilustracao_pos_resultado_sem_mudar_catalogo'))
    pd.DataFrame(examples).to_csv(OUT/'exemplos.csv',index=False)
    summary=dict(ancoras=len(cat),minimos=int((cat.tipo=='minimo').sum()),maximos=int((cat.tipo=='maximo').sum()),
        aquisicoes=cat.aquisicao.nunique(),pares_ancoras=len(selected),comparacoes_planejadas=len(all_ev),
        comparacoes_calculadas=len(ev),comparacoes_nao_estimaveis=len(all_ev)-len(ev),
        comparacoes_fixos=len(fixed),ancoras_em_mistura=int(cat.ancora_em_zona_mistura.sum()),
        ancoras_com_sobreposicao=int((cat.segmentos_sobrepostos>0).sum()),
        ancoras_com_outros_marcos=int((cat.vizinhos_selecionados_no_contexto>0).sum()),
        minimos_com_sobreposicao=int(((cat.tipo=='minimo')&(cat.segmentos_sobrepostos>0)).sum()),
        minimos_com_outros_marcos=int(((cat.tipo=='minimo')&(cat.vizinhos_selecionados_no_contexto>0)).sum()),
        minimos_em_mistura=int(((cat.tipo=='minimo')&cat.ancora_em_zona_mistura).sum()),
        contextos_com_deslocamentos_nao_positivos=int((cat.deslocamentos_nao_positivos_contexto>0).sum()),
        costuras_contexto_min=int(cat.costuras_contexto.min()),costuras_contexto_max=int(cat.costuras_contexto.max()),
        invalidas_denominador=int(ev.escalas_invalidas.sum()),erro_equacao_maximo=float(ev.erro_equacao.max()),
        maximo_wco=float(ev.maximo_integral.max()),minimo_wco=float(ev.minimo_integral.min()),
        maximo_wco_morlet=float(ev[ev.variante!='v5_mesmos_segmentos'].maximo_integral.max()),
        maximo_wco_v5=float(ev[ev.variante=='v5_mesmos_segmentos'].maximo_integral.max()),
        mapas_ancorados_v5_acima_1=int(((ev.variante=='v5_mesmos_segmentos')&(ev.maximo_integral>1+1e-8)).sum()),
        mapas_fixos_v5_acima_1=int(((fixed.variante=='v5_mesmos_segmentos')&(fixed.maximo_integral>1+1e-8)).sum()),
        ausencias_por_variante=all_ev[all_ev.status!='calculado'].groupby('variante').size().to_dict(),
        criterios_descritivos=dict(margem=CFG['efeito_descritivo_minimo'],cobertura=CFG['cobertura_minima_fracao']),
        status='PASS')
    (OUT/'resumo.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(summary,indent=2,ensure_ascii=False));print(cs[cs.variante=='principal'].to_string(index=False))

if __name__=='__main__':main()
