"""Agrega por par de aquisições, preservando cobertura e dependência."""
import json
import numpy as np
import pandas as pd
from analysis_common import ROOT,BASE

def aggregate(df):
    keys=['variante','grupo','run_x','run_y','pattern_x','pattern_y']
    out=df.groupby(keys,dropna=False).agg(escore=('escore','mean'),n_segmentos=('escore','size'),
        conservador=('escore_conservador','mean'),zero=('media_lag_zero','mean'),
        maximo=('maximo_integral','mean')).reset_index()
    out['padrao_combinacao']=[' × '.join(sorted({x,y})) for x,y in zip(out.pattern_x,out.pattern_y)]
    return out

def main():
    out=ROOT/'results'
    e=pd.read_csv(out/'resultados_eventos.csv');f=pd.read_csv(out/'resultados_fixos.csv')
    assert len(e)==165*8 and len(f)==975+135
    a=aggregate(e);b=aggregate(f)
    a.to_csv(out/'agregados_eventos.csv',index=False);b.to_csv(out/'agregados_fixos.csv',index=False)
    rows=[]
    for var in sorted(e.variante.unique()):
        refvar='principal' if var.startswith(('bordas_','ancoras_')) else var
        for scope in ['intra','entre']:
            eq=a[(a.variante==var)&(a.grupo==scope+'_equivalentes')]
            ne=a[(a.variante==var)&(a.grupo==scope+'_distintos')]
            fi=b[(b.variante==refvar)&(b.grupo==scope+'_fixo')]
            for name,comp in [('distintos',ne),('fixos',fi)]:
                merged=eq.merge(comp,on=['run_x','run_y'],suffixes=('_eq','_ref'),validate='one_to_one')
                for r in merged.itertuples(index=False):
                    rows.append(dict(variante=var,escopo=scope,comparador=name,
                        run_x=r.run_x,run_y=r.run_y,padrao=r.padrao_combinacao_eq,
                        equivalentes=r.escore_eq,controle=r.escore_ref,
                        diferenca=r.escore_eq-r.escore_ref,
                        diferenca_conservador=r.conservador_eq-r.conservador_ref,
                        n_eventos_eq=r.n_segmentos_eq,n_controles=r.n_segmentos_ref))
    contrasts=pd.DataFrame(rows);contrasts.to_csv(out/'contrastes_pareados.csv',index=False)
    by=contrasts.groupby(['variante','escopo','comparador','padrao']).agg(
        equivalentes=('equivalentes','mean'),controle=('controle','mean'),diferenca=('diferenca','mean'),
        diferenca_conservador=('diferenca_conservador','mean'),
        pares_aquisicoes=('diferenca','size'),minimo_par=('diferenca','min'),maximo_par=('diferenca','max')).reset_index()
    by.to_csv(out/'contrastes_por_padrao.csv',index=False)
    balance=by.groupby(['variante','escopo','comparador']).agg(
        equivalentes=('equivalentes','mean'),controle=('controle','mean'),diferenca=('diferenca','mean'),
        padroes=('padrao','size'),pares_aquisicoes=('pares_aquisicoes','sum')).reset_index()
    balance.to_csv(out/'contrastes_balanceados.csv',index=False)
    primary=contrasts[contrasts.variante=='principal'];loo=[]
    for (scope,comp,pattern),g in primary.groupby(['escopo','comparador','padrao']):
        runs=sorted(set(g.run_x)|set(g.run_y))
        for run in runs:
            sub=g[(g.run_x!=run)&(g.run_y!=run)]
            loo.append(dict(escopo=scope,comparador=comp,padrao=pattern,retirada=run,
                            pares_restantes=len(sub),diferenca=sub.diferenca.mean()))
    pd.DataFrame(loo).to_csv(out/'retirada_aquisicao.csv',index=False)
    cat=pd.read_csv(out/'catalogo_eventos.csv');pairs=pd.read_csv(out/'pares_eventos.csv')
    inv=pd.read_csv(out/'auditoria_entradas.csv').query('incluida');cover=[]
    for pat,g in inv.groupby('padrao'):
        eligible=cat[(cat.padrao==pat)&cat.incluida]
        eq=pairs[(pairs.pattern_x==pat)&(pairs.grupo=='intra_equivalentes')]
        eqevents=set(eq.event_x)|set(eq.event_y);eqruns=set(eq.run_x)|set(eq.run_y)
        cand=cat[cat.padrao==pat]
        cover.append(dict(padrao=pat,aquisicoes=len(g),candidatos=len(cand),eventos_aceitos=len(eligible),
            aquisicoes_com_evento=eligible.aquisicao.nunique(),eventos_em_equivalentes_intra=len(eqevents),
            aquisicoes_em_equivalentes_intra=len(eqruns),pares_eventos_equivalentes_intra=len(eq),
            pares_aquisicoes_equivalentes_intra=len(eq[['run_x','run_y']].drop_duplicates()),
            fracao_aquisicoes_equivalentes_intra=len(eqruns)/len(g)))
    coverage=pd.DataFrame(cover);coverage.to_csv(out/'cobertura.csv',index=False)
    # Peso por evento e por par de aquisições são intencionalmente separados.
    dist= a[a.variante=='principal'].groupby(['grupo','padrao_combinacao']).agg(
        media=('escore','mean'),mediana=('escore','median'),minimo=('escore','min'),maximo=('escore','max'),
        pares_aquisicoes=('escore','size'),pares_eventos=('n_segmentos','sum')).reset_index()
    dist.to_csv(out/'distribuicoes.csv',index=False)
    base=json.loads((BASE/'results/preprocessing_conclusion.json').read_text())
    (out/'referencia_historica.json').write_text(json.dumps(dict(
        origem='../v16_spatial2_18_71/results/preprocessing_conclusion.json',
        dominio='50–3.000 px; todos os lags; máximos por escala ponderados por WPS da referência 18–71',
        comparacao_direta_com_media_eventos_valida=False,dados=base),indent=2,ensure_ascii=False)+'\n')
    # Todos os exemplos extremos são somente ilustrações posteriores, sem alterar o catálogo.
    p=e[(e.variante=='principal')&(e.grupo=='intra_equivalentes')]
    best=p.loc[p.escore.idxmax()];worst=p.loc[p.escore.idxmin()]
    examples={'favoravel':best.pair_id,'desfavoravel':worst.pair_id,
              'regra':'Maior e menor escore intrapadrão, somente após o lote; não são exemplos representativos.'}
    (out/'exemplos.json').write_text(json.dumps(examples,ensure_ascii=False,indent=2)+'\n')
    accepted=cat[cat.incluida]
    conclusion=dict(
        status='PASS',eventos_aceitos=len(accepted),pares_eventos=len(pairs),
        aquisicoes_elegiveis=26,aquisicoes_com_eventos=int(accepted.aquisicao.nunique()),
        aquisicoes_com_equivalentes_intra=int(coverage.aquisicoes_em_equivalentes_intra.sum()),
        principais=by[by.variante=='principal'].to_dict('records'),
        balanceado=balance[balance.variante=='principal'].to_dict('records'),
        eventos_sem_parceiro=int((accepted.parceiros_equivalentes_outra_aquisicao==0).sum()),
        eventos_com_vizinhos_aceitos=int((accepted.outros_eventos_no_contexto>0).sum()),
        max_mudanca_fronteira=float(accepted[['mudanca_fronteira_limiar55','mudanca_fronteira_limiar65']].max().max()),
        escalas_invalidas_morlet=int(e[e.variante!='v5_mesmos_segmentos'].escalas_invalidas.sum()+f[f.variante!='v5_mesmos_segmentos'].escalas_invalidas.sum()),
        erro_equacao_maximo=float(max(e.erro_equacao.max(),f.erro_equacao.max())),
        maior_wco_morlet=float(max(e[e.variante!='v5_mesmos_segmentos'].maximo_integral.max(),f[f.variante!='v5_mesmos_segmentos'].maximo_integral.max())),
        maior_wco_v5=float(max(e[e.variante=='v5_mesmos_segmentos'].maximo_integral.max(),f[f.variante=='v5_mesmos_segmentos'].maximo_integral.max())),
        slug_sensibilidade_faixa30_49=float(by.query('variante=="faixa30_49" and escopo=="intra" and comparador=="distintos" and padrao=="slug"').diferenca.iloc[0]),
        frase='Não conseguimos melhorar de forma convincente a alternative WCO com este protocolo. Plug teve melhora marginal; slug dependeu de um único par e o contraste se inverteu na faixa 30–49.')
    (out/'conclusao.json').write_text(json.dumps(conclusion,ensure_ascii=False,indent=2)+'\n')
    print(by[by.variante=='principal'].to_string(index=False));print(coverage.to_string(index=False))

if __name__=='__main__':main()
