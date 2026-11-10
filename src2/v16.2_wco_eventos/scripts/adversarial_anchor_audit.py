"""Distingue correspondência local de repetibilidade da seleção completa."""
from datetime import datetime, timezone
import json
import numpy as np
import pandas as pd
from analysis_common import ROOT, load_signals
from catalog_texture_anchors import detect, select

OUT = ROOT/'results/auditoria_adversarial'
EXT = ROOT/'results/extensao_churn_anular'


def main():
    _, signals = load_signals()
    cat = pd.read_csv(EXT/'catalogo_ancoras.csv').query('selecionada')
    records = []
    for run, group in cat.groupby('aquisicao'):
        for name, band, sigma in [('sigma12','principal',12), ('sigma20','principal',20), ('faixa30_49','faixa30_49',16)]:
            candidates, _ = detect(signals[run][band], sigma)
            selected = candidates.loc[sorted(select(candidates))]
            for anchor in group.itertuples():
                g = candidates[(candidates.tipo==anchor.tipo)&candidates.elegivel_contexto]
                h = selected[selected.tipo==anchor.tipo]
                nearest_candidate = min(abs(g.ancora-anchor.ancora)) if len(g) else np.inf
                nearest_selected = min(abs(h.ancora-anchor.ancora)) if len(h) else np.inf
                records.append(dict(anchor_id=anchor.anchor_id, aquisicao=run, padrao=anchor.padrao,
                    tipo=anchor.tipo, alternativa=name, distancia_candidato=float(nearest_candidate),
                    distancia_selecionado=float(nearest_selected), candidato_em_64=nearest_candidate<=64,
                    selecionado_em_64=nearest_selected<=64, selecionado_em_16=nearest_selected<=16,
                    numero_selecionados_tipo=len(h)))
    df = pd.DataFrame(records)
    df.to_csv(OUT/'repetibilidade_selecao_detalhes.csv', index=False)
    summary = df.groupby(['padrao','tipo','alternativa'], as_index=False).agg(
        total=('anchor_id','size'), candidatos_em_64=('candidato_em_64','sum'),
        selecionados_em_64=('selecionado_em_64','sum'), selecionados_em_16=('selecionado_em_16','sum'),
        maior_distancia_ao_selecionado=('distancia_selecionado','max'))
    summary.to_csv(OUT/'repetibilidade_selecao.csv', index=False)
    coverage = []
    cp = pd.read_csv(EXT/'contrastes_pareados.csv').query("variante=='principal'")
    for (pattern, control), g in cp.groupby(['classe','comparador']):
        total = {'churn':10,'wavy annular':9,'entre':19}[pattern]
        nr = len(set(g.run_x)|set(g.run_y))
        pairs = {'churn':45,'wavy annular':36,'entre':90}[pattern]
        coverage.append(dict(padrao=pattern,controle=control,aquisicoes=nr,aquisicoes_possiveis=total,
            fracao_aquisicoes=nr/total,pares=len(g),pares_possiveis=pairs,fracao_pares=len(g)/pairs))
    pd.DataFrame(coverage).to_csv(OUT/'cobertura_denominadores.csv',index=False)
    candidates = pd.read_csv(EXT/'pares_candidatos.csv')
    available = candidates[(candidates.escopo=='intra')&candidates.grupo.isin(['minimos_compativeis','minimos_distintos'])]
    report = dict(status='PASS', data_utc=datetime.now(timezone.utc).isoformat(),
        testes='Reexecução do detector e da seleção por estrato; nenhuma nova WCO nem alteração do catálogo.',
        resumo=summary.to_dict('records'), cobertura=coverage,
        pares_intra_admissiveis_por_grupo=available.groupby(['pattern_x','grupo']).size().to_dict())
    report['pares_intra_admissiveis_por_grupo']={' / '.join(k):int(v) for k,v in report['pares_intra_admissiveis_por_grupo'].items()}
    (OUT/'ancoras.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    print(summary.to_string(index=False))
    print(report['pares_intra_admissiveis_por_grupo'])


if __name__=='__main__':main()
