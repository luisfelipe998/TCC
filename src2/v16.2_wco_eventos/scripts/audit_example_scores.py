"""Confere extremos ilustrados, recortes e distinção entre intensidade e WCO."""
from datetime import datetime, timezone
import hashlib
import json
import numpy as np
import pandas as pd
from PIL import Image
from scipy.ndimage import gaussian_filter1d
from analysis_common import ROOT, DATA, LAM
from texture_common import OUT, MAPS
from wavelet_core import transform_correlation


def main():
    examples = pd.read_csv(OUT/'exemplos.csv')
    results = pd.read_csv(OUT/'resultados_ancoras.csv').query("variante=='principal'")
    anchors = pd.read_csv(OUT/'catalogo_ancoras.csv').query('selecionada').set_index('anchor_id')
    all_path = ROOT/'results/auditoria_adversarial/todos_pares_admissiveis.csv'
    admissible = pd.read_csv(all_path)
    records, universes, inputs = [], [], [OUT/'exemplos.csv', OUT/'resultados_ancoras.csv', OUT/'catalogo_ancoras.csv', all_path]
    for pattern in ['churn', 'wavy annular']:
        g = results[(results.pattern_x==pattern)&(results.pattern_y==pattern)&(results.grupo=='minimos_compativeis')]
        ordered = g.sort_values(['escore','pair_id'])
        expanded = admissible[(admissible.etapa=='extensao')&(admissible.pattern_x==pattern)&(admissible.grupo=='minimos_compativeis')].sort_values(['escore','pair_id'])
        universes.append(dict(padrao=pattern, pares_selecionados=len(g), pares_admissiveis=len(expanded),
            minimo_selecionados=float(ordered.escore.iloc[0]), maximo_selecionados=float(ordered.escore.iloc[-1]),
            minimo_admissiveis=float(expanded.escore.iloc[0]), maximo_admissiveis=float(expanded.escore.iloc[-1]),
            par_minimo_admissiveis=expanded.pair_id.iloc[0], par_maximo_admissiveis=expanded.pair_id.iloc[-1]))
        expected = {'menor_escore':ordered.iloc[0], 'maior_escore':ordered.iloc[-1],
                    'predeterminado':g.sort_values('pair_id').iloc[0]}
        for typ, row in expected.items():
            example = examples[(examples.classe==pattern)&(examples.tipo==typ)]
            assert len(example)==1 and example.iloc[0].pair_id==row.pair_id
            prepared, optical = [], []
            for key in ['anchor_x','anchor_y']:
                a = anchors.loc[row[key]]
                assert a.tipo=='minimo' and a.padrao==pattern
                path = DATA/a.aquisicao/'__panorama_v11.png'; inputs.append(path)
                image = np.asarray(Image.open(path).convert('L'))
                raw = image[18:72].sum(axis=0,dtype=np.int64)/54
                c, lo, hi = int(a.ancora), int(a.segmento_inicio), int(a.segmento_fim)
                assert (lo,hi)==(c-1024,c+1024)
                sm = gaussian_filter1d(raw,16,mode='reflect')
                assert sm[c]<=sm[c-1] and sm[c]<=sm[c+1]
                # Preparo separado: regressão por mínimos quadrados e desvio populacional.
                v = 255-raw[lo:hi]
                design = np.column_stack([np.arange(2048,dtype=float),np.ones(2048)])
                v = v-design@np.linalg.lstsq(design,v,rcond=None)[0]
                v -= v.mean(); v /= v.std(); prepared.append(v)
                optical.append(dict(anchor_id=a.name, centro=c, intensidade_suavizada=float(sm[c]), minimo_local=True))
            x,y = prepared
            wxx,_ = transform_correlation(np.correlate(x,x,'full')/2048,LAM)
            wyy,_ = transform_correlation(np.correlate(y,y,'full')/2048,LAM)
            wxy,_ = transform_correlation(np.correlate(x,y,'full')/2048,LAM)
            co = abs(wxy)**2/(abs(wxx[:,2047])*abs(wyy[:,2047]))[:,None]
            recalculated = float(co[:,2015:2080].mean())
            path = MAPS/f'{row.pair_id}.npz'; inputs.append(path)
            with np.load(path) as z:
                meta = json.loads(str(z['meta_json'])); assert meta['pair_id']==row.pair_id
                for key in ['anchor_x','anchor_y']: assert meta[key]==row[key]
                np.testing.assert_allclose(z['x'],x,atol=1e-11,rtol=1e-11)
                np.testing.assert_allclose(z['y'],y,atol=1e-11,rtol=1e-11)
                saved = float(z['wco'][:,abs(z['lags'])<=32].mean())
                map_error = float(np.max(abs(co-z['wco'])))
                assert z['denominator_valid'].all() and map_error<1e-10
            np.testing.assert_allclose([saved,recalculated,example.iloc[0].escore],row.escore,atol=1e-11,rtol=1e-11)
            records.append(dict(padrao=pattern,tipo=typ,pair_id=row.pair_id,
                escore_tabela=float(row.escore),escore_mapa=saved,escore_recalculado=recalculated,
                erro_maximo_mapa=map_error,ancoras=optical))
    report = dict(status='PASS',data_utc=datetime.now(timezone.utc).isoformat(),
        definicao='Média de WCO em 200 escalas logarítmicas (50–1024 px) e 65 lags (−32…+32).',
        metodo='Preparo por regressão e correlação direta independentes; CWT pelo núcleo principal preservado.',
        interpretacao='As figuras ordenam pares selecionados de mínimos compatíveis pela WCO do contexto, não pela escuridão da âncora. Os extremos não são necessariamente os de todos os pares admissíveis.',
        universos=universes, exemplos=records,
        entradas_sha256={str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(set(inputs))})
    (OUT/'verificacao_exemplos.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    return report


if __name__=='__main__':
    r=main()
    print(json.dumps({'status':r['status'],'universos':r['universos'],'exemplos_conferidos':len(r['exemplos']),
        'erro_maximo_mapa':max(e['erro_maximo_mapa'] for e in r['exemplos'])},indent=2,ensure_ascii=False))
