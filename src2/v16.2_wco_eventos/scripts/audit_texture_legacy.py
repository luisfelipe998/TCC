"""Confere ultrapassagens unitárias do núcleo discreto v5 sem clipping."""
import json
from datetime import datetime,timezone
import numpy as np
import pandas as pd
import pywt
from analysis_common import ROOT,load_signals,crop_prepare,evaluate,LAM
from texture_common import OUT,verify_freeze
from wavelet_core import LEGACY_WAVELET

def main():
    verify_freeze();_,signals=load_signals()
    cat=pd.read_csv(OUT/'catalogo_ancoras.csv').set_index('anchor_id')
    segments=pd.read_csv(ROOT/'results/segmentos_fixos.csv').set_index('segment_id')
    cases=[];rows=[]
    for name,kind in [('resultados_ancoras.csv','ancoras'),('resultados_fixos.csv','fixos')]:
        d=pd.read_csv(OUT/name).query("status=='calculado'")
        assert d[d.variante!='v5_mesmos_segmentos'].maximo_integral.max()<=1+1e-8
        for _,r in d[d.maximo_integral>1+1e-8].iterrows():
            assert r.variante=='v5_mesmos_segmentos';cases.append((kind,r))
    for kind,r in cases:
        if kind=='ancoras':
            a,b=cat.loc[r.anchor_x],cat.loc[r.anchor_y];sx,sy=int(a.segmento_inicio),int(b.segmento_inicio)
        else:
            ix,iy=r.pair_id.split('__');sx,sy=int(segments.loc[ix].inicio),int(segments.loc[iy].inicio)
        x=crop_prepare(signals[r.run_x]['principal'],sx);y=crop_prepare(signals[r.run_y]['principal'],sy)
        m,s=evaluate(x,y,engine='legacy')
        np.testing.assert_allclose([s['escore'],s['maximo_integral']],[r.escore,r.maximo_integral],atol=1e-11,rtol=1e-11)
        mask=m['wco']>1+1e-8;js=np.flatnonzero(mask.any(axis=1));scales=m['scales'][js]
        # Caminho independente: somas diretas + CWT por convolução temporal,
        # preservando a convenção cruzada efetiva do código legado e refletindo o lag.
        rx=np.correlate(x,x,'full')/len(x);ry=np.correlate(y,y,'full')/len(y)
        rxy=np.correlate(y,x,'full')/len(x)
        wx=pywt.cwt(rx,scales,LEGACY_WAVELET,method='conv',precision=12)[0]
        wy=pywt.cwt(ry,scales,LEGACY_WAVELET,method='conv',precision=12)[0]
        wxy=pywt.cwt(rxy,scales,LEGACY_WAVELET,method='conv',precision=12)[0][:,::-1]
        co=abs(wxy)**2/(abs(wx[:,2047])*abs(wy[:,2047]))[:,None]
        error=float(np.max(abs(co-m['wco'][js])));assert error<1e-9
        j,k=np.unravel_index(np.nanargmax(m['wco']),m['wco'].shape)
        rows.append(dict(origem=kind,pair_id=r.pair_id,variante=r.variante,escore=r.escore,
            maximo_fft=float(m['wco'][j,k]),maximo_convolucao=float(co.max()),
            lambda_pico_px=float(LAM[j]),lag_pico_px=int(m['lags'][k]),
            celulas_acima_1=int(mask.sum()),celulas_acima_1_lags_principais=int(mask[:,abs(m['lags'])<=32].sum()),
            escalas_com_ultrapassagem=int(len(js)),erro_fft_vs_convolucao=error))
    pd.DataFrame(rows).to_csv(OUT/'auditoria_ultrapassagens_v5.csv',index=False)
    report=dict(status='PASS',data_utc=datetime.now(timezone.utc).isoformat(),casos=len(rows),
        maior_valor=max(r['maximo_fft'] for r in rows),erro_maximo=max(r['erro_fft_vs_convolucao'] for r in rows),
        interpretacao='Ultrapassagens do núcleo discreto legado reproduzidas por somas diretas e convolução temporal. Não são erro de FFT nem motivo para clipping. O limite unitário é verificado separadamente na Morlet principal; a variante v5 é sensibilidade numérica.',
        casos_detalhados=rows)
    (OUT/'auditoria_ultrapassagens_v5.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(report,indent=2,ensure_ascii=False))

if __name__=='__main__':main()
