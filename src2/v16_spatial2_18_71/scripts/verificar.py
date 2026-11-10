"""Verificação final de fontes, produtos, cobertura, apresentação e arquivo histórico."""
from itertools import combinations
from pathlib import Path
import json
import re
import sys
from urllib.parse import unquote
sys.dont_write_bytecode = True
from common import *
import wavelet_core as core


def main():
    acquisitions=records()
    computed=pairs()
    runs=[r["ensaio"] for r in acquisitions]
    expected={frozenset(pair) for pair in combinations(runs,2)}
    actual=[frozenset((r["run_x"],r["run_y"])) for r in computed]
    assert len(runs)==len(set(runs))==26 and len(actual)==len(set(actual))==325
    assert set(actual)==expected
    assert sum(r["tipo"]=="intra" for r in computed)==87
    assert sum(r["tipo"]=="inter" for r in computed)==238
    assert sha(ROOT/"scripts/wavelet_core.py")==sha(OLD/"scripts/wavelet_core.py")
    sources=read_csv(INTER/"auditoria/entradas.csv")
    assert len(sources)==31
    for row in sources:
        run=row["ensaio"]
        assert sha(DATA/run/"__panorama_v11.png")==row["panorama_sha256"]
        assert sha(DATA/run/"__dados_cropped_v11.xlsx")==row["planilha_sha256"]
        raw=load_raw(run)
        assert hashlib.sha256(raw.astype("<f8",copy=False).tobytes()).hexdigest()==row["raw_sha256"]
        np.testing.assert_array_equal(raw,load_image(run)[18:72].mean(axis=0))
    print("Fontes: 31 entradas e 62 arquivos de origem verificados.",flush=True)
    equation_error=0.
    expected_titles={"Sinais normalizados","Espectro de potência wavelet (WPS) · A",
                     "Espectro de potência wavelet (WPS) · B","Espectro cruzado wavelet (WCS) · módulo",
                     "Coerência wavelet (WCO)"}
    for i,row in enumerate(computed,1):
        a,b=row["run_x"],row["run_y"]
        stem=INTER/"pares"/f"{a}__{b}"
        assert sha(stem.with_suffix(".npz"))==row["npz_sha256"]
        assert sha(ROOT/row["arquivo"])==row["png_sha256"]
        assert row["assinaturas"]["core"]==sha(ROOT/"scripts/wavelet_core.py")
        assert row["assinaturas"]["analysis"]==sha(ROOT/"scripts/analise.py")
        assert row["figuras_sha256"]==sha(ROOT/"scripts/figuras.py")
        assert json.loads(stem.with_suffix(".json").read_text())==row
        assert row["erro_metricas_vs_original"]<1e-11
        layout=row["apresentacao"]
        assert expected_titles.issubset(layout["titulos_paineis"])
        assert "Extensão do panorama (px)" in layout["rotulos_x"]
        assert "Altura do\npanorama (px)" in layout["rotulos_y"]
        with np.load(stem.with_suffix(".npz")) as saved:
            np.testing.assert_array_equal(saved["wavelengths"],WAVES)
            n=row["comprimento"]
            lags=saved["lags"]
            assert lags[0]==-n+1 and lags[-1]==n-1 and np.all(np.diff(lags)>0)
            assert int(saved["full_lag_count"])==2*n-1
            assert len(np.flatnonzero(lags==0))==1
            j0=int(np.flatnonzero(lags==0)[0])
            np.testing.assert_allclose(saved["x"],core.prepare(load_raw(a)[:n]),atol=1e-13,rtol=1e-13)
            np.testing.assert_allclose(saved["y"],core.prepare(load_raw(b)[:n]),atol=1e-13,rtol=1e-13)
            zx,zy=saved["wps_zero_x"],saved["wps_zero_y"]
            np.testing.assert_allclose(zx,abs(saved["wxx"][:,j0]),atol=1e-13,rtol=1e-13)
            np.testing.assert_allclose(zy,abs(saved["wyy"][:,j0]),atol=1e-13,rtol=1e-13)
            wco=saved["wco"]
            reconstructed=abs(saved["wxy"])**2/(zx*zy)[:,None]
            error=float(abs(wco-reconstructed).max())
            assert error<1e-12 and np.isfinite(wco).all()
            equation_error=max(equation_error,error)
            assert wco.min()>=0 and wco.max()<=1+1e-8
            assert wco.max()<=row["wco_maxima"]+1e-12
            peaks=saved["scale_peak_wco"]
            assert abs(peaks.max()-row["wco_maxima"])<1e-12
            assert (wco.max(axis=1)<=peaks+1e-12).all()
            for cut in (.01,.05,.1):
                active=(zx>=cut*zx.max()) & (zy>=cut*zy.max())
                weights=np.where(active,np.sqrt(zx*zy),0)
                score=np.dot(weights,peaks)/weights.sum()
                assert abs(score-row[f"escore_piso{int(cut*100):02d}"])<1e-12
            for field,max_field in (("mean_x_all","max_wxx"),("mean_y_all","max_wyy")):
                assert saved[field].shape==(200,) and np.isfinite(saved[field]).all()
                assert np.min(saved[field])>=0 and np.max(saved[field])<=float(saved[max_field])+1e-12
        if i%50==0 or i==325:print(f"Mapas e escores: {i}/325",flush=True)

    for folder,count in (("01_panoramas",26),("02_wavelets",87),("03_wavelets_inter",238),("04_estatistica",4)):
        pngs=list((ROOT/folder).glob("*.png"))
        assert len(pngs)==count,(folder,len(pngs))
        for path in pngs:
            with Image.open(path) as image:image.verify()
    layouts=json.loads((INTER/"apresentacao/figuras.json").read_text())+json.loads((INTER/"apresentacao/estatistica.json").read_text())
    assert len(layouts)==355
    assert all(r["fonte_minima_inserida_pt"]>=8 and r["altura_inserida_cm"]<=24 for r in layouts)
    short_names=[run.removesuffix("_001") for run in runs]
    for name in ("matriz_pearson_todos.png","matriz_escore_wco.png"):
        layout=next(v for v in layouts if v['arquivo']==f'04_estatistica/{name}')
        assert layout['ensaios_x']==layout['ensaios_y']==short_names
        assert layout['rotulos_sem_sobreposicao']
    with np.load(INTER/"descritores.npz") as saved:
        assert saved["runs"].tolist()==runs
        matrix=saved["r"]
        np.testing.assert_allclose(matrix,np.corrcoef(saved["residuals"]),atol=1e-13,rtol=1e-13)
        np.testing.assert_array_equal(matrix,matrix.T)
        np.testing.assert_array_equal(matrix.diagonal(),np.ones(26))
        csv_rows=read_csv(TABLES/"matriz_pearson_todos.csv")
        assert [r["ensaio"] for r in csv_rows]==runs
        np.testing.assert_allclose(matrix,[[float(row[run]) for run in runs] for row in csv_rows],atol=1e-14)
        averages=read_csv(TABLES/"medias_por_ensaio.csv")
        assert len(averages)==26
        for i,row in enumerate(averages):
            assert row["ensaio"]==runs[i]
            assert abs(float(row["pearson_media_todos"])-(matrix[i].sum()-1)/25)<1e-13
    assert len(read_csv(TABLES/"pares_pearson_wco.csv"))==325
    assert len(read_csv(TABLES/"perfis_espectrais.csv"))==5200
    adversarial=json.loads((INTER/'auditoria/adversarial/resultado.json').read_text())
    assert adversarial['status']=='PASS_COM_RESSALVAS_DE_ENTRADA_E_METODO'
    assert adversarial['quadros']==15531 and adversarial['panoramas_reconstruidos']==31
    assert adversarial['equacoes']['status']==adversarial['estatistica']['status']=='PASS'
    assert adversarial['equacoes']['pdf_sha256']==sha(ROOT.parents[1]/'assets/schuck-bardo-wavelet-coherence.pdf')
    for name,digest in adversarial['scripts_sha256'].items():
        assert sha(ROOT/'scripts'/name)==digest,name
    assert adversarial['metadados']['status']=='PASS'
    metadata=read_csv(TABLES/'metadados_aquisicao.csv')
    assert len(metadata)==31
    for row in metadata:
        assert sha(DATA/row['ensaio']/'_config.xsv')==row['config_xsv_sha256']
        assert sha(DATA/row['ensaio']/'_config.mii')==row['config_mii_sha256']
    audited_inputs=read_csv(TABLES/'auditoria_entradas.csv')
    assert {r['ensaio'] for r in audited_inputs}=={r['ensaio'] for r in sources}
    for row in audited_inputs:
        assert float(row['erro_reconstrucao'])==float(row['erro_extracao'])==0
        assert float(row['erro_planilha'])<1e-10
        if row['ensaio'] in runs:assert row['rotulo_condicao_concordante']=='True'
    manifest=json.loads((INTER/'auditoria/adversarial/quadros_sha256.json').read_text())
    assert len(manifest)==31
    for record in manifest:
        assert len(record['arquivos'])==len(record['sha256'])==501
        for name,digest in zip(record['arquivos'],record['sha256']):
            assert sha(DATA/record['ensaio']/name)==digest,(record['ensaio'],name)
    # Recompute the exclusion summaries directly from the original pair table.
    pair_table=read_csv(TABLES/'pares_pearson_wco.csv')
    sensitivity=read_csv(TABLES/'sensibilidade_exclusao_m03q40.csv')
    assert len(sensitivity)==4
    for row in sensitivity:
        selected=[p for p in pair_table if row['conjunto']=='26_originais'
                  or 'm03q40_001' not in (p['ensaio_x'],p['ensaio_y'])]
        metric='pearson_r' if row['metrica']=='pearson' else 'escore_piso05'
        grouped={}
        means={}
        for p in selected:
            key=tuple(sorted((p['padrao_x'],p['padrao_y'])))
            grouped.setdefault(key,[]).append(float(p[metric]))
        for kind in ('intra','inter'):
            values=[float(p[metric]) for p in selected if p['tipo']==kind]
            assert len(values)==int(row[f'pares_{kind}'])
            means[kind]=float(np.mean(values))
            assert abs(means[kind]-float(row[kind]))<1e-13
            balanced=np.mean([np.mean(v) for (a,b),v in grouped.items() if (a==b)==(kind=='intra')])
            assert abs(balanced-float(row[f'{kind}_balanceado']))<1e-13
        assert abs(means['intra']-means['inter']-float(row['contraste']))<1e-13
    revision=OLD/'revisao_adversarial_20261027'
    for row in json.loads((revision/'manifesto.json').read_text()):
        assert sha(revision/row['arquivo'])==row['sha256']
    print('Auditoria adversarial: 15.531 quadros íntegros, rótulos e sensibilidade conferidos.',flush=True)
    original=json.loads((INTER/"auditoria/arquivamento_original.json").read_text())
    for row in original["files"]:
        assert sha(OLD/row["path"])==row["sha256"],row["path"]
    pilot=json.loads((INTER/"auditoria/arquivamento_piloto.json").read_text())
    for row in pilot:
        assert sha(OLD/"piloto_aprovado"/row["path"])==row["sha256"],row["path"]
    print("Arquivo: 583 arquivos originais e 17 do piloto preservados.",flush=True)
    for doc in [*STATS.glob("*.md"),*(ROOT/"05_conclusao").glob("*.md"),*(INTER/"galerias").glob("*.html")]:
        text=doc.read_text()
        targets=re.findall(r'\]\(([^)]+)\)',text) if doc.suffix==".md" else re.findall(r'(?:href|src)="([^"]+)"',text)
        for target in targets:
            if target.startswith(("https:","http:","#")):continue
            path=doc.parent/unquote(target.split("#")[0].strip("<>"))
            assert path.exists(),(doc,target)
    allowed={"01_panoramas","02_wavelets","03_wavelets_inter","04_estatistica","05_conclusao","intermediario","scripts","_OLD"}
    assert {p.name for p in ROOT.iterdir()}==allowed
    assert not list((ROOT/"scripts").glob("__pycache__"))
    assert not list(ROOT.glob("*.log"))
    artifacts=[]
    for path in sorted(ROOT.rglob("*")):
        relative=path.relative_to(ROOT)
        if path.is_file() and "_OLD" not in relative.parts and path.name not in ("manifesto_artefatos.json","verificacao_final.json"):
            artifacts.append(dict(arquivo=str(relative),bytes=path.stat().st_size,sha256=sha(path)))
    write_json(INTER/"auditoria/manifesto_artefatos.json",artifacts)
    result=dict(status="PASS",entradas_auditadas=31,ensaios=26,pares=325,pares_intra=87,pares_inter=238,
                figuras=355,figuras_panoramas=26,figuras_intra=87,figuras_inter=238,figuras_estatisticas=4,
                maior_erro_equacao_wco_salva=equation_error,matriz_pearson_conferida=True,
                medias_sem_diagonal=True,links_locais_validos=True,
                arquivos_originais_preservados=len(original["files"]),arquivos_piloto_preservados=len(pilot),
                menor_fonte_inserida_pt=min(r["fonte_minima_inserida_pt"] for r in layouts),
                maior_altura_inserida_cm=max(r["altura_inserida_cm"] for r in layouts),
                artefatos_ativos=len(artifacts),tamanho_ativo_bytes=sum(r["bytes"] for r in artifacts),
                core_igual_ao_original=True,raiz_sem_arquivos_avulsos=True)
    result.update(matrizes_identificadas_por_ensaio=True,
        auditoria_adversarial=adversarial['status'],
        qualidade_entrada='Ressalva em m03q40: 60/500 deslocamentos substituídos por 1 px e 25,43% dos pixels no teto 255; conjunto original preservado e sensibilidade com 25 disponível.',
        quadros_originais_conferidos=15531,sensibilidade_exclusao_conferida=True)
    write_json(INTER/"auditoria/verificacao_final.json",result)
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)


if __name__=="__main__":
    main()
