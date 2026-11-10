"""Auditoria independente dos resultados, da entrega e da preservação das fontes."""
from pathlib import Path
from itertools import combinations
import hashlib
import importlib.metadata
import json
import re
import sys
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from PIL import Image
from profile_core import ROOT, BASE, wc
from preflight import sha


def main():
    out=ROOT/'results'
    final=out/'final_verification.json'
    final.write_text(json.dumps(dict(status='RUNNING'))+'\n')
    manifest_path=out/'artifact_manifest.json'
    manifest_path.write_text(json.dumps(dict(status='RUNNING'))+'\n')
    cfg=json.loads((ROOT/'config.json').read_text())
    assert sha(ROOT/'PROTOCOL.md')==cfg['protocol_sha256']
    snapshot=json.loads((out/'source_snapshot.json').read_text())
    for rel,expected in snapshot.items():assert sha(BASE/rel)==expected, f'Fonte modificada: {rel}'
    inputs=json.loads((out/'input_hashes.json').read_text())
    for rel,expected in inputs.items():assert sha(ROOT.parent/rel)==expected, f'Entrada modificada: {rel}'
    previous=json.loads((BASE/'results/source_v16_sha256.json').read_text())
    for rel,expected in previous.items():assert sha(ROOT.parent/'v16_spatial2'/rel)==expected, f'Antecessor modificado: {rel}'
    # Verificação de fontes externas utilizadas para justificar a formulação.
    inherited_manifest=json.loads((BASE/'results/artifact_manifest.json').read_text())
    external={}
    for name,expected in inherited_manifest['source_sha256'].items():
        path=Path(name)
        assert path.exists() and sha(path)==expected, f'Fonte externa alterada: {path}'
        external[name]=expected
    assert sha(BASE/'scripts/wavelet_core.py')==sha(ROOT.parent/'v16_spatial2/scripts/wavelet_core.py')
    inv=pd.read_csv(BASE/'results/input_audit.csv');selected=inv[inv.included]
    assert len(inv)==31 and len(selected)==26
    assert sorted(selected.run)==cfg['expected_included_runs']
    audit=pd.read_csv(out/'input_verification.csv')
    assert len(audit)==31 and audit.erro_extracao.eq(0).all()

    normalized=pd.read_csv(out/'profiles_normalized.csv')
    maxmass=maxarea=0.
    for _,g in normalized.groupby(['aquisicao','variante','dominio']):
        g=g.sort_values('comprimento_onda_px')
        x=np.log(g.comprimento_onda_px.to_numpy());p=g.perfil.to_numpy()
        area=float(np.trapezoid(p,x=x))
        maxarea=max(maxarea,abs(area-g.area.iloc[0]))
        assert np.allclose(g.area,area,rtol=1e-11,atol=1e-12)
        step=(x[-1]-x[0])/(len(x)-1)
        assert np.allclose(np.diff(x),step,rtol=1e-10,atol=1e-13)
        weights=np.full(len(x),step);weights[[0,-1]]/=2
        independent=p*weights/area
        maxmass=max(maxmass,float(np.max(abs(independent-g.massa))))
        assert np.allclose(independent,g.massa,rtol=1e-10,atol=1e-12)
        assert np.isclose(g.massa.sum(),1,atol=1e-12)
        assert (g.massa>=0).all()
    main=normalized[(normalized.variante=='principal')&(normalized.dominio=='principal')]
    density=main.pivot(index='aquisicao',columns='comprimento_onda_px',values='densidade_log')
    values=pd.read_csv(out/'all_pairs.csv')
    expected=set(combinations(sorted(selected.run),2))
    summary=pd.read_csv(out/'ranking_summary.csv').set_index('metodo')
    maxdistance=maxranking=maxbalanced=0.
    for method,g in values.groupby('metodo'):
        assert len(g)==325 and not g[['aquisicao_x','aquisicao_y']].duplicated().any()
        assert set(zip(g.aquisicao_x,g.aquisicao_y))==expected
        assert g.mesmo_padrao.sum()==87
        assert np.isfinite(g[['valor','similaridade']]).all().all()
        if method.startswith('wco_'):
            assert np.array_equal(g.valor.to_numpy(),g.similaridade.to_numpy())
        else:
            assert np.array_equal(-g.valor.to_numpy(),g.similaridade.to_numpy())
            assert g.valor.between(-1e-14,1+1e-14).all()
        same=g[g.mesmo_padrao];cross=g[~g.mesmo_padrao]
        auc=float(mannwhitneyu(same.similaridade,cross.similaridade,method='asymptotic').statistic/(len(same)*len(cross)))
        maxranking=max(maxranking,abs(auc-summary.loc[method,'A_pares']))
        assert np.isclose(auc,summary.loc[method,'A_pares'],atol=1e-13)
        # Mesma estatística balanceada por outra decomposição: 4×10 comparações de grupos.
        group_auc=[]
        for _,s in same.groupby('grupo'):
            for _,c in cross.groupby('grupo'):
                group_auc.append(float(mannwhitneyu(s.similaridade,c.similaridade).statistic/(len(s)*len(c))))
        bal=float(np.mean(group_auc))
        maxbalanced=max(maxbalanced,abs(bal-summary.loc[method,'A_balanceada']))
        assert np.isclose(bal,summary.loc[method,'A_balanceada'],atol=1e-13)
    for row in values[values.metodo=='wps_principal'].itertuples():
        px=density.loc[row.aquisicao_x].to_numpy();py=density.loc[row.aquisicao_y].to_numpy()
        # Integral independente das densidades, sem chamar a implementação de Hellinger.
        dist=float(np.sqrt(np.trapezoid((np.sqrt(px)-np.sqrt(py))**2,x=np.log(density.columns.to_numpy()))/2))
        maxdistance=max(maxdistance,abs(dist-row.valor))
        assert np.isclose(dist,row.valor,atol=1e-12)
    deletion=pd.read_csv(out/'acquisition_deletion.csv')
    assert len(deletion)==26*values.metodo.nunique()
    assert (deletion.pares_intra+deletion.pares_entre).eq(300).all()
    assert deletion[deletion.aquisicao_retirada=='m07q01_001'].grupos_entre.eq(6).all()
    # Conferir que a influência publicada corresponde à retirada de todos os pares incidentes.
    for m,g in values.groupby('metodo'):
        for row in deletion[deletion.metodo==m].itertuples():
            sub=g[(g.aquisicao_x!=row.aquisicao_retirada)&(g.aquisicao_y!=row.aquisicao_retirada)]
            s=sub[sub.mesmo_padrao].similaridade;c=sub[~sub.mesmo_padrao].similaridade
            auc=float(mannwhitneyu(s,c).statistic/(len(s)*len(c)))
            assert np.isclose(auc,row.A_pares,atol=1e-13)
    prefix=pd.read_csv(out/'pair_prefix_profiles.csv')
    assert prefix.groupby(['aquisicao','n']).size().eq(200).all()
    assert len(prefix.groupby(['aquisicao','n']))==350
    maps=pd.read_csv(out/'illustrated_maps.csv')
    assert len(maps)==8 and maps.erro_equacao.max()<1e-10
    assert len(list((ROOT/'figures').glob('*.png')))==17
    pngs=list((ROOT/'figures').rglob('*.png'))
    for path in pngs:
        with Image.open(path) as im:im.verify()
        assert path.with_suffix('.pdf').read_bytes().startswith(b'%PDF-')
    broken=[]
    for path in ROOT.glob('*.md'):
        for target in re.findall(r'\]\(([^)]+)\)',path.read_text()):
            if target.startswith(('http:','https:','#','mailto:')):continue
            target=target.split('#')[0]
            if not (path.parent/target).exists():broken.append(f'{path.name}: {target}')
    assert not broken, broken
    html=(ROOT/'RESUMO_EXECUTIVO.html').read_text()
    embedded_figures=re.findall(r'^!\[[^\]]*\]\(([^)]+)\)',(ROOT/'EXECUTIVE_SUMMARY.md').read_text(),re.MULTILINE)
    assert len(embedded_figures)==5
    assert '<html lang="pt-BR">' in html and html.count('data:image/png;base64,')==len(embedded_figures)
    guide=json.loads((out/'reading_guide.json').read_text())
    for rel,expected_hash in {**guide['source_sha256'],**guide['output_sha256']}.items():
        assert sha(ROOT/rel)==expected_hash, f'Fonte/figura didática alterada: {rel}'
    same=values[(values.metodo=='wps_principal') & values.mesmo_padrao].valor.to_numpy()
    cross=values[(values.metodo=='wps_principal') & ~values.mesmo_padrao].valor.to_numpy()
    counts=dict(favorable=int((same[:,None]<cross).sum()),ties=int((same[:,None]==cross).sum()),
                unfavorable=int((same[:,None]>cross).sum()),total=len(same)*len(cross))
    assert guide['order_comparisons']==counts
    histogram=next(x for x in guide['figures'] if x['name']=='02_distribuicao_distancias')
    for group,distances in zip(histogram['groups'],[same,cross]):
        observed,_=np.histogram(distances,bins=histogram['bin_edges'])
        assert observed.sum()==len(distances)==group['n']
        assert observed.tolist()==group['counts']
        assert np.allclose(observed/len(distances),group['proportions'],rtol=0,atol=1e-14)
    for path in (ROOT/'scripts').glob('*.py'):compile(path.read_text(),str(path),'exec')
    environment=dict(python=sys.version,
        distribution_versions={p:importlib.metadata.version(p) for p in ['numpy','scipy','pandas','PyWavelets','matplotlib','Pillow','openpyxl']},
        pywavelets_module=wc.pywt.__version__,
        pywavelets_source_sha256=sha(Path(wc.pywt.__file__).parent/'_cwt.py'))
    result=dict(status='PASS',verified_utc=datetime.now(timezone.utc).isoformat(),
        base_files_unchanged=len(snapshot),predecessor_files_unchanged=len(previous),input_hashes_verified=len(inputs),
        external_sources_verified=len(external),audited_inputs=31,eligible_acquisitions=26,
        methods=int(values.metodo.nunique()),pairs_per_method=325,same_pairs=87,cross_pairs=238,
        independently_recomputed_main_distances=325,independent_ranking_methods=int(values.metodo.nunique()),
        distance_max_error=maxdistance,ranking_max_error=maxranking,balanced_ranking_max_error=maxbalanced,
        mass_max_error=maxmass,area_max_error=maxarea,acquisition_deletions=len(deletion),
        unique_pair_prefixes=350,illustrated_maps=8,map_equation_max_error=float(maps.erro_equacao.max()),
        illustrated_invalid_denominators=int(maps.denominadores_invalidos.sum()),
        png_verified=len(pngs),pdf_verified=len(pngs),broken_links=broken,
        guide_embedded_figures=len(embedded_figures),guide_order_comparisons=counts,
        protocol_hash_verified=True,core_unchanged=True,
        canonical_folder=ROOT.name,historical_alias_present=(ROOT.parent/'v16.1').is_symlink(),
        limitations='Verificação computacional e de procedência; sem nova reconstrução dos panoramas, validação física, inferência independente ou teste externo.')
    final.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    manifest=dict(created_utc=datetime.now(timezone.utc).isoformat(),environment=environment,
        base_directory=BASE.name,base_core_sha256=sha(BASE/'scripts/wavelet_core.py'),
        protocol_sha256=sha(ROOT/'PROTOCOL.md'),configuration_sha256=sha(ROOT/'config.json'),
        source_snapshot_sha256=sha(out/'source_snapshot.json'),external_sources=external,
        files={str(p.relative_to(ROOT)):sha(p) for p in sorted(ROOT.rglob('*'))
               if p.is_file() and p!=manifest_path and '__pycache__' not in p.parts})
    manifest_path.write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)


if __name__=='__main__':
    try:
        main()
    except Exception as e:
        (ROOT/'results/final_verification.json').write_text(json.dumps(dict(status='FAIL',error=str(e)),ensure_ascii=False,indent=2)+'\n')
        raise
