"""Final integration checks: cohort, pair joins, exported data and report links.

Numerical recomputation is performed by run_full_inter verify and
statistics_anchors calculate; this step checks the delivered artifacts.
"""
import sys
sys.dont_write_bytecode=True
import json
import re
from collections import Counter
from itertools import combinations
from urllib.parse import unquote

import numpy as np
from PIL import Image
from common import INTER, ROOT, SCRIPTS, read_csv, sha, write_json
from pair_identity import pair_key
from run_full_inter import filename_base


def main():
    out=INTER/"verificacao_entrega"
    out.mkdir(exist_ok=True)
    intercheck=json.loads((INTER/"lote_inter/verificacao.json").read_text())
    statcheck=json.loads((INTER/"estatistica/verificacao.json").read_text())
    assert intercheck["status"]==statcheck["status"]=="PASS"
    assert statcheck["script_sha256"]==sha(SCRIPTS/"statistics_anchors.py")
    protocol=json.loads((INTER/"lote_inter/protocolo.json").read_text())
    assert protocol["script_sha256"]==sha(SCRIPTS/"run_full_inter.py")
    index=read_csv(ROOT/"04_estatistica/dados/indice_ensaios.csv")
    runs=[v["ensaio"] for v in index]
    assert len(set(runs))==25 and not any(v.endswith("_002") for v in runs)
    rows=read_csv(ROOT/"04_estatistica/dados/pares_pearson_wco.csv")
    assert len(rows)==300
    assert {frozenset((v["ensaio_x"],v["ensaio_y"])) for v in rows}=={
        frozenset(pair) for pair in combinations(runs,2)}
    waves=[]
    for group in ("lote_intra","lote_inter"):
        waves+=json.loads((INTER/group/"metricas_wavelets.json").read_text())
    lookup={frozenset((v["run_x"],v["run_y"])):v for v in waves}
    with np.load(INTER/"estatistica/descritores.npz") as a:
        assert list(a["runs"])==runs
        for row in rows:
            i,j=int(row["numero_x"])-1,int(row["numero_y"])-1
            assert runs[i]==row["ensaio_x"] and runs[j]==row["ensaio_y"] and i<j
            assert float(row["pearson_r"])==a["r"][i,j]
            assert float(row["pearson_r_completo"])==a["r_full"][i,j]
            source=lookup[frozenset((runs[i],runs[j]))]
            assert row["tipo"]==source["escopo"]
            assert float(row["wco_media_central"])==source["wco_media"]
            assert row["mapa_sha256"]==source["mapa_sha256"]
        for name,key in (("matriz_pearson_todos","r"),("matriz_pearson_completo","r_full")):
            csv=read_csv(ROOT/"04_estatistica/dados"/(name+".csv"))
            assert [r["ensaio"] for r in csv]==runs
            np.testing.assert_array_equal(np.array([[float(row[r]) for r in runs] for row in csv]),a[key])
    images=json.loads((INTER/"lote_inter/figuras_wavelets.json").read_text())
    labels=json.loads((INTER/"lote_inter/rotulos_visuais.json").read_text())
    interrows=json.loads((INTER/"lote_inter/metricas_wavelets.json").read_text())
    expected={filename_base(v)+"_"+labels[pair_key(v)]["rotulo"]+".png" for v in interrows}
    assert expected=={v["arquivo"] for v in images}=={p.name for p in (ROOT/"03_wavelets_inter").iterdir()}
    inspected={(v["evidencia"],v["posicao_na_galeria"]) for v in labels.values()}
    assert len(inspected)==213
    for value in labels.values():
        assert sha(ROOT/value["evidencia"])==value["evidencia_sha256"]
    deliverables=[p for folder in ("03_wavelets_inter","04_estatistica","05_conclusao")
                  for p in (ROOT/folder).rglob("*") if p.is_file()]
    assert not any(p.suffix.lower()==".pdf" for p in deliverables)
    broken=[]
    link_count=0
    for path in deliverables:
        if path.suffix==".md":
            for target in re.findall(r"\]\(([^)]+)\)",path.read_text()):
                if target.startswith(("http:","https:","#")):
                    continue
                target=unquote(target.strip("<>").split("#")[0])
                link_count+=1
                if not (path.parent/target).exists():
                    broken.append([str(path),target])
        elif path.suffix==".png":
            with Image.open(path) as im:
                assert im.format=="PNG"
                im.verify()
    assert not broken,broken
    manifest={str(p.relative_to(ROOT)):sha(p) for p in deliverables}
    manifest.update({str(p.relative_to(ROOT)):sha(p) for p in
        (SCRIPTS/"run_full_inter.py",SCRIPTS/"record_inter_review.py",SCRIPTS/"statistics_anchors.py",
         SCRIPTS/"write_executive.py",SCRIPTS/"validate_inter_statistics.py",INTER/"estatistica/descritores.npz")})
    write_json(out/"manifesto_sha256.json",manifest)
    result=dict(status="PASS",ensaios=25,pares_intra=87,pares_inter=213,
        celulas_pearson=625,pares_pearson=300,figuras_inter=213,figuras_estatisticas=5,
        rotulos_inter=dict(Counter(v["rotulo"] for v in labels.values())),
        links_locais_validos=link_count,arquivos_finais=len(deliverables),
        coeficientes_csv_iguais_ao_npz=True,pares_duplicados=0,arquivos_pdf=0,
        auditoria_inter="../lote_inter/verificacao.json",auditoria_estatistica="../estatistica/verificacao.json")
    write_json(out/"verificacao.json",result)
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
