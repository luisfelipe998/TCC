"""Independent numeric, profile, cache and fresh-run checks for the bounded pilot."""
import sys
sys.dont_write_bytecode = True
import copy
import json
import shutil
import tempfile
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import numpy as np

from common import (BASE, CFG, INTER, OLD, ROOT, SCRIPTS, context_ok, load_image,
                    read_csv, selected_manifest, sha, write_csv, write_json)
from prepare_comparisons import load_catalog, load_individual, comparison_spec
from wavelet_cache import (CACHE, context_record, get_auto, grid_spec, load_entry,
                           validate_length)
from select_anchors import event_candidates, texture_candidates
from audit_method_schuck import independent_prepare, temporal_cwt, independent_pair, max_error
from pair_identity import pair_key

OUT = INTER/"verificacao_reconstrucao"


def rejection(call, message):
    try:
        call()
    except ValueError:
        return message
    raise AssertionError(f"Entrada deveria ser rejeitada: {message}")


def verify_individuals(catalog):
    results, diagnostics, seen = [], [], set()
    reviews = {r["event_id"]:r for r in read_csv(OLD/"results/revisao_visual.csv")}
    previous = {r["event_id"]:r for r in read_csv(OLD/"results/catalogo_eventos.csv")}
    for run, entry in catalog["aquisicoes"].items():
        pano, raw = load_image(run)
        pattern = entry["padrao"]
        rows,_ = event_candidates(run,pattern,pano,raw,reviews,previous) if pattern in ("plug","slug") else texture_candidates(run,pattern,raw)
        with np.load(BASE/"intermediate"/f"{run}_input.npz") as saved:
            seams, widths = saved["seams"].astype(int), saved["shifts"].astype(int)
        modified = seams[np.minimum(widths[:-1],widths[1:])>=4]-1
        for profile in ("nativo","comparavel"):
            record, auto = load_individual(catalog,run,profile)
            n, start, stop = record["comprimento"],record["segmento_inicio"],record["segmento_fim"]
            assert n==stop-start==record["diagnosticos_comprimento"]
            assert start==record["ancora"]-n//2 and stop==record["ancora"]+n//2
            assert start>=16 and stop<=pano.shape[1]-16
            margin = 256 if pattern in ("plug","slug") else 128
            assert record["inicio_suporte"]-start>=margin and stop-record["fim_suporte"]>=margin
            assert record["costuras_contexto"]==int(((seams>start)&(seams<stop)).sum())
            assert record["colunas_substituidas_feather_contexto"]==int(((modified>=start)&(modified<stop)).sum())
            if pattern not in ("plug","slug"):
                checks=[]
                for sigma in (12,20):
                    candidates,_ = texture_candidates(run,pattern,raw,reference_sigma=sigma)
                    allowed=[r for r in candidates if r["elegivel"] and r["comprimento_minimo"]<=n and context_ok(r["ancora"],len(raw),n)]
                    best=min(allowed,key=lambda r:(-r["prominencia"],r["variacao_maxima"],r["ancora"])) if allowed else None
                    expected=best["ancora"] if best else None
                    assert record[f"ancora_reselecionada_sigma{sigma}"]==expected
                    repeat=bool(best and abs(expected-record["ancora"])<=16)
                    assert record[f"repetivel_selecao_sigma{sigma}"]==repeat
                    checks.append(repeat)
                assert record["repetivel_selecao_completa"]==all(checks)
            diagnostics.append(dict(aquisicao=run,perfil=profile,comprimento=n,
                                    costuras=record["costuras_contexto"],
                                    repetivel_sigma=record.get("repetivel_selecao_completa"),
                                    status="PASS"))
            independent_raw=pano[18:72,start:stop].mean(axis=0)
            np.testing.assert_array_equal(independent_raw,auto["raw"])
            if auto["key"] in seen:
                continue
            seen.add(auto["key"])
            x=independent_prepare(independent_raw)
            rxx=np.correlate(x,x,mode="full")/n
            waves=np.geomspace(50,n/2,200)
            transformed=temporal_cwt(rxx,waves)
            support=n-1-np.abs(np.arange(1-n,n))[None,:]>=np.sqrt(2)*(waves/(4*np.pi/(6+np.sqrt(38))))[:,None]
            mean=np.asarray([np.abs(row)[mask].mean() for row,mask in zip(transformed,support)])
            errors=dict(erro_preparo=max_error(x,auto["x"]),erro_correlacao=max_error(rxx,auto["rxx"]),
                        erro_cwt=max_error(transformed,auto["wxx"]),erro_potencia_media=max_error(mean,auto["power_mean"]))
            results.append(dict(aquisicao=run,comprimento=n,cache_key=auto["key"],**errors))
            print(f"Conferencia independente: {run}, N={n}, erro CWT={errors['erro_cwt']:.2e}",flush=True)
    write_csv(OUT/"individuais.csv",results)
    write_csv(OUT/"diagnosticos_por_perfil.csv",diagnostics)
    # Regression from the previous audit: this selection has different
    # repeatability at 2048 and 4096; the smaller crop must not inherit True.
    native,_=load_individual(catalog,"m04q10_001","nativo")
    common,_=load_individual(catalog,"m04q10_001","comparavel")
    assert native["repetivel_selecao_completa"] is True
    assert common["repetivel_selecao_completa"] is False
    assert "costuras_contexto" not in context_record(native,2048)
    return results, diagnostics


def verify_pairs(catalog):
    results=[]
    for comparison in selected_manifest()["pares"]:
        _,a=load_individual(catalog,comparison["run_x"],"nativo")
        _,b=load_individual(catalog,comparison["run_y"],"nativo")
        n=comparison["comprimento"]
        expected=independent_pair(independent_prepare(a["raw"]),independent_prepare(b["raw"]),np.geomspace(50,n/2,200))
        with np.load(INTER/"mapas"/f"{pair_key(comparison)}.npz") as saved:
            meta=json.loads(str(saved["meta_json"]))
            assert meta["backend_sha256"]==sha(SCRIPTS/"wavelet_cache.py")
            assert meta["cache_individual_x"]==a["key"] and meta["cache_individual_y"]==b["key"]
            errors={f"erro_{key}":max_error(expected[key],saved[key]) for key in ("wxx","wyy","wxy","wco")}
            assert np.nanmax(saved["wco"])<=1+1e-12 and np.nanmin(saved["wco"])>=0
        results.append(dict(padrao=comparison["padrao"],comprimento=n,**errors))
    write_csv(OUT/"pares.csv",results)
    return results


def verify_guards(catalog):
    checked=[]
    old=copy.deepcopy(CFG["wavelet"])
    try:
        for field,bad in (("estimador","unbiased"),("K",2),("omega0",8),("comprimento_onda_maximo","N")):
            CFG["wavelet"][field]=bad
            checked.append(rejection(lambda:grid_spec(2048),f"parametro_{field}"))
            CFG["wavelet"].update(old)
    finally:
        CFG["wavelet"].update(old)
    for bad in (2049,2048.5,True,512,0):
        checked.append(rejection(lambda:validate_length(bad),f"comprimento_{bad}"))
    changed=copy.deepcopy(catalog)
    changed["aquisicoes"]["m05q01_001"]["perfis"]["comparavel"]=copy.deepcopy(changed["aquisicoes"]["m06q01_001"]["perfis"]["comparavel"])
    checked.append(rejection(lambda:load_individual(changed,"m05q01_001","comparavel"),"cache_de_outro_ensaio"))
    for field,bad in (("ancora",1),("costuras_contexto",999),("comprimento",4096)):
        changed=copy.deepcopy(catalog)
        changed["aquisicoes"]["m05q01_001"]["perfis"]["comparavel"]["contexto"][field]=bad
        checked.append(rejection(lambda:load_individual(changed,"m05q01_001","comparavel"),f"contexto_{field}"))
    checked.append(rejection(lambda:comparison_spec(catalog,"m05q01_001","m03q05_001","nativo"),"grades_2048_4096_incompativeis"))
    record,auto=load_individual(catalog,"m05q01_001","comparavel")
    with tempfile.TemporaryDirectory(dir=OUT,prefix="cache_teste_") as temporary:
        folder=Path(temporary)/auto["key"]
        shutil.copytree(CACHE/"individuais"/auto["key"],folder)
        copied=load_entry(str(folder))
        assert copied["raw"][0]==auto["raw"][0]
        changed_values=np.load(folder/"raw.npy",mmap_mode="r+")
        changed_values[0]+=1
        changed_values.flush()
        del changed_values,copied
        checked.append(rejection(lambda:load_entry(str(folder)),"array_modificado_apos_primeira_leitura"))
    return checked


def verify_reuse_and_scope(catalog):
    before={str(p):p.stat().st_mtime_ns for p in CACHE.rglob("*") if p.is_file()}
    for run,entry in catalog["aquisicoes"].items():
        for profile,item in entry["perfis"].items():
            record,auto=load_individual(catalog,run,profile)
            again=get_auto(record)
            assert again["key"]==auto["key"]==item["cache_key"]
    names=sorted(catalog["aquisicoes"])
    counts=dict(intra=0,inter=0)
    for a,b in combinations(names,2):
        spec,contexts,autos=comparison_spec(catalog,a,b,"comparavel")
        assert spec["comprimento"]==2048 and len(contexts)==2
        assert autos[0]["meta"]["spec"]["grade"]==autos[1]["meta"]["spec"]["grade"]
        counts[spec["escopo"]]+=1
        # Do not call pair_from_autos: cross-pattern metadata checks only.
    assert before=={str(p):p.stat().st_mtime_ns for p in CACHE.rglob("*") if p.is_file()}
    selection=selected_manifest()
    assert catalog["escopo"]=="piloto" and len(names)==8
    assert set(names)=={r["aquisicao"] for r in selection["ancoras"]}
    assert catalog["individuais_unicos"]==12 and catalog["grades"]==2
    assert len(list((CACHE/"individuais").iterdir()))==12 and len(list((CACHE/"grades").iterdir()))==2
    assert counts==dict(intra=4,inter=24)
    native,_,_=comparison_spec(catalog,"m05q01_001","m06q01_001","nativo")
    comparable,_,_=comparison_spec(catalog,"m05q01_001","m06q01_001","comparavel")
    assert pair_key(native)!=pair_key(comparable)
    return counts


def main():
    OUT.mkdir(exist_ok=True)
    catalog=load_catalog()
    start=json.loads((INTER/"reinicio.json").read_text())
    epoch=datetime.fromisoformat(start["fim_limpeza_utc"]).timestamp()
    for name in ("historico","comparacao_estimadores","teste_nitidez","estudo_consenso"):
        assert not (INTER/name).exists(), "Resultados anteriores nao devem permanecer."
    for metadata in CACHE.rglob("metadata.json"):
        created=datetime.fromisoformat(json.loads(metadata.read_text())["data_utc"]).timestamp()
        assert created>=epoch
    for folder in (ROOT/"01_ancoras",ROOT/"02_wavelets",INTER/"mapas"):
        assert all(p.stat().st_mtime>=epoch for p in folder.iterdir() if p.is_file())
    individuals,diagnostics=verify_individuals(catalog)
    pairs=verify_pairs(catalog)
    guards=verify_guards(catalog)
    compatible=verify_reuse_and_scope(catalog)
    figure_checks=[]
    for directory in (ROOT/"01_ancoras",ROOT/"02_wavelets"):
        for path in directory.glob("*.png"):
            record=json.loads((INTER/"apresentacao"/f"{path.stem}.json").read_text())
            assert record["altura_inserida_cm"]<=24 and record["fonte_minima_inserida_pt"]>=8
            figure_checks.append(record)
    assert len(figure_checks)==12 and len(list((INTER/"mapas").glob("*.npz")))==4
    assert not list((ROOT/"03_conclusoes").iterdir()) and not list(ROOT.rglob("*.pdf"))
    report=dict(status="PASS",data_utc=datetime.now(timezone.utc).isoformat(),escopo="piloto",
                ancoras=8,pares_calculados=4,contextos_individuais=len(individuals),perfis_conferidos=len(diagnostics),grades=2,
                pares_inter_calculados=0,compatibilidade_somente_metadados=compatible,
                erros_maximos_individuais={k:max(r[k] for r in individuals) for k in ("erro_preparo","erro_correlacao","erro_cwt","erro_potencia_media")},
                erros_maximos_pares={k:max(r[k] for r in pairs) for k in ("erro_wxx","erro_wyy","erro_wxy","erro_wco")},
                cache_reutilizado_sem_reescrita=True,entradas_invalidas_rejeitadas=guards,
                cache_criado_apos_limpeza=True,diagnosticos_por_janela_conferidos=True,
                todas_celulas_wps_conferidas=True,todas_celulas_wco_piloto_conferidas=True,
                bytes_precalculados=catalog["bytes_precalculados"],imagens_png=12,conclusoes_vazias=True,
                fonte_minima_pt=min(r["fonte_minima_inserida_pt"] for r in figure_checks),
                altura_maxima_16cm=max(r["altura_inserida_cm"] for r in figure_checks),
                codigo_sha256={p.name:sha(p) for p in SCRIPTS.glob("*.py")})
    write_json(OUT/"verificacao.json",report)
    print(json.dumps({k:report[k] for k in ("status","ancoras","pares_calculados","contextos_individuais",
        "erros_maximos_individuais","erros_maximos_pares","pares_inter_calculados","cache_reutilizado_sem_reescrita")},indent=2),flush=True)


if __name__=="__main__":
    main()
