"""Independent audit of the complete batch; production artifacts are read only.

Evidence is written exclusively to intermediario/auditoria_lote_completo.
Uses original frames, independent one-dimensional event morphology, explicit
Gaussian convolution, temporal Morlet convolution and regenerated figure pixels.
Run stages: ancoras, numerica, controles, compatibilidade, figuras, concluir.
"""
import sys
sys.dont_write_bytecode = True
import argparse
import copy
import json
import math
import time
from datetime import datetime, timezone
from itertools import combinations

import numpy as np
from PIL import Image
from scipy.signal import find_peaks

from common import (ANCHORS, BASE, CFG, DATA, INTER, OLD, ROOT, SCRIPTS, WAVELETS,
                    context_ok, load_image, read_csv, sha, slug, write_csv, write_json)
from prepare_comparisons import load_catalog, load_individual, comparison_spec
from wavelet_cache import CACHE, get_grid
from audit_method_schuck import independent_prepare, independent_pair, explicit_cwt, max_error
from pair_identity import pair_key
import run_full_intra as full

OUT = INTER / "auditoria_lote_completo"


def protected():
    folders = (ANCHORS, WAVELETS, CACHE, INTER/"preparacao_comparacoes", full.MAPS)
    paths = [p for folder in folders for p in folder.rglob("*") if p.is_file()]
    paths += [p for p in SCRIPTS.iterdir() if p.is_file() and p.name != "audit_full_analysis.py"]
    paths += list(full.OUT.glob("*.json")) + list(full.OUT.glob("*.csv"))
    return {str(p.relative_to(ROOT)):dict(sha256=sha(p), mtime_ns=p.stat().st_mtime_ns) for p in paths}


def independent_intervals(pano, threshold):
    fraction = np.count_nonzero(pano[18:72] < threshold, axis=0)/54
    med = np.median(np.lib.stride_tricks.sliding_window_view(np.pad(fraction, 2, mode="edge"), 5), axis=1)
    mask = med >= .25
    ones = np.flatnonzero(mask)
    if len(ones) == 0:
        return []
    # Closing by a 13-column element fills internal gaps of at most 12 columns.
    breaks = np.flatnonzero(np.diff(ones) > 1)
    for index in breaks:
        if ones[index+1]-ones[index]-1 <= 12:
            mask[ones[index]+1:ones[index+1]] = True
    columns = np.flatnonzero(mask)
    breaks = np.flatnonzero(np.diff(columns) > 1)
    left = np.r_[columns[0], columns[breaks+1]]
    right = np.r_[columns[breaks]+1, columns[-1]+1]
    return [(int(a),int(b)) for a,b in zip(left,right) if b-a >= 96]


def required(center, left, right, margin):
    return max(1024, int(math.ceil(2*(max(center-left,right-center)+margin)/256))*256)


def event_options(run, pano, raw):
    reviews = {r["event_id"]:r for r in read_csv(OLD/"results/revisao_visual.csv")}
    old = {r["event_id"]:r for r in read_csv(OLD/"results/catalogo_eventos.csv")}
    variants = {q:independent_intervals(pano,q) for q in (55,60,65)}
    results = []
    for i,(left,right) in enumerate(variants[60],1):
        ident = f"{run}_E{i:03d}"
        center = (left+right)//2
        approved = ident in reviews and ident in old and reviews[ident]["aceito"] == "True"
        if ident in old:
            assert (left,right,center) == tuple(int(old[ident][k]) for k in ("inicio","fim","ancora"))
        bounds, shifts = [(left,right)], []
        for threshold in (55,65):
            matches = [(a,b) for a,b in variants[threshold]
                       if max(0,min(b,right)-max(a,left))/(max(b,right)-min(a,left)) >= .5]
            if len(matches) != 1:
                approved = False
                continue
            a,b = matches[0]
            bounds.append((a,b))
            shifts.append(max(abs(left-a),abs(right-b)))
        support_left,support_right = min(a for a,b in bounds),max(b for a,b in bounds)
        length = required(center,support_left,support_right,256)
        allowed = approved and 96 <= right-left <= 1536 and left>0 and right<len(raw) and context_ok(center,len(raw),length)
        if not allowed:
            continue
        contrast = min(np.median(raw[left-32:left])-np.median(raw[left:left+32]),
                       np.median(raw[right:right+32])-np.median(raw[right-32:right]))
        assert np.std(raw[center-length//2:center+length//2]) > 1e-10
        results.append(dict(anchor_id=ident,ancora=center,inicio=left,fim=right,
                            comprimento_minimo=length,inicio_suporte=support_left,fim_suporte=support_right,
                            variacao_maxima=max(shifts),contraste_minimo=float(contrast),
                            score=(max(shifts),-float(contrast),center)))
    return results


def independent_minima(raw, sigma):
    radius = int(4*sigma+.5)
    k = np.arange(-radius,radius+1)
    weights = np.exp(-k*k/(2*sigma*sigma)); weights /= weights.sum()
    smooth = np.convolve(np.pad(raw,radius,mode="symmetric"),weights,mode="valid")
    threshold = max(3., .25*float(np.diff(np.percentile(smooth,[5,95]))[0]))
    centers, props = find_peaks(-smooth, prominence=threshold, width=(16,512), distance=128, wlen=1025)
    return [dict(ancora=int(c),base_esquerda=int(props["left_bases"][i]),
                 base_direita=int(props["right_bases"][i]),prominencia=float(props["prominences"][i]),
                 largura=float(props["widths"][i]),inicio=float(props["left_ips"][i]),fim=float(props["right_ips"][i]),
                 limiar_prominencia=threshold) for i,c in enumerate(centers)]


def texture_options(run, raw):
    variants = {s:independent_minima(raw,s) for s in (12,16,20)}
    results = []
    for candidate in variants[16]:
        c = candidate["ancora"]
        bounds = [(candidate["base_esquerda"],candidate["base_direita"])]
        shifts = []
        for sigma in (12,20):
            matches = [r for r in variants[sigma] if abs(r["ancora"]-c)<=16 and
                       context_ok(r["ancora"],len(raw),required(r["ancora"],r["base_esquerda"],r["base_direita"],128))]
            if len(matches)!=1:
                break
            other = matches[0]
            shifts.append(abs(other["ancora"]-c))
            bounds.append((other["base_esquerda"],other["base_direita"]))
        if len(shifts)!=2:
            continue
        left,right = min(a for a,b in bounds),max(b for a,b in bounds)
        length = required(c,left,right,128)
        if context_ok(c,len(raw),length):
            results.append(dict(candidate,anchor_id=f"{run}_MIN{c:05d}",comprimento_minimo=length,
                                inicio_suporte=left,fim_suporte=right,variacao_maxima=max(shifts),
                                score=(-candidate["prominencia"],max(shifts),c)))
    return results


def anchors():
    catalog = load_catalog()
    candidate_rows = read_csv(INTER/"preparacao_comparacoes/candidatos.csv")
    results, source_hashes = [], {}
    for run,entry in catalog["aquisicoes"].items():
        pano,raw = load_image(run)
        with np.load(BASE/"intermediate"/f"{run}_input.npz") as saved:
            widths = saved["shifts"].astype(int)
            shifts = saved["raw_shifts"]
            np.testing.assert_array_equal(widths,np.maximum(np.r_[shifts,shifts[-1]],1))
        files = sorted((DATA/run).glob("ImgA*.jpg"))
        assert [f.name for f in files] == [f"ImgA{i:06d}.jpg" for i in range(501)]
        pieces = []
        for frame,width in zip(files,widths):
            with Image.open(frame) as im:
                values = np.asarray(im.convert("L"))
            first = values.shape[1]//2-width//2
            pieces.append(values[:,first:first+width])
            source_hashes[str(frame.relative_to(ROOT.parent))] = sha(frame)
        nominal = np.concatenate(pieces,axis=1)
        rebuilt = nominal.copy()
        seams = np.cumsum(widths)[:-1]
        for i,boundary in enumerate(seams):
            if min(widths[i],widths[i+1])>=4:
                rebuilt[:,boundary-1] = pieces[i+1][:,1]
        np.testing.assert_array_equal(pano,rebuilt)
        np.testing.assert_array_equal(raw,pano[18:72].sum(axis=0,dtype=np.uint64)/54)
        n = entry["perfis"]["nativo"]["contexto"]
        c = entry["perfis"]["comparavel"]["contexto"]
        options = event_options(run,pano,raw) if entry["padrao"] in ("plug","slug") else texture_options(run,raw)
        saved_eligible = {r["anchor_id"] for r in candidate_rows if r["aquisicao"]==run and r["elegivel_individual"]=="True"}
        assert {r["anchor_id"] for r in options} == saved_eligible, run
        rankings = {}
        for profile,record in (("nativo",n),("comparavel",c)):
            length = record["comprimento"]
            eligible = sorted((r for r in options if r["comprimento_minimo"]<=length and context_ok(r["ancora"],len(raw),length)),key=lambda r:r["score"])
            expected_rank = next(i for i,r in enumerate(eligible,1) if r["anchor_id"]==record["anchor_id"])
            assert record["ranking"] == expected_rank
            assert record["primeiro_candidato_no_contexto"] == eligible[0]["anchor_id"]
            actual = eligible[expected_rank-1]
            for field in ("ancora","inicio","fim","comprimento_minimo","inicio_suporte","fim_suporte","variacao_maxima"):
                np.testing.assert_allclose(record[field],actual[field],atol=1e-10,rtol=1e-12)
            assert record["segmento_inicio"] == record["ancora"]-length//2
            assert record["segmento_fim"] == record["ancora"]+length//2
            rankings[profile] = expected_rank
        assert rankings["nativo"] == 1
        assert n["ancora"] == c["ancora"] and n["anchor_id"] == c["anchor_id"]
        assert n["comprimento"] == CFG["janelas_por_padrao"][entry["padrao"]]
        if run=="m07q05_001":
            center = n["ancora"]
            changed = pano.copy(); changed[:,center] = nominal[:,center]
            alt = event_options(run,changed,changed[18:72].mean(axis=0))
            alt = sorted((r for r in alt if r["comprimento_minimo"]<=2048 and context_ok(r["ancora"],len(raw),2048)),key=lambda r:r["score"])
            assert alt[0]["ancora"]==center and alt[0]["anchor_id"]==n["anchor_id"]
            np.savez_compressed(OUT/"teste_coluna_costura.npz",raw=raw,raw_nominal=nominal[18:72].mean(axis=0),ancora=center)
        results.append(dict(aquisicao=run,padrao=entry["padrao"],ancora=n["ancora"],
                            panorama_reconstruido_exatamente=True,candidatos_conferidos=len(options),
                            ranking_nativo=rankings["nativo"],ranking_2048=rankings["comparavel"],
                            melhor_no_contexto_2048=c["primeiro_candidato_no_contexto"],
                            repetivel_ranking_nativo=n.get("repetivel_selecao_completa"),
                            repetivel_ranking_2048=c.get("repetivel_selecao_completa"),
                            deslocamento_local_maximo=n["variacao_maxima"],
                            ancora_em_coluna_substituida=n["ancora_na_coluna_substituida_feather"],
                            deslocamentos_nao_positivos_contexto=n["deslocamentos_nao_positivos_contexto"]))
        print(f"Fonte e ranking independentes: {run} / {n['ancora']}",flush=True)
    write_csv(OUT/"ancoras.csv",results)
    write_json(OUT/"quadros_fontes_sha256.json",source_hashes)
    return dict(ancoras=len(results),quadros_reconstruidos=len(source_hashes),
                sensiveis_ranking_nativo=[r["aquisicao"] for r in results if r["repetivel_ranking_nativo"] is False],
                sensiveis_ranking_2048=[r["aquisicao"] for r in results if r["repetivel_ranking_2048"] is False],
                ranking_2048_diferente=[r["aquisicao"] for r in results if r["ranking_2048"]>1],
                coluna_substituida=[r["aquisicao"] for r in results if r["ancora_em_coluna_substituida"]])


def numerical():
    import verify_rebuild
    catalog = load_catalog()
    verify_rebuild.OUT = OUT
    individuals,diagnostics = verify_rebuild.verify_individuals(catalog)
    guards = verify_rebuild.verify_guards(catalog)
    rows = json.loads((full.OUT/"metricas_wavelets.json").read_text())
    results = []
    maximum, minimum = 0., 1.
    for i,row in enumerate(rows,1):
        spec,records,maps = full.read_pair(catalog,row)
        x,y = [independent_prepare(load_image(spec[role])[1][records[spec[anchor]]["segmento_inicio"]:records[spec[anchor]]["segmento_fim"]])
               for role,anchor in (("run_x","anchor_x"),("run_y","anchor_y"))]
        expected = independent_pair(x,y,np.geomspace(50,1024,200))
        errors = {f"erro_{key}":max_error(expected[key],maps[key]) for key in ("wxx","wyy","wxy","wco")}
        np.testing.assert_array_equal(expected["valid"],maps["denominator_valid"])
        actual_max,actual_min = float(np.nanmax(maps["wco"])),float(np.nanmin(maps["wco"]))
        assert 0<=actual_min and actual_max<=1+1e-12
        maximum,minimum = max(maximum,actual_max),min(minimum,actual_min)
        # Direct scalar products check the FFT convolution indexing and phase.
        indices = np.array([0,49,99,149,199]); positions = np.array([0,2015,2047,2079,4094])
        explicit = explicit_cwt(np.correlate(x,y,"full")/2048,maps["wavelengths"][indices],positions)
        errors["erro_soma_direta"] = max_error(explicit,maps["wxy"][np.ix_(indices,positions)])
        assert spec["pattern_x"] == spec["pattern_y"] and spec["comprimento"] == 2048
        assert np.all(maps["support"][:,np.abs(maps["lags"])<=32])
        results.append(dict(pair_id=pair_key(spec),**errors,wco_minima=actual_min,wco_maxima=actual_max))
        if i%10==0:
            print(f"Auditoria numerica independente: {i}/87",flush=True)
    assert len(rows)==87 and len({(r["run_x"],r["run_y"]) for r in rows})==87
    write_csv(OUT/"pares.csv",results)
    return dict(individuais=len(individuals),perfis=len(diagnostics),pares=len(results),
                erros_individuais={k:max(r[k] for r in individuals) for k in ("erro_preparo","erro_correlacao","erro_cwt","erro_potencia_media")},
                erros_pares={k:max(r[k] for r in results) for k in errors},
                wco_minima=minimum,wco_maxima=maximum,protecoes_aprovadas=guards)


def controls():
    from wavelet_core import pair
    n = 2048; waves = np.geomspace(50,1024,200)
    rng = np.random.default_rng(20261025)
    x,y = rng.normal(size=(2,n))
    auto = pair(x,x,waves)
    err_self = max_error(auto["wco"][:,n-1],np.ones(200))
    a,b = pair(x,y,waves),pair(y,x,waves)
    err_reverse = max_error(b["wxy"],np.conj(a["wxy"][:,::-1]))
    err_reverse_coh = max_error(b["wco"],a["wco"][:,::-1])
    gain = pair(-3.7*x,.2*y,waves)
    err_gain = max_error(gain["wco"],a["wco"])
    u,v = np.zeros(n),np.zeros(n); u[700]=1; v[1200]=1
    delay = pair(u,v,waves)
    np.testing.assert_array_equal(delay["lags"][np.argmax(delay["wco"],axis=1)],np.full(200,-500))
    err_delay = max_error(delay["wco"][:,n-1-500],np.ones(200))
    unbiased = pair(u,v,waves,estimator="unbiased")
    expected_unbiased = (n/(n-500))**2
    max_error(unbiased["wco"][:,n-1-500],np.full(200,expected_unbiased))
    sinusoid = np.cos(2*np.pi*np.arange(n)/256)
    sinusoid -= sinusoid.mean()
    spectral = pair(sinusoid,sinusoid,waves)
    peak = float(waves[np.argmax(np.abs(spectral["wxx"][:,n-1]))])
    assert abs(peak/256-1)<.03
    catalogue = load_catalog()
    seam = np.load(OUT/"teste_coluna_costura.npz")
    record,a = load_individual(catalogue,"m07q05_001","comparavel")
    raw = a["raw"].copy(); index = record["ancora"]-record["segmento_inicio"]
    nominal = float(seam["raw_nominal"][record["ancora"]]); change = nominal-float(raw[index]); raw[index]=nominal
    corrected = independent_prepare(raw)
    effects = []
    for run in ("m05q01_001","m06q01_001"):
        _,b = load_individual(catalogue,run,"comparavel")
        original = pair(b["x"],a["x"],waves)
        variant = pair(b["x"],corrected,waves)
        delta = variant["wco"]-original["wco"]
        effects.append(dict(run_x=run,run_y="m07q05_001",diferenca_maxima_wco=float(np.nanmax(np.abs(delta))),
                            diferenca_media_central=float(np.mean(delta[:,n-1-32:n+32]))))
    write_csv(OUT/"sensibilidade_coluna_costura.csv",effects)
    return dict(erro_auto_coerencia=err_self,erro_troca_wcs=err_reverse,erro_troca_wco=err_reverse_coh,
                erro_invariancia_ganho=err_gain,erro_atraso=err_delay,atraso_esperado=-500,
                pico_lambda_seno_256=peak,wco_nao_enviesada_impulso=expected_unbiased,
                efeito_coluna_costura_intensidade=change,efeito_coluna_costura=effects,
                teste_costura_nao_altera_ancora=True)


def figures():
    import plots
    catalog = load_catalog()
    output = OUT/"regeneracao_conferencia"
    output.mkdir(exist_ok=True)
    old_paths = plots.INTER,plots.ANCHORS,plots.WAVELETS
    plots.INTER,plots.ANCHORS,plots.WAVELETS = output,output/"01_ancoras",output/"02_wavelets"
    results = []
    def compare(original,recreated):
        with Image.open(original) as a,Image.open(recreated) as b:
            np.testing.assert_array_equal(np.asarray(a),np.asarray(b))
            assert a.width>=1800 and a.height>=1900 and 16*a.height/a.width<=24
        results.append(dict(arquivo=str(original.relative_to(ROOT)),pixels_identicos=True,
                            bytes_identicos=sha(original)==sha(recreated),sha256=sha(original)))
    try:
        for run,item in catalog["aquisicoes"].items():
            record = item["perfis"]["nativo"]["contexto"]
            plots.anchor_figure(record,*load_image(run))
            name = f"{slug(record['padrao'])}_{record['anchor_id']}.png"
            compare(ANCHORS/name,plots.ANCHORS/name)
        rows = json.loads((full.OUT/"metricas_wavelets.json").read_text())
        labels = json.loads((full.OUT/"rotulos_visuais.json").read_text())
        for i,row in enumerate(rows,1):
            spec,records,maps = full.read_pair(catalog,row)
            label = labels[pair_key(spec)]
            assert label["mapa_sha256"] == row["mapa_sha256"]
            spec.pop("perfil")
            plots.pair_figure(spec,records,maps,row,label=label["rotulo"])
            name = f"{pair_key(spec)}_{label['rotulo']}.png"
            compare(WAVELETS/name,plots.WAVELETS/name)
            if i%10==0:
                print(f"Pixels das figuras conferidos: {i}/87 pares",flush=True)
    finally:
        plots.INTER,plots.ANCHORS,plots.WAVELETS = old_paths
    write_csv(OUT/"figuras.csv",results)
    return dict(imagens=len(results),pixels_identicos=sum(r["pixels_identicos"] for r in results),
                bytes_identicos=sum(r["bytes_identicos"] for r in results))


def compatibility():
    from collections import Counter
    from wavelet_cache import get_auto
    catalog = load_catalog()
    before = {str(p):p.stat().st_mtime_ns for p in CACHE.rglob("*") if p.is_file()}
    counts = Counter()
    for a,b in combinations(sorted(catalog["aquisicoes"]),2):
        spec,contexts,autos = comparison_spec(catalog,a,b,"comparavel")
        assert spec["comprimento"] == 2048
        assert autos[0]["meta"]["spec"]["grade"] == autos[1]["meta"]["spec"]["grade"]
        assert all(r["diagnosticos_comprimento"]==2048 for r in contexts.values())
        counts[spec["escopo"]] += 1
    for entry in catalog["aquisicoes"].values():
        for item in entry["perfis"].values():
            assert get_auto(item["contexto"])["key"] == item["cache_key"]
    assert before == {str(p):p.stat().st_mtime_ns for p in CACHE.rglob("*") if p.is_file()}
    assert counts == Counter(intra=87,inter=213)
    result = dict(pares_validados_por_metadados=dict(counts),inter_calculados=0,
                  reutilizacao_sem_reescrita=True,individuais=44,perfis=50,grades=2)
    write_json(OUT/"compatibilidade_futura.json",result)
    return result


def finish():
    baseline = json.loads((OUT/"arquivos_protegidos_inicio.json").read_text())
    assert protected()==baseline, "A auditoria nao deve modificar dados, imagens ou scripts de producao."
    reports = {name:json.loads((OUT/f"{name}.json").read_text()) for name in ("ancoras","numerica","controles","figuras")}
    assert reports["ancoras"]["ancoras"]==25 and reports["numerica"]["pares"]==87 and reports["figuras"]["pixels_identicos"]==112
    compatibility = json.loads((OUT/"compatibilidade_futura.json").read_text())
    assert compatibility["pares_validados_por_metadados"] == dict(intra=87,inter=213)
    assert compatibility["inter_calculados"] == 0 and compatibility["reutilizacao_sem_reescrita"]
    report = dict(status="PASS_COM_LIMITES_METODOLOGICOS",producao_preservada=True,
                  referencia_pdf_sha256=sha(ROOT.parents[1]/"assets/schuck-bardo-wavelet-coherence.pdf"),
                  paginas_impressas_conferidas=[411,413,414],
                  formulas="CWT 32.4, Morlet 32.6, lambda 32.8, alternativas 32.14-17 com K=1",
                  adaptacoes="Divisor N (biased), em lugar de N-|lag| das Eqs.32.21-23; dx=1 px. Potencia media e resumo adicional sobre lags com suporte.",
                  observacao_32_20="A Eq.32.20 impressa usa modulo, ao contrario da definicao complexa 32.16. Com K=1, modulo para exibir WCS e WCO sao os mesmos; nao generalizar a K>1 sem definir a media.",
                  limites=["Ancora escolhida no contexto nativo e preservada em 2048; quatro rankings locais diferentes.",
                           "Correspondencia local em sigma 12/16/20 nao garante repeticao do vencedor global: cinco casos nativos sensiveis.",
                           "Um centro coincide com coluna substituida pelo feather herdado; teste pontual registrado, sem validar fisicamente toda a montagem.",
                           "K=1 e rotulos favoravel/desfavoravel sao descritivos; nao representam significancia ou acuracia de classificacao.",
                           "Zoom das 19 ancoras de churn/wavy annular mostra 4096 px; comparacoes usam 2048 px e os mesmos centros."],
                  compatibilidade_futura=compatibility,
                  codigo_auditoria_sha256=sha(__file__),resultados=reports)
    write_json(OUT/"verificacao.json",report)
    return dict(status=report["status"],producao_preservada=True,imagens=112,pares=87,precomputacoes=44)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("etapa",choices=["ancoras","numerica","controles","compatibilidade","figuras","concluir"])
    args = parser.parse_args()
    OUT.mkdir(exist_ok=True)
    baseline = OUT/"arquivos_protegidos_inicio.json"
    if not baseline.exists():
        write_json(baseline,protected())
    start = time.perf_counter()
    result = {"ancoras":anchors,"numerica":numerical,"controles":controls,
              "compatibilidade":compatibility,"figuras":figures,"concluir":finish}[args.etapa]()
    result.update(data_utc=datetime.now(timezone.utc).isoformat(),segundos=time.perf_counter()-start)
    write_json(OUT/f"{args.etapa}.json",result)
    print(json.dumps(result,indent=2,ensure_ascii=False),flush=True)


if __name__=="__main__":
    main()
