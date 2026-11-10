"""Recalcula a v16 integral e gera as figuras intra e interpadrões aprovadas."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
import json
import os
import sys
import time

sys.dont_write_bytecode = True
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
from common import *
import wavelet_core as core
import figuras

_reference = None
_intra_reference = None


def zero_time(rxx, wavelengths=WAVES):
    lags = np.arange(-(len(rxx)//2), len(rxx)//2+1)
    values = []
    for scale in wavelengths/core.FOURIER_FACTOR:
        u = lags/scale
        kernel = np.pi**(-.25)/np.sqrt(scale)*np.exp(-u*u/2)*np.cos(6*u)
        values.append(abs(np.sum(rxx*kernel)))
    return np.asarray(values)


def audit_inputs():
    config = json.loads((OLD/"results/analysis_config.json").read_text())
    assert config["extraction_rows_inclusive"] == [18, 71]
    assert config["wavelengths_px"] == [50, 3000] and config["primary_scale_count"] == 200
    assert config["primary_correlation"] == "biased" and not config["coi"] and not config["conservative_mask"]
    assert sha(ROOT/"scripts/wavelet_core.py") == sha(OLD/"scripts/wavelet_core.py")
    inventory = read_csv(OLD/"results/input_audit.csv")
    audit, selected = [], []
    for row in inventory:
        run = row["run"]
        image_path = DATA/run/"__panorama_v11.png"
        workbook_path = DATA/run/"__dados_cropped_v11.xlsx"
        assert sha(image_path) == row["panorama_sha256"], run
        assert sha(workbook_path) == row["xlsx_sha256"], run
        pano = load_image(run)
        raw = pano[ROW_START:ROW_STOP].sum(axis=0, dtype=np.uint64)/ROW_COUNT
        assert pano.shape == (int(row["height"]), int(row["n"]))
        np.testing.assert_array_equal(raw, pano[ROW_START:ROW_STOP].mean(axis=0))
        with np.load(OLD/"intermediate"/f"{run}_input.npz") as old:
            np.testing.assert_array_equal(raw, old["raw"])
        raw_sha = hashlib.sha256(raw.astype("<f8", copy=False).tobytes()).hexdigest()
        assert raw_sha == row["raw_array_sha256"]
        np.savez_compressed(INTER/"entradas"/f"{run}.npz", raw=raw)
        current = dict(ensaio=run, padrao=pattern_name(row["pattern"]),
                       incluido=row["included"] == "True", comprimento=len(raw),
                       linhas_inicio=18, linhas_fim=71, erro_extracao=0.,
                       panorama_sha256=row["panorama_sha256"], planilha_sha256=row["xlsx_sha256"],
                       raw_sha256=raw_sha)
        audit.append(current)
        if current["incluido"]:
            selected.append(dict(current))
    selected.sort(key=lambda row: (PATTERNS.index(row["padrao"]), row["ensaio"]))
    assert len(inventory) == 31 and len(selected) == 26
    for i, row in enumerate(selected, 1):
        row["numero"] = i
    write_csv(INTER/"auditoria/entradas.csv", audit)
    write_json(INTER/"aquisicoes.json", selected)
    write_csv(TABLES/"indice_ensaios.csv", [dict(numero=r["numero"], ensaio=r["ensaio"],
              padrao=r["padrao"], comprimento_px=r["comprimento"]) for r in selected])
    write_json(INTER/"configuracao.json", dict(ensaios=26, auditados=31, pares=325, intra=87, inter=238,
        linhas_inclusivas=[18,71], wavelengths_px=[50,3000], escalas=200,
        autocorrelacao="enviesada, divisor N", omega0=6, amostragem_px=1,
        selecao_par="prefixo de comprimento min(Na,Nb); preparo após o corte",
        dominio="todos os lags, sem COI ou máscara geométrica",
        wps="abs(CWT(Rxx))", wcs="CWT(Rxy)", wco="abs(WCS)^2 / (WPSa(0)*WPSb(0))",
        potencia_media="média da WPS em todos os lags",
        escore_principal="máximos da WCO por escala, ponderados pela média geométrica das WPS no lag zero; piso 5%",
        pearson="resíduos de log(WPS(0)+1e-12) após ajuste contra 1 e log(lambda), panoramas completos",
        core_sha256=sha(ROOT/"scripts/wavelet_core.py")))
    print("Entradas: 31 panoramas e planilhas conferidos; 26 ensaios elegíveis.", flush=True)
    return selected


def individual(row):
    run, pattern = row["ensaio"], row["padrao"]
    raw = load_raw(run)
    x = core.prepare(raw)
    rxx = core.correlation(x, x)
    wxx, scales = core.transform_correlation(rxx, WAVES)
    zero = np.abs(wxx[:, len(x)-1])
    mean = np.abs(wxx).mean(axis=1)
    independent = zero_time(rxx)
    np.testing.assert_allclose(zero, independent, atol=1e-11, rtol=1e-11)
    prior = [r for r in read_csv(OLD/"results/profiles.csv") if r["run"] == run and r["method"] == "v16_biased"]
    assert len(prior) == 200
    np.testing.assert_allclose(WAVES, [float(r["wavelength_px"]) for r in prior], atol=1e-12)
    old_zero = np.array([float(r["zero_lag_wps"]) for r in prior])
    old_mean = np.array([float(r["lag_mean_all"]) for r in prior])
    np.testing.assert_allclose(zero, old_zero, atol=1e-11, rtol=1e-11)
    np.testing.assert_allclose(mean, old_mean, atol=1e-11, rtol=1e-11)
    np.savez_compressed(INTER/"individuais"/f"{run}.npz", wavelengths=WAVES, prepared=x,
                        wps_zero=zero, wps_mean_all=mean)
    layout = figuras.panorama_figure(run, load_image(run), raw, x, pattern=pattern)
    result = dict(ensaio=run, padrao=pattern, comprimento=len(raw),
                  pico_wps_zero_lambda_px=float(WAVES[np.argmax(zero)]),
                  pico_potencia_media_lambda_px=float(WAVES[np.argmax(mean)]),
                  erro_wps_independente=float(abs(zero-independent).max()),
                  erro_perfil_vs_original=float(max(abs(zero-old_zero).max(), abs(mean-old_mean).max())),
                  npz_sha256=sha(INTER/"individuais"/f"{run}.npz"), apresentacao=layout)
    return result


def specifications(selected):
    tasks = []
    for a, b in combinations(sorted(selected, key=lambda r: r["ensaio"]), 2):
        same = a["padrao"] == b["padrao"]
        pattern = a["padrao"] if same else a["padrao"]+" × "+b["padrao"]
        prefix = slug(a["padrao"]) if same else slug(a["padrao"])+"_"+slug(b["padrao"])
        name = prefix+"_"+a["ensaio"]+"_"+b["ensaio"]
        folder = "02_wavelets" if same else "03_wavelets_inter"
        tasks.append(dict(run_x=a["ensaio"], run_y=b["ensaio"], padrao_x=a["padrao"],
                          padrao_y=b["padrao"], tipo="intra" if same else "inter",
                          padrao=pattern, arquivo=f"{folder}/{name}.png",
                          input_x_sha256=a["raw_sha256"], input_y_sha256=b["raw_sha256"]))
    assert len(tasks) == 325 and sum(t["tipo"] == "intra" for t in tasks) == 87
    return tasks


def references():
    global _reference, _intra_reference
    if _reference is None:
        _reference = {(r["run_x"],r["run_y"]): r for r in read_csv(OLD/"results/preprocessing_confirm_pairs.csv")
                      if r["variant"] == "baseline"}
        _intra_reference = {(r["run_x"],r["run_y"]): r for r in read_csv(OLD/"results/pair_methods.csv")
                            if r["method"] == "v16_biased"}
        assert len(_reference) == 325 and len(_intra_reference) == 87
    return _reference, _intra_reference


def calculate_pair(spec, force=False):
    a, b = spec["run_x"], spec["run_y"]
    stem = INTER/"pares"/f"{a}__{b}"
    marker = stem.with_suffix(".json")
    signatures = dict(core=sha(ROOT/"scripts/wavelet_core.py"), analysis=sha(__file__),
                      input_x=spec["input_x_sha256"], input_y=spec["input_y_sha256"])
    if not force and marker.exists() and stem.with_suffix(".npz").exists():
        saved = json.loads(marker.read_text())
        if saved["assinaturas"] == signatures and saved["figuras_sha256"] == sha(ROOT/"scripts/figuras.py") and (ROOT/spec["arquivo"]).exists():
            return saved
    raw_a, raw_b = load_raw(a), load_raw(b)
    n = min(len(raw_a), len(raw_b))
    x, y = core.prepare(raw_a[:n]), core.prepare(raw_b[:n])
    maps = core.pair(x, y, WAVES)
    co, j0 = maps["wco"], n-1
    assert maps["denominator_valid"].all() and np.isfinite(co).all()
    assert float(co.min()) >= 0 and float(co.max()) <= 1+1e-8
    zx, zy = np.abs(maps["wxx"][:, j0]), np.abs(maps["wyy"][:, j0])
    independent_errors = []
    indices = [0, 31, 79, 127, 165, 199]
    for prepared, zero in ((x,zx),(y,zy)):
        independent = zero_time(core.correlation(prepared, prepared), WAVES[indices])
        independent_errors.append(float(abs(independent-zero[indices]).max()))
        np.testing.assert_allclose(independent, zero[indices], atol=1e-11, rtol=1e-11)
    rxy = core.correlation(x,y)
    direct_errors = []
    for lag in (-n+1, -101, -1, 0, 1, 101, n-1):
        expected = np.dot(x[lag:], y[:n-lag])/n if lag >= 0 else np.dot(x[:n+lag], y[-lag:])/n
        direct_errors.append(float(abs(rxy[j0+lag]-expected)))
    assert max(direct_errors) < 1e-12
    weight = np.sqrt(zx*zy)
    peak = co.max(axis=1)
    scores = {}
    for cut in (.01,.05,.10):
        active = (zx >= cut*zx.max()) & (zy >= cut*zy.max())
        weighted = np.where(active,weight,0)
        scores[f"escore_piso{int(cut*100):02d}"] = float(np.sum(weighted*peak)/weighted.sum())
        if cut == .05:
            powered_count = int(active.sum())
    peak_scale, peak_lag = np.unravel_index(np.argmax(co),co.shape)
    metrics = dict(comprimento=n, comprimento_x=len(raw_a), comprimento_y=len(raw_b),
                   descarte_x=len(raw_a)-n, descarte_y=len(raw_b)-n,
                   escalas_validas=int(maps["denominator_valid"].sum()), escalas_piso05=powered_count,
                   wco_maxima=float(co.max()), wco_media_integral=float(co.mean()),
                   wco_maxima_zero=float(co[:,j0].max()), wco_mediana_zero=float(np.median(co[:,j0])),
                   wco_ponderada_zero=float(np.sum(weight*co[:,j0])/weight.sum()),
                   pico_lag_px=int(maps["lags"][peak_lag]), pico_lambda_px=float(WAVES[peak_scale]), **scores)
    reference, intra_reference = references()
    old = reference[(a,b)]
    old_errors = [abs(metrics["wco_maxima"]-float(old["max_wco_all"])),
                  abs(metrics["wco_ponderada_zero"]-float(old["weighted_zero_lag_wco"]))]
    for suffix in ("01","05","10"):
        old_errors.append(abs(metrics[f"escore_piso{suffix}"]-float(old[f"weighted_scale_peak_power{suffix}"])))
    assert powered_count == int(old["eligible_scales_power05"])
    assert max(old_errors) < 1e-11, (a,b,old_errors)
    map_error = None
    if (a,b) in intra_reference:
        baseline = intra_reference[(a,b)]
        for key, old_key in (("wco_maxima_zero","max_zero_lag_wco"),("wco_mediana_zero","median_zero_lag_wco")):
            assert abs(metrics[key]-float(baseline[old_key])) < 1e-11
        with np.load(OLD/"intermediate"/f"{a}__{b}_maps.npz") as archive:
            idx = archive["lags"].astype(int)+j0
            errors = [float(abs(maps[key][:,idx]-archive[key]).max()) for key in ("wxx","wyy","wxy","wco")]
            assert max(errors) < 1e-11
            map_error = max(errors)
    step = max(1,int(np.ceil(len(maps["lags"])/1800)))
    idx = np.unique(np.r_[np.arange(0,len(maps["lags"]),step),j0,len(maps["lags"])-1])
    compact = {key:maps[key][:,idx] for key in ("wxx","wyy","wxy","wco")}
    compact.update(lags=maps["lags"][idx], wavelengths=WAVES, scales=maps["scales"], x=x,y=y,
                   mean_x_all=np.abs(maps["wxx"]).mean(axis=1),
                   mean_y_all=np.abs(maps["wyy"]).mean(axis=1),
                   max_wxx=np.abs(maps["wxx"]).max(), max_wyy=np.abs(maps["wyy"]).max(),
                   max_wxy=np.abs(maps["wxy"]).max(), wps_zero_x=zx, wps_zero_y=zy,
                   scale_peak_wco=peak, full_lag_count=len(maps["lags"]))
    np.savez_compressed(stem.with_suffix(".npz"), **compact)
    del maps,co
    layout = figuras.pair_figure([load_image(a),load_image(b)],x,y,compact,
                                 runs=(a,b),pattern=spec["padrao"],destination=ROOT/spec["arquivo"])
    result = dict(**spec, **metrics, erro_metricas_vs_original=max(old_errors),
                  erro_mapas_intra_vs_original=map_error, erro_soma_direta=max(direct_errors),
                  erro_wps_independente=max(independent_errors),
                  npz_sha256=sha(stem.with_suffix(".npz")), png_sha256=sha(ROOT/spec["arquivo"]),
                  assinaturas=signatures, figuras_sha256=sha(ROOT/"scripts/figuras.py"), apresentacao=layout)
    write_json(marker,result)
    return result


def run(workers, force=False):
    started = time.perf_counter()
    selected = audit_inputs()
    individual_rows = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(individual,row) for row in selected]
        for future in as_completed(futures):
            individual_rows.append(future.result())
            print(f"Panoramas e perfis: {len(individual_rows)}/26",flush=True)
    individual_rows.sort(key=lambda r: next(v["numero"] for v in selected if v["ensaio"]==r["ensaio"]))
    write_json(INTER/"individuais.json",individual_rows)
    tasks = specifications(selected)
    # Long titles and the smallest common length exercise both figure layouts first.
    preflight = [next(t for t in tasks if (t["run_x"],t["run_y"]) == ("m03q20_001","m03q40_001")),
                 next(t for t in tasks if (t["run_x"],t["run_y"]) == ("m03q20_001","m07q01_001"))]
    output = [calculate_pair(t, force=force) for t in preflight]
    done = {(t["run_x"],t["run_y"]) for t in preflight}
    print("Prévia de lote: par longo intrapadrão e par entre padrões verificados.",flush=True)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(calculate_pair,t,force) for t in tasks if (t["run_x"],t["run_y"]) not in done]
        for future in as_completed(futures):
            output.append(future.result())
            if len(output)%10 == 0 or len(output)==325:
                progress = dict(pares_concluidos=len(output),total=325,decorrido_s=round(time.perf_counter()-started,1))
                write_json(INTER/"execucao.json",progress)
                print(f"Pares: {len(output)}/325 · {progress['decorrido_s']} s",flush=True)
    output.sort(key=lambda r:(r["run_x"],r["run_y"]))
    write_json(INTER/"pares.json",output)
    layout = [r["apresentacao"] for r in individual_rows+output]
    write_json(INTER/"apresentacao/figuras.json",layout)
    flat_keys = [k for k in output[0] if k not in ("apresentacao","assinaturas")]
    write_csv(TABLES/"metricas_wavelets.csv",[{k:r[k] for k in flat_keys} for r in output])
    result = dict(status="PASS",ensaios=26,pares=325,intra=87,inter=238,
                  segundos=time.perf_counter()-started,
                  maior_erro_metricas_vs_original=max(r["erro_metricas_vs_original"] for r in output),
                  maior_erro_mapas_intra=max(r["erro_mapas_intra_vs_original"] for r in output if r["tipo"]=="intra"),
                  maior_erro_perfis=max(r["erro_perfil_vs_original"] for r in individual_rows),
                  maior_erro_soma_direta=max(r["erro_soma_direta"] for r in output),
                  maior_erro_wps_independente=max(r["erro_wps_independente"] for r in output),
                  menor_fonte_pt=min(r["fonte_minima_inserida_pt"] for r in layout),
                  maior_altura_cm=max(r["altura_inserida_cm"] for r in layout),
                  core_sha256=sha(ROOT/"scripts/wavelet_core.py"), script_sha256=sha(__file__))
    write_json(INTER/"auditoria/calculo.json",result)
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers",type=int,default=2)
    parser.add_argument("--force",action="store_true")
    args=parser.parse_args()
    assert 1<=args.workers<=4
    run(args.workers,args.force)
