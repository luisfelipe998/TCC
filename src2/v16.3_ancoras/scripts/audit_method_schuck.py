"""Read-only numerical/method audit; all new evidence stays in intermediario.

Reference: Schuck Jr. and Bodmann (2019), printed pp. 411, 413-414.
The independent path uses least squares, np.correlate and the temporal Morlet
convolution (Eq. 32.4), never the production CWT or its cached Fourier kernels.
Only the four pilot pairs and synthetic controls are evaluated. No real inter
pattern comparison is calculated, and no production configuration is changed.
"""
import sys
sys.dont_write_bytecode = True
import copy
import json
from datetime import datetime, timezone

import numpy as np
from PIL import Image
from scipy.signal import fftconvolve

from common import BASE, CFG, INTER, ROOT, SCRIPTS, image_path, read_csv, sha, write_csv, write_json
from prepare_comparisons import load_catalog, load_individual
from wavelet_cache import CACHE, context_record, grid_spec, signature

OUT = INTER / "auditoria_metodo_schuck"
FACTOR = 4*np.pi/(6+np.sqrt(38))


def independent_prepare(raw):
    signal = 255-np.asarray(raw, dtype=float)
    design = np.column_stack((np.ones(len(signal)), np.linspace(-1, 1, len(signal))))
    signal -= design @ np.linalg.lstsq(design, signal, rcond=None)[0]
    signal -= np.mean(signal)
    return signal / np.sqrt(np.dot(signal, signal)/len(signal))


def direct_correlation(x, y, unbiased=False):
    n = len(x)
    lags = np.arange(1-n, n)
    return np.correlate(x, y, "full") / (n-np.abs(lags) if unbiased else n)


def temporal_cwt(values, waves):
    """Full linear convolution with sampled Eq. 32.6; dx=1, Eq. 32.4."""
    offsets = np.arange(1-len(values), len(values))
    output = []
    for wavelength in waves:
        scale = wavelength / FACTOR
        u = offsets / scale
        kernel = np.pi**(-.25) / np.sqrt(scale) * np.exp(-u*u/2) * np.exp(6j*u)
        output.append(fftconvolve(values, kernel, mode="same"))
    return np.asarray(output)


def explicit_cwt(values, waves, positions):
    """Scalar products also check the independent convolution's indexing/sign."""
    index = np.arange(len(values))
    output = []
    for wavelength in waves:
        scale = wavelength / FACTOR
        row = []
        for position in positions:
            u = (index-position)/scale
            row.append(np.dot(values, np.pi**(-.25)/np.sqrt(scale)*np.exp(-u*u/2)*np.exp(-6j*u)))
        output.append(row)
    return np.asarray(output)


def independent_pair(x, y, waves, unbiased=False):
    transforms = [temporal_cwt(direct_correlation(a, b, unbiased), waves)
                  for a, b in ((x,x), (y,y), (x,y))]
    wxx, wyy, wxy = transforms
    center = len(x)-1
    a, b = np.abs(wxx[:, center]), np.abs(wyy[:, center])
    valid = (a > 1e-10*a.max()) & (b > 1e-10*b.max())
    coherence = np.full(wxy.shape, np.nan)
    coherence[valid] = np.abs(wxy[valid])**2/(a[valid]*b[valid])[:, None]
    return dict(wxx=wxx, wyy=wyy, wxy=wxy, wco=coherence, valid=valid)


def score(maps, n):
    return float(maps["wco"][maps["valid"], n-1-32:n+32].mean())


def max_error(a, b):
    np.testing.assert_allclose(a, b, atol=1e-10, rtol=1e-10, equal_nan=True)
    return float(np.nanmax(np.abs(a-b)))


def protected_files():
    files = [p for folder in (ROOT/"01_ancoras", ROOT/"02_wavelets", SCRIPTS, INTER/"mapas",
                              INTER/"preparacao_comparacoes") for p in folder.iterdir() if p.is_file()]
    files += [INTER/name for name in ("selecao.json", "congelamento.json", "metricas_wavelets.json")]
    files += list(CACHE.rglob("*.npy")) + list(CACHE.rglob("metadata.json"))
    return {str(p.relative_to(ROOT)): sha(p) for p in files}


def audit_individuals(catalog):
    results, metadata = [], []
    visited, images = {}, {}
    candidate_rows = read_csv(INTER/"preparacao_comparacoes/candidatos.csv")
    for run, entry in catalog["aquisicoes"].items():
        with Image.open(image_path(run)) as im:
            images[run] = np.asarray(im.convert("L"))
        with np.load(BASE/"intermediate"/f"{run}_input.npz") as data:
            widths = data["shifts"].astype(int)
            seams = data["seams"].astype(int)
            raw_shifts = np.r_[data["raw_shifts"], data["raw_shifts"][-1]]
        bounds = np.r_[0, np.cumsum(widths)]
        modified = seams[np.minimum(widths[:-1], widths[1:]) >= 4]-1
        for profile in ("nativo", "comparavel"):
            record, auto = load_individual(catalog, run, profile)
            spec = auto["meta"]["spec"]
            n, start, stop = record["comprimento"], record["segmento_inicio"], record["segmento_fim"]
            assert spec["aquisicao"] == run and spec["ancora"] == record["ancora"]
            assert spec["inicio"] == start and spec["fim"] == stop and spec["comprimento"] == n == stop-start
            assert spec["linhas"] == [18,71]
            assert spec["panorama_sha256"] == sha(image_path(run))
            # Independent geometric checks use the actual union of detector bounds.
            margin = 256 if record["tipo"] == "evento_delimitado" else 128
            assert start >= 16 and stop <= images[run].shape[1]-16
            assert record["inicio_suporte"]-start >= margin and stop-record["fim_suporte"] >= margin
            raw = images[run][18:72, start:stop].mean(axis=0)
            np.testing.assert_array_equal(raw, auto["raw"])
            first = np.searchsorted(bounds, start, side="right")-1
            last = np.searchsorted(bounds, stop-1, side="right")-1
            neighbors = [r for r in candidate_rows if r["aquisicao"] == run
                         and r["anchor_id"] != record["anchor_id"] and r["elegivel_individual"] == "True"
                         and (float(r["inicio"]) < stop and float(r["fim"]) > start if record["tipo"] == "evento_delimitado"
                              else start <= int(r["ancora"]) < stop)]
            fresh = dict(costuras_contexto=int(((seams>start)&(seams<stop)).sum()),
                         colunas_substituidas_feather_contexto=int(((modified>=start)&(modified<stop)).sum()),
                         deslocamentos_nao_positivos_contexto=int((raw_shifts[first:last+1]<=0).sum()),
                         outros_candidatos_estaveis_no_contexto=len(neighbors))
            for name, actual in fresh.items():
                if record[name] != actual:
                    metadata.append(dict(aquisicao=run, perfil=profile, comprimento=n, campo=name,
                                         armazenado=record[name], recalculado=actual))
            if auto["key"] in visited:
                continue
            prepared = independent_prepare(raw)
            error_prepare = max_error(prepared, auto["x"])
            corr = direct_correlation(prepared, prepared)
            error_corr = max_error(corr, auto["rxx"])
            waves = np.geomspace(50, n/2, 200)
            transformed = temporal_cwt(corr, waves)
            error_cwt = max_error(transformed, auto["wxx"])
            lags = np.arange(1-n, n)
            support = n-1-np.abs(lags)[None,:] >= np.sqrt(2)*(waves/FACTOR)[:,None]
            mean = np.asarray([np.mean(np.abs(t)[s]) for t,s in zip(transformed,support)])
            error_mean = max_error(mean, auto["power_mean"])
            center = np.abs(transformed[:,n-1])
            row = dict(aquisicao=run, comprimento=n, cache_key=auto["key"],
                       erro_preparo=error_prepare, erro_correlacao=error_corr,
                       erro_cwt_temporal=error_cwt, erro_potencia_media=error_mean,
                       escalas_denominador_valido=int((center>1e-10*center.max()).sum()),
                       fracao_variancia_retirada_pela_tendencia=float(1-np.var(255-raw-np.polyval(np.polyfit(np.arange(n),255-raw,1),np.arange(n)))/np.var(raw)),
                       **fresh)
            visited[auto["key"]] = True
            results.append(row)
            print(f"Individual conferido: {run}, N={n}, erro CWT={error_cwt:.2e}", flush=True)
    write_csv(OUT/"individuais.csv", results)
    write_csv(OUT/"diagnosticos_desatualizados.csv", metadata)
    return results, metadata


def audit_pilot(catalog):
    selection = json.loads((INTER/"selecao.json").read_text())
    results, estimators, profiles = [], [], []
    for comparison in selection["pares"]:
        n = comparison["comprimento"]
        waves = np.geomspace(50, n/2, 200)
        _, a = load_individual(catalog, comparison["run_x"], "nativo")
        _, b = load_individual(catalog, comparison["run_y"], "nativo")
        # Prepare each signal independently from the stored raw values already
        # checked against the original panorama in audit_individuals.
        x, y = independent_prepare(a["raw"]), independent_prepare(b["raw"])
        direct = independent_pair(x, y, waves)
        filename = f"{comparison['padrao'].replace(' ','_')}_{comparison['run_x']}_{comparison['run_y']}.npz"
        with np.load(INTER/"mapas"/filename) as stored:
            errors = {key:max_error(direct[key], stored[key]) for key in ("wxx", "wyy", "wxy", "wco")}
            np.testing.assert_array_equal(direct["valid"], stored["denominator_valid"])
        indexes = np.array([0,43,99,149,199])
        positions = np.array([0,n-1-32,n-1,n-1+32,2*n-2])
        explicit = explicit_cwt(direct_correlation(x,y), waves[indexes], positions)
        error_explicit = max_error(explicit, direct["wxy"][np.ix_(indexes,positions)])
        reverse = temporal_cwt(direct_correlation(y,x),waves[indexes])
        error_reverse = max_error(reverse, np.conj(direct["wxy"][indexes,::-1]))
        biased_score = score(direct,n)
        assert np.nanmax(direct["wco"]) <= 1+1e-10
        results.append(dict(padrao=comparison["padrao"], comprimento=n,
                            **{f"erro_{k}":v for k,v in errors.items()},
                            erro_produto_escalar=error_explicit, erro_troca_xy=error_reverse,
                            wco_media=biased_score, max_wco=float(np.nanmax(direct["wco"]))))
        literal = independent_pair(x,y,waves,unbiased=True)
        lags = np.arange(1-n,n)
        support = n-1-np.abs(lags)[None,:] >= np.sqrt(2)*(waves/FACTOR)[:,None]
        estimators.append(dict(padrao=comparison["padrao"], comprimento=n, media_divisor_N=biased_score,
                               media_divisor_N_menos_abs_lag=score(literal,n),
                               diferenca=score(literal,n)-biased_score,
                               maximo_literal=float(np.nanmax(literal["wco"])),
                               maximo_literal_com_suporte=float(np.nanmax(np.where(support,literal["wco"],np.nan))),
                               celulas_literal_acima_1=int(np.sum(literal["wco"]>1+1e-10))))
        _, ac = load_individual(catalog, comparison["run_x"], "comparavel")
        _, bc = load_individual(catalog, comparison["run_y"], "comparavel")
        common_n = len(ac["x"])
        common_waves = np.geomspace(50, common_n/2,200)
        if common_n == n:
            common_score = restricted_score = biased_score
        else:
            reduced = independent_pair(independent_prepare(ac["raw"]),independent_prepare(bc["raw"]),common_waves)
            common_score = score(reduced,common_n)
            same_band = independent_pair(x,y,common_waves)
            restricted_score = score(same_band,n)
        profiles.append(dict(padrao=comparison["padrao"], janela_nativa=n, janela_comum=common_n,
                             media_nativa=biased_score, media_nativa_banda_50_1024=restricted_score,
                             media_comparavel=common_score,
                             diferenca_somente_janela_mesma_banda=common_score-restricted_score))
        print(f"Par piloto: {comparison['padrao']}, erro independente WCO={errors['wco']:.2e}; "
              f"N={biased_score:.4f}, N-|m|={score(literal,n):.4f}, perfil comum={common_score:.4f}",flush=True)
    write_csv(OUT/"pares_piloto.csv",results)
    write_csv(OUT/"sensibilidade_estimador.csv",estimators)
    write_csv(OUT/"sensibilidade_perfil.csv",profiles)
    return results, estimators, profiles


def audit_guards(catalog):
    baseline = signature(grid_spec(2048))
    old = copy.deepcopy(CFG["wavelet"])
    changes = []
    try:
        for name,value in (("estimador","unbiased"),("K",2),("omega0",8)):
            CFG["wavelet"][name] = value
            changes.append(dict(parametro=name, valor_testado=value,
                                assinatura_permanece_igual=signature(grid_spec(2048))==baseline))
            CFG["wavelet"].update(old)
    finally:
        CFG["wavelet"].update(old)
    record,_ = load_individual(catalog,"m05q01_001","comparavel")
    odd = context_record(record,2049)
    changed = copy.deepcopy(catalog)
    changed["aquisicoes"]["m05q01_001"]["perfis"]["comparavel"] = copy.deepcopy(
        changed["aquisicoes"]["m06q01_001"]["perfis"]["comparavel"])
    bad_record,bad_auto = load_individual(changed,"m05q01_001","comparavel")
    # This deliberately altered dictionary is never written back to disk.
    return dict(parametros_nao_aplicados=changes, janela_impar_solicitada=2049,
                amostras_janela_impar=odd["segmento_fim"]-odd["segmento_inicio"],
                vinculo_errado_no_indice_aceito=bad_auto["meta"]["spec"]["aquisicao"]!="m05q01_001",
                arquivo_indice_original_preservado=True)


def synthetic_controls():
    rng = np.random.default_rng(20261025)
    n = 512
    waves = np.geomspace(50,256,40)
    x = rng.normal(size=n)
    x -= x.mean()
    identity = independent_pair(x,x,waves)
    inverted = independent_pair(x,-x,waves)
    identity_error = float(np.max(abs(identity["wco"][:,n-1]-1)))
    inversion_error = max_error(identity["wco"],inverted["wco"])
    assert identity_error < 1e-12
    impulse_x, impulse_y = np.zeros(n), np.zeros(n)
    impulse_x[180], impulse_y[217] = 1,1
    delay = independent_pair(impulse_x,impulse_y,np.array([80.]))
    measured = int(np.nanargmax(delay["wco"][0])-(n-1))
    assert measured == -37
    # Correlation/WCO cannot locate an event after a joint translation of both
    # records when the support remains entirely inside the finite window.
    shifted_x, shifted_y = np.zeros(n), np.zeros(n)
    shifted_x[280], shifted_y[317] = 1,1
    shifted = independent_pair(shifted_x,shifted_y,np.array([80.]))
    shift_error = max_error(delay["wco"],shifted["wco"])
    # Illustration of finite-record coherence under independent noise. It is
    # not a calibrated null for the flow records or their anchor selection.
    waves_null = np.array([50.,100.,250.,500.,1024.])
    n_null = 2048
    values = []
    for trial in range(128):
        a,b = [independent_prepare(128+rng.normal(size=n_null)) for _ in range(2)]
        aa,bb,ab = [explicit_cwt(direct_correlation(u,v),waves_null,[n_null-1])[:,0]
                    for u,v in ((a,a),(b,b),(a,b))]
        values.append(abs(ab)**2/(abs(aa)*abs(bb)))
    values = np.asarray(values)
    null_rows = [dict(comprimento_onda=float(w), media=float(v.mean()), p95=float(np.quantile(v,.95)),
                     fracao_acima_05=float(np.mean(v>.5))) for w,v in zip(waves_null,values.T)]
    write_csv(OUT/"ruido_independente_ilustrativo.csv",null_rows)
    return dict(erro_identidade_lag0=identity_error, erro_inversao_global=inversion_error,
                lag_impulso_atrasado_37=measured, erro_translacao_conjunta=shift_error,
                nulo_ilustrativo=dict(semente=20261025, realizacoes=128, comprimento=2048,
                    dominio="Somente lag zero; sem selecao de ancora. Nao e teste de significancia dos dados reais.",
                    resultados=null_rows),
                ensemble_fases_opostas=dict(modulo_media_complexa=0.,media_modulos=float(abs(1+2j))))


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    protected = protected_files()
    write_json(OUT/"estado_antes.json", protected)
    catalog = load_catalog()
    individuals, metadata = audit_individuals(catalog)
    pairs, estimators, profiles = audit_pilot(catalog)
    controls = synthetic_controls()
    guards = audit_guards(catalog)
    assert protected == protected_files(), "A auditoria modificou um artefato protegido."
    assert not list((ROOT/"03_conclusoes").iterdir())
    report = dict(data_utc=datetime.now(timezone.utc).isoformat(),
        status_numerico="PASS", status_metodo="RESSALVAS_E_CORRECOES_NECESSARIAS",
        individuais_conferidos=len(individuals), aquisicoes=len(catalog["aquisicoes"]),
        pares_piloto_conferidos=len(pairs), pares_inter_calculados=0,
        todas_as_celulas_wps_conferidas=True, todas_as_celulas_wco_piloto_conferidas=True,
        erros_maximos_individuais={k:max(r[k] for r in individuals) for k in
            ("erro_preparo","erro_correlacao","erro_cwt_temporal","erro_potencia_media")},
        erros_maximos_pares={k:max(r[k] for r in pairs) for k in
            ("erro_wxx","erro_wyy","erro_wxy","erro_wco","erro_produto_escalar","erro_troca_xy")},
        diagnosticos_desatualizados=len(metadata),
        perfis_com_diagnosticos_desatualizados=len({(r['aquisicao'],r['perfil']) for r in metadata}),
        sensibilidade_estimador=estimators, sensibilidade_perfil=profiles,
        controles=controls, protecoes=guards,
        artefatos_atuais_preservados=True,
        fonte=dict(arquivo="assets/schuck-bardo-wavelet-coherence.pdf",
                   sha256=sha(ROOT.parent.parent/"assets/schuck-bardo-wavelet-coherence.pdf"),
                   paginas_impressas_conferidas_visualmente=[411,413,414],
                   equacoes=["32.4","32.6","32.8","32.14-32.23"]),
        script_sha256=sha(SCRIPTS/"audit_method_schuck.py"))
    write_json(OUT/"auditoria.json",report)
    print(json.dumps({k:report[k] for k in ("status_numerico","status_metodo","individuais_conferidos",
        "erros_maximos_individuais","erros_maximos_pares","perfis_com_diagnosticos_desatualizados",
        "pares_inter_calculados")},indent=2),flush=True)


if __name__ == "__main__":
    main()
