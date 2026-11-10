"""Compute only the four frozen pilot pairs using the unchanged v16.2 core."""
import json
from datetime import datetime, timezone
import numpy as np

from common import (CFG, INTER, OLD, SCRIPTS, load_image, selected_manifest,
                    sha, slug, write_csv, write_json)
from wavelet_cache import pair_from_autos, PREPARATION
from prepare_comparisons import load_catalog, load_individual
from pair_identity import pair_key


def calculate():
    selection = selected_manifest()
    assert sha(SCRIPTS / "wavelet_core.py") == sha(OLD / "scripts/wavelet_core.py")
    assert len(selection["pares"]) <= 4
    catalog = load_catalog()
    inspected = json.loads((INTER/"inspecao_visual_piloto.json").read_text())
    assert inspected["antes_das_wavelets"]
    assert {r["anchor_id"] for r in inspected["ancoras"] if r["verificada_visualmente"]} == {r["anchor_id"] for r in selection["ancoras"]}
    anchors = {r["anchor_id"]: r for r in selection["ancoras"]}
    metrics = []
    for comparison in selection["pares"]:
        records = [anchors[comparison[k]] for k in ("anchor_x", "anchor_y")]
        length = comparison["comprimento"]
        autos = []
        for record in records:
            context, auto = load_individual(catalog, record["aquisicao"], "nativo")
            for field in ("anchor_id", "ancora", "comprimento", "segmento_inicio", "segmento_fim"):
                if context[field] != record[field]:
                    raise ValueError("O perfil individual nao corresponde ao par congelado do piloto.")
            autos.append(auto)
        maps = pair_from_autos(*autos)
        wavelengths = maps["wavelengths"]
        lags, scales = maps["lags"], maps["scales"]
        support = (length - 1 - np.abs(lags)[None, :]) >= np.sqrt(2) * scales[:, None]
        domain = np.abs(lags) <= CFG["wavelet"]["lag_resumo"]
        assert support[:, domain].all()
        valid = maps["denominator_valid"]
        assert valid.any(), comparison
        values = maps["wco"][valid][:, domain]
        summary = dict(padrao=comparison["padrao"], run_x=comparison["run_x"], run_y=comparison["run_y"],
                       anchor_x=comparison["anchor_x"], anchor_y=comparison["anchor_y"], comprimento=length,
                       lambda_min=float(wavelengths[0]), lambda_max=float(wavelengths[-1]),
                       escalas_validas=int(valid.sum()), lags_resumo=int(domain.sum()),
                       wco_media=float(np.mean(np.mean(values, axis=1))),
                       wco_mediana=float(np.median(values)),
                       fracao_wco_acima_05=float(np.mean(values >= .5)),
                       fracao_wco_acima_07=float(np.mean(values >= .7)),
                       wco_maxima_suporte=float(np.nanmax(np.where(support, maps["wco"], np.nan))),
                       wco_maxima_integral=float(np.nanmax(maps["wco"])),
                       razao_larguras=max(r["largura"] for r in records) / min(r["largura"] for r in records))
        meta = dict(comparacao=comparison, metricas=summary,
                    congelamento_sha256=sha(INTER / "congelamento.json"),
                    nucleo_sha256=sha(SCRIPTS / "wavelet_core.py"),
                    preparo=PREPARATION, cache_individual_x=autos[0]["key"], cache_individual_y=autos[1]["key"],
                    backend_sha256=sha(SCRIPTS / "wavelet_cache.py"),
                    preparacao_sha256=catalog["congelamento_sha256"],
                    calculo_sha256=sha(SCRIPTS/"compute_wavelets.py"),
                    potencia_media=CFG["wavelet"]["potencia_media"],
                    data_utc=datetime.now(timezone.utc).isoformat())
        np.savez_compressed(INTER / "mapas" / f"{pair_key(comparison)}.npz", **maps,
                            meta_json=np.array(json.dumps(meta, ensure_ascii=False)))
        metrics.append(summary)
        from plots import pair_figure
        pair_figure(comparison, anchors, maps, summary)
        print(f"{comparison['padrao']}: N={length}; WCO media={summary['wco_media']:.4f}; "
              f"fracao >= 0,5: {summary['fracao_wco_acima_05']:.1%}", flush=True)
    write_csv(INTER / "metricas_wavelets.csv", metrics)
    write_json(INTER / "metricas_wavelets.json", metrics)
    return metrics


def render_final():
    from plots import pair_figure
    selection = selected_manifest()
    labels = json.loads((INTER / "rotulos_visuais.json").read_text())
    metrics = {r["padrao"]: r for r in json.loads((INTER / "metricas_wavelets.json").read_text())}
    anchors = {r["anchor_id"]: r for r in selection["ancoras"]}
    for comparison in selection["pares"]:
        pattern = comparison["padrao"]
        label = labels[pattern]
        assert label["rotulo"] in ("favoravel", "desfavoravel")
        assert label["run_x"] == comparison["run_x"] and label["run_y"] == comparison["run_y"]
        path = INTER / "mapas" / f"{pair_key(comparison)}.npz"
        if "mapa_sha256" in label:
            assert sha(path) == label["mapa_sha256"]
        with np.load(path) as stored:
            maps = {key: stored[key].copy() for key in stored.files if key != "meta_json"}
        pair_figure(comparison, anchors, maps, metrics[pattern], label=label["rotulo"])
        print(f"Figura final: {pattern} / {label['rotulo']}", flush=True)
