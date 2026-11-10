"""Expand the frozen individual catalogue to all intra-pattern pairs only.

Anchor figures use the native context. All pair figures use the shared 2048 px
profile, also prepared for future inter-pattern work. Pair archives reference
individual WPS arrays rather than duplicating them. Visual labels are reviewed
separately and never affect anchor selection or numerical results.
"""
import sys
sys.dont_write_bytecode = True
import argparse
import json
import math
import shutil
import time
from collections import Counter
from datetime import datetime, timezone
from itertools import combinations

import numpy as np

from common import (ANCHORS, CFG, CONCLUSIONS, INTER, ROOT, SCRIPTS, WAVELETS,
                    load_image, sha, slug, write_csv, write_json)
from pair_identity import pair_key
from prepare_comparisons import load_catalog, comparison_spec, load_individual
from wavelet_cache import CACHE, get_grid, pair_from_autos

OUT = INTER / "lote_intra"
MAPS = OUT / "mapas"
PROFILE = "comparavel"


def now():
    return datetime.now(timezone.utc).isoformat()


def comparisons(catalog):
    if catalog["escopo"] != "completo":
        raise ValueError("Prepare o catalogo completo antes do lote intra padrao.")
    for pattern in CFG["padroes"]:
        runs = sorted(run for run, item in catalog["aquisicoes"].items() if item["padrao"] == pattern)
        for a, b in combinations(runs, 2):
            spec, records, autos = comparison_spec(catalog, a, b, PROFILE)
            assert spec["escopo"] == "intra" and spec["comprimento"] == 2048
            yield spec, records, autos


def setup_plots():
    import plots
    plots.INTER = OUT
    return plots


def prepare():
    from prepare_comparisons import prepare_all
    if (OUT/"inicio.json").exists():
        raise ValueError("Lote ja iniciado; continue pela etapa pendente para preservar seu congelamento.")
    archive = OUT/"piloto_anterior"
    archive.mkdir(exist_ok=True)
    for source in (INTER/"preparacao_comparacoes", ANCHORS, WAVELETS):
        shutil.copytree(source, archive/source.name)
    write_json(OUT/"inicio.json", dict(data_utc=now(), cache_preexistente={
        str(p.relative_to(CACHE)):dict(mtime_ns=p.stat().st_mtime_ns, sha256=sha(p))
        for p in CACHE.rglob("*") if p.is_file()}))
    catalog = prepare_all("completo")
    return dict(aquisicoes=len(catalog["aquisicoes"]), intra=catalog["pares_intra_planejados"], inter_calculados=0)


def anchors():
    catalog = load_catalog()
    plots = setup_plots()
    records = []
    for pattern in CFG["padroes"]:
        for run, item in sorted(catalog["aquisicoes"].items()):
            if item["padrao"] != pattern:
                continue
            record = item["perfis"]["nativo"]["contexto"]
            pano, raw = load_image(run)
            plots.anchor_figure(record, pano, raw)
            path = ANCHORS / f"{slug(pattern)}_{record['anchor_id']}.png"
            records.append(dict(padrao=pattern, aquisicao=run, anchor_id=record["anchor_id"],
                                ancora=record["ancora"], comprimento=record["comprimento"],
                                arquivo=str(path.relative_to(ROOT)), png_sha256=sha(path)))
            print(f"Ancora {len(records)}: {pattern} / {run} / coluna {record['ancora']} px", flush=True)
    write_json(OUT / "figuras_ancoras.json", records)
    # Compact scientific contact sheets for inspection, kept out of final folders.
    gallery = OUT / "galerias_ancoras"
    gallery.mkdir(exist_ok=True)
    import matplotlib.pyplot as plt
    for page, offset in enumerate(range(0, len(records), 3), 1):
        subset = records[offset:offset+3]
        fig, axs = plt.subplots(1, 3, figsize=(15, 7.5))
        for ax in axs:
            ax.set_axis_off()
        for ax, record in zip(axs, subset):
            path = OUT / "apresentacao" / "previas_16cm" / (ROOT / record["arquivo"]).name
            ax.imshow(plt.imread(path))
        fig.subplots_adjust(left=.003, right=.997, bottom=.003, top=.997, wspace=.02)
        fig.savefig(gallery / f"ancoras_{page:02d}.png", dpi=140)
        plt.close(fig)
    return dict(ancoras=len(records))


def summary(spec, records, maps):
    domain = np.abs(maps["lags"]) <= CFG["wavelet"]["lag_resumo"]
    valid = maps["denominator_valid"]
    assert maps["support"][:, domain].all() and valid.any()
    values = maps["wco"][valid][:, domain]
    assert np.nanmin(maps["wco"]) >= 0 and np.nanmax(maps["wco"]) <= 1+1e-12
    return dict(**spec, lambda_min=float(maps["wavelengths"][0]),
                lambda_max=float(maps["wavelengths"][-1]), escalas_validas=int(valid.sum()),
                lags_resumo=int(domain.sum()), wco_media=float(values.mean()),
                wco_mediana=float(np.median(values)), fracao_wco_acima_05=float((values>=.5).mean()),
                fracao_wco_acima_07=float((values>=.7).mean()),
                wco_maxima_suporte=float(np.nanmax(np.where(maps["support"], maps["wco"], np.nan))),
                wco_maxima_integral=float(np.nanmax(maps["wco"])),
                razao_larguras=max(r["largura"] for r in records.values())/min(r["largura"] for r in records.values()))


def calculate():
    catalog = load_catalog()
    inspected = json.loads((OUT / "inspecao_ancoras.json").read_text())
    assert inspected["antes_das_comparacoes"]
    assert inspected["congelamento_sha256"] == catalog["congelamento_sha256"]
    expected = {x["anchor_id"] for x in catalog["aquisicoes"].values()}
    assert {r["anchor_id"] for r in inspected["ancoras"] if r["verificada_visualmente"]} == expected
    for record in inspected["ancoras"]:
        assert sha(ROOT / record["arquivo"]) == record["png_sha256"]
    MAPS.mkdir(exist_ok=True)
    rows = []
    for spec, records, autos in comparisons(catalog):
        maps = pair_from_autos(*autos)
        stats = summary(spec, records, maps)
        meta = dict(schema=1, armazenamento="WCS/WCO do par; WPS e sinais referenciados no cache individual",
                    comparacao=spec, metricas=stats, cache_individual_x=autos[0]["key"],
                    cache_individual_y=autos[1]["key"], preparacao_sha256=catalog["congelamento_sha256"],
                    nucleo_sha256=sha(SCRIPTS/"wavelet_core.py"), backend_sha256=sha(SCRIPTS/"wavelet_cache.py"),
                    calculo_sha256=sha(__file__), data_utc=now())
        path = MAPS / f"{pair_key(spec)}.npz"
        temp = path.with_suffix(".tmp.npz")
        np.savez_compressed(temp, wxy=maps["wxy"], wco=maps["wco"],
                            denominator_valid=maps["denominator_valid"],
                            meta_json=np.array(json.dumps(meta, ensure_ascii=False)))
        temp.replace(path)
        rows.append(dict(**stats, arquivo=str(path.relative_to(ROOT)), mapa_sha256=sha(path)))
        print(f"Par {len(rows)}/{catalog['pares_intra_planejados']}: {pair_key(spec)}", flush=True)
    assert len(rows) == catalog["pares_intra_planejados"]
    write_json(OUT / "metricas_wavelets.json", rows)
    write_csv(OUT / "metricas_wavelets.csv", rows)
    return dict(pares=len(rows), pares_inter=0)


def read_pair(catalog, row):
    spec, records, (a, b) = comparison_spec(catalog, row["run_x"], row["run_y"], PROFILE)
    path = ROOT / row["arquivo"]
    assert sha(path) == row["mapa_sha256"]
    grid = get_grid(spec["comprimento"])
    with np.load(path) as stored:
        meta = json.loads(str(stored["meta_json"]))
        assert meta["comparacao"] == spec
        assert meta["preparacao_sha256"] == catalog["congelamento_sha256"]
        assert meta["cache_individual_x"] == a["key"] and meta["cache_individual_y"] == b["key"]
        assert meta["nucleo_sha256"] == sha(SCRIPTS/"wavelet_core.py")
        assert meta["backend_sha256"] == sha(SCRIPTS/"wavelet_cache.py")
        assert meta["calculo_sha256"] == sha(__file__)
        maps = {key:stored[key].copy() for key in ("wxy", "wco", "denominator_valid")}
    maps.update(wxx=a["wxx"], wyy=b["wxx"], x=a["x"], y=b["x"], raw_x=a["raw"], raw_y=b["raw"],
                power_mean_x=a["power_mean"], power_mean_y=b["power_mean"],
                **{key:grid[key] for key in ("lags", "scales", "wavelengths", "support")})
    return spec, records, maps


def galleries():
    import matplotlib.pyplot as plt
    catalog = load_catalog()
    rows = json.loads((OUT / "metricas_wavelets.json").read_text())
    folder = OUT / "galerias_wco"
    folder.mkdir(exist_ok=True)
    for pattern in CFG["padroes"]:
        subset = [r for r in rows if r["padrao"] == pattern]
        for page, offset in enumerate(range(0, len(subset), 9), 1):
            group = subset[offset:offset+9]
            fig, axs = plt.subplots(3, 3, figsize=(15, 9), layout="constrained")
            for ax in axs.flat:
                ax.set_visible(False)
            for position, (ax, row) in enumerate(zip(axs.flat, group), offset+1):
                _, _, maps = read_pair(catalog, row)
                ax.set_visible(True)
                # Consistent full lag domain, logarithmic wavelength and 0..1 scale.
                indexes = np.unique(np.r_[np.arange(0, len(maps["lags"]), 4), len(maps["x"])-1])
                mesh = ax.pcolormesh(maps["lags"][indexes], maps["wavelengths"], maps["wco"][:, indexes],
                                     shading="auto", cmap="magma", vmin=0, vmax=1)
                ax.set_yscale("log")
                ax.set_ylim(50, 1024)
                ax.set_yticks([64, 256, 1024], ["64", "256", "1024"])
                for sign in (-1, 1):
                    ax.axvline(sign*32, color="#7f939e", lw=.7, ls="--")
                    ax.plot(sign*(2047-np.sqrt(2)*maps["scales"]), maps["wavelengths"], color="#888888", lw=.7, ls=":")
                ax.set_title(f"{position:02d} · {row['run_x']} × {row['run_y']}", fontsize=9)
                ax.set_xlabel("Lag (px)")
                ax.set_ylabel("λ (px)")
                fig.colorbar(mesh, ax=ax, ticks=[0, .5, 1], fraction=.025)
            fig.suptitle(f"{pattern} · WCO · 2.048 px · inspeção {page}", fontsize=13)
            fig.savefig(folder / f"{slug(pattern)}_{page:02d}.png", dpi=135)
            plt.close(fig)
    return dict(pares=len(rows))


def render():
    catalog = load_catalog()
    rows = json.loads((OUT / "metricas_wavelets.json").read_text())
    labels = json.loads((OUT / "rotulos_visuais.json").read_text())
    assert set(labels) == {pair_key(row) for row in rows}
    plots = setup_plots()
    destination = OUT / "figuras_prontas"
    destination.mkdir(exist_ok=True)
    plots.WAVELETS = destination
    figures = []
    for i, row in enumerate(rows, 1):
        spec, records, maps = read_pair(catalog, row)
        key = pair_key(spec)
        label = labels[key]
        assert label["rotulo"] in ("favoravel", "desfavoravel")
        assert label["mapa_sha256"] == row["mapa_sha256"]
        # Final filenames follow the requested pattern/runs/label. The explicit
        # profile remains in the numerical archive key and figure manifest.
        visible_spec = dict(spec)
        visible_spec.pop("perfil")
        plots.pair_figure(visible_spec, records, maps, row, label=label["rotulo"])
        name = f"{pair_key(visible_spec)}_{label['rotulo']}.png"
        figures.append(dict(pair_id=key, arquivo=name, perfil=PROFILE, comprimento=spec["comprimento"],
                            mapa_sha256=row["mapa_sha256"], png_sha256=sha(destination/name)))
        print(f"PNG {i}/{len(rows)}: {name}", flush=True)
    assert set(p.name for p in destination.glob("*.png")) == {r["arquivo"] for r in figures}
    # Install only after the entire new batch rendered successfully. Keep the
    # previous native pilot figures in intermediario, including its 4096 px pairs.
    archive = OUT / "piloto_anterior" / "02_wavelets"
    archive.mkdir(parents=True, exist_ok=True)
    for path in WAVELETS.glob("*.png"):
        if not (archive/path.name).exists():
            shutil.copy2(path, archive/path.name)
        path.unlink()
    for row in figures:
        (destination/row["arquivo"]).replace(WAVELETS/row["arquivo"])
    write_json(OUT / "figuras_wavelets.json", figures)
    return dict(imagens=len(figures))


def verify():
    from audit_method_schuck import independent_prepare, independent_pair, max_error
    import verify_rebuild
    catalog = load_catalog()
    verification = OUT / "verificacao"
    verification.mkdir(exist_ok=True)
    verify_rebuild.OUT = verification
    individuals, diagnostics = verify_rebuild.verify_individuals(catalog)
    rows = json.loads((OUT / "metricas_wavelets.json").read_text())
    pairs = []
    expected_pairs = {(spec["run_x"], spec["run_y"]) for spec, _, _ in comparisons(catalog)}
    actual_pairs = {(row["run_x"], row["run_y"]) for row in rows}
    assert actual_pairs == expected_pairs and len(rows) == len(actual_pairs) == 87
    assert len(list(MAPS.glob("*.npz"))) == len(rows)
    for i, row in enumerate(rows, 1):
        spec, records, maps = read_pair(catalog, row)
        x, y = independent_prepare(maps["raw_x"]), independent_prepare(maps["raw_y"])
        expected = independent_pair(x, y, np.geomspace(50, 1024, 200))
        errors = {f"erro_{key}":max_error(expected[key], maps[key]) for key in ("wxx", "wyy", "wxy", "wco")}
        np.testing.assert_array_equal(expected["valid"], maps["denominator_valid"])
        repeated = summary(spec, records, maps)
        for key, value in repeated.items():
            assert row[key] == value, (key, row[key], value)
        pairs.append(dict(pair_id=pair_key(spec), **errors))
        if i % 10 == 0:
            print(f"Conferencia independente: {i}/{len(rows)} pares", flush=True)
    write_csv(verification / "pares.csv", pairs)
    baseline = json.loads((OUT/"inicio.json").read_text())
    for relative, state in baseline["cache_preexistente"].items():
        path = CACHE / relative
        assert path.stat().st_mtime_ns == state["mtime_ns"] and sha(path) == state["sha256"]
    before = {str(p):p.stat().st_mtime_ns for p in CACHE.rglob("*") if p.is_file()}
    # All future pair bindings are checked, but only intra WCS/WCO were calculated.
    from wavelet_cache import get_auto
    for run, item in catalog["aquisicoes"].items():
        for profile in ("nativo", "comparavel"):
            record, auto = load_individual(catalog, run, profile)
            assert get_auto(record)["key"] == auto["key"]
    compatibility = Counter()
    for a, b in combinations(sorted(catalog["aquisicoes"]), 2):
        spec, _, autos = comparison_spec(catalog, a, b, PROFILE)
        compatibility[spec["escopo"]] += 1
        assert autos[0]["meta"]["spec"]["grade"] == autos[1]["meta"]["spec"]["grade"]
    assert dict(compatibility) == dict(intra=87, inter=213)
    assert before == {str(p):p.stat().st_mtime_ns for p in CACHE.rglob("*") if p.is_file()}
    figures = json.loads((OUT/"figuras_wavelets.json").read_text())
    anchor_figures = json.loads((OUT/"figuras_ancoras.json").read_text())
    assert len(figures) == 87 and len(anchor_figures) == 25
    assert {p.name for p in WAVELETS.glob("*.png")} == {r["arquivo"] for r in figures}
    assert {str(p.relative_to(ROOT)) for p in ANCHORS.glob("*.png")} == {r["arquivo"] for r in anchor_figures}
    presentation = []
    for folder, entries in ((WAVELETS, figures), (ROOT, anchor_figures)):
        for row in entries:
            path = folder / row["arquivo"]
            assert sha(path) == row["png_sha256"]
            record = json.loads((OUT/"apresentacao"/f"{path.stem}.json").read_text())
            assert record["altura_inserida_cm"] <= 24 and record["fonte_minima_inserida_pt"] >= 8
            assert "favor" not in record["titulo"].lower()
            assert "\n" not in record["titulo"]
            presentation.append(record)
    assert not list(CONCLUSIONS.iterdir()) and not list(ROOT.rglob("*.pdf"))
    counts = Counter(item["padrao"] for item in catalog["aquisicoes"].values())
    report = dict(status="PASS", data_utc=now(), ancoras=25, pares_intra=87, pares_inter_calculados=0,
                  por_padrao={p:dict(ensaios=n, pares=math.comb(n, 2)) for p,n in counts.items()},
                  perfil_comparacoes=PROFILE, comprimento_comparacoes=2048,
                  individuais_conferidos=len(individuals), perfis_conferidos=len(diagnostics),
                  grades=catalog["grades"], compatibilidade_futura_somente_metadados=dict(compatibility),
                  cache_do_piloto_preservado=True, cache_reutilizado_sem_reescrita=True,
                  erros_maximos_individuais={k:max(r[k] for r in individuals) for k in ("erro_preparo","erro_correlacao","erro_cwt","erro_potencia_media")},
                  erros_maximos_pares={k:max(r[k] for r in pairs) for k in ("erro_wxx","erro_wyy","erro_wxy","erro_wco")},
                  imagens_png=len(presentation), fonte_minima_pt=min(r["fonte_minima_inserida_pt"] for r in presentation),
                  altura_maxima_16cm=max(r["altura_inserida_cm"] for r in presentation), conclusoes_vazias=True)
    write_json(verification/"verificacao.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("etapa", choices=["preparar", "ancoras", "calcular", "galerias", "figuras", "verificar"])
    args = parser.parse_args()
    OUT.mkdir(exist_ok=True)
    start = time.perf_counter()
    result = {"preparar":prepare, "ancoras":anchors, "calcular":calculate, "galerias":galleries, "figuras":render, "verificar":verify}[args.etapa]()
    result.update(etapa=args.etapa, segundos=time.perf_counter()-start, data_utc=now())
    write_json(OUT/f"execucao_{args.etapa}.json", result)
    print(json.dumps(result, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
