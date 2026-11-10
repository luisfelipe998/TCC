"""All cross-pattern pairs, with the frozen common context and existing WPS.

No resizing, padding of measured signals, or pair-dependent anchor selection.
Numerical archives retain the profile in their identity. Final PNG names follow
the requested pattern/pattern/run/run/visual-label convention.
"""
import sys
sys.dont_write_bytecode = True
import argparse
import json
from collections import Counter
from itertools import combinations, product

import numpy as np

from common import CFG, INTER, ROOT, SCRIPTS, image_path, sha, slug, write_csv, write_json
from pair_identity import pair_key
from prepare_comparisons import comparison_spec, load_catalog, load_individual
from wavelet_cache import CACHE, get_grid, pair_from_autos
from run_full_intra import now, summary

OUT = INTER / "lote_inter"
MAPS = OUT / "mapas"
DEST = ROOT / "03_wavelets_inter"
PROFILE = "comparavel"


def comparisons(catalog):
    assert catalog["escopo"] == "completo"
    for px, py in combinations(CFG["padroes"], 2):
        xs = sorted(r for r, v in catalog["aquisicoes"].items() if v["padrao"] == px)
        ys = sorted(r for r, v in catalog["aquisicoes"].items() if v["padrao"] == py)
        for x, y in product(xs, ys):
            spec, contexts, autos = comparison_spec(catalog, x, y, PROFILE)
            assert spec["escopo"] == "inter" and spec["comprimento"] == 2048
            assert (spec["pattern_x"], spec["pattern_y"]) == (px, py)
            for run, auto in zip((x, y), autos):
                c = contexts[catalog["aquisicoes"][run]["anchor_id"]]
                native = catalog["aquisicoes"][run]["perfis"]["nativo"]["contexto"]
                assert c["ancora"] == native["ancora"]
                assert c["segmento_fim"]-c["segmento_inicio"] == len(auto["x"]) == 2048
            yield spec, contexts, autos


def filename_base(spec):
    return f"{slug(spec['pattern_x'])}_{slug(spec['pattern_y'])}_{spec['run_x']}_{spec['run_y']}"


def protected(catalog):
    folders = [ROOT/"01_ancoras", ROOT/"02_wavelets", CACHE,
               INTER/"preparacao_comparacoes", INTER/"lote_intra/mapas"]
    files = {p for folder in folders for p in folder.rglob("*") if p.is_file()}
    files |= {SCRIPTS/name for name in ("common.py", "config.json", "select_anchors.py",
              "anchor_diagnostics.py", "wavelet_core.py", "wavelet_cache.py",
              "prepare_comparisons.py", "run_full_intra.py", "plots.py", "plot_config.json")}
    files |= {image_path(run) for run in catalog["aquisicoes"]}
    return {str(p):sha(p) for p in sorted(files)}


def calculate():
    catalog = load_catalog()
    OUT.mkdir(exist_ok=True)
    MAPS.mkdir(exist_ok=True)
    protocol = OUT/"protocolo.json"
    if not protocol.exists():
        write_json(protocol, dict(data_utc=now(), script_sha256=sha(__file__),
            preparacao_sha256=catalog["congelamento_sha256"], perfil=PROFILE,
            comprimento=2048, lambda_px=[50, 1024], escalas=200, pares=213,
            ancora="Centro nativo congelado; recorte comum antes do preparo; sem reamostragem.",
            rotulo="Orientacao visual: WCO baixa favorece distinguir padroes diferentes; sem teste inferencial.",
            preservados_sha256=protected(catalog)))
    old = json.loads(protocol.read_text())
    assert old["script_sha256"] == sha(__file__)
    assert old["preparacao_sha256"] == catalog["congelamento_sha256"]
    rows = []
    for spec, records, autos in comparisons(catalog):
        maps = pair_from_autos(*autos)
        stats = summary(spec, records, maps)
        meta = dict(schema=1, comparacao=spec, metricas=stats,
            cache_individual_x=autos[0]["key"], cache_individual_y=autos[1]["key"],
            preparacao_sha256=catalog["congelamento_sha256"],
            nucleo_sha256=sha(SCRIPTS/"wavelet_core.py"),
            backend_sha256=sha(SCRIPTS/"wavelet_cache.py"), calculo_sha256=sha(__file__),
            armazenamento="Somente WCS/WCO do par; WPS individuais referenciadas.", data_utc=now())
        path = MAPS / f"{pair_key(spec)}.npz"
        if path.exists():
            raise FileExistsError(f"Nao sobrescrever silenciosamente um mapa existente: {path}")
        temp = path.with_suffix(".tmp.npz")
        np.savez_compressed(temp, wxy=maps["wxy"], wco=maps["wco"],
                            denominator_valid=maps["denominator_valid"],
                            meta_json=np.array(json.dumps(meta, ensure_ascii=False)))
        temp.replace(path)
        rows.append(dict(**stats, arquivo=str(path.relative_to(ROOT)), mapa_sha256=sha(path)))
        print(f"Inter {len(rows)}/213: {filename_base(spec)}", flush=True)
    assert len(rows) == catalog["pares_inter_planejados"] == 213
    assert len({frozenset((r["run_x"], r["run_y"])) for r in rows}) == 213
    write_json(OUT/"metricas_wavelets.json", rows)
    write_csv(OUT/"metricas_wavelets.csv", rows)
    return dict(pares=len(rows), grupos=dict(Counter(r["padrao"] for r in rows)))


def read_pair(catalog, row):
    spec, records, (a, b) = comparison_spec(catalog, row["run_x"], row["run_y"], PROFILE)
    path = ROOT/row["arquivo"]
    assert sha(path) == row["mapa_sha256"]
    with np.load(path, allow_pickle=False) as stored:
        meta = json.loads(str(stored["meta_json"]))
        assert meta["comparacao"] == spec
        assert meta["preparacao_sha256"] == catalog["congelamento_sha256"]
        assert (meta["cache_individual_x"], meta["cache_individual_y"]) == (a["key"], b["key"])
        assert meta["nucleo_sha256"] == sha(SCRIPTS/"wavelet_core.py")
        assert meta["backend_sha256"] == sha(SCRIPTS/"wavelet_cache.py")
        assert meta["calculo_sha256"] == sha(__file__)
        assert meta["metricas"] == {k:row[k] for k in meta["metricas"]}
        maps = {key:stored[key].copy() for key in ("wxy", "wco", "denominator_valid")}
    grid = get_grid(spec["comprimento"])
    maps.update(wxx=a["wxx"], wyy=b["wxx"], x=a["x"], y=b["x"], raw_x=a["raw"], raw_y=b["raw"],
        power_mean_x=a["power_mean"], power_mean_y=b["power_mean"],
        **{k:grid[k] for k in ("lags", "scales", "wavelengths", "support")})
    return spec, records, maps


def galleries():
    import plots  # Configure Agg and the same palettes before pyplot.
    import matplotlib.pyplot as plt
    catalog = load_catalog()
    rows = json.loads((OUT/"metricas_wavelets.json").read_text())
    folder = OUT/"galerias_wco"
    folder.mkdir(exist_ok=True)
    index = []
    for px, py in combinations(CFG["padroes"], 2):
        subset = [r for r in rows if (r["pattern_x"], r["pattern_y"]) == (px, py)]
        for page, offset in enumerate(range(0, len(subset), 9), 1):
            group = subset[offset:offset+9]
            fig, axs = plt.subplots(3, 3, figsize=(15, 9), layout="constrained")
            for ax in axs.flat:
                ax.set_visible(False)
            path = folder/f"{slug(px)}_{slug(py)}_{page:02d}.png"
            for position, (ax, row) in enumerate(zip(axs.flat, group), offset+1):
                _, _, maps = read_pair(catalog, row)
                ax.set_visible(True)
                ix = np.unique(np.r_[np.arange(0, len(maps["lags"]), 4), 2047])
                mesh = ax.pcolormesh(maps["lags"][ix], maps["wavelengths"], maps["wco"][:, ix],
                                    shading="auto", cmap="magma", vmin=0, vmax=1)
                ax.set_yscale("log")
                ax.set_yticks([64,256,1024], ["64","256","1024"])
                ax.set_ylim(50,1024)
                for sign in (-1, 1):
                    ax.axvline(sign*32, color="#80909d", lw=.7, ls="--")
                    ax.plot(sign*(2047-np.sqrt(2)*maps["scales"]), maps["wavelengths"],
                            color="#888888", lw=.7, ls=":")
                ax.set_title(f"{position:02d} · {row['run_x']} × {row['run_y']}\n"
                             f"média central {row['wco_media']:.3f}", fontsize=9)
                ax.set_xlabel("Lag (px)"); ax.set_ylabel("λ (px)")
                fig.colorbar(mesh, ax=ax, ticks=[0,.5,1], fraction=.025)
                index.append(dict(pair_id=pair_key(row), grupo=f"{px} × {py}", posicao=position,
                                  galeria=str(path.relative_to(ROOT))))
            fig.suptitle(f"{px} × {py} · WCO · 2.048 px · inspeção {page}", fontsize=13)
            fig.savefig(path, dpi=135)
            plt.close(fig)
    write_json(OUT/"indice_galerias.json", index)
    return dict(pares=len(index))


def render():
    import plots
    catalog = load_catalog()
    rows = json.loads((OUT/"metricas_wavelets.json").read_text())
    labels = json.loads((OUT/"rotulos_visuais.json").read_text())
    assert set(labels) == {pair_key(row) for row in rows}
    DEST.mkdir(exist_ok=True)
    plots.INTER, plots.WAVELETS = OUT, DEST
    # The established layout is unchanged; only the output filename is adapted.
    plots.pair_key = filename_base
    figures = []
    for i, row in enumerate(rows, 1):
        spec, records, maps = read_pair(catalog, row)
        label = labels[pair_key(spec)]
        assert label["rotulo"] in ("favoravel", "desfavoravel")
        assert label["mapa_sha256"] == row["mapa_sha256"]
        plots.pair_figure(spec, records, maps, row, label=label["rotulo"])
        name = filename_base(spec)+"_"+label["rotulo"]+".png"
        figures.append(dict(pair_id=pair_key(spec), arquivo=name, perfil=PROFILE,
                            comprimento=2048, mapa_sha256=row["mapa_sha256"], png_sha256=sha(DEST/name)))
        print(f"PNG inter {i}/{len(rows)}: {name}", flush=True)
    assert {p.name for p in DEST.iterdir()} == {r["arquivo"] for r in figures}
    write_json(OUT/"figuras_wavelets.json", figures)
    return dict(imagens=len(figures))


def verify():
    from audit_method_schuck import independent_prepare, explicit_cwt, temporal_cwt, max_error
    from PIL import Image
    catalog = load_catalog()
    rows = json.loads((OUT/"metricas_wavelets.json").read_text())
    grid = get_grid(2048)
    waves = grid["wavelengths"]
    independently_prepared = {}
    for run in catalog["aquisicoes"]:
        _, a = load_individual(catalog, run, PROFILE)
        x = independent_prepare(a["raw"])
        max_error(x, a["x"])
        rxx = np.correlate(x, x, "full")/len(x)
        p = abs(explicit_cwt(rxx, waves, [2047])[:, 0])
        max_error(p, abs(a["wxx"][:, 2047]))
        independently_prepared[run] = (x, p)
    checks = []
    for i, row in enumerate(rows, 1):
        spec, records, maps = read_pair(catalog, row)
        x, px = independently_prepared[row["run_x"]]
        y, py = independently_prepared[row["run_y"]]
        corr = np.correlate(x, y, "full")/2048
        wxy = temporal_cwt(corr, waves)
        valid = (px>1e-10*px.max()) & (py>1e-10*py.max())
        np.testing.assert_array_equal(valid, maps["denominator_valid"])
        wco = np.full(wxy.shape, np.nan)
        wco[valid] = abs(wxy[valid])**2/(px[valid]*py[valid])[:, None]
        error_wxy = max_error(wxy, maps["wxy"])
        error_wco = max_error(wco, maps["wco"])
        positions = [0,2015,2047,2079,4094]
        error_dot = max_error(explicit_cwt(corr, waves[::25], positions), wxy[::25][:,positions])
        fresh = summary(spec, records, maps)
        assert all(fresh[k] == row[k] for k in fresh)
        assert np.nanmin(wco) >= 0 and np.nanmax(wco) <= 1+1e-10
        checks.append(dict(pair_id=pair_key(row), erro_wcs=error_wxy, erro_wco=error_wco,
                           erro_soma_temporal=error_dot, minimo=float(np.nanmin(wco)),
                           maximo=float(np.nanmax(wco))))
        if i%20 == 0 or i==213:
            print(f"Verificacao independente inter: {i}/213", flush=True)
    old = json.loads((OUT/"protocolo.json").read_text())
    assert old["preservados_sha256"] == protected(catalog), "Arquivos anteriores alterados."
    figures = json.loads((OUT/"figuras_wavelets.json").read_text())
    assert len(figures) == len(rows) == 213
    assert set(p.name for p in DEST.iterdir()) == {r["arquivo"] for r in figures}
    layout = []
    for figure in figures:
        path = DEST/figure["arquivo"]
        assert sha(path) == figure["png_sha256"]
        with Image.open(path) as image:
            assert image.format == "PNG" and image.height > image.width
            image.verify()
        diagnostic = json.loads((OUT/"apresentacao"/(path.stem+".json")).read_text())
        assert "favoravel" not in diagnostic["titulo"].lower()
        assert diagnostic["altura_inserida_cm"]<=24 and diagnostic["fonte_minima_inserida_pt"]>=8
        layout.append(diagnostic)
    write_csv(OUT/"verificacao_pares.csv", checks)
    result = dict(status="PASS", pares=213, imagens=213, preservados=len(old["preservados_sha256"]),
        maior_erro_wcs=max(r["erro_wcs"] for r in checks),
        maior_erro_wco=max(r["erro_wco"] for r in checks),
        menor_wco=min(r["minimo"] for r in checks), maior_wco=max(r["maximo"] for r in checks),
        altura_maxima_cm=max(d["altura_inserida_cm"] for d in layout),
        fonte_minima_pt=min(d["fonte_minima_inserida_pt"] for d in layout),
        fonte="Mesmos panoramas, ancoras, WPS e intra preservados; soma temporal independente para todos os pares.")
    write_json(OUT/"verificacao.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["calculate", "galleries", "render", "verify"])
    args = parser.parse_args()
    print(json.dumps(globals()[args.stage](), ensure_ascii=False, indent=2), flush=True)
