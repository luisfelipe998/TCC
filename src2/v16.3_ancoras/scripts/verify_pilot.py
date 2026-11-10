"""Check selection, provenance, spatial domains, equations and final artifacts."""
import json
from datetime import datetime
import numpy as np
from PIL import Image

from common import (ANCHORS, CFG, CONCLUSIONS, INTER, OLD, ROOT, SCRIPTS, WAVELETS,
                    context_ok, load_image, read_csv, selected_manifest, sha, slug, write_json)
from compute_wavelets import pair_key
from wavelet_core import correlation


def check_image(stem):
    with Image.open(stem.with_suffix(".png")) as im:
        width, height = im.size
        assert .65 <= width / height <= 1.0, "PNG deve caber em retrato sem reducao excessiva."
        assert width >= 1800 and height >= 1900
    assert not stem.with_suffix(".pdf").exists()
    return dict(arquivo=str(stem.relative_to(ROOT)), pixels=[width, height],
                altura_ao_inserir_com_16cm=16*height/width)


def verify():
    selection = selected_manifest()
    anchors = selection["ancoras"]
    assert len(selection["pares"]) == 4 and len(anchors) == 8
    assert {r["padrao"] for r in selection["pares"]} == set(CFG["padroes"])
    assert len({r["aquisicao"] for r in anchors}) == 8
    assert not list(CONCLUSIONS.iterdir()), "A pasta de conclusoes deve continuar vazia."
    assert sha(SCRIPTS / "wavelet_core.py") == sha(OLD / "scripts/wavelet_core.py")
    catalog = read_csv(INTER / "catalogo_candidatos.csv")
    reviews = {r["event_id"]: r for r in read_csv(INTER / "revisao_visual_herdada.csv")}
    visual = json.loads((INTER / "inspecao_visual_piloto.json").read_text())
    assert visual["antes_das_wavelets"] and len(visual["ancoras"]) == 8
    visual_ids = {r["anchor_id"] for r in visual["ancoras"] if r["verificada_visualmente"]}
    assert visual_ids == {r["anchor_id"] for r in anchors}
    anchor_by_id = {r["anchor_id"]: r for r in anchors}
    expected_anchors, expected_pairs, figures = set(), set(), []
    for record in anchors:
        run, center, n = record["aquisicao"], record["ancora"], record["comprimento"]
        assert not run.endswith("_002") and record["padrao"] != "dispersed bubbles"
        _, raw = load_image(run)
        assert context_ok(center, len(raw), n)
        assert record["segmento_inicio"] == center - n//2
        assert record["segmento_fim"] == center + n//2
        assert n >= record["comprimento_minimo"] and n % CFG["contexto"]["passo"] == 0
        margin = CFG["contexto"]["margem_evento" if record["tipo"] == "evento_delimitado" else "margem_vale"]
        assert record["inicio_suporte"] - record["segmento_inicio"] >= margin
        assert record["segmento_fim"] - record["fim_suporte"] >= margin
        source = ROOT.parent / record["arquivo_quadro_origem"]
        assert sha(source) == record["quadro_origem_sha256"] and record["coluna_origem_conferida"]
        pano, _ = load_image(run)
        with Image.open(source) as frame:
            original = np.asarray(frame.convert("L"))
        np.testing.assert_array_equal(pano[:, center], original[:, record["coluna_quadro_ancora"]])
        eligible = [r for r in catalog if r["aquisicao"] == run and r["elegivel"] == "True"]
        if record["tipo"] == "evento_delimitado":
            assert reviews[record["anchor_id"]]["aceito"] == "True"
            assert center == (record["inicio"] + record["fim"]) // 2
            margin = CFG["contexto"]["margem_evento"]
            assert record["inicio"] - record["segmento_inicio"] >= margin
            assert record["segmento_fim"] - record["fim"] >= margin
            shifts = [abs(record[k] - record[f"{k}_limiar{q}"])
                      for k in ("inicio", "fim") for q in (55, 65)]
            assert max(shifts) == record["variacao_maxima"]
            w = CFG["eventos"]["nitidez_faixa_px"]
            a, b = record["inicio"], record["fim"]
            contrast = min(np.median(raw[a-w:a])-np.median(raw[a:a+w]),
                           np.median(raw[b:b+w])-np.median(raw[b-w:b]))
            assert abs(contrast-record["contraste_minimo"]) < 1e-12
            winner = min(eligible, key=lambda r: (float(r["variacao_maxima"]),
                         -float(r["contraste_minimo"]), int(r["ancora"])))
        else:
            shifts = [abs(record[f"ancora_sigma{s}"] - center) for s in (12, 20)]
            assert max(shifts) == record["variacao_maxima"] <= CFG["textura"]["tolerancia"]
            assert 16 <= record["largura"] <= 512
            assert record["prominencia"] >= record["limiar_prominencia"]
            margin = CFG["contexto"]["margem_vale"]
            assert record["base_esquerda"] - record["segmento_inicio"] >= margin
            assert record["segmento_fim"] - record["base_direita"] >= margin
            winner = min(eligible, key=lambda r: (-float(r["prominencia"]), float(r["variacao_maxima"]), int(r["ancora"])))
        assert winner["anchor_id"] == record["anchor_id"]
        stem = ANCHORS / f"{slug(record['padrao'])}_{record['anchor_id']}"
        expected_anchors.add(stem.with_suffix(".png").name)
        figures.append(check_image(stem))

    for trial in selection["ensaios"]:
        rows = [r for r in catalog if r["aquisicao"] == trial["aquisicao"]]
        assert trial["elegiveis"] == sum(r["elegivel"] == "True" for r in rows)
    diagnostics = read_csv(INTER / "diagnosticos_ancoras.csv")
    assert {r["anchor_id"] for r in diagnostics} == {r["anchor_id"] for r in anchors}

    frozen = json.loads((INTER / "congelamento.json").read_text())
    labels = json.loads((INTER / "rotulos_visuais.json").read_text())
    summaries = {r["padrao"]: r for r in json.loads((INTER / "metricas_wavelets.json").read_text())}
    numeric = []
    for comparison in selection["pares"]:
        n = comparison["comprimento"]
        left, right = [anchor_by_id[comparison[k]] for k in ("anchor_x", "anchor_y")]
        assert left["padrao"] == right["padrao"] == comparison["padrao"]
        assert left["comprimento"] == right["comprimento"] == n
        assert left["aquisicao"] != right["aquisicao"]
        with np.load(INTER / "mapas" / f"{pair_key(comparison)}.npz") as stored:
            maps = dict(stored)
        meta = json.loads(str(maps["meta_json"]))
        assert meta["congelamento_sha256"] == sha(INTER / "congelamento.json")
        assert datetime.fromisoformat(meta["data_utc"]) > datetime.fromisoformat(frozen["data_utc"])
        assert datetime.fromisoformat(meta["data_utc"]) > datetime.fromisoformat(visual["data_utc"])
        for key, record in (("x", left), ("y", right)):
            _, raw = load_image(record["aquisicao"])
            segment = raw[record["segmento_inicio"]:record["segmento_fim"]]
            np.testing.assert_array_equal(segment, maps["raw_"+key])
            # Independent linear least-squares preparation, outside the copied core.
            signal = 255 - segment
            axis = np.arange(n, dtype=float)
            prepared = signal - np.polyval(np.polyfit(axis, signal, 1), axis)
            prepared -= prepared.mean()
            prepared /= prepared.std()
            np.testing.assert_allclose(prepared, maps[key], atol=1e-12, rtol=1e-12)
        direct = np.correlate(maps["x"], maps["y"], mode="full") / n
        np.testing.assert_allclose(correlation(maps["x"], maps["y"]), direct, atol=1e-12, rtol=1e-12)
        assert maps["wco"].shape == (200, 2*n-1)
        assert maps["lags"][n-1] == 0
        assert maps["wavelengths"][0] == 50 and maps["wavelengths"][-1] == n/2
        valid = maps["denominator_valid"]
        denominator = np.abs(maps["wxx"][:, n-1]) * np.abs(maps["wyy"][:, n-1])
        independent = np.abs(maps["wxy"][valid])**2 / denominator[valid, None]
        error = float(np.max(abs(independent - maps["wco"][valid])))
        assert error < 1e-12
        support = n - 1 - abs(maps["lags"])[None, :] >= np.sqrt(2)*maps["scales"][:, None]
        np.testing.assert_array_equal(support, maps["support"])
        for key, mean in (("wxx", "power_mean_x"), ("wyy", "power_mean_y")):
            check = np.array([np.abs(row)[mask].mean() for row, mask in zip(maps[key], support)])
            np.testing.assert_allclose(check, maps[mean], atol=1e-12, rtol=1e-12)
        domain = abs(maps["lags"]) <= CFG["wavelet"]["lag_resumo"]
        assert support[:, domain].all() and domain.sum() == 65
        mean = float(maps["wco"][valid][:, domain].mean(axis=1).mean())
        assert abs(mean - summaries[comparison["padrao"]]["wco_media"]) < 1e-12
        label = labels[comparison["padrao"]]["rotulo"]
        stem = WAVELETS / f"{pair_key(comparison)}_{label}"
        expected_pairs.add(stem.with_suffix(".png").name)
        figures.append(check_image(stem))
        numeric.append(dict(padrao=comparison["padrao"], comprimento=n, erro_equacao=error,
                            escalas_validas=int(valid.sum()), wco_media=mean))

    # Boundary cases directly relevant to the requested absence handling.
    from select_anchors import event_candidates, match_minimum, texture_candidates
    flat = np.full(4096, 128.)
    rows, chosen = texture_candidates("verificacao", "churn", flat)
    assert not rows and chosen is None
    assert match_minimum(100, [dict(ancora=99, contexto_valido=True),
                                dict(ancora=101, contexto_valido=True)]) is None
    pano = np.full((80, 4096), 200, dtype=np.uint8)
    pano[18:72, 1500:1800] = 20
    _, chosen = event_candidates("verificacao", "plug", pano, pano[18:72].mean(axis=0), {}, {})
    assert chosen is None, "Sem revisao visual, o script nao pode aprovar a morfologia."
    x, y = np.zeros(128), np.zeros(128)
    x[30], y[47] = 1, 1
    assert np.arange(-127, 128)[np.argmax(correlation(x, y))] == -17

    # A larger requested context must back off rather than eliminate the only slug anchor.
    from select_anchors import common_window
    examples = [(1112, 1536, 5071), (2264, 1024, 6323)]
    groups = []
    for center, required, panorama in examples:
        row = dict(ancora=center, comprimento_minimo=required, tamanho_panorama=panorama,
                   elegivel=True, variacao_maxima=1, contraste_minimo=0.)
        groups.append(dict(pattern="slug", rows=[row], best=row))
    original_window = CFG["janelas_por_padrao"]["slug"]
    try:
        CFG["janelas_por_padrao"]["slug"] = 3072
        assert common_window(groups) == 2048
        assert all(g["best"]["selecionada"] for g in groups)
    finally:
        CFG["janelas_por_padrao"]["slug"] = original_window

    assert {p.name for p in ANCHORS.iterdir()} == expected_anchors
    assert {p.name for p in WAVELETS.iterdir()} == expected_pairs
    assert len(list((INTER / "mapas").glob("*.npz"))) == 4
    assert not list(ROOT.rglob("*.pdf")), "Esta analise deve conter somente imagens PNG."
    assert CFG["figuras"]["paleta_wps_wcs"] == "viridis"
    assert CFG["figuras"]["paleta_wco"] == "magma"
    artifacts = {str(p.relative_to(ROOT)): sha(p) for folder in (ANCHORS, WAVELETS, SCRIPTS)
                 for p in sorted(folder.iterdir()) if p.is_file()}
    report = dict(status="PASS", ancoras=8, pares=4, imagens=len(figures),
                  conclusoes_vazias=True, numero_de_candidatos=len(catalog),
                  verificacoes_numericas=numeric, imagens_png=figures, artefatos_sha256=artifacts,
                  casos_limite=["sinal_constante_sem_ancora", "correspondencia_ambigua_rejeitada",
                                "morfologia_sem_revisao_rejeitada", "convencao_de_lag",
                                "janela_maior_inviavel_reduzida_para_2048"])
    write_json(INTER / "verificacao.json", report)
    print(f"PASS: {len(anchors)} ancoras, 4 pares, {len(figures)} PNGs compactos; sem PDFs; conclusoes vazias.")
