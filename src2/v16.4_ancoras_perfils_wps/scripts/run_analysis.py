"""Compare frozen anchor crops with WPS profiles and matched spatial controls.

Only v16.4 and a temporary Matplotlib configuration directory may be written.
The original v16.1 operations and v16.3 read-only catalogue API are reused.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

sys.dont_write_bytecode = True
os.environ.setdefault("MPLCONFIGDIR", "/tmp/v16_4_ancoras_wps_mpl")
ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT.parent
V161 = SRC / "v16.1_perfis_wps"
V163 = SRC / "v16.3_ancoras"
OUT = ROOT / "results"
CFG = json.loads((ROOT / "config.json").read_text())

sys.path.insert(0, str(V163 / "scripts"))
from common import image_path, load_image  # noqa: E402
from prepare_comparisons import load_catalog, load_individual  # noqa: E402
from wavelet_cache import CACHE, load_entry  # noqa: E402
import wavelet_core as wc  # noqa: E402

sys.path.insert(0, str(V161 / "scripts"))
from profile_core import distribution, quadrature, hellinger, cosine_distance, zero_profile  # noqa: E402
from evaluation import summarize, summarize_pattern  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy  # noqa: E402
from scipy.spatial.distance import euclidean  # noqa: E402
from scipy.stats import rankdata  # noqa: E402

PATTERNS = ["plug", "slug", "churn", "wavy annular"]


def sha(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2,
                                    allow_nan=False) + "\n", encoding="utf-8")


def original_state():
    return {str(p.relative_to(SRC)): (p.stat().st_size, p.stat().st_mtime_ns)
            for root in (V161, V163) for p in root.rglob("*") if p.is_file()}


def independent_prepare(raw):
    y = 255 - np.asarray(raw, float)
    t = np.arange(len(y), dtype=float)
    t -= t.mean()
    y = y - y.mean() - t * (t @ y) / (t @ t)
    y -= y.mean()
    return y / np.sqrt(np.mean(y * y))


def temporal_kernels(n, wavelengths):
    lag = np.arange(-n + 1, n)
    a = wavelengths[:, None] / wc.FOURIER_FACTOR
    u = lag[None, :] / a
    return np.pi**(-.25) / np.sqrt(a) * np.exp(-u*u/2) * np.cos(6*u)


def independent_A(same, cross):
    """Rank-sum formulation, including half credit for exact ties."""
    n, m = len(same), len(cross)
    ranks = rankdata(np.r_[same, cross], method="average")
    return float((ranks[:n].sum() - n*(n+1)/2) / (n*m))


def independent_balanced(frame):
    intra = [v.similaridade.to_numpy() for _, v in frame[frame.mesmo_padrao].groupby("grupo")]
    inter = [v.similaridade.to_numpy() for _, v in frame[~frame.mesmo_padrao].groupby("grupo")]
    return float(np.mean([independent_A(a, b) for a in intra for b in inter]))


def manifest_inputs(catalog):
    paths = {ROOT / "config.json", ROOT / "PROTOCOLO.md",
             V161 / "scripts/profile_core.py", V161 / "scripts/evaluation.py",
             SRC / "v16_spatial2_18_71/scripts/wavelet_core.py",
             V161 / "results/all_pairs.csv", V161 / "results/ranking_summary.csv",
             V163 / "intermediario/preparacao_comparacoes/indice.json",
             V163 / "intermediario/preparacao_comparacoes/congelamento_individual.json",
             V163 / "intermediario/estatistica/descritores.npz",
             V163 / "intermediario/estatistica/verificacao.json",
             V163 / "04_estatistica/dados/indice_ensaios.csv",
             V163 / "04_estatistica/dados/pares_pearson_wco.csv"}
    frozen = json.loads((V163 / "intermediario/preparacao_comparacoes/congelamento_individual.json").read_text())
    paths.update(V163 / p for p in frozen["codigo_sha256"])
    paths.update(SRC / p for p in frozen["fontes_sha256"])
    for batch in ("lote_intra", "lote_inter"):
        paths.add(V163 / "intermediario" / batch / "metricas_wavelets.json")
    for entry in catalog["aquisicoes"].values():
        item = entry["perfis"]["comparavel"]
        for folder in (CACHE / "individuais" / item["cache_key"], CACHE / "grades" / item["grade"]):
            paths.update(p for p in folder.iterdir() if p.is_file())
    current = {str(p.relative_to(SRC)): sha(p) for p in sorted(paths)}
    target = OUT / "entradas_congeladas.json"
    if target.exists():
        previous = json.loads(target.read_text())
        if previous["sha256"] != current:
            raise ValueError("Entradas/protocolo mudaram; não reutilizar silenciosamente este estudo.")
    else:
        write_json(target, dict(data_utc=datetime.now(timezone.utc).isoformat(), sha256=current))
    return current


def run():
    started = time.perf_counter()
    before = original_state()
    OUT.mkdir(parents=True, exist_ok=True)
    assert CFG["janela_px"] == 2048 and CFG["linhas_inclusivas"] == [18, 71]
    assert (CFG["omega0"], CFG["estimador"], CFG["K"]) == (6, "biased", 1)
    assert sha(V163 / "scripts/wavelet_core.py") == sha(SRC / "v16_spatial2_18_71/scripts/wavelet_core.py")
    catalog = load_catalog()  # Checks frozen sources, code, labels and crop bindings; no writes.
    inputs = manifest_inputs(catalog)
    old_checks = json.loads((V163 / "intermediario/estatistica/verificacao.json").read_text())
    archive_path = V163 / "intermediario/estatistica/descritores.npz"
    assert old_checks["status"] == "PASS" and sha(archive_path) == old_checks["descritores_sha256"]
    old = dict(np.load(archive_path, allow_pickle=False))
    index = pd.read_csv(V163 / "04_estatistica/dados/indice_ensaios.csv").set_index("ensaio")
    runs = sorted(catalog["aquisicoes"], key=lambda r: (PATTERNS.index(catalog["aquisicoes"][r]["padrao"]), r))
    assert len(runs) == len(set(runs)) == 25
    assert set(runs) == set(index.index) == set(old["runs"])
    positions = {r: i for i, r in enumerate(old["runs"].tolist())}
    labels = [catalog["aquisicoes"][r]["padrao"] for r in runs]
    assert pd.Series(labels).value_counts().to_dict() == {"churn": 10, "wavy annular": 9, "plug": 3, "slug": 3}
    waves = np.geomspace(CFG["lambda_min_px"], CFG["lambda_max_px"], CFG["numero_escalas"])
    np.testing.assert_array_equal(waves, old["wavelengths"])
    n = CFG["janela_px"]
    assert np.all(n-1 >= np.sqrt(2) * waves/wc.FOURIER_FACTOR)
    conservative = 3*waves/wc.FOURIER_FACTOR <= n/2
    kernels = temporal_kernels(n, waves)
    conditions = CFG["condicoes"]
    powers = {key: [] for key in conditions}
    records, windows, checks = [], [], []
    for i, run_id in enumerate(runs):
        record, auto = load_individual(catalog, run_id, "comparavel")
        pano, raw = load_image(run_id)
        np.testing.assert_array_equal(raw, pano[18:72].sum(axis=0, dtype=np.uint64)/54)
        row = index.loc[run_id]
        assert row.padrao == labels[i] and row.panorama_sha256 == sha(image_path(run_id))
        center = record["ancora"]
        start = record["segmento_inicio"]
        assert (start, record["segmento_fim"], row.ancora_px, row.inicio_px, row.fim_exclusivo_px) == (
            center-n//2, center+n//2, center, center-n//2, center+n//2)
        grid = load_entry(str(CACHE / "grades" / auto["meta"]["spec"]["grade"]))
        np.testing.assert_array_equal(grid["wavelengths"], waves)
        np.testing.assert_array_equal(auto["raw"], raw[start:start+n])
        shift = CFG["deslocamento_sensibilidade_px"]
        cuts = {"ancora_2048": (start, start+n), "inicio_2048": (0, n),
                "centro_2048": ((len(raw)-n)//2, (len(raw)-n)//2+n),
                "fim_2048": (len(raw)-n, len(raw)), "panorama_completo": (0, len(raw)),
                "ancora_menos16": (start-shift, start+n-shift),
                "ancora_mais16": (start+shift, start+n+shift)}
        assert set(cuts) == set(conditions)
        for condition in conditions:
            lo, hi = cuts[condition]
            assert 0 <= lo < hi <= len(raw), (run_id, condition, lo, hi)
            segment = raw[lo:hi]
            x = wc.prepare(segment)
            independent = independent_prepare(segment)
            np.testing.assert_allclose(x, independent, atol=1e-11, rtol=1e-11)
            p_fft = zero_profile(x, waves)
            if condition == "ancora_2048":
                p = abs(auto["wxx"][:, n-1])
                np.testing.assert_allclose(p, old["p"][positions[run_id]], atol=1e-11, rtol=1e-11)
                np.testing.assert_allclose(x, auto["x"], atol=1e-11, rtol=1e-11)
            elif condition == "panorama_completo":
                p = old["p_full"][positions[run_id]]
            else:
                p = p_fft
            np.testing.assert_allclose(p, p_fft, atol=1e-9, rtol=1e-10)
            if len(segment) == n:
                direct_rxx = np.correlate(independent, independent, "full") / n
                p_time = abs(kernels @ direct_rxx)
                np.testing.assert_allclose(p, p_time, atol=1e-9, rtol=1e-10)
                error_time = float(np.max(abs(p-p_time)))
                verification = "correlacao direta + Morlet temporal em todas as escalas"
            else:
                # Existing p_full was calculated by temporal sums and audited in v16.3.
                error_time = float(np.max(abs(p-p_fft)))
                verification = "perfil temporal auditado v16.3 confrontado com FFT v16.1"
            powers[condition].append(p)
            windows.append(dict(ensaio=run_id, padrao=labels[i], condicao=condition,
                inicio_px=lo, fim_exclusivo_px=hi, comprimento_px=hi-lo, ancora_px=center,
                fracao_panorama=(hi-lo)/len(raw)))
            checks.append(dict(ensaio=run_id, condicao=condition,
                erro_preparo=float(np.max(abs(x-independent))),
                erro_wps_fft=float(np.max(abs(p-p_fft))), erro_wps_temporal=error_time,
                verificacao=verification))
        records.append(dict(numero=i+1, ensaio=run_id, padrao=labels[i], ancora_px=center,
            inicio_px=start, fim_exclusivo_px=start+n, comprimento_panorama_px=len(raw),
            tipo_ancora=record["tipo"], ancora_id=record["anchor_id"],
            cache_individual=auto["key"], panorama_sha256=sha(image_path(run_id)),
            costuras_contexto=record["costuras_contexto"],
            outros_candidatos=record["outros_candidatos_estaveis_no_contexto"],
            alertas_opticos=record.get("alertas_opticos", "")))
        print(f"Perfis e fontes conferidos: {i+1}/25 — {run_id}", flush=True)
    powers = {key: np.asarray(value) for key, value in powers.items()}
    pd.DataFrame(records).to_csv(OUT / "indice_ensaios.csv", index=False)
    pd.DataFrame(windows).to_csv(OUT / "recortes.csv", index=False)
    pd.DataFrame(checks).to_csv(OUT / "verificacao_perfis.csv", index=False)
    lengths = {(r["ensaio"], r["condicao"]): r["comprimento_px"] for r in windows}
    rows, normalized, methods, masses = [], [], {}, {}

    def pair_row(method, i, j, value, higher, nx, ny):
        return dict(metodo=method, aquisicao_x=runs[i], aquisicao_y=runs[j],
            padrao_x=labels[i], padrao_y=labels[j], mesmo_padrao=labels[i] == labels[j],
            grupo=" / ".join(sorted([labels[i], labels[j]])), valor=float(value),
            similaridade=float(value if higher else -value), n_x=nx, n_y=ny)

    max_distance_error = 0.
    for condition in conditions:
        for domain, mask in (("principal", np.ones(len(waves), dtype=bool)), ("conservador", conservative)):
            axis = waves[mask]
            values = powers[condition][:, mask]
            mass = np.asarray([distribution(p, axis) for p in values])
            np.testing.assert_allclose(mass.sum(axis=1), 1., atol=1e-14)
            masses[condition + "__" + domain] = mass
            method = condition + "__hellinger_" + domain
            methods[method] = dict(condicao=condition, dominio=domain, metrica="Hellinger",
                lambda_min_px=float(axis.min()), lambda_max_px=float(axis.max()), escalas=len(axis),
                maior_significa_mais_semelhante=False)
            for i, p in enumerate(values):
                weights = quadrature(axis)
                area = float(weights @ p)
                np.testing.assert_allclose(distribution(7.3*p, axis), mass[i], atol=1e-14)
                for k, wave in enumerate(axis):
                    normalized.append(dict(ensaio=runs[i], padrao=labels[i], condicao=condition,
                        dominio=domain, lambda_px=wave, wps_lag_zero=p[k], peso_log=weights[k],
                        area_log=area, densidade_log=p[k]/area, massa=mass[i, k]))
            for i, j in combinations(range(25), 2):
                h = hellinger(mass[i], mass[j])
                alt = euclidean(np.sqrt(mass[i]), np.sqrt(mass[j])) / np.sqrt(2)
                max_distance_error = max(max_distance_error, abs(h-alt))
                assert 0 <= h <= 1
                rows.append(pair_row(method, i, j, h, False, lengths[runs[i], condition], lengths[runs[j], condition]))
        if condition in ("ancora_2048", "panorama_completo"):
            method = condition + "__cosseno_principal"
            methods[method] = dict(condicao=condition, dominio="principal", metrica="Cosseno",
                lambda_min_px=50, lambda_max_px=1024, escalas=200, maior_significa_mais_semelhante=False)
            for i, j in combinations(range(25), 2):
                value = cosine_distance(powers[condition][i], powers[condition][j], waves)
                rows.append(pair_row(method, i, j, value, False, lengths[runs[i], condition], lengths[runs[j], condition]))
    assert hellinger(np.array([1., 0]), np.array([0., 1])) == 1
    assert hellinger(np.array([.2, .8]), np.array([.2, .8])) == 0
    assert independent_A([0, 0], [0, 0]) == .5
    assert independent_A([1, 2], [-1, 0]) == 1
    assert independent_A([-1, 0], [1, 2]) == 0
    pd.DataFrame(normalized).to_csv(OUT / "perfis_normalizados.csv", index=False)
    np.savez_compressed(OUT / "perfis.npz", runs=np.asarray(runs), patterns=np.asarray(labels),
        wavelengths=waves, conditions=np.asarray(conditions), powers=np.asarray([powers[c] for c in conditions]),
        conservative_mask=conservative, **masses)

    references = pd.read_csv(V163 / "04_estatistica/dados/pares_pearson_wco.csv")
    ref_lookup = {frozenset((r.ensaio_x, r.ensaio_y)): r for r in references.itertuples()}
    wco_rows = []
    for batch in ("lote_intra", "lote_inter"):
        wco_rows += json.loads((V163 / "intermediario" / batch / "metricas_wavelets.json").read_text())
    wco_lookup = {frozenset((r["run_x"], r["run_y"])): r for r in wco_rows}
    assert len(ref_lookup) == len(references) == len(wco_lookup) == len(wco_rows) == 300
    old_pairs = pd.read_csv(V161 / "results/all_pairs.csv")
    old_pairs = old_pairs[(old_pairs.metodo == "wps_principal") & old_pairs.aquisicao_x.isin(runs) & old_pairs.aquisicao_y.isin(runs)]
    history = {frozenset((r.aquisicao_x, r.aquisicao_y)): r for r in old_pairs.itertuples()}
    assert len(history) == len(old_pairs) == 300
    for method, metric, low, high, count in (
        ("v163_pearson_ancora", "Pearson residual", 50, 1024, 200),
        ("v163_wco_central", "WCO media central", 50, 1024, 200),
        ("v161_principal_25", "Hellinger histórico — faixa diferente", 50, float(np.geomspace(50, 3000, 200)[186]), 187)):
        methods[method] = dict(condicao="ancora_2048" if method.startswith("v163") else "panorama_completo",
            dominio="principal" if method.startswith("v163") else "historico", metrica=metric,
            lambda_min_px=low, lambda_max_px=high, escalas=count,
            maior_significa_mais_semelhante=method.startswith("v163"))
    max_pearson_error, max_wco_error = 0., 0.
    for i, j in combinations(range(25), 2):
        key = frozenset((runs[i], runs[j]))
        ref, wco = ref_lookup[key], wco_lookup[key]
        assert {ref.padrao_x, ref.padrao_y} == {labels[i], labels[j]}
        assert ref.tipo == ("intra" if labels[i] == labels[j] else "inter")
        assert wco["perfil"] == "comparavel" and wco["comprimento"] == 2048
        assert wco["escalas_validas"] == 200 and wco["lags_resumo"] == 65
        r = np.corrcoef(old["z"][positions[runs[i]]], old["z"][positions[runs[j]]])[0, 1]
        max_pearson_error = max(max_pearson_error, abs(r-ref.pearson_r))
        max_wco_error = max(max_wco_error, abs(wco["wco_media"]-ref.wco_media_central))
        rows.append(pair_row("v163_pearson_ancora", i, j, ref.pearson_r, True, n, n))
        rows.append(pair_row("v163_wco_central", i, j, ref.wco_media_central, True, n, n))
        previous = history[key]
        assert bool(previous.mesmo_padrao) == (labels[i] == labels[j])
        rows.append(pair_row("v161_principal_25", i, j, previous.valor, False,
                             lengths[runs[i], "panorama_completo"], lengths[runs[j], "panorama_completo"]))
    assert max_pearson_error < 1e-12 and max_wco_error < 1e-12
    frame = pd.DataFrame(rows)
    assert not frame.duplicated(["metodo", "aquisicao_x", "aquisicao_y"]).any()
    summary, by_pattern, deletion, groups = [], [], [], []
    max_A_error = 0.
    for method, part in frame.groupby("metodo", sort=False):
        result = summarize(part)
        assert result["pares_intra"] == 87 and result["pares_entre"] == 213
        assert result["grupos_intra"] == 4 and result["grupos_entre"] == 6
        same, cross = part[part.mesmo_padrao], part[~part.mesmo_padrao]
        a_error = abs(result["A_pares"] - independent_A(same.similaridade, cross.similaridade))
        balanced_error = abs(result["A_balanceada"] - independent_balanced(part))
        max_A_error = max(max_A_error, a_error, balanced_error)
        comparisons = same.similaridade.to_numpy()[:, None] - cross.similaridade.to_numpy()[None, :]
        summary.append(dict(metodo=method, **methods[method], **result,
            mediana_intra=float(same.valor.median()), mediana_entre=float(cross.valor.median()),
            ordens_favoraveis=int((comparisons > 0).sum()), empates=int((comparisons == 0).sum()),
            ordens_contrarias=int((comparisons < 0).sum()), comparacoes=87*213))
        for pattern in PATTERNS:
            by_pattern.append(dict(metodo=method, padrao=pattern, **summarize_pattern(part, pattern)))
        for run_id in runs:
            subset = part[(part.aquisicao_x != run_id) & (part.aquisicao_y != run_id)]
            deletion.append(dict(metodo=method, aquisicao_retirada=run_id,
                padrao_retirado=labels[runs.index(run_id)], **summarize(subset)))
        for group, group_frame in part.groupby("grupo"):
            groups.append(dict(metodo=method, grupo=group, pares=len(group_frame),
                mesmo_padrao=bool(group_frame.mesmo_padrao.iloc[0]),
                media=float(group_frame.valor.mean()), mediana=float(group_frame.valor.median()),
                q25=float(group_frame.valor.quantile(.25)), q75=float(group_frame.valor.quantile(.75))))
    assert max_A_error < 1e-12 and max_distance_error < 1e-12
    sums, deletions = pd.DataFrame(summary).set_index("metodo"), pd.DataFrame(deletion)
    main = "ancora_2048__hellinger_principal"
    effects = []
    for other in sums.index:
        if other == main:
            continue
        paired = deletions[deletions.metodo == main].set_index("aquisicao_retirada")[
            ["A_pares", "A_balanceada"]] - deletions[deletions.metodo == other].set_index("aquisicao_retirada")[["A_pares", "A_balanceada"]]
        effects.append(dict(referencia=other,
            delta_A_pares=float(sums.loc[main, "A_pares"] - sums.loc[other, "A_pares"]),
            delta_A_balanceada=float(sums.loc[main, "A_balanceada"] - sums.loc[other, "A_balanceada"]),
            delta_retiradas_min=float(paired.A_pares.min()), delta_retiradas_max=float(paired.A_pares.max()),
            delta_retiradas_balanceada_min=float(paired.A_balanceada.min()),
            delta_retiradas_balanceada_max=float(paired.A_balanceada.max()),
            retiradas_ancora_maior=int((paired.A_pares > 0).sum()),
            retiradas_ancora_maior_balanceada=int((paired.A_balanceada > 0).sum())))
    sensitivity = []
    for condition in conditions:
        for i, run_id in enumerate(runs):
            sensitivity.append(dict(ensaio=run_id, padrao=labels[i], controle=condition,
                H_em_relacao_ancora=hellinger(masses["ancora_2048__principal"][i], masses[condition+"__principal"][i])))
    for name, table in (("pares", frame), ("concordancia", sums.reset_index()),
                        ("por_padrao", pd.DataFrame(by_pattern)), ("influencia_aquisicoes", deletions),
                        ("grupos_pares", pd.DataFrame(groups)), ("efeitos_pareados", pd.DataFrame(effects)),
                        ("sensibilidade_perfis", pd.DataFrame(sensitivity))):
        table.to_csv(OUT / f"{name}.csv", index=False)
    for relative, expected in inputs.items():
        assert sha(SRC / relative) == expected, f"Entrada mudou durante a execução: {relative}"
    assert before == original_state(), "Arquivos anteriores mudaram durante a análise."
    write_json(OUT / "verificacao.json", dict(status="PASS", ensaios=25, condicoes=len(conditions),
        metodos=len(methods), pares_por_metodo=300, perfis_verificados=len(checks),
        escalas_principais=200, escalas_conservadoras=int(conservative.sum()),
        lambda_conservadora_max_px=float(waves[conservative].max()),
        maior_erro_preparo=max(v["erro_preparo"] for v in checks),
        maior_erro_wps=max(max(v["erro_wps_temporal"], v["erro_wps_fft"]) for v in checks),
        maior_erro_H=max_distance_error, maior_erro_A=max_A_error,
        maior_erro_referencia_Pearson=max_pearson_error, maior_erro_referencia_WCO=max_wco_error,
        arquivos_fontes_congelados=len(inputs), arquivos_anteriores_preservados=len(before),
        status_arquivos_anteriores="mesmos nomes, tamanhos e mtimes antes/depois; hashes das entradas consumidas iguais",
        verificacoes=["fontes, catálogo e recortes", "preparo por formulação independente",
            "150 recortes: autocorrelação direta e Morlet temporal em todas as escalas",
            "25 panoramas completos: perfil temporal auditado confrontado com FFT",
            "Hellinger contra distância euclidiana e casos limite", "invariância de normalização à amplitude",
            "A por soma de postos e macro média entre grupos", "referências vinculadas por IDs, não por posição"],
        limites="Não reaudita quadros originais, escolha morfológica ou os 300 mapas WCO; reutiliza a auditoria v16.3.",
        duracao_segundos=time.perf_counter()-started, data_utc=datetime.now(timezone.utc).isoformat(),
        versoes=dict(python=sys.version.split()[0], numpy=np.__version__, scipy=scipy.__version__, pandas=pd.__version__),
        codigo_sha256={"scripts/run_analysis.py": sha(__file__)}))
    print(sums[["A_pares", "A_balanceada"]].to_string(), flush=True)
    print(f"PASS — {time.perf_counter()-started:.1f} s", flush=True)


if __name__ == "__main__":
    run()
