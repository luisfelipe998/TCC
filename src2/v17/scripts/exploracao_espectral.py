"""Exploração descritiva independente; nenhuma predição ou aprendizagem.

Execute da raiz: python3 src2/v17/scripts/exploracao_espectral.py
Parâmetros fixados em ../PROTOCOLO.txt. Somente escreve dentro de v17.
"""
from pathlib import Path
import os
import sys
import hashlib
import json
from itertools import combinations

sys.dont_write_bytecode = True
BASE = Path(__file__).resolve().parents[1]
PROJECT = BASE.parents[1]
OUT = BASE / "resultados/exploracao"
OUT.mkdir(parents=True, exist_ok=True)
os.environ["MPLCONFIGDIR"] = str(BASE / "resultados/cache_matplotlib")
for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"]:
    os.environ[name] = "1"
import numpy as np
import pandas as pd
import scipy
from scipy.signal import correlate, detrend, welch
from scipy.spatial.distance import pdist, squareform
from scipy.stats import wasserstein_distance
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SOURCE = PROJECT / "src2"
FACTOR = 4 * np.pi / (6 + np.sqrt(38))
CLASSES = ["churn", "wavy annular", "plug", "slug"]
SCOPES = {"quatro_padroes": CLASSES, "churn_annular": CLASSES[:2], "plug_slug": CLASSES[2:]}
INPUTS = {}


def track(path):
    path = Path(path)
    INPUTS[str(path.relative_to(PROJECT))] = hashlib.file_digest(path.open("rb"), "sha256").hexdigest()
    return path


def save_json(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def prepare(raw):
    x = detrend(255 - np.asarray(raw, float), type="linear")
    x -= x.mean()
    assert x.std() > 1e-12
    return x / x.std()


def kernel(lam, n):
    lag = np.arange(-n + 1, n)
    a = lam / FACTOR
    tau = lag[None, :] / a[:, None]
    return np.pi ** (-0.25) / np.sqrt(a[:, None]) * np.exp(-0.5 * tau**2 - 6j * tau)


def power_zero(x, k):
    # Soma discreta direta, não chama o código wavelet das fontes.
    r = correlate(x, x, mode="full", method="fft") / len(x)
    return np.abs(k @ r)


def mass(profile, lam):
    t = np.log(lam)
    w = np.r_[np.diff(t)[0] / 2, (t[2:] - t[:-2]) / 2, np.diff(t)[-1] / 2]
    p = profile * w
    return p / p.sum(axis=-1, keepdims=True)


def residual(profile, lam):
    design = np.column_stack([np.ones(len(lam)), np.log(lam)])
    y = np.log(profile + 1e-12)
    z = y - (np.linalg.lstsq(design, y.T, rcond=None)[0].T @ design.T)
    z -= z.mean(axis=-1, keepdims=True)
    norms = np.linalg.norm(z, axis=-1, keepdims=True)
    tolerance = 1e-12 * np.maximum(1.0, np.linalg.norm(y, axis=-1, keepdims=True))
    if np.any(norms <= tolerance):
        raise ValueError("Perfil sem forma residual identificável; Pearson é indefinido.")
    return z / norms


def width(p, lam):
    # Convenção contínua explícita de interpolação da CDF discretizada.
    cdf = p.cumsum()
    t = np.log(lam)
    q = np.interp([0.25, 0.75], np.r_[0, cdf], np.r_[t[0], t])
    return q[1] - q[0]


def distances(profiles, lam):
    p = mass(profiles, lam)
    cdf = p.cumsum(axis=1)
    h = squareform(pdist(np.sqrt(p))) / np.sqrt(2)
    w = np.sum(np.abs(cdf[:, None, :-1] - cdf[None, :, :-1]) * np.diff(np.log(lam)), axis=-1)
    z = residual(profiles, lam)
    r = squareform(pdist(z))  # chordal distance sqrt(2*(1-r)).
    widths = np.array([width(v, lam) for v in p])
    # Independent identity checks against SciPy, not tests mirroring implementation.
    assert abs(w[0, 1] - wasserstein_distance(np.log(lam), np.log(lam), p[0], p[1])) < 1e-12
    for d in [h, w, r]:
        assert np.isfinite(d).all() and np.allclose(d, d.T) and np.allclose(np.diag(d), 0)
    return {"hellinger": h, "wasserstein": w, "pearson_residual": r,
            "largura": squareform(pdist(widths[:, None]))}, widths


def separation(d, labels):
    classes = np.unique(labels)
    sil = np.zeros(len(labels))
    within = []
    between = []
    for c in classes:
        ix = np.where(labels == c)[0]
        assert len(ix) > 1, "Singleton não entra na análise de silhueta."
        a = d[np.ix_(ix, ix)].sum(axis=1) / (len(ix) - 1)
        b = np.min([d[np.ix_(ix, np.where(labels == other)[0])].mean(axis=1)
                    for other in classes if other != c], axis=0)
        sil[ix] = np.divide(b - a, np.maximum(a, b), out=np.zeros(len(ix)), where=np.maximum(a, b) > 0)
        within.append(d[np.ix_(ix, ix)][np.triu_indices(len(ix), 1)].mean())
    for a, b in combinations(classes, 2):
        between.append(d[np.ix_(labels == a, labels == b)].mean())
    return {"silhueta_macro": float(np.mean([sil[labels == c].mean() for c in classes])),
            "distancia_intra_macro": float(np.mean(within)),
            "distancia_inter_macro": float(np.mean(between)),
            "contraste_distancia_macro": float(np.mean(between) - np.mean(within)),
            "n_silhueta_negativa": int((sil < 0).sum())}, sil


def composition():
    p = SOURCE / "v16_spatial2_18_71/04_estatistica/dados"
    pairs = pd.read_csv(track(p / "pares_pearson_wco.csv"))
    output, by_group = [], []
    for name, exclusions in {"26_todos": [], "25_sem_bolha": ["m07q01_001"],
                             "25_sem_m03q40": ["m03q40_001"],
                             "24_sem_ambos": ["m07q01_001", "m03q40_001"]}.items():
        df = pairs[~pairs.ensaio_x.isin(exclusions) & ~pairs.ensaio_y.isin(exclusions)].copy()
        canonical = [tuple(sorted((a, b))) for a, b in zip(df.padrao_x, df.padrao_y)]
        df["grupo_a"] = [a for a, b in canonical]
        df["grupo_b"] = [b for a, b in canonical]
        for metric in ["pearson_r", "escore_piso05"]:
            m = {}
            for kind in ["intra", "inter"]:
                v = df[df.tipo == kind]
                groups = v.groupby(["grupo_a", "grupo_b"])[metric].mean()
                m[kind] = float(v[metric].mean())
                m[kind + "_macro"] = float(groups.mean())
                output.append(dict(conjunto=name, metrica=metric, tipo=kind, pares=len(v),
                                   grupos=len(groups), media=m[kind], media_balanceada=m[kind + "_macro"]))
                for (a, b), val in groups.items():
                    by_group.append(dict(conjunto=name, metrica=metric, padrao_x=a, padrao_y=b,
                                         tipo=kind, media=val))
            output.append(dict(conjunto=name, metrica=metric, tipo="intra_menos_inter", pares=len(df),
                               grupos=None, media=m["intra"] - m["inter"],
                               media_balanceada=m["intra_macro"] - m["inter_macro"]))
    pd.DataFrame(output).to_csv(OUT / "efeito_composicao_tcc.csv", index=False)
    pd.DataFrame(by_group).to_csv(OUT / "efeito_composicao_grupos.csv", index=False)


def main():
    index = pd.read_csv(track(SOURCE / "estatistica8/dados/inventario.csv"))
    runs = index.ensaio.to_numpy(str)
    labels = index.padrao.to_numpy(str)
    assert len(runs) == 25 and index.groupby("padrao").size().to_dict() == {
        "churn": 10, "plug": 3, "slug": 3, "wavy annular": 9}
    powers = {}
    features, windows, checks = [], [], []
    grids = {cap: np.geomspace(50, cap, 200) for cap in [512, 1024]}
    kernels = {cap: kernel(lam, 2048) for cap, lam in grids.items()}
    frozen = np.load(track(SOURCE / "estatistica8/dados/perfis_ancora.npz"))
    for i, row in enumerate(index.itertuples()):
        path = track(SOURCE / "data_processed" / row.ensaio / "__panorama_v11.png")
        with Image.open(path) as im:
            raw = np.asarray(im.convert("L"))[18:72].mean(axis=0)
        assert len(raw) == row.panorama_px
        starts = {"prefixo": 0, "centro": (len(raw) - 2048) // 2,
                  "sufixo": len(raw) - 2048, "ancora": row.inicio_px}
        for window, start in starts.items():
            x = prepare(raw[start:start + 2048])
            windows.append(dict(ensaio=row.ensaio, padrao=row.padrao, janela=window,
                                inicio=start, fim_exclusivo=start + 2048, panorama_px=len(raw)))
            f, psd = welch(x, fs=1.0, window="hann", nperseg=1024, noverlap=512,
                           nfft=8192, detrend=False, scaling="density")
            for cap, lam in grids.items():
                alt = power_zero(x, kernels[cap])
                if cap == 1024 and window == "ancora":
                    err = float(np.max(np.abs(alt - frozen["p"][i])))
                    checks.append(dict(ensaio=row.ensaio, erro_wps_ancora=err))
                    assert err < 1e-10
                profiles = {"wps": alt, "wps_ponderada": alt / np.sqrt(lam),
                            "welch_logenergia": np.interp(1 / lam, f, psd) / lam}
                for method, profile in profiles.items():
                    assert np.isfinite(profile).all() and np.all(profile > 0)
                    key = (cap, window, method)
                    powers.setdefault(key, []).append(profile)
                    p = mass(profile, lam)
                    features.append(dict(ensaio=row.ensaio, padrao=row.padrao, m=int(row.ensaio[1:3]),
                                         q=int(row.ensaio[4:6]), teto_lambda=cap, janela=window, metodo=method,
                                         largura=width(p, lam), lambda_pico=lam[profile.argmax()],
                                         massa_ultimo_quarto_log=float(p[150:].sum()),
                                         pico_no_limite=bool(profile.argmax() in [0, 199])))
        print(f"Espectros {i + 1}/25: {row.ensaio}", flush=True)
    rows, individual, summaries, width_rows = [], [], [], []
    archives = {"runs": runs, "labels": labels}
    for (cap, window, method), profiles in powers.items():
        profiles = np.array(profiles)
        archives[f"{cap}_{window}_{method}"] = profiles
        # Also report arithmetic mean of normalized window distributions, per acquisition.
        scenarios = [(window, profiles)]
        if window == "prefixo":
            pp = np.mean([mass(np.array(powers[(cap, w, method)]), grids[cap])
                          for w in ["prefixo", "centro", "sufixo"]], axis=0)
            # Undo quadrature weights so distances() normalizes precisely this mixture.
            quadrature = mass(np.ones(200), grids[cap])
            scenarios.append(("media_tres_janelas", pp / quadrature))
        for scenario, profile in scenarios:
            ds, widths = distances(profile, grids[cap])
            for run, label, val in zip(runs, labels, widths):
                width_rows.append(dict(teto_lambda=cap, janela=scenario, metodo=method,
                                       ensaio=run, padrao=label, largura=val,
                                       m=int(run[1:3]), q=int(run[4:6])))
            for scope, classes in SCOPES.items():
                ix = np.where(np.isin(labels, classes))[0]
                for metric, d in ds.items():
                    summary, sil = separation(d[np.ix_(ix, ix)], labels[ix])
                    rows.append(dict(teto_lambda=cap, janela=scenario, metodo=method, metrica=metric,
                                     escopo=scope, n=len(ix), **summary))
                    for j, val in zip(ix, sil):
                        individual.append(dict(teto_lambda=cap, janela=scenario, metodo=method,
                                               metrica=metric, escopo=scope, ensaio=runs[j],
                                               padrao=labels[j], silhueta=val))
            for c in CLASSES:
                vals = widths[labels == c]
                summaries.append(dict(teto_lambda=cap, janela=scenario, metodo=method,
                                      padrao=c, n=len(vals), mediana_largura=np.median(vals),
                                      minimo=vals.min(), maximo=vals.max()))
    feat = pd.DataFrame(features)
    stats = pd.DataFrame(rows)
    feat.to_csv(OUT / "descritores_por_janela.csv", index=False)
    pd.DataFrame(windows).to_csv(OUT / "janelas.csv", index=False)
    stats.to_csv(OUT / "separabilidade_todas_variantes.csv", index=False)
    pd.DataFrame(individual).to_csv(OUT / "silhueta_por_aquisicao.csv", index=False)
    pd.DataFrame(summaries).to_csv(OUT / "largura_por_classe.csv", index=False)
    widths_frame = pd.DataFrame(width_rows)
    widths_frame.to_csv(OUT / "largura_todas_representacoes.csv", index=False)
    width_controls(widths_frame)
    pd.DataFrame(checks).to_csv(OUT / "verificacao_wps_independente.csv", index=False)
    np.savez_compressed(OUT / "espectros_janelas.npz", **archives)
    variability = feat.groupby(["teto_lambda", "metodo", "ensaio", "padrao"]).apply(
        lambda g: pd.Series({"amplitude_largura_tres_janelas": g[g.janela != "ancora"].largura.max()
                            - g[g.janela != "ancora"].largura.min(),
                            "largura_ancora": g[g.janela == "ancora"].largura.iloc[0]}), include_groups=False).reset_index()
    variability.to_csv(OUT / "variacao_intra_aquisicao.csv", index=False)
    composition()
    plots(stats, feat)
    save_json(OUT / "fontes_sha256.json", INPUTS)
    save_json(OUT / "execucao.json", {"status": "PASS", "aquisicoes": 25, "janelas_por_aquisicao": 4,
              "bandas": [[50, 512], [50, 1024]], "n_avaliacoes_descritivas": len(stats),
              "n_classificadores": 0, "n_testes_significancia": 0,
              "max_erro_wps_ancora": max(v["erro_wps_ancora"] for v in checks),
              "python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__,
              "script_sha256": hashlib.file_digest(Path(__file__).open("rb"), "sha256").hexdigest()})
    print(stats[(stats.janela == "media_tres_janelas") & (stats.teto_lambda == 512)].to_string(index=False))


def width_controls(frame):
    """Adendo exploratório de robustez, registrado após a triagem inicial.

    Sem p-valores; índices descritivos de todos os métodos/bandas/janelas.
    Não confundir retirada de aquisição com validação preditiva.
    """
    strata, influence = [], []
    for (cap, window, method), group in frame.groupby(["teto_lambda", "janela", "metodo"]):
        group = group[group.padrao.isin(CLASSES[:2])]
        subsets = [("todos", group)]
        for name in ["m", "q"]:
            subsets += [(f"{name}={val}", d) for val, d in group.groupby(name)
                        if len(d.padrao.unique()) == 2]
        for stratum, d in subsets:
            a = d[d.padrao == "churn"].largura.to_numpy()
            b = d[d.padrao == "wavy annular"].largura.to_numpy()
            probability = float(((a[:, None] > b) + 0.5 * (a[:, None] == b)).mean())
            strata.append(dict(teto_lambda=cap, janela=window, metodo=method, estrato=stratum,
                               n_churn=len(a), n_anular=len(b), mediana_churn=np.median(a),
                               mediana_anular=np.median(b), probabilidade_churn_maior=probability,
                               cliff_delta=2 * probability - 1))
        for remove in [None] + group.ensaio.tolist():
            d = group if remove is None else group[group.ensaio != remove]
            matrix = squareform(pdist(d.largura.to_numpy()[:, None]))
            result, _ = separation(matrix, d.padrao.to_numpy())
            influence.append(dict(teto_lambda=cap, janela=window, metodo=method,
                                  removido=remove or "nenhum", n=len(d), **result))
    pd.DataFrame(strata).to_csv(OUT / "largura_estratos_operacionais.csv", index=False)
    pd.DataFrame(influence).to_csv(OUT / "largura_influencia_aquisicoes.csv", index=False)


def plots(stats, feat):
    plt.rcParams.update({"font.size": 10})
    fig, axes = plt.subplots(1, 3, figsize=(15, 6), constrained_layout=True)
    for ax, (scope, classes) in zip(axes, SCOPES.items()):
        d = stats[(stats.teto_lambda == 512) & (stats.escopo == scope)].copy()
        d["representacao"] = d.metodo + " / " + d.metrica
        tab = d.pivot(index="representacao", columns="janela", values="silhueta_macro")
        tab = tab[["prefixo", "centro", "sufixo", "media_tres_janelas", "ancora"]]
        im = ax.imshow(tab, cmap="RdBu", vmin=-0.6, vmax=0.6, aspect="auto")
        ax.set_xticks(range(len(tab.columns)), ["início", "centro", "fim", "média 3", "âncora"], rotation=45)
        ax.set_yticks(range(len(tab)), tab.index if ax == axes[0] else [])
        ax.set_title(scope.replace("_", " "))
        for i in range(len(tab)):
            for j in range(len(tab.columns)):
                ax.text(j, i, f"{tab.iloc[i,j]:.2f}", ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=axes, label="Silhueta macro (maior: mais separação descritiva)", shrink=0.8)
    fig.suptitle("50–512 px; mesmas aquisições e regra de janela; nenhuma classificação")
    fig.savefig(OUT / "separabilidade_janelas.png", dpi=180)
    fig.savefig(OUT / "separabilidade_janelas.pdf")
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), constrained_layout=True)
    for ax, method in zip(axes, ["wps", "wps_ponderada", "welch_logenergia"]):
        d = feat[(feat.teto_lambda == 512) & (feat.metodo == method) & (feat.janela != "ancora")]
        for i, c in enumerate(CLASSES):
            for _, g in d[d.padrao == c].groupby("ensaio"):
                ax.plot([i - .18, i, i + .18], g.set_index("janela").loc[["prefixo", "centro", "sufixo"], "largura"],
                        marker="o", alpha=.65, linewidth=.8, markersize=3)
        ax.set_xticks(range(4), ["churn", "anular", "plug", "slug"])
        ax.set_title(method)
        ax.set_ylabel("Largura interquartil em ln(λ)")
        ax.grid(alpha=.2)
    fig.suptitle("Cada linha liga início, centro e fim de uma aquisição — sensibilidades dependentes")
    fig.savefig(OUT / "largura_variacao_por_aquisicao.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
