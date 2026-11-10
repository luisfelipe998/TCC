"""Recalcula somente o par piloto da v16 e apresenta os resultados como na v16.3.

Execução a partir da raiz do projeto:
    python3 src2/v16_spatial2_18_71/scripts/gerar_piloto.py
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
os.environ.setdefault("MPLCONFIGDIR", "/tmp/pearson-v16-comparacao-mpl")

import numpy as np
from PIL import Image
from scipy.stats import pearsonr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.text import Text
from matplotlib.ticker import FuncFormatter, MaxNLocator, NullLocator

import wavelet_core as core
from extraction_config import ROW_START, ROW_END, ROW_STOP, ROW_COUNT

ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT / "_OLD"
INTER = ROOT / "intermediario"
DATA = ROOT.parent / "data_processed"
RUNS = ("m05q01_001", "m06q01_001")
WAVES = np.geomspace(50, 3000, 200)
TEAL, ORANGE, INK = "#087f8c", "#a05f30", "#263c48"
WAVELET_TITLES = {
    "wxx": "Espectro de potência wavelet (WPS) · A",
    "wyy": "Espectro de potência wavelet (WPS) · B",
    "wxy": "Espectro cruzado wavelet (WCS) · módulo",
    "wco": "Coerência wavelet (WCO)",
}

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10,
    "axes.labelsize": 9, "xtick.labelsize": 9, "ytick.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": .65, "savefig.facecolor": "white",
    "axes.titleweight": "medium", "axes.titlepad": 6,
})


def sha(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def read_csv(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def write_csv(path, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def decimal(value, _=None):
    return f"{value:g}".replace(".", ",")


def save_figure(fig, path):
    fig.canvas.draw()
    bbox = fig.get_tightbbox(fig.canvas.get_renderer()).padded(.06)
    scale = 16 / (bbox.width*2.54)
    texts = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text()]
    info = dict(arquivo=str(path.relative_to(ROOT)), largura_insercao_cm=16,
                altura_inserida_cm=bbox.height*2.54*scale,
                fonte_minima_inserida_pt=min(t.get_fontsize() for t in texts)*scale,
                titulo=fig._suptitle.get_text(),
                titulos_paineis=[ax.get_title(loc="left") or ax.get_title() for ax in fig.axes],
                rotulos_x=[ax.get_xlabel() for ax in fig.axes],
                rotulos_y=[ax.get_ylabel() for ax in fig.axes])
    assert info["altura_inserida_cm"] <= 24, info
    assert info["fonte_minima_inserida_pt"] >= 8, info
    fig.savefig(path, dpi=300, bbox_inches="tight", pad_inches=.06)
    plt.close(fig)
    return info


def image_panel(ax, image):
    ax.imshow(image, cmap="gray", vmin=0, vmax=255, aspect="auto",
              extent=[-.5, image.shape[1]-.5, image.shape[0]-.5, -.5], interpolation="nearest")
    ax.set_xlim(-.5, image.shape[1]-.5)
    ax.set_yticks([0, 79], ["0", "79"])
    ax.set_ylabel("Altura do\npanorama (px)")
    ax.set_xlabel("Extensão do panorama (px)", labelpad=2)
    ax.xaxis.set_major_locator(MaxNLocator(5, integer=True))


def panorama_figure(run, image, raw, prepared, common_length):
    fig = plt.figure(figsize=(6.3, 6.4), facecolor="white")
    fig.suptitle(f"plug · {run} · panorama de {len(raw):,} px".replace(",", "."), y=.99, fontsize=10)
    grid = fig.add_gridspec(3, 1, left=.13, right=.98, top=.91, bottom=.115,
                           height_ratios=[.75, 1.3, 1.3], hspace=.6)
    ax = fig.add_subplot(grid[0])
    image_panel(ax, image)
    for boundary in (ROW_START-.5, ROW_END+.5):
        ax.axhline(boundary, color=TEAL, lw=1)
    ax.set_title("Panorama completo · faixa 18–71", loc="left")
    for slot, values, title, ylabel in [
        (grid[1], raw, "Média por coluna · linhas 18–71", "Intensidade média"),
        (grid[2], prepared, "Sinal preparado · panorama completo", "Sinal (σ)"),
    ]:
        ax = fig.add_subplot(slot)
        ax.plot(np.arange(len(raw)), values, color=INK if ylabel.startswith("Intensidade") else TEAL, lw=.8)
        ax.set(xlim=(0, len(raw)-1), xlabel="Extensão do panorama (px)", ylabel=ylabel)
        ax.set_title(title, loc="left")
        ax.xaxis.set_major_locator(MaxNLocator(5, integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_formatter(FuncFormatter(decimal))
        ax.grid(axis="y", color="#e8edef", lw=.6)
        if common_length < len(raw):
            ax.axvspan(common_length-.5, len(raw)-.5, color="#a8afb5", alpha=.18)
            ax.axvline(common_length-.5, color=ORANGE, lw=.85, ls="--")
    if common_length < len(raw):
        fig.text(.56, .018, "Faixa cinza: trecho fora do prefixo usado no par.", ha="center", fontsize=8.5)
    return save_figure(fig, ROOT / "01_panoramas" / f"plug_{run}.png")


def wavelength_axis(ax, label=True):
    ax.set_yscale("log")
    ax.set_ylim(50, 3000)
    ax.set_yticks([50, 250, 1000, 3000], ["50", "250", "1000", "3000"])
    ax.yaxis.set_minor_locator(NullLocator())
    if label:
        ax.set_ylabel("Comprimento de onda\nλ (px)")
    else:
        ax.tick_params(axis="y", labelleft=False, length=0)


def heatmap(fig, slot, maps, key, title, upper, means=None, mean_upper=None):
    sub = slot.subgridspec(1, 3 if means is not None else 2,
                          width_ratios=[4.65, .14, 1.40] if means is not None else [6.58, .14],
                          wspace=.38 if means is not None else .065)
    ax, bar = fig.add_subplot(sub[0, 0]), fig.add_subplot(sub[0, 1])
    values = maps[key] if key == "wco" else np.abs(maps[key])
    lags = maps["lags"]
    step = max(1, int(np.ceil(len(lags)/1800)))
    index = np.unique(np.r_[np.arange(0, len(lags), step), len(lags)//2, len(lags)-1])
    picture = ax.pcolormesh(lags[index], WAVES, values[:, index], shading="auto",
                            cmap="magma" if key == "wco" else "viridis", norm=Normalize(0, upper))
    ax.set_xlim(lags[0], lags[-1])
    ax.set_title(title, loc="left", fontsize=9)
    ax.set_xlabel("Lag (px)", labelpad=2)
    ax.xaxis.set_major_locator(MaxNLocator(4, integer=True))
    wavelength_axis(ax)
    cb = fig.colorbar(picture, cax=bar, format=FuncFormatter(decimal))
    cb.ax.tick_params(labelsize=8.5, length=2, pad=2)
    cb.locator = MaxNLocator(4)
    cb.update_ticks()
    if means is not None:
        ax = fig.add_subplot(sub[0, 2])
        ax.plot(means, WAVES, color=TEAL, lw=1.1)
        ax.fill_betweenx(WAVES, 0, means, color=TEAL, alpha=.08)
        wavelength_axis(ax, label=False)
        ax.set_xlim(0, mean_upper)
        ax.xaxis.set_major_locator(MaxNLocator(2))
        ax.xaxis.set_major_formatter(FuncFormatter(decimal))
        ax.set_xlabel("u.a.", labelpad=2)
        ax.set_title("Potência média", fontsize=8.5)
        ax.grid(axis="x", color="#edf0f0", lw=.5)


def pair_figure(images, x, y, maps, *, runs=RUNS, pattern="plug", destination=None):
    """Apresentação comum aos pares intra e entre padrões."""
    n = len(x)
    fig = plt.figure(figsize=(6.3, 9.25), facecolor="white")
    fig.suptitle(f"{pattern} · {runs[0]} × {runs[1]} · {n:,} px".replace(",", "."), y=.995, fontsize=10)
    grid = fig.add_gridspec(7, 1, left=.095, right=.945, top=.945, bottom=.035,
                           height_ratios=[.52, .52, .78, 1.16, 1.16, 1.12, 1.23], hspace=.87)
    for i, (role, run, image) in enumerate(zip(("A", "B"), runs, images)):
        ax = fig.add_subplot(grid[i])
        image_panel(ax, image[:, :n])
        ax.set_title(f"{role} - {run} · prefixo 0–{n-1} px", loc="left", fontsize=9)
    ax = fig.add_subplot(grid[2])
    ax.plot(np.arange(n), x, color=TEAL, lw=.8, label="A")
    ax.plot(np.arange(n), y, color=ORANGE, lw=.8, ls="--", label="B")
    ax.set_title("Sinais normalizados", loc="left", fontsize=9)
    ax.set(xlim=(0, n-1), ylabel="Sinal (σ)")
    ax.set_xlabel("Extensão do panorama (px)", labelpad=2)
    ax.xaxis.set_major_locator(MaxNLocator(5, integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(3))
    ax.yaxis.set_major_formatter(FuncFormatter(decimal))
    ax.legend(loc="upper right", ncol=2, frameon=False, fontsize=8.5)
    power_max = max(float(np.abs(maps[k]).max()) for k in ("wxx", "wyy"))
    means = [np.abs(maps[k]).mean(axis=1) for k in ("wxx", "wyy")]
    mean_max = max(float(v.max()) for v in means)*1.05
    heatmap(fig, grid[3], maps, "wxx", WAVELET_TITLES["wxx"], power_max, means[0], mean_max)
    heatmap(fig, grid[4], maps, "wyy", WAVELET_TITLES["wyy"], power_max, means[1], mean_max)
    heatmap(fig, grid[5], maps, "wxy", WAVELET_TITLES["wxy"], float(np.abs(maps["wxy"]).max()))
    heatmap(fig, grid[6], maps, "wco", WAVELET_TITLES["wco"], 1.)
    if destination is None:
        destination = ROOT / "02_wavelets" / f"plug_{runs[0]}_{runs[1]}.png"
    return save_figure(fig, destination)


def render_saved():
    """Atualiza os PNGs a partir dos dados salvos, sem repetir os cálculos."""
    pair_path = INTER/"pares"/f"{RUNS[0]}__{RUNS[1]}.npz"
    protected = [pair_path, INTER/"individuais/descritores_piloto.npz",
                 *(INTER/"entradas"/f"{run}.npz" for run in RUNS),
                 *sorted((ROOT/"04_estatistica/dados").glob("*.csv"))]
    before = {str(path.relative_to(ROOT)): sha(path) for path in protected}
    audit_path = INTER/"auditoria/piloto.json"
    audit = json.loads(audit_path.read_text())
    sources = {row["ensaio"]: row for row in audit["entradas"]}
    with np.load(pair_path) as archive:
        maps = {key: archive[key] for key in archive.files}
    x, y = maps["x"], maps["y"]
    images, presentation = [], []
    for run in RUNS:
        image_path = DATA/run/"__panorama_v11.png"
        assert sha(image_path) == sources[run]["panorama_sha256"]
        with Image.open(image_path) as image:
            pano = np.asarray(image.convert("L"))
        with np.load(INTER/"entradas"/f"{run}.npz") as archive:
            presentation.append(panorama_figure(run, pano, archive["raw"], archive["prepared"], len(x)))
        images.append(pano)
    presentation.append(pair_figure(images, x, y, maps))
    assert before == {str(path.relative_to(ROOT)): sha(path) for path in protected}
    write_json(INTER/"apresentacao/piloto.json", presentation)
    audit["apresentacao"] = presentation
    audit["atualizacao_rotulos"] = dict(script_sha256=sha(__file__), dados_preservados_sha256=before,
                                         calculos_reexecutados=False)
    write_json(audit_path, audit)
    print(json.dumps(dict(status="PASS", figuras=len(presentation), dados_numericos_preservados=True), indent=2))


def zero_time(rxx):
    lags = np.arange(-(len(rxx)//2), len(rxx)//2+1)
    values = []
    for scale in WAVES/core.FOURIER_FACTOR:
        u = lags/scale
        kernel = np.pi**(-.25)/np.sqrt(scale)*np.exp(-u*u/2)*np.cos(6*u)
        values.append(abs(np.sum(rxx*kernel)))
    return np.asarray(values)


def main():
    for folder in ("01_panoramas", "02_wavelets", "04_estatistica/dados", "05_conclusao",
                   "intermediario/entradas", "intermediario/individuais", "intermediario/pares",
                   "intermediario/auditoria", "intermediario/apresentacao"):
        (ROOT/folder).mkdir(parents=True, exist_ok=True)
    assert sha(ROOT/"scripts/wavelet_core.py") == sha(OLD/"scripts/wavelet_core.py")
    inventory = {row["run"]: row for row in read_csv(OLD/"results/input_audit.csv")}
    images, raws, prepared, full_powers, checks = [], [], [], [], []
    old_profiles = [r for r in read_csv(OLD/"results/profiles.csv") if r["method"] == "v16_biased"]
    for run in RUNS:
        row = inventory[run]
        assert row["included"] == "True" and row["pattern"] == "plug"
        path = DATA/run/"__panorama_v11.png"
        assert sha(path) == row["panorama_sha256"]
        with Image.open(path) as image:
            pano = np.asarray(image.convert("L"))
        raw = pano[ROW_START:ROW_STOP].sum(axis=0, dtype=np.uint64)/ROW_COUNT
        with np.load(OLD/"intermediate"/f"{run}_input.npz") as archive:
            np.testing.assert_array_equal(raw, archive["raw"])
        x = core.prepare(raw)
        rxx = core.correlation(x, x, estimator="biased")
        p = zero_time(rxx)
        prior = [r for r in old_profiles if r["run"] == run]
        np.testing.assert_allclose(WAVES, [float(r["wavelength_px"]) for r in prior], atol=1e-12)
        expected = np.array([float(r["zero_lag_wps"]) for r in prior])
        np.testing.assert_allclose(p, expected, atol=1e-11, rtol=1e-11)
        np.savez_compressed(INTER/"entradas"/f"{run}.npz", raw=raw, prepared=x, wavelengths=WAVES, wps_zero=p)
        checks.append(dict(ensaio=run, panorama_sha256=sha(path), comprimento=len(raw),
                           linhas=[ROW_START, ROW_END], erro_raw=0.,
                           erro_wps_zero_completo=float(np.max(abs(p-expected)))))
        images.append(pano); raws.append(raw); prepared.append(x); full_powers.append(p)

    n = min(map(len, raws))
    x, y = [core.prepare(raw[:n]) for raw in raws]
    maps = core.pair(x, y, WAVES, estimator="biased", engine="morlet")
    maps["wavelengths"] = WAVES
    assert maps["denominator_valid"].all()
    assert np.isfinite(maps["wco"]).all()
    assert float(maps["wco"].max()) <= 1+1e-8
    j0 = n-1
    zx, zy = [np.abs(maps[key][:, j0]) for key in ("wxx", "wyy")]
    np.testing.assert_allclose(zx, zero_time(core.correlation(x, x)), atol=1e-11, rtol=1e-11)
    np.testing.assert_allclose(zy, zero_time(core.correlation(y, y)), atol=1e-11, rtol=1e-11)
    # Independent direct sums verify lag orientation and the divisor N.
    rxy = core.correlation(x, y)
    for lag in (-n+1, -101, -1, 0, 1, 101, n-1):
        expected = np.dot(x[lag:], y[:n-lag])/n if lag >= 0 else np.dot(x[:n+lag], y[-lag:])/n
        assert abs(rxy[j0+lag]-expected) < 1e-12

    old_errors = {}
    with np.load(OLD/"intermediate"/f"{RUNS[0]}__{RUNS[1]}_maps.npz") as archive:
        indexes = archive["lags"].astype(int)+j0
        np.testing.assert_array_equal(archive["wavelengths"], WAVES)
        for key in ("wxx", "wyy", "wxy", "wco"):
            np.testing.assert_allclose(maps[key][:, indexes], archive[key], atol=1e-12, rtol=1e-12)
            old_errors[key] = float(np.max(abs(maps[key][:, indexes]-archive[key])))
        for key, old_key in (("wxx", "mean_x_all"), ("wyy", "mean_y_all")):
            np.testing.assert_allclose(np.abs(maps[key]).mean(axis=1), archive[old_key], atol=1e-12, rtol=1e-12)

    weight = np.sqrt(zx*zy)
    active = (zx >= .05*zx.max()) & (zy >= .05*zy.max())
    weighted = np.where(active, weight, 0)
    peak = maps["wco"].max(axis=1)
    score = float(np.sum(weighted*peak)/weighted.sum())
    reference = next(r for r in read_csv(OLD/"results/preprocessing_confirm_pairs.csv")
                     if r["run_x"] == RUNS[0] and r["run_y"] == RUNS[1] and r["variant"] == "baseline")
    assert abs(score-float(reference["weighted_scale_peak_power05"])) < 1e-12
    assert abs(float(maps["wco"].max())-float(reference["max_wco_all"])) < 1e-12

    log_power = np.log(np.asarray(full_powers)+1e-12)
    design = np.column_stack([np.ones(len(WAVES)), np.log(WAVES)])
    coefficients = np.linalg.lstsq(design, log_power.T, rcond=None)[0]
    z = log_power-(design@coefficients).T
    r = float(pearsonr(z[0], z[1]).statistic)
    np.testing.assert_allclose(r, np.corrcoef(z)[0, 1], atol=1e-13)
    np.savez_compressed(INTER/"individuais"/"descritores_piloto.npz", runs=np.asarray(RUNS),
                        wavelengths=WAVES, power=np.asarray(full_powers), residuals=z, pearson=np.corrcoef(z))
    metrics = dict(ensaio_a=RUNS[0], ensaio_b=RUNS[1], padrao="plug", n_a=len(raws[0]), n_b=len(raws[1]),
                   comprimento_comum=n, descarte_a=len(raws[0])-n, descarte_b=len(raws[1])-n,
                   lambda_min=50, lambda_max=3000, escalas=200, lags=2*n-1,
                   wco_maxima_integral=float(maps["wco"].max()),
                   wco_maxima_lag_zero=float(maps["wco"][:, j0].max()),
                   wco_mediana_lag_zero=float(np.median(maps["wco"][:, j0])),
                   wco_media_integral=float(maps["wco"].mean()),
                   escore_ponderado_picos_piso05=score, escalas_piso05=int(active.sum()),
                   pearson_perfis_completos=r)
    write_csv(ROOT/"04_estatistica/dados/metricas_piloto.csv", [metrics])
    write_csv(ROOT/"04_estatistica/dados/perfis_piloto.csv", [dict(ensaio=run, lambda_px=float(lam),
        wps_zero=float(full_powers[i][k]), residuo=float(z[i, k]))
        for i, run in enumerate(RUNS) for k, lam in enumerate(WAVES)])

    # Retain every lag for this single pilot, so its numerical artifact is complete.
    np.savez_compressed(INTER/"pares"/f"{RUNS[0]}__{RUNS[1]}.npz", x=x, y=y, **maps)
    presentation = [panorama_figure(run, image, raw, full, n)
                    for run, image, raw, full in zip(RUNS, images, raws, prepared)]
    presentation.append(pair_figure(images, x, y, maps))
    write_json(INTER/"apresentacao/piloto.json", presentation)

    # Verify that archiving preserved every byte, including old code and figures.
    archived = json.loads((INTER/"auditoria/arquivamento.json").read_text())
    for item in archived["files"]:
        assert sha(OLD/item["path"]) == item["sha256"], item["path"]
    audit = dict(status="PASS", escopo="somente um par piloto; lote completo não executado",
                 entradas=checks, par=list(RUNS), metodo="v16, prefixo comum antes do preparo",
                 nucleo_identico_ao_original=True, erros_mapas_vs_original=old_errors,
                 escore_ponderado_conferido=True, pearson_conferido_scipy_numpy=True,
                 arquivos_antigos_preservados=len(archived["files"]), apresentacao=presentation,
                 script_sha256=sha(__file__), core_sha256=sha(ROOT/"scripts/wavelet_core.py"))
    write_json(INTER/"auditoria/piloto.json", audit)
    number = lambda value: f"{value:.6f}".replace(".", ",")
    report = f'''# Piloto da análise v16

Foi refeito somente o par **m05q01_001 × m06q01_001**, do padrão *plug*, para revisão do formato antes do processamento dos demais ensaios. A apresentação seguiu a organização da v16.3, com os parâmetros e procedimentos da v16 original. Os 583 arquivos anteriores permaneceram em `_OLD`; sua integridade foi conferida por SHA-256.

## Figuras

- [Panorama e sinal de m05q01_001](../01_panoramas/plug_m05q01_001.png).
- [Panorama e sinal de m06q01_001](../01_panoramas/plug_m06q01_001.png).
- [Par piloto: panoramas, sinais, WPS, WCS e WCO](../02_wavelets/plug_m05q01_001_m06q01_001.png).

![Análise wavelet do par piloto](../02_wavelets/plug_m05q01_001_m06q01_001.png)

*Os dois primeiros painéis mostraram os prefixos empregados no par; o terceiro apresentou os sinais normalizados no segmento comum. Os painéis seguintes mostraram os espectros de potência wavelet individuais e suas médias em todos os lags, o módulo do espectro cruzado wavelet e a coerência wavelet. O eixo dos sinais representou a extensão do panorama; o eixo dos mapas representou o lag da correlação.*

## Metodologia executada

Os sinais foram novamente extraídos dos panoramas originais pela média das linhas 18–71, inclusive, com indexação iniciada em zero. O primeiro panorama continha 5.940 colunas; o segundo, 7.798. A comparação empregou as primeiras 5.940 colunas de ambos, seguindo o prefixo comum da v16. O trecho final de 1.858 colunas do segundo panorama ficou fora desta comparação e foi assinalado em cinza na apresentação individual.

Após o corte comum, aplicaram-se a inversão `255 − intensidade`, a remoção da tendência linear, a centralização e a divisão pelo desvio padrão populacional. Não houve deslocamento de eventos entre os sinais. Os sinais individuais de `01_panoramas` foram preparados em seus comprimentos completos; os sinais do par foram preparados novamente após o corte.

Utilizou-se a correlação bilateral com divisor N e a CWT Morlet de energia unitária, com ω₀ = 6, amostragem de 1 px e 200 comprimentos de onda espaçados em escala logarítmica entre 50 e 3.000 px. Foram calculados todos os 11.879 lags, de −5.939 a +5.939 px, sem COI nem máscara geométrica. Preservaram-se as definições alternativas da v16:

- WPS: módulo da CWT da autocorrelação, sem elevar o módulo novamente ao quadrado.
- WCS: CWT complexa da correlação cruzada, apresentada em módulo.
- WCO: módulo quadrático da WCS dividido pelo produto das WPS individuais no lag zero.
- Potência média: média da WPS sobre todos os lags, com peso igual por lag.

As duas WPS compartilharam a escala de cores e os limites da potência média. A WCO foi exibida entre 0 e 1, sem truncamento dos dados. As métricas foram calculadas em todos os lags; a redução horizontal dos mapas ocorreu somente para apresentação. Os dados numéricos completos do piloto foram preservados em `intermediario/pares`.

## Resultados do piloto

| Medida | Valor |
| --- | ---: |
| WCO máxima no plano completo | {number(metrics['wco_maxima_integral'])} |
| WCO máxima no lag zero | {number(metrics['wco_maxima_lag_zero'])} |
| Mediana da WCO no lag zero | {number(metrics['wco_mediana_lag_zero'])} |
| Média da WCO no plano completo | {number(metrics['wco_media_integral'])} |
| Escore ponderado dos máximos por escala, piso de 5% | {number(score)} |
| Pearson dos perfis dos panoramas completos | {number(r)} |

O escore ponderado empregou o mesmo resumo da v16: o máximo da WCO em todos os lags de cada escala, ponderado pela média geométrica das WPS no lag zero. Foram retidas as 154 escalas em que ambas as WPS atingiram pelo menos 5% de seus respectivos máximos. Esse piso foi aplicado ao resumo ponderado; todos os mapas conservaram as 200 escalas.

O Pearson foi calculado entre os resíduos de `log(WPS no lag zero + 10⁻¹²)` após ajuste linear contra `log(λ)`, seguindo a definição estatística da v16.3. Cada descritor utilizou seu panorama completo e a grade nativa da v16, de 50 a 3.000 px. Portanto, seu domínio diferiu daquele da comparação anterior entre versões, que utilizou 50–1.024 px e 25 ensaios comuns. O valor alto deste par indicou semelhança da forma residual dos perfis espectrais; não implicou WCO elevada em todas as escalas e lags.

## Verificação e alcance

Os sinais extraídos coincidiram exatamente com as entradas arquivadas. Os mapas recalculados coincidiram com os valores salvos da v16 nas mesmas posições de lag, com erro máximo zero para WPS, WCS e WCO. O escore ponderado e o máximo integral também coincidiram. A implementação da correlação foi conferida por somas diretas; a WPS no lag zero, por uma implementação temporal independente; o Pearson, por SciPy e NumPy.

Este piloto verificou um par e o formato de apresentação. Não foram calculados os demais pares nem produzidas conclusões sobre separação global entre padrões. A execução completa ficou pendente da revisão deste piloto pelo usuário.

Dados: [métricas](../04_estatistica/dados/metricas_piloto.csv), [perfis espectrais](../04_estatistica/dados/perfis_piloto.csv), [auditoria](../intermediario/auditoria/piloto.json) e [reprodução](../scripts/README.txt).
'''
    (ROOT/"05_conclusao/Piloto.md").write_text(report, encoding="utf-8")
    print(json.dumps(dict(status="PASS", metricas=metrics, apresentacao=presentation,
                         arquivos_antigos_preservados=len(archived["files"])), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--somente-figuras", action="store_true", help="Regenera os PNGs usando os dados existentes.")
    args = parser.parse_args()
    render_saved() if args.somente_figuras else main()
