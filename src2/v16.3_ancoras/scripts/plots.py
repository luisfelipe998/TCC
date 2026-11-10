"""Compact portrait PNGs for insertion into a document; no page or PDF export."""
from common import ANCHORS, CFG, INTER, SCRIPTS, WAVELETS, load_image, slug, write_json
from pair_identity import pair_key
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, to_rgba
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from matplotlib.ticker import FuncFormatter, MaxNLocator, NullLocator

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10,
    "axes.labelsize": 9, "xtick.labelsize": 9, "ytick.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": .65, "pdf.fonttype": 42, "savefig.facecolor": "white",
    "axes.titleweight": "medium", "axes.titlepad": 6,
})
TEAL = "#087f8c"
ORANGE = "#d17c1e"
INK = "#263c48"
POWER_CMAP = plt.get_cmap(CFG["figuras"]["paleta_wps_wcs"])
COHERENCE_CMAP = plt.get_cmap(CFG["figuras"]["paleta_wco"])
STYLE = CFG["figuras"] | json.loads((SCRIPTS / "plot_config.json").read_text())


def decimal(value, _=None):
    return f"{value:g}".replace(".", ",")


def figure(kind):
    return plt.figure(figsize=STYLE[kind + "_inches"], facecolor="white")


def save_figure(fig, stem):
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text()]
    bbox = fig.get_tightbbox(renderer).padded(.06)
    scale = STYLE["largura_insercao_cm"] / (bbox.width * 2.54)
    diagnostics = dict(
        arquivo=stem.name + ".png", largura_insercao_cm=STYLE["largura_insercao_cm"],
        altura_inserida_cm=bbox.height * 2.54 * scale,
        fonte_minima_inserida_pt=min(t.get_fontsize() for t in texts) * scale,
        titulo=fig._suptitle.get_text(),
        titulos_paineis=[ax.get_title(loc="left") or ax.get_title() for ax in fig.axes],
        rotulos=[ax.get_ylabel() for ax in fig.axes])
    assert diagnostics["altura_inserida_cm"] <= STYLE["altura_maxima_insercao_cm"], diagnostics
    assert diagnostics["fonte_minima_inserida_pt"] >= STYLE["fonte_minima_inserida_pt"], diagnostics
    diagnostic_dir = INTER / "apresentacao"
    diagnostic_dir.mkdir(exist_ok=True)
    write_json(diagnostic_dir / f"{stem.name}.json", diagnostics)
    fig.savefig(stem.with_suffix(".png"), dpi=STYLE["dpi"],
                bbox_inches="tight", pad_inches=.06)
    # At 96 dpi this preview approximates the text size when inserted at 16 cm.
    preview_dir = diagnostic_dir / "previas_16cm"
    preview_dir.mkdir(exist_ok=True)
    fig.savefig(preview_dir / f"{stem.name}.png", dpi=96*scale,
                bbox_inches="tight", pad_inches=.06)
    plt.close(fig)


def show_image(ax, pano, left, right):
    ax.imshow(pano, cmap="gray", vmin=0, vmax=255, aspect="auto",
              extent=[left-.5, right-.5, pano.shape[0]-.5, -.5], interpolation="nearest")
    ax.set_xlim(left-.5, right-.5)
    ax.set_yticks([0, 40, 79], ["0", "40", "79"])
    ax.set_ylabel("Linha")
    ax.tick_params(length=3)


def mark_anchor(ax, record, origin=0):
    ax.axvline(record["ancora"]-origin, color=ORANGE, lw=1.9, ls="--", zorder=5)


def anchor_figure(record, pano, raw):
    fig = figure("ancoras")
    pattern, run, center = record["padrao"], record["aquisicao"], record["ancora"]
    start, stop = record["segmento_inicio"], record["segmento_fim"]
    length = stop - start
    fig.suptitle(f"{pattern} · Âncora · {run} · coluna {center} px",
                 y=.99, fontsize=10, fontweight="medium")
    grid = fig.add_gridspec(3, 1, left=.11, right=.98, top=.915, bottom=.055,
                           height_ratios=[.72, 1.25, 1.40], hspace=.43)
    overview = fig.add_subplot(grid[0])
    show_image(overview, pano, 0, len(raw))
    overview.add_patch(Rectangle((start-.5, -.5), length, pano.shape[0],
                                 facecolor=to_rgba(TEAL, .18),
                                 edgecolor=to_rgba(TEAL, .70), lw=1.3))
    overview.axvline(center, color=ORANGE, lw=1.9, ls="--")
    overview.set_title("Panorama completo", loc="left")
    overview.set_xlabel("Coluna (px)")
    overview.xaxis.set_major_locator(MaxNLocator(6, integer=True))

    zoom = fig.add_subplot(grid[1])
    show_image(zoom, pano[:, start:stop], start, stop)
    first_row, last_row = CFG["linhas"]
    # Pixel edges delimit the inclusive rows used to calculate the mean.
    for boundary in (first_row-.5, last_row+.5):
        zoom.axhline(boundary, color=TEAL, lw=1.1, zorder=4)
    mark_anchor(zoom, record)
    zoom.set_title(f"Contexto · {length:,} px".replace(",", "."), loc="left")
    zoom.set_title(f"{first_row}–{last_row} px", loc="right", fontsize=9)
    zoom.set_xlabel("Coluna (px)")

    profile = fig.add_subplot(grid[2], sharex=zoom)
    columns = np.arange(start, stop)
    profile.plot(columns, raw[start:stop], color=INK, lw=.8)
    mark_anchor(profile, record)
    profile.set_title(f"Média por coluna · {first_row}–{last_row} px", loc="left")
    profile.set_ylabel("Intensidade média")
    profile.set_xlabel("Coluna (px)")
    profile.xaxis.set_major_locator(MaxNLocator(5, integer=True))
    profile.yaxis.set_major_locator(MaxNLocator(5))
    profile.grid(axis="y", color="#e8edef", lw=.6)
    save_figure(fig, ANCHORS / f"{slug(pattern)}_{record['anchor_id']}")


def wavelength_axis(ax, wavelengths, label=True):
    ax.set_yscale("log")
    ax.set_ylim(wavelengths[0], wavelengths[-1])
    ticks = [x for x in (64, 256, 1024, 4096) if x <= wavelengths[-1]]
    if wavelengths[-1] / ticks[-1] >= 1.75:
        ticks.append(int(wavelengths[-1]))
    ax.set_yticks(ticks, [str(t) for t in ticks])
    ax.yaxis.set_minor_locator(NullLocator())
    if label:
        ax.set_ylabel("Comprimento de onda\nλ (px)")
    else:
        ax.tick_params(axis="y", labelleft=False, length=0)


def heatmap(fig, slot, maps, key, title, upper, mean=None, power_upper=None):
    grid = slot.subgridspec(1, 3 if mean is not None else 2,
                           width_ratios=[4.65, .14, 1.40] if mean is not None else [6.58, .14],
                           wspace=.38 if mean is not None else .065)
    ax, bar = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])
    values = maps[key] if key == "wco" else np.abs(maps[key])
    lags, wavelengths = maps["lags"], maps["wavelengths"]
    zero = len(maps["x"]) - 1
    # Decimation is display-only; all metrics use every integer lag.
    indexes = np.unique(np.r_[np.arange(0, len(lags), 3), len(lags)-1, zero,
                               zero-CFG["wavelet"]["lag_resumo"], zero+CFG["wavelet"]["lag_resumo"]])
    cmap = COHERENCE_CMAP if key == "wco" else POWER_CMAP
    picture = ax.pcolormesh(lags[indexes], wavelengths, values[:, indexes], shading="auto",
                            cmap=cmap, norm=Normalize(0, upper), rasterized=True)
    edge = zero - np.sqrt(2) * maps["scales"]
    for sign in (-1, 1):
        ax.plot(sign*edge, wavelengths, color="#70777a", ls=":", lw=.7, alpha=.85)
    if key == "wco":
        for sign in (-1, 1):
            ax.axvline(sign*CFG["wavelet"]["lag_resumo"], color="#6c8190", lw=.6, ls="--")
    ax.set_xlim(lags[0], lags[-1])
    ax.set_title(title, loc="left", fontsize=9)
    ax.set_xlabel("Lag (px)", labelpad=2)
    ax.xaxis.set_major_locator(MaxNLocator(4, integer=True))
    wavelength_axis(ax, wavelengths)
    cbar = fig.colorbar(picture, cax=bar, format=FuncFormatter(decimal))
    cbar.ax.tick_params(labelsize=8.5, length=2, pad=2)
    cbar.locator = MaxNLocator(4)
    cbar.update_ticks()
    if mean is not None:
        power = fig.add_subplot(grid[0, 2])
        power.plot(mean, wavelengths, color=TEAL, lw=1.1)
        power.fill_betweenx(wavelengths, 0, mean, color=TEAL, alpha=.08)
        wavelength_axis(power, wavelengths, label=False)
        power.set_xlim(0, power_upper)
        power.xaxis.set_major_locator(MaxNLocator(2))
        power.xaxis.set_major_formatter(FuncFormatter(decimal))
        power.set_xlabel("u.a.", labelpad=2)
        power.set_title("Potência média", fontsize=8.5)
        power.grid(axis="x", color="#edf0f0", lw=.5)


def pair_figure(pair, anchors, maps, metrics, label=None):
    fig = figure("pares")
    pattern = pair["padrao"]
    length_label = f"{pair['comprimento']:,}".replace(",", ".")
    fig.suptitle(f"{pattern} · {pair['run_x']} × {pair['run_y']} · {length_label} px",
                 y=.995, fontsize=10)
    grid = fig.add_gridspec(7, 1, left=.095, right=.945, top=.945, bottom=.035,
                           height_ratios=[.52, .52, .78, 1.16, 1.16, 1.12, 1.23], hspace=.87)
    for index, (role, anchor_id) in enumerate((("A", pair["anchor_x"]), ("B", pair["anchor_y"]))):
        record = anchors[anchor_id]
        start, stop, center = record["segmento_inicio"], record["segmento_fim"], record["ancora"]
        pano, _ = load_image(record["aquisicao"])
        ax = fig.add_subplot(grid[index])
        show_image(ax, pano[:, start:stop], start-center, stop-center)
        ax.set_yticks([0, 79], ["0", "79"])
        mark_anchor(ax, record, origin=center)
        ax.set_title(f"{role} - {record['aquisicao']}", loc="left", fontsize=9)
        ax.set_xlabel("Posição relativa (px)", labelpad=2)
        ax.xaxis.set_major_locator(MaxNLocator(5, integer=True))
    ax = fig.add_subplot(grid[2])
    positions = np.arange(len(maps["x"])) - len(maps["x"])//2
    ax.plot(positions, maps["x"], color=TEAL, lw=.8, label="A")
    ax.plot(positions, maps["y"], color="#a05f30", lw=.8, ls="--", label="B")
    ax.axvline(0, color="#939b9d", lw=.55, ls=":")
    ax.set_title("Sinais alinhados", loc="left", fontsize=9)
    ax.set_xlabel("Posição relativa (px)", labelpad=2)
    ax.set_ylabel("Sinal (σ)")
    ax.set_xlim(positions[0], positions[-1])
    ax.yaxis.set_major_locator(MaxNLocator(3))
    ax.yaxis.set_major_formatter(FuncFormatter(decimal))
    ax.legend(loc="upper right", ncol=2, frameon=False, fontsize=8.5)
    power_max = max(float(np.abs(maps["wxx"]).max()), float(np.abs(maps["wyy"]).max()))
    mean_max = 1.05 * max(float(maps["power_mean_x"].max()), float(maps["power_mean_y"].max()))
    heatmap(fig, grid[3], maps, "wxx", "WPS · A", power_max, maps["power_mean_x"], mean_max)
    heatmap(fig, grid[4], maps, "wyy", "WPS · B", power_max, maps["power_mean_y"], mean_max)
    heatmap(fig, grid[5], maps, "wxy", "WCS · módulo", float(np.abs(maps["wxy"]).max()))
    heatmap(fig, grid[6], maps, "wco", "WCO",
            max(1.0, float(np.nanmax(maps["wco"]))))
    name = f"{pair_key(pair)}_{label or 'previa'}"
    save_figure(fig, (WAVELETS if label else INTER / "previas") / name)
