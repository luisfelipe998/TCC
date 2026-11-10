"""Apresentação local dos panoramas e campos wavelet da análise vigente."""
from runtime import OUTPUT_ROOT, reference
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.text import Text
from matplotlib.ticker import FuncFormatter, MaxNLocator, NullLocator

ROOT = OUTPUT_ROOT / "04_analise_simplificada"
WAVES = np.geomspace(50, 3000, 200)
ROW_START, ROW_END = 18, 71

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


def decimal(value, _=None):
    return f"{value:g}".replace(".", ",")


def save_figure(fig, path):
    fig.canvas.draw()
    bbox = fig.get_tightbbox(fig.canvas.get_renderer()).padded(.06)
    scale = 16 / (bbox.width*2.54)
    texts = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text()]
    info = dict(arquivo=reference(path), largura_insercao_cm=16,
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


def panorama_figure(run, image, raw, prepared, *, pattern):
    common_length = len(raw)
    fig = plt.figure(figsize=(6.3, 6.4), facecolor="white")
    fig.suptitle(f"{pattern} · {run} · panorama de {len(raw):,} px".replace(",", "."), y=.99, fontsize=10)
    grid = fig.add_gridspec(3, 1, left=.13, right=.98, top=.91, bottom=.115,
                           height_ratios=[.75, 1.3, 1.3], hspace=.6)
    ax = fig.add_subplot(grid[0])
    image_panel(ax, image)
    for boundary in (ROW_START-.5, ROW_END+.5):
        ax.axhline(boundary, color=TEAL, lw=1)
    ax.set_title("Panorama completo · faixa 18–71", loc="left")
    for slot, values, title, ylabel in [
        (grid[1], raw, "Média por coluna · linhas 18–71", "Intensidade média"),
        (grid[2], prepared, "Sinal normalizado", "Amplitude\npadronizada"),
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
    return save_figure(fig, ROOT / "01_panoramas" / f"{pattern.replace(chr(32), chr(95))}_{run}.png")


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
    step = 1 if "mean_x_all" in maps else max(1, int(np.ceil(len(lags)/1800)))
    index = np.unique(np.r_[np.arange(0, len(lags), step), len(lags)//2, len(lags)-1])
    picture = ax.pcolormesh(lags[index], WAVES, values[:, index], shading="auto",
                            cmap="magma" if key == "wco" else "viridis", norm=Normalize(0, upper))
    ax.set_xlim(lags[0], lags[-1])
    ax.set_title(title, loc="left", fontsize=9)
    ax.set_xlabel("Atraso (pixels)", labelpad=2)
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


def pair_figure(images, x, y, maps, *, runs, pattern, destination):
    """Apresentação comum aos pares intra e entre padrões."""
    n = len(x)
    fig = plt.figure(figsize=(6.3, 9.25), facecolor="white")
    title = f"{pattern} · {runs[0]} × {runs[1]} · {n:,} px".replace(",", ".")
    if len(title) > 65:
        title = f"{pattern}\n{runs[0]} × {runs[1]} · {n:,} px".replace(",", ".")
    fig.suptitle(title, y=.995, fontsize=10)
    grid = fig.add_gridspec(7, 1, left=.095, right=.945, top=.925 if "\n" in title else .945, bottom=.035,
                           height_ratios=[.52, .52, .78, 1.16, 1.16, 1.12, 1.23], hspace=.87)
    for i, (role, run, image) in enumerate(zip(("A", "B"), runs, images)):
        ax = fig.add_subplot(grid[i])
        image_panel(ax, image[:, :n])
        ax.set_title(f"{role} - {run} · prefixo 0–{n-1} px", loc="left", fontsize=9)
    ax = fig.add_subplot(grid[2])
    ax.plot(np.arange(n), x, color=TEAL, lw=.8, label="A")
    ax.plot(np.arange(n), y, color=ORANGE, lw=.8, ls="--", label="B")
    ax.set_title("Sinais preparados", loc="left", fontsize=9)
    ax.set(xlim=(0, n-1), ylabel="Amplitude\npadronizada")
    ax.set_xlabel("Extensão do panorama (px)", labelpad=2)
    ax.xaxis.set_major_locator(MaxNLocator(5, integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(3))
    ax.yaxis.set_major_formatter(FuncFormatter(decimal))
    ax.legend(loc="upper right", ncol=2, frameon=False, fontsize=8.5)
    power_max = max(float(maps["max_wxx"]), float(maps["max_wyy"]))
    means = [maps["mean_x_all"], maps["mean_y_all"]]
    mean_max = max(float(v.max()) for v in means)*1.05
    heatmap(fig, grid[3], maps, "wxx", WAVELET_TITLES["wxx"], power_max, means[0], mean_max)
    heatmap(fig, grid[4], maps, "wyy", WAVELET_TITLES["wyy"], power_max, means[1], mean_max)
    heatmap(fig, grid[5], maps, "wxy", WAVELET_TITLES["wxy"], float(maps["max_wxy"]))
    heatmap(fig, grid[6], maps, "wco", WAVELET_TITLES["wco"], 1.)
    return save_figure(fig, destination)
