"""Ilustração sintética da extração dos dois perfis de potência wavelet.

Reprodução, a partir da raiz do repositório:
    python3 -B TCC/figuras/metodologia/gerar_perfis_wavelet.py

Não utiliza ensaios experimentais. A potência é o módulo da CWT da
autocorrelação, conforme a implementação efetivamente usada no trabalho.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import sys

sys.dont_write_bytecode = True
os.environ.setdefault("MPLCONFIGDIR", "/tmp/tcc-perfis-wavelet-mpl")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
import scipy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, NullLocator

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CORE_PATH = ROOT / "src2/v16_spatial2_18_71/scripts/wavelet_core.py"
PREFIX = HERE / "perfis_wavelet"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    core_hash = digest(CORE_PATH)
    spec = importlib.util.spec_from_file_location("perfis_didaticos_core", CORE_PATH)
    core = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(core)

    n = 8192
    positions = np.arange(n, dtype=float)
    envelope = np.exp(-0.5 * ((positions - 4200.0) / 650.0) ** 2)
    synthetic_raw = 127.5 + 30.0 * (
        0.7 * np.sin(2 * np.pi * positions / 760.0)
        + 2.2 * envelope * np.sin(2 * np.pi * positions / 210.0 + 0.4)
        + 0.3 * np.cos(2 * np.pi * positions / 1450.0)
    )
    signal = core.prepare(synthetic_raw)
    wavelengths = np.geomspace(50.0, 3000.0, 200)
    correlation = core.correlation(signal, signal, estimator="biased")
    coefficients, scales = core.transform_correlation(
        correlation, wavelengths, engine="morlet"
    )
    power = np.abs(coefficients)
    lags = np.arange(-n + 1, n)
    zero_index = n - 1
    profile_zero = power[:, zero_index].copy()
    profile_mean = power.mean(axis=1)
    assert np.isfinite(power).all() and np.all(power >= 0)
    assert lags[zero_index] == 0
    np.testing.assert_array_equal(profile_zero, np.abs(coefficients[:, zero_index]))
    np.testing.assert_array_equal(profile_mean, np.sum(power, axis=1) / len(lags))

    data_path = HERE / "perfis_wavelet_dados.npz"
    np.savez_compressed(
        data_path, positions=positions, synthetic_raw=synthetic_raw,
        signal=signal, envelope=envelope, autocorrelation=correlation,
        lags=lags, wavelengths=wavelengths, scales=scales, power=power,
        profile_zero=profile_zero, profile_mean=profile_mean,
        zero_index=zero_index,
    )

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.titlesize": 9, "axes.labelsize": 9,
        "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False,
    })
    width_cm, height_cm = 16.0, 7.8
    fig = plt.figure(figsize=(width_cm / 2.54, height_cm / 2.54))
    grid = fig.add_gridspec(
        1, 3, width_ratios=[2.45, 1, 1], left=0.105, right=0.985,
        bottom=0.28, top=0.73, wspace=0.30,
    )
    ax_map = fig.add_subplot(grid[0, 0])
    ax_zero = fig.add_subplot(grid[0, 1], sharey=ax_map)
    ax_mean = fig.add_subplot(grid[0, 2], sharey=ax_map)
    # Only the drawing is reduced in the lag direction. All saved data and
    # both profiles use every one of the 16383 lag samples.
    step = 8
    drawing_indices = np.unique(np.r_[np.arange(0, len(lags), step), zero_index, len(lags)-1])
    mesh = ax_map.pcolormesh(
        lags[drawing_indices], wavelengths, power[:, drawing_indices],
        shading="nearest", cmap="viridis", rasterized=True, vmin=0,
        vmax=float(power.max()),
    )
    zero_color = "#c34a25"
    mean_color = "#087f8c"
    ax_map.axvline(0, color="white", linewidth=2.6, alpha=0.9)
    ax_map.axvline(0, color=zero_color, linewidth=1.25)
    ax_map.set_title("Potência wavelet\nda autocorrelação", pad=9)
    ax_map.set_xlabel("Atraso (pixels)", labelpad=3)
    ax_map.set_ylabel("Comprimento de onda (pixels)", labelpad=3)
    ax_map.set_xticks([-8000, 0, 8000], ["−8000", "0", "8000"])
    ax_map.set_xlim(lags[0], lags[-1])
    ax_zero.plot(profile_zero, wavelengths, color=zero_color, linewidth=1.6)
    ax_zero.fill_betweenx(wavelengths, 0, profile_zero, color=zero_color, alpha=0.08)
    ax_zero.set_title("Potência no\natraso zero\nComparação por\nPearson", color=zero_color, fontsize=8.5, pad=9)
    ax_mean.plot(profile_mean, wavelengths, color=mean_color, linewidth=1.6)
    ax_mean.fill_betweenx(wavelengths, 0, profile_mean, color=mean_color, alpha=0.08)
    ax_mean.set_title("Potência média\nsobre os atrasos\nCurva lateral\ndos mapas", color=mean_color, fontsize=8.5, pad=9)
    for axis in (ax_map, ax_zero, ax_mean):
        axis.set_yscale("log")
        axis.set_ylim(50, 3000)
        axis.set_yticks([50, 250, 1000, 3000], ["50", "250", "1000", "3000"])
        axis.yaxis.set_minor_locator(NullLocator())
    for axis in (ax_zero, ax_mean):
        axis.set_xlabel("Potência (u.a.)", labelpad=3)
        axis.tick_params(axis="y", labelleft=False)
        axis.grid(axis="y", color="#dce4e8", linewidth=0.5)
        axis.set_xlim(left=0)
        axis.xaxis.set_major_locator(plt.MaxNLocator(3))
        axis.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}".replace(".", ",")))
    # A common horizontal scale makes the reduction caused by averaging
    # visible without changing either of the numerical profiles.
    xmax = 1.04 * float(max(profile_zero.max(), profile_mean.max()))
    ax_zero.set_xlim(0, xmax)
    ax_mean.set_xlim(0, xmax)
    color_axis = fig.add_axes([0.105, 0.12, 0.41, 0.025])
    colorbar = fig.colorbar(mesh, cax=color_axis, orientation="horizontal")
    colorbar.set_label("Potência (u.a.)", fontsize=8.5, labelpad=1)
    colorbar.ax.tick_params(labelsize=8.5, pad=1)
    colorbar.ax.xaxis.set_major_locator(plt.MaxNLocator(4))
    colorbar.ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}".replace(".", ",")))

    pdf_path = PREFIX.with_suffix(".pdf")
    png_path = PREFIX.with_suffix(".png")
    fig.savefig(pdf_path, dpi=360)
    fig.savefig(png_path, dpi=220)
    plt.close(fig)
    assert digest(CORE_PATH) == core_hash
    metadata = {
        "purpose": "Ilustração didática sintética dos dois perfis extraídos de um mesmo campo de potência. Não representa um ensaio ou padrão experimental.",
        "signal": {
            "samples": n, "sample_spacing_pixels": 1,
            "construction": "127.5 + 30 * (0.7 sin(2πn/760) + 2.2 exp(-0.5((n-4200)/650)^2) sin(2πn/210+0.4) + 0.3 cos(2πn/1450))",
            "randomness": "Nenhuma. Construção determinística.",
            "preparation": "core.prepare: inversão de intensidade, retirada de tendência linear, centralização e divisão pelo desvio padrão populacional.",
        },
        "wavelet": {
            "engine": "morlet", "omega0": core.OMEGA0,
            "fourier_factor": core.FOURIER_FACTOR,
            "wavelengths": {"count": len(wavelengths), "minimum": 50, "maximum": 3000, "spacing": "geomspace"},
            "autocorrelation_estimator": "biased, divisão pelo número de amostras do sinal",
            "lag_count": len(lags), "zero_index": zero_index,
            "power_definition": "abs(CWT(autocorrelation)), sem elevar ao quadrado",
            "cwt_padding_sigma": 8,
        },
        "profiles": {
            "zero": "power[:, zero_index]",
            "mean": "power.mean(axis=1), incluindo todos os atrasos positivos, negativos e zero",
            "zero_max_error_against_complex_modulus": float(np.max(np.abs(profile_zero-np.abs(coefficients[:,zero_index])))),
            "mean_max_error_against_explicit_sum": float(np.max(np.abs(profile_mean-np.sum(power,axis=1)/len(lags)))),
            "normalization_for_display": "Nenhuma. Os dois perfis e o mapa conservam a potência calculada, em unidades arbitrárias.",
        },
        "figure": {
            "width_cm": width_cm, "height_cm": height_cm,
            "minimum_font_points": 8.5, "pdf_raster_dpi": 360,
            "preview_dpi": 220, "shared_y_axis": "logarítmico, 50 a 3000 pixels",
            "profile_x_axes": "Mesmos limites em unidades arbitrárias",
            "drawing_lag_indices": drawing_indices.tolist(),
            "drawing_only_reduction": "A cada 8 atrasos, incluindo zero e o último. Perfis calculados do campo completo.",
        },
        "versions": {"numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
        "source_core": str(CORE_PATH.relative_to(ROOT)), "source_core_sha256": core_hash,
        "script_sha256": digest(Path(__file__)),
        "files_sha256": {path.name: digest(path) for path in (pdf_path,png_path,data_path)},
        "status": "PASS",
    }
    (HERE / "perfis_wavelet_origem.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps({"status": "PASS", "outputs": [p.name for p in (pdf_path,png_path,data_path)],
                      "power_shape": list(power.shape), "profile_zero_peak": float(profile_zero.max()),
                      "profile_mean_peak": float(profile_mean.max())}, ensure_ascii=False))


if __name__ == "__main__":
    main()
