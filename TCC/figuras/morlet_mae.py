"""Reproduce the Morlet mother-wavelet illustration used in the thesis.

Run from any working directory with ``python3 TCC/figuras/morlet_mae.py``.
The figure evaluates eq:morlet, without a zero-mean correction term.
It contains no experimental measurements and does not change the analysis runs.
"""

from pathlib import Path
import hashlib
import json
import os
import re
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "tcc-morlet-mpl"))

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OMEGA0 = 6.0


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    # A unit-scale, zero-position mother wavelet in dimensionless coordinates.
    t = np.linspace(-4.0, 4.0, 4001)
    envelope = np.pi ** (-0.25) * np.exp(-0.5 * t**2)
    psi = envelope * np.exp(1j * OMEGA0 * t)
    np.testing.assert_allclose(np.abs(psi), envelope, rtol=3e-16, atol=2e-16)
    np.testing.assert_allclose(psi.real, envelope * np.cos(OMEGA0 * t), atol=2e-16)
    np.testing.assert_allclose(psi.imag, envelope * np.sin(OMEGA0 * t), atol=2e-16)

    # Check the normalization over a wider integration domain than the plot.
    check_t = np.linspace(-8.0, 8.0, 32001)
    energy = float(np.trapezoid(np.pi ** (-0.5) * np.exp(-check_t**2), check_t))
    assert abs(energy - 1.0) < 1e-12

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.7,
        "pdf.fonttype": 42,
        "savefig.facecolor": "white",
    })
    fig, ax = plt.subplots(figsize=(7.0, 3.25), layout="constrained")
    ax.axhline(0, color="#87929a", linewidth=0.65, zorder=0)
    ax.plot(t, psi.real, color="#087f8c", linewidth=1.8, label=r"Parte real de $\psi(t)$")
    ax.plot(t, psi.imag, color="#b55b27", linewidth=1.55,
            linestyle=(0, (5, 2)), label=r"Parte imaginária de $\psi(t)$")
    ax.set(xlim=(-4, 4), ylim=(-0.83, 0.85),
           xlabel=r"Coordenada $t$ (adimensional)", ylabel="Amplitude (adimensional)")
    ax.set_xticks(np.arange(-4, 5))
    ax.set_yticks([-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}".replace(".", ",")))
    ax.grid(axis="x", color="#d5dce0", linewidth=0.45, alpha=0.7)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.20), ncol=2,
              frameon=False, fontsize=9, handlelength=2.7, columnspacing=1.3)

    pdf = HERE / "morlet_mae.pdf"
    png = HERE / "morlet_mae.png"
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.05,
                metadata={"Title": "Wavelet mãe de Morlet com omega0 = 6",
                          "Author": "TCC", "CreationDate": None, "ModDate": None})
    fig.savefig(png, dpi=240, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)

    core = ROOT / "src2/v16_spatial2_18_71/scripts/wavelet_core.py"
    tex = ROOT / "TCC/fundamentacao.tex"
    equation = re.search(r"\\begin\{equation\}\\label\{eq:morlet\}.*?\\end\{equation\}",
                         tex.read_text(), flags=re.S).group(0)
    manifest = {
        "kind": "Analytical illustration, not experimental data",
        "equation_label": "eq:morlet",
        "formula": "psi(t) = pi^(-1/4) * exp(-t^2/2) * exp(1j*6*t)",
        "omega0": OMEGA0,
        "zero_mean_correction": False,
        "scale": 1.0,
        "position": 0.0,
        "coordinate": "dimensionless t",
        "plot_domain": [-4.0, 4.0],
        "samples": len(t),
        "plotted_components": ["real part of psi(t)", "imaginary part of psi(t)"],
        "plotted_arrays_sha256_float64_little_endian": {
            "coordinate": hashlib.sha256(t.astype("<f8").tobytes()).hexdigest(),
            "real": hashlib.sha256(psi.real.astype("<f8").tobytes()).hexdigest(),
            "imaginary": hashlib.sha256(psi.imag.astype("<f8").tobytes()).hexdigest(),
        },
        "energy_integral_on_minus8_plus8": energy,
        "source_equation_sha256": hashlib.sha256(equation.encode()).hexdigest(),
        "analysis_implementation": str(core.relative_to(ROOT)),
        "analysis_implementation_sha256": sha256(core),
        "script": str(Path(__file__).resolve().relative_to(ROOT)),
        "script_sha256": sha256(Path(__file__)),
        "outputs": {str(p.relative_to(ROOT)): sha256(p) for p in (pdf, png)},
    }
    (HERE / "morlet_mae_origem.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"pdf": str(pdf), "png": str(png), "energy": energy}))


if __name__ == "__main__":
    main()
