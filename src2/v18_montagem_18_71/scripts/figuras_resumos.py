#!/usr/bin/env python3
"""Duas distribuições de pares, com caixa interquartil, mediana e média."""
from pathlib import Path
import csv
import hashlib
import json
import os
import sys

sys.dont_write_bytecode = True
os.environ.setdefault("MPLCONFIGDIR", "/tmp/tcc-v18-resumos-mpl")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, MaxNLocator

from runtime import OUTPUT_ROOT, reference

OUT = OUTPUT_ROOT / "06_resumos"
COLORS = ["#27818B", "#CC9142"]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate():
    source = OUT / "dados/pares_pearson_coerencia.csv"
    with source.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 325
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.titlesize": 9.5, "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(16/2.54, 8.2/2.54))
    fig.subplots_adjust(left=.095, right=.99, bottom=.23, top=.89, wspace=.31)
    rng = np.random.default_rng(20261030)
    plotted = []
    for panel, (ax, metric, title) in enumerate(zip(axes,
            ["pearson_r", "wco_media_integral"],
            ["(a) Coeficiente de Pearson", "(b) Coerência média por par"])):
        groups = [np.asarray([float(row[metric]) for row in rows if row["tipo"] == kind])
                  for kind in ("intra", "inter")]
        assert [len(group) for group in groups] == [87, 238]
        for position, (values, color) in enumerate(zip(groups, COLORS), 1):
            jitter = rng.uniform(-.19, .19, len(values))
            ax.scatter(position+jitter, values, s=8, color=color, alpha=.30,
                       edgecolors="none", zorder=2, rasterized=False)
            ax.boxplot([values], positions=[position], widths=.44,
                whis=(0, 100), showfliers=False, patch_artist=True, manage_ticks=False,
                boxprops=dict(facecolor=color+"38", edgecolor=color, linewidth=1),
                medianprops=dict(color="#263642", linewidth=1.5),
                whiskerprops=dict(color=color, linewidth=.8),
                capprops=dict(color=color, linewidth=.8))
            ax.scatter([position], [values.mean()], marker="D", s=30,
                       color="#263642", edgecolors="white", linewidths=.6, zorder=5)
            q25, q75 = np.quantile(values, [.25, .75], method="linear")
            plotted.append(dict(painel=panel+1, metrica=metric, tipo="intra" if position==1 else "inter",
                pares=len(values), minimo=float(values.min()), media=float(values.mean()),
                mediana=float(np.median(values)), maximo=float(values.max()), q25=float(q25), q75=float(q75)))
        ax.set_title(title, loc="left", pad=9)
        ax.set_xlim(.48, 2.52)
        ax.set_xticks([1, 2], ["Mesmo padrão\n87 pares", "Padrões diferentes\n238 pares"])
        digits = 1 if metric == "pearson_r" else 2
        ax.yaxis.set_major_formatter(FuncFormatter(
            lambda value, _, digits=digits: f"{value:.{digits}f}".replace(".", ",")))
        ax.yaxis.set_major_locator(MaxNLocator(5))
        ax.grid(axis="y", color="#D8DDE0", linewidth=.5, zorder=0)
        ax.spines[["top", "right"]].set_visible(False)
        if metric == "pearson_r":
            ax.set_ylim(-1.05, 1.05)
            ax.set_ylabel("Coeficiente de Pearson r")
        else:
            maximum = max(float(values.max()) for values in groups)
            ax.set_ylim(0, np.ceil(maximum*1.07/.05)*.05)
            ax.set_ylabel("Coerência média do campo")
    legend = [Patch(facecolor="#C9D6DA", edgecolor="#526971", label="25º–75º percentis"),
              Line2D([], [], color="#263642", linewidth=1.5, label="Mediana"),
              Line2D([], [], color="#263642", marker="D", linestyle="none", markersize=5, label="Média")]
    fig.legend(handles=legend, loc="lower center", bbox_to_anchor=(.53, .015),
               ncol=3, frameon=False, fontsize=8.5, handlelength=1.4, columnspacing=1.8)
    destination = OUT / "figuras"
    destination.mkdir(parents=True, exist_ok=True)
    for extension in ("pdf", "png"):
        fig.savefig(destination / f"distribuicoes_pearson_coerencia.{extension}", dpi=360,
                    facecolor="white", metadata={"Creator": "v18 figuras_resumos.py"} if extension=="pdf" else None)
    plt.close(fig)
    metadata = dict(tipo="Distribuições descritivas de pares que compartilham 26 aquisições",
        entrada=dict(arquivo=reference(source), sha256=sha(source)),
        script=dict(arquivo=reference(Path(__file__)), sha256=sha(Path(__file__))),
        dimensoes_cm=[16, 8.2], fonte_pt=9, semente_posicoes_graficas=20261030,
        pontos="Um ponto por par, deslocamento horizontal aleatório apenas para legibilidade.",
        caixa="Percentis 25 e 75, mediana na linha, média no losango, hastes do mínimo ao máximo.",
        inferencia=False, valores=plotted,
        arquivos_sha256={f"distribuicoes_pearson_coerencia.{extension}":sha(destination / f"distribuicoes_pearson_coerencia.{extension}")
                         for extension in ("pdf", "png")})
    (destination / "origem.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2)+"\n")


if __name__ == "__main__":
    generate()
