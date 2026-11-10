#!/usr/bin/env python3
"""
v13 — Mapa fração de vazio (α) × velocidade superficial do líquido (Jl)
=======================================================================

Combina as duas análises anteriores:
  - α  (fração de vazio cordal)       : reutiliza void_fraction() de v11
  - Jl (velocidade superficial líq.)  : reutiliza compute_runs() de v12

Produz o diagrama α × Jl análogo ao gráfico de referência da literatura
(Baker-style, mas no plano void × Jl em vez de Jg × Jl), colorido por
padrão de escoamento.

Saídas em src2/v13/:
  void_vs_jl_v13.png    — figura principal (scatter + clusters)
  void_vs_jl_table_v13.md — tabela por ensaio + resumo por padrão

Uso:
  cd src2/v13 && python3 void_vs_jl_v13.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

# Adiciona caminhos de v11 e v12 ao sys.path para importação direta.
HERE = Path(__file__).resolve().parent
SRC2 = HERE.parent
sys.path.insert(0, str(SRC2 / "v11"))
sys.path.insert(0, str(SRC2 / "v12"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.ticker import MultipleLocator

try:
    from scipy.spatial import ConvexHull
    _HAVE_SCIPY = True
except ImportError:
    _HAVE_SCIPY = False

# Importa funções de v11 e v12 (sem re-implementar).
from void_fraction_v11 import (
    list_scenarios, folder_flow_type, void_fraction, parse_m_q,
    GROUP_MAP, GROUP_COLOR, GROUP_ORDER, SAT_LIMIT, PANORAMA, DATA_ROOT,
)
from velocity_map_v12 import (
    fluid_properties, compute_runs,
    PATTERN_ORDER, PATTERN_PT, PALETTE, MARKERS,
)
from PIL import Image

# Mapeamento: rótulo fino (v11) → chave canônica (v12)
_V11_TO_V12 = {
    "plug":         "plug",
    "slug":         "slug",
    "churn":        "churn",
    "wavy annular": "wavy annular",
}


# ---------------------------------------------------------------------------
# Coleta α por ensaio (v11)
# ---------------------------------------------------------------------------
def collect_alpha() -> dict[str, dict]:
    """Retorna {run_key: {alpha, alpha_std, group, m, q, dropped}}."""
    result = {}
    for name in list_scenarios():
        folder = os.path.join(DATA_ROOT, name)
        pano = np.array(Image.open(os.path.join(folder, PANORAMA)).convert("L"))
        alpha, meta = void_fraction(pano)
        label = folder_flow_type(name)
        group = GROUP_MAP.get(label, label)
        m, q = parse_m_q(name)
        dropped = meta["sat_frac"] > SAT_LIMIT
        # chave sem _001/_002 para cruzar com v12, mais o sufixo exato
        result[name] = dict(
            alpha=alpha, alpha_std=meta["alpha_std"],
            group=group, m=m, q=q, dropped=dropped,
        )
    return result


# ---------------------------------------------------------------------------
# Cruzamento α (v11) × Jl (v12) por ensaio
# ---------------------------------------------------------------------------
def build_combined(alpha_map: dict, v12_rows: list) -> list[dict]:
    """
    Cruza por (m, q, rep): ensaio presente em AMBOS os datasets.

    v11 usa 26 ensaios de padrão puro; v12 usa os mesmos 28 (mais bolha
    dispersa e smooth annular mapeados como wavy annular).
    Só entram os ensaios que existem em v11 (têm panorama + rótulo puro).
    """
    # Índice de v12: (m, q, rep) -> row
    v12_idx: dict[tuple, dict] = {}
    for r in v12_rows:
        if r["pattern"] is None:
            continue
        rep = r["run"].split("_")[1] if "_" in r["run"] else "001"
        v12_idx[(r["m"], r["q"], rep)] = r

    combined = []
    for name, av in alpha_map.items():
        rep = name.split("_")[1] if "_" in name else "001"
        key = (av["m"], av["q"], rep)
        if key not in v12_idx:
            continue
        v12r = v12_idx[key]
        combined.append({
            "run": name,
            "m": av["m"], "q": av["q"],
            "G": v12r["G"],
            "alpha": av["alpha"],
            "alpha_std": av["alpha_std"],
            "Jl": v12r["Jl"],
            "Jg": v12r["Jg"],
            "x": v12r["x"],
            "group": av["group"],          # rótulo v11 (plug/slug/churn/annular)
            "pattern": v12r["pattern"],    # rótulo v12 canônico
            "dropped": av["dropped"],
        })
    return combined


# ---------------------------------------------------------------------------
# Figura: α × Jl
# ---------------------------------------------------------------------------
def _draw_cluster(ax, pts, color):
    P = np.asarray(pts, float)
    n = len(P)
    degenerate = n < 3
    if not degenerate and _HAVE_SCIPY:
        try:
            area = ConvexHull(P).volume
            degenerate = area < 1e-9
        except Exception:
            degenerate = True
    if not degenerate and _HAVE_SCIPY:
        hull = ConvexHull(P)
        poly = P[hull.vertices]
        ax.add_patch(Polygon(poly, closed=True, facecolor=color, alpha=0.14,
                             edgecolor=color, linewidth=1.3, zorder=1))
    else:
        mu = P.mean(axis=0)
        std = P.std(axis=0) if n >= 2 else np.array([0.02, 0.01])
        half = np.maximum(std * 1.5, [0.02, 0.01])
        from matplotlib.patches import Ellipse
        ax.add_patch(Ellipse(mu, width=2 * half[0], height=2 * half[1],
                             facecolor=color, alpha=0.12, edgecolor=color,
                             linewidth=1.1, linestyle="--", zorder=1))


def plot_void_vs_jl(rows: list[dict], outpath: Path) -> None:
    kept = [r for r in rows if not r["dropped"]]

    fig, ax = plt.subplots(figsize=(8.5, 6.5), dpi=140)

    # grupos na ordem de v12 (bubble, plug, slug, churn, wavy annular)
    # mas filtra só os que v11 tem (sem bubble aqui)
    order = ["plug", "slug", "churn", "wavy annular"]
    for pat in order:
        pr = [r for r in kept if r["pattern"] == pat]
        if not pr:
            continue
        xs = np.array([r["Jl"] for r in pr])
        ys = np.array([r["alpha"] for r in pr])
        c = PALETTE[pat]
        _draw_cluster(ax, np.column_stack([xs, ys]), c)
        ax.scatter(xs, ys, s=75, marker=MARKERS[pat], facecolor=c,
                   edgecolor="white", linewidth=1.0, zorder=3,
                   label=PATTERN_PT[pat])
        # centróide
        ax.scatter([xs.mean()], [ys.mean()], s=170, marker=MARKERS[pat],
                   facecolor="none", edgecolor=c, linewidth=2.0, zorder=4)

    # saturado (descartado) — marcador vazado cinza
    sat = [r for r in rows if r["dropped"]]
    if sat:
        ax.scatter([r["Jl"] for r in sat], [r["alpha"] for r in sat],
                   s=80, facecolors="none", edgecolors="0.5", linewidth=1.3,
                   marker="o", zorder=3, label="saturado (fora da estatística)")

    ax.set_xlabel(r"$J_l$ — velocidade superficial de líquido  [m/s]", fontsize=11)
    ax.set_ylabel(r"Fração de vazio cordal  $\alpha$  [–]", fontsize=11)
    ax.set_title("Fração de vazio × velocidade superficial de líquido\n"
                 "R600a, Tsat 27 °C, D 1 mm  |  α cordal (v11) × Jl modelo equilíbrio (v12)",
                 fontsize=11)

    ax.xaxis.set_major_locator(MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.05))
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.set_xlim(left=0.0)
    ax.set_ylim(bottom=0.0, top=1.0)

    ax.grid(True, which="major", color="#dddddd", linewidth=0.6, zorder=0)
    ax.grid(True, which="minor", color="#eeeeee", linewidth=0.3, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#888888")

    leg = ax.legend(title="Padrão", frameon=True, framealpha=0.95,
                    loc="upper right", fontsize=10)
    leg.get_frame().set_edgecolor("#cccccc")
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)
    print(f"figura: {outpath}")


# ---------------------------------------------------------------------------
# Tabela markdown
# ---------------------------------------------------------------------------
def write_table(rows: list[dict], outpath: Path) -> None:
    kept = [r for r in rows if not r["dropped"]]
    lines = []
    ap = lines.append

    ap("# v13 — Fração de vazio × velocidade superficial de líquido\n")
    ap("Combina α cordal (v11, Otsu sobre panorama) com Jl modelo de equilíbrio (v12).\n")

    ap("## Por ensaio")
    ap("| Ensaio | Padrão | G [kg/m²s] | q″ [kW/m²] | α [-] | Jl [m/s] | Jg [m/s] | sat |")
    ap("|--------|--------|----------:|-----------:|------:|---------:|---------:|-----|")
    for r in sorted(kept, key=lambda r: (PATTERN_ORDER.index(r["pattern"])
                                         if r["pattern"] in PATTERN_ORDER else 99,
                                         r["G"], r["q"])):
        sat_flag = "✗" if r["dropped"] else ""
        ap(f"| {r['run']} | {PATTERN_PT.get(r['pattern'], r['pattern'])} | {r['G']} | {r['q']} | "
           f"{r['alpha']:.3f} | {r['Jl']:.4f} | {r['Jg']:.3f} | {sat_flag} |")
    ap("")

    ap("## Resumo por padrão (α e Jl)")
    ap("| Padrão | n | α mín | α máx | α média | Jl mín | Jl máx | Jl média |")
    ap("|--------|--:|------:|------:|--------:|-------:|-------:|---------:|")
    for pat in PATTERN_ORDER:
        pr = [r for r in kept if r["pattern"] == pat]
        if not pr:
            continue
        alphas = np.array([r["alpha"] for r in pr])
        jls = np.array([r["Jl"] for r in pr])
        ap(f"| {PATTERN_PT[pat]} | {len(pr)} | "
           f"{alphas.min():.3f} | {alphas.max():.3f} | {alphas.mean():.3f} | "
           f"{jls.min():.4f} | {jls.max():.4f} | {jls.mean():.4f} |")
    ap("")

    ap("## Notas")
    ap("- **α cordal**: fração de vazio calculada por Otsu no panorama espaço-tempo "
       "(v11); reflete uma corda vertical central, não a seção transversal. "
       "Viés cordal > areal para padrões com núcleo gasoso central.")
    ap("- **Jl modelo**: calculado via balanço de energia (modelo homogêneo, "
       "x_in=0, Tsat=27 °C, R600a) — estimativa de equilíbrio, não medido.")
    ap("- m03q40_001 excluído da estatística (saturação 34% do núcleo → α artificialmente baixo).")

    outpath.write_text("\n".join(lines), encoding="utf-8")
    print(f"tabela: {outpath}")


# ---------------------------------------------------------------------------
def main():
    props = fluid_properties()
    v12_rows = compute_runs(props)

    print("Coletando α de v11 ...")
    alpha_map = collect_alpha()

    combined = build_combined(alpha_map, v12_rows)
    n_kept = sum(1 for r in combined if not r["dropped"])
    print(f"{len(combined)} ensaios cruzados ({n_kept} na estatística)")

    out = HERE
    plot_void_vs_jl(combined, out / "void_vs_jl_v13.png")
    write_table(combined, out / "void_vs_jl_table_v13.md")


if __name__ == "__main__":
    main()
