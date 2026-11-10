#!/usr/bin/env python3
"""
v12 — Mapa de velocidades superficiais (Jg × Jl) por padrão de escoamento
===========================================================================

Objetivo
--------
Posicionar cada ensaio (mXXqYY) no plano velocidade superficial de gás (Jg) ×
velocidade superficial de líquido (Jl), agrupado por padrão de escoamento, com
clusters por padrão e uma tabela min/máx/média/mediana de Jg e Jl por padrão.

Considera-se plug, slug, churn, bolha dispersa e anular ondulado (wavy annular
+ smooth annular unificados como "wavy annular"). Descartam-se somente os
ensaios de transição com dois padrões simultâneos (plug-slug, slug-churn).

Por que não o v10
-----------------
O v10 tentou obter velocidade a partir do deslocamento óptico das imagens
(px/quadro × FPS), mas NÃO existe calibração px→mm nos dados (ScaleFactor=1.0),
então o v10 ficou preso a uma "celeridade de estrutura" em px/s — não eram Jg/Jl
físicos. O v12 abandona a rota óptica e calcula velocidades FÍSICAS a partir dos
parâmetros do ensaio (G, q″) via balanço de energia — sem precisar de calibração.

Física
------
Para cada ensaio, com A = πD²/4 e área aquecida = πD·L:

    x_out = 4·q″·L / (G·D·i_lv)          # título de vapor (equilíbrio; entrada
                                          # líquido saturado) — idêntico a Bo·4L/D
    x     = clip(x_out, 0, 1)
    Jg    = G·x       / ρ_v               # velocidade superficial de gás    [m/s]
    Jl    = G·(1 − x) / ρ_l               # velocidade superficial de líquido [m/s]

Propriedades do R600a (isobutano) a Tsat = 27 °C obtidas do CoolProp.

Limitações (declaradas honestamente — ver bloco no .md gerado)
--------------------------------------------------------------
1. x_in = 0 (entrada líquido saturado): ignora o sub-resfriamento imposto pelo
   pré-aquecedor. Logo x é uma estimativa de EQUILÍBRIO (limite superior), não
   medida. Não há temperatura de entrada nos dados brutos para corrigir.
2. G se cancela em Jg: como x ∝ q″/G, Jg = (4L/(D·i_lv·ρ_v))·q″ depende SÓ de q″.
   No plano Jg×Jl os 5 ensaios de mesmo q″ ficam numa linha vertical; a informação
   de G está em Jl. É propriedade conhecida do modelo de equilíbrio, não um erro.
3. x = Bo·4L/D: o título é um reescalonamento do número de ebulição do próprio
   trabalho-base; o escorregamento é ignorado (modelo homogêneo). Jg/Jl são, assim,
   ilustrativos do modelo.

Saídas (em src2/v12/)
---------------------
- velocity_map_v12.png     : PRIMÁRIO — log-log Jg×Jl, clusters por padrão
- velocity_table_v12.md    : tabela por ensaio + resumo por padrão + limitações
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Polygon, Ellipse
    from matplotlib.ticker import (FixedLocator, FixedFormatter, LogLocator,
                                   MultipleLocator, NullFormatter)
except Exception as e:  # pragma: no cover
    sys.exit(f"matplotlib é necessário: {e}")

try:
    from scipy.spatial import ConvexHull
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False

try:
    import CoolProp
    from CoolProp.CoolProp import PropsSI
except Exception as e:  # pragma: no cover
    sys.exit(
        "CoolProp é necessário para as propriedades do R600a.\n"
        "Instale com:  pip3 install CoolProp\n"
        f"(erro: {e})"
    )

# --------------------------------------------------------------------------- #
# Constantes físicas / geometria (fundamentacao.tex, metodologia.tex)
# --------------------------------------------------------------------------- #
FLUID = "R600a"            # isobutano
T_SAT = 300.15             # 27 °C em K
D = 1.0e-3                 # diâmetro interno da TS [m]
L_HEATED = 0.265           # comprimento aquecido da TS [m]
A_CROSS = np.pi * D**2 / 4.0           # área da seção [m²]
A_HEATED = np.pi * D * L_HEATED        # área aquecida (perímetro × L) [m²]

# m → fluxo mássico G [kg/(m²·s)]  (fundamentacao.tex:291; v10/void_velocity_v10.py:77)
G_MAP = {3: 380, 4: 510, 5: 640, 6: 760, 7: 890}

# --------------------------------------------------------------------------- #
# Rótulos de padrão — fonte de verdade: CLASSIFICACAO_PADRAOES_MANUAL.md
# Inclui: plug, slug, churn, bolha dispersa e anular (wavy+smooth → "wavy annular").
# Descartam-se apenas transições com dois padrões (plug-slug, slug-churn).
# --------------------------------------------------------------------------- #
RAW_LABELS = {
    "m03q01": "plug-slug",       "m03q05": "churn",
    "m03q10": "wavy annular",    "m03q20": "wavy annular",
    "m03q40": "wavy annular",    "m03q80": "smooth annular",
    "m04q01": "slug",            "m04q05": "slug-churn",
    "m04q10": "churn",           "m04q20": "churn",
    "m04q40": "wavy annular",    "m04q80": "wavy annular",
    "m05q01": "plug",            "m05q05": "slug",
    "m05q10": "churn",           "m05q20": "churn",
    "m05q40": "churn",           "m05q80": "wavy annular",
    "m06q01": "plug",            "m06q05": "slug",
    "m06q10": "churn",           "m06q20": "churn",
    "m06q40": "wavy annular",    "m06q80": "wavy annular",
    "m07q01": "dispersed bubbles", "m07q05": "plug",
    "m07q10": "slug-churn",      "m07q20": "churn",
    "m07q40": "churn",           "m07q80": "wavy annular",
}

# Normalização para os 5 padrões-alvo; None = descartar (somente transições).
def canonical_pattern(label: str):
    if label in ("plug", "slug", "churn"):
        return label
    if "annular" in label:          # wavy annular + smooth annular → wavy annular
        return "wavy annular"
    if label == "dispersed bubbles":
        return "bubble"
    return None                     # plug-slug, slug-churn

PATTERN_ORDER = ["bubble", "plug", "slug", "churn", "wavy annular"]
PATTERN_PT = {
    "bubble":       "bolha dispersa",
    "plug":         "plug",
    "slug":         "slug",
    "churn":        "churn",
    "wavy annular": "anular ondulado",
}

# Paleta categórica validada (dataviz validate_palette: ALL CHECKS PASS, light).
# Encoding secundário por marcador (CVD/impressão): nunca cor sozinha.
PALETTE = {
    "bubble":       "#e6ab02",   # âmbar
    "plug":         "#4269d0",   # azul
    "slug":         "#d55e00",   # vermelhão
    "churn":        "#009e73",   # verde
    "wavy annular": "#9c6ade",   # roxo
}
MARKERS = {"bubble": "P", "plug": "o", "slug": "s", "churn": "^", "wavy annular": "D"}

RUN_RE = re.compile(r"^m(\d{2})q(\d{2})(?:_(\d{3}))?$")


def parse_run(name: str):
    """'m05q10' ou 'm05q05_002' -> (m, q, rep).  (parse_m_q, v10:294, estendido)"""
    m = RUN_RE.match(name)
    if not m:
        return None
    return int(m.group(1)), int(m.group(2)), (m.group(3) or "001")


# --------------------------------------------------------------------------- #
# Propriedades do fluido
# --------------------------------------------------------------------------- #
def fluid_properties():
    rho_l = PropsSI("D", "T", T_SAT, "Q", 0, FLUID)
    rho_v = PropsSI("D", "T", T_SAT, "Q", 1, FLUID)
    h_l = PropsSI("H", "T", T_SAT, "Q", 0, FLUID)
    h_v = PropsSI("H", "T", T_SAT, "Q", 1, FLUID)
    i_lv = h_v - h_l
    p_sat = PropsSI("P", "T", T_SAT, "Q", 0, FLUID)
    return dict(rho_l=rho_l, rho_v=rho_v, i_lv=i_lv, p_sat=p_sat,
                coolprop_version=CoolProp.__version__)


# --------------------------------------------------------------------------- #
# Cálculo por ensaio
# --------------------------------------------------------------------------- #
def compute_runs(props):
    rho_l, rho_v, i_lv = props["rho_l"], props["rho_v"], props["i_lv"]
    rows = []
    # Enumera todas as condições mXXqYY (001) + a repetição m05q05_002.
    conditions = list(RAW_LABELS.keys()) + ["m05q05_002"]
    for cond in conditions:
        base = cond.split("_")[0]
        parsed = parse_run(cond)
        if parsed is None:
            continue
        m, q, rep = parsed
        label = RAW_LABELS[base]
        pat = canonical_pattern(label)
        if pat is None:
            rows.append(dict(run=cond, m=m, q=q, G=G_MAP[m], label=label,
                             pattern=None))           # descartado
            continue
        # m05q05_002: marcador __plug.txt conflita com a classificação (slug).
        # Resolve-se pela classificação (fonte de verdade) → slug.
        G = G_MAP[m]
        q_wm2 = q * 1e3                               # kW/m² → W/m²
        x_out = 4.0 * q_wm2 * L_HEATED / (G * D * i_lv)   # = Bo·4L/D
        x = float(np.clip(x_out, 0.0, 1.0))
        Jg = G * x / rho_v
        Jl = G * (1.0 - x) / rho_l
        mdot = G * A_CROSS                            # vazão mássica [kg/s]
        rows.append(dict(run=cond, m=m, q=q, G=G, label=label, pattern=pat,
                         x_out=x_out, x=x, Jg=Jg, Jl=Jl, mdot=mdot))
    return rows


# --------------------------------------------------------------------------- #
# Verificações internas (falha ruidosa se violadas)
# --------------------------------------------------------------------------- #
def sanity_checks(rows, props):
    kept = [r for r in rows if r["pattern"] is not None]
    msgs = []

    # Área e vazão mássica
    assert abs(A_CROSS - 7.854e-7) < 1e-9, f"A_CROSS={A_CROSS}"
    mdot_m05 = G_MAP[5] * A_CROSS * 1e3               # g/s
    assert 0.30 <= mdot_m05 <= 0.70, f"mdot(m05)={mdot_m05:.3f} g/s fora de 0.3-0.7"
    msgs.append(f"A={A_CROSS:.3e} m²; mdot(m05)={mdot_m05:.3f} g/s  [ok]")

    # x_out extremos (bate com a verificação adversarial)
    xd = {r["run"]: r for r in rows}
    x_hi = 4 * 80e3 * L_HEATED / (G_MAP[3] * D * props["i_lv"])
    x_lo = 4 * 1e3 * L_HEATED / (G_MAP[7] * D * props["i_lv"])
    assert 0.65 <= x_hi <= 0.72, f"x_out(m03q80)={x_hi:.3f}"
    assert x_lo < 0.01, f"x_out(m07q01)={x_lo:.4f}"
    msgs.append(f"x_out: m03q80={x_hi:.3f} (máx), m07q01={x_lo:.4f} (mín)  [ok]")

    # Monotonicidade: x sobe com q (G fixo) e desce com G (q fixo)
    for m in G_MAP:
        xs = [(r["q"], r["x_out"]) for r in kept if r["m"] == m]
        xs.sort()
        vals = [v for _, v in xs]
        assert all(a <= b + 1e-12 for a, b in zip(vals, vals[1:])), \
            f"x não-monotônico em q para m={m}: {xs}"
    msgs.append("x monotônico ↑ com q″ (G fixo)  [ok]")

    # G cancela em Jg: mesmo q″ → mesmo Jg (dentro do clip)
    byq = {}
    for r in kept:
        if r["x"] < 1.0:                              # fora do clip
            byq.setdefault(r["q"], []).append(r["Jg"])
    for q, jgs in byq.items():
        if len(jgs) > 1:
            assert max(jgs) - min(jgs) < 1e-6 * max(jgs), \
                f"Jg varia em q={q}: {jgs}"
    msgs.append("Jg idêntico entre runs de mesmo q″ (cancelamento de G)  [ok]")

    # Contagem
    n_by = {p: sum(1 for r in kept if r["pattern"] == p) for p in PATTERN_ORDER}
    assert n_by == {"bubble": 1, "plug": 3, "slug": 4, "churn": 10, "wavy annular": 10}, n_by
    assert len(kept) == 28, f"n mantidos={len(kept)}"
    msgs.append(f"ensaios: {n_by} (total {len(kept)}); descartados "
                f"{len(rows)-len(kept)}  [ok]")
    return msgs


# --------------------------------------------------------------------------- #
# Clusters: casca convexa com fallback (elipse 1.5σ / segmento)
# --------------------------------------------------------------------------- #
def _draw_cluster(ax, pts, color, logx=False, logy=False):
    """Desenha a região do cluster em pts (array Nx2, em coords de dados)."""
    P = np.asarray(pts, float)
    # trabalha no espaço em que os eixos vivem (log se aplicável) p/ hull/elipse
    Q = P.copy()
    if logx:
        Q[:, 0] = np.log10(Q[:, 0])
    if logy:
        Q[:, 1] = np.log10(Q[:, 1])

    def to_data(xy):
        xy = np.asarray(xy, float).copy()
        if logx:
            xy[..., 0] = 10 ** xy[..., 0]
        if logy:
            xy[..., 1] = 10 ** xy[..., 1]
        return xy

    n = len(Q)
    # degenerado: <3 pontos OU todos colineares → elipse 1.5σ (ou segmento)
    degenerate = n < 3
    if not degenerate and _HAVE_SCIPY:
        area = 0.0
        try:
            area = ConvexHull(Q).volume   # 2D → área
        except Exception:
            degenerate = True
        if area < 1e-9:
            degenerate = True

    if not degenerate and _HAVE_SCIPY:
        hull = ConvexHull(Q)
        poly = to_data(Q[hull.vertices])
        ax.add_patch(Polygon(poly, closed=True, facecolor=color, alpha=0.14,
                             edgecolor=color, linewidth=1.3, zorder=1))
    else:
        # fallback: elipse de covariância 1.5σ no espaço (log-)transformado
        mu = Q.mean(axis=0)
        if n >= 2:
            cov = np.cov(Q.T)
            cov = np.atleast_2d(cov)
            if cov.shape == (2, 2) and np.all(np.isfinite(cov)):
                vals, vecs = np.linalg.eigh(cov)
                vals = np.clip(vals, 1e-6, None)
                order = vals.argsort()[::-1]
                vals, vecs = vals[order], vecs[:, order]
                ang = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
                w, h = 2 * 1.5 * np.sqrt(vals)
                # elipse precisa de centro em coords de dados; aproxima via patch
                # desenhado no espaço transformado não é trivial → usa transform
                cx, cy = to_data(mu)
                # Para log, uma elipse "visual" é desenhada em torno do centro
                # usando a transformação de eixos; aqui aproximamos por um
                # polígono amostrado da elipse no espaço transformado.
                t = np.linspace(0, 2 * np.pi, 60)
                ell = np.column_stack([w / 2 * np.cos(t), h / 2 * np.sin(t)])
                R = np.array([[np.cos(np.radians(ang)), -np.sin(np.radians(ang))],
                              [np.sin(np.radians(ang)),  np.cos(np.radians(ang))]])
                ell = (ell @ R.T) + mu
                ax.add_patch(Polygon(to_data(ell), closed=True, facecolor=color,
                                     alpha=0.12, edgecolor=color, linewidth=1.1,
                                     linestyle="--", zorder=1))


def _scatter_by_pattern(ax, rows, xkey, ykey, logx, logy):
    for pat in PATTERN_ORDER:
        pr = [r for r in rows if r["pattern"] == pat]
        if not pr:
            continue
        xs = np.array([r[xkey] for r in pr], float)
        ys = np.array([r[ykey] for r in pr], float)
        c = PALETTE[pat]
        _draw_cluster(ax, np.column_stack([xs, ys]), c, logx=logx, logy=logy)
        ax.scatter(xs, ys, s=70, marker=MARKERS[pat], facecolor=c,
                   edgecolor="white", linewidth=1.1, zorder=3,
                   label=PATTERN_PT[pat])
        # centroide (média geométrica em log)
        cx = 10 ** np.mean(np.log10(xs)) if logx else xs.mean()
        cy = 10 ** np.mean(np.log10(ys)) if logy else ys.mean()
        ax.scatter([cx], [cy], s=160, marker=MARKERS[pat], facecolor="none",
                   edgecolor=c, linewidth=2.0, zorder=4)


def _style_axes(ax):
    ax.grid(True, which="both", color="#dddddd", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#888888")


def plot_velocity_map(rows, props, outpath):
    kept = [r for r in rows if r["pattern"] is not None]
    # piso em Jg para segurança no eixo log (runs q01 ~ monofásico)
    jg_pos = [r["Jg"] for r in kept if r["Jg"] > 0]
    floor = min(jg_pos) * 0.5
    for r in kept:
        r["Jg_plot"] = max(r["Jg"], floor)

    fig, ax = plt.subplots(figsize=(8.2, 6.4), dpi=140)
    _scatter_by_pattern(ax, kept, "Jg_plot", "Jl", logx=True, logy=True)
    ax.set_xscale("log")
    ax.set_yscale("log")

    # Eixos log com MAIS marcações: ticks maiores rotulados nos valores físicos
    # relevantes + ticks menores (minor) em toda a década para leitura fina.
    jg_major = [0.3, 0.5, 1, 2, 3, 5, 10, 20, 30]
    jl_major = [0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0]
    ax.xaxis.set_major_locator(FixedLocator(jg_major))
    ax.xaxis.set_major_formatter(FixedFormatter([f"{v:g}" for v in jg_major]))
    ax.yaxis.set_major_locator(FixedLocator(jl_major))
    ax.yaxis.set_major_formatter(FixedFormatter([f"{v:g}" for v in jl_major]))
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1,
                                          numticks=100))
        axis.set_minor_formatter(NullFormatter())
    ax.tick_params(axis="both", which="major", labelsize=9)

    ax.set_xlabel(r"$J_g$ — velocidade superficial de gás  [m/s]", fontsize=11)
    ax.set_ylabel(r"$J_l$ — velocidade superficial de líquido  [m/s]", fontsize=11)
    ax.set_title("Mapa de padrões: velocidade superficial de gás × de líquido\n"
                 "R600a, Tsat 27 °C, D 1 mm — modelo de título de equilíbrio",
                 fontsize=11.5)
    _style_axes(ax)
    leg = ax.legend(title="Padrão", frameon=True, framealpha=0.95,
                    loc="lower left", fontsize=10)
    leg.get_frame().set_edgecolor("#cccccc")
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)


def plot_flux_quality_map(rows, outpath):
    kept = [r for r in rows if r["pattern"] is not None]
    fig, ax = plt.subplots(figsize=(8.2, 6.4), dpi=140)
    _scatter_by_pattern(ax, kept, "x", "G", logx=False, logy=False)

    # Mais marcações: x a cada 0.05 (minor 0.025); G nos cinco valores reais de
    # fluxo mássico, com ticks menores a cada 20.
    ax.xaxis.set_major_locator(MultipleLocator(0.05))
    ax.xaxis.set_minor_locator(MultipleLocator(0.025))
    g_vals = sorted(G_MAP.values())           # 380,510,640,760,890
    ax.yaxis.set_major_locator(FixedLocator(g_vals))
    ax.yaxis.set_major_formatter(FixedFormatter([f"{v:g}" for v in g_vals]))
    ax.yaxis.set_minor_locator(MultipleLocator(20))
    ax.set_xlim(left=-0.01)
    ax.tick_params(axis="both", which="major", labelsize=9)

    ax.set_xlabel(r"Título de vapor de equilíbrio  $x$  [-]", fontsize=11)
    ax.set_ylabel(r"Fluxo mássico  $G$  [kg/(m²·s)]", fontsize=11)
    ax.set_title("Cross-check honesto: fluxo mássico × título de vapor\n"
                 "(G é eixo independente — sem o cancelamento do mapa Jg×Jl)",
                 fontsize=11.5)
    _style_axes(ax)
    leg = ax.legend(title="Padrão", frameon=True, framealpha=0.95,
                    loc="upper right", fontsize=10)
    leg.get_frame().set_edgecolor("#cccccc")
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Tabela / relatório
# --------------------------------------------------------------------------- #
def _stats(vals):
    a = np.asarray(vals, float)
    return a.min(), a.max(), a.mean(), np.median(a)


def write_report(rows, props, msgs, outpath):
    kept = [r for r in rows if r["pattern"] is not None]
    dropped = [r for r in rows if r["pattern"] is None]

    L = []
    ap = L.append
    ap("# v12 — Velocidades superficiais (Jg, Jl) por padrão de escoamento\n")
    ap("Mapa de padrões de escoamento no plano velocidade superficial de gás "
       "(Jg) × velocidade superficial de líquido (Jl), com clusters por padrão. "
       "Padrões incluídos: **bolha dispersa, plug, slug, churn e anular ondulado**.\n")

    ap("## Propriedades do R600a (isobutano) a Tsat = 27 °C")
    ap(f"- Fonte: **CoolProp {props['coolprop_version']}** (`PropsSI`, fluido `R600a`)")
    ap(f"- ρ_l = **{props['rho_l']:.2f} kg/m³**  ·  ρ_v = **{props['rho_v']:.3f} kg/m³**")
    ap(f"- i_lv = **{props['i_lv']/1e3:.2f} kJ/kg**  ·  P_sat = {props['p_sat']/1e5:.3f} bar")
    ap(f"- Geometria: D = 1,0 mm · L_aquecido = 0,265 m · A = {A_CROSS:.3e} m²\n")

    ap("## Modelo")
    ap("```")
    ap("x_out = 4·q″·L / (G·D·i_lv)     (título de equilíbrio; entrada líquido saturado)")
    ap("Jg    = G·x / ρ_v              Jl = G·(1−x) / ρ_l")
    ap("```\n")

    # ---- tabela por ensaio ----
    ap("## Por ensaio")
    ap("| Ensaio | Padrão | G [kg/m²s] | q″ [kW/m²] | x [-] | Jg [m/s] | Jl [m/s] |")
    ap("|--------|--------|-----------:|-----------:|------:|---------:|---------:|")
    for r in sorted(kept, key=lambda r: (PATTERN_ORDER.index(r["pattern"]),
                                          r["G"], r["q"])):
        ap(f"| {r['run']} | {PATTERN_PT[r['pattern']]} | {r['G']} | {r['q']} | "
           f"{r['x']:.4f} | {r['Jg']:.3f} | {r['Jl']:.4f} |")
    ap("")

    # ---- resumo por padrão ----
    ap("## Resumo por padrão (min / máx / média / mediana)")
    ap("| Padrão | n | Jg min | Jg máx | Jg média | Jg mediana | "
       "Jl min | Jl máx | Jl média | Jl mediana |")
    ap("|--------|--:|-------:|-------:|---------:|-----------:|"
       "-------:|-------:|---------:|-----------:|")
    for pat in PATTERN_ORDER:
        pr = [r for r in kept if r["pattern"] == pat]
        jg = _stats([r["Jg"] for r in pr])
        jl = _stats([r["Jl"] for r in pr])
        ap(f"| {PATTERN_PT[pat]} | {len(pr)} | "
           f"{jg[0]:.3f} | {jg[1]:.3f} | {jg[2]:.3f} | {jg[3]:.3f} | "
           f"{jl[0]:.4f} | {jl[1]:.4f} | {jl[2]:.4f} | {jl[3]:.4f} |")
    ap("")

    # ---- descartados ----
    ap("## Ensaios descartados (transições com dois padrões simultâneos)")
    for r in sorted(dropped, key=lambda r: r["run"]):
        ap(f"- `{r['run']}` — {r['label']}")
    ap("")
    ap("> `m05q05_002`: o arquivo-marcador diz `__plug.txt`, mas a classificação "
       "manual (fonte de verdade) rotula m05q05 como **slug**; adotou-se **slug** "
       "(G e q″ idênticos a m05q05_001, logo mesma velocidade).\n")

    # ---- limitações ----
    ap("## Limitações (honestas)")
    ap("1. **x_in = 0 (entrada líquido saturado):** ignora o sub-resfriamento do "
       "pré-aquecedor. x é estimativa de **equilíbrio (limite superior)**, não medida; "
       "não há temperatura de entrada nos dados brutos para corrigir. Consequência: "
       "ensaios anulares saem com x≈0,3–0,68 — abaixo do esperado (alto vazio) para o rótulo.")
    ap("2. **G se cancela em Jg:** como x ∝ q″/G, Jg = (4L/(D·i_lv·ρ_v))·q″ depende "
       "**só de q″**. No plano Jg×Jl os ensaios de mesmo q″ alinham-se numa vertical; "
       "a informação de G está em Jl. Ensaio m07q01 (bolha) tem q″=1 kW/m², logo "
       "Jg=0,336 m/s — mesma posição horizontal que m05q01 e m06q01 (plug), separados "
       "somente por Jl (G diferente).")
    ap("3. **x = Bo·4L/D:** o título é um reescalonamento do número de ebulição do "
       "trabalho-base (de Oliveira et al. 2021); o **escorregamento é ignorado** (modelo "
       "homogêneo). Logo Jg/Jl são ilustrativos do modelo, não velocidades medidas.")
    ap("4. Ensaios q01 (q″=1 kW/m²) dão x≈0,004–0,006 (≈ monofásico); mantidos com "
       "piso em Jg para o eixo log.\n")

    # ---- verificações ----
    ap("## Verificações internas executadas")
    for m in msgs:
        ap(f"- {m}")
    ap("")

    Path(outpath).write_text("\n".join(L), encoding="utf-8")


# --------------------------------------------------------------------------- #
def main():
    here = Path(__file__).resolve().parent
    props = fluid_properties()
    print(f"CoolProp {props['coolprop_version']}: rho_l={props['rho_l']:.2f} "
          f"rho_v={props['rho_v']:.3f} i_lv={props['i_lv']/1e3:.2f} kJ/kg "
          f"Psat={props['p_sat']/1e5:.3f} bar")

    rows = compute_runs(props)
    msgs = sanity_checks(rows, props)
    for m in msgs:
        print("  ✔", m)

    plot_velocity_map(rows, props, here / "velocity_map_v12.png")
    write_report(rows, props, msgs, here / "velocity_table_v12.md")
    print(f"\nEscrito em {here}:")
    print("  - velocity_map_v12.png")
    print("  - velocity_table_v12.md")


if __name__ == "__main__":
    main()
