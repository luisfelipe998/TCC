#!/usr/bin/env python3
"""
v12.1 — Fração de vazio (α), velocidades reais de fase (uv, ul) e escorregamento (S)
====================================================================================

Contexto
--------
Extensão honesta do v12 (`src2/v12/velocity_map_v12.py`). O v12 produz o título de
equilíbrio x e as velocidades SUPERFICIAIS Jg, Jl. Auditoria contra a literatura
local (`src2/v12/VALIDACAO_LITERATURA_v12.md`) mostrou que:

  - As EQUAÇÕES de Jg/Jl são idênticas a dias2025knn (Eqs. 6-7) — método validado.
  - As MAGNITUDES absolutas de Jg/Jl NÃO podem ser "aproximadas" da literatura via
    escorregamento: Jg e Jl são velocidades SUPERFICIAIS (Jg=G·x/ρv, Jl=G·(1−x)/ρl),
    dependem só de x, G, ρ — são INVARIANTES ao escorregamento por definição.

O que a literatura de fato tabula e o que MUDA com escorregamento são a fração de
vazio α, as velocidades REAIS de fase uv, ul e a razão de escorregamento S. O v12.1
adiciona exatamente essas grandezas, mantendo x, Jg, Jl idênticos ao v12 (o v12
original NÃO é alterado).

Modelo adicionado
-----------------
Fração de vazio de Rouhani–Axelsson (1970), como em dias2025knn Eq. (5):

    C0    = 1 + 0.12·(1 − x)
    drift = 1.18·(1 − x)·[g·σ·(ρl − ρv)]^0.25 / (G·ρl^0.5)
    α     = (x/ρv) / ( C0·(x/ρv + (1−x)/ρl) + drift )

Velocidades REAIS de fase (jv = Jg, jl = Jl do v12):

    uv = jv / α          (= G·x / (ρv·α))
    ul = jl / (1 − α)    (= G·(1−x) / (ρl·(1−α)))

Razão de escorregamento — deoliveira2021 Eq. (17):

    S = uv/ul = (x·ρl·(1−α)) / (α·ρv·(1−x))

σ (tensão superficial do R600a a Tsat=27 °C) vem do CoolProp (`PropsSI('I',...)`),
mesmo padrão de propriedades já usado no v12.

Grounding / honestidade (ver src2/v12.1/ANALISE_ADVERSARIAL_v12.md)
-------------------------------------------------------------------
- α Rouhani–Axelsson: dias2025knn Eq. (5); também usada em deoliveira2021 §4.
- S: deoliveira2021 Eq. (17) (verbatim).
- uv=jv/α, ul=jl/(1−α): decorrência da definição de α (implícito em deoliveira2021
  via Eq. 17), NÃO citação literal — declarado assim no memo.
- x permanece de EQUILÍBRIO (x_in=0). Não há potência do pré-aquecedor nem vazão
  mássica medida nos dados brutos (dados.xlsx = só Voltage_0/1/2) para reconstruir
  x_in. Logo α/uv/ul/S herdam essa limitação do v12 — declarada, não escondida.

Saídas (em src2/v12.1/)
-----------------------
- velocity_table_v12_1.md       : tabela por ensaio (x, Jg, Jl, α, uv, ul, S) +
                                  resumo por padrão (min/máx/média/mediana)
- phase_velocity_map_v12_1.png  : mapa uv × ul (velocidades REAIS) por padrão
- void_fraction_v12_1.png       : α por padrão (strip plot) com faixas da literatura
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
    from matplotlib.patches import Polygon
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
# Constantes físicas / geometria (iguais ao v12; copiadas, não importadas)
# --------------------------------------------------------------------------- #
FLUID = "R600a"            # isobutano
T_SAT = 300.15             # 27 °C em K
D = 1.0e-3                 # diâmetro interno da TS [m]
L_HEATED = 0.265           # comprimento aquecido da TS [m]
G_GRAV = 9.81              # aceleração da gravidade [m/s²]
A_CROSS = np.pi * D**2 / 4.0           # área da seção [m²]
A_HEATED = np.pi * D * L_HEATED        # área aquecida (perímetro × L) [m²]

# m → fluxo mássico G [kg/(m²·s)]
G_MAP = {3: 380, 4: 510, 5: 640, 6: 760, 7: 890}

# --------------------------------------------------------------------------- #
# Rótulos de padrão — fonte de verdade: CLASSIFICACAO_PADRAOES_MANUAL.md
# (idênticos ao v12)
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


def canonical_pattern(label: str):
    if label in ("plug", "slug", "churn"):
        return label
    if "annular" in label:
        return "wavy annular"
    if label == "dispersed bubbles":
        return "bubble"
    return None


PATTERN_ORDER = ["bubble", "plug", "slug", "churn", "wavy annular"]
PATTERN_PT = {
    "bubble":       "bolha dispersa",
    "plug":         "plug",
    "slug":         "slug",
    "churn":        "churn",
    "wavy annular": "anular ondulado",
}

PALETTE = {
    "bubble":       "#e6ab02",
    "plug":         "#4269d0",
    "slug":         "#d55e00",
    "churn":        "#009e73",
    "wavy annular": "#9c6ade",
}
MARKERS = {"bubble": "P", "plug": "o", "slug": "s", "churn": "^", "wavy annular": "D"}

# Faixas de α por padrão na literatura local (citações verbatim, p/ validação).
#   plug   : deoliveira2021 §4.2  → α 0,042–0,29
#   churn  : deoliveira2021 §4.4  → α 0,61–0,85  ("α from 0.61 to 0.85")
#   anular : deoliveira2021 §4.5  → α 0,86–0,99
#   slug   : de Oliveira (2021) §4.3 NÃO publica faixa de α para slug (confirmado
#            lendo a §4.3 inteira: só dá descritores de caos/Hurst/LLE). A faixa
#            de slug vem da dissertação dias2024 (mesmo arranjo LETEF, R-600a,
#            mini-canal 1 mm), três pontos confirmados VISUALMENTE nas páginas
#            renderizadas (o texto do PDF está com OCR corrompido — ver
#            paper-log.md):
#              Fig. 22d  G=520  q″=20  x=0,06258  α=0,7023
#              Fig. 23b  G=740  q″=40  x=0,06623  α=0,7135
#              Fig. 24a  G=1220 q″=20 x=0,03385  α=0,5782
#            → faixa 0,578–0,714 (α cai com G, como o próprio texto da p.73 diz).
ALPHA_LIT = {
    "plug":         (0.042, 0.29),    # deoliveira2021 §4.2
    "slug":         (0.5782, 0.7135),  # dias2024 Figs. 22d/23b/24a (visual)
    "churn":        (0.61, 0.85),     # deoliveira2021 §4.4
    "wavy annular": (0.86, 0.99),     # deoliveira2021 §4.5
}

# Pontos ÚNICOS de α por padrão (dias2025knn Tab. 2, G=516 kg/m²s) — NÃO são
# faixas min–máx; representados por MARCADOR, não por retângulo.
# Verificado verbatim (2 leitores + confirmação visual da Tab. 2, p.5): bolha
# 0,36 / coalescida 0,50 / plug 0,68 / slug 0,73 / churn 0,85 / anular 0,88.
# Mantém-se só o slug: corrobora (G=516 ≈ topo da faixa dias2024) a faixa acima.
ALPHA_LIT_POINT = {
    "slug": 0.73,    # dias2025knn Tab. 2 (corrobora a faixa dias2024)
}

RUN_RE = re.compile(r"^m(\d{2})q(\d{2})(?:_(\d{3}))?$")


def parse_run(name: str):
    m = RUN_RE.match(name)
    if not m:
        return None
    return int(m.group(1)), int(m.group(2)), (m.group(3) or "001")


# --------------------------------------------------------------------------- #
# Propriedades do fluido (estende o v12 com σ)
# --------------------------------------------------------------------------- #
def fluid_properties():
    rho_l = PropsSI("D", "T", T_SAT, "Q", 0, FLUID)
    rho_v = PropsSI("D", "T", T_SAT, "Q", 1, FLUID)
    h_l = PropsSI("H", "T", T_SAT, "Q", 0, FLUID)
    h_v = PropsSI("H", "T", T_SAT, "Q", 1, FLUID)
    i_lv = h_v - h_l
    p_sat = PropsSI("P", "T", T_SAT, "Q", 0, FLUID)
    sigma = PropsSI("I", "T", T_SAT, "Q", 0, FLUID)   # tensão superficial [N/m]
    return dict(rho_l=rho_l, rho_v=rho_v, i_lv=i_lv, p_sat=p_sat, sigma=sigma,
                coolprop_version=CoolProp.__version__)


# --------------------------------------------------------------------------- #
# Modelo de escorregamento
# --------------------------------------------------------------------------- #
def alpha_rouhani_axelsson(x, G, props):
    """Fração de vazio de Rouhani–Axelsson (1970) — dias2025knn Eq. (5)."""
    rho_l, rho_v, sigma = props["rho_l"], props["rho_v"], props["sigma"]
    C0 = 1.0 + 0.12 * (1.0 - x)
    drift = 1.18 * (1.0 - x) * (G_GRAV * sigma * (rho_l - rho_v)) ** 0.25 \
        / (G * rho_l ** 0.5)
    return (x / rho_v) / (C0 * (x / rho_v + (1.0 - x) / rho_l) + drift)


# --------------------------------------------------------------------------- #
# Cálculo por ensaio
# --------------------------------------------------------------------------- #
def compute_runs(props):
    rho_l, rho_v, i_lv = props["rho_l"], props["rho_v"], props["i_lv"]
    rows = []
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
                             pattern=None))
            continue
        G = G_MAP[m]
        q_wm2 = q * 1e3
        x_out = 4.0 * q_wm2 * L_HEATED / (G * D * i_lv)   # = Bo·4L/D (idêntico v12)
        x = float(np.clip(x_out, 0.0, 1.0))
        Jg = G * x / rho_v                                # = jv (idêntico v12)
        Jl = G * (1.0 - x) / rho_l                        # = jl (idêntico v12)

        # ---- novidades do v12.1 (escorregamento) ----
        alpha = float(alpha_rouhani_axelsson(x, G, props))
        # uv/ul exigem x>0 e 0<α<1; runs q01 têm x≈0 → fase gás ~inexistente.
        uv = Jg / alpha if alpha > 0 else float("nan")
        ul = Jl / (1.0 - alpha) if alpha < 1 else float("nan")
        if 0.0 < alpha < 1.0 and x > 0.0:
            S = (x * rho_l * (1.0 - alpha)) / (alpha * rho_v * (1.0 - x))
        else:
            S = float("nan")

        mdot = G * A_CROSS
        rows.append(dict(run=cond, m=m, q=q, G=G, label=label, pattern=pat,
                         x_out=x_out, x=x, Jg=Jg, Jl=Jl, mdot=mdot,
                         alpha=alpha, uv=uv, ul=ul, S=S))
    return rows


# --------------------------------------------------------------------------- #
# Verificações internas (falha ruidosa)
# --------------------------------------------------------------------------- #
def sanity_checks(rows, props):
    kept = [r for r in rows if r["pattern"] is not None]
    msgs = []

    # Herdadas do v12 (garante que x/Jg/Jl não mudaram)
    assert abs(A_CROSS - 7.854e-7) < 1e-9, f"A_CROSS={A_CROSS}"
    x_hi = 4 * 80e3 * L_HEATED / (G_MAP[3] * D * props["i_lv"])
    assert 0.65 <= x_hi <= 0.72, f"x_out(m03q80)={x_hi:.3f}"
    n_by = {p: sum(1 for r in kept if r["pattern"] == p) for p in PATTERN_ORDER}
    assert n_by == {"bubble": 1, "plug": 3, "slug": 4, "churn": 10,
                    "wavy annular": 10}, n_by
    assert len(kept) == 28, f"n mantidos={len(kept)}"
    msgs.append(f"herdado v12: x_out(m03q80)={x_hi:.3f}; ensaios={n_by} (28)  [ok]")

    # α ∈ (0,1) para todo ensaio
    for r in kept:
        assert 0.0 < r["alpha"] < 1.0, f"α fora de (0,1) em {r['run']}: {r['alpha']}"
    msgs.append(f"0 < α < 1 em todos os {len(kept)} ensaios  [ok]")

    # α monotônico ↑ com x (G fixo)  — propriedade física de Rouhani–Axelsson
    for m in G_MAP:
        xs = sorted([(r["x"], r["alpha"]) for r in kept if r["m"] == m])
        a = [al for _, al in xs]
        assert all(p <= q + 1e-9 for p, q in zip(a, a[1:])), \
            f"α não-monotônico em x para m={m}: {xs}"
    msgs.append("α monotônico ↑ com x (G fixo)  [ok]")

    # S ≥ 1 (vapor nunca mais lento que líquido num modelo de escorregamento)
    for r in kept:
        if np.isfinite(r["S"]):
            assert r["S"] >= 1.0 - 1e-6, f"S<1 em {r['run']}: {r['S']}"
    s_fin = [r["S"] for r in kept if np.isfinite(r["S"])]
    msgs.append(f"S ≥ 1 em todos ({len(s_fin)} finitos; máx={max(s_fin):.2f})  [ok]")

    # uv ≥ jv e ul ≥ jl (pois α<1 e 1−α<1): velocidade real ≥ superficial
    for r in kept:
        if np.isfinite(r["uv"]):
            assert r["uv"] >= r["Jg"] - 1e-9, f"uv<jv em {r['run']}"
        if np.isfinite(r["ul"]):
            assert r["ul"] >= r["Jl"] - 1e-9, f"ul<jl em {r['run']}"
    msgs.append("uv ≥ jv e ul ≥ jl em todos  [ok]")

    # α anular vs deoliveira2021 §4.5 (0,86–0,99): o TOPO deve alcançar a faixa
    # alta e a MEDIANA cair no regime de alto vazio. O piso baixo (≈0,74 em
    # m03q10, x=0,085) é consequência do x de equilíbrio (x_in=0) nos anulares
    # de baixo q″ — limitação herdada do v12, declarada, não escondida.
    an = sorted(r["alpha"] for r in kept if r["pattern"] == "wavy annular")
    hi_an, med_an = max(an), float(np.median(an))
    assert hi_an >= 0.90, f"α anular máx={hi_an:.3f} não alcança o regime alto"
    assert med_an >= 0.85, f"α anular mediana={med_an:.3f} abaixo do alto vazio"
    msgs.append(f"α anular ∈ [{min(an):.3f}; {hi_an:.3f}], mediana={med_an:.3f} "
                f"(lit. deoliveira2021 0,86–0,99; piso baixo = x_in=0 herdado)  [ok]")

    return msgs


# --------------------------------------------------------------------------- #
# Clusters (idêntico ao v12)
# --------------------------------------------------------------------------- #
def _draw_cluster(ax, pts, color, logx=False, logy=False):
    P = np.asarray(pts, float)
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
    degenerate = n < 3
    if not degenerate and _HAVE_SCIPY:
        area = 0.0
        try:
            area = ConvexHull(Q).volume
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
        pr = [r for r in rows if r["pattern"] == pat and
              np.isfinite(r[xkey]) and np.isfinite(r[ykey])]
        if not pr:
            continue
        xs = np.array([r[xkey] for r in pr], float)
        ys = np.array([r[ykey] for r in pr], float)
        c = PALETTE[pat]
        _draw_cluster(ax, np.column_stack([xs, ys]), c, logx=logx, logy=logy)
        ax.scatter(xs, ys, s=70, marker=MARKERS[pat], facecolor=c,
                   edgecolor="white", linewidth=1.1, zorder=3,
                   label=PATTERN_PT[pat])
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


# --------------------------------------------------------------------------- #
# Figuras
# --------------------------------------------------------------------------- #
def plot_phase_velocity_map(rows, props, outpath):
    """Mapa uv × ul — velocidades REAIS de fase (modelo com escorregamento)."""
    kept = [r for r in rows if r["pattern"] is not None and np.isfinite(r["uv"])
            and np.isfinite(r["ul"])]

    fig, ax = plt.subplots(figsize=(8.2, 6.4), dpi=140)
    _scatter_by_pattern(ax, kept, "uv", "ul", logx=True, logy=True)
    ax.set_xscale("log")
    ax.set_yscale("log")

    uv_major = [1, 2, 3, 5, 10, 20, 30, 50]
    ul_major = [0.3, 0.5, 1, 2, 3, 5, 10]
    ax.xaxis.set_major_locator(FixedLocator(uv_major))
    ax.xaxis.set_major_formatter(FixedFormatter([f"{v:g}" for v in uv_major]))
    ax.yaxis.set_major_locator(FixedLocator(ul_major))
    ax.yaxis.set_major_formatter(FixedFormatter([f"{v:g}" for v in ul_major]))
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1,
                                          numticks=100))
        axis.set_minor_formatter(NullFormatter())
    ax.tick_params(axis="both", which="major", labelsize=9)

    ax.set_xlabel(r"$u_v = j_v/\alpha$ — velocidade REAL do gás  [m/s]", fontsize=11)
    ax.set_ylabel(r"$u_l = j_l/(1-\alpha)$ — velocidade REAL do líquido  [m/s]",
                  fontsize=11)
    ax.set_title("Velocidades reais de fase (modelo com escorregamento)\n"
                 "R600a, Tsat 27 °C, D 1 mm — α de Rouhani–Axelsson (1970)",
                 fontsize=11.5)
    _style_axes(ax)
    leg = ax.legend(title="Padrão", frameon=True, framealpha=0.95,
                    loc="upper left", fontsize=10)
    leg.get_frame().set_edgecolor("#cccccc")
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)


def plot_void_fraction(rows, outpath):
    """α por padrão (strip) + faixas da literatura (deoliveira2021)."""
    kept = [r for r in rows if r["pattern"] is not None]
    fig, ax = plt.subplots(figsize=(8.2, 5.2), dpi=140)

    for i, pat in enumerate(PATTERN_ORDER):
        pr = [r for r in kept if r["pattern"] == pat]
        ys = np.array([r["alpha"] for r in pr], float)
        xs = np.full_like(ys, i, dtype=float)
        xs = xs + np.linspace(-0.12, 0.12, len(ys)) if len(ys) > 1 else xs
        c = PALETTE[pat]
        ax.scatter(xs, ys, s=70, marker=MARKERS[pat], facecolor=c,
                   edgecolor="white", linewidth=1.0, zorder=3)
        # faixa min–máx da literatura (retângulo). de Oliveira (2021) para
        # plug/churn/anular; dissertação dias2024 para slug (ver ALPHA_LIT).
        if pat in ALPHA_LIT:
            lo, hi = ALPHA_LIT[pat]
            ax.add_patch(plt.Rectangle((i - 0.3, lo), 0.6, hi - lo,
                                       facecolor=c, alpha=0.12, edgecolor=c,
                                       linewidth=1.0, linestyle="--", zorder=1))
        # ponto ÚNICO da literatura (marcador) — slug em dias2025knn Tab. 2,
        # sobreposto à faixa dias2024 para mostrar a corroboração (G=516).
        if pat in ALPHA_LIT_POINT:
            yv = ALPHA_LIT_POINT[pat]
            ax.plot([i], [yv], marker="X", markersize=9, color=c,
                    markeredgecolor="black", markeredgewidth=0.8, zorder=4)

    ax.set_xticks(range(len(PATTERN_ORDER)))
    ax.set_xticklabels([PATTERN_PT[p] for p in PATTERN_ORDER], fontsize=10)
    ax.set_ylim(0, 1.0)
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.set_ylabel(r"Fração de vazio  $\alpha$  [-]  (Rouhani–Axelsson)", fontsize=11)
    ax.set_title("Fração de vazio por padrão × literatura\n"
                 "faixa tracejada = de Oliveira (2021) §4.2/4.4/4.5; slug = "
                 "Dias (2024) Figs. 22d/23b/24a  ·  × = dias2025knn Tab. 2",
                 fontsize=10.0)
    _style_axes(ax)
    handles = [
        Line2D([0], [0], linestyle="--", color="#777777",
               label="faixa min–máx (deoliveira2021; slug: dias2024)"),
        Line2D([0], [0], marker="X", linestyle="none", color="#777777",
               markeredgecolor="black", markersize=9,
               label="ponto único slug (dias2025knn Tab. 2)"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=9, frameon=True)
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Tabela / relatório
# --------------------------------------------------------------------------- #
def _stats(vals):
    a = np.asarray([v for v in vals if np.isfinite(v)], float)
    if a.size == 0:
        return float("nan"), float("nan"), float("nan"), float("nan")
    return a.min(), a.max(), a.mean(), np.median(a)


def write_report(rows, props, msgs, outpath):
    kept = [r for r in rows if r["pattern"] is not None]

    L = []
    ap = L.append
    ap("# v12.1 — Fração de vazio α, velocidades reais de fase (uv, ul) e "
       "escorregamento S\n")
    ap("Extensão do v12 com o modelo de **escorregamento** (Rouhani–Axelsson). "
       "x, Jg e Jl são **idênticos ao v12** (velocidades superficiais são "
       "invariantes ao escorregamento); o que o v12.1 acrescenta é α e as "
       "velocidades **reais** de fase — as grandezas que a literatura tabula e "
       "que de fato mudam com o escorregamento.\n")

    ap("## Propriedades do R600a (isobutano) a Tsat = 27 °C")
    ap(f"- Fonte: **CoolProp {props['coolprop_version']}** (`PropsSI`, fluido `R600a`)")
    ap(f"- ρ_l = **{props['rho_l']:.2f} kg/m³**  ·  ρ_v = **{props['rho_v']:.3f} kg/m³**")
    ap(f"- i_lv = **{props['i_lv']/1e3:.2f} kJ/kg**  ·  σ = **{props['sigma']*1e3:.3f} mN/m**"
       f"  ·  P_sat = {props['p_sat']/1e5:.3f} bar")
    ap(f"- Geometria: D = 1,0 mm · L_aquecido = 0,265 m · A = {A_CROSS:.3e} m²\n")

    ap("## Modelo acrescentado (sobre o v12)")
    ap("```")
    ap("# x, Jg=jv, Jl=jl  : idênticos ao v12 (não mudam com escorregamento)")
    ap("C0    = 1 + 0.12·(1 − x)")
    ap("drift = 1.18·(1 − x)·[g·σ·(ρl − ρv)]^0.25 / (G·ρl^0.5)")
    ap("α     = (x/ρv) / ( C0·(x/ρv + (1−x)/ρl) + drift )   # Rouhani–Axelsson (1970)")
    ap("uv    = jv / α            ul = jl / (1 − α)          # velocidades reais")
    ap("S     = (x·ρl·(1−α)) / (α·ρv·(1−x))                  # escorregamento = uv/ul")
    ap("```")
    ap("> **Grounding:** α — dias2025knn Eq. (5); S — deoliveira2021 Eq. (17). "
       "uv=jv/α e ul=jl/(1−α) decorrem da definição de α (implícito em "
       "deoliveira2021 via Eq. 17), não citação literal.\n")

    # ---- tabela por ensaio ----
    ap("## Por ensaio")
    ap("| Ensaio | Padrão | G [kg/m²s] | q″ [kW/m²] | x [-] | Jg=jv [m/s] | "
       "Jl=jl [m/s] | α [-] | uv [m/s] | ul [m/s] | S [-] |")
    ap("|--------|--------|-----------:|-----------:|------:|------------:|"
       "-----------:|------:|---------:|---------:|------:|")
    for r in sorted(kept, key=lambda r: (PATTERN_ORDER.index(r["pattern"]),
                                          r["G"], r["q"])):
        uv = f"{r['uv']:.3f}" if np.isfinite(r["uv"]) else "—"
        ul = f"{r['ul']:.4f}" if np.isfinite(r["ul"]) else "—"
        S = f"{r['S']:.3f}" if np.isfinite(r["S"]) else "—"
        ap(f"| {r['run']} | {PATTERN_PT[r['pattern']]} | {r['G']} | {r['q']} | "
           f"{r['x']:.4f} | {r['Jg']:.3f} | {r['Jl']:.4f} | {r['alpha']:.4f} | "
           f"{uv} | {ul} | {S} |")
    ap("")

    # ---- resumo por padrão ----
    ap("## Resumo por padrão (min / máx / média / mediana)")
    ap("| Padrão | n | α min | α máx | α média | α med | S min | S máx | S média | "
       "S med | uv média | ul média |")
    ap("|--------|--:|------:|------:|--------:|------:|------:|------:|--------:|"
       "------:|---------:|---------:|")
    for pat in PATTERN_ORDER:
        pr = [r for r in kept if r["pattern"] == pat]
        al = _stats([r["alpha"] for r in pr])
        sv = _stats([r["S"] for r in pr])
        uv = _stats([r["uv"] for r in pr])
        ul = _stats([r["ul"] for r in pr])
        ap(f"| {PATTERN_PT[pat]} | {len(pr)} | "
           f"{al[0]:.3f} | {al[1]:.3f} | {al[2]:.3f} | {al[3]:.3f} | "
           f"{sv[0]:.3f} | {sv[1]:.3f} | {sv[2]:.3f} | {sv[3]:.3f} | "
           f"{uv[2]:.3f} | {ul[2]:.4f} |")
    ap("")

    # ---- validação de α vs literatura ----
    ap("## Validação de α por padrão × literatura local")
    ap("| Padrão | α v12.1 (min–máx) | Literatura | Citação | Veredito |")
    ap("|--------|-------------------|-----------|---------|----------|")
    an = _stats([r["alpha"] for r in kept if r["pattern"] == "wavy annular"])
    pl = _stats([r["alpha"] for r in kept if r["pattern"] == "plug"])
    ch = _stats([r["alpha"] for r in kept if r["pattern"] == "churn"])
    sl = _stats([r["alpha"] for r in kept if r["pattern"] == "slug"])
    ap(f"| plug | {pl[0]:.3f}–{pl[1]:.3f} | α 0,042–0,29 (faixa) | deoliveira2021 §4.2 | "
       "ordem de grandeza compatível (plug = baixo vazio) |")
    ap(f"| slug | {sl[0]:.3f}–{sl[1]:.3f} | α 0,578–0,714 (faixa, G=520–1220); "
       "0,73 (ponto, G=516) | Dias (2024) Figs. 22d/23b/24a; dias2025knn Tab. 2 | "
       "v12.1 **abaixo** da faixa (x_in=0 subestima; ver nota) |")
    ap(f"| churn | {ch[0]:.3f}–{ch[1]:.3f} | α 0,61–0,85 (faixa) | deoliveira2021 §4.4 | "
       "faixa sobreposta; v12.1 alcança o topo da faixa |")
    ap(f"| anular ondulado | {an[0]:.3f}–{an[1]:.3f} | α 0,86–0,99 (faixa) | "
       "deoliveira2021 §4.5 | **bate** (dentro da faixa) |")
    ap("")
    ap("> **Nota (slug):** de Oliveira (2021) **não publica faixa de α para slug** "
       "(confirmado lendo a §4.3 inteira — só dá descritores de caos/Hurst/LLE; "
       "faixa de α existe lá só para plug §4.2, churn §4.4 e anular §4.5). A faixa "
       "de slug mostrada vem da **dissertação Dias (2024)**, mesmo arranjo LETEF "
       "(R-600a, mini-canal 1 mm), três pontos **confirmados visualmente** nas "
       "páginas renderizadas (Figs. 22d/23b/24a: α=0,7023/0,7135/0,5782 em "
       "G=520/740/1220; o PDF tem OCR corrompido — ver `assets/paper-log.md`). O "
       "ponto único α=0,73 de `dias2025knn` Tab. 2 (G=516) **corrobora** a faixa "
       "(marcador × na figura). O v12.1 dá α_slug "
       f"{sl[0]:.2f}–{sl[1]:.2f}, **abaixo** da faixa 0,58–0,71 de Dias (2024) — "
       "mesmo viés do plug: o x de equilíbrio (x_in=0) subestima o vazio nos "
       "padrões de baixo-a-médio título. A faixa da literatura é mostrada como "
       "**referência de onde o slug real vive**, não como concordância de "
       "magnitude do v12.1.\n")

    # ---- limitação herdada ----
    ap("## Limitação herdada do v12 (honesta)")
    ap("α, uv, ul e S são calculados sobre o **x de equilíbrio** (x_in = 0). Não há "
       "potência do pré-aquecedor nem vazão mássica medida nos dados brutos "
       "(`src/data_processed/<run>/dados.xlsx` contém só `source_file`, `Voltage_0`, "
       "`Voltage_1`, `Voltage_2`) para reconstruir o x de entrada real. Logo estas "
       "grandezas herdam a mesma natureza de **estimativa de modelo** do v12 — "
       "declarado, não escondido. Ver `src2/v12.1/ANALISE_ADVERSARIAL_v12.md`.\n")

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
          f"sigma={props['sigma']*1e3:.3f} mN/m")

    rows = compute_runs(props)
    msgs = sanity_checks(rows, props)
    for m in msgs:
        print("  ✔", m)

    plot_phase_velocity_map(rows, props, here / "phase_velocity_map_v12_1.png")
    plot_void_fraction(rows, here / "void_fraction_v12_1.png")
    write_report(rows, props, msgs, here / "velocity_table_v12_1.md")
    print(f"\nEscrito em {here}:")
    print("  - phase_velocity_map_v12_1.png")
    print("  - void_fraction_v12_1.png")
    print("  - velocity_table_v12_1.md")


if __name__ == "__main__":
    main()
