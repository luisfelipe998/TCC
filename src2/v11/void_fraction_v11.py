"""
v11 — Fração de vazio por padrão de escoamento (tabela + mapa de regiões).

Para cada ensaio de PADRÃO PURO em ``src2/data_processed/mXXqYY_*`` calcula a
fração de vazio α diretamente do panorama reconstruído (``__panorama_v11.png``)
e produz:

  1. Fração de vazio **cordal média no tempo** α por limiar de Otsu no núcleo do
     panorama (linhas 20–59, paredes excluídas; gás = lado escuro, ``valor <= t``,
     consistente com a classe que ``manual_otsu`` acumula). Para cada coluna j
     (um instante) α_j = fração cordal de gás; α = média temporal de α_j. Como o
     panorama é uma reconstrução espaço-tempo (altura = corda vertical; colunas =
     tempo), esta é uma fração de vazio CORDAL (uma linha), NÃO uma fração areal
     de seção transversal; mantém o viés cordal>areal conhecido para núcleo
     gasoso central (Hewitt; sondas óticas de feixe único, Kataoka–Ishii–
     Serizawa). Em pixels, sem calibração px→mm. A comparação ENTRE padrões é
     legítima; o valor absoluto não é a fração de área da seção.

  2. ``alpha_std`` = intermitência temporal (desvio de α_j ao longo do tempo).

Escopo: apenas os quatro padrões PUROS pedidos — plug, slug, churn e wavy
annular. Excluídos: ``synthetic1/2``; *smooth annular* (``m03q80``); *dispersed
bubbles* (``m07q01``); e as transições com DOIS padrões ao mesmo tempo
*plug-slug* (``m03q01``) e *slug-churn* (``m04q05``, ``m07q10``). ⇒ 26 ensaios
(plug 4, slug 3, churn 10, wavy annular 9).

Decisões de rotulagem/qualidade:
  - m05q05 (G=640, q=5) tem DOIS ensaios com rótulos distintos nos marcadores
    ``__*.txt`` (fonte de verdade): ``_001`` = slug, ``_002`` = plug (outlier de
    plug confirmado pelo dono). Ambos são mantidos; é a única condição (m, q) que
    aparece em duas classes — anotado na tabela e no mapa.
  - m03q40 é descartado da ESTATÍSTICA por saturação: 34,3 % do núcleo satura em
    255 (vs ≤14 % nos demais) e core-mean ~198 (≈2× o restante). O recorte do
    sinal destrói o contraste gás/líquido no lado claro e o Otsu classifica gás
    saturado como líquido ⇒ α=0,275 artificialmente baixo (menor do conjunto).
    Correções de de-trend/saturação/limiar local já foram testadas e descartadas
    (v10): recorte é perda de informação, não recuperável por re-limiarização.
    A linha é MANTIDA na tabela por ensaio (sinalizada), apenas excluída das
    estatísticas por padrão (annular passa a n=8).

Saídas (gravadas em ``src2/v11/``):
  void_fraction_table_v11.md  — tabela por PADRÃO (mín/máx/média/mediana) +
                                tabela por ensaio (transparência)
  alpha_regions_map_v11.png   — mapa α × intermitência, regiões por padrão
                                (rótulo CONHECIDO — não é clustering)

Uso:
  cd src2/v11 && python3 void_fraction_v11.py
"""

import os
import re
import glob

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from PIL import Image
from scipy.spatial import ConvexHull, QhullError

# ── Constantes ────────────────────────────────────────────────────────────────
HERE      = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.join(HERE, "..", "data_processed")

PANORAMA  = "__panorama_v11.png"      # nome do panorama (pipeline_v11.py)
CORE_ROWS = slice(20, 60)             # núcleo do escoamento (paredes excluídas)

G_MAP = {3: 380, 4: 510, 5: 640, 6: 760, 7: 890}   # fluxo mássico [kg/(m²·s)]

# Ensaios reais excluídos por decisão (não são padrão puro dos quatro pedidos).
EXCLUDE_RUNS = {"m03q80_001",   # smooth annular
                "m07q01_001"}   # dispersed bubbles

RUN_RE = re.compile(r"^m0\dq\d+_\d+$")

# Fração de pixels saturados (==255) no núcleo acima da qual o ensaio sai da
# estatística por padrão (Otsu deixa de separar gás/líquido). Só m03q40 passa.
SAT_LIMIT = 0.20

# Agrupamento para cor/região da FIGURA (a TABELA mantém o rótulo fino).
GROUP_MAP = {
    "plug":          "plug",
    "slug":          "slug",
    "churn":         "churn",
    "wavy annular":  "annular",
}
# Rótulos puros aceitos; os demais (plug-slug, slug-churn, smooth annular,
# dispersed bubbles) são transições/casos excluídos da análise.
PURE_LABELS = set(GROUP_MAP)
GROUP_ORDER = ["plug", "slug", "churn", "annular"]

# Paleta qualitativa segura para daltônicos (Okabe–Ito).
GROUP_COLOR = {
    "plug":    "#0072B2",
    "slug":    "#009E73",
    "churn":   "#D55E00",
    "annular": "#CC79A7",
}


# ── Fração de vazio (portado de v10, verbatim) ──────────────────────────────

def manual_otsu(vals):
    """Limiar de Otsu (256 níveis) sobre valores 0–255. Só numpy.

    Devolve o nível ``t`` que maximiza a variância entre classes; a classe
    "escura" é ``valor <= t`` (o próprio bin ``t`` é acumulado na classe escura
    por ``wB += hist[t]`` antes do teste, então o corte é inclusivo).
    """
    hist = np.histogram(vals, bins=256, range=(0, 256))[0].astype(float)
    total = hist.sum()
    sum_all = np.dot(np.arange(256), hist)

    wB = 0.0
    sumB = 0.0
    best = -1.0
    thr = 0
    for t in range(256):
        wB += hist[t]
        if wB == 0:
            continue
        wF = total - wB
        if wF == 0:
            break
        sumB += t * hist[t]
        mB = sumB / wB
        mF = (sum_all - sumB) / wF
        between = wB * wF * (mB - mF) ** 2
        if between > best:
            best = between
            thr = t
    return thr


def void_fraction(panorama):
    """Fração de vazio cordal média no tempo (gás = lado escuro, ``valor <= t``).

    panorama : array (80, N) uint8 em escala de cinza.

    O limiar de Otsu ``t`` é escolhido por imagem sobre todo o núcleo
    (linhas 20–59). Para cada coluna j (um instante no tempo) a fração de gás é
    ``alpha_j = mean(core[:, j] <= t)`` (fração cordal da seção), e a fração de
    vazio do ensaio é a média temporal ``alpha = mean_j(alpha_j)``.

    A comparação é ``<= t`` (não ``< t``): ``manual_otsu`` acumula o próprio bin
    ``t`` na classe escura, então a classe de gás é ``intensidade <= t``.

    Devolve (alpha, meta) com alpha em [0, 1]; meta inclui ``alpha_std``
    (intermitência temporal), ``core_mean``, ``otsu_thr`` e ``sat_frac``
    (fração de pixels saturados em 255 no núcleo).
    """
    core = panorama[CORE_ROWS, :].astype(float)   # 40 × N, paredes excluídas
    thr = manual_otsu(core.ravel())
    aj = (core <= thr).mean(axis=0)                # fração cordal por coluna (tempo)
    alpha = float(aj.mean())                       # média temporal
    return alpha, {
        "otsu_thr": int(thr),
        "core_mean": float(core.mean()),
        "alpha_std": float(aj.std()),
        "sat_frac": float((core >= 255).mean()),
    }


# ── Metadados: rótulo, (m, q) (portado de v10) ──────────────────────────────

def parse_m_q(run):
    """mXXqYY_ZZZ -> (m, q). Ex.: m05q10_001 -> (5, 10)."""
    mobj = re.match(r"m0?(\d+)q0?(\d+)", run)
    return int(mobj.group(1)), int(mobj.group(2))


def folder_flow_type(name):
    """Rótulo autoritativo pelo marcador __<rótulo>.txt.

    O ``sorted`` seleciona 'wavy annular' (espaço) antes de 'wavy_annular'
    (sublinhado), deduplicando os marcadores redundantes de m04q80/m05q80.
    """
    folder = os.path.join(DATA_ROOT, name)
    for p in sorted(glob.glob(os.path.join(folder, "*.txt"))):
        stem = os.path.splitext(os.path.basename(p))[0].lstrip("_")
        if stem != "dados":
            return stem
    return name


def list_scenarios():
    """Ensaios reais de PADRÃO PURO (plug, slug, churn, wavy annular).

    Exclui synthetic*, EXCLUDE_RUNS e qualquer ensaio cujo rótulo seja de
    transição/caso especial (plug-slug, slug-churn, smooth annular, dispersed
    bubbles) — identificado pelo marcador __<rótulo>.txt não estar em
    PURE_LABELS.
    """
    runs = []
    for d in sorted(os.listdir(DATA_ROOT)):
        if not RUN_RE.match(d) or d in EXCLUDE_RUNS:
            continue
        if not os.path.exists(os.path.join(DATA_ROOT, d, PANORAMA)):
            continue
        if folder_flow_type(d) not in PURE_LABELS:
            continue
        runs.append(d)
    return runs


# ── Saídas: tabelas ──────────────────────────────────────────────────────────

PER_RUN_COLS = ["scenario", "pattern", "G [kg/m2 s]", "q [kW/m2]",
                "alpha", "alpha_std", "sat %", "N [px]", "flag"]
STATS_COLS   = ["padrão", "n", "α mín", "α máx", "α média", "α mediana"]


def _per_run_values(r):
    return [r["scenario"], r["pattern"], r["G"], r["q"],
            round(r["alpha"], 3), round(r["alpha_std"], 3),
            round(r["sat_frac"] * 100, 1), r["N"], r["flag"]]


def _md_table(fh, cols, rows_values):
    fh.write("| " + " | ".join(cols) + " |\n")
    fh.write("|" + "|".join(["---"] * len(cols)) + "|\n")
    for vals in rows_values:
        fh.write("| " + " | ".join(str(v) for v in vals) + " |\n")


def per_pattern_stats(rows):
    """Mín/máx/média/mediana de α por padrão, sobre os ensaios NÃO descartados.

    Devolve lista de dicts na ordem GROUP_ORDER (só grupos presentes).
    """
    stats = []
    for g in GROUP_ORDER:
        vals = [r["alpha"] for r in rows if r["group"] == g and not r["dropped"]]
        if not vals:
            continue
        a = np.array(vals)
        stats.append({
            "group": g, "n": len(vals),
            "min": float(a.min()), "max": float(a.max()),
            "avg": float(a.mean()), "median": float(np.median(a)),
        })
    return stats


def write_table(rows, stats):
    """Grava a tabela por PADRÃO (destaque) + a tabela por ensaio (transparência)."""
    out = os.path.join(HERE, "void_fraction_table_v11.md")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("# Fração de vazio por padrão de escoamento (v11)\n\n")
        fh.write(
            "Fração de vazio **cordal média no tempo** α por limiar de Otsu no "
            "núcleo do panorama (linhas 20–59, paredes excluídas; gás = `valor "
            "<= t`). Para cada coluna (instante) α_j = fração cordal de gás; α = "
            "média temporal de α_j. Por ser **cordal** (uma linha, não a seção), "
            "α carrega o viés cordal>areal conhecido para núcleo gasoso central "
            "(Hewitt; Kataoka–Ishii–Serizawa); em pixels, sem calibração px→mm. "
            "A comparação ENTRE padrões é legítima; o valor absoluto **não** é a "
            "fração de área da seção transversal.\n\n")

        n_analisados = sum(1 for r in rows if not r["dropped"])
        fh.write(f"Padrões: plug, slug, churn, wavy annular — {len(rows)} ensaios "
                 f"({n_analisados} na estatística). Transições (plug-slug, "
                 "slug-churn), smooth annular, dispersed bubbles e synthetic "
                 "foram desconsiderados.\n\n")

        fh.write("## Estatística por padrão\n\n")
        _md_table(fh, STATS_COLS, [
            [s["group"], s["n"], round(s["min"], 3), round(s["max"], 3),
             round(s["avg"], 3), round(s["median"], 3)]
            for s in stats])

        fh.write("\n")
        fh.write(
            "> **Amostra pequena:** slug (n=3) e plug (n=4) têm poucos ensaios — "
            "a mediana é praticamente o valor central e `α mín`/`α máx` são os "
            "próprios extremos amostrais, não estimativas. Ler com cautela; sem "
            "teste estatístico.\n\n")
        fh.write(
            "> **m05q05 (G=640, q=5):** única condição (m, q) presente em DUAS "
            "classes — `m05q05_001` = slug e `m05q05_002` = plug, conforme os "
            "marcadores `__*.txt` (fonte de verdade; `_002` é outlier de plug "
            "confirmado pelo dono). Ambos mantidos.\n\n")

        dropped = [r for r in rows if r["dropped"]]
        if dropped:
            fh.write(
                "> **Descartados da estatística (mantidos na tabela por ensaio):** "
                + "; ".join(f"{r['scenario']} ({r['flag']})" for r in dropped)
                + ". Saturação do panorama (núcleo recortado em 255) impede o Otsu "
                  "de separar gás/líquido no lado claro — α fica artificialmente "
                  "baixo. Correções de de-trend/saturação/limiar local foram "
                  "testadas e descartadas (perda de informação por recorte).\n\n")

        fh.write("## Fração de vazio por ensaio (transparência)\n\n")
        fh.write("Ordenado por α crescente. `sat %` = fração do núcleo saturada "
                 "em 255; `flag` sinaliza ensaio fora da estatística.\n\n")
        rows_sorted = sorted(rows, key=lambda r: r["alpha"])
        _md_table(fh, PER_RUN_COLS, [_per_run_values(r) for r in rows_sorted])
    print(f"tabela: {out}")


# ── Saídas: mapa (regiões por rótulo conhecido — NÃO é clustering) ───────────

def _draw_region(ax, pts, color, label, label_shift=(0.0, 0.0)):
    """Desenha a região (fecho convexo) de um grupo + rótulo, em coordenadas de
    dados (ambos os eixos lineares neste mapa).

    label_shift : deslocamento (dx, dy) do rótulo em unidades de dados.
    """
    pts = np.asarray(pts, dtype=float)
    centroid = pts.mean(axis=0)

    drew_hull = False
    if len(pts) >= 3:
        try:
            hull = ConvexHull(pts)
            verts = pts[hull.vertices]
            ax.fill(verts[:, 0], verts[:, 1], color=color, alpha=0.15,
                    zorder=1, lw=0)
            closed = np.vstack([verts, verts[0]])
            ax.plot(closed[:, 0], closed[:, 1], color=color, lw=1.2,
                    alpha=0.6, zorder=2)
            drew_hull = True
        except QhullError:
            drew_hull = False

    if not drew_hull:
        # fallback: elipse compacta (grupos com <3 pontos ou colineares).
        # Pisos pequenos p/ eixos lineares em [0,~0.6]×[0,~0.3].
        half = pts.std(axis=0) if len(pts) >= 2 else np.array([0.0, 0.0])
        half = np.maximum(half, [0.015, 0.010]) + [0.010, 0.006]
        ell = mpatches.Ellipse(
            centroid, width=2 * half[0], height=2 * half[1],
            facecolor=color, alpha=0.15, edgecolor=color, lw=1.2, zorder=1)
        ax.add_patch(ell)

    lx = centroid[0] + label_shift[0]
    ly = centroid[1] + label_shift[1]
    txt = ax.text(lx, ly, label, color=color, fontsize=11, fontweight="bold",
                  ha="center", va="center", zorder=5)
    txt.set_path_effects([pe.withStroke(linewidth=2.5, foreground="white")])


def make_map(rows):
    """Mapa α × intermitência std(α_j), colorido por padrão, com regiões por
    rótulo CONHECIDO (não é clustering — os grupos vêm do marcador __*.txt).

    Ensaios descartados da estatística (saturados) aparecem como marcador vazado
    cinza, fora das regiões.
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    LABEL_SHIFT = {
        "plug":    (0.0,   0.030),
        "slug":    (0.0,  -0.035),
        "churn":   (0.0,  -0.035),
        "annular": (0.035, 0.025),
    }

    by_group = {}
    for r in rows:
        if r["dropped"]:
            continue
        by_group.setdefault(r["group"], []).append(r)

    # regiões primeiro (zorder baixo), só sobre ensaios NÃO descartados
    for g in GROUP_ORDER:
        pts = by_group.get(g)
        if not pts:
            continue
        coords = [(p["alpha"], p["alpha_std"]) for p in pts]
        _draw_region(ax, coords, GROUP_COLOR[g], g,
                     label_shift=LABEL_SHIFT.get(g, (0.0, 0.0)))

    # pontos por cima
    for g in GROUP_ORDER:
        pts = by_group.get(g)
        if not pts:
            continue
        xs = [p["alpha"] for p in pts]
        ys = [p["alpha_std"] for p in pts]
        ax.scatter(xs, ys, s=55, color=GROUP_COLOR[g], edgecolor="black",
                   linewidth=0.5, label=f"{g} (n={len(pts)})", zorder=4)

    # ensaios descartados (saturados): marcador vazado cinza
    drop = [r for r in rows if r["dropped"]]
    if drop:
        ax.scatter([r["alpha"] for r in drop], [r["alpha_std"] for r in drop],
                   s=70, facecolors="none", edgecolors="0.4", linewidth=1.3,
                   marker="o", zorder=4,
                   label="saturado (fora da estatística)")
        for r in drop:
            ax.annotate(r["scenario"], (r["alpha"], r["alpha_std"]),
                        textcoords="offset points", xytext=(6, 4),
                        fontsize=7, color="0.4")

    ax.set_xlabel("Fração de vazio cordal média  $\\alpha$  [–]")
    ax.set_ylabel("Intermitência temporal  $\\mathrm{std}(\\alpha_j)$  [–]")
    n_in = sum(len(v) for v in by_group.values())
    ax.set_title("Fração de vazio × intermitência — regiões por padrão "
                 f"(rótulo conhecido) — {n_in} ensaios")
    ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.4)
    ax.legend(loc="best", fontsize=8, framealpha=0.9)
    out = os.path.join(HERE, "alpha_regions_map_v11.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"figura: {out}")


# ── Orquestração ─────────────────────────────────────────────────────────────

def main():
    scenarios = list_scenarios()
    print(f"{len(scenarios)} ensaios de padrão puro candidatos "
          "(transições, synthetic*, smooth annular e dispersed bubbles "
          "excluídos)\n")

    rows = []
    for name in scenarios:
        folder = os.path.join(DATA_ROOT, name)
        panorama = np.array(Image.open(os.path.join(folder, PANORAMA)).convert("L"))
        alpha, meta = void_fraction(panorama)

        label = folder_flow_type(name)
        group = GROUP_MAP.get(label, label)
        m, q = parse_m_q(name)

        dropped = meta["sat_frac"] > SAT_LIMIT
        flag = (f"saturado (sat={meta['sat_frac'] * 100:.0f}%) — "
                "descartado da estatística") if dropped else ""

        rows.append({
            "scenario": name, "pattern": label, "group": group,
            "m": m, "q": q, "G": G_MAP[m], "alpha": alpha,
            "alpha_std": meta["alpha_std"], "sat_frac": meta["sat_frac"],
            "otsu_thr": meta["otsu_thr"], "core_mean": meta["core_mean"],
            "N": int(panorama.shape[1]), "dropped": dropped, "flag": flag,
        })

    stats = per_pattern_stats(rows)

    write_table(rows, stats)
    make_map(rows)

    # fail-loud: contagem analisada deve bater com a agrupada
    n = len(rows)
    counts = {g: sum(1 for r in rows if r["group"] == g) for g in GROUP_ORDER}
    print(f"\n{n} ensaios ANALISADOS "
          + "(" + " + ".join(f"{g} {counts[g]}" for g in GROUP_ORDER) + ")")
    if n != sum(counts.values()):
        raise RuntimeError(
            f"contagem inconsistente: {n} linhas vs {sum(counts.values())} agrupadas.")

    dropped = [r for r in rows if r["dropped"]]
    if dropped:
        print(f"{len(dropped)} descartado(s) da estatística: "
              + "; ".join(f"{r['scenario']} (sat={r['sat_frac'] * 100:.0f}%)"
                          for r in dropped))

    print("\nEstatística por padrão (α):")
    print(f"  {'padrão':9} {'n':>2}  {'mín':>6} {'máx':>6} {'média':>6} {'mediana':>7}")
    for s in stats:
        print(f"  {s['group']:9} {s['n']:2d}  {s['min']:6.3f} {s['max']:6.3f} "
              f"{s['avg']:6.3f} {s['median']:7.3f}")


if __name__ == "__main__":
    main()
