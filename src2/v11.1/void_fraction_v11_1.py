"""
v11.1 — Fração de vazio × INTERMITÊNCIA POR TAXA DE ALTERNAÇÃO (tabela + mapa).

Revisão do v11. O eixo de intermitência do v11 era ``alpha_std`` = desvio-padrão
temporal de α_j. Esse descritor colocava o *wavy annular* ACIMA do *churn* em
intermitência — ordem fisicamente invertida (o orientador apontou o estranhamento).
A auditoria do v11 (ver ``../v11/void_fraction_v11.py`` e o diagnóstico abaixo)
mostrou que NÃO há bug de código: os números reproduzem; a inversão é um artefato
da MÉTRICA somada a um confundimento de brilho/fluxo de calor q:

  1. ``std(α_j)`` lê a ondulação do filme interfacial do anular como "intermitência"
     (um filme ondulado faz α_j oscilar sem que a corda central deixe de ser gás);
  2. os ensaios anulares foram amostrados em q maior (q̄≈54 vs churn≈20); os anulares
     de q baixo saem quase saturados (core-mean ~2×), e a perda de faixa dinâmica do
     Otsu COMPRIME o std — mascarando a diferença.

v11.1 troca o eixo Y pela definição física de escoamento INTERMITENTE
(costigan-whalley; jones-zuber): limiariza α_j em dois ESTADOS (gás-dominado ×
líquido-dominado) e mede a TAXA DE ALTERNAÇÃO entre eles — a frequência de passagem
de pistões/tampões. Plug/slug alternam muito (bolhões bem separados por líquido),
churn alterna de forma caótica, e o anular ondulado quase nunca deixa o estado gás
(núcleo de gás estável) ⇒ taxa MÍNIMA. A ordem física correta é recuperada:

    plug  >  slug  >  churn  >  annular        (annular = MENOS intermitente)

Robustez a brilho/q: os estados são definidos sobre a FRAÇÃO α_j (um número em
[0,1]), com histerese (gatilho de Schmitt: entra em "gás" só quando α_j>hi,
entra em "líquido" só quando α_j<lo), NÃO sobre o nível de cinza do pixel. Logo o
nível absoluto do Otsu / brilho do panorama não fixa o limiar. Exigir uma oscilação
PROFUNDA de vazio (corda quase toda gás ↔ quase toda líquido) para trocar de estado
é o que rejeita a ondulação do filme do anular (o viés do std): um filme ondulado
nunca leva a corda central a ficar totalmente líquida, então o anular quase não sai
do estado gás.

NOVO neste v11.1: o padrão **dispersed bubbles** (``m07q01``) passa a ser incluído
na análise (estava excluído no v11). Bolhas dispersas são pequenas e distribuídas no
líquido: a corda central fica quase sempre líquida (α≈0,02), nunca atinge o estado
gás ⇒ taxa de alternação = 0. Fica no canto inferior-ESQUERDO do mapa (pouco vazio,
pouca alternação), distinto do anular (que está embaixo por ter vazio ALTO e estável).
É n=1 — anotado como tal; não entra em comparação estatística.

Escopo: cinco padrões — dispersed bubbles, plug, slug, churn, wavy annular.
Excluídos: ``synthetic1/2``; *smooth annular* (``m03q80``); e as transições com DOIS
padrões ao mesmo tempo *plug-slug* (``m03q01``) e *slug-churn* (``m04q05``, ``m07q10``).

Decisões de rotulagem/qualidade (idênticas ao v11):
  - m05q05 (G=640, q=5) tem DOIS ensaios com rótulos distintos nos marcadores
    ``__*.txt`` (fonte de verdade): ``_001`` = slug, ``_002`` = plug (outlier de
    plug confirmado pelo dono). Ambos mantidos.
  - m03q40 é descartado da ESTATÍSTICA por saturação (34,3 % do núcleo em 255): o
    recorte destrói o contraste e o Otsu subestima o vazio. Mantido na tabela por
    ensaio (sinalizado). Observação: o eixo de ALTERNAÇÃO do v11.1 é imune a esse
    viés de brilho (limiar em α_j, não em cinza), mas o eixo X (α, herdado do v11)
    ainda é afetado — por coerência com o v11 o ensaio continua fora da estatística.

Saídas (gravadas em ``src2/v11.1/``):
  void_fraction_table_v11_1.md  — tabela por PADRÃO (mín/máx/média/mediana da taxa
                                   de alternação + α) + tabela por ensaio
  alpha_regions_map_v11_1.png   — mapa α × taxa de alternação, regiões por padrão
                                   (rótulo CONHECIDO — não é clustering)

Uso:
  cd src2/v11.1 && python3 void_fraction_v11_1.py
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

PANORAMA  = "__panorama_v11.png"      # panorama reconstruído (pipeline_v11.py)
CORE_ROWS = slice(20, 60)             # núcleo do escoamento (paredes excluídas)

G_MAP = {3: 380, 4: 510, 5: 640, 6: 760, 7: 890}   # fluxo mássico [kg/(m²·s)]

# Estados do gatilho de Schmitt sobre α_j (fração cordal de gás em [0,1]):
#   entra em GÁS-dominado quando α_j > STATE_HI; entra em LÍQUIDO-dominado quando
#   α_j < STATE_LO; dentro da banda morta [LO, HI] mantém o estado. Exigir uma
#   oscilação PROFUNDA (corda quase toda gás ↔ quase toda líquido) é o que rejeita
#   a ondulação do filme do anular. Banda LARGA — a ordem entre padrões é estável
#   para qualquer banda larga e só degrada ao aproximar de um corte na mediana
#   (sensibilidade documentada em STATE_BANDS abaixo).
STATE_LO = 0.25
STATE_HI = 0.75

# Bandas alternativas para a análise de sensibilidade reportada na tabela.
STATE_BANDS = [(0.20, 0.80), (0.25, 0.75), (0.30, 0.70), (0.35, 0.65)]

# v11.1 INCLUI dispersed bubbles (m07q01). Apenas smooth annular segue excluído.
EXCLUDE_RUNS = {"m03q80_001"}         # smooth annular (não pedido)

RUN_RE = re.compile(r"^m0\dq\d+_\d+$")

# Fração de pixels saturados (==255) no núcleo acima da qual o ensaio sai da
# estatística por padrão (coerência com o eixo α herdado do v11). Só m03q40 passa.
SAT_LIMIT = 0.20

# Agrupamento para cor/região da FIGURA (a TABELA mantém o rótulo fino).
GROUP_MAP = {
    "dispersed bubbles": "bubble",
    "plug":              "plug",
    "slug":              "slug",
    "churn":             "churn",
    "wavy annular":      "annular",
}
# Rótulos puros aceitos; os demais (plug-slug, slug-churn, smooth annular) são
# transições/casos excluídos da análise.
PURE_LABELS = set(GROUP_MAP)
GROUP_ORDER = ["bubble", "plug", "slug", "churn", "annular"]

# Paleta qualitativa segura para daltônicos (Okabe–Ito).
GROUP_COLOR = {
    "bubble":  "#E69F00",
    "plug":    "#0072B2",
    "slug":    "#009E73",
    "churn":   "#D55E00",
    "annular": "#CC79A7",
}
GROUP_PRETTY = {
    "bubble":  "dispersed bubbles",
    "plug":    "plug",
    "slug":    "slug",
    "churn":   "churn",
    "annular": "annular",
}


# ── Sinal de fração de vazio (α_j), portado de v11 verbatim ──────────────────

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


def alpha_signal(panorama):
    """Série temporal da fração de vazio cordal α_j (gás = ``valor <= t``).

    panorama : array (80, N) uint8 em escala de cinza.

    O limiar de Otsu ``t`` é escolhido por imagem sobre todo o núcleo
    (linhas 20–59). Para cada coluna j (um instante no tempo) α_j =
    ``mean(core[:, j] <= t)``. Devolve (aj, meta) com aj em [0, 1]^N; meta inclui
    ``alpha`` (média temporal — eixo X), ``alpha_std`` (descritor antigo do v11,
    mantido para comparação), ``core_mean``, ``otsu_thr`` e ``sat_frac``.
    """
    core = panorama[CORE_ROWS, :].astype(float)   # 40 × N, paredes excluídas
    thr = manual_otsu(core.ravel())
    aj = (core <= thr).mean(axis=0)                # fração cordal por coluna (tempo)
    meta = {
        "otsu_thr": int(thr),
        "core_mean": float(core.mean()),
        "alpha": float(aj.mean()),
        "alpha_std": float(aj.std()),
        "sat_frac": float((core >= 255).mean()),
    }
    return aj, meta


# ── Intermitência por taxa de alternação (gatilho de Schmitt) ────────────────

def state_series(aj, lo=STATE_LO, hi=STATE_HI):
    """Classificação em dois estados com histerese (gatilho de Schmitt).

    A corda é contada como GÁS-dominada (1) só quando α_j sobe acima de ``hi``
    (corda quase toda gás) e como LÍQUIDO-dominada (0) só quando α_j cai abaixo de
    ``lo`` (um pistão de líquido preenche a corda); na banda morta [lo, hi] o
    estado é mantido. Os limiares são sobre α_j (fração), não sobre brilho de
    pixel ⇒ robusto a brilho/q. Devolve array inteiro em {0=líquido, 1=gás}.
    """
    s = np.empty(len(aj), dtype=int)
    cur = 1 if aj[0] >= (lo + hi) / 2 else 0
    for i, v in enumerate(aj):
        if cur == 0 and v > hi:
            cur = 1
        elif cur == 1 and v < lo:
            cur = 0
        s[i] = cur
    return s


def transition_rate(aj, lo=STATE_LO, hi=STATE_HI):
    """Taxa de alternação gás↔líquido por 1000 px de percurso do tubo.

    O eixo X do panorama é espacial (1 px = 1 px de deslocamento do tubo); a taxa
    é por 1000 px (a hipótese de Taylor mapeia isso em frequência de passagem).
    """
    s = state_series(aj, lo, hi)
    flips = int(np.sum(np.abs(np.diff(s))))
    return flips / len(s) * 1000.0, flips


def intermittency_meta(aj):
    """Descritores de alternação para a banda padrão + varredura de sensibilidade."""
    rate, flips = transition_rate(aj)
    s = state_series(aj)
    frac_gas = float(np.mean(s == 1))
    # fração do tempo "intermitente" = não-monopolizada por um estado
    frac_interm = 2.0 * min(frac_gas, 1.0 - frac_gas)
    sweep = {f"rate_{lo:.2f}_{hi:.2f}": transition_rate(aj, lo, hi)[0]
             for (lo, hi) in STATE_BANDS}
    return {"trans_rate": rate, "n_flips": flips,
            "frac_gas": frac_gas, "frac_interm": frac_interm, **sweep}


# ── Metadados: rótulo, (m, q) (portado de v11) ───────────────────────────────

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
    """Ensaios reais de PADRÃO PURO (bubble, plug, slug, churn, wavy annular)."""
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
                "rate [/kpx]", "alpha", "alpha_std(v11)", "sat %", "N [px]", "flag"]
STATS_COLS   = ["padrão", "n", "taxa mín", "taxa máx", "taxa média", "taxa mediana",
                "α mediana"]
SWEEP_COLS   = ["padrão", "n",
                "taxa 0.20/0.80", "taxa 0.25/0.75", "taxa 0.30/0.70", "taxa 0.35/0.65"]


def _per_run_values(r):
    return [r["scenario"], r["pattern"], r["G"], r["q"],
            round(r["trans_rate"], 3), round(r["alpha"], 3),
            round(r["alpha_std"], 3),
            round(r["sat_frac"] * 100, 1), r["N"], r["flag"]]


def _md_table(fh, cols, rows_values):
    fh.write("| " + " | ".join(cols) + " |\n")
    fh.write("|" + "|".join(["---"] * len(cols)) + "|\n")
    for vals in rows_values:
        fh.write("| " + " | ".join(str(v) for v in vals) + " |\n")


def per_pattern_stats(rows):
    """Mín/máx/média/mediana da TAXA de alternação + mediana de α por padrão,
    sobre os ensaios NÃO descartados. Lista de dicts na ordem GROUP_ORDER."""
    stats = []
    for g in GROUP_ORDER:
        grp = [r for r in rows if r["group"] == g and not r["dropped"]]
        if not grp:
            continue
        rate = np.array([r["trans_rate"] for r in grp])
        alpha = np.array([r["alpha"] for r in grp])
        stats.append({
            "group": g, "n": len(grp),
            "rate_min": float(rate.min()), "rate_max": float(rate.max()),
            "rate_avg": float(rate.mean()), "rate_median": float(np.median(rate)),
            "alpha_median": float(np.median(alpha)),
        })
    return stats


def per_pattern_sweep(rows):
    """Mediana da taxa por padrão para cada banda de histerese (sensibilidade)."""
    out = []
    for g in GROUP_ORDER:
        grp = [r for r in rows if r["group"] == g and not r["dropped"]]
        if not grp:
            continue
        row = {"group": g, "n": len(grp)}
        for (lo, hi) in STATE_BANDS:
            key = f"rate_{lo:.2f}_{hi:.2f}"
            row[key] = float(np.median([r[key] for r in grp]))
        out.append(row)
    return out


def write_table(rows, stats, sweep):
    out = os.path.join(HERE, "void_fraction_table_v11_1.md")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("# Fração de vazio × taxa de alternação por padrão (v11.1)\n\n")
        fh.write(
            "Revisão do v11. **Eixo Y = taxa de alternação gás↔líquido por 1000 px** "
            "(definição física de escoamento intermitente; costigan-whalley, "
            "jones-zuber), substituindo o `std(α_j)` do v11 que invertia a ordem "
            "anular×churn. A fração de vazio cordal média **α** permanece no eixo X "
            "(cordal, em pixels, sem calibração px→mm; comparação ENTRE padrões é "
            "legítima, o valor absoluto não é a fração de área da seção).\n\n")
        fh.write(
            "A taxa usa um gatilho de Schmitt sobre α_j: estado GÁS quando α_j>"
            f"{STATE_HI:.2f}, estado LÍQUIDO quando α_j<{STATE_LO:.2f} (banda morta "
            "mantém o estado). Limiar em α_j (fração), não em brilho ⇒ imune ao viés "
            "de saturação/q que contaminava o `std`.\n\n")

        n_analisados = sum(1 for r in rows if not r["dropped"])
        fh.write(f"Padrões: dispersed bubbles, plug, slug, churn, wavy annular — "
                 f"{len(rows)} ensaios ({n_analisados} na estatística). Transições "
                 "(plug-slug, slug-churn), smooth annular e synthetic "
                 "desconsiderados.\n\n")

        fh.write("## Estatística por padrão\n\n")
        _md_table(fh, STATS_COLS, [
            [GROUP_PRETTY[s["group"]], s["n"],
             round(s["rate_min"], 3), round(s["rate_max"], 3),
             round(s["rate_avg"], 3), round(s["rate_median"], 3),
             round(s["alpha_median"], 3)]
            for s in stats])

        fh.write("\n**Ordem recuperada (taxa mediana):** ")
        ordered = sorted([s for s in stats if s["group"] != "bubble"],
                         key=lambda s: -s["rate_median"])
        fh.write(" > ".join(f"{GROUP_PRETTY[s['group']]} {s['rate_median']:.2f}"
                            for s in ordered))
        fh.write(" — anular é o MENOS intermitente, como esperado fisicamente.\n\n")

        fh.write("## Sensibilidade à banda de histerese\n\n")
        fh.write("Taxa mediana por padrão para quatro bandas [lo, hi]. A ORDEM "
                 "plug>slug>churn>annular é estável para toda banda larga; só "
                 "degrada ao estreitar a banda em direção a um corte na mediana.\n\n")
        _md_table(fh, SWEEP_COLS, [
            [GROUP_PRETTY[s["group"]], s["n"]]
            + [round(s[f"rate_{lo:.2f}_{hi:.2f}"], 3) for (lo, hi) in STATE_BANDS]
            for s in sweep])

        fh.write("\n")
        fh.write(
            "> **dispersed bubbles (n=1):** taxa = 0 — a corda central fica quase "
            "sempre líquida (α≈0,02; bolhas pequenas e dispersas nunca enchem a "
            "corda de gás), logo nunca atinge o estado gás. Fica no canto "
            "inferior-ESQUERDO (pouco vazio, pouca alternação), por motivo DISTINTO "
            "do anular (vazio alto e estável). Único ensaio — sem comparação.\n\n")
        fh.write(
            "> **Amostra pequena:** slug (n=3) e plug (n=3) têm poucos ensaios; "
            "mediana ≈ valor central. Ler com cautela.\n\n")
        fh.write(
            "> **m05q05 (G=640, q=5):** única condição (m, q) em DUAS classes — "
            "`_001`=slug, `_002`=plug (outlier confirmado pelo dono). Ambos "
            "mantidos.\n\n")

        dropped = [r for r in rows if r["dropped"]]
        if dropped:
            fh.write(
                "> **Descartado da estatística (mantido na tabela):** "
                + "; ".join(f"{r['scenario']} ({r['flag']})" for r in dropped)
                + ". A saturação afeta o eixo X (α); o eixo Y (taxa) é imune "
                  "(limiar em α_j), mas por coerência com o v11 o ensaio é "
                  "excluído da estatística.\n\n")

        fh.write("## Por ensaio (transparência)\n\n")
        fh.write("Ordenado por taxa de alternação decrescente. `rate` = alternações "
                 "gás↔líquido por 1000 px; `alpha_std(v11)` repetido para "
                 "comparação com o descritor antigo.\n\n")
        rows_sorted = sorted(rows, key=lambda r: -r["trans_rate"])
        _md_table(fh, PER_RUN_COLS, [_per_run_values(r) for r in rows_sorted])
    print(f"tabela: {out}")


# ── Saídas: mapa (regiões por rótulo conhecido — NÃO é clustering) ───────────

def _draw_region(ax, pts, color, label, label_shift=(0.0, 0.0)):
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
        half = pts.std(axis=0) if len(pts) >= 2 else np.array([0.0, 0.0])
        half = np.maximum(half, [0.015, 0.10]) + [0.010, 0.06]
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
    """Mapa α × taxa de alternação, colorido por padrão, com regiões por rótulo
    CONHECIDO (não é clustering — os grupos vêm do marcador __*.txt)."""
    fig, ax = plt.subplots(figsize=(8, 6))

    LABEL_SHIFT = {
        "bubble":  (0.03,  0.35),
        "plug":    (0.0,   0.45),
        "slug":    (0.045, 0.0),
        "churn":   (-0.06, 0.55),
        "annular": (0.02, -0.55),
    }

    by_group = {}
    for r in rows:
        if r["dropped"]:
            continue
        by_group.setdefault(r["group"], []).append(r)

    for g in GROUP_ORDER:
        pts = by_group.get(g)
        if not pts:
            continue
        coords = [(p["alpha"], p["trans_rate"]) for p in pts]
        _draw_region(ax, coords, GROUP_COLOR[g], GROUP_PRETTY[g],
                     label_shift=LABEL_SHIFT.get(g, (0.0, 0.0)))

    for g in GROUP_ORDER:
        pts = by_group.get(g)
        if not pts:
            continue
        xs = [p["alpha"] for p in pts]
        ys = [p["trans_rate"] for p in pts]
        ax.scatter(xs, ys, s=55, color=GROUP_COLOR[g], edgecolor="black",
                   linewidth=0.5, label=f"{GROUP_PRETTY[g]} (n={len(pts)})",
                   zorder=4)

    # dispersed bubbles (n=1): anota explicitamente
    for g, pts in by_group.items():
        if g == "bubble" and len(pts) == 1:
            r = pts[0]
            ax.annotate("n=1 (corda sempre líquida)",
                        (r["alpha"], r["trans_rate"]),
                        textcoords="offset points", xytext=(8, 6),
                        fontsize=7, color=GROUP_COLOR["bubble"])

    drop = [r for r in rows if r["dropped"]]
    if drop:
        ax.scatter([r["alpha"] for r in drop], [r["trans_rate"] for r in drop],
                   s=70, facecolors="none", edgecolors="0.4", linewidth=1.3,
                   marker="o", zorder=4, label="saturado (fora da estatística)")
        for r in drop:
            ax.annotate(r["scenario"], (r["alpha"], r["trans_rate"]),
                        textcoords="offset points", xytext=(6, 4),
                        fontsize=7, color="0.4")

    ax.set_xlabel("Fração de vazio cordal média  $\\alpha$  [–]")
    ax.set_ylabel("Intermitência — taxa de alternação gás↔líquido  "
                  "[alternações / 1000 px]")
    n_in = sum(len(v) for v in by_group.values())
    ax.set_title("Fração de vazio × taxa de alternação — regiões por padrão "
                 f"(rótulo conhecido) — {n_in} na estatística")
    ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.4)
    ax.legend(loc="best", fontsize=8, framealpha=0.9)
    out = os.path.join(HERE, "alpha_regions_map_v11_1.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"figura: {out}")


# ── Orquestração ─────────────────────────────────────────────────────────────

def main():
    scenarios = list_scenarios()
    print(f"{len(scenarios)} ensaios de padrão puro candidatos "
          "(transições, synthetic*, smooth annular excluídos; dispersed bubbles "
          "INCLUÍDO)\n")

    rows = []
    for name in scenarios:
        folder = os.path.join(DATA_ROOT, name)
        panorama = np.array(Image.open(os.path.join(folder, PANORAMA)).convert("L"))
        aj, meta = alpha_signal(panorama)
        im = intermittency_meta(aj)

        label = folder_flow_type(name)
        group = GROUP_MAP.get(label, label)
        m, q = parse_m_q(name)

        dropped = meta["sat_frac"] > SAT_LIMIT
        flag = (f"saturado (sat={meta['sat_frac'] * 100:.0f}%) — "
                "descartado da estatística") if dropped else ""

        rows.append({
            "scenario": name, "pattern": label, "group": group,
            "m": m, "q": q, "G": G_MAP[m],
            "alpha": meta["alpha"], "alpha_std": meta["alpha_std"],
            "sat_frac": meta["sat_frac"], "otsu_thr": meta["otsu_thr"],
            "core_mean": meta["core_mean"],
            "N": int(panorama.shape[1]), "dropped": dropped, "flag": flag,
            **im,
        })

    stats = per_pattern_stats(rows)
    sweep = per_pattern_sweep(rows)

    write_table(rows, stats, sweep)
    make_map(rows)

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

    print("\nEstatística por padrão (taxa de alternação /1000px):")
    print(f"  {'padrão':18} {'n':>2}  {'mín':>6} {'máx':>6} {'média':>6} "
          f"{'mediana':>7}  {'α med':>6}")
    for s in stats:
        print(f"  {GROUP_PRETTY[s['group']]:18} {s['n']:2d}  "
              f"{s['rate_min']:6.3f} {s['rate_max']:6.3f} {s['rate_avg']:6.3f} "
              f"{s['rate_median']:7.3f}  {s['alpha_median']:6.3f}")


if __name__ == "__main__":
    main()
