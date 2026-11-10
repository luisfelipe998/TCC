"""
v11.3 — núcleo compartilhado das duas análises de fração de vazio × intermitência.

Reúne o que as duas análises irmãs do v11.3 reusam, portado VERBATIM de
``../v11.2/nonchordal_void_v11_2.py`` (I/O, Otsu, bright-core, flip-rate,
metadados de rótulo/condição, estatística não-paramétrica e desenho de regiões).
Não altera o v11.1 nem o v11.2.

Contexto óptico (a razão de tudo): a óptica do microtubo de 1 mm retroiluminado é
**W/M** — o núcleo de gás aparece CLARO no centro (linha ~40), as interfaces/
meniscos aparecem ESCUROS fora do centro (linhas ~20/60), paredes saturam
(linhas 0–7). Logo "gás = pixel escuro" (α cordal do v11/v11.1) INVERTE o sinal.
As métricas aqui (núcleo CLARO) respeitam a física.
"""

import os
import re
import glob
from math import erf, sqrt

import numpy as np
from PIL import Image

# ── Caminhos e constantes compartilhadas (idênticas ao v11.2) ────────────────
HERE      = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.join(HERE, "..", "data_processed")

PANORAMA  = "__panorama_v11.png"          # panorama cordal (herdado do v11)
RAW_GLOB  = "ImgA*.jpg"                    # quadros brutos (80 linhas × 1016 px)

CORE_ROWS     = slice(20, 60)             # núcleo cordal (flip-rate e α de ref.)
INTERIOR_ROWS = slice(10, 70)            # janela interior (bright-core / areal)
CENTER_ROW    = 40
AXIAL_STEP    = 8                        # subamostragem axial (custo)
MAX_FRAMES    = 120                      # subamostra de quadros por ensaio

G_MAP = {3: 380, 4: 510, 5: 640, 6: 760, 7: 890}   # fluxo mássico [kg/(m²·s)]

EXCLUDE_RUNS = {"m03q80_001"}            # smooth annular (fora do escopo)
SAT_LIMIT    = 0.20                      # m03q40 sai da estatística (sat 34%)

RUN_RE = re.compile(r"^m0\dq\d+_\d+$")

# Nota: os marcadores de dispersed bubbles existem com DUAS grafias no disco —
# "__dispersed bubbles.txt" (espaço, m07q01_001) e "__dispersed_bubbles.txt"
# (sublinhado, m05q05_002). As duas são o MESMO padrão (o marcador manda; decisão
# do dono: m05q05_002 É dispersed bubbles). Ambas mapeadas para o grupo "bubble"
# ⇒ n=2. (O docstring do v11/v11.1 que chamava m05q05_002 de "plug" é obsoleto.)
GROUP_MAP = {
    "dispersed bubbles": "bubble",
    "dispersed_bubbles": "bubble",
    "plug":              "plug",
    "slug":              "slug",
    "churn":             "churn",
    "wavy annular":      "annular",
}
PURE_LABELS  = set(GROUP_MAP)
GROUP_ORDER  = ["bubble", "plug", "slug", "churn", "annular"]
GROUP_COLOR  = {"bubble": "#E69F00", "plug": "#0072B2", "slug": "#009E73",
                "churn": "#D55E00", "annular": "#CC79A7"}
GROUP_PRETTY = {"bubble": "dispersed bubbles", "plug": "plug", "slug": "slug",
                "churn": "churn", "annular": "annular"}


# ── Otsu (portado de v11.2 verbatim) ─────────────────────────────────────────

def manual_otsu(vals):
    """Limiar de Otsu (256 níveis) sobre valores 0–255. Classe escura = ``<= t``."""
    hist = np.histogram(vals, bins=256, range=(0, 256))[0].astype(float)
    total = hist.sum()
    sum_all = np.dot(np.arange(256), hist)
    wB = 0.0; sumB = 0.0; best = -1.0; thr = 0
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
            best = between; thr = t
    return thr


# ── Leitura de dados (portado de v11.2 verbatim) ─────────────────────────────

def load_panorama(name):
    p = os.path.join(DATA_ROOT, name, PANORAMA)
    return np.array(Image.open(p).convert("L"))


def load_raw_frames(name, max_frames=MAX_FRAMES):
    """Subamostra até ``max_frames`` quadros brutos (80×1016), uniformemente."""
    files = sorted(glob.glob(os.path.join(DATA_ROOT, name, RAW_GLOB)))
    if not files:
        return None
    if len(files) > max_frames:
        idx = np.linspace(0, len(files) - 1, max_frames).round().astype(int)
        files = [files[i] for i in idx]
    return np.stack([np.array(Image.open(f).convert("L")) for f in files])


# ── MÉTODO 1: central-bright-core (CONFIRMADO no v11.2) ──────────────────────

def _central_bright_run(col, thr):
    """Largura (fração) da maior sequência contígua de pixels CLAROS (``>thr``)
    que CONTÉM o centro da janela. 0 se o centro não for claro."""
    bright = col > thr
    c = len(col) // 2
    if not bright[c]:
        return 0.0
    lo = c
    while lo - 1 >= 0 and bright[lo - 1]:
        lo -= 1
    hi = c
    while hi + 1 < len(col) and bright[hi + 1]:
        hi += 1
    return (hi - lo + 1) / len(col)


def bright_core_void(frames, rows=INTERIOR_ROWS, use_midpoint=False):
    """Vazio 'bright-core': por quadro, limiar de Otsu (ou ponto-médio por coluna
    sem nível absoluto, se ``use_midpoint``); por coluna axial toma o núcleo
    CLARO central; mediana sobre colunas+quadros. Portado verbatim do v11.2."""
    vals = []
    for fr in frames:
        win = fr[rows, :].astype(float)
        if use_midpoint:
            thr_cols = 0.5 * (win.min(axis=0) + win.max(axis=0))
        else:
            thr = manual_otsu(win.ravel())
        for j in range(0, win.shape[1], AXIAL_STEP):
            col = win[:, j]
            t = thr_cols[j] if use_midpoint else thr
            vals.append(_central_bright_run(col, t))
    return float(np.median(vals)) if vals else np.nan


def axial_continuity(frames, rows=INTERIOR_ROWS):
    """Correlação de continuidade axial do núcleo claro (estrutura, não magnitude):
    anular (filme liso, núcleo contínuo) > churn (esfarrapado). Verbatim do v11.2."""
    corrs = []
    for fr in frames:
        win = fr[rows, :].astype(float)
        thr = manual_otsu(win.ravel())
        widths = np.array([_central_bright_run(win[:, j], thr)
                           for j in range(win.shape[1])])
        if widths.std() < 1e-9:
            continue
        corrs.append(np.corrcoef(widths[:-1], widths[1:])[0, 1])
    return float(np.nanmedian(corrs)) if corrs else np.nan


# ── α cordal de referência (eixo X do v11.1, para comparação) ────────────────

def alpha_chordal(panorama):
    """α cordal médio (gás = pixel escuro — INVERTIDO; só para comparação)."""
    core = panorama[CORE_ROWS, :].astype(float)
    thr = manual_otsu(core.ravel())
    aj = (core <= thr).mean(axis=0)
    return {"alpha_chordal": float(aj.mean()), "otsu_thr": int(thr),
            "core_mean": float(core.mean()),
            "sat_frac": float((core >= 255).mean())}


# ── Eixo Y comum: taxa de troca temporal (flip-rate, CONFIRMADO no v11.2) ─────

def temporal_flip_rate(panorama, thr=None):
    """Taxa de troca gás↔líquido por PIXEL no núcleo 2-D (espaço×tempo).
    Churn (caótico) troca ~2× mais que anular (núcleo estável). Verbatim do v11.2.
    Eixo ORTOGONAL de ESTABILIDADE — NÃO a magnitude do vazio."""
    core = panorama[CORE_ROWS, :].astype(float)
    t = manual_otsu(core.ravel()) if thr is None else thr
    gas = (core <= t).astype(int)
    flips = np.abs(np.diff(gas, axis=1)).sum()
    denom = gas.shape[0] * (gas.shape[1] - 1)
    rate = flips / denom
    alpha = float(gas.mean())
    null = 2.0 * alpha * (1.0 - alpha)
    rate_norm = rate / null if null > 1e-9 else np.nan
    return {"flip_rate": float(rate), "flip_rate_bnull": float(rate_norm),
            "alpha_core": alpha}


# ── Areal flat-field por baseline de líquido (variante B validada) ───────────

LIQ_BASELINE_PCT = 20       # percentil por pixel = fundo dominado por líquido
AREAL_GAS_RATIO  = 1.3      # Fn > 1.3 ⇒ excesso claro sobre o baseline = gás


def liquid_baseline(frames, rows=INTERIOR_ROWS, pct=LIQ_BASELINE_PCT):
    """Fundo por pixel = percentil ``pct`` no tempo (dominado por líquido; NÃO
    acompanha o núcleo de gás móvel). Base do flat-field que remove o reflexo W/M
    sem autonormalizar o sinal (a mediana temporal autonormaliza e destrói)."""
    win = frames[:, rows, :].astype(float)              # (T, nrows, 1016)
    return np.percentile(win, pct, axis=0)              # (nrows, 1016)


def areal_flatfield_void(frames, rows=INTERIOR_ROWS, pct=LIQ_BASELINE_PCT,
                         ratio=AREAL_GAS_RATIO):
    """Vazio AREAL 2-D sobre o quadro flat-field (variante B). Divide cada quadro
    pelo baseline de líquido, classifica gás = Fn>ratio (excesso claro), toma a
    fração 2-D de pixels-gás por quadro, mediana sobre quadros.

    NÃO é cordal — é a fração de ÁREA (schepperle2024). Recupera
    dispersed<plug<slug<churn, mas o annular cai ABAIXO do churn (bandas de
    interface escuras do anular persistem e derrotam o limiar de área clara).
    Proxy óptico RELATIVO não-calibrado. Ver MEMO_v11_3.md."""
    base = liquid_baseline(frames, rows, pct)
    base = np.where(base < 1.0, 1.0, base)              # evita divisão por ~0
    fracs = []
    for fr in frames:
        win = fr[rows, :].astype(float)
        fn = win / base
        fracs.append(float((fn > ratio).mean()))
    return float(np.median(fracs)) if fracs else np.nan


def radial_std_rawvsflat(frames, rows=INTERIOR_ROWS, pct=LIQ_BASELINE_PCT):
    """Prova da remoção do reflexo: std radial do col-mean por linha (média sobre
    quadros), bruto vs flat-field. O W/M bruto tem std radial grande (~60–70);
    após o flat-field colapsa ~100–300×."""
    win = frames[:, rows, :].astype(float)              # (T, nrows, 1016)
    raw_prof = win.mean(axis=2)                         # (T, nrows) col-mean/linha
    raw_std = float(raw_prof.mean(axis=0).std())        # std radial, média no tempo
    base = liquid_baseline(frames, rows, pct)
    base = np.where(base < 1.0, 1.0, base)
    fn = win / base[None, :, :]
    flat_prof = fn.mean(axis=2)
    flat_std = float(flat_prof.mean(axis=0).std())
    return raw_std, flat_std


# ── Metadados de rótulo/condição (portado de v11.2 verbatim) ─────────────────

def parse_m_q(run):
    mobj = re.match(r"m0?(\d+)q0?(\d+)", run)
    return int(mobj.group(1)), int(mobj.group(2))


def folder_flow_type(name):
    """Rótulo autoritativo pelo marcador ``__<rótulo>.txt`` (a fonte de verdade).
    m05q05_002 → 'dispersed bubbles' (marcador __dispersed_bubbles.txt)."""
    folder = os.path.join(DATA_ROOT, name)
    for p in sorted(glob.glob(os.path.join(folder, "*.txt"))):
        stem = os.path.splitext(os.path.basename(p))[0].lstrip("_")
        if stem != "dados":
            return stem
    return name


def list_scenarios():
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


# ── Estatística não-paramétrica (portado de v11.2 verbatim) ──────────────────

def mannwhitney_u_greater(x, y):
    """U e p unilateral (H1: x estocasticamente > y), aproximação normal com
    correção de empates. Sem scipy.stats."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    nx, ny = len(x), len(y)
    allv = np.concatenate([x, y])
    order = allv.argsort()
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(allv) + 1)
    _, inv, counts = np.unique(allv, return_inverse=True, return_counts=True)
    sums = np.zeros(len(counts)); np.add.at(sums, inv, ranks)
    avg = sums / counts
    ranks = avg[inv]
    Rx = ranks[:nx].sum()
    Ux = Rx - nx * (nx + 1) / 2.0
    mu = nx * ny / 2.0
    n = nx + ny
    tie = np.sum(counts ** 3 - counts)
    sigma = np.sqrt(nx * ny / 12.0 * ((n + 1) - tie / (n * (n - 1))))
    if sigma == 0:
        return Ux, 1.0
    z = (Ux - mu - 0.5) / sigma
    p = 0.5 * (1 - erf(z / sqrt(2)))
    return float(Ux), float(p)


def auc_greater(x, y):
    """AUC = P(x > y) + 0.5 P(x=y)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    gt = sum((xi > y).sum() for xi in x)
    eq = sum((xi == y).sum() for xi in x)
    return (gt + 0.5 * eq) / (len(x) * len(y))


# ── Utilidades de tabela/figura (portado de v11.2) ───────────────────────────

def md_table(fh, cols, rows_values):
    fh.write("| " + " | ".join(cols) + " |\n")
    fh.write("|" + "|".join(["---"] * len(cols)) + "|\n")
    for v in rows_values:
        fh.write("| " + " | ".join(str(x) for x in v) + " |\n")


def grp(rows, g, key, drop_sat=True):
    return [r[key] for r in rows if r["group"] == g
            and (not drop_sat or not r["dropped"]) and not np.isnan(r[key])]


def per_pattern_median(rows, key, drop_sat=True):
    out = []
    for g in GROUP_ORDER:
        vals = grp(rows, g, key, drop_sat)
        out.append((g, len(vals), (float(np.median(vals)) if vals else None)))
    return out


def draw_region(ax, pts, color, label, shift=(0.0, 0.0)):
    import matplotlib.patches as mpatches
    import matplotlib.patheffects as pe
    from scipy.spatial import ConvexHull, QhullError
    pts = np.asarray(pts, float)
    centroid = pts.mean(axis=0)
    if len(pts) >= 3:
        try:
            hull = ConvexHull(pts)
            v = pts[hull.vertices]
            ax.fill(v[:, 0], v[:, 1], color=color, alpha=0.15, zorder=1, lw=0)
            cl = np.vstack([v, v[0]])
            ax.plot(cl[:, 0], cl[:, 1], color=color, lw=1.2, alpha=0.6, zorder=2)
        except QhullError:
            pass
    else:
        half = pts.std(axis=0) if len(pts) >= 2 else np.array([0.0, 0.0])
        half = np.maximum(half, [0.01, 0.002]) + [0.008, 0.001]
        ax.add_patch(mpatches.Ellipse(centroid, 2 * half[0], 2 * half[1],
                     facecolor=color, alpha=0.15, edgecolor=color, lw=1.2,
                     zorder=1))
    t = ax.text(centroid[0] + shift[0], centroid[1] + shift[1], label,
                color=color, fontsize=11, fontweight="bold", ha="center",
                va="center", zorder=5)
    t.set_path_effects([pe.withStroke(linewidth=2.5, foreground="white")])


def make_void_flip_map(rows, xkey, xlabel, title, out_png, shift=None):
    """Mapa X=vazio × Y=flip-rate, regiões por padrão (rótulo CONHECIDO, não é
    clustering). dispersed bubbles (n=2) anotado explicitamente."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    shift = shift or {}
    fig, ax = plt.subplots(figsize=(8, 6))
    by = {}
    for r in rows:
        if r["dropped"] or np.isnan(r[xkey]):
            continue
        by.setdefault(r["group"], []).append(r)
    for g in GROUP_ORDER:
        pts = by.get(g)
        if not pts:
            continue
        coords = [(p[xkey], p["flip_rate"]) for p in pts]
        draw_region(ax, coords, GROUP_COLOR[g], GROUP_PRETTY[g],
                    shift.get(g, (0.0, 0.0)))
    for g in GROUP_ORDER:
        pts = by.get(g)
        if not pts:
            continue
        ax.scatter([p[xkey] for p in pts], [p["flip_rate"] for p in pts],
                   s=55, color=GROUP_COLOR[g], edgecolor="black", linewidth=0.5,
                   label=f"{GROUP_PRETTY[g]} (n={len(pts)})", zorder=4)
    # dispersed bubbles: anota n=2 (não entra em comparação estatística)
    bub = by.get("bubble")
    if bub:
        r = max(bub, key=lambda p: p[xkey])
        ax.annotate("n=2 (espúrio/instável)", (r[xkey], r["flip_rate"]),
                    textcoords="offset points", xytext=(8, 6), fontsize=7,
                    color=GROUP_COLOR["bubble"])
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Intermitência — taxa de troca temporal gás↔líquido  [por pixel]")
    ax.set_title(title)
    ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.4)
    ax.legend(loc="best", fontsize=8, framealpha=0.9)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"figura: {out_png}")
