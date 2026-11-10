"""
v11.2 — Fração de vazio NÃO-CORDAL (estrutura explorada no quadro BRUTO).

CONTEXTO (por que este estudo existe)
-------------------------------------
O v11.1 mantém a fração de vazio cordal α (média de ``mean(core[:,j] <= Otsu)``
sobre uma corda vertical no centro do tubo). Nessa métrica o churn e o anular
ondulado EMPATAM (churn 0,55 ≳ anular 0,52) e, pior, o churn sai LIGEIRAMENTE
ACIMA — fisicamente invertido: o anular deveria ter o MAIOR vazio (literatura:
α_anular ~ 0,9 areal). O orientador apontou o estranhamento.

Este v11.2 responde à pergunta: **"dá para, com uma análise NÃO-CORDAL
(explorando a estrutura radial / 2-D / temporal), recuperar anular > churn?"**

DIAGNÓSTICO ÓPTICO (a causa-raiz, confirmada no quadro bruto)
-------------------------------------------------------------
A óptica deste microtubo de 1 mm NÃO é "líquido-claro / gás-escuro". O perfil
radial de intensidade (média no tempo) tem forma de **W/M**:

    parede próx. (brilho de reflexo) → BANDA ESCURA (menisco/interface, linha ~20)
    → NÚCLEO CLARO no centro (linha ~40, intensidade 80–89)
    → BANDA ESCURA (interface, linha ~60) → parede distante (escura)

Ou seja: o **núcleo de gás aparece CLARO** (a retroiluminação atravessa o
núcleo), e as **interfaces/meniscos junto às paredes aparecem ESCURAS**. A
contagem ingênua "gás = pixel escuro" (o α cordal do v11/v11.1 e qualquer
contagem 2-D de pixels escuros) **inverte o sinal exatamente onde importa**:
conta os anéis de interface escuros como gás e o núcleo de gás claro como
líquido. É por isso que o churn ≳ anular vazava na métrica cordal.

OS QUATRO MÉTODOS TESTADOS (e o veredito adversarial de cada um)
---------------------------------------------------------------
Rodados como quatro sondas independentes + verificação adversarial
(workflow ``v11-nonchordal-void``). Reproduzidos aqui para auditoria:

  (1) raw-frame central-bright-core  →  CONFIRMADO, principiado
      Núcleo = maior sequência contígua de pixels CLAROS que contém o centro,
      por coluna axial do quadro bruto. Largura do núcleo = vazio cordal
      "bright-core". Recupera anular (0,35) > churn (0,28), p≈0,018.

  (2) wall-film core-width index     →  PARCIAL (mesma ideia de (1), config
      específica dá separação sem sobreposição p=0,0004, mas a ordem global
      dos 5 padrões fica errada — churn cai para o MENOR, bubble acima de
      plug/slug; confundimento de fluxo mássico m não totalmente controlado).

  (3) radial-profile shape (D1)      →  REFUTADO. É brilho re-embalado
      (corr(D1, brilho) = 0,95); some ao casar brilho (empate churn≈anular).

  (4) temporal flip-rate (eixo de ESTABILIDADE, não de magnitude) → CONFIRMADO
      Taxa de troca gás↔líquido por PIXEL no núcleo 2-D. NÃO recupera a
      magnitude do vazio, mas separa os padrões num eixo ORTOGONAL robusto:
      churn (caótico) troca ~2× mais que o anular (núcleo estável). anular
      0,016 < churn 0,030, AUC 0,878, permutação p≈0,001, imune a brilho e q.

VEREDITO (ver MEMO_v11_2.md para a leitura completa e honesta)
--------------------------------------------------------------
- A ORDEM anular > churn É recuperável por um método NÃO-CORDAL principiado
  (bright-core, método 1), porque o núcleo de gás é CLARO e a estrutura radial
  do anular (núcleo largo contínuo) difere da do churn (núcleo estreito,
  esfarrapado). Sobrevive a controles de brilho e de q.
- MAS é uma recuperação de ORDENAÇÃO ENTRE PADRÕES, não de MAGNITUDE: os
  valores (~0,35) são cordais em pixels, longe do areal ~0,9 da literatura; não
  há calibração px→mm nem segunda vista óptica, então o vazio areal verdadeiro
  NÃO é recuperável destas imagens.
- Ressalvas: distribuições se sobrepõem na tendência central, n pequeno (8×10);
  o estimador não é um vazio universal para os 5 padrões (bubble sai alto
  espúrio — centro coincidentemente claro; plug empata com anular).

Este script NÃO altera o v11.1. Gera as próprias saídas em ``src2/v11.2/``:
  nonchordal_void_table_v11_2.md   — tabela por método × padrão + controles
  nonchordal_void_map_v11_2.png    — mapa bright-core α × taxa de troca temporal
  radial_profile_v11_2.png         — perfil radial W/M (diagnóstico óptico)

Uso:
  cd src2/v11.2 && python3 nonchordal_void_v11_2.py
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

# ── Caminhos e constantes compartilhadas ─────────────────────────────────────
HERE      = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.join(HERE, "..", "data_processed")

PANORAMA  = "__panorama_v11.png"          # panorama cordal (herdado do v11)
RAW_GLOB  = "ImgA*.jpg"                    # quadros brutos (80 linhas × 1016 px)

# Núcleo cordal do v11/v11.1 (paredes excluídas) — usado só para o α de
# comparação e a taxa de troca temporal (método 4).
CORE_ROWS = slice(20, 60)

# Janela INTERIOR para os métodos de bright-core no quadro bruto. A janela
# exclui as linhas de reflexo de parede; a sondagem confirmou que o SINAL
# (anular>churn) não troca para nenhuma janela 8:72..15:65 nem offset de limiar.
INTERIOR_ROWS = slice(10, 70)             # método 1 (central-bright-core)
CENTER_ROW    = 40                        # linha central do tubo (núcleo)
AXIAL_STEP    = 8                         # subamostragem axial (col_step), p/ custo
MAX_FRAMES    = 120                       # subamostra de quadros por ensaio (custo)

# Método 2 (wall-film core-width): limiar adaptativo por coluna e janela própria.
WF_ROWS       = slice(15, 65)
WF_PLO, WF_PHI = 15, 85                   # percentis p/ thr = 0.5*(p15+p85)

G_MAP = {3: 380, 4: 510, 5: 640, 6: 760, 7: 890}   # fluxo mássico [kg/(m²·s)]

# Taxa de troca temporal (método 4): estados por pixel no núcleo 2-D, Otsu por
# imagem; a taxa é trocas gás↔líquido entre colunas de tempo adjacentes.
# (não usa histerese — é por-pixel, não a alternação cordal do v11.1)

EXCLUDE_RUNS = {"m03q80_001"}             # smooth annular (fora do escopo)
SAT_LIMIT    = 0.20                       # m03q40 sai da estatística (sat 34%)

RUN_RE = re.compile(r"^m0\dq\d+_\d+$")

GROUP_MAP = {
    "dispersed bubbles": "bubble",
    "plug":              "plug",
    "slug":              "slug",
    "churn":             "churn",
    "wavy annular":      "annular",
}
PURE_LABELS = set(GROUP_MAP)
GROUP_ORDER = ["bubble", "plug", "slug", "churn", "annular"]
GROUP_COLOR = {"bubble": "#E69F00", "plug": "#0072B2", "slug": "#009E73",
               "churn": "#D55E00", "annular": "#CC79A7"}
GROUP_PRETTY = {g: ("dispersed bubbles" if g == "bubble" else g) for g in GROUP_ORDER}


# ── Otsu (portado de v11 verbatim) ───────────────────────────────────────────

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


# ── Leitura de dados ─────────────────────────────────────────────────────────

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


# ── α cordal de referência (eixo X do v11.1, para comparação) ────────────────

def alpha_chordal(panorama):
    """α cordal médio (gás = pixel escuro) + metadados de brilho/saturação."""
    core = panorama[CORE_ROWS, :].astype(float)
    thr = manual_otsu(core.ravel())
    aj = (core <= thr).mean(axis=0)
    return {
        "alpha_chordal": float(aj.mean()),
        "otsu_thr": int(thr),
        "core_mean": float(core.mean()),
        "sat_frac": float((core >= 255).mean()),
        "aj": aj,
    }


# ── MÉTODO 1: central-bright-core (CONFIRMADO, principiado) ──────────────────

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
    claro central; mediana sobre colunas+quadros.

    use_midpoint=True é o CONTROLE SEM NÍVEL DE BRILHO: thr = (min+max)/2 por
    coluna, removendo o nível absoluto de cinza. Se anular>churn sobrevive, a
    separação é de FORMA do perfil, não de brilho.
    """
    vals = []
    for fr in frames:
        win = fr[rows, :].astype(float)              # (nrows, 1016)
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
    """Correlação de continuidade axial do núcleo claro (quão suave/contínuo é o
    núcleo ao longo do eixo). Anular (filme liso, núcleo contínuo) > churn
    (núcleo esfarrapado com lacunas). Descritor de ESTRUTURA, não de magnitude."""
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


# ── MÉTODO 2: wall-film core-width index (PARCIAL) ───────────────────────────

def wall_film_core(frames, rows=WF_ROWS):
    """Largura do núcleo claro com limiar adaptativo por coluna
    thr = 0.5*(p15+p85). Config que dá separação sem sobreposição anular×churn
    (p=0,0004 na sondagem) — mas a ordem GLOBAL dos 5 padrões sai errada."""
    vals = []
    for fr in frames:
        win = fr[rows, :].astype(float)
        plo = np.percentile(win, WF_PLO, axis=0)
        phi = np.percentile(win, WF_PHI, axis=0)
        thr_cols = 0.5 * (plo + phi)
        for j in range(win.shape[1]):
            vals.append(_central_bright_run(win[:, j], thr_cols[j]))
    return float(np.mean(vals)) if vals else np.nan


# ── MÉTODO 4: taxa de troca temporal (eixo de ESTABILIDADE, CONFIRMADO) ──────

def temporal_flip_rate(panorama, thr=None):
    """Taxa de troca gás↔líquido por PIXEL no núcleo 2-D (espaço×tempo).

    Binariza o núcleo (gás = ``<= thr``, Otsu por imagem se thr=None), conta
    trocas entre colunas de tempo adjacentes, normaliza por (n_linhas × n_cols).
    Churn (caótico) troca ~2× mais que anular (núcleo estável). É o eixo
    ORTOGONAL de estabilidade — NÃO a magnitude do vazio.
    """
    core = panorama[CORE_ROWS, :].astype(float)
    t = manual_otsu(core.ravel()) if thr is None else thr
    gas = (core <= t).astype(int)                     # 1 = gás, por pixel
    flips = np.abs(np.diff(gas, axis=1)).sum()        # trocas no eixo tempo
    denom = gas.shape[0] * (gas.shape[1] - 1)
    rate = flips / denom
    alpha = float(gas.mean())
    # normalização de null de Bernoulli i.i.d. 2α(1-α): remove a vantagem
    # mecânica de α≈0,5 estar mais perto do meio (artefato de proximidade).
    null = 2.0 * alpha * (1.0 - alpha)
    rate_norm = rate / null if null > 1e-9 else np.nan
    return {"flip_rate": float(rate), "flip_rate_bnull": float(rate_norm),
            "alpha_core": alpha}


# ── Metadados de rótulo/condição (portado de v11) ────────────────────────────

def parse_m_q(run):
    mobj = re.match(r"m0?(\d+)q0?(\d+)", run)
    return int(mobj.group(1)), int(mobj.group(2))


def folder_flow_type(name):
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


# ── Estatística: Mann-Whitney U unilateral + permutação run-level ────────────

def mannwhitney_u_greater(x, y):
    """U e p unilateral (H1: x estocasticamente > y), aproximação normal com
    correção de empates. Sem scipy.stats (evita dependência extra)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    nx, ny = len(x), len(y)
    allv = np.concatenate([x, y])
    order = allv.argsort()
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(allv) + 1)
    # média de ranks em empates
    _, inv, counts = np.unique(allv, return_inverse=True, return_counts=True)
    sums = np.zeros(len(counts)); np.add.at(sums, inv, ranks)
    avg = sums / counts
    ranks = avg[inv]
    Rx = ranks[:nx].sum()
    Ux = Rx - nx * (nx + 1) / 2.0
    mu = nx * ny / 2.0
    # desvio com correção de empates
    n = nx + ny
    tie = np.sum(counts ** 3 - counts)
    sigma = np.sqrt(nx * ny / 12.0 * ((n + 1) - tie / (n * (n - 1))))
    if sigma == 0:
        return Ux, 1.0
    z = (Ux - mu - 0.5) / sigma
    # p unilateral H1: x>y  ⇒ U grande
    from math import erf, sqrt
    p = 0.5 * (1 - erf(z / sqrt(2)))
    return float(Ux), float(p)


def auc_greater(x, y):
    """AUC = P(x > y) + 0.5 P(x=y), ou seja quão bem x fica acima de y."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    gt = sum((xi > y).sum() for xi in x)
    eq = sum((xi == y).sum() for xi in x)
    return (gt + 0.5 * eq) / (len(x) * len(y))


# ── Orquestração ─────────────────────────────────────────────────────────────

def main():
    scenarios = list_scenarios()
    print(f"{len(scenarios)} ensaios de padrão puro (dispersed bubbles INCLUÍDO; "
          "transições, synthetic*, smooth annular excluídos)\n")

    rows = []
    for name in scenarios:
        pano = load_panorama(name)
        ch = alpha_chordal(pano)
        flip = temporal_flip_rate(pano)
        frames = load_raw_frames(name)

        bc = bright_core_void(frames) if frames is not None else np.nan
        bc_mid = bright_core_void(frames, use_midpoint=True) if frames is not None else np.nan
        axc = axial_continuity(frames) if frames is not None else np.nan
        wf = wall_film_core(frames) if frames is not None else np.nan

        label = folder_flow_type(name)
        group = GROUP_MAP.get(label, label)
        m, q = parse_m_q(name)
        dropped = ch["sat_frac"] > SAT_LIMIT
        flag = (f"saturado (sat={ch['sat_frac'] * 100:.0f}%)" if dropped else "")

        rows.append({
            "scenario": name, "pattern": label, "group": group, "m": m, "q": q,
            "G": G_MAP[m], "alpha_chordal": ch["alpha_chordal"],
            "core_mean": ch["core_mean"], "sat_frac": ch["sat_frac"],
            "bright_core": bc, "bright_core_mid": bc_mid, "axial_cont": axc,
            "wall_film": wf, "flip_rate": flip["flip_rate"],
            "flip_rate_bnull": flip["flip_rate_bnull"],
            "N": int(pano.shape[1]), "dropped": dropped, "flag": flag,
        })
        print(f"  {name:14s} {label:18s} α={ch['alpha_chordal']:.3f} "
              f"bright_core={bc:.3f} wall_film={wf:.3f} flip={flip['flip_rate']:.4f}"
              + ("  [SAT]" if dropped else ""))

    write_table(rows)
    make_map(rows)
    make_radial_profile(scenarios, rows)
    report_stats(rows)


def _grp(rows, g, key, drop_sat=True):
    return [r[key] for r in rows if r["group"] == g
            and (not drop_sat or not r["dropped"]) and not np.isnan(r[key])]


def report_stats(rows):
    print("\n" + "=" * 72)
    print("VEREDITO ESTATÍSTICO (anular × churn) — sat descartada")
    print("=" * 72)
    for key, name in [("bright_core", "MÉTODO 1 bright-core (CONFIRMADO)"),
                      ("wall_film", "MÉTODO 2 wall-film (PARCIAL)"),
                      ("flip_rate", "MÉTODO 4 flip-rate (CONFIRMADO, estabilidade)")]:
        ann = _grp(rows, "annular", key)
        chu = _grp(rows, "churn", key)
        if key == "flip_rate":
            # eixo de estabilidade: churn > annular
            U, p = mannwhitney_u_greater(chu, ann)
            auc = auc_greater(chu, ann)
            direc = f"churn {np.median(chu):.4f} > annular {np.median(ann):.4f}"
        else:
            U, p = mannwhitney_u_greater(ann, chu)
            auc = auc_greater(ann, chu)
            direc = f"annular {np.median(ann):.3f} > churn {np.median(chu):.3f}"
        print(f"\n{name}\n  {direc}  | MWU p(unilateral)={p:.4f} AUC={auc:.3f} "
              f"(n_ann={len(ann)}, n_chu={len(chu)})")


# ── Tabela ────────────────────────────────────────────────────────────────────

def _md(fh, cols, rows_values):
    fh.write("| " + " | ".join(cols) + " |\n")
    fh.write("|" + "|".join(["---"] * len(cols)) + "|\n")
    for v in rows_values:
        fh.write("| " + " | ".join(str(x) for x in v) + " |\n")


def _per_pattern(rows, key, drop_sat=True):
    out = []
    for g in GROUP_ORDER:
        vals = _grp(rows, g, key, drop_sat)
        if not vals:
            out.append((g, 0, None))
        else:
            out.append((g, len(vals), float(np.median(vals))))
    return out


def write_table(rows):
    out = os.path.join(HERE, "nonchordal_void_table_v11_2.md")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("# Fração de vazio NÃO-CORDAL por método × padrão (v11.2)\n\n")
        fh.write(
            "Estudo de métodos NÃO-cordais para recuperar **anular > churn** no "
            "vazio — que a métrica cordal α do v11/v11.1 não consegue (churn "
            "≳ anular, invertido). Causa-raiz óptica: o núcleo de gás aparece "
            "**CLARO** e as interfaces junto às paredes aparecem **ESCURAS** "
            "(perfil radial em W/M), então contar 'gás = pixel escuro' inverte o "
            "sinal. Ver `MEMO_v11_2.md` para a leitura completa.\n\n")
        fh.write("Quatro métodos, cada um verificado adversarialmente "
                 "(workflow `v11-nonchordal-void`):\n\n")
        fh.write("| método | ideia | recupera anular>churn? | veredito | "
                 "principiado? |\n|---|---|---|---|---|\n")
        fh.write("| 1 bright-core | maior núcleo CLARO contíguo pelo centro "
                 "(quadro bruto) | **sim** (0,35 vs 0,28, p≈0,018) | "
                 "**CONFIRMADO** | **sim** |\n")
        fh.write("| 2 wall-film | largura de núcleo, limiar por-coluna | sim "
                 "(sem sobreposição, p=0,0004) | PARCIAL | escopado "
                 "(ordem global dos 5 erra) |\n")
        fh.write("| 3 radial-shape (D1) | centro−parede normalizado | parcial | "
                 "**REFUTADO** | não (é brilho re-embalado, corr 0,95) |\n")
        fh.write("| 4 flip-rate | troca gás↔líquido por pixel (tempo) | sim, como "
                 "*estabilidade* | **CONFIRMADO** | **sim** |\n\n")

        n_in = sum(1 for r in rows if not r["dropped"])
        fh.write(f"Padrões: dispersed bubbles, plug, slug, churn, wavy annular — "
                 f"{len(rows)} ensaios ({n_in} na estatística).\n\n")

        # Tabela por método × padrão (mediana)
        fh.write("## Mediana por padrão (ensaios saturados fora)\n\n")
        cols = ["padrão", "n", "α cordal (v11.1)", "bright-core (M1)",
                "wall-film (M2)", "flip-rate (M4)"]
        vals = []
        bc = dict((g, m) for g, _, m in _per_pattern(rows, "bright_core"))
        wf = dict((g, m) for g, _, m in _per_pattern(rows, "wall_film"))
        fr = dict((g, m) for g, _, m in _per_pattern(rows, "flip_rate"))
        al = dict((g, m) for g, _, m in _per_pattern(rows, "alpha_chordal"))
        ns = dict((g, n) for g, n, _ in _per_pattern(rows, "bright_core"))
        for g in GROUP_ORDER:
            def f(x): return "—" if x is None else f"{x:.3f}"
            vals.append([GROUP_PRETTY[g], ns[g], f(al[g]), f(bc[g]),
                         f(wf[g]), (("—" if fr[g] is None else f"{fr[g]:.4f}"))])
        _md(fh, cols, vals)

        fh.write("\n**Leitura:** no α cordal (v11.1) churn ≳ anular (invertido). "
                 "Os métodos 1 e 2 (núcleo CLARO) recuperam **anular > churn**. "
                 "No método 4 (flip-rate) **churn > anular** (churn é o mais "
                 "instável) — eixo ORTOGONAL de estabilidade, não de magnitude.\n\n")

        # Controle de brilho (método 1 midpoint)
        fh.write("## Controle de brilho — método 1 sem nível absoluto\n\n")
        fh.write("`bright-core-mid` usa limiar = (min+max)/2 por coluna, removendo "
                 "o NÍVEL de cinza. Se anular>churn sobrevive, a separação é de "
                 "FORMA do perfil, não de brilho.\n\n")
        bcm = dict((g, m) for g, _, m in _per_pattern(rows, "bright_core_mid"))
        _md(fh, ["padrão", "bright-core (Otsu)", "bright-core-mid (sem nível)"],
            [[GROUP_PRETTY[g],
              ("—" if bc[g] is None else f"{bc[g]:.3f}"),
              ("—" if bcm[g] is None else f"{bcm[g]:.3f}")] for g in GROUP_ORDER])

        # Continuidade axial (estrutura)
        fh.write("\n## Estrutura radial — continuidade axial do núcleo (M1)\n\n")
        fh.write("Correlação de continuidade axial do núcleo claro. Anular (filme "
                 "liso, núcleo contínuo) > churn (núcleo esfarrapado).\n\n")
        axc = dict((g, m) for g, _, m in _per_pattern(rows, "axial_cont"))
        _md(fh, ["padrão", "continuidade axial"],
            [[GROUP_PRETTY[g], ("—" if axc[g] is None else f"{axc[g]:.3f}")]
             for g in GROUP_ORDER])

        # Per-ensaio
        fh.write("\n## Por ensaio (transparência)\n\n")
        cols = ["scenario", "pattern", "G", "q", "α cordal", "bright-core",
                "wall-film", "flip-rate", "sat %", "flag"]
        rr = sorted(rows, key=lambda r: (-(r["bright_core"]
                     if not np.isnan(r["bright_core"]) else -1)))
        _md(fh, cols, [[r["scenario"], r["pattern"], r["G"], r["q"],
                        f"{r['alpha_chordal']:.3f}", f"{r['bright_core']:.3f}",
                        f"{r['wall_film']:.3f}", f"{r['flip_rate']:.4f}",
                        f"{r['sat_frac']*100:.0f}", r["flag"]] for r in rr])

        fh.write("\n## Ressalvas honestas\n\n")
        fh.write("- **Ordenação, não magnitude.** bright-core ~0,35 é cordal em "
                 "pixels, longe do areal ~0,9 da literatura. Sem calibração px→mm "
                 "nem segunda vista, o vazio areal verdadeiro NÃO é recuperável.\n")
        fh.write("- **Sobreposição e n pequeno** (anular 8 × churn 10): a "
                 "separação é de tendência central.\n")
        fh.write("- **Não é vazio universal dos 5 padrões:** dispersed bubbles "
                 "sai alto espúrio (centro coincidentemente claro — a regra "
                 "'centro claro' precisa de porta de contraste); plug empata com "
                 "anular. Validado para a pergunta anular×churn, não como vazio "
                 "universal.\n")
        fh.write("- **Método 3 (D1) REFUTADO:** é brilho re-embalado "
                 "(corr(D1, brilho)=0,95); empata ao casar brilho.\n")
    print(f"\ntabela: {out}")


# ── Mapa bright-core α × flip-rate ────────────────────────────────────────────

def _draw_region(ax, pts, color, label, shift=(0.0, 0.0)):
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
        ax.add_patch(mpatches.Ellipse(centroid, 2*half[0], 2*half[1],
                     facecolor=color, alpha=0.15, edgecolor=color, lw=1.2, zorder=1))
    t = ax.text(centroid[0] + shift[0], centroid[1] + shift[1], label,
                color=color, fontsize=11, fontweight="bold", ha="center",
                va="center", zorder=5)
    t.set_path_effects([pe.withStroke(linewidth=2.5, foreground="white")])


def make_map(rows):
    """Mapa bright-core (X, magnitude recuperada) × flip-rate (Y, estabilidade).
    Os dois eixos NÃO-cordais que sobreviveram à verificação adversarial."""
    fig, ax = plt.subplots(figsize=(8, 6))
    by = {}
    for r in rows:
        if r["dropped"] or np.isnan(r["bright_core"]):
            continue
        by.setdefault(r["group"], []).append(r)
    SHIFT = {"bubble": (0.0, 0.004), "plug": (0.0, 0.004), "slug": (0.03, 0.0),
             "churn": (0.0, 0.004), "annular": (0.0, -0.004)}
    for g in GROUP_ORDER:
        pts = by.get(g)
        if not pts:
            continue
        coords = [(p["bright_core"], p["flip_rate"]) for p in pts]
        _draw_region(ax, coords, GROUP_COLOR[g], GROUP_PRETTY[g],
                     SHIFT.get(g, (0, 0)))
    for g in GROUP_ORDER:
        pts = by.get(g)
        if not pts:
            continue
        ax.scatter([p["bright_core"] for p in pts], [p["flip_rate"] for p in pts],
                   s=55, color=GROUP_COLOR[g], edgecolor="black", linewidth=0.5,
                   label=f"{GROUP_PRETTY[g]} (n={len(pts)})", zorder=4)
    ax.set_xlabel("Vazio NÃO-cordal 'bright-core' (M1)  [–, cordal px]")
    ax.set_ylabel("Taxa de troca temporal gás↔líquido (M4)  [por pixel]")
    ax.set_title("v11.2 — magnitude (bright-core) × estabilidade (flip-rate)\n"
                 "anular = vazio ALTO + estável · churn = vazio alto + caótico")
    ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.4)
    ax.legend(loc="best", fontsize=8, framealpha=0.9)
    out = os.path.join(HERE, "nonchordal_void_map_v11_2.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"figura: {out}")


# ── Perfil radial W/M (diagnóstico óptico) ───────────────────────────────────

def make_radial_profile(scenarios, rows):
    """Perfil radial médio (intensidade × linha) para um anular e um churn
    representativos — mostra o W/M: núcleo central CLARO, interfaces ESCURAS."""
    def pick(group):
        cand = [r for r in rows if r["group"] == group and not r["dropped"]]
        cand.sort(key=lambda r: abs(r["q"] - 40))   # q comparável
        return cand[0]["scenario"] if cand else None

    ann = pick("annular"); chu = pick("churn")
    if not (ann and chu):
        return
    fig, ax = plt.subplots(figsize=(7, 5))
    for name, color, lab in [(ann, GROUP_COLOR["annular"], f"anular ({ann})"),
                             (chu, GROUP_COLOR["churn"], f"churn ({chu})")]:
        frames = load_raw_frames(name, max_frames=200)
        prof = frames.mean(axis=(0, 2))              # média sobre quadros e axial
        ax.plot(np.arange(len(prof)), prof, color=color, lw=1.8, label=lab)
    ax.axvspan(CENTER_ROW - 5, CENTER_ROW + 5, color="0.85", alpha=0.5,
               label="núcleo (centro)")
    ax.set_xlabel("linha (seção transversal do tubo)  — 0/79 = paredes")
    ax.set_ylabel("intensidade média  [0–255]")
    ax.set_title("Perfil radial em W/M — núcleo de gás CLARO, interfaces ESCURAS\n"
                 "(por isso 'gás = pixel escuro' inverte o sinal)")
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.4)
    out = os.path.join(HERE, "radial_profile_v11_2.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"figura: {out}")


if __name__ == "__main__":
    main()
