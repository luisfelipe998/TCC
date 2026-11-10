"""
run_stats.py — estatística dos QUATRO padrões (churn, anular ondulado, plug, slug).

Evolução de `../estatistica2/run_stats.py`. Mesma filosofia: TUDO sai de UMA matriz
D = 1 − Pearson(descritor) (26×26). Sem scikit-bio, sem ML, sem PERMANOVA/PERMDISP.

O que muda em relação à estatistica2 (2 padrões → 4 padrões):
  A0  Figura 17 generalizada: heatmap 26×26 da matriz de distância, ordenado por padrão
      (churn, anular, plug, slug) com linhas separadoras e rótulos de bloco → v8-estat-dmatrix.
  A1  ANOSIM R GLOBAL multigrupo (R de Clarke, 4 grupos; within = união dos pares intra de todos
      os grupos; between = pares entre grupos distintos). p por:
        · enumeração EXATA q-estratificada (18 000 arranjos; mal condicionada — ver A4),
        · Monte Carlo m-estratificado (melhor condicionada; nenhum estrato puro),
        · Monte Carlo livre (contexto).
      Permutação é da ETIQUETA de 4 vias dentro de cada estrato (multiset por estrato fixo).
  A1b MATRIZ de R par-a-par entre os 6 pares de padrões (quais padrões se separam de quais).
      plug/slug são minúsculos (4/3) → R par-a-par com <30 distâncias leva ressalva de n pequeno.
  A2  ECDF intra vs inter + correlações pairwise corr(D,|Δq|), corr(D,|Δm|), corr(D,mesmo padrão).
  A4  diagnóstico de confound: pureza dos estratos q vs m (por que m é a referência).
  A3' robustez piso50 vs piso10 (sonda), idêntico em espírito à estatistica2.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from itertools import combinations, product
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))

# ordem de exibição dos padrões na Figura 17 (regime: intermitente lento → rápido disperso)
PATTERN_ORDER = ["churn", "wavy annular", "plug", "slug"]
PATTERN_PT = {"churn": "churn", "wavy annular": "anular ondulado", "plug": "plug", "slug": "slug"}
PATTERN_ABBR = {"churn": "CH", "wavy annular": "AO", "plug": "PL", "slug": "SL"}

RNG_SEED = 0
N_MC = 50000          # sorteios Monte Carlo (m-estratificado e livre); p_min = 1/(N_MC+1)


# ── matriz de distância ──────────────────────────────────────────────────────

def dist_matrix(Z):
    """D = 1 - Pearson entre linhas (runs). Z: (n_runs, n_scales)."""
    C = np.corrcoef(Z)
    return 1.0 - C


# ── ANOSIM R multigrupo (R de Clarke) ────────────────────────────────────────

def anosim_R(D, groups):
    """R de Clarke, VÁLIDO PARA QUALQUER nº de grupos.

    groups: vetor inteiro de rótulo por run. within = todos os pares cujos dois runs têm o
    MESMO rótulo (de qualquer grupo); between = pares de rótulos distintos. Ranks sobre TODAS
    as M = C(n,2) distâncias. R = (r_between − r_within)/(M/2).
    """
    n = D.shape[0]
    iu = np.triu_indices(n, 1)
    d = D[iu]
    same = groups[iu[0]] == groups[iu[1]]
    ranks = d.argsort().argsort().astype(float) + 1.0     # rank 1..M (ties: ordem estável)
    r_between = ranks[~same].mean()
    r_within  = ranks[same].mean()
    M = len(d)
    return (r_between - r_within) / (M / 2.0)


# ── permutação estratificada da etiqueta de 4 vias ───────────────────────────

def _strata_index(strata):
    """dict estrato -> lista de índices de run."""
    idx = {}
    for i, s in enumerate(strata):
        idx.setdefault(s, []).append(i)
    return idx


def stratified_arrangements_exact(strata, labels):
    """Enumera TODOS os arranjos da etiqueta (4 vias) permutando dentro de cada estrato.

    Mantém fixo o multiset de rótulos de cada estrato. Produto cartesiano das permutações
    distintas por estrato. Usável só quando o total é pequeno (q-estratificado: 18 000).
    """
    idx_by_s = _strata_index(strata)
    per_s_options = []
    for s in sorted(idx_by_s, key=lambda x: (float("inf") if x is None else x)):
        members = idx_by_s[s]
        lab = [labels[i] for i in members]
        seen = set()
        opts = []
        # permutações distintas do multiset de rótulos deste estrato
        for perm in set(_distinct_perms(tuple(lab))):
            assign = {members[k]: perm[k] for k in range(len(members))}
            opts.append(assign)
        per_s_options.append(opts)

    arrangements = []
    for combo in product(*per_s_options):
        vec = np.empty(len(strata), dtype=object)
        for assign in combo:
            for i, c in assign.items():
                vec[i] = c
        arrangements.append(vec)
    return arrangements


def _distinct_perms(seq):
    """Permutações distintas de um multiset (evita explosão por rótulos repetidos)."""
    from itertools import permutations
    return set(permutations(seq))


def count_arrangements(strata, labels):
    """Multinomial por estrato (produto) — nº total de arranjos estratificados."""
    import math
    idx_by_s = _strata_index(strata)
    tot = 1
    for s, members in idx_by_s.items():
        c = Counter(labels[i] for i in members)
        num = math.factorial(len(members))
        for v in c.values():
            num //= math.factorial(v)
        tot *= num
    return tot


def mc_stratified_labels(strata, labels, rng):
    """Um sorteio: embaralha os rótulos DENTRO de cada estrato (multiset preservado)."""
    out = np.array(labels, dtype=object).copy()
    idx_by_s = _strata_index(strata)
    for members in idx_by_s.values():
        perm = rng.permutation(len(members))
        orig = [labels[i] for i in members]
        for k, i in enumerate(members):
            out[i] = orig[perm[k]]
    return out


def to_int(labels):
    """Mapeia rótulos (strings) para inteiros para o anosim_R."""
    uniq = {p: k for k, p in enumerate(PATTERN_ORDER)}
    return np.array([uniq[str(x)] for x in labels], dtype=int)


# ── A0: Figura 17 (matriz de distância 26×26 ordenada por padrão) ─────────────

def figure_dmatrix(D, labels, out_png):
    order = []
    counts = []
    for p in PATTERN_ORDER:
        members = [i for i in range(len(labels)) if labels[i] == p]
        order.extend(members)
        counts.append(len(members))
    order = np.array(order)
    Dord = D[np.ix_(order, order)]
    n = len(order)

    fig, ax = plt.subplots(figsize=(8.4, 7.2))
    im = ax.imshow(Dord, cmap="viridis_r", origin="upper")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(r"distância $1 - \rho$ de Pearson")

    # linhas separadoras entre blocos de padrão
    boundaries = np.cumsum(counts)[:-1]
    for b in boundaries:
        ax.axhline(b - 0.5, color="white", lw=1.6)
        ax.axvline(b - 0.5, color="white", lw=1.6)

    # rótulos de bloco acima da matriz (eixo secundário no topo, sem colidir com o título)
    starts = np.concatenate([[0], np.cumsum(counts)])
    secx = ax.secondary_xaxis("top")
    centers = [(starts[k] + starts[k + 1] - 1) / 2.0 for k in range(len(PATTERN_ORDER))]
    secx.set_xticks(centers)
    secx.set_xticklabels(["%s (%d)" % (PATTERN_PT[p], counts[k])
                          for k, p in enumerate(PATTERN_ORDER)], fontsize=8)
    secx.tick_params(length=0)

    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(range(1, n + 1), fontsize=6)
    ax.set_yticklabels(range(1, n + 1), fontsize=6)
    ax.set_xlabel("ensaio"); ax.set_ylabel("ensaio")
    ax.set_title("Matriz de distância espectral entre os 26 ensaios\n"
                 "(descritor WPS no atraso nulo, branqueado)", fontsize=10, pad=28)
    fig.savefig(out_png, dpi=150, bbox_inches="tight"); plt.close(fig)
    return order, counts


# ── A1: ANOSIM R global multigrupo + p ────────────────────────────────────────

def a1_global(Z, labels, q, m, fh):
    D = dist_matrix(Z)
    g = to_int(labels)
    R_obs = anosim_R(D, g)

    # p EXATO q-estratificado (18 000) — enumeração
    arr_q = stratified_arrangements_exact(q, labels)
    Rs_q = np.array([anosim_R(D, to_int(a)) for a in arr_q])
    p_q = float(np.mean(Rs_q >= R_obs))
    rank_q = int(np.sum(Rs_q > R_obs))

    # p Monte Carlo m-estratificado (108.9M arranjos → não enumerável)
    rng = np.random.default_rng(RNG_SEED)
    Rs_m = np.empty(N_MC)
    for t in range(N_MC):
        Rs_m[t] = anosim_R(D, to_int(mc_stratified_labels(m, labels, rng)))
    p_m = (1 + int(np.sum(Rs_m >= R_obs))) / (1 + N_MC)   # p MC inclui o observado
    n_m_arr = count_arrangements(m, labels)

    # p livre (contexto): embaralho irrestrito do rótulo de 4 vias
    Rs_free = np.empty(N_MC)
    base = np.array(labels, dtype=object)
    for t in range(N_MC):
        Rs_free[t] = anosim_R(D, to_int(rng.permutation(base)))
    p_free = (1 + int(np.sum(Rs_free >= R_obs))) / (1 + N_MC)

    # percentil 95 do nulo m-estratificado (régua do acaso)
    p95_m = float(np.quantile(Rs_m, 0.95))

    out = (
        "[A1] ANOSIM R GLOBAL (4 grupos: churn, anular, plug, slug)\n"
        "  R de Clarke (observado)              = %.4f\n"
        "  p q-estratificado (EXATO, %d arr.)  = %.4f  (R_max=%.4f; rank=%d; p_min=%.2e)\n"
        "  p m-estratificado (MC %d; %s arr.)  = %.4f  (R_max=%.4f; p95_nulo=%.3f)\n"
        "  p não-estratificado (MC %d)          = %.4f  (R_max=%.4f)\n"
        "  -> q-estrat é MAL condicionada (estratos quase puros — ver A4); a m-estrat é a\n"
        "     referência (nenhum estrato puro). Efeito-tamanho DESCRITIVO, não teste forte.\n"
        % (R_obs, len(arr_q), p_q, Rs_q.max(), rank_q, 1.0 / len(arr_q),
           N_MC, "{:,}".format(n_m_arr), p_m, Rs_m.max(), p95_m,
           N_MC, p_free, Rs_free.max())
    )
    print(out, end=""); fh.write(out)
    return R_obs, p_q, p_m, p_free, p95_m


# ── A1b: matriz R par-a-par entre padrões ─────────────────────────────────────

def a1b_pairwise(Z, labels, fh):
    D = dist_matrix(Z)
    out = ["[A1b] R de Clarke PAR-A-PAR (submatriz de cada par de padrões)\n"
           "      quais padrões se separam de quais; n_pares = C(n_i+n_j, 2)\n"]
    counts = {p: int((labels == p).sum()) for p in PATTERN_ORDER}
    rows = []
    for a, b in combinations(PATTERN_ORDER, 2):
        sel = np.where((labels == a) | (labels == b))[0]
        Dsub = D[np.ix_(sel, sel)]
        g = (labels[sel] == a).astype(int)
        R = anosim_R(Dsub, g)
        npairs = len(sel) * (len(sel) - 1) // 2
        caveat = "  (n pequeno)" if npairs < 30 else ""
        rows.append((a, b, R, npairs, caveat))
        out.append("  %-14s × %-14s  R=%+.3f  (n_pares=%d)%s\n"
                   % (PATTERN_PT[a], PATTERN_PT[b], R, npairs, caveat))
    out.append("  NOTA: R>0 ⇒ os dois padrões separam; R≈0 ⇒ indistinguíveis pelo descritor.\n"
               "        Pares com plug(4)/slug(3) repousam em <30 distâncias → ler com cautela\n"
               "        (R par-a-par instável para grupos minúsculos; sem correção de múltiplas\n"
               "         comparações — 6 pares são exploratórios, não família de testes).\n")
    txt = "".join(out)
    print(txt, end=""); fh.write(txt)
    return rows


# ── A2: ECDF intra vs inter + correlações pairwise ───────────────────────────

def a2(Z, labels, q, m, fh):
    D = dist_matrix(Z)
    n = D.shape[0]
    intra, inter = [], []
    dq_all, dm_all, D_all, same_all = [], [], [], []
    for i, j in combinations(range(n), 2):
        d = D[i, j]
        D_all.append(d); dq_all.append(abs(q[i] - q[j])); dm_all.append(abs(m[i] - m[j]))
        same = labels[i] == labels[j]
        same_all.append(same)
        (intra if same else inter).append(d)
    intra = np.array(intra); inter = np.array(inter)
    D_all = np.array(D_all); dq_all = np.array(dq_all); dm_all = np.array(dm_all)
    same_all = np.array(same_all, float)

    def ecdf(x):
        xs = np.sort(x)
        return xs, np.arange(1, len(xs) + 1) / len(xs)

    fig, ax = plt.subplots(figsize=(8, 5))
    for arr_, lbl in [(intra, "intra (mesmo padrão)"), (inter, "inter (padrões distintos)")]:
        xs, ys = ecdf(arr_)
        ax.step(xs, ys, where="post", label="%s (n=%d)" % (lbl, len(arr_)))
    ax.set_xlabel("distância 1 − Pearson"); ax.set_ylabel("ECDF")
    ax.set_title("ECDF das distâncias intra vs inter — 4 padrões "
                 "(descritor WPS b=0 branqueado, piso 50 px)", fontsize=10)
    ax.legend(); ax.grid(True, alpha=0.3)
    out_png = os.path.join(HERE, "v8-estat-ecdf.png")
    fig.savefig(out_png, dpi=140, bbox_inches="tight"); plt.close(fig)

    cq = np.corrcoef(D_all, dq_all)[0, 1]
    cm = np.corrcoef(D_all, dm_all)[0, 1]
    cs = np.corrcoef(D_all, same_all)[0, 1]
    out = (
        "[A2] ECDF distâncias (figura: v8-estat-ecdf.png)\n"
        "  grupo        n     mediana   média\n"
        "  intra      %3d     %.3f     %.3f\n"
        "  inter      %3d     %.3f     %.3f\n"
        "  corr(D, |Δq|)        = %+.3f   (proximidade de ponto operacional em q)\n"
        "  corr(D, |Δm|)        = %+.3f   (efeito/confound de fluxo mássico)\n"
        "  corr(D, mesmo padrão)= %+.3f   (efeito de padrão)\n"
        % (len(intra), np.median(intra), np.mean(intra),
           len(inter), np.median(inter), np.mean(inter), cq, cm, cs)
    )
    print(out, end=""); fh.write(out)
    return cq, cm, cs


# ── A4: diagnóstico de confound (pureza dos estratos) ─────────────────────────

def a4_confound(labels, q, m, fh):
    def purity(strata, name):
        idx_by_s = _strata_index(strata)
        lines = ["  estratos por %s:\n" % name]
        n_pure = 0
        for s in sorted(idx_by_s, key=lambda x: (float("inf") if x is None else x)):
            c = Counter(labels[i] for i in idx_by_s[s])
            top = max(c.values()); tot = sum(c.values())
            pure = (len(c) == 1)
            near = (top / tot >= 0.75) and not pure
            tag = "PURO" if pure else ("quase-puro" if near else "misto")
            if pure or near:
                n_pure += 1
            comp = ", ".join("%s=%d" % (PATTERN_ABBR[k], v) for k, v in c.items())
            lines.append("    %-4s {%s}  [%s]\n" % (str(s), comp, tag))
        lines.append("    -> %d de %d estratos puros/quase-puros\n" % (n_pure, len(idx_by_s)))
        return "".join(lines), n_pure, len(idx_by_s)

    tq, nq, totq = purity(q, "fluxo de calor q")
    tm, nm, totm = purity(m, "fluxo mássico m")
    out = (
        "[A4] Condicionamento dos estratos (por que m é a referência)\n"
        + tq + tm +
        "  VEREDITO: q ordena quase perfeitamente os padrões (plug/slug em q baixo, churn médio,\n"
        "            anular alto) → %d/%d estratos de q puros/quase-puros ⇒ a estratificação por q\n"
        "            remove quase todo o sinal de padrão (MAL condicionada). m mistura ≥3 padrões\n"
        "            em cada estrato (%d/%d puros) ⇒ m-estratificação é a leitura confiável.\n"
        % (nq, totq, nm, totm)
    )
    print(out, end=""); fh.write(out)


# ── A3': robustez ao artefato de costura (piso50 vs piso10) ──────────────────

def _stitch_band_energy_fraction(d):
    import pywt
    cf = pywt.central_frequency("cmor2.0-0.9549")
    wl = d["scales"] / cf
    sm = d["shift_mean"].astype(float)
    fracs = []
    for i in range(d["wps_b0"].shape[0]):
        if np.isnan(sm[i]):
            continue
        band = (wl >= 0.7 * sm[i]) & (wl <= 1.4 * sm[i])
        tot = d["wps_b0"][i].sum()
        if tot > 0:
            fracs.append(d["wps_b0"][i][band].sum() / tot)
    return float(np.mean(fracs)) if fracs else np.nan


def a3prime(data, fh):
    out = ["[A3'] Robustez ao artefato de costura — R GLOBAL em piso50 (λ≥50px, principal) vs\n"
           "      piso10 (λ≥10px, inclui a banda λ≈shift que o piso50 exclui):\n"]
    R = {}
    for variant in ("cropped", "center"):
        R[variant] = {}
        for floor in ("piso50", "piso10"):
            d = data[(variant, floor)]
            R[variant][floor] = anosim_R(dist_matrix(d["Z"]), to_int(d["labels"]))
        out.append("  %-8s  R(piso50)=%.4f   R(piso10)=%.4f   Δ=%+.4f\n"
                   % (variant, R[variant]["piso50"], R[variant]["piso10"],
                      R[variant]["piso10"] - R[variant]["piso50"]))
    ef = _stitch_band_energy_fraction(data[("cropped", "piso10")])
    out.append("  energia média na banda de costura (cropped, piso10) = %.1f%% do |WPS b=0|\n"
               % (100 * ef))
    cr50, cr10 = R["cropped"]["piso50"], R["cropped"]["piso10"]
    rises = cr10 > cr50 + 1e-4
    if rises and ef < 0.05:
        out.append(
            "  VEREDITO: R sobe (%.3f→%.3f) ao incluir pequena escala, MAS a banda de costura\n"
            "            carrega energia ínfima (%.1f%%) e não é o pico → a subida NÃO é o artefato\n"
            "            de costura; vem de conteúdo difuso de pequena escala colinear com\n"
            "            velocidade/comprimento (confound). piso50 = leitura CONSERVADORA.\n"
            % (cr50, cr10, 100 * ef))
    elif rises:
        out.append("  VEREDITO: R sobe (%.3f→%.3f) e a banda carrega %.1f%% → parte pode ser\n"
                   "            costura; piso50 conservador.\n" % (cr50, cr10, 100 * ef))
    else:
        out.append("  VEREDITO: R não sobe → A1 (piso50) robusto à inclusão de pequena escala.\n")
    txt = "".join(out)
    print(txt, end=""); fh.write(txt)


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    data = {}
    for variant in ("cropped", "center"):
        for floor in ("piso50", "piso10"):
            p = os.path.join(HERE, "descriptors_%s_%s.npz" % (variant, floor))
            d = np.load(p, allow_pickle=True)
            data[(variant, floor)] = {k: d[k] for k in d.files}

    res_path = os.path.join(HERE, "stats_results3.txt")
    with open(res_path, "w") as fh:
        d0 = data[("cropped", "piso50")]       # principal
        labels = d0["labels"]

        corr = np.corrcoef(d0["shift_mean"].astype(float), d0["n_original"].astype(float))[0, 1]
        warn = ("[AVISO] corr(shift_mean, n_samples) = %.4f — mesma variável de velocidade.\n\n"
                % corr)
        print(warn, end=""); fh.write(warn)

        for p in PATTERN_ORDER:
            sel = labels == p
            line = ("  %-13s n=%2d  m̄=%.2f  q̄=%.2f  n_samples̄=%.0f\n"
                    % (p, sel.sum(), np.nanmean(d0["m"][sel]),
                       np.nanmean(d0["q"][sel]), np.nanmean(d0["n_original"][sel])))
            print(line, end=""); fh.write(line)
        fh.write("\n"); print()

        # A0 — Figura 17 (matriz de distância)
        D = dist_matrix(d0["Z"])
        figure_dmatrix(D, labels, os.path.join(HERE, "v8-estat-dmatrix.png"))
        fh.write("[A0] matriz de distância 26x26 salva em v8-estat-dmatrix.png "
                 "(ordem: churn, anular, plug, slug)\n\n")

        a1_global(d0["Z"], labels, d0["q"], d0["m"], fh); fh.write("\n"); print()
        a1b_pairwise(d0["Z"], labels, fh); fh.write("\n"); print()
        a2(d0["Z"], labels, d0["q"], d0["m"], fh); fh.write("\n"); print()
        a4_confound(labels, d0["q"], d0["m"], fh); fh.write("\n"); print()
        a3prime(data, fh)

    print("\nsalvo: %s + v8-estat-dmatrix.png + v8-estat-ecdf.png" % res_path)


if __name__ == "__main__":
    main()
