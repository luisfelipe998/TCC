"""
run_stats.py — núcleo mínimo A1/A2/A3 sobre os descritores por run.

Tudo sai de UMA matriz D = 1 - Pearson(descritor) (19x19). Sem scikit-bio, sem ML, sem
PERMANOVA, sem PERMDISP (redundantes: ANOSIM-R corr 0,98 com pseudo-F sobre 19 pontos).

  A1  ANOSIM R (Clarke) + p por ENUMERAÇÃO EXATA dos 200 arranjos q-estratificados.
      Reporta R estratificado E não-estratificado.
  A2  ECDF das distâncias INTRA (mesmo padrão) vs INTER (padrões diferentes), inter separado
      em q-vizinho vs q-distante. Figura. A coesão intra é a cauda esquerda (subsome C1).
  A3  Robustez: o MESMO R de A1 em cropped/center x banda-de-costura completa/mascarada.

CAVEAT (impresso): shift_mean tem corr≈1 com n_samples — mascarar a banda de costura remove
artefato E sinal físico colinear juntos; o center (mesma costura, colapso diferente) é a
desambiguação. Nunca concluir "artefato isolado" só disto.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from itertools import combinations, product

HERE = os.path.dirname(os.path.abspath(__file__))


# ── matriz de distância ──────────────────────────────────────────────────────

def dist_matrix(Z):
    """D = 1 - Pearson entre linhas (runs). Z: (n_runs, n_scales)."""
    C = np.corrcoef(Z)
    return 1.0 - C


# ── A1: ANOSIM R + p exato ───────────────────────────────────────────────────

def anosim_R(D, groups):
    """R de Clarke. groups: array de rótulos (0/1) por run."""
    n = D.shape[0]
    iu = np.triu_indices(n, 1)
    d = D[iu]
    same = groups[iu[0]] == groups[iu[1]]
    ranks = d.argsort().argsort().astype(float) + 1.0     # rank 1..M
    r_between = ranks[~same].mean()
    r_within  = ranks[same].mean()
    M = len(d)
    return (r_between - r_within) / (M / 2.0)


def stratified_arrangements(strata_labels, q, is_churn):
    """Enumera TODOS os arranjos de rótulo churn/annular permutando só DENTRO de cada estrato q.

    Devolve lista de vetores booleanos (True=churn). Estratos puros (q=5, q=80) não movem.
    """
    idx_by_q = {}
    for i, qi in enumerate(q):
        idx_by_q.setdefault(qi, []).append(i)

    per_q_options = []
    for qi in sorted(idx_by_q):
        members = idx_by_q[qi]
        k = int(sum(is_churn[i] for i in members))        # nº de churn a colocar neste estrato
        opts = []
        for churn_positions in combinations(members, k):
            assign = {i: False for i in members}
            for i in churn_positions:
                assign[i] = True
            opts.append(assign)
        per_q_options.append(opts)

    arrangements = []
    for combo in product(*per_q_options):
        vec = np.zeros(len(q), bool)
        for assign in combo:
            for i, c in assign.items():
                vec[i] = c
        arrangements.append(vec)
    return arrangements


def a1(Z, labels, q, tag, mask=None, fh=None):
    is_churn = labels == "churn"
    groups = is_churn.astype(int)
    Zm = Z[:, ~mask] if mask is not None else Z
    D = dist_matrix(Zm)

    R_obs = anosim_R(D, groups)

    # p estratificado por enumeração exata
    arr = stratified_arrangements(None, q, is_churn)
    Rs = np.array([anosim_R(D, a.astype(int)) for a in arr])
    p_strat = float(np.mean(Rs >= R_obs))     # fração de arranjos com R >= observado

    # R não-estratificado: R observado é o mesmo; para contexto, distribuição sob permutação livre
    n = len(groups)
    rng = np.random.default_rng(0)
    R_free = []
    for _ in range(5000):
        g = rng.permutation(groups)
        R_free.append(anosim_R(D, g))
    R_free = np.array(R_free)
    p_free = float(np.mean(R_free >= R_obs))

    out = (
        "[A1] %s\n"
        "  ANOSIM R (observado)          = %.4f\n"
        "  p estratificado por q (exato) = %.4f  (n_arranjos=%d; p mínimo=%.4f)\n"
        "  p não-estratificado (5000 MC) = %.4f\n"
        "  -> efeito-tamanho DESCRITIVO, NÃO teste inferencial forte (q prediz padrão)\n"
        % (tag, R_obs, p_strat, len(arr), 1.0/len(arr), p_free)
    )
    print(out, end="")
    if fh: fh.write(out)
    return R_obs, p_strat, p_free


# ── A2: ECDF intra vs inter ──────────────────────────────────────────────────

def a2(Z, labels, q, fh):
    D = dist_matrix(Z)
    n = D.shape[0]
    intra, inter_near, inter_far = [], [], []
    for i, j in combinations(range(n), 2):
        d = D[i, j]
        if labels[i] == labels[j]:
            intra.append(d)
        else:
            (inter_near if abs(q[i] - q[j]) <= 20 else inter_far).append(d)
    intra = np.array(intra); inter_near = np.array(inter_near); inter_far = np.array(inter_far)

    def ecdf(x):
        xs = np.sort(x)
        return xs, np.arange(1, len(xs) + 1) / len(xs)

    fig, ax = plt.subplots(figsize=(8, 5))
    for arr_, lbl in [(intra, "intra (mesmo padrão)"),
                      (inter_near, "inter q-vizinho (|Δq|≤20)"),
                      (inter_far, "inter q-distante")]:
        if len(arr_):
            xs, ys = ecdf(arr_)
            ax.step(xs, ys, where="post", label="%s (n=%d)" % (lbl, len(arr_)))
    ax.set_xlabel("distância 1 − Pearson"); ax.set_ylabel("ECDF")
    ax.set_title("A2 — distâncias intra vs inter (descritor WPS b=0 whitened)")
    ax.legend(); ax.grid(True, alpha=0.3)
    out_png = os.path.join(HERE, "a2_ecdf.png")
    fig.savefig(out_png, dpi=140, bbox_inches="tight"); plt.close(fig)

    out = (
        "[A2] ECDF distâncias (figura: a2_ecdf.png)\n"
        "  intra      n=%3d  mediana=%.3f\n"
        "  inter viz  n=%3d  mediana=%s\n"
        "  inter dist n=%3d  mediana=%s\n"
        "  -> a coesão intra é a cauda esquerda deste ECDF (subsome C1)\n"
        % (len(intra), np.median(intra),
           len(inter_near), ("%.3f" % np.median(inter_near)) if len(inter_near) else "—",
           len(inter_far),  ("%.3f" % np.median(inter_far))  if len(inter_far)  else "—")
    )
    print(out, end=""); fh.write(out)


# ── A3: robustez (4 valores de R) ────────────────────────────────────────────

def stitch_band_mask(scales, shift_mean):
    """Máscara de escalas na banda de costura λ∈[0.7,1.4]×shift_mean médio dos runs.

    A CWT indexa por ESCALA (scales = center_freq/freq); o comprimento de onda físico é
    λ = 1/freq = scales/center_freq. A banda de costura é dada em px de comprimento de onda,
    então converte-se para comprimento de onda ANTES de comparar.
    """
    import pywt
    center_freq = pywt.central_frequency("cmor2.0-0.9549")
    wavelength = scales / center_freq            # px
    sm = shift_mean.astype(float)
    if np.all(np.isnan(sm)):
        return np.zeros(scales.shape, bool)
    lam = float(np.nanmean(sm))
    return (wavelength >= 0.7 * lam) & (wavelength <= 1.4 * lam)


def a3(data_by_variant, fh):
    out = ["[A3] Robustez — mesmo ANOSIM R em 4 condições:\n"]
    for variant, d in data_by_variant.items():
        for band, mask in [("banda-completa", None),
                            ("banda-mascarada", stitch_band_mask(d["scales"], d["shift_mean"]))]:
            R, _, _ = a1(d["Z"], d["labels"], d["q"],
                         tag="(A3 %s / %s)" % (variant, band), mask=mask)
            nmask = 0 if mask is None else int(mask.sum())
            out.append("  %-8s %-16s R=%.4f  (escalas mascaradas=%d)\n"
                       % (variant, band, R, nmask))
    txt = "".join(out)
    fh.write(txt)
    print(txt, end="")


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    data = {}
    for variant in ("cropped", "center"):
        p = os.path.join(HERE, "descriptors_%s.npz" % variant)
        d = np.load(p, allow_pickle=True)
        data[variant] = {k: d[k] for k in d.files}

    res_path = os.path.join(HERE, "stats_results.txt")
    with open(res_path, "w") as fh:
        # aviso de colinearidade shift ⟂ n
        d0 = data["cropped"]
        corr = np.corrcoef(d0["shift_mean"].astype(float), d0["n_original"].astype(float))[0, 1]
        warn = ("[AVISO] corr(shift_mean, n_samples) = %.4f — são a MESMA variável de "
                "velocidade;\n        a banda de costura (A3) é colinear com comprimento/"
                "velocidade/padrão.\n\n" % corr)
        print(warn, end=""); fh.write(warn)

        # confounds descritivos
        for lbl in ("churn", "wavy annular"):
            sel = d0["labels"] == lbl
            line = ("  %-13s n=%2d  m̄=%.2f  q̄=%.2f  n_samples̄=%.0f\n"
                    % (lbl, sel.sum(), np.nanmean(d0["m"][sel]),
                       np.nanmean(d0["q"][sel]), np.nanmean(d0["n_original"][sel])))
            print(line, end=""); fh.write(line)
        fh.write("\n")

        # A1 principal (cropped, banda completa)
        print("\n"); fh.write("\n")
        a1(d0["Z"], d0["labels"], d0["q"], tag="(principal: cropped, banda-completa)", fh=fh)
        fh.write("\n")

        # A2 (cropped)
        print(); a2(d0["Z"], d0["labels"], d0["q"], fh)
        fh.write("\n")

        # A3 (4 R)
        print(); a3(data, fh)

    print("\nsalvo: %s  e  a2_ecdf.png" % res_path)


if __name__ == "__main__":
    main()
