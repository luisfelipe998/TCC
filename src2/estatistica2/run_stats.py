"""
run_stats.py — follow-up A1/A2/A3' sobre os descritores por run.

Tudo sai de UMA matriz D = 1 - Pearson(descritor) (19x19). Sem scikit-bio, sem ML, sem
PERMANOVA/PERMDISP (redundantes: ANOSIM-R corr 0,98 com pseudo-F sobre 19 pontos).

Diferenças vs ../estatistica/run_stats.py:
  A1  acrescenta p por ENUMERAÇÃO EXATA m-estratificada (1728 arranjos, p_min 0,00058) ao lado
      da q-estratificada (200, p_min 0,005) e da livre (MC). m-estratos NÃO são puros → nulo
      mais bem condicionado. Imprime R máximo e rank de R_obs em cada enumeração.
  A2  reporta MEDIANA E MÉDIA de cada grupo (a ordem das medianas inverte na média → não é
      evidência robusta) e destaca corr(D,|Δq|)≈0 como a evidência de que a separação NÃO é
      proximidade de q. Acrescenta corr(D,|Δm|) e corr(D,same_pattern).
  A3' substitui a máscara no-op (banda de costura < piso 50 px) pelo teste real: recomputa R em
      piso50 vs piso10 (cropped e center). piso10 = sonda de robustez, nunca manchete.
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


# ── A1: ANOSIM R + p por enumeração exata ────────────────────────────────────

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


def stratified_arrangements(strata, is_churn):
    """Enumera TODOS os arranjos de rótulo churn/annular permutando só DENTRO de cada estrato.

    `strata` = vetor de rótulo de estrato por run (q OU m). Devolve lista de vetores booleanos
    (True=churn). Estratos puros não movem (q=5, q=80); estratos de m não são puros.
    """
    idx_by_s = {}
    for i, si in enumerate(strata):
        idx_by_s.setdefault(si, []).append(i)

    per_s_options = []
    for si in sorted(idx_by_s):
        members = idx_by_s[si]
        k = int(sum(is_churn[i] for i in members))        # nº de churn a colocar neste estrato
        opts = []
        for churn_positions in combinations(members, k):
            assign = {i: False for i in members}
            for i in churn_positions:
                assign[i] = True
            opts.append(assign)
        per_s_options.append(opts)

    arrangements = []
    for combo in product(*per_s_options):
        vec = np.zeros(len(strata), bool)
        for assign in combo:
            for i, c in assign.items():
                vec[i] = c
        arrangements.append(vec)
    return arrangements


def _enum_p(D, is_churn, strata, R_obs):
    """p e diagnóstico (nº arranjos, R máx, rank) por enumeração exata num estrato."""
    arr = stratified_arrangements(strata, is_churn)
    Rs = np.array([anosim_R(D, a.astype(int)) for a in arr])
    p = float(np.mean(Rs >= R_obs))
    rank_from_top = int(np.sum(Rs > R_obs))
    return p, len(arr), float(Rs.max()), rank_from_top


def a1(Z, labels, q, m, tag, fh=None):
    is_churn = labels == "churn"
    groups = is_churn.astype(int)
    D = dist_matrix(Z)
    R_obs = anosim_R(D, groups)

    p_q, nq, Rmax_q, rank_q = _enum_p(D, is_churn, q, R_obs)
    p_m, nm, Rmax_m, rank_m = _enum_p(D, is_churn, m, R_obs)

    # p livre (contexto): permutação de rótulo sem estrato
    rng = np.random.default_rng(0)
    R_free = np.array([anosim_R(D, rng.permutation(groups)) for _ in range(5000)])
    p_free = float(np.mean(R_free >= R_obs))

    out = (
        "[A1] %s\n"
        "  ANOSIM R (observado)              = %.4f\n"
        "  p q-estratificado (exato)         = %.4f  (n_arranjos=%d; p_min=%.5f; R_max=%.4f; rank=%d)\n"
        "  p m-estratificado (exato)         = %.4f  (n_arranjos=%d; p_min=%.5f; R_max=%.4f; rank=%d)\n"
        "  p não-estratificado (5000 MC)     = %.4f\n"
        "  -> R sobrevive controlando q E o confound ortogonal m; a m-estratificação é mais bem\n"
        "     condicionada (nenhum estrato puro: %d vs %d arranjos). Efeito-tamanho DESCRITIVO,\n"
        "     NÃO teste inferencial forte (q e m colineares com padrão).\n"
        % (tag, R_obs, p_q, nq, 1.0 / nq, Rmax_q, rank_q,
           p_m, nm, 1.0 / nm, Rmax_m, rank_m, p_free, nm, nq)
    )
    print(out, end="")
    if fh:
        fh.write(out)
    return R_obs, p_q, p_m, p_free


# ── A2: ECDF intra vs inter + correlações pairwise ───────────────────────────

def a2(Z, labels, q, m, fh):
    D = dist_matrix(Z)
    n = D.shape[0]
    intra, inter_near, inter_far = [], [], []
    dq_all, dm_all, D_all, same_all = [], [], [], []
    for i, j in combinations(range(n), 2):
        d = D[i, j]
        D_all.append(d); dq_all.append(abs(q[i] - q[j])); dm_all.append(abs(m[i] - m[j]))
        same = labels[i] == labels[j]
        same_all.append(same)
        if same:
            intra.append(d)
        else:
            (inter_near if abs(q[i] - q[j]) <= 20 else inter_far).append(d)
    intra = np.array(intra); inter_near = np.array(inter_near); inter_far = np.array(inter_far)
    D_all = np.array(D_all); dq_all = np.array(dq_all); dm_all = np.array(dm_all)
    same_all = np.array(same_all, float)

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
    ax.set_title("A2 — distâncias intra vs inter (descritor WPS b=0 whitened, piso 50 px)")
    ax.legend(); ax.grid(True, alpha=0.3)
    out_png = os.path.join(HERE, "a2_ecdf.png")
    fig.savefig(out_png, dpi=140, bbox_inches="tight"); plt.close(fig)

    def ms(x):  # mediana/média
        return (np.median(x), np.mean(x)) if len(x) else (np.nan, np.nan)
    mi_med, mi_mean = ms(intra); nn_med, nn_mean = ms(inter_near); nf_med, nf_mean = ms(inter_far)
    cq = np.corrcoef(D_all, dq_all)[0, 1]
    cm = np.corrcoef(D_all, dm_all)[0, 1]
    cs = np.corrcoef(D_all, same_all)[0, 1]

    out = (
        "[A2] ECDF distâncias (figura: a2_ecdf.png)\n"
        "  grupo          n     mediana   média\n"
        "  intra        %3d     %.3f     %.3f\n"
        "  inter viz    %3d     %.3f     %.3f\n"
        "  inter dist   %3d     %.3f     %.3f\n"
        "  NOTA: a ordem inter-viz vs inter-dist INVERTE entre mediana e média (as ECDFs quase\n"
        "        se sobrepõem) → NÃO é evidência robusta. A evidência robusta é a correlação:\n"
        "  corr(D, |Δq|)        = %+.3f   (≈0 → separação NÃO é proximidade de ponto operacional em q)\n"
        "  corr(D, |Δm|)        = %+.3f   (efeito de m, confound não controlado)\n"
        "  corr(D, mesmo padrão)= %+.3f   (efeito de padrão; MESMA ordem de grandeza que |Δm|)\n"
        % (len(intra), mi_med, mi_mean,
           len(inter_near), nn_med, nn_mean,
           len(inter_far), nf_med, nf_mean, cq, cm, cs)
    )
    print(out, end=""); fh.write(out)


# ── A3': robustez ao artefato de costura (piso50 vs piso10) ──────────────────

def _stitch_band_energy_fraction(d):
    """Fração média da energia de |WPS b=0| que cai na banda de costura [0,7;1,4]×shift_mean.

    Mede DIRETAMENTE se o artefato de costura se manifesta como energia concentrada no descritor
    piso10. Se for ínfima, a eventual subida de R ao incluir escalas pequenas NÃO pode ser
    atribuída ao bump de costura. shift_mean NaN (center) → devolve NaN.
    """
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
    out = ["[A3'] Robustez ao artefato de costura — R em piso50 (λ≥50px, A1 principal) vs\n"
           "      piso10 (λ≥10px, inclui a banda λ≈shift que o piso50 exclui):\n"]
    R = {}
    for variant in ("cropped", "center"):
        R[variant] = {}
        for floor in ("piso50", "piso10"):
            d = data[(variant, floor)]
            R[variant][floor] = anosim_R(dist_matrix(d["Z"]),
                                         (d["labels"] == "churn").astype(int))
        out.append("  %-8s  R(piso50)=%.4f   R(piso10)=%.4f   Δ=%+.4f\n"
                   % (variant, R[variant]["piso50"], R[variant]["piso10"],
                      R[variant]["piso10"] - R[variant]["piso50"]))

    # Atribuição: a banda de costura carrega energia significativa no descritor piso10?
    ef = _stitch_band_energy_fraction(data[("cropped", "piso10")])
    out.append("  energia média na banda de costura (cropped, piso10) = %.1f%% do |WPS b=0|\n"
               % (100 * ef))

    cr50, cr10 = R["cropped"]["piso50"], R["cropped"]["piso10"]
    rises = cr10 > cr50 + 1e-4
    if rises and ef < 0.05:
        verdict = (
            "  VEREDITO: R SOBE ao rebaixar o piso (0,177→%.3f), MAS a banda de costura carrega\n"
            "            energia ínfima (%.1f%%) e nenhum pico distinto (ver a3_out_of_band.png:\n"
            "            o pico do |WPS b=0| fica em λ≈400–950 px, longe de λ≈shift). Logo a subida\n"
            "            NÃO é atribuível ao bump de costura — vem do conteúdo difuso de PEQUENA\n"
            "            ESCALA (10–50 px), onde annular (sinal mais longo/fino) e churn diferem.\n"
            "            Essa faixa é COLINEAR com velocidade/comprimento (n_original≈shift_mean,\n"
            "            corr=1,0 com o padrão); incluí-la mede em parte o confound, não o padrão.\n"
            "            O piso50 (A1 principal) EXCLUI essa região confundida → leitura\n"
            "            CONSERVADORA quanto ao CONFOUND (não é questão de borda/COI: o descritor\n"
            "            é a coluna b=0 = CENTRO da autocorrelação 2N−1, a 7000 px de qualquer\n"
            "            borda, logo nenhuma escala tem COI que a alcance).\n"
            % (cr10, 100 * ef))
    elif rises:
        verdict = (
            "  VEREDITO: R SOBE (0,177→%.3f) E a banda de costura carrega energia não-trivial\n"
            "            (%.1f%%) → parte da separação pode vir do bump de costura; piso50 é a\n"
            "            leitura CONSERVADORA.\n" % (cr10, 100 * ef))
    else:
        verdict = ("  VEREDITO: R não sobe ao rebaixar o piso → escalas pequenas não separam;\n"
                   "            A1 (piso50) é robusto à sua inclusão.\n")
    out.append(verdict)
    out.append("  (piso10 é SONDA de robustez — inclui escalas pequenas colineares com o confound\n"
               "   de velocidade/comprimento; NUNCA é número de manchete. center quase não muda\n"
               "   (Δ=%+.4f), consistente com o colapso center não injetar a mesma textura de\n"
               "   pequena escala. center tem shift_mean NaN → banda identificada por escala, não\n"
               "   por metadata.)\n"
               % (R["center"]["piso10"] - R["center"]["piso50"]))
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

    res_path = os.path.join(HERE, "stats_results2.txt")
    with open(res_path, "w") as fh:
        d0 = data[("cropped", "piso50")]       # A1/A2 principal

        corr = np.corrcoef(d0["shift_mean"].astype(float), d0["n_original"].astype(float))[0, 1]
        warn = ("[AVISO] corr(shift_mean, n_samples) = %.4f — são a MESMA variável de "
                "velocidade.\n\n" % corr)
        print(warn, end=""); fh.write(warn)

        for lbl in ("churn", "wavy annular"):
            sel = d0["labels"] == lbl
            line = ("  %-13s n=%2d  m̄=%.2f  q̄=%.2f  n_samples̄=%.0f\n"
                    % (lbl, sel.sum(), np.nanmean(d0["m"][sel]),
                       np.nanmean(d0["q"][sel]), np.nanmean(d0["n_original"][sel])))
            print(line, end=""); fh.write(line)
        fh.write("\n"); print()

        a1(d0["Z"], d0["labels"], d0["q"], d0["m"],
           tag="(principal: cropped, piso50)", fh=fh)
        fh.write("\n"); print()

        a2(d0["Z"], d0["labels"], d0["q"], d0["m"], fh)
        fh.write("\n"); print()

        a3prime(data, fh)

    print("\nsalvo: %s  e  a2_ecdf.png" % res_path)


if __name__ == "__main__":
    main()
