"""
plot_dmatrix.py — mapa de calor da matriz de distância D = 1 - Pearson entre os 19 ensaios.

Figura puramente ILUSTRATIVA para o TCC. NÃO recomputa estatística: lê o descritor já extraído
(descriptors_cropped_piso50.npz, o mesmo da análise A1 principal) e reusa a MESMA definição de
distância de run_stats.dist_matrix (D = 1 - np.corrcoef(Z)). Ordena os ensaios por padrão
(churn primeiro, anular ondulado depois) para que os blocos intra-padrão apareçam na diagonal.

Saída: ../../TCC/figuras/v8-estat-dmatrix.png
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from run_stats import dist_matrix        # reusa a definição exata (1 - Pearson)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT  = os.path.join(HERE, "..", "..", "TCC", "figuras", "v8-estat-dmatrix.png")

LABEL_PT = {"churn": "churn", "wavy annular": "anular ondulado"}


def main():
    d = np.load(os.path.join(HERE, "descriptors_cropped_piso50.npz"), allow_pickle=True)
    Z = d["Z"]
    labels = d["labels"].astype(str)

    # ordem: churn primeiro, depois anular ondulado (blocos intra na diagonal)
    order = np.concatenate([np.where(labels == "churn")[0],
                            np.where(labels == "wavy annular")[0]])
    Zo = Z[order]
    labo = labels[order]
    n_churn = int(np.sum(labo == "churn"))

    D = dist_matrix(Zo)

    fig, ax = plt.subplots(figsize=(6.6, 5.6))
    im = ax.imshow(D, cmap="viridis_r", vmin=0.0, vmax=float(D.max()))
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(r"distância $1-\rho$ de Pearson")

    # linhas separando os dois blocos de padrão
    ax.axhline(n_churn - 0.5, color="white", lw=1.5)
    ax.axvline(n_churn - 0.5, color="white", lw=1.5)

    ticks = np.arange(len(labo))
    tick_lbl = ["%s" % (i + 1) for i in range(len(labo))]
    ax.set_xticks(ticks); ax.set_xticklabels(tick_lbl, fontsize=7)
    ax.set_yticks(ticks); ax.set_yticklabels(tick_lbl, fontsize=7)
    ax.set_xlabel("ensaio")
    ax.set_ylabel("ensaio")

    # rótulos de bloco
    ax.text((n_churn - 1) / 2, -0.9, "churn (10)", ha="center", va="bottom", fontsize=9)
    ax.text(n_churn + (len(labo) - n_churn - 1) / 2, -0.9, "anular ondulado (9)",
            ha="center", va="bottom", fontsize=9)

    ax.set_title("Matriz de distância espectral entre os 19 ensaios\n"
                 "(descritor WPS no atraso nulo, branqueado)", fontsize=10, pad=18)
    fig.tight_layout()
    fig.savefig(OUT, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("salvo:", os.path.abspath(OUT))
    print("D min/median/max = %.3f / %.3f / %.3f"
          % (D[np.triu_indices(len(labo), 1)].min(),
             np.median(D[np.triu_indices(len(labo), 1)]),
             D[np.triu_indices(len(labo), 1)].max()))


if __name__ == "__main__":
    main()
