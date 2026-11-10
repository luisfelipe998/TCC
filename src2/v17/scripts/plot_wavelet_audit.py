"""Plots from audit_wavelets.py tables; outputs restricted to v17."""
from pathlib import Path
import os
import sys
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR",str(ROOT/"logs/mpl_wavelets"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

DATA=ROOT/"resultados/wavelets"
null=pd.read_csv(DATA/"resumo_nulos_brancos.csv")
windows=pd.read_csv(DATA/"sensibilidade_janelas.csv")
fig,ax=plt.subplots(1,2,figsize=(11,4.6),layout="constrained")
colors=["#16697a","#e3853f","#725b9f"]
for (n,g),color in zip(null[null.band=="fixed_pixels"].groupby("n"),colors):
    ax[0].plot(g.lambda_px,g.peak_median,"o-",label=f"N = {n:,}".replace(",","."),color=color)
    ax[0].fill_between(g.lambda_px,g.peak_q05,g.peak_q95,alpha=.10,color=color)
ax[0].set(xscale="log",xlabel="Comprimento de onda λ (px)",ylabel="Máximo da WCO em atraso",
          ylim=(0,1.02),title="Ruídos independentes com a mesma distribuição")
ax[0].legend(fontsize=9)
ax[0].text(.03,.98,"Mediana e faixa 5–95% de 48 pares por N",transform=ax[0].transAxes,va="top",fontsize=9)
names={"full_50_3000":"Completo · 50–3.000","full_50_1024":"Completo · 50–1.024",
       "prefix_3171_50_1024":"Prefixo 3.171 · 50–1.024","center_3171_50_1024":"Centro 3.171 · 50–1.024",
       "suffix_3171_50_1024":"Sufixo 3.171 · 50–1.024","prefix_3171_50_3000":"Prefixo 3.171 · 50–3.000",
       "center_3171_50_3000":"Centro 3.171 · 50–3.000","suffix_3171_50_3000":"Sufixo 3.171 · 50–3.000"}
y=np.arange(len(windows))
ax[1].barh(y,windows.delta,color=["#e3853f" if n.endswith("3000") else "#16697a" for n in windows.config])
ax[1].set(yticks=y,yticklabels=[names[n] for n in windows.config],xlabel="Δr = média intra − média inter",
          title="Sensibilidade do contraste de Pearson",xlim=(0,.31))
ax[1].invert_yaxis()
for i,v in enumerate(windows.delta):ax[1].text(v+.003,i,f"{v:.3f}".replace(".",","),va="center",fontsize=9)
fig.suptitle("Auditoria v16 · resultados exploratórios, sem aprendizado de máquina",fontsize=13)
for a in ax:a.grid(axis="x",alpha=.2)
fig.savefig(DATA/"diagnosticos_wavelets.png",dpi=180)
fig.savefig(DATA/"diagnosticos_wavelets.pdf")
plt.close(fig)
print("Saved diagnosticos_wavelets.png/pdf")
