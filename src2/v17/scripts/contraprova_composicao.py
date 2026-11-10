"""Independent review of composition() using matrices and pair JSON, not its CSV source."""
from pathlib import Path
import csv
import json
import sys
sys.dont_write_bytecode=True
import numpy as np

BASE=Path(__file__).resolve().parents[1]
SOURCE=BASE.parent/"v16_spatial2_18_71"
OUT=BASE/"resultados/wavelets"
records=json.loads((SOURCE/"intermediario/aquisicoes.json").read_text())
runs=[r["ensaio"] for r in records]
labels=[r["padrao"] for r in records]
with np.load(SOURCE/"intermediario/descritores.npz") as a:
    assert a["runs"].tolist()==runs
    # Independent Pearson reconstruction from saved power, with scalar OLS.
    t=np.log(a["wavelengths"]);t-=t.mean()
    logp=np.log(a["power"]+1e-12);logp-=logp.mean(axis=1,keepdims=True)
    z=logp-(logp@t/np.dot(t,t))[:,None]*t
    pearson=np.corrcoef(z)
score=np.full((len(runs),len(runs)),np.nan)
for p in json.loads((SOURCE/"intermediario/pares.json").read_text()):
    i,j=runs.index(p["run_x"]),runs.index(p["run_y"])
    score[i,j]=score[j,i]=p["escore_piso05"]

reference_path=BASE/"resultados/exploracao/efeito_composicao_tcc.csv"
with reference_path.open(newline="") as f: reference=list(csv.DictReader(f))
lookup={(r["conjunto"],r["metrica"],r["tipo"]):r for r in reference}
rows=[];max_error=0.
for name,exclusions in [("26_todos",[]),("25_sem_bolha",["m07q01_001"]),
                       ("25_sem_m03q40",["m03q40_001"]),
                       ("24_sem_ambos",["m03q40_001","m07q01_001"])]:
    keep=[i for i,r in enumerate(runs) if r not in exclusions]
    for metric,matrix in [("pearson_r",pearson),("escore_piso05",score)]:
        grouped={};wi=[];be=[]
        for a,i in enumerate(keep):
            for j in keep[a+1:]:
                key=tuple(sorted((labels[i],labels[j])))
                value=float(matrix[i,j])
                grouped.setdefault(key,[]).append(value)
                (wi if labels[i]==labels[j] else be).append(value)
        wm=np.mean([np.mean(v) for (a,b),v in grouped.items() if a==b])
        bm=np.mean([np.mean(v) for (a,b),v in grouped.items() if a!=b])
        row=dict(conjunto=name,metrica=metric,n=len(keep),pares_intra=len(wi),pares_inter=len(be),
                 delta=float(np.mean(wi)-np.mean(be)),delta_macro=float(wm-bm))
        ref=lookup[name,metric,"intra_menos_inter"]
        error=max(abs(row["delta"]-float(ref["media"])),abs(row["delta_macro"]-float(ref["media_balanceada"])))
        max_error=max(max_error,error)
        row["erro_contra_root"]=error
        rows.append(row)
with (OUT/"contraprova_composicao.csv").open("w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
assert max_error<1e-12
all_row=next(r for r in rows if r["conjunto"]=="26_todos" and r["metrica"]=="pearson_r")
without_row=next(r for r in rows if r["conjunto"]=="25_sem_bolha" and r["metrica"]=="pearson_r")
decomposition=dict(status="PASS",maximum_error=max_error,
    delta_reduction_when_singleton_removed=all_row["delta"]-without_row["delta"],
    relative_delta_reduction=(all_row["delta"]-without_row["delta"])/all_row["delta"],
    macro_reduction_when_singleton_removed=all_row["delta_macro"]-without_row["delta_macro"],
    relative_macro_reduction=(all_row["delta_macro"]-without_row["delta_macro"])/all_row["delta_macro"],
    between_group_weight_with_singleton=4/10,between_pair_weight_with_singleton=25/238,
    method="Rebuilt Pearson from fixed profiles; WCO from pair JSON. Independent loops over acquisition indices; no composition() call.")
(OUT/"contraprova_composicao.json").write_text(json.dumps(decomposition,indent=2,ensure_ascii=False)+"\n")
print(json.dumps(decomposition,indent=2))
