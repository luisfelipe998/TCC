"""Orquestra wavelet_alternativa_v2.py (v7) sobre data_processed.

v7: para cada padrao de escoamento (flow pattern), roda uma analise por run,
usando aquela run como HOLDOUT (sinal 2 = run individual crua) e a media das
autocorrelacoes das DEMAIS runs do padrao como sinal 1 (assinatura media da
classe, hold-out para evitar vazamento). Padroes com < 2 runs sao ignorados
(sem hold-in para formar a media).

Chama o script no modo default (salva graficos em subpasta
<PATTERN>_mean_K<n>_<PATTERN>__x__<HOLDOUT>_<FLOW2>).
"""
import os
import glob
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.join(HERE, "..", "data_processed")
SCRIPT = os.path.join(HERE, "wavelet_alternativa_v2.py")


def folder_flow_type(name):
    """Padrao de escoamento = nome do .txt marcador da pasta (mesma regra do script)."""
    folder = os.path.join(DATA_ROOT, name)
    for p in sorted(glob.glob(os.path.join(folder, "*.txt"))):
        stem = os.path.splitext(os.path.basename(p))[0]
        if stem != "dados":
            return stem
    return name


def list_runs():
    """Pastas de data_processed que tem dados.xlsx (uma run cada)."""
    out = []
    for d in sorted(os.listdir(DATA_ROOT)):
        if os.path.exists(os.path.join(DATA_ROOT, d, "dados.xlsx")):
            out.append(d)
    return out


def group_by_pattern(runs):
    groups = {}
    for r in runs:
        groups.setdefault(folder_flow_type(r), []).append(r)
    return groups


def main():
    runs = list_runs()
    groups = group_by_pattern(runs)

    # Uma analise por (padrao, run-holdout). Precisa >= 2 runs para ter hold-in.
    jobs = []
    for pattern, members in sorted(groups.items()):
        if len(members) < 2:
            print("skip '%s': apenas %d run (sem hold-in para a media)" % (pattern, len(members)))
            continue
        for holdout in sorted(members):
            jobs.append((pattern, holdout))

    print("total de analises (media da classe x run holdout): %d" % len(jobs))
    for i, (pattern, holdout) in enumerate(jobs, 1):
        print("[%d/%d] %s: mean(K=%d) x %s" % (i, len(jobs), pattern, len(groups[pattern]) - 1, holdout))
        subprocess.run([sys.executable, SCRIPT, pattern, holdout], check=True)


if __name__ == "__main__":
    main()
