"""Orchestrates wavelet_v8.py over all INTER-pattern run pairs in data_processed.

Diferença para o v8 (intra-padrão): aqui cada par (A, B) tem A e B de padrões
DIFERENTES (ex.: churn1 x annular1, churn1 x annular2, churn2 x annular1, ...).
Restringe-se aos padrões que também apareciam no v8 (>= 2 runs): churn,
wavy annular, plug, slug, slug-churn. Padrões degenerados (n=1: plug-slug,
dispersed bubbles, smooth annular) são ignorados, exatamente como no v8.

Para cada par de padrões distintos gera o produto cartesiano completo dos runs.
A ORDEM dentro do par (quem é a base NAME1 passada ao wavelet_v8.py) segue a
prioridade: churn > wavy annular > demais (lexicográfico). Assim churn e wavy
annular ficam sempre como primeiro sinal quando participam; em churn x wavy
annular a base é churn.
"""
import os
import glob
import subprocess
import sys
from itertools import combinations, product

HERE      = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.join(HERE, "../", "data_processed")
SCRIPT    = os.path.join(HERE, "wavelet_v8.py")
PYTHON    = sys.executable

# Padrões considerados = os mesmos do v8 (>= 2 runs). Igual ao critério do v8,
# que pulava qualquer padrão com < 2 runs.
MIN_RUNS = 2

# Prioridade de base (NAME1): índice menor = fica como primeiro sinal.
BASE_PRIORITY = {"churn": 0, "wavy annular": 1}


def base_rank(pattern):
    return BASE_PRIORITY.get(pattern, 2)


def folder_flow_type(name):
    folder = os.path.join(DATA_ROOT, name)
    for p in sorted(glob.glob(os.path.join(folder, "*.txt"))):
        stem = os.path.splitext(os.path.basename(p))[0].lstrip("_")
        if stem != "dados":
            return stem
    return name


def list_runs():
    runs = []
    for d in sorted(os.listdir(DATA_ROOT)):
        path = os.path.join(DATA_ROOT, d)
        if (os.path.isdir(path)
                and os.path.exists(os.path.join(path, "__dados_cropped_v11.xlsx"))
                and os.path.exists(os.path.join(path, "__dados_center_v13.xlsx"))):
            runs.append(d)
    return runs


def main():
    runs   = list_runs()
    groups = {}
    for r in runs:
        groups.setdefault(folder_flow_type(r), []).append(r)

    # Mantém só os padrões com >= MIN_RUNS runs (os que o v8 processou).
    eligible = {p: sorted(members) for p, members in groups.items()
                if len(members) >= MIN_RUNS}
    for p, members in sorted(groups.items()):
        if len(members) < MIN_RUNS:
            print(f"skip '{p}': only {len(members)} run")

    # Pares de PADRÕES distintos; para cada, produto cartesiano dos runs.
    # A base (primeiro padrão) é a de maior prioridade: churn > wavy annular > resto.
    jobs = []
    for pat_x, pat_y in combinations(sorted(eligible), 2):
        pat_a, pat_b = sorted((pat_x, pat_y), key=lambda p: (base_rank(p), p))
        for a, b in product(eligible[pat_a], eligible[pat_b]):
            jobs.append((pat_a, pat_b, a, b))

    print(f"total inter-pattern pairs: {len(jobs)}")
    expected = ["signals_cropped.png", "signals_center.png",
                "wavelet_cropped.png", "wavelet_center.png"]
    for i, (pat_a, pat_b, a, b) in enumerate(jobs, 1):
        out = os.path.join(HERE, f"{a}_{pat_a}__x__{b}_{pat_b}")
        if os.path.isdir(out) and all(
                os.path.exists(os.path.join(out, p)) for p in expected):
            print(f"[{i}/{len(jobs)}] skip (done) {pat_a} x {pat_b}: {a} x {b}")
            continue
        print(f"[{i}/{len(jobs)}] {pat_a} x {pat_b}: {a} x {b}")
        subprocess.run([PYTHON, SCRIPT, a, b], check=True)


if __name__ == "__main__":
    main()
