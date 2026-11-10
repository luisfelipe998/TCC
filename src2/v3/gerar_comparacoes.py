"""Orquestra wavelet_v3.py sobre data_processed."""
import os
import glob
import subprocess
import sys
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.join(HERE, "..", "data_processed")
SCRIPT = os.path.join(HERE, "wavelet_v3.py")


def folder_flow_type(name):
    folder = os.path.join(DATA_ROOT, name)
    for p in sorted(glob.glob(os.path.join(folder, "*.txt"))):
        stem = os.path.splitext(os.path.basename(p))[0].lstrip("_")
        if stem != "dados":
            return stem
    return name


def list_runs():
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

    pairs = []
    for pattern, members in sorted(groups.items()):
        if len(members) < 2:
            print("skip '%s': apenas %d run (sem par)" % (pattern, len(members)))
            continue
        for a, b in combinations(sorted(members), 2):
            pairs.append((pattern, a, b))

    print("total de pares (mesmo padrao): %d" % len(pairs))
    for i, (pattern, a, b) in enumerate(pairs, 1):
        print("[%d/%d] %s: %s x %s" % (i, len(pairs), pattern, a, b))
        subprocess.run([sys.executable, SCRIPT, a, b], check=True)


if __name__ == "__main__":
    main()
