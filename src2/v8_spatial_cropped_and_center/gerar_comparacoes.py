"""Orchestrates wavelet_v8.py over all same-pattern pairs in data_processed.

For each flow pattern with >= 2 runs, generates every unique pair (A, B) where
A < B (lexicographic) and calls wavelet_v8.py. Skips patterns with < 2 runs.
"""
import os
import glob
import subprocess
import sys
from itertools import combinations

HERE      = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.join(HERE, "../", "data_processed")
SCRIPT    = os.path.join(HERE, "wavelet_v8.py")
PYTHON    = sys.executable


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

    jobs = []
    for pattern, members in sorted(groups.items()):
        if len(members) < 2:
            print(f"skip '{pattern}': only {len(members)} run")
            continue
        for a, b in combinations(sorted(members), 2):
            jobs.append((pattern, a, b))

    print(f"total pairs: {len(jobs)}")
    for i, (pattern, a, b) in enumerate(jobs, 1):
        print(f"[{i}/{len(jobs)}] {pattern}: {a} x {b}")
        subprocess.run([PYTHON, SCRIPT, a, b], check=True)


if __name__ == "__main__":
    main()
