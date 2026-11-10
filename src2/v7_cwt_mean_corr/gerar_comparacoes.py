"""Orquestra wavelet_v7.py (v7_cwt_mean_corr) sobre data_processed.

Para cada padrao com K>=2 runs, faz leave-one-out:
  probe = run k, refs = todos os outros K-1 runs do mesmo padrao.
"""
import os
import glob
import subprocess
import sys

HERE      = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.join(HERE, "..", "data_processed")
SCRIPT    = os.path.join(HERE, "wavelet_v7.py")


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
        if os.path.exists(os.path.join(DATA_ROOT, d, "__dados_cropped_v11.xlsx")):
            out.append(d)
    return out


def group_by_pattern(runs):
    groups = {}
    for r in runs:
        groups.setdefault(folder_flow_type(r), []).append(r)
    return groups


def main():
    runs   = list_runs()
    groups = group_by_pattern(runs)

    jobs = []
    for pattern, members in sorted(groups.items()):
        members = sorted(members)
        K = len(members)
        if K < 2:
            print("skip '%s': apenas %d run (sem par)" % (pattern, K))
            continue
        for k in range(K):
            probe = members[k]
            refs  = members[:k] + members[k + 1:]
            jobs.append((pattern, probe, refs))

    print("total de jobs (LOO): %d" % len(jobs))
    for i, (pattern, probe, refs) in enumerate(jobs, 1):
        print("[%d/%d] %s: probe=%s refs=%s" % (i, len(jobs), pattern, probe, refs))
        subprocess.run([sys.executable, SCRIPT, probe] + refs, check=True)


if __name__ == "__main__":
    main()
