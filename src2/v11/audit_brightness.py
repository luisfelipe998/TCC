import os, re, glob
import numpy as np
from PIL import Image

HERE = "/Users/i521862/Documents/personal/unisinos/TCC/src2/v11"
DATA_ROOT = "/Users/i521862/Documents/personal/unisinos/TCC/src2/data_processed"
PANORAMA = "__panorama_v11.png"
CORE = slice(20, 60)
RUN_RE = re.compile(r"^m0\dq\d+_\d+$")
EXCLUDE = {"m03q80_001", "m07q01_001"}
PURE = {"plug", "slug", "churn", "wavy annular"}
GROUP = {"plug": "plug", "slug": "slug", "churn": "churn", "wavy annular": "annular"}

def manual_otsu(vals):
    hist = np.histogram(vals, bins=256, range=(0, 256))[0].astype(float)
    total = hist.sum(); sum_all = np.dot(np.arange(256), hist)
    wB = 0.0; sumB = 0.0; best = -1.0; thr = 0
    for t in range(256):
        wB += hist[t]
        if wB == 0: continue
        wF = total - wB
        if wF == 0: break
        sumB += t * hist[t]
        mB = sumB / wB; mF = (sum_all - sumB) / wF
        between = wB * wF * (mB - mF) ** 2
        if between > best: best = between; thr = t
    return thr

def flow_type(name):
    for p in sorted(glob.glob(os.path.join(DATA_ROOT, name, "*.txt"))):
        stem = os.path.splitext(os.path.basename(p))[0].lstrip("_")
        if stem != "dados": return stem
    return name

def metrics_from_core(core):
    thr = manual_otsu(core.ravel())
    aj = (core <= thr).mean(axis=0)
    return thr, float(aj.mean()), float(aj.std())

runs = []
for d in sorted(os.listdir(DATA_ROOT)):
    if not RUN_RE.match(d) or d in EXCLUDE: continue
    if not os.path.exists(os.path.join(DATA_ROOT, d, PANORAMA)): continue
    lab = flow_type(d)
    if lab not in PURE: continue
    runs.append((d, lab))

rows = []
for name, lab in runs:
    pan = np.array(Image.open(os.path.join(DATA_ROOT, name, PANORAMA)).convert("L"))
    core = pan[CORE, :].astype(float)
    thr, alpha, astd = metrics_from_core(core)
    cm = float(core.mean()); cstd = float(core.std())
    sat = float((core >= 255).mean())
    # normalization 1: z-score then rescale to common mean/std (mu=63,sigma=40 target typical)
    z = (core - cm) / (cstd + 1e-9)
    core_z = np.clip(z * 40.0 + 63.0, 0, 255)
    thr_z, alpha_z, astd_z = metrics_from_core(core_z)
    rows.append(dict(name=name, group=GROUP[lab], cm=cm, cstd=cstd, sat=sat,
                     thr=thr, alpha=alpha, astd=astd, astd_z=astd_z, alpha_z=alpha_z,
                     N=pan.shape[1]))

# histogram matching to a reference: pick the median-brightness normal run's core as reference
def hist_match(src, ref_sorted):
    s = src.ravel()
    order = np.argsort(s)
    ranks = np.empty_like(order); ranks[order] = np.arange(len(s))
    q = ranks / (len(s) - 1)
    idx = np.clip((q * (len(ref_sorted) - 1)).round().astype(int), 0, len(ref_sorted) - 1)
    return ref_sorted[idx].reshape(src.shape)

# reference = pooled core of all non-saturated (sat<0.20) runs
ref_pool = []
for name, lab in runs:
    pan = np.array(Image.open(os.path.join(DATA_ROOT, name, PANORAMA)).convert("L"))
    core = pan[CORE, :].astype(float)
    if (core >= 255).mean() <= 0.20:
        ref_pool.append(core.ravel())
ref_sorted = np.sort(np.concatenate(ref_pool))
# subsample ref for speed
ref_sub = ref_sorted[np.linspace(0, len(ref_sorted) - 1, 200000).astype(int)]
ref_sub = np.sort(ref_sub)

for r in rows:
    pan = np.array(Image.open(os.path.join(DATA_ROOT, r["name"], PANORAMA)).convert("L"))
    core = pan[CORE, :].astype(float)
    matched = hist_match(core, ref_sub)
    thr_h, alpha_h, astd_h = metrics_from_core(matched)
    r["astd_h"] = astd_h; r["alpha_h"] = alpha_h

def med(vals): return float(np.median(vals)) if vals else float("nan")

print("=== PER-RUN ===")
print(f"{'run':14}{'grp':8}{'cm':>7}{'cstd':>7}{'sat%':>7}{'thr':>5}{'alpha':>7}{'astd':>7}{'astd_z':>8}{'astd_h':>8}")
for r in sorted(rows, key=lambda x: (x["group"], x["name"])):
    print(f"{r['name']:14}{r['group']:8}{r['cm']:7.1f}{r['cstd']:7.1f}{r['sat']*100:7.1f}{r['thr']:5d}"
          f"{r['alpha']:7.3f}{r['astd']:7.3f}{r['astd_z']:8.3f}{r['astd_h']:8.3f}")

def group_med(key, g):
    return med([r[key] for r in rows if r["group"] == g])

print("\n=== PER-PATTERN MEDIAN astd ===")
for g in ["plug", "slug", "churn", "annular"]:
    n = sum(1 for r in rows if r["group"] == g)
    print(f"{g:8} n={n:2d}  astd={group_med('astd',g):.3f}  astd_z={group_med('astd_z',g):.3f}  astd_h={group_med('astd_h',g):.3f}  cm={group_med('cm',g):.1f}")

# correlations across all runs
import numpy as np
def corr(a, b):
    a = np.array(a); b = np.array(b)
    return float(np.corrcoef(a, b)[0, 1])
allast = [r["astd"] for r in rows]; allcm = [r["cm"] for r in rows]; allsat = [r["sat"] for r in rows]
print("\n=== CORRELATIONS (all runs) ===")
print(f"astd vs core_mean : r={corr(allast, allcm):.3f}")
print(f"astd vs sat_frac  : r={corr(allast, allsat):.3f}")
# within annular
ann = [r for r in rows if r["group"] == "annular"]
print(f"\n=== WITHIN ANNULAR (n={len(ann)}) ===")
print(f"astd vs core_mean : r={corr([r['astd'] for r in ann], [r['cm'] for r in ann]):.3f}")
print(f"astd vs sat_frac  : r={corr([r['astd'] for r in ann], [r['sat'] for r in ann]):.3f}")
chu = [r for r in rows if r["group"] == "churn"]
print(f"\n=== WITHIN CHURN (n={len(chu)}) ===")
print(f"astd vs core_mean : r={corr([r['astd'] for r in chu], [r['cm'] for r in chu]):.3f}")

# bright vs normal split within annular
print("\n=== ANNULAR bright(cm>100) vs normal ===")
bright = [r for r in ann if r["cm"] > 100]
normal = [r for r in ann if r["cm"] <= 100]
print(f"bright (n={len(bright)}): " + ", ".join(f"{r['name']}(cm{r['cm']:.0f},astd{r['astd']:.3f})" for r in bright))
print(f"normal (n={len(normal)}): " + ", ".join(f"{r['name']}(cm{r['cm']:.0f},astd{r['astd']:.3f})" for r in normal))
print(f"bright astd median = {med([r['astd'] for r in bright]):.3f}")
print(f"normal astd median = {med([r['astd'] for r in normal]):.3f}")
print(f"churn astd median  = {group_med('astd','churn'):.3f}")
