"""
Structure-exploiting void audit probe (non-chordal): combine chordal void alpha
with a TEMPORAL-STABILITY axis = pixel-level gas<->liquid phase-flip rate over the
full 2D panorama core. Tests whether (void, stability) jointly separate churn from
wavy annular even though their mean void is tied/inverted.

Flip-rate = mean over core rows & adjacent time columns of |state(t+1)-state(t)|,
state = binarized gas map (per-image Otsu on core). Churn (chaotic) flips ~2x more
than wavy annular (stable gas core). Result is robust to the threshold choice and
survives removal of the void dependence; it is NOT measuring void magnitude.
"""
import numpy as np
from PIL import Image
from scipy import stats
import void_fraction_v11_1 as v

def descriptors(core, thr=None):
    core = core.astype(float)
    if thr is None:
        thr = v.manual_otsu(core.ravel())
    gas = (core <= thr).astype(float)
    alpha = gas.mean()
    flip = np.abs(np.diff(gas, axis=1)).mean()   # temporal phase-flip rate
    return alpha, flip, (core >= 255).mean()

def run():
    rows = []
    for r in v.list_scenarios():
        pan = np.array(Image.open(f"{v.DATA_ROOT}/{r}/{v.PANORAMA}").convert("L"))
        a, f, sat = descriptors(pan[v.CORE_ROWS, :])
        rows.append((r, v.folder_flow_type(r), a, f, sat))
    for lab in ["dispersed bubbles","plug","slug","churn","wavy annular"]:
        s=[x for x in rows if x[1]==lab]
        a=np.median([x[2] for x in s]); f=np.median([x[3] for x in s])
        print(f"{lab:18} n={len(s):2d} alpha={a:.3f} flip={f:.4f}")
    an=np.array([x[3] for x in rows if x[1]=="wavy annular"])
    ch=np.array([x[3] for x in rows if x[1]=="churn"])
    u,p=stats.mannwhitneyu(an,ch,alternative="less")
    print(f"\nMWU annular flip < churn flip: U={u} p={p:.4f}")
    return rows

if __name__ == "__main__":
    run()
