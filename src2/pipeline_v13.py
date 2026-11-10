"""
Flow signal pipeline v13 — two-signal comparison from v11 panorama.

Reads the existing __panorama_v11.png from a run folder and extracts two
1-D grayscale signals:
  - Signal A: column mean across rows 20–59 (flow core average)
  - Signal B: single column values from row 40 (tube center only)

Both signals share the same x-axis (panorama column = tube position in px).

Usage:
  python3 pipeline_v13.py <folder>

Outputs (written inside <folder>):
  __signals_v13.png        — two-panel comparison plot
  __dados_center_v13.xlsx  — center-row signal (row 40) as xlsx
"""

import os
import sys
import numpy as np
from PIL import Image
import openpyxl
import matplotlib.pyplot as plt

PANORAMA_V11      = "__panorama_v11.png"
SIGNALS_V13       = "__signals_v13.png"
DADOS_CENTER_V13  = "__dados_center_v13.xlsx"

CROP_ROWS  = slice(20, 60)
CENTER_ROW = 40


def extract_signals(folder):
    path = os.path.join(folder, PANORAMA_V11)
    if not os.path.exists(path):
        print(f"Error: {path} not found. Run pipeline_v11.py first.")
        sys.exit(1)
    panorama    = np.array(Image.open(path).convert("L"), dtype=float)
    sig_cropped = panorama[CROP_ROWS, :].mean(axis=0)
    sig_center  = panorama[CENTER_ROW, :]
    return sig_cropped, sig_center


def save_center_xlsx(folder, sig_center):
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Sheet"
    ws.append(["position_px", "Voltage_0", "Voltage_1", "Voltage_2"])
    for i, val in enumerate(sig_center):
        ws.append([i, 0.0, 0.0, float(val)])

    ws_meta = wb.create_sheet("metadata")
    ws_meta.append(["key",              "value"])
    ws_meta.append(["domain",           "spatial"])
    ws_meta.append(["spatial_sampling", "1 px/sample"])
    ws_meta.append(["n_samples",        len(sig_center)])
    ws_meta.append(["source_row",       CENTER_ROW])
    ws_meta.append(["source",          f"panorama row {CENTER_ROW} (no averaging)"])

    path = os.path.join(folder, DADOS_CENTER_V13)
    wb.save(path)
    print(f"Saved: {path}  ({len(sig_center)} samples, row {CENTER_ROW})")


def plot(folder, sig_cropped, sig_center):
    name = os.path.basename(folder.rstrip("/"))
    x    = np.arange(len(sig_cropped))

    fig, axes = plt.subplots(2, 1, figsize=(14, 4), sharex=True)

    axes[0].plot(x, sig_cropped, color="tab:orange", linewidth=0.6)
    axes[0].set_title(f"Raw — rows {CROP_ROWS.start}–{CROP_ROWS.stop - 1}  ({name})")
    axes[0].set_ylabel("Grayscale (0–255)")
    axes[0].grid(True, linestyle="--", linewidth=0.4, alpha=0.6)

    axes[1].plot(x, sig_center, color="tab:blue", linewidth=0.6)
    axes[1].set_title(f"Center row (px {CENTER_ROW})  ({name})")
    axes[1].set_ylabel("Grayscale (0–255)")
    axes[1].set_xlabel("Position (px)")
    axes[1].grid(True, linestyle="--", linewidth=0.4, alpha=0.6)

    fig.tight_layout()
    out = os.path.join(folder, SIGNALS_V13)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    folder = sys.argv[1] if len(sys.argv) > 1 else "."
    sig_cropped, sig_center = extract_signals(folder)
    save_center_xlsx(folder, sig_center)
    plot(folder, sig_cropped, sig_center)
