"""
Flow signal pipeline v11 — spatial-domain panorama and preprocessed signals.

Pipeline steps:
  1. Load raw camera frames and measure the pixel displacement between each
     consecutive pair (inter-frame shift) via cross-correlation.
  2. Reconstruct a panoramic image by stitching the center strip of every
     frame at its measured position (center-slice approach).
  3. Extract two 1-D grayscale signals from the panorama:
       - full-height : column mean across all 80 rows
       - cropped     : column mean across rows 20–59 (flow core only)
  4. Apply a cumulative preprocessing chain to the cropped signal and plot
     every stage so the effect of each step is visible.

All signals are in the spatial domain: x-axis = tube position in pixels,
1 sample = 1 px of tube travel (no time axis, no sampling frequency).

Usage:
  python3 pipeline_v11.py <folder> [shift_px]

Outputs (written inside <folder>):
  __panorama_v11.png        — reconstructed panoramic image
  __dados_v11.xlsx          — full-height signal  (all 80 rows averaged)
  __dados_cropped_v11.xlsx  — cropped signal      (rows 20–59 averaged)
  __signals_v11.png         — multi-panel plot of all signal stages
"""

import os
import sys
import numpy as np
from PIL import Image
import openpyxl
import matplotlib.pyplot as plt
from scipy.signal import detrend

# ── Output file names ─────────────────────────────────────────────────────────
PANORAMA_V11        = "__panorama_v11.png"
DADOS_V11           = "__dados_v11.xlsx"
DADOS_CROPPED_V11   = "__dados_cropped_v11.xlsx"
SIGNALS_V11         = "__signals_v11.png"

# Rows of the panorama that represent the flow core (excludes tube walls).
CROP_ROWS = slice(20, 60)


# ── Step 1a: load frame file list ─────────────────────────────────────────────

def load_frames(folder):
    """Return sorted list of raw camera frame filenames found in folder."""
    files = sorted(
        f for f in os.listdir(folder)
        if f.lower().endswith((".jpg", ".jpeg", ".png"))
        and not f.startswith("__")    # skip pipeline outputs
        and not f.startswith("dados")
        and not f.startswith("panorama")
        and not f.startswith("signals")
    )
    if not files:
        print("No image frames found.")
        sys.exit(1)
    return files


# ── Step 1b: measure inter-frame pixel shifts ─────────────────────────────────

def _best_shift(row_a, row_b, max_shift):
    """Return shift s (px) that maximises cross-correlation between two 1-D row profiles.

    The convention is a[n+s] versus b[n]. If frame B is translated to the right
    relative to A, the returned lag is negative. This is a correlation lag,
    not a validated physical displacement or velocity.
    """
    n = len(row_a)
    best_s, best_r = 0, -np.inf
    for s in range(-max_shift, max_shift + 1):
        if s >= 0:
            r = np.corrcoef(row_a[s:], row_b[:n - s])[0, 1]
        else:
            r = np.corrcoef(row_a[:n + s], row_b[-s:])[0, 1]
        if r > best_r:
            best_r, best_s = r, s
    return best_s


def measure_all_shifts(folder, files, row_slice=slice(20, 60)):
    """Measure inter-frame pixel displacement for every consecutive frame pair.

    Uses cross-correlation on the horizontally-averaged row profile of the
    center band (row_slice) — robust against vertical variation (tube walls).
    Returns a list of length len(files); the last entry duplicates its predecessor
    because there is no frame N+1 to compare against.
    """
    shifts = []
    print(f"Measuring per-frame shifts for {len(files)} frames …")
    for i in range(len(files) - 1):
        a = np.array(Image.open(os.path.join(folder, files[i])).convert("L"), dtype=float)
        b = np.array(Image.open(os.path.join(folder, files[i + 1])).convert("L"), dtype=float)
        row_a = a[row_slice, :].mean(axis=0)
        row_b = b[row_slice, :].mean(axis=0)
        shifts.append(_best_shift(row_a, row_b, min(200, len(row_a) // 4)))

    shifts.append(shifts[-1])          # last frame has no successor — reuse previous shift
    # Historical v11 policy: nonpositive lags become one-pixel strips. Such
    # replacements are quality flags, not measured one-pixel displacements.
    # See v16_spatial2_18_71/05_conclusao/Analise adversarial.md.
    shifts = [max(s, 1) for s in shifts]

    arr = np.array(shifts)
    print(f"Per-frame shifts — mean: {arr.mean():.2f} px  median: {np.median(arr):.1f} px  "
          f"std: {arr.std():.2f} px  min: {arr.min()}  max: {arr.max()}")
    return shifts


# ── Step 2: reconstruct panorama ─────────────────────────────────────────────

def _build_canvas(files, shifts, folder):
    """Stitch center strips from all frames into a single panorama array.

    For each frame i, cut a strip of width = shifts[i] centred on the frame's
    horizontal midpoint. Concatenating these strips side-by-side reconstructs
    the tube as the flow moves through the camera field of view.

    A 2-pixel linear feather blend is applied at each seam to suppress hard edges
    caused by slight brightness differences between adjacent frames.
    """
    first = Image.open(os.path.join(folder, files[0]))
    W, H  = first.size
    cx    = W // 2

    # Cut the center strip for each frame
    slices = []
    for i, f in enumerate(files):
        s         = shifts[i]
        col_start = cx - s // 2
        col_end   = col_start + s
        strip = np.array(
            Image.open(os.path.join(folder, f)).convert("L"), dtype=np.float32
        )[:, col_start:col_end]
        slices.append(strip)

    # Concatenate all strips into one canvas
    canvas = np.concatenate(slices, axis=1).copy()

    # Feather blend at each seam to smooth brightness transitions
    seams = np.concatenate([[0], np.cumsum([sl.shape[1] for sl in slices])])
    for i in range(1, len(slices)):
        x       = int(seams[i])
        feather = min(2, min(slices[i - 1].shape[1], slices[i].shape[1]) // 2)
        if feather == 0:
            continue
        ramp = np.linspace(0.0, 1.0, feather)
        prev = slices[i - 1][:, -feather:]
        curr = slices[i][:,  :feather]
        w    = ramp[np.newaxis, :]
        canvas[:, x - feather:x] = prev * (1 - w) + curr * w

    return np.clip(canvas, 0, 255).astype(np.uint8), H


def reconstruct(folder, shift_override=None):
    """Orchestrate shift measurement → canvas build → panorama save.

    Returns (panorama_array, shifts_list).
    """
    files = load_frames(folder)
    first = Image.open(os.path.join(folder, files[0]))
    W, H  = first.size
    print(f"{len(files)} frames | {W}x{H} px")

    if shift_override is not None:
        shifts = [max(shift_override, 1)] * len(files)
    else:
        shifts = measure_all_shifts(folder, files)

    panorama, H = _build_canvas(files, shifts, folder)

    out = os.path.join(folder, PANORAMA_V11)
    Image.fromarray(panorama).save(out)
    print(f"Panorama saved: {out}  ({panorama.shape[1]}x{H} px)")
    return panorama, shifts


# ── Step 3: extract signals from panorama ────────────────────────────────────

def _save_signal_xlsx(path, signal, shift_mean, row_desc):
    """Write a 1-D spatial signal to an xlsx file with a metadata sheet."""
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Sheet"
    ws.append(["position_px", "Voltage_0", "Voltage_1", "Voltage_2"])
    for i, val in enumerate(signal):
        ws.append([i, 0.0, 0.0, float(val)])

    ws_meta = wb.create_sheet("metadata")
    ws_meta.append(["key",              "value"])
    ws_meta.append(["domain",           "spatial"])
    ws_meta.append(["spatial_sampling", "1 px/sample"])
    ws_meta.append(["n_samples",        len(signal)])
    ws_meta.append(["shift_mean_px",    shift_mean])
    ws_meta.append(["rows_averaged",    row_desc])

    wb.save(path)
    print(f"xlsx saved: {path}  ({len(signal)} samples, rows {row_desc})")


def preprocess(raw_signal):
    """Apply preprocessing steps cumulatively to the cropped signal.

    Each step is applied on top of the previous one so the effect of adding
    each operation is visible when the stages are plotted separately.

    Returns a list of (label, signal) tuples — one per stage including raw.
    """
    # Invert so liquid regions are high and gas bubbles are low
    inverted  = 255.0 - raw_signal

    # Remove linear trend (slow drift along the tube caused by illumination gradient)
    detrended = detrend(inverted, type="linear")

    # Subtract the mean so the signal is centred around zero
    zero_mean = detrended - detrended.mean()

    # Divide by the standard deviation so amplitude is in units of σ (z-score)
    z_scored  = zero_mean / zero_mean.std()

    return [
        (f"Raw — rows {CROP_ROWS.start}–{CROP_ROWS.stop - 1}",  raw_signal),
        ("+ Invert (255 − x)",                                    inverted),
        ("+ Linear detrend",                                      detrended),
        ("+ Zero-mean",                                           zero_mean),
        ("+ Z-score (amplitude in σ)",                            z_scored),
    ]


def extract_signals(folder, panorama, shifts):
    """Extract signals from the panorama, save xlsx files, and build preprocessing chain.

    Returns (sig_full, preprocessing_steps) where preprocessing_steps is a list of
    (label, signal) tuples produced by preprocess().
    """
    shift_mean = float(np.mean(shifts))

    # 1. Extract raw signals: column mean collapses the 2-D panorama to 1-D
    sig_full    = panorama.mean(axis=0).astype(float)           # all 80 rows
    sig_cropped = panorama[CROP_ROWS, :].mean(axis=0).astype(float)  # flow core only

    # 2. Save both raw signals to xlsx
    _save_signal_xlsx(
        os.path.join(folder, DADOS_V11),
        sig_full, shift_mean, f"0–{panorama.shape[0] - 1} (all)"
    )
    _save_signal_xlsx(
        os.path.join(folder, DADOS_CROPPED_V11),
        sig_cropped, shift_mean, f"{CROP_ROWS.start}–{CROP_ROWS.stop - 1}"
    )

    # 3. Build the preprocessing chain on the cropped signal
    preprocessing_steps = preprocess(sig_cropped)

    return sig_full, preprocessing_steps


# ── Step 4: plot all signal stages ───────────────────────────────────────────

def plot(folder, sig_full, preprocessing_steps):
    """Plot the full-height signal and each preprocessing stage in separate panels."""
    name   = os.path.basename(folder.rstrip("/"))
    n_rows = 1 + len(preprocessing_steps)
    fig, axes = plt.subplots(n_rows, 1, figsize=(14, 3 * n_rows), sharex=True)

    # Panel 0: full-height raw signal — shows the complete tube cross-section
    axes[0].plot(np.arange(len(sig_full)), sig_full, color="tab:blue", linewidth=0.6)
    axes[0].set_title(f"Full-height signal (all rows averaged) — {name}")
    axes[0].set_ylabel("Grayscale (0–255)")
    axes[0].grid(True, linestyle="--", linewidth=0.4, alpha=0.6)

    # Panels 1…N: each preprocessing stage applied on the cropped signal
    colors = ["tab:orange", "tab:green", "tab:red", "tab:purple"]
    for idx, (label, sig) in enumerate(preprocessing_steps):
        ax = axes[1 + idx]
        ax.plot(np.arange(len(sig)), sig, color=colors[idx % len(colors)], linewidth=0.6)
        ax.set_title(label)
        ax.set_ylabel("Amplitude")
        ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.6)

    axes[-1].set_xlabel("Position (px)")
    fig.tight_layout()
    out = os.path.join(folder, SIGNALS_V11)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Plot saved: {out}")


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    folder = sys.argv[1] if len(sys.argv) > 1 else "."
    shift  = int(sys.argv[2]) if len(sys.argv) > 2 else None

    # Step 1+2: reconstruct the panoramic image from raw frames
    panorama, shifts = reconstruct(folder, shift)

    # Step 3: extract signals from the panorama and save xlsx files
    sig_full, preprocessing_steps = extract_signals(folder, panorama, shifts)

    # Step 4: plot all signal stages
    plot(folder, sig_full, preprocessing_steps)
