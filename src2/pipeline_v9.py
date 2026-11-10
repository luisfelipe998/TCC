"""
Flow signal pipeline v9 — time-domain per-frame signal (precompute stage).

Counterpart to pipeline_v11.py. Where v11 stitches the 501 camera frames into a
*spatial* panorama (x-axis = tube position in px), v9 collapses each frame to a
single scalar and stacks them into a *time* series (x-axis = time in seconds,
1 sample = 1 camera frame). The time sampling rate is the camera frame rate.

Pipeline steps (per run folder in data_processed):
  1. Load the raw camera frames (ImgA*.jpg).
  2. For each frame, take the mean grayscale of the flow-core rows (CROP_ROWS,
     rows 20-59 — same core band v11 uses). One scalar per frame.
  3. Read fps / frame count from _config.xsv; derive dt, T, fs, f_min, Nyquist.
  4. Apply a cumulative preprocessing chain (linear detrend -> zero-mean ->
     z-score) to the raw per-frame series. No brightness inversion.
  5. Store the preprocessed time series plus timing metadata to
     __timeseries_v9.xlsx inside the run folder, so the wavelet analysis reads
     xlsx only and decodes the JPGs exactly once per run.

The stored Voltage_2 column is the analysis signal (matching v11/v13, where
Voltage_2 is the optical-sensor-analogue channel the wavelet scripts read).
Voltage_0/Voltage_1 are written as zeros to keep the column layout identical to
the v11/v13 xlsx files.

Usage:
  python3 pipeline_v9.py <folder>        # one run
  python3 pipeline_v9.py --all           # every run under data_processed

Outputs (written inside each <folder>):
  __timeseries_v9.xlsx   — time-domain per-frame signal + timing metadata
  __signals_v9.png       — multi-panel plot of each preprocessing stage
"""

import os
import sys
import glob
import numpy as np
from PIL import Image
import openpyxl
import matplotlib.pyplot as plt
from scipy.signal import detrend, resample_poly
from math import gcd

# ── Paths ─────────────────────────────────────────────────────────────────────
HERE      = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.join(HERE, "data_processed")

# ── Output file names ─────────────────────────────────────────────────────────
TIMESERIES_V9 = "__timeseries_v9.xlsx"
SIGNALS_V9    = "__signals_v9.png"

# Rows of each frame that represent the flow core (excludes tube walls).
# Same band pipeline_v11.py averages spatially; here we average it per frame.
CROP_ROWS = slice(20, 60)

# Common sampling rate all runs are resampled to. The dataset was filmed at
# 4000/5000/7000 fps; cross-pattern pairs mixed rates on one axis, so every run
# is anti-alias downsampled to the SLOWEST rate (4000 fps). Downsampling to the
# common floor discards no analyzed signal (<1.2% of power sits above 1000 Hz,
# and the axis caps at 1000 Hz anyway), keeps every sample real (no interpolated
# invention), and preserves clip DURATION T — so the per-run COI floor, which is
# set by T not by fps, is unchanged (10.6/13.2/18.5 Hz for the 4000/5000/7000 fps
# groups). Upsampling would instead fabricate samples and cannot add real content
# above each run's original Nyquist.
TARGET_FPS = 4000.0


def resample_to_target(raw_signal, native_fps, target_fps=TARGET_FPS):
    """Anti-alias resample a per-frame series from native_fps to target_fps.

    Uses scipy.signal.resample_poly (polyphase FIR, built-in anti-alias low-pass),
    with the up/down ratio reduced to lowest terms. A no-op when the rates match.
    Duration is preserved; the sample count changes (e.g. 501 -> 287 for 7000->4000).
    """
    if native_fps == target_fps:
        return raw_signal
    up, down = int(round(target_fps)), int(round(native_fps))
    g = gcd(up, down)
    return resample_poly(raw_signal, up // g, down // g)



# ── Step 1: load frame file list ──────────────────────────────────────────────

def load_frames(folder):
    """Return sorted list of raw camera frame filenames found in folder.

    Mirrors pipeline_v11.load_frames: skips pipeline output files so only the
    raw ImgA*.jpg acquisition frames are returned.
    """
    files = sorted(
        f for f in os.listdir(folder)
        if f.lower().endswith((".jpg", ".jpeg", ".png"))
        and not f.startswith("__")
        and not f.startswith("dados")
        and not f.startswith("panorama")
        and not f.startswith("signals")
    )
    if not files:
        print("No image frames found in %s." % folder)
        sys.exit(1)
    return files


# ── Step 2: collapse each frame to one scalar ─────────────────────────────────

def build_raw_series(folder, files, row_slice=CROP_ROWS):
    """Collapse every frame to the mean grayscale of its flow-core rows.

    Returns a 1-D array of length len(files): one sample per camera frame, in
    acquisition order. This is the time-domain analogue of v11's spatial column
    mean — there we average columns of a panorama, here we average the flow-core
    rows *and* all columns of each individual frame to get its instantaneous
    brightness.
    """
    series = np.empty(len(files), dtype=float)
    for i, f in enumerate(files):
        img = np.array(Image.open(os.path.join(folder, f)).convert("L"), dtype=float)
        series[i] = img[row_slice, :].mean()
    return series


# ── Step 3: timing from the camera config ─────────────────────────────────────

def read_config(folder):
    """Parse _config.xsv into a flat dict of key -> value (strings)."""
    path = os.path.join(folder, "_config.xsv")
    cfg = {}
    with open(path, encoding="latin-1") as fh:
        for line in fh:
            if "=" in line and not line.lstrip().startswith("["):
                k, _, v = line.partition("=")
                cfg[k.strip()] = v.strip()
    return cfg


def timing_from_config(cfg, n_frames, fps_override=None):
    """Derive time-domain timing parameters.

    fs is the camera frame rate 'Rate' (fps) unless fps_override is given (used
    after resampling to a common rate). Everything else follows: dt = 1/fs,
    T = n_frames/fs, f_min = 1/T, Nyquist = fs/2. When called with the RESAMPLED
    length and the resample target rate, T reproduces the original clip duration.
    """
    fps = float(cfg["Rate"])
    fs  = float(fps_override) if fps_override is not None else fps
    dt  = 1.0 / fs
    T   = n_frames / fs
    return {
        "fps":         fps,
        "fs_hz":       fs,
        "dt_s":        dt,
        "T_s":         T,
        "f_min_run_hz": 1.0 / T,
        "nyquist_hz":  fs / 2.0,
        "n_frames":    n_frames,
    }


# ── Step 4: preprocessing chain (identical operations to pipeline_v11) ────────

def preprocess(raw_signal):
    """Apply preprocessing steps cumulatively.

    No brightness inversion (v11 inverted so liquid reads high; the wavelet works
    on magnitude and the sign is irrelevant, so v9 keeps the raw grayscale).
    Returns a list of (label, signal) tuples — one per stage including raw — and
    the final z-scored signal is the last entry.
    """
    detrended = detrend(raw_signal, type="linear")   # remove slow drift
    zero_mean = detrended - detrended.mean()          # centre at zero
    z_scored  = zero_mean / zero_mean.std(ddof=1)     # amplitude in sigma

    return [
        ("Raw — flow-core mean per frame", raw_signal),
        ("+ Linear detrend",               detrended),
        ("+ Zero-mean",                    zero_mean),
        ("+ Z-score (amplitude in sigma)", z_scored),
    ]


# ── Step 5: persist the time series ───────────────────────────────────────────

def save_timeseries_xlsx(path, z_signal, timing):
    """Write the preprocessed time series to xlsx with a timing metadata sheet.

    Column layout mirrors the v11/v13 xlsx (position/Voltage_0/1/2) so downstream
    loaders stay uniform; here the first column is time_s and Voltage_2 holds the
    analysis signal.
    """
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Sheet"
    ws.append(["time_s", "Voltage_0", "Voltage_1", "Voltage_2"])
    dt = timing["dt_s"]
    for i, val in enumerate(z_signal):
        ws.append([i * dt, 0.0, 0.0, float(val)])

    ws_meta = wb.create_sheet("metadata")
    ws_meta.append(["key", "value"])
    ws_meta.append(["domain",       "time"])
    ws_meta.append(["fs_hz",        timing["fs_hz"]])
    ws_meta.append(["dt_s",         timing["dt_s"]])
    ws_meta.append(["T_s",          timing["T_s"]])
    ws_meta.append(["n_frames",     timing["n_frames"]])
    ws_meta.append(["f_min_run_hz", timing["f_min_run_hz"]])
    ws_meta.append(["nyquist_hz",   timing["nyquist_hz"]])
    ws_meta.append(["rows_averaged", "%d-%d" % (CROP_ROWS.start, CROP_ROWS.stop - 1)])

    wb.save(path)
    print("xlsx saved: %s  (%d samples, fs=%g Hz, T=%.4f s, f_min=%.2f Hz, Nyq=%g Hz)"
          % (path, len(z_signal), timing["fs_hz"], timing["T_s"],
             timing["f_min_run_hz"], timing["nyquist_hz"]))


def plot_stages(folder, steps, dt):
    """Plot each preprocessing stage against time (seconds)."""
    name   = os.path.basename(folder.rstrip("/"))
    n_rows = len(steps)
    t      = np.arange(len(steps[0][1])) * dt
    fig, axes = plt.subplots(n_rows, 1, figsize=(14, 2.6 * n_rows), sharex=True)
    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]
    for idx, (label, sig) in enumerate(steps):
        ax = axes[idx]
        ax.plot(t, sig, color=colors[idx % len(colors)], linewidth=0.8)
        ax.set_title(label)
        ax.set_ylabel("Amplitude")
        ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.6)
    axes[-1].set_xlabel("Time (s)")
    fig.suptitle("Time-domain per-frame signal stages — %s" % name, y=1.001)
    fig.tight_layout()
    out = os.path.join(folder, SIGNALS_V9)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("plot saved: %s" % out)


# ── Orchestration ─────────────────────────────────────────────────────────────

def process_run(folder):
    """Precompute and store the time-domain signal for one run folder."""
    name  = os.path.basename(folder.rstrip("/"))
    files = load_frames(folder)
    cfg   = read_config(folder)

    n_cfg = int(cfg.get("Frames", len(files)))
    if n_cfg != len(files):
        print("  note: config Frames=%d but found %d image files; using file count"
              % (n_cfg, len(files)))

    native_fps = float(cfg["Rate"])
    raw_native = build_raw_series(folder, files)
    # Unify sampling rate: anti-alias downsample every run to TARGET_FPS so all
    # pairs (including cross-fps ones) share one time/frequency axis honestly.
    raw = resample_to_target(raw_native, native_fps)

    # Timing reflects the RESAMPLED signal: fs = TARGET_FPS, length changed,
    # duration T preserved (so the per-run COI floor is unchanged).
    timing = timing_from_config(cfg, len(raw), fps_override=TARGET_FPS)

    print("%s: %d frames @ %g fps -> %d samples @ %g fps | T=%.4fs"
          % (name, len(files), native_fps, len(raw), timing["fs_hz"], timing["T_s"]))

    steps = preprocess(raw)
    z     = steps[-1][1]

    save_timeseries_xlsx(os.path.join(folder, TIMESERIES_V9), z, timing)
    plot_stages(folder, steps, timing["dt_s"])


def list_run_folders():
    """Every data_processed subfolder that contains a camera config."""
    runs = []
    for d in sorted(os.listdir(DATA_ROOT)):
        path = os.path.join(DATA_ROOT, d)
        if os.path.isdir(path) and os.path.exists(os.path.join(path, "_config.xsv")):
            runs.append(path)
    return runs


if __name__ == "__main__":
    args = sys.argv[1:]
    if not args:
        print(__doc__)
        sys.exit(1)

    if args[0] == "--all":
        folders = list_run_folders()
        print("processing %d runs under %s" % (len(folders), DATA_ROOT))
        for f in folders:
            process_run(f)
    else:
        folder = args[0]
        if not os.path.isdir(folder):
            folder = os.path.join(DATA_ROOT, folder)
        process_run(folder)
