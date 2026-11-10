import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, MaxNLocator, EngFormatter
from matplotlib.gridspec import GridSpec
import pandas as pd
import pywt
from scipy.signal import detrend as scipy_detrend

# ══════════════════════════════════════════════════════════════════════════════
# wavelet_v9 — TIME-DOMAIN wavelet cross-spectrum / coherence.
#
# Time-domain counterpart of v8_spatial_cropped_and_center/wavelet_v8.py. The
# wavelet machinery (biased auto/cross-correlograms -> CWT -> WPS/WCS/WCO) is
# identical; only the domain changes:
#   - signal      : per-frame time series from __timeseries_v9.xlsx (built by
#                   pipeline_v9.py), one sample per camera frame.
#   - sampling    : sampling_period = dt (seconds), read from run metadata.
#   - freq axis   : Hz, on a SHARED baseline [FREQ_MIN .. FREQ_MAX] common to
#                   every run (least-restrictive f_min .. smallest Nyquist), so
#                   all pairs are plotted on one comparable scale.
#   - lag axis    : seconds (lag_samples * dt).
#   - COI         : time-domain cone of influence; the low-frequency band a run
#                   cannot resolve (its own 1/T above a plotted scale, or the
#                   correlogram-edge e-folding) is shaded and excluded.
# ══════════════════════════════════════════════════════════════════════════════

# ── Signal source ─────────────────────────────────────────────────────────────
HERE      = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.join(HERE, "../", "data_processed")
OS_COL    = "Voltage_2"
TIMESERIES = "__timeseries_v9.xlsx"

# ── Shared frequency axis (Hz) ────────────────────────────────────────────────
# FREQ_MIN = least-restrictive per-run f_min = 1/T_max (4000 fps, T=0.1253 s).
#   Runs with shorter T (5000/7000 fps) cannot reach this floor -> COI-masked.
# FREQ_MAX = 1000 Hz. Well below the smallest per-run Nyquist (2000 Hz, 4000 fps)
#   and above all observed flow energy: across every pattern measured, ~0% of
#   WPS power lies above 1000 Hz (peaks 22-59 Hz), so the 1000-2000 Hz band was
#   empty and is cropped for readability. No data is filtered — the axis stops at
#   1000 Hz; scales above it are simply not computed.
FREQ_MIN         = 1.0 / (501 / 4000.0)          # 7.984 Hz
FREQ_MAX         = 1000.0                          # Hz (axis cap; < min Nyquist 2000)
FREQ_TICKS       = [7.984, 10, 20, 50, 100, 200, 500, 1000]
FREQ_TICK_LABELS = ["8 Hz", "10 Hz", "20 Hz", "50 Hz", "100 Hz",
                    "200 Hz", "500 Hz", "1000 Hz"]

N_SCALES = 200

# ── COI parameters ────────────────────────────────────────────────────────────
# e-folding factor for the cmor cone of influence (same convention as v8).
COI_EFOLD_FACTOR = np.sqrt(2)
COI_FILL_COLOR   = "0.35"
COI_FILL_ALPHA   = 0.72
COI_LINE_COLOR   = "white"
COI_LINE_WIDTH   = 1.0
COI_LINE_STYLE   = "--"

WAVELET = "cmor2.0-0.9549"


# ── Helpers ───────────────────────────────────────────────────────────────────

def _read_meta(path):
    """Return the timing metadata dict from a __timeseries_v9.xlsx file."""
    md = pd.read_excel(path, sheet_name="metadata")
    d  = dict(zip(md["key"].astype(str), md["value"]))
    return {
        "fs_hz":        float(d["fs_hz"]),
        "dt_s":         float(d["dt_s"]),
        "T_s":          float(d["T_s"]),
        "f_min_run_hz": float(d["f_min_run_hz"]),
        "nyquist_hz":   float(d["nyquist_hz"]),
    }


def _load_and_preprocess(path):
    """Load the Voltage_2 time series and re-standardise (z-score).

    pipeline_v9 already z-scores; we detrend + z-score again here so the wavelet
    input is guaranteed zero-mean unit-variance regardless of the stored stage,
    matching wavelet_v8._load_and_preprocess exactly.
    """
    df = pd.read_excel(path, engine="openpyxl", sheet_name="Sheet")
    if OS_COL not in df.columns:
        raise KeyError("Canal '%s' ausente em %s. Colunas: %s"
                       % (OS_COL, path, list(df.columns)))
    sig = pd.to_numeric(df[OS_COL], errors="coerce").to_numpy(float)
    sig = sig[~np.isnan(sig)]
    sig = scipy_detrend(sig, type="linear")
    sig = sig - sig.mean()
    sig = sig / sig.std(ddof=1)
    return sig


def load_signal(name):
    """Load the precomputed time series and its timing metadata for a run."""
    path = os.path.join(DATA_ROOT, name, TIMESERIES)
    if not os.path.exists(path):
        raise FileNotFoundError(
            "%s missing — run pipeline_v9.py first (precompute stage)." % path)
    return _load_and_preprocess(path), _read_meta(path)


def folder_flow_type(name):
    folder = os.path.join(DATA_ROOT, name)
    for p in sorted(glob.glob(os.path.join(folder, "*.txt"))):
        stem = os.path.splitext(os.path.basename(p))[0].lstrip("_")
        if stem != "dados":
            return stem
    return name


def coi_freq_floor(N, dt, f_min_run):
    """Lowest frequency this run can trust, as a Hz floor for COI shading.

    Two limits combine (the more restrictive wins, i.e. the higher frequency):
      1. Record length: a run resolves at best its own f_min_run = 1/T. On the
         shared axis, everything below f_min_run is unreachable for this run.
      2. Correlogram edge e-folding: the cmor cone limits the largest usable
         period near the signal edge. The bilateral correlogram has half-length
         (N-1); the e-folding period there is  P_edge = 2*(N-1)*dt / (efold*fc*2)
         expressed via the wavelet centre frequency fc, giving a frequency floor
         f_edge = 1 / P_edge_max.

    Returns a frequency (Hz); the band below it is shaded/excluded.
    """
    fc = pywt.central_frequency(WAVELET)
    # Largest period the correlogram edge supports (samples -> seconds).
    lam_samples = (2 * N - 1) / (2 * COI_EFOLD_FACTOR * fc)   # in samples
    period_edge = lam_samples * dt                             # seconds
    f_edge      = 1.0 / period_edge if period_edge > 0 else FREQ_MIN
    # Record-length limit and edge limit: take the higher (more restrictive).
    return max(f_min_run, f_edge)


def apply_freq_axis(ax, coi_freq_floor_hz=None):
    """Log Hz y-axis on the shared range; shade the sub-floor COI band."""
    ax.set_yscale("log")
    ax.set_ylim(FREQ_MIN, FREQ_MAX)
    ticks, labels = list(FREQ_TICKS), list(FREQ_TICK_LABELS)

    if coi_freq_floor_hz is not None and FREQ_MIN < coi_freq_floor_hz < FREQ_MAX:
        ax.axhspan(FREQ_MIN, coi_freq_floor_hz,
                   color=COI_FILL_COLOR, alpha=COI_FILL_ALPHA, zorder=5)
        ax.axhline(coi_freq_floor_hz, color=COI_LINE_COLOR, linewidth=COI_LINE_WIDTH,
                   linestyle=COI_LINE_STYLE, zorder=6)
        # Drop any fixed tick within 8% of the COI tick so labels don't overprint
        # on the log axis (e.g. a 19 Hz COI floor colliding with the 20 Hz tick).
        keep = [(t, l) for t, l in zip(ticks, labels)
                if abs(np.log10(t) - np.log10(coi_freq_floor_hz)) > np.log10(1.08)]
        ticks  = [t for t, _ in keep] + [coi_freq_floor_hz]
        labels = [l for _, l in keep] + ["%.0f Hz (COI)" % coi_freq_floor_hz]

    ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.yaxis.set_major_formatter(FixedFormatter(labels))
    ax.yaxis.set_minor_locator(FixedLocator([]))


# ── Signals figure ────────────────────────────────────────────────────────────

def plot_signals(signal1, signal2, dt, name1, flow1, name2, flow2):
    """3-panel time-domain signals figure (x-axis in seconds)."""
    N = min(len(signal1), len(signal2))
    s1, s2 = signal1[:N], signal2[:N]
    t = np.arange(N) * dt

    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    fig.suptitle("%s (%s)  vs  %s (%s)  [time domain]"
                 % (name1, flow1, name2, flow2), fontsize=11)

    axes[0].plot(t, s1, color="tab:blue", linewidth=0.8)
    axes[0].set_title("Signal 1 — %s (%s)" % (name1, flow1))
    axes[0].set_ylabel("Amplitude (sigma)")
    axes[0].grid(True)

    axes[1].plot(t, s2, color="tab:orange", linewidth=0.8)
    axes[1].set_title("Signal 2 — %s (%s)" % (name2, flow2))
    axes[1].set_ylabel("Amplitude (sigma)")
    axes[1].grid(True)

    axes[2].plot(t, s1, color="tab:blue",   linewidth=0.8, label="%s (%s)" % (name1, flow1))
    axes[2].plot(t, s2, color="tab:orange", linewidth=0.8, alpha=0.7, label="%s (%s)" % (name2, flow2))
    axes[2].set_title("Both signals (z-scored)")
    axes[2].set_ylabel("Amplitude (sigma)")
    axes[2].set_xlabel("Time (s)")
    axes[2].legend(fontsize=8)
    axes[2].grid(True)

    plt.tight_layout()
    return fig


# ── Wavelet computation ───────────────────────────────────────────────────────

def compute_wavelet(signal1, signal2, dt):
    """Compute WPS1, WPS2, WCS, WCO from two time-domain signals.

    Identical estimator to wavelet_v8.compute_wavelet (biased auto/cross-
    correlograms, bilateral construction, cmor CWT), with sampling_period = dt
    (seconds) and a Hz frequency grid.
    """
    N = min(len(signal1), len(signal2))
    s1, s2 = signal1[:N], signal2[:N]

    def biased_autocorr(x):
        return np.array([np.sum(x[:N - m] * x[m:]) / N for m in range(N)])

    def biased_crosscorr(x, y):
        return np.array([np.sum(x[:N - m] * y[m:]) / N for m in range(N)])

    Rxx = biased_autocorr(s1)
    Ryy = biased_autocorr(s2)
    Rxy = biased_crosscorr(s1, s2)
    Ryx = biased_crosscorr(s2, s1)

    Rxx_full = np.concatenate([Rxx[N - 1:0:-1], Rxx])
    Ryy_full = np.concatenate([Ryy[N - 1:0:-1], Ryy])
    Rxy_full = np.concatenate([Ryx[N - 1:0:-1], Rxy])
    lags     = np.arange(-(N - 1), N)
    assert Rxx_full.shape[0] == 2 * N - 1  # N-1 is zero-lag column

    center_freq = pywt.central_frequency(WAVELET)
    freqs       = np.logspace(np.log10(FREQ_MIN), np.log10(FREQ_MAX), N_SCALES)
    # scale = center_freq / (freq * dt)  == center_freq / (freq * sampling_period)
    scales      = center_freq / (freqs * dt)

    WPS1, _ = pywt.cwt(Rxx_full.real, scales, WAVELET, sampling_period=dt)
    WPS2, _ = pywt.cwt(Ryy_full.real, scales, WAVELET, sampling_period=dt)
    WCS,  _ = pywt.cwt(Rxy_full.real, scales, WAVELET, sampling_period=dt)
    WCO     = np.abs(WCS) ** 2 / (np.abs(WPS1[:, [N - 1]]) * np.abs(WPS2[:, [N - 1]]))

    lags_s = lags * dt
    return WPS1, WPS2, WCS, WCO, lags_s, freqs, N


# ── Wavelet figure ────────────────────────────────────────────────────────────

def plot_wavelet(WPS1, WPS2, WCS, WCO, lags_s, freqs, N,
                 name1, flow1, name2, flow2, coi1, coi2):
    """3x4 wavelet figure with Hz frequency axis and seconds lag axis."""
    mean_power1 = np.mean(np.abs(WPS1), axis=1)
    mean_power2 = np.mean(np.abs(WPS2), axis=1)
    power_xlim  = max(mean_power1.max(), mean_power2.max()) * 1.1

    LL, FF = np.meshgrid(lags_s, freqs)

    fig = plt.figure(figsize=(20, 13))
    gs  = GridSpec(3, 4, figure=fig, width_ratios=[9, 2.5, 9, 2.5],
                   hspace=0.35, wspace=0.35)

    ax_wps1 = fig.add_subplot(gs[0, 0])
    im0 = ax_wps1.pcolormesh(LL, FF, np.abs(WPS1), cmap="viridis", shading="auto")
    fig.colorbar(im0, ax=ax_wps1, label="Magnitude", fraction=0.046, pad=0.04)
    ax_wps1.set_title(r"WPS — %s (%s) from $R_{xx}$" % (name1, flow1))
    ax_wps1.set_xlabel("Lag (s)")
    ax_wps1.set_ylabel("Frequency (Hz)")
    apply_freq_axis(ax_wps1, coi1)

    ax_mp1 = fig.add_subplot(gs[0, 1], sharey=ax_wps1)
    ax_mp1.plot(mean_power1, freqs, color="tab:blue")
    ax_mp1.set_xlim(0, power_xlim)
    ax_mp1.xaxis.set_major_locator(MaxNLocator(3))
    ax_mp1.xaxis.set_major_formatter(EngFormatter())
    ax_mp1.set_title("Mean power")
    ax_mp1.set_xlabel("Mean |WPS|")
    apply_freq_axis(ax_mp1, coi1)
    plt.setp(ax_mp1.get_yticklabels(), visible=False)

    ax_wps2 = fig.add_subplot(gs[0, 2])
    im1 = ax_wps2.pcolormesh(LL, FF, np.abs(WPS2), cmap="viridis", shading="auto")
    fig.colorbar(im1, ax=ax_wps2, label="Magnitude", fraction=0.046, pad=0.04)
    ax_wps2.set_title(r"WPS — %s (%s) from $R_{yy}$" % (name2, flow2))
    ax_wps2.set_xlabel("Lag (s)")
    ax_wps2.set_ylabel("Frequency (Hz)")
    apply_freq_axis(ax_wps2, coi2)

    ax_mp2 = fig.add_subplot(gs[0, 3], sharey=ax_wps2)
    ax_mp2.plot(mean_power2, freqs, color="tab:blue")
    ax_mp2.set_xlim(0, power_xlim)
    ax_mp2.xaxis.set_major_locator(MaxNLocator(3))
    ax_mp2.xaxis.set_major_formatter(EngFormatter())
    ax_mp2.set_title("Mean power")
    ax_mp2.set_xlabel("Mean |WPS|")
    apply_freq_axis(ax_mp2, coi2)
    plt.setp(ax_mp2.get_yticklabels(), visible=False)

    ax_wcs = fig.add_subplot(gs[1, 0:4])
    im2 = ax_wcs.pcolormesh(LL, FF, np.abs(WCS), cmap="coolwarm", shading="auto")
    fig.colorbar(im2, ax=ax_wcs, label="Cross Power", fraction=0.023, pad=0.02)
    ax_wcs.set_title(
        "Wavelet Cross Spectrum — %s (%s) x %s (%s)  "
        "[lag < 0: %s leads | lag > 0: %s leads]"
        % (name1, flow1, name2, flow2, name2, name1))
    ax_wcs.axvline(0, color="white", linewidth=0.8, linestyle="--")
    ax_wcs.set_xlabel("Lag (s)")
    ax_wcs.set_ylabel("Frequency (Hz)")
    apply_freq_axis(ax_wcs, max(coi1, coi2))

    ax_wco = fig.add_subplot(gs[2, 0:4])
    im3 = ax_wco.pcolormesh(LL, FF, WCO, cmap="plasma", shading="auto", vmin=0, vmax=1)
    fig.colorbar(im3, ax=ax_wco, label="Coherence", fraction=0.023, pad=0.02)
    ax_wco.set_title(
        "Wavelet Coherence — %s (%s) x %s (%s)  "
        "[lag < 0: %s leads | lag > 0: %s leads]"
        % (name1, flow1, name2, flow2, name2, name1))
    ax_wco.axvline(0, color="white", linewidth=0.8, linestyle="--")
    ax_wco.set_xlabel("Lag (s)")
    ax_wco.set_ylabel("Frequency (Hz)")
    apply_freq_axis(ax_wco, max(coi1, coi2))

    return fig


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    SHOW  = "--show" in sys.argv
    _pos  = [a for a in sys.argv[1:] if not a.startswith("--")]
    NAME1 = _pos[0] if len(_pos) > 0 else "m05q10_001"
    NAME2 = _pos[1] if len(_pos) > 1 else "m05q20_001"
    FLOW1 = folder_flow_type(NAME1)
    FLOW2 = folder_flow_type(NAME2)

    sig1, meta1 = load_signal(NAME1)
    sig2, meta2 = load_signal(NAME2)

    # Both runs of a same-pattern pair share fps in this dataset; guard anyway.
    if abs(meta1["dt_s"] - meta2["dt_s"]) > 1e-12:
        print("WARNING: differing dt (%g vs %g s); using signal-1 dt for the axis."
              % (meta1["dt_s"], meta2["dt_s"]))
    dt = meta1["dt_s"]

    fig_sig = plot_signals(sig1, sig2, dt, NAME1, FLOW1, NAME2, FLOW2)

    WPS1, WPS2, WCS, WCO, lags_s, freqs, N = compute_wavelet(sig1, sig2, dt)
    coi1 = coi_freq_floor(N, dt, meta1["f_min_run_hz"])
    coi2 = coi_freq_floor(N, dt, meta2["f_min_run_hz"])
    fig_wav = plot_wavelet(WPS1, WPS2, WCS, WCO, lags_s, freqs, N,
                           NAME1, FLOW1, NAME2, FLOW2, coi1, coi2)

    if SHOW:
        plt.show()
    else:
        out_dir = os.path.join(HERE, "%s_%s__x__%s_%s" % (NAME1, FLOW1, NAME2, FLOW2))
        os.makedirs(out_dir, exist_ok=True)
        fig_sig.savefig(os.path.join(out_dir, "signals.png"), dpi=120, bbox_inches="tight")
        fig_wav.savefig(os.path.join(out_dir, "wavelet.png"), dpi=120, bbox_inches="tight")
        plt.close("all")
        print("saved: %s" % out_dir)
