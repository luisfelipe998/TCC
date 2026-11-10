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

# ── Signal source ─────────────────────────────────────────────────────────────
HERE      = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.join(HERE, "../", "data_processed")
OS_COL    = "Voltage_2"

# ── Spatial frequency / wavelength axis ──────────────────────────────────────
FREQ_MIN         = 1 / 3000
FREQ_MAX         = 1 / 50
FREQ_TICKS       = [1/3000, 1/1000, 1/500, 1/200, 1/100, 1/50]
FREQ_TICK_LABELS = ["3000 px", "1000 px", "500 px", "200 px", "100 px", "50 px"]

# ── COI parameters ────────────────────────────────────────────────────────────
COI_EFOLD_FACTOR = np.sqrt(2)
COI_CAP_BY_PATTERN = {
    "wavy annular":      3000,
    "churn":             3000,
    "plug":              3000,
    "slug":              3000,
    "plug-slug":         3000,
    "slug-churn":        3000,
    "dispersed bubbles": 2400,
}
COI_CAP_DEFAULT = 3000
COI_FILL_COLOR  = "0.75"
COI_FILL_ALPHA  = 0.45
COI_LINE_COLOR  = "white"
COI_LINE_WIDTH  = 1.0
COI_LINE_STYLE  = "--"


# ── Helpers ───────────────────────────────────────────────────────────────────

def _load_and_preprocess(path):
    df  = pd.read_excel(path, engine="openpyxl", sheet_name="Sheet")
    if OS_COL not in df.columns:
        raise KeyError("Canal '%s' ausente em %s. Colunas: %s" % (OS_COL, path, list(df.columns)))
    sig = pd.to_numeric(df[OS_COL], errors="coerce").to_numpy(float)
    sig = sig[~np.isnan(sig)]
    sig = scipy_detrend(sig, type="linear")
    sig = sig - sig.mean()
    sig = sig / sig.std(ddof=1)
    return sig


def load_signal(name):
    """Load from __dados_cropped_v11.xlsx (rows 20–59 mean)."""
    path = os.path.join(DATA_ROOT, name, "__dados_cropped_v11.xlsx")
    return _load_and_preprocess(path), 1.0


def load_signal_center(name):
    """Load from __dados_center_v13.xlsx (single center row 40)."""
    path = os.path.join(DATA_ROOT, name, "__dados_center_v13.xlsx")
    return _load_and_preprocess(path), 1.0


def folder_flow_type(name):
    folder = os.path.join(DATA_ROOT, name)
    for p in sorted(glob.glob(os.path.join(folder, "*.txt"))):
        stem = os.path.splitext(os.path.basename(p))[0].lstrip("_")
        if stem != "dados":
            return stem
    return name


def coi_boundary(N, pattern):
    wavelet_center_freq = pywt.central_frequency("cmor2.0-0.9549")
    lam_from_n  = min((2 * N - 1) / (2 * COI_EFOLD_FACTOR * wavelet_center_freq), 3000)
    pattern_cap = COI_CAP_BY_PATTERN.get(pattern, COI_CAP_DEFAULT)
    return min(lam_from_n, pattern_cap)


def apply_freq_axis(ax, coi_wavelengths=()):
    ax.set_yscale("log")
    ax.set_ylim(FREQ_MIN, FREQ_MAX)
    ax.yaxis.set_major_locator(FixedLocator(FREQ_TICKS))
    ax.yaxis.set_major_formatter(FixedFormatter(FREQ_TICK_LABELS))
    ax.yaxis.set_minor_locator(FixedLocator([]))

    extra_ticks, extra_labels = [], []
    for lam_coi in coi_wavelengths:
        freq_coi = 1.0 / lam_coi
        if FREQ_MIN < freq_coi < FREQ_MAX:
            ax.axhspan(FREQ_MIN, freq_coi, color=COI_FILL_COLOR, alpha=COI_FILL_ALPHA, zorder=5)
            ax.axhline(freq_coi, color=COI_LINE_COLOR, linewidth=COI_LINE_WIDTH,
                       linestyle=COI_LINE_STYLE, zorder=6)
            extra_ticks.append(freq_coi)
            extra_labels.append(f"{lam_coi:.0f} px")

    if extra_ticks:
        ax.yaxis.set_major_locator(FixedLocator(FREQ_TICKS + extra_ticks))
        ax.yaxis.set_major_formatter(FixedFormatter(FREQ_TICK_LABELS + extra_labels))


# ── Signals figure ───────────────────────────────────────────────────────────

def plot_signals(signal1, signal2, name1, flow1, name2, flow2, label):
    """Build and return the 3-panel signals figure."""
    N = min(len(signal1), len(signal2))
    s1, s2 = signal1[:N], signal2[:N]
    n = np.arange(N)

    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    fig.suptitle(f"{name1} ({flow1})  vs  {name2} ({flow2})  [{label}]", fontsize=11)

    axes[0].plot(n, s1, color="tab:blue", linewidth=0.6)
    axes[0].set_title(f"Signal 1 — {name1} ({flow1})")
    axes[0].set_ylabel("Amplitude (σ)")
    axes[0].grid(True)

    axes[1].plot(n, s2, color="tab:orange", linewidth=0.6)
    axes[1].set_title(f"Signal 2 — {name2} ({flow2})")
    axes[1].set_ylabel("Amplitude (σ)")
    axes[1].grid(True)

    axes[2].plot(n, s1, color="tab:blue",   linewidth=0.6, label=f"{name1} ({flow1})")
    axes[2].plot(n, s2, color="tab:orange", linewidth=0.6, alpha=0.7, label=f"{name2} ({flow2})")
    axes[2].set_title(f"Both signals (z-scored) — {label}")
    axes[2].set_ylabel("Amplitude (σ)")
    axes[2].set_xlabel("Position (px)")
    axes[2].legend(fontsize=8)
    axes[2].grid(True)

    plt.tight_layout()
    return fig


# ── Wavelet computation ───────────────────────────────────────────────────────

def compute_wavelet(signal1, signal2):
    """Compute WPS1, WPS2, WCS, WCO from two spatial signals."""
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
    assert Rxx_full.shape[0] == 2 * N - 1  # N-1 is zero-lag column; fails fast if bilateral construction changes

    wavelet     = "cmor2.0-0.9549"
    center_freq = pywt.central_frequency(wavelet)
    freqs       = np.logspace(np.log10(FREQ_MIN), np.log10(FREQ_MAX), 200)
    scales      = center_freq / freqs   # sampling_rate = 1.0

    WPS1, _ = pywt.cwt(Rxx_full.real, scales, wavelet, sampling_period=1.0)
    WPS2, _ = pywt.cwt(Ryy_full.real, scales, wavelet, sampling_period=1.0)
    WCS,  _ = pywt.cwt(Rxy_full.real, scales, wavelet, sampling_period=1.0)
    WCO     = np.abs(WCS) ** 2 / (np.abs(WPS1[:, [N - 1]]) * np.abs(WPS2[:, [N - 1]]))

    return WPS1, WPS2, WCS, WCO, lags, freqs, N


# ── Wavelet figure ────────────────────────────────────────────────────────────

def plot_wavelet(WPS1, WPS2, WCS, WCO, lags, freqs, N,
                 name1, flow1, name2, flow2, coi1, coi2):
    """Build and return the 3×4 wavelet figure."""
    mean_power1 = np.mean(np.abs(WPS1), axis=1)
    mean_power2 = np.mean(np.abs(WPS2), axis=1)
    power_xlim  = max(mean_power1.max(), mean_power2.max()) * 1.1

    LL, FF = np.meshgrid(lags, freqs)

    fig = plt.figure(figsize=(20, 13))
    gs  = GridSpec(3, 4, figure=fig, width_ratios=[9, 2.5, 9, 2.5],
                   hspace=0.35, wspace=0.35)

    ax_wps1 = fig.add_subplot(gs[0, 0])
    im0 = ax_wps1.pcolormesh(LL, FF, np.abs(WPS1), cmap="viridis", shading="auto")
    fig.colorbar(im0, ax=ax_wps1, label="Magnitude", fraction=0.046, pad=0.04)
    ax_wps1.set_title(f"WPS — {name1} ({flow1}) from $R_{{xx}}$")
    ax_wps1.set_xlabel("Lag (px)")
    ax_wps1.set_ylabel("Wavelength (px)")
    apply_freq_axis(ax_wps1, coi_wavelengths=[coi1])

    ax_mp1 = fig.add_subplot(gs[0, 1], sharey=ax_wps1)
    ax_mp1.plot(mean_power1, freqs, color="tab:blue")
    ax_mp1.set_xlim(0, power_xlim)
    ax_mp1.xaxis.set_major_locator(MaxNLocator(3))
    ax_mp1.xaxis.set_major_formatter(EngFormatter())
    ax_mp1.set_title("Mean power")
    ax_mp1.set_xlabel("Mean |WPS|")
    apply_freq_axis(ax_mp1, coi_wavelengths=[coi1])
    plt.setp(ax_mp1.get_yticklabels(), visible=False)

    ax_wps2 = fig.add_subplot(gs[0, 2])
    im1 = ax_wps2.pcolormesh(LL, FF, np.abs(WPS2), cmap="viridis", shading="auto")
    fig.colorbar(im1, ax=ax_wps2, label="Magnitude", fraction=0.046, pad=0.04)
    ax_wps2.set_title(f"WPS — {name2} ({flow2}) from $R_{{yy}}$")
    ax_wps2.set_xlabel("Lag (px)")
    ax_wps2.set_ylabel("Wavelength (px)")
    apply_freq_axis(ax_wps2, coi_wavelengths=[coi2])

    ax_mp2 = fig.add_subplot(gs[0, 3], sharey=ax_wps2)
    ax_mp2.plot(mean_power2, freqs, color="tab:blue")
    ax_mp2.set_xlim(0, power_xlim)
    ax_mp2.xaxis.set_major_locator(MaxNLocator(3))
    ax_mp2.xaxis.set_major_formatter(EngFormatter())
    ax_mp2.set_title("Mean power")
    ax_mp2.set_xlabel("Mean |WPS|")
    apply_freq_axis(ax_mp2, coi_wavelengths=[coi2])
    plt.setp(ax_mp2.get_yticklabels(), visible=False)

    ax_wcs = fig.add_subplot(gs[1, 0:4])
    im2 = ax_wcs.pcolormesh(LL, FF, np.abs(WCS), cmap="coolwarm", shading="auto")
    fig.colorbar(im2, ax=ax_wcs, label="Cross Power", fraction=0.023, pad=0.02)
    ax_wcs.set_title(
        f"Wavelet Cross Spectrum — {name1} ({flow1}) x {name2} ({flow2})  "
        f"[lag < 0: {name2} leads | lag > 0: {name1} leads]"
    )
    ax_wcs.axvline(0, color="white", linewidth=0.8, linestyle="--")
    ax_wcs.set_xlabel("Lag (px)")
    ax_wcs.set_ylabel("Wavelength (px)")
    apply_freq_axis(ax_wcs, coi_wavelengths=[coi1, coi2])

    ax_wco = fig.add_subplot(gs[2, 0:4])
    im3 = ax_wco.pcolormesh(LL, FF, WCO, cmap="plasma", shading="auto", vmin=0, vmax=1)
    fig.colorbar(im3, ax=ax_wco, label="Coherence", fraction=0.023, pad=0.02)
    ax_wco.set_title(
        f"Wavelet Coherence — {name1} ({flow1}) x {name2} ({flow2})  "
        f"[lag < 0: {name2} leads | lag > 0: {name1} leads]"
    )
    ax_wco.axvline(0, color="white", linewidth=0.8, linestyle="--")
    ax_wco.set_xlabel("Lag (px)")
    ax_wco.set_ylabel("Wavelength (px)")
    apply_freq_axis(ax_wco, coi_wavelengths=[coi1, coi2])

    return fig


# ── Entry point ───────────────────────────────────────────────────────────────

SHOW  = "--show" in sys.argv
_pos  = [a for a in sys.argv[1:] if not a.startswith("--")]
NAME1 = _pos[0] if len(_pos) > 0 else "synthetic1"
NAME2 = _pos[1] if len(_pos) > 1 else "synthetic2"
FLOW1 = folder_flow_type(NAME1)
FLOW2 = folder_flow_type(NAME2)

# ── Load cropped signals ──────────────────────────────────────────────────────
sig1_cropped, _ = load_signal(NAME1)
sig2_cropped, _ = load_signal(NAME2)

# ── Load center signals ───────────────────────────────────────────────────────
sig1_center, _ = load_signal_center(NAME1)
sig2_center, _ = load_signal_center(NAME2)

# ── Signals figures ───────────────────────────────────────────────────────────
fig_sig_cropped = plot_signals(sig1_cropped, sig2_cropped, NAME1, FLOW1, NAME2, FLOW2, "cropped rows 20–59")
fig_sig_center  = plot_signals(sig1_center,  sig2_center,  NAME1, FLOW1, NAME2, FLOW2, "center row 40")

# ── Wavelet: cropped ──────────────────────────────────────────────────────────
WPS1_c, WPS2_c, WCS_c, WCO_c, lags_c, freqs_c, N_c = compute_wavelet(sig1_cropped, sig2_cropped)
coi1_c = coi_boundary(N_c, FLOW1)
coi2_c = coi_boundary(N_c, FLOW2)
fig_wav_cropped = plot_wavelet(WPS1_c, WPS2_c, WCS_c, WCO_c, lags_c, freqs_c, N_c,
                               NAME1, FLOW1, NAME2, FLOW2, coi1_c, coi2_c)

# ── Wavelet: center ───────────────────────────────────────────────────────────
WPS1_k, WPS2_k, WCS_k, WCO_k, lags_k, freqs_k, N_k = compute_wavelet(sig1_center, sig2_center)
coi1_k = coi_boundary(N_k, FLOW1)
coi2_k = coi_boundary(N_k, FLOW2)
fig_wav_center = plot_wavelet(WPS1_k, WPS2_k, WCS_k, WCO_k, lags_k, freqs_k, N_k,
                              NAME1, FLOW1, NAME2, FLOW2, coi1_k, coi2_k)

# ── Save or show ──────────────────────────────────────────────────────────────
if SHOW:
    plt.show()
else:
    out_dir = os.path.join(HERE, "%s_%s__x__%s_%s" % (NAME1, FLOW1, NAME2, FLOW2))
    os.makedirs(out_dir, exist_ok=True)
    stale = os.path.join(out_dir, "signals.png")
    if os.path.exists(stale):
        os.remove(stale)
    fig_sig_cropped.savefig(os.path.join(out_dir, "signals_cropped.png"), dpi=120, bbox_inches="tight")
    fig_sig_center.savefig( os.path.join(out_dir, "signals_center.png"),  dpi=120, bbox_inches="tight")
    fig_wav_cropped.savefig(os.path.join(out_dir, "wavelet_cropped.png"), dpi=120, bbox_inches="tight")
    fig_wav_center.savefig( os.path.join(out_dir, "wavelet_center.png"),  dpi=120, bbox_inches="tight")
    plt.close("all")
    print("saved: %s" % out_dir)
