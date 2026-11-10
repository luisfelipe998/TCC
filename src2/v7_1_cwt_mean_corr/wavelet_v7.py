"""v7_1_cwt_mean_corr — LOO multi-realization WCO (CWT each realization first, then mean magnitudes).
UNBIASED correlation estimator: divides by (N-m) instead of N.

Usage:
    python wavelet_v7.py <probe_run> <ref1> <ref2> ... [--show]

For probe run x and reference runs y1..yK-1 (all same pattern):
  1. Truncate all signals to min(N)
  2. Compute unbiased R_xx (probe), R_yy^j and R_xy^j for each ref j
  3. CWT each R_yy^j and R_xy^j individually
  4. WPS_x = |CWT(R_xx)| at lag=0
     WPS_y = mean(|CWT(R_yy^j)|) at lag=0 over j
     WCS   = mean(|CWT(R_xy^j)|) over j
  5. WCO = WCS^2 / (WPS_x * WPS_y)
"""
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

HERE      = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.join(HERE, "../", "data_processed")
OS_COL    = "Voltage_2"

FREQ_MIN         = 1 / 3000
FREQ_MAX         = 1 / 50
FREQ_TICKS       = [1/3000, 1/1000, 1/500, 1/200, 1/100, 1/50]
FREQ_TICK_LABELS = ["3000 px", "1000 px", "500 px", "200 px", "100 px", "50 px"]

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

COI_FILL_COLOR = "0.75"
COI_FILL_ALPHA = 0.45
COI_LINE_COLOR = "white"
COI_LINE_WIDTH = 1.0
COI_LINE_STYLE = "--"


def load_signal(name):
    path = os.path.join(DATA_ROOT, name, "__dados_cropped_v11.xlsx")
    df   = pd.read_excel(path, engine="openpyxl", sheet_name="Sheet")
    if OS_COL not in df.columns:
        raise KeyError("Canal '%s' ausente em %s. Colunas: %s" % (OS_COL, path, list(df.columns)))
    sig = pd.to_numeric(df[OS_COL], errors="coerce").to_numpy(float)
    sig = sig[~np.isnan(sig)]
    sig = 255.0 - sig
    sig = scipy_detrend(sig, type="linear")
    sig = sig - sig.mean()
    sig = sig / sig.std()
    return sig


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

    extra_ticks  = []
    extra_labels = []
    for lam_coi in coi_wavelengths:
        freq_coi = 1.0 / lam_coi
        if FREQ_MIN < freq_coi < FREQ_MAX:
            ax.axhspan(FREQ_MIN, freq_coi,
                       color=COI_FILL_COLOR, alpha=COI_FILL_ALPHA, zorder=5)
            ax.axhline(freq_coi,
                       color=COI_LINE_COLOR, linewidth=COI_LINE_WIDTH,
                       linestyle=COI_LINE_STYLE, zorder=6)
            extra_ticks.append(freq_coi)
            extra_labels.append(f"{lam_coi:.0f} px")

    if extra_ticks:
        combined_ticks  = FREQ_TICKS  + extra_ticks
        combined_labels = FREQ_TICK_LABELS + extra_labels
        ax.yaxis.set_major_locator(FixedLocator(combined_ticks))
        ax.yaxis.set_major_formatter(FixedFormatter(combined_labels))


# ── Entry point ───────────────────────────────────────────────────────────────

SHOW  = "--show" in sys.argv
_pos  = [a for a in sys.argv[1:] if not a.startswith("--")]
if len(_pos) < 2:
    print("Usage: wavelet_v7.py <probe_run> <ref1> [ref2 ...] [--show]")
    sys.exit(1)

PROBE_NAME = _pos[0]
REF_NAMES  = _pos[1:]
K_REFS     = len(REF_NAMES)
FLOW       = folder_flow_type(PROBE_NAME)

# ── Load all signals ──────────────────────────────────────────────────────────
probe_raw = load_signal(PROBE_NAME)
refs_raw  = [load_signal(r) for r in REF_NAMES]

N = min(len(probe_raw), min(len(r) for r in refs_raw))
probe_sig = probe_raw[:N]
ref_sigs  = [r[:N] for r in refs_raw]

sampling_rate = 1.0

# ── Biased correlations (divide by N) ────────────────────────────────────────
def biased_autocorr(x):
    return np.array([np.sum(x[:N - m] * x[m:]) / (N - m) for m in range(N)])

def biased_crosscorr(x, y):
    return np.array([np.sum(x[:N - m] * y[m:]) / (N - m) for m in range(N)])

# Probe autocorrelation (WPS_x)
Rxx = biased_autocorr(probe_sig)
Rxx_full = np.concatenate([Rxx[N - 1:0:-1], Rxx])

lags = np.arange(-(N - 1), N)

# ── CWT setup ────────────────────────────────────────────────────────────────
wavelet     = "cmor2.0-0.9549"
center_freq = pywt.central_frequency(wavelet)
freqs       = np.logspace(np.log10(FREQ_MIN), np.log10(FREQ_MAX), 200)
scales      = center_freq * sampling_rate / freqs
n_scales    = len(scales)
n_lags      = 2 * N - 1

# WPS_x from probe
CWT_xx, _ = pywt.cwt(Rxx_full.real, scales, wavelet, sampling_period=1 / sampling_rate)
WPS_x = np.abs(CWT_xx[:, [N - 1]])   # (n_scales, 1)

# Per-reference CWTs, accumulate for mean
WPS_y_accum = np.zeros((n_scales, 1))
WCS_accum   = np.zeros((n_scales, n_lags))
CWT_yy_mean = np.zeros((n_scales, n_lags))   # for plotting WPS_y panel

for ref in ref_sigs:
    ryy = biased_autocorr(ref)
    ryy_full = np.concatenate([ryy[N - 1:0:-1], ryy])

    rxy     = biased_crosscorr(probe_sig, ref)
    ryx     = biased_crosscorr(ref, probe_sig)
    rxy_full = np.concatenate([ryx[N - 1:0:-1], rxy])

    cwt_yy, _ = pywt.cwt(ryy_full.real, scales, wavelet, sampling_period=1 / sampling_rate)
    cwt_xy, _ = pywt.cwt(rxy_full.real, scales, wavelet, sampling_period=1 / sampling_rate)

    WPS_y_accum += np.abs(cwt_yy[:, [N - 1]])
    WCS_accum   += np.abs(cwt_xy)
    CWT_yy_mean += np.abs(cwt_yy)

WPS_y     = WPS_y_accum / K_REFS   # (n_scales, 1)
WCS       = WCS_accum   / K_REFS   # (n_scales, n_lags)
CWT_yy_mean /= K_REFS

WCO = WCS ** 2 / (WPS_x * WPS_y)

# ── COI boundary ──────────────────────────────────────────────────────────────
coi_lam = coi_boundary(N, FLOW)

# ── Mean power profiles ───────────────────────────────────────────────────────
mean_power_x = np.mean(np.abs(CWT_xx), axis=1)
mean_power_y = np.mean(CWT_yy_mean, axis=1)
POWER_MARGIN = 1.1
power_xlim   = max(mean_power_x.max(), mean_power_y.max()) * POWER_MARGIN

# ── Meshgrid ──────────────────────────────────────────────────────────────────
LL, FF = np.meshgrid(lags, freqs)

# ── Output path ───────────────────────────────────────────────────────────────
out_name = "%s_mean_K%d_%s__x__%s_%s" % (FLOW, K_REFS, FLOW, PROBE_NAME, FLOW)
out_dir  = os.path.join(HERE, out_name)
os.makedirs(out_dir, exist_ok=True)

# ── Figure 1: signals ─────────────────────────────────────────────────────────
fig_sig, axes_sig = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
fig_sig.suptitle(
    f"{PROBE_NAME} ({FLOW})  vs  mean of {K_REFS} other {FLOW} runs\n"
    f"[v7_cwt_mean_corr: CWT each ref then mean]",
    fontsize=10
)

axes_sig[0].plot(np.arange(N), probe_sig, color="tab:blue", linewidth=0.6)
axes_sig[0].set_title(f"Probe — {PROBE_NAME} ({FLOW})")
axes_sig[0].set_ylabel("Amplitude (σ)")
axes_sig[0].grid(True)

axes_sig[1].plot(np.arange(N), probe_sig, color="tab:blue", linewidth=0.6, label=PROBE_NAME)
for ref_name, ref in zip(REF_NAMES, ref_sigs):
    axes_sig[1].plot(np.arange(N), ref, linewidth=0.4, alpha=0.5, label=ref_name)
axes_sig[1].set_title(f"Probe + {K_REFS} reference runs ({FLOW})")
axes_sig[1].set_ylabel("Amplitude (σ)")
axes_sig[1].set_xlabel("Position (px)")
axes_sig[1].legend(fontsize=7, ncol=min(K_REFS + 1, 5))
axes_sig[1].grid(True)

plt.tight_layout()

# ── Figure 2: wavelet panels ──────────────────────────────────────────────────
fig = plt.figure(figsize=(20, 13))
gs  = GridSpec(3, 4, figure=fig, width_ratios=[9, 2.5, 9, 2.5],
               hspace=0.35, wspace=0.35)

# WPS_x (probe)
ax_wps1 = fig.add_subplot(gs[0, 0])
im0 = ax_wps1.pcolormesh(LL, FF, np.abs(CWT_xx), cmap="viridis", shading="auto")
fig.colorbar(im0, ax=ax_wps1, label="Magnitude", fraction=0.046, pad=0.04)
ax_wps1.set_title(f"WPS — probe {PROBE_NAME} ({FLOW})")
ax_wps1.set_xlabel("Lag (px)")
ax_wps1.set_ylabel("Wavelength (px)")
apply_freq_axis(ax_wps1, coi_wavelengths=[coi_lam])

ax_mp1 = fig.add_subplot(gs[0, 1], sharey=ax_wps1)
ax_mp1.plot(mean_power_x, freqs, color="tab:blue")
ax_mp1.set_xlim(0, power_xlim)
ax_mp1.xaxis.set_major_locator(MaxNLocator(3))
ax_mp1.xaxis.set_major_formatter(EngFormatter())
ax_mp1.set_title("Mean power")
ax_mp1.set_xlabel("Mean |WPS|")
apply_freq_axis(ax_mp1, coi_wavelengths=[coi_lam])
plt.setp(ax_mp1.get_yticklabels(), visible=False)

# WPS_y (mean of ref CWTs)
ax_wps2 = fig.add_subplot(gs[0, 2])
im1 = ax_wps2.pcolormesh(LL, FF, CWT_yy_mean, cmap="viridis", shading="auto")
fig.colorbar(im1, ax=ax_wps2, label="Magnitude", fraction=0.046, pad=0.04)
ax_wps2.set_title(f"WPS — mean of {K_REFS} refs ({FLOW}) [CWT each then mean]")
ax_wps2.set_xlabel("Lag (px)")
ax_wps2.set_ylabel("Wavelength (px)")
apply_freq_axis(ax_wps2, coi_wavelengths=[coi_lam])

ax_mp2 = fig.add_subplot(gs[0, 3], sharey=ax_wps2)
ax_mp2.plot(mean_power_y, freqs, color="tab:blue")
ax_mp2.set_xlim(0, power_xlim)
ax_mp2.xaxis.set_major_locator(MaxNLocator(3))
ax_mp2.xaxis.set_major_formatter(EngFormatter())
ax_mp2.set_title("Mean power")
ax_mp2.set_xlabel("Mean |WPS|")
apply_freq_axis(ax_mp2, coi_wavelengths=[coi_lam])
plt.setp(ax_mp2.get_yticklabels(), visible=False)

# WCS
ax_wcs = fig.add_subplot(gs[1, 0:4])
im2 = ax_wcs.pcolormesh(LL, FF, WCS, cmap="coolwarm", shading="auto")
fig.colorbar(im2, ax=ax_wcs, label="Cross Power", fraction=0.023, pad=0.02)
ax_wcs.set_title(
    f"WCS — {PROBE_NAME} ({FLOW}) x mean of {K_REFS} refs  "
    f"[lag<0: refs lead | lag>0: probe leads]"
)
ax_wcs.axvline(0, color="white", linewidth=0.8, linestyle="--")
ax_wcs.set_xlabel("Lag (px)")
ax_wcs.set_ylabel("Wavelength (px)")
apply_freq_axis(ax_wcs, coi_wavelengths=[coi_lam])

# WCO
ax_wco = fig.add_subplot(gs[2, 0:4])
im3 = ax_wco.pcolormesh(LL, FF, WCO, cmap="plasma", shading="auto", vmin=0, vmax=1)
fig.colorbar(im3, ax=ax_wco, label="Coherence", fraction=0.023, pad=0.02)
ax_wco.set_title(
    f"WCO — {PROBE_NAME} ({FLOW}) x mean of {K_REFS} refs  "
    f"[v7_cwt_mean_corr | lag<0: refs lead | lag>0: probe leads]"
)
ax_wco.axvline(0, color="white", linewidth=0.8, linestyle="--")
ax_wco.set_xlabel("Lag (px)")
ax_wco.set_ylabel("Wavelength (px)")
apply_freq_axis(ax_wco, coi_wavelengths=[coi_lam])

# ── Save or show ──────────────────────────────────────────────────────────────
if SHOW:
    plt.show()
else:
    fig_sig.savefig(os.path.join(out_dir, "signals.png"), dpi=120, bbox_inches="tight")
    fig.savefig(os.path.join(out_dir, "wavelet.png"),     dpi=120, bbox_inches="tight")
    plt.close("all")
    print("saved: %s" % out_dir)
