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
# Full axis always shown 50–3000 px; COI overlay marks the unreliable region.
FREQ_MIN         = 1 / 3000   # cycles/px  ↔  3000 px wavelength
FREQ_MAX         = 1 / 50     # cycles/px  ↔  50 px wavelength
FREQ_TICKS       = [1/3000, 1/1000, 1/500, 1/200, 1/100, 1/50]
FREQ_TICK_LABELS = ["3000 px", "1000 px", "500 px", "200 px", "100 px", "50 px"]

# ── COI (Cone of Influence) parameters ───────────────────────────────────────
# The COI boundary is the MORE RESTRICTIVE of two limits:
#   1. λ_max_from_N  — Morlet e-folding limit derived from signal length N
#                      λ = center_freq * sqrt(2) * N/2
#   2. COI_CAP_BY_PATTERN — per-pattern upper bound from power-distribution
#                           analysis in v5_spatial/INSIGHTS.md.
#
# Anything above the boundary (longer wavelengths) is shaded as unreliable.

# Torrence & Compo (1998) eq. 19: for a complex Morlet wavelet the COI
# e-folding distance at scale s is sqrt(2)*s samples from each edge.
# At the center of a signal of length N (farthest point from either edge),
# the maximum reliable scale is s_max = N/(2*sqrt(2)).
# pywt converts scale to wavelength as:  lambda = scale / center_freq
# Therefore:  lambda_max = N / (2 * sqrt(2) * center_freq)
COI_EFOLD_FACTOR = np.sqrt(2)   # Morlet e-folding factor (T&C 1998 eq. 19)

# Per-pattern caps (px) — floor-to-100 of lambda_max for the shortest run in
# each pattern, derived from INSIGHTS.md power-distribution analysis.
# Key = exact string returned by folder_flow_type().
COI_CAP_BY_PATTERN = {
    "wavy annular":      3000,   # all runs lambda_max > 5300 px; full axis reliable
    "churn":             2600,   # shortest run m05q10_001 N=7001 → lambda_max=2640
    "plug":              2200,   # shortest run m05q01_001 N=5940 → lambda_max=2240
    "slug":              1900,   # shortest run m04q01_001 N=5071 → lambda_max=1912
    "plug-slug":         1900,   # only run  m03q01_001 N=5148  → lambda_max=1941
    "slug-churn":        2000,   # shortest run m04q05_001 N=5460 → lambda_max=2059
    "dispersed bubbles": 1100,   # only run  m07q01_001 N=3171  → lambda_max=1196
    # "smooth annular" omitted — degenerate run, no reliable analysis possible
}
COI_CAP_DEFAULT = 3000   # fallback for unknown / synthetic patterns

# Visual style for the unreliable region overlay
COI_FILL_COLOR = "0.75"   # light grey
COI_FILL_ALPHA = 0.45
COI_LINE_COLOR = "white"
COI_LINE_WIDTH = 1.0
COI_LINE_STYLE = "--"


# ── Helpers ───────────────────────────────────────────────────────────────────

def load_signal(name):
    """Load Voltage_2 from __dados_cropped_v11.xlsx and apply the pipeline_v11
    preprocessing chain (invert → detrend → zero-mean → z-score).

    Returns (signal, sampling_rate) where sampling_rate = 1.0 px/sample.
    """
    path = os.path.join(DATA_ROOT, name, "__dados_cropped_v11.xlsx")
    df   = pd.read_excel(path, engine="openpyxl", sheet_name="Sheet")
    if OS_COL not in df.columns:
        raise KeyError("Canal '%s' ausente em %s. Colunas: %s" % (OS_COL, path, list(df.columns)))
    sig = pd.to_numeric(df[OS_COL], errors="coerce").to_numpy(float)
    sig = sig[~np.isnan(sig)]
    sig = 255.0 - sig                        # invert: liquid → high, gas → low
    sig = scipy_detrend(sig, type="linear")  # remove illumination gradient
    sig = sig - sig.mean()                   # zero-mean
    sig = sig / sig.std()                    # z-score
    return sig, 1.0                          # 1 px/sample, spatial domain


def folder_flow_type(name):
    """Return flow pattern label from the __<pattern>.txt marker file."""
    folder = os.path.join(DATA_ROOT, name)
    for p in sorted(glob.glob(os.path.join(folder, "*.txt"))):
        stem = os.path.splitext(os.path.basename(p))[0].lstrip("_")
        if stem != "dados":
            return stem
    return name


def coi_boundary(N, pattern):
    """Effective COI boundary wavelength (px) for a signal of length N and flow pattern.

    Returns the more restrictive of:
      - λ_max_from_N: T&C (1998) eq.19 Morlet e-folding limit
                      = N / (2 * sqrt(2) * center_freq)
      - COI_CAP_BY_PATTERN[pattern]: per-pattern cap from INSIGHTS.md analysis
    """
    wavelet_center_freq = pywt.central_frequency("cmor2.0-0.9549")
    lam_from_n  = N / (2 * COI_EFOLD_FACTOR * wavelet_center_freq)
    pattern_cap = COI_CAP_BY_PATTERN.get(pattern, COI_CAP_DEFAULT)
    return min(lam_from_n, pattern_cap)


def apply_freq_axis(ax, coi_wavelengths=()):
    """Apply log wavelength axis ticks and optional COI shading.

    coi_wavelengths: sequence of boundary wavelengths (px). For each value that
    falls within the visible range, fills the unreliable region above it with a
    semi-transparent grey band and draws a dashed line at the boundary.
    """
    ax.set_yscale("log")
    ax.set_ylim(FREQ_MIN, FREQ_MAX)
    ax.yaxis.set_major_locator(FixedLocator(FREQ_TICKS))
    ax.yaxis.set_major_formatter(FixedFormatter(FREQ_TICK_LABELS))
    ax.yaxis.set_minor_locator(FixedLocator([]))

    extra_ticks  = []
    extra_labels = []
    for lam_coi in coi_wavelengths:
        freq_coi = 1.0 / lam_coi   # convert wavelength → frequency for the y-axis
        if FREQ_MIN < freq_coi < FREQ_MAX:
            # Fill unreliable region: wavelengths LONGER than COI boundary = LOWER frequencies
            ax.axhspan(FREQ_MIN, freq_coi,
                       color=COI_FILL_COLOR, alpha=COI_FILL_ALPHA, zorder=5)
            # Dashed boundary line
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
NAME1 = _pos[0] if len(_pos) > 0 else "synthetic1"
NAME2 = _pos[1] if len(_pos) > 1 else "synthetic2"
FLOW1 = folder_flow_type(NAME1)
FLOW2 = folder_flow_type(NAME2)

# ── Load and truncate signals ─────────────────────────────────────────────────
signal1, fs1 = load_signal(NAME1)
signal2, fs2 = load_signal(NAME2)
if fs1 != fs2:
    print("Warning: fs mismatch %g vs %g — using fs1" % (fs1, fs2))
sampling_rate = fs1   # always 1.0 for spatial-domain signals

N       = min(len(signal1), len(signal2))
signal1 = signal1[:N]
signal2 = signal2[:N]

# ── Correlations (biased estimator, divides by N) ────────────────────────────
# Biased: always divides by N, suppressing noise amplification at large lags.
#   R_xx(m) = (1/N) * sum_{n=0}^{N-1-m} x[n]*x[n+m],  0 <= m <= N-1
def biased_autocorrelation(x):
    return np.array([np.sum(x[:N - m] * x[m:]) / N for m in range(N)])

def biased_crosscorrelation(x, y):
    return np.array([np.sum(x[:N - m] * y[m:]) / N for m in range(N)])

Rxx = biased_autocorrelation(signal1)
Ryy = biased_autocorrelation(signal2)

# Bilateral forms (length 2N-1, lag=0 at index N-1)
# Autocorrelations are symmetric: R_xx(-m) = R_xx(m)
Rxx_full = np.concatenate([Rxx[N - 1:0:-1], Rxx])       # [Rxx[N-1..1] | Rxx[0..N-1]]
Ryy_full = np.concatenate([Ryy[N - 1:0:-1], Ryy])

# Bilateral cross-correlation: negative lags (y leads) | zero | positive lags (x leads)
Rxy      = biased_crosscorrelation(signal1, signal2)
Ryx      = biased_crosscorrelation(signal2, signal1)
lags     = np.arange(-(N - 1), N)                       # length 2N-1
Rxy_full = np.concatenate([Ryx[N - 1:0:-1], Rxy])       # [Ryx[N-1..1] | Rxy[0..N-1]]

# ── CWT (cmor2.0-0.9549, ω₀ ≈ 6.0) ──────────────────────────────────────────
# scale = center_freq * fs / f  →  with fs=1.0: scale = center_freq / f
wavelet     = "cmor2.0-0.9549"
center_freq = pywt.central_frequency(wavelet)
freqs       = np.logspace(np.log10(FREQ_MIN), np.log10(FREQ_MAX), 200)
scales      = center_freq * sampling_rate / freqs

# WPS from bilateral autocorrelations (length 2N-1), WCS from bilateral cross-correlation.
# All three share the same lag axis; lag=0 is at column index N-1 in every array.
# Using bilateral R_xx/R_yy (symmetric) means the CWT at b=0 (center column) sees
# the full symmetric context on both sides — this is the correct WPS_x(0,a) value
# per Schuck & Bodmann eq. 32.14/32.17 (continuous theory, R_xx defined for all tau).
WPS1, _ = pywt.cwt(Rxx_full.real, scales, wavelet, sampling_period=1 / sampling_rate)
WPS2, _ = pywt.cwt(Ryy_full.real, scales, wavelet, sampling_period=1 / sampling_rate)
WCS,  _ = pywt.cwt(Rxy_full.real, scales, wavelet, sampling_period=1 / sampling_rate)

# WCO (Schuck & Bodmann eq. 32.17) — normalised at lag zero (center column N-1)
WCO = np.abs(WCS) ** 2 / (np.abs(WPS1[:, [N - 1]]) * np.abs(WPS2[:, [N - 1]]))

# ── COI boundaries ────────────────────────────────────────────────────────────
coi_boundary1 = coi_boundary(N, FLOW1)   # effective limit for signal 1
coi_boundary2 = coi_boundary(N, FLOW2)   # effective limit for signal 2

# ── Mean power profiles (for side panels) ────────────────────────────────────
mean_power1 = np.mean(np.abs(WPS1), axis=1)
mean_power2 = np.mean(np.abs(WPS2), axis=1)
POWER_MARGIN = 1.1
power_xlim   = max(mean_power1.max(), mean_power2.max()) * POWER_MARGIN

# ── Meshgrids for pcolormesh ──────────────────────────────────────────────────
# WPS1, WPS2, WCS, WCO all have the same bilateral lag axis (length 2N-1).
LL, FF = np.meshgrid(lags, freqs)   # shared grid for all four panels

# ── Figure 1: signals ─────────────────────────────────────────────────────────
n_auto = np.arange(N)
fig_sig = plt.figure(figsize=(10, 4))
plt.plot(n_auto, signal1, label=f"{NAME1} ({FLOW1})")
plt.plot(n_auto, signal2, label=f"{NAME2} ({FLOW2})", alpha=0.7)
plt.title(f"Signals (z-scored) — {NAME1} ({FLOW1}) vs {NAME2} ({FLOW2})")
plt.xlabel("Position (px)")
plt.ylabel("Amplitude (σ)")
plt.legend()
plt.grid(True)

# ── Figure 2: wavelet panels ──────────────────────────────────────────────────
# Row 0: WPS1 + mean-power sidebar | WPS2 + mean-power sidebar
# Row 1: WCS bilateral
# Row 2: WCO bilateral
fig = plt.figure(figsize=(20, 13))
gs  = GridSpec(3, 4, figure=fig, width_ratios=[9, 2.5, 9, 2.5],
               hspace=0.35, wspace=0.35)

# — WPS1 —
ax_wps1 = fig.add_subplot(gs[0, 0])
im0 = ax_wps1.pcolormesh(LL, FF, np.abs(WPS1), cmap="viridis", shading="auto")
fig.colorbar(im0, ax=ax_wps1, label="Magnitude", fraction=0.046, pad=0.04)
ax_wps1.set_title(f"WPS — {NAME1} ({FLOW1}) from $R_{{xx}}$")
ax_wps1.set_xlabel("Lag (px)")
ax_wps1.set_ylabel("Wavelength (px)")
apply_freq_axis(ax_wps1, coi_wavelengths=[coi_boundary1])

ax_mp1 = fig.add_subplot(gs[0, 1], sharey=ax_wps1)
ax_mp1.plot(mean_power1, freqs, color="tab:blue")
ax_mp1.set_xlim(0, power_xlim)
ax_mp1.xaxis.set_major_locator(MaxNLocator(3))
ax_mp1.xaxis.set_major_formatter(EngFormatter())
ax_mp1.set_title("Mean power")
ax_mp1.set_xlabel("Mean |WPS|")
apply_freq_axis(ax_mp1, coi_wavelengths=[coi_boundary1])
plt.setp(ax_mp1.get_yticklabels(), visible=False)

# — WPS2 —
ax_wps2 = fig.add_subplot(gs[0, 2])
im1 = ax_wps2.pcolormesh(LL, FF, np.abs(WPS2), cmap="viridis", shading="auto")
fig.colorbar(im1, ax=ax_wps2, label="Magnitude", fraction=0.046, pad=0.04)
ax_wps2.set_title(f"WPS — {NAME2} ({FLOW2}) from $R_{{yy}}$")
ax_wps2.set_xlabel("Lag (px)")
ax_wps2.set_ylabel("Wavelength (px)")
apply_freq_axis(ax_wps2, coi_wavelengths=[coi_boundary2])

ax_mp2 = fig.add_subplot(gs[0, 3], sharey=ax_wps2)
ax_mp2.plot(mean_power2, freqs, color="tab:blue")
ax_mp2.set_xlim(0, power_xlim)
ax_mp2.xaxis.set_major_locator(MaxNLocator(3))
ax_mp2.xaxis.set_major_formatter(EngFormatter())
ax_mp2.set_title("Mean power")
ax_mp2.set_xlabel("Mean |WPS|")
apply_freq_axis(ax_mp2, coi_wavelengths=[coi_boundary2])
plt.setp(ax_mp2.get_yticklabels(), visible=False)

# — WCS bilateral —
ax_wcs = fig.add_subplot(gs[1, 0:4])
im2 = ax_wcs.pcolormesh(LL, FF, np.abs(WCS), cmap="coolwarm", shading="auto")
fig.colorbar(im2, ax=ax_wcs, label="Cross Power", fraction=0.023, pad=0.02)
ax_wcs.set_title(
    f"Wavelet Cross Spectrum — {NAME1} ({FLOW1}) x {NAME2} ({FLOW2})  "
    f"[lag < 0: {NAME2} leads | lag > 0: {NAME1} leads]"
)
ax_wcs.axvline(0, color="white", linewidth=0.8, linestyle="--")
ax_wcs.set_xlabel("Lag (px)")
ax_wcs.set_ylabel("Wavelength (px)")
apply_freq_axis(ax_wcs, coi_wavelengths=[coi_boundary1, coi_boundary2])

# — WCO bilateral —
ax_wco = fig.add_subplot(gs[2, 0:4])
im3 = ax_wco.pcolormesh(LL, FF, WCO, cmap="plasma", shading="auto", vmin=0, vmax=1)
fig.colorbar(im3, ax=ax_wco, label="Coherence", fraction=0.023, pad=0.02)
ax_wco.set_title(
    f"Wavelet Coherence — {NAME1} ({FLOW1}) x {NAME2} ({FLOW2})  "
    f"[lag < 0: {NAME2} leads | lag > 0: {NAME1} leads]"
)
ax_wco.axvline(0, color="white", linewidth=0.8, linestyle="--")
ax_wco.set_xlabel("Lag (px)")
ax_wco.set_ylabel("Wavelength (px)")
apply_freq_axis(ax_wco, coi_wavelengths=[coi_boundary1, coi_boundary2])

# ── Save or show ──────────────────────────────────────────────────────────────
if SHOW:
    plt.show()
else:
    out_dir = os.path.join(HERE, "%s_%s__x__%s_%s" % (NAME1, FLOW1, NAME2, FLOW2))
    os.makedirs(out_dir, exist_ok=True)
    fig_sig.savefig(os.path.join(out_dir, "signals.png"), dpi=120, bbox_inches="tight")
    fig.savefig(os.path.join(out_dir, "wavelet.png"),     dpi=120, bbox_inches="tight")
    plt.close("all")
    print("saved: %s" % out_dir)
