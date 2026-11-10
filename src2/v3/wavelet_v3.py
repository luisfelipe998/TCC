import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, MaxNLocator, EngFormatter
from matplotlib.gridspec import GridSpec
from scipy.signal import butter, sosfiltfilt
import pandas as pd
import pywt

# --- Signal source ---
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.join(HERE, "../", "data_processed")
OS_COL = "Voltage_2"
sampling_rate = 1000  # Hz

FREQ_MIN, FREQ_MAX = 0.5, 60
FREQ_TICKS = [0.5, 1, 1.5, 3, 5, 10, 20, 30, 40, 60]

# --- Preprocessing parameters ---
# Step 1: Butterworth bandpass [FREQ_MIN, FREQ_MAX] Hz, order 4.
# Removes camera noise above 60 Hz and slow drift below 0.5 Hz before correlation.
# sosfiltfilt (zero-phase) avoids phase distortion.
BUTTER_ORDER = 4

# Step 2: Normalize by std after DC removal and bandpass.
# Makes WCO insensitive to amplitude differences between runs (illumination, contrast).


def bandpass(sig, fs, lo, hi, order=BUTTER_ORDER):
    nyq = fs / 2.0
    sos = butter(order, [lo / nyq, hi / nyq], btype='band', output='sos')
    return sosfiltfilt(sos, sig)


def load_signal(name, fs):
    """Le Voltage_2, remove NaN, remove DC, aplica bandpass e normaliza por std."""
    path = os.path.join(DATA_ROOT, name, "dados.xlsx")
    df = pd.read_excel(path, engine="openpyxl")
    if OS_COL not in df.columns:
        raise KeyError("Canal '%s' ausente em %s. Colunas: %s" % (OS_COL, path, list(df.columns)))
    sig = pd.to_numeric(df[OS_COL], errors="coerce").to_numpy(float)
    sig = sig[~np.isnan(sig)]
    sig = sig - np.mean(sig)
    sig = bandpass(sig, fs, FREQ_MIN, FREQ_MAX)
    std = np.std(sig)
    if std > 0:
        sig = sig / std
    return sig


def load_fs(name):
    """Le effective_fs da sheet metadata de dados.xlsx. Fallback: sampling_rate global."""
    path = os.path.join(DATA_ROOT, name, "dados.xlsx")
    try:
        df_meta = pd.read_excel(path, sheet_name="metadata", engine="openpyxl")
        row = df_meta[df_meta["key"] == "effective_fs"]
        if not row.empty:
            return float(row["value"].iloc[0])
    except Exception:
        pass
    return sampling_rate


def folder_flow_type(name):
    folder = os.path.join(DATA_ROOT, name)
    for p in sorted(glob.glob(os.path.join(folder, "*.txt"))):
        stem = os.path.splitext(os.path.basename(p))[0].lstrip("_")
        if stem != "dados":
            return stem
    return name


SHOW = "--show" in sys.argv
_pos = [a for a in sys.argv[1:] if not a.startswith("--")]
NAME1 = _pos[0] if len(_pos) > 0 else "synthetic1"
NAME2 = _pos[1] if len(_pos) > 1 else "synthetic2"
FLOW1 = folder_flow_type(NAME1)
FLOW2 = folder_flow_type(NAME2)

# Use effective_fs from metadata if available; both runs should have the same fs.
# If they differ (mixed synthetic + real), use the average and warn.
fs1 = load_fs(NAME1)
fs2 = load_fs(NAME2)
if fs1 != fs2:
    print(f"WARNING: fs mismatch — {NAME1}={fs1:.2f} Hz, {NAME2}={fs2:.2f} Hz. Using average.")
fs = (fs1 + fs2) / 2.0
print(f"Sampling rate: {fs:.2f} Hz")

signal1 = load_signal(NAME1, fs)
signal2 = load_signal(NAME2, fs)

N = min(len(signal1), len(signal2))
signal1 = signal1[:N]
signal2 = signal2[:N]

# --- Correlacoes (estimador enviesado, divide por N) ---
# Biased: divide por N, suprimindo amplificacao de ruido em lags grandes.
#   R_xx(m) = (1/N) * sum_{n=0}^{N-1-m} x[n]*x[n+m],  0 <= m <= N-1
def biased_autocorrelation(x):
    return np.array([np.sum(x[:N - m] * x[m:]) / N for m in range(N)])

def biased_crosscorrelation(x, y):
    return np.array([np.sum(x[:N - m] * y[m:]) / N for m in range(N)])

Rxx = biased_autocorrelation(signal1)
Ryy = biased_autocorrelation(signal2)
Rxy = biased_crosscorrelation(signal1, signal2)
Ryx = biased_crosscorrelation(signal2, signal1)

# Correlacao cruzada bilateral: lags -(N-1)..0..+(N-1).
# Lags negativos (NAME2 lidera): Ryx[N-1..1] revertido.
# Lags positivos (NAME1 lidera): Rxy[0..N-1].
lags = np.arange(-(N - 1), N)
Rxy_full = np.concatenate([Ryx[N - 1:0:-1], Rxy])

# --- Escalas CWT (cmor2.0-0.9549, omega_0 ~ 6.0) ---
wavelet = 'cmor2.0-0.9549'
center_freq = pywt.central_frequency(wavelet)
freqs = np.logspace(np.log10(FREQ_MIN), np.log10(FREQ_MAX), 200)
scales = center_freq * fs / freqs

# --- CWT ---
WPS1, _ = pywt.cwt(Rxx.real, scales, wavelet, sampling_period=1 / fs)
WPS2, _ = pywt.cwt(Ryy.real, scales, wavelet, sampling_period=1 / fs)
WCS,  _ = pywt.cwt(Rxy_full.real, scales, wavelet, sampling_period=1 / fs)

# --- WCO (Schuck & Bodmann, eq. 32.17) ---
# Normaliza pelo WPS no lag zero (b=0, pico da autocorrelacao).
WCO = np.abs(WCS) ** 2 / (np.abs(WPS1[:, [0]]) * np.abs(WPS2[:, [0]]))

# --- Potencia media por frequencia ---
mean_power1 = np.mean(np.abs(WPS1), axis=1)
mean_power2 = np.mean(np.abs(WPS2), axis=1)
POWER_MARGIN = 1.1
power_xlim = max(mean_power1.max(), mean_power2.max()) * POWER_MARGIN

# Malhas para pcolormesh.
n_auto = np.arange(N)
NN_auto, F_auto = np.meshgrid(n_auto, freqs)
LL, F_cross = np.meshgrid(lags, freqs)


def apply_freq_axis(ax):
    ax.set_yscale('log')
    ax.set_ylim(FREQ_MIN, FREQ_MAX)
    ax.yaxis.set_major_locator(FixedLocator(FREQ_TICKS))
    ax.yaxis.set_major_formatter(FixedFormatter([str(x) for x in FREQ_TICKS]))
    ax.yaxis.set_minor_locator(FixedLocator([]))


# --- Sinais preprocessados (figura separada) ---
fig_sig = plt.figure(figsize=(10, 4))
plt.plot(n_auto, signal1, label=f'{NAME1} ({FLOW1})')
plt.plot(n_auto, signal2, label=f'{NAME2} ({FLOW2})', alpha=0.7)
plt.title(f"Signals (bandpass + normalized) — {NAME1} ({FLOW1}) vs {NAME2} ({FLOW2})")
plt.xlabel("Sample")
plt.ylabel("Amplitude (normalized)")
plt.legend()
plt.grid(True)

# --- Figura principal: 3 linhas ---
# Linha 0: WPS1 + painel lateral | WPS2 + painel lateral
# Linha 1: WCS bilateral
# Linha 2: WCO bilateral
fig = plt.figure(figsize=(20, 13))
gs = GridSpec(3, 4, figure=fig, width_ratios=[9, 2.5, 9, 2.5],
              hspace=0.35, wspace=0.35)

# WPS1
ax_wps1 = fig.add_subplot(gs[0, 0])
im0 = ax_wps1.pcolormesh(NN_auto, F_auto, np.abs(WPS1), cmap='viridis', shading='auto')
fig.colorbar(im0, ax=ax_wps1, label='Magnitude', fraction=0.046, pad=0.04)
ax_wps1.set_title(f"WPS — {NAME1} ({FLOW1}) from $R_{{xx}}$")
ax_wps1.set_xlabel("Lag (samples)")
ax_wps1.set_ylabel("Frequency [Hz]")
apply_freq_axis(ax_wps1)

ax_mp1 = fig.add_subplot(gs[0, 1], sharey=ax_wps1)
ax_mp1.plot(mean_power1, freqs, color="tab:blue")
ax_mp1.set_xlim(0, power_xlim)
ax_mp1.xaxis.set_major_locator(MaxNLocator(3))
ax_mp1.xaxis.set_major_formatter(EngFormatter())
ax_mp1.set_title("Mean power")
ax_mp1.set_xlabel("Mean |WPS|")
apply_freq_axis(ax_mp1)
plt.setp(ax_mp1.get_yticklabels(), visible=False)

# WPS2
ax_wps2 = fig.add_subplot(gs[0, 2])
im1 = ax_wps2.pcolormesh(NN_auto, F_auto, np.abs(WPS2), cmap='viridis', shading='auto')
fig.colorbar(im1, ax=ax_wps2, label='Magnitude', fraction=0.046, pad=0.04)
ax_wps2.set_title(f"WPS — {NAME2} ({FLOW2}) from $R_{{yy}}$")
ax_wps2.set_xlabel("Lag (samples)")
ax_wps2.set_ylabel("Frequency [Hz]")
apply_freq_axis(ax_wps2)

ax_mp2 = fig.add_subplot(gs[0, 3], sharey=ax_wps2)
ax_mp2.plot(mean_power2, freqs, color="tab:blue")
ax_mp2.set_xlim(0, power_xlim)
ax_mp2.xaxis.set_major_locator(MaxNLocator(3))
ax_mp2.xaxis.set_major_formatter(EngFormatter())
ax_mp2.set_title("Mean power")
ax_mp2.set_xlabel("Mean |WPS|")
apply_freq_axis(ax_mp2)
plt.setp(ax_mp2.get_yticklabels(), visible=False)

# WCS bilateral
ax_wcs = fig.add_subplot(gs[1, 0:4])
im2 = ax_wcs.pcolormesh(LL, F_cross, np.abs(WCS), cmap='coolwarm', shading='auto')
fig.colorbar(im2, ax=ax_wcs, label='Cross Power', fraction=0.023, pad=0.02)
ax_wcs.set_title(
    f"Wavelet Cross Spectrum — {NAME1} ({FLOW1}) x {NAME2} ({FLOW2})  "
    f"[lag < 0: {NAME2} leads | lag > 0: {NAME1} leads]"
)
ax_wcs.axvline(0, color='white', linewidth=0.8, linestyle='--')
ax_wcs.set_xlabel("Lag (samples)")
ax_wcs.set_ylabel("Frequency [Hz]")
apply_freq_axis(ax_wcs)

# WCO bilateral
ax_wco = fig.add_subplot(gs[2, 0:4])
im3 = ax_wco.pcolormesh(LL, F_cross, WCO, cmap='plasma', shading='auto', vmin=0, vmax=1)
fig.colorbar(im3, ax=ax_wco, label='Coherence', fraction=0.023, pad=0.02)
ax_wco.set_title(
    f"Wavelet Coherence — {NAME1} ({FLOW1}) x {NAME2} ({FLOW2})  "
    f"[lag < 0: {NAME2} leads | lag > 0: {NAME1} leads]"
)
ax_wco.axvline(0, color='white', linewidth=0.8, linestyle='--')
ax_wco.set_xlabel("Lag (samples)")
ax_wco.set_ylabel("Frequency [Hz]")
apply_freq_axis(ax_wco)


if SHOW:
    plt.show()
else:
    out_dir = os.path.join(HERE, "%s_%s__x__%s_%s" % (NAME1, FLOW1, NAME2, FLOW2))
    os.makedirs(out_dir, exist_ok=True)
    fig_sig.savefig(os.path.join(out_dir, "signals.png"), dpi=120, bbox_inches="tight")
    fig.savefig(os.path.join(out_dir, "wavelet.png"), dpi=120, bbox_inches="tight")
    print("saved: %s" % out_dir)
