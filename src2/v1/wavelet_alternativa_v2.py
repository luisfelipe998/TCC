import os
import sys

import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, MaxNLocator, EngFormatter
from matplotlib.gridspec import GridSpec
import pandas as pd
import pywt

# --- Signal source ---
# Le dois sinais reais de src/data_processed: canal optico Voltage_2 da aba
# Sheet1, amostrados a fs = 1 kHz.
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.join(HERE, "../", "data_processed")
OS_COL = "Voltage_2"
sampling_rate = 1000  # Hz

# Faixa de frequencia de interesse (Hz) e ticks a exibir.
FREQ_MIN, FREQ_MAX = 0.5, 60
FREQ_TICKS = [0.5, 1, 1.5, 3, 5, 10, 20, 30, 40, 60]


def load_signal(name):
    """Le o canal optico (Voltage_2) de data_processed/<name>/dados.xlsx."""
    path = os.path.join(DATA_ROOT, name, "dados.xlsx")
    df = pd.read_excel(path, engine="openpyxl")
    if OS_COL not in df.columns:
        raise KeyError("Canal '%s' ausente em %s. Colunas: %s" % (OS_COL, path, list(df.columns)))
    sig = pd.to_numeric(df[OS_COL], errors="coerce").to_numpy(float)
    sig = sig[~np.isnan(sig)]
    # Remove o offset DC. O OS (Voltage_2) tem media grande (~-1) e flutuacao
    # pequena; sem remover a media a autocorrelacao preserva o DC, que domina a
    # CWT e faz toda WPS colapsar num borrao de baixa frequencia (todo padrao
    # picaria em ~1.7 Hz). Nao normaliza pela amplitude (mantem std original).
    return sig - np.mean(sig)


def folder_flow_type(name):
    """Tipo de escoamento (ou 'synthetic') = nome do .txt marcador da pasta."""
    folder = os.path.join(DATA_ROOT, name)
    for p in sorted(glob.glob(os.path.join(folder, "*.txt"))):
        stem = os.path.splitext(os.path.basename(p))[0].lstrip("_")
        if stem != "dados":
            return stem
    return name


# Pastas-alvo em data_processed: padrao synthetic1/synthetic2; se passados
# argumentos na linha de comando, usa-os (python wavelet_alternativa_v2.py <f1> <f2>).
# Flag --show: reproduz o comportamento interativo (plt.show). Sem ela, salva os
# graficos em uma subpasta NAME1_FLOW1__x__NAME2_FLOW2 ao lado do script.
SHOW = "--show" in sys.argv
_pos = [a for a in sys.argv[1:] if not a.startswith("--")]
NAME1 = _pos[0] if len(_pos) > 0 else "synthetic1"
NAME2 = _pos[1] if len(_pos) > 1 else "synthetic2"
FLOW1 = folder_flow_type(NAME1)
FLOW2 = folder_flow_type(NAME2)

signal1 = load_signal(NAME1)
signal2 = load_signal(NAME2)

# Quando os sinais tem tamanhos diferentes, trunca o maior ao tamanho do menor.
num_samples = min(len(signal1), len(signal2))
signal1 = signal1[:num_samples]
signal2 = signal2[:num_samples]

n = np.arange(num_samples)  # eixo de amostras

# --- Autocorrelation and Cross-correlation (estimador enviesado) ---
# Estimador enviesado: divide sempre por N (tamanho total), nao por (N-m).
# Isso suprime a amplificacao de ruido em lags grandes — o estimador nao-enviesado
# divide por (N-m), o que amplifica o produto de um unico par de amostras em m=N-1
# por um fator N, gerando a mancha brilhante no final da WPS. O estimador enviesado
# e consistente (converge com N->inf), tem menor variancia e e o padrao em numpy/scipy.
#
# Autocorrelacao enviesada:
#   R_xx(m) = (1/N) * sum_{n=0}^{N-1-m} x[n] * conj(x[n+m]),  0 <= m <= N-1
#
# Cross-correlacao enviesada:
#   R_xy(m) = (1/N) * sum_{n=0}^{N-1-m} x[n] * conj(y[n+m]),  0 <= m <= N-1
def biased_autocorrelation(x):
    N = len(x)
    return np.array([
        np.sum(x[:N - m] * np.conj(x[m:])) / N
        for m in range(N)
    ])

def biased_crosscorrelation(x, y):
    N = len(x)
    return np.array([
        np.sum(x[:N - m] * np.conj(y[m:])) / N
        for m in range(N)
    ])

Rxx = biased_autocorrelation(signal1)
Ryy = biased_autocorrelation(signal2)
Rxy = biased_crosscorrelation(signal1, signal2)
Ryx = biased_crosscorrelation(signal2, signal1)

# --- Scales chosen so the CWT covers exactly FREQ_MIN..FREQ_MAX Hz ---
# Morlet complexa 'cmorB-C': a frequencia (em unidades de fs) de uma escala a e
# C/a, logo f[Hz] = C*fs/a  ->  a = C*fs/f. Escalas log-espacadas amostram a
# faixa de frequencia de forma uniforme em log.
wavelet = 'cmor2.0-0.9549'
center_freq = pywt.central_frequency(wavelet)          # C em unidades de fs
freqs = np.logspace(np.log10(FREQ_MIN), np.log10(FREQ_MAX), 200)  # Hz, crescente
scales = center_freq * sampling_rate / freqs           # escala correspondente

# Transformada wavelet continua aplicada as correlacoes (Schuck & Bodmann, eq. 32.14-32.20).
# WPS = |CWT{R_xx}|: espectro de potencia wavelet do sinal x via autocorrelacao.
# WCS = CWT{R_xy}: espectro cruzado wavelet entre x e y via correlacao cruzada.
# Para sinais reais, R_xx e R_xy sao reais; .real e no-op mas explicita a intencao.
#   W(a,b) = (1/sqrt(a)) * sum_t R(t) * conj(psi((t-b)/a))
WPS1_alt, _ = pywt.cwt(Rxx.real, scales, wavelet, sampling_period=1 / sampling_rate)
WPS2_alt, _ = pywt.cwt(Ryy.real, scales, wavelet, sampling_period=1 / sampling_rate)
WCS_alt, _ = pywt.cwt(Rxy.real, scales, wavelet, sampling_period=1 / sampling_rate)
WCS_rev, _ = pywt.cwt(Ryx.real, scales, wavelet, sampling_period=1 / sampling_rate)

# --- Coherence (WCO — Schuck & Bodmann, eq. 32.17) ---
# Coerencia wavelet normalizada pela potencia no lag zero (b=0, primeira coluna):
#   WCO(a,b) = |WCS_xy(a,b)|^2 / ( |WPS_x(a,0)| * |WPS_y(a,0)| )
# A normalizacao em b=0 usa o pico da autocorrelacao (lag zero = maxima correlacao),
# garantindo WCO in [0,1] no ponto de normalizacao. Em outros lags pode exceder 1
# marginalmente; vmax=1 no pcolormesh satura esses valores.
WCO_alt = np.abs(WCS_alt) ** 2 / (np.abs(WPS1_alt[:, [0]]) * np.abs(WPS2_alt[:, [0]]))
WCO_rev = np.abs(WCS_rev) ** 2 / (np.abs(WPS2_alt[:, [0]]) * np.abs(WPS1_alt[:, [0]]))

# Potencia media (media ao longo das amostras) de cada WPS, por frequencia.
mean_power1 = np.mean(np.abs(WPS1_alt), axis=1)
mean_power2 = np.mean(np.abs(WPS2_alt), axis=1)

# Limite comum do eixo x dos paineis laterais de potencia media: maior pico entre
# as duas WPS com folga de 10%. Ambos os paineis usam este limite para que as
# potencias sejam comparaveis lado a lado.
POWER_MARGIN = 1.1
power_xlim = max(mean_power1.max(), mean_power2.max()) * POWER_MARGIN

# Malha lag x frequencia para pcolormesh (eixo x = lag em amostras, eixo y = Hz).
NN, F = np.meshgrid(n, freqs)
# Malha com lag negativo para os paineis reversed: lag negativo indica sinal2 liderando.
NN_rev, _ = np.meshgrid(-n, freqs)


def apply_freq_axis(ax):
    """Eixo y em frequencia (log), limitado a FREQ_MIN..FREQ_MAX, com ticks fixos."""
    ax.set_yscale('log')
    ax.set_ylim(FREQ_MIN, FREQ_MAX)
    ax.yaxis.set_major_locator(FixedLocator(FREQ_TICKS))
    ax.yaxis.set_major_formatter(FixedFormatter([str(x) for x in FREQ_TICKS]))
    ax.yaxis.set_minor_locator(FixedLocator([]))


# --- Signals (separate figure) ---
fig_sig = plt.figure(figsize=(10, 4))
plt.plot(n, signal1, label='Original Signal')
plt.plot(n, signal2, label='Delayed Signal', alpha=0.7)
plt.title(f"Signals — {NAME1} ({FLOW1}) vs {NAME2} ({FLOW2})")
plt.xlabel("Sample")
plt.ylabel("Amplitude")
plt.legend()
plt.grid(True)

# --- Tres linhas de graficos wavelet numa unica figura ---
# Linha 0: WPS1 (+ painel lateral de potencia media) | WPS2 (+ painel lateral).
# Linha 1: WCS forward (signal1 lidera) | WCO forward.
# Linha 2: WCS reversed (signal2 lidera, eixo x negativo) | WCO reversed.
# Cada WPS tem um eixo estreito a direita compartilhando o eixo de frequencia,
# com o mesmo limite de potencia (power_xlim, auto) nos dois paineis.
fig = plt.figure(figsize=(20, 13))
gs = GridSpec(3, 4, figure=fig, width_ratios=[9, 2.5, 9, 2.5],
              hspace=0.3, wspace=0.35)

# WPS1 + painel lateral
ax_wps1 = fig.add_subplot(gs[0, 0])
im0 = ax_wps1.pcolormesh(NN, F, np.abs(WPS1_alt), cmap='viridis', shading='auto')
fig.colorbar(im0, ax=ax_wps1, label='Magnitude', fraction=0.046, pad=0.04)
ax_wps1.set_title(f"WPS — {NAME1} ({FLOW1}) from $R_{{xx}}(t)$")
ax_wps1.set_xlabel("Sample")
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

# WPS2 + painel lateral
ax_wps2 = fig.add_subplot(gs[0, 2])
im1 = ax_wps2.pcolormesh(NN, F, np.abs(WPS2_alt), cmap='viridis', shading='auto')
fig.colorbar(im1, ax=ax_wps2, label='Magnitude', fraction=0.046, pad=0.04)
ax_wps2.set_title(f"WPS — {NAME2} ({FLOW2}) from $R_{{yy}}(t)$")
ax_wps2.set_xlabel("Sample")
ax_wps2.set_ylabel("Frequency [Hz]")
apply_freq_axis(ax_wps2)

ax_mp2 = fig.add_subplot(gs[0, 3], sharey=ax_wps2)
ax_mp2.plot(mean_power2, freqs, color="tab:blue")
ax_mp2.set_xlim(0, power_xlim)     # mesmo limite comum da WPS1
ax_mp2.xaxis.set_major_locator(MaxNLocator(3))
ax_mp2.xaxis.set_major_formatter(EngFormatter())
ax_mp2.set_title("Mean power")
ax_mp2.set_xlabel("Mean |WPS|")
apply_freq_axis(ax_mp2)
plt.setp(ax_mp2.get_yticklabels(), visible=False)

# WCS
ax_wcs = fig.add_subplot(gs[1, 0:2])
im2 = ax_wcs.pcolormesh(NN, F, np.abs(WCS_alt), cmap='coolwarm', shading='auto')
fig.colorbar(im2, ax=ax_wcs, label='Cross Power', fraction=0.046, pad=0.04)
ax_wcs.set_title(f"Wavelet Cross Spectrum — {NAME1} ({FLOW1}) x {NAME2} ({FLOW2}) from $R_{{xy}}(t)$")
ax_wcs.set_xlabel("Lag (samples, positive = signal1 leads)")
ax_wcs.set_ylabel("Frequency [Hz]")
apply_freq_axis(ax_wcs)

# WCO
ax_wco = fig.add_subplot(gs[1, 2:4])
im3 = ax_wco.pcolormesh(NN, F, WCO_alt, cmap='plasma', shading='auto', vmin=0, vmax=1)
fig.colorbar(im3, ax=ax_wco, label='Coherence', fraction=0.046, pad=0.04)
ax_wco.set_title(f"Wavelet Coherence — {NAME1} ({FLOW1}) x {NAME2} ({FLOW2})")
ax_wco.set_xlabel("Lag (samples, positive = signal1 leads)")
ax_wco.set_ylabel("Frequency [Hz]")
apply_freq_axis(ax_wco)

# WCS reversed
ax_wcs_rev = fig.add_subplot(gs[2, 0:2])
im4 = ax_wcs_rev.pcolormesh(NN_rev, F, np.abs(WCS_rev), cmap='coolwarm', shading='auto')
fig.colorbar(im4, ax=ax_wcs_rev, label='Cross Power', fraction=0.046, pad=0.04)
ax_wcs_rev.set_title(f"Wavelet Cross Spectrum (reversed) — {NAME2} ({FLOW2}) x {NAME1} ({FLOW1}) from $R_{{yx}}(t)$")
ax_wcs_rev.set_xlabel("Lag (samples, negative = signal2 leads)")
ax_wcs_rev.set_ylabel("Frequency [Hz]")
apply_freq_axis(ax_wcs_rev)

# WCO reversed
ax_wco_rev = fig.add_subplot(gs[2, 2:4])
im5 = ax_wco_rev.pcolormesh(NN_rev, F, WCO_rev, cmap='plasma', shading='auto', vmin=0, vmax=1)
fig.colorbar(im5, ax=ax_wco_rev, label='Coherence', fraction=0.046, pad=0.04)
ax_wco_rev.set_title(f"Wavelet Coherence (reversed) — {NAME2} ({FLOW2}) x {NAME1} ({FLOW1})")
ax_wco_rev.set_xlabel("Lag (samples, negative = signal2 leads)")
ax_wco_rev.set_ylabel("Frequency [Hz]")
apply_freq_axis(ax_wco_rev)


if SHOW:
    plt.show()
else:
    out_dir = os.path.join(HERE, "%s_%s__x__%s_%s" % (NAME1, FLOW1, NAME2, FLOW2))
    os.makedirs(out_dir, exist_ok=True)
    fig_sig.savefig(os.path.join(out_dir, "signals.png"), dpi=120, bbox_inches="tight")
    fig.savefig(os.path.join(out_dir, "wavelet.png"), dpi=120, bbox_inches="tight")
    print("saved: %s" % out_dir)
