import os
import sys

import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, MaxNLocator, EngFormatter
from matplotlib.gridspec import GridSpec
import pandas as pd
import pywt
from scipy.signal import butter, filtfilt

# --- Signal source ---
# Le dois sinais reais de src/data_processed: canal optico Voltage_0 da aba
# Sheet1, amostrados a fs = 1 kHz.
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.join(HERE, "../", "data_processed")
OS_COL = "Voltage_0"
sampling_rate = 1000  # Hz

# Frequencia de corte do filtro passa-baixas (Hz) aplicado aos sinais antes da
# analise. Butterworth ordem 4, aplicado com fase zero (filtfilt).
LOWPASS_HZ = 50.0
LOWPASS_ORDER = 4

# Faixa de frequencia de interesse (Hz) e ticks a exibir.
FREQ_MIN, FREQ_MAX = 0.5, 128.0
FREQ_TICKS = [0.5, 1, 1.5, 3, 5, 10, 20, 30, 50, 70, 100, 128]


def load_signal(name):
    """Le o canal optico (Voltage_0) de data_processed/<name>/dados.xlsx."""
    path = os.path.join(DATA_ROOT, name, "dados.xlsx")
    df = pd.read_excel(path, engine="openpyxl")
    if OS_COL not in df.columns:
        raise KeyError("Canal '%s' ausente em %s. Colunas: %s" % (OS_COL, path, list(df.columns)))
    sig = pd.to_numeric(df[OS_COL], errors="coerce").to_numpy(float)
    # Remove o offset DC. O OS (Voltage_2) tem media grande (~-1) e flutuacao
    # pequena; sem remover a media a autocorrelacao preserva o DC, que domina a
    # CWT e faz toda WPS colapsar num borrao de baixa frequencia (todo padrao
    # picaria em ~1.7 Hz). Nao normaliza pela amplitude (mantem std original).
    return sig


def folder_flow_type(name):
    """Tipo de escoamento (ou 'synthetic') = nome do .txt marcador da pasta."""
    folder = os.path.join(DATA_ROOT, name)
    for p in sorted(glob.glob(os.path.join(folder, "*.txt"))):
        stem = os.path.splitext(os.path.basename(p))[0]
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

# --- Filtro passa-baixas (corte em LOWPASS_HZ) ---
# Butterworth ordem P: |H(f)|^2 = 1 / (1 + (f/fc)^(2P)), fc = LOWPASS_HZ.
# Aplicado com filtfilt (fase zero: passa nos dois sentidos, resposta |H(f)|^2 e
# sem atraso de fase entre os dois sinais). Remove o conteudo acima de fc antes
# de calcular correlacoes e CWT.
_b, _a = butter(LOWPASS_ORDER, LOWPASS_HZ / (sampling_rate / 2.0), btype="low")
signal1 = filtfilt(_b, _a, signal1)
signal2 = filtfilt(_b, _a, signal2)

n = np.arange(num_samples)  # eixo de amostras


# --- Autocorrelation and Cross-correlation (estimador ENVIESADO, divisor 1/N) ---
# Estimador enviesado (divisor fixo 1/N no lugar do nao-enviesado 1/(N-m)):
#   R_xx(m) = (1/N) * sum_n x[n]*conj(x[n+m]) ; idem R_xy com y no lugar de x.
# Como a soma em m tem so N-m termos, dividir por N (e nao por N-m) equivale a
# aplicar um peso triangular (N-m)/N que vai a zero quando m -> N. Assim a cauda
# de lag alto decai suavemente em vez de estourar em variancia (no estimador
# nao-enviesado, m ~ N deixa 1-2 termos sobre N-m ~ 1, gerando ruido que a CWT
# espalha na banda baixa e satura a WCO). O preco e um vies que subestima a
# magnitude nos lags altos; os lags pequenos/medios (onde vive o sinal) ficam ~iguais.
#
# O estimador so cobre lags m >= 0, ou seja, o "lado" em que signal1 lidera
# signal2. Para inspecionar o outro lado do eixo de lag basta trocar os sinais:
# R_yx(m) = conj(R_xy(-m)) e a metade de lag negativo refletida. Por isso a
# linha "reversed" da figura recalcula a WCS/WCO com signal1 e signal2 trocados,
# mantendo cada painel unilateral (o script nunca monta um eixo de lag bilateral).
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


# --- Scales chosen so the CWT covers exactly FREQ_MIN..FREQ_MAX Hz ---
# Morlet complexa 'cmorB-C': a frequencia (em unidades de fs) de uma escala a e
# C/a, logo f[Hz] = C*fs/a  ->  a = C*fs/f. Escalas log-espacadas amostram a
# faixa de frequencia de forma uniforme em log.
wavelet = 'cmor2.0-0.9549'
center_freq = pywt.central_frequency(wavelet)          # C em unidades de fs
freqs = np.logspace(np.log10(FREQ_MIN), np.log10(FREQ_MAX), 200)  # Hz, crescente
scales = center_freq * sampling_rate / freqs           # escala correspondente


def cwt_of(seq):
    """CWT da parte real de uma sequencia (auto/cross-correlacao)."""
    W, _ = pywt.cwt(seq.real, scales, wavelet, sampling_period=1 / sampling_rate)
    return W


def cross_terms(sig_x, sig_y):
    """WCS e WCO (definicao alternativa, normalizada em b0=0) do par (x, y).

    Ordem importa: R_xy(m) usa sig_x como o sinal que lidera (lag >= 0). Trocar
    x e y da o outro lado do eixo de lag (a metade de lag negativo refletida).
    """
    Rxy = biased_crosscorrelation(sig_x, sig_y)
    WCS = cwt_of(Rxy)
    WPSx = cwt_of(biased_autocorrelation(sig_x))
    WPSy = cwt_of(biased_autocorrelation(sig_y))
    # Coerencia wavelet normalizada pela potencia em b0=0 (primeira coluna):
    #   WCO(a,b) = |W_xy(a,b)|^2 / ( |W_xx(a,b0)| * |W_yy(a,b0)| )
    WCO = np.abs(WCS) ** 2 / (np.abs(WPSx[:, [0]]) * np.abs(WPSy[:, [0]]))
    return WCS, WCO


# Auto-espectros (WPS) das duas runs, na orientacao normal.
WPS1_alt = cwt_of(biased_autocorrelation(signal1))
WPS2_alt = cwt_of(biased_autocorrelation(signal2))

# Termos cruzados nos dois sentidos do eixo de lag.
WCS_alt, WCO_alt = cross_terms(signal1, signal2)          # signal1 lidera (lag >= 0)
WCS_rev, WCO_rev = cross_terms(signal2, signal1)          # signal2 lidera (lado reverso)

# Potencia media (media ao longo das amostras) de cada WPS, por frequencia.
mean_power1 = np.mean(np.abs(WPS1_alt), axis=1)
mean_power2 = np.mean(np.abs(WPS2_alt), axis=1)

# Limite comum do eixo x dos dois paineis laterais de "Mean power": pega o maior
# pico entre as duas curvas e aplica POWER_MARGIN (folga de 10%) para o pico nao
# encostar na borda direita. Ambos os paineis usam este mesmo limite, deixando as
# potencias das duas WPS comparaveis lado a lado. So afeta a escala do grafico.
POWER_MARGIN = 1.1
power_xlim = max(mean_power1.max(), mean_power2.max()) * POWER_MARGIN

# Limite comum de cor da WCS nos dois sentidos, para que a linha normal e a
# reversa fiquem comparaveis na mesma escala.
wcs_vmax = max(np.abs(WCS_alt).max(), np.abs(WCS_rev).max())

# Malha amostra x frequencia para pcolormesh (posicionamento correto em Hz).
NN, F = np.meshgrid(n, freqs)


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

# --- Wavelet charts on one image ---
# Linha 0: WPS1 (+ painel lateral de potencia media) | WPS2 (+ painel lateral).
# Linha 1: WCS | WCO na orientacao normal (signal1 lidera, lag >= 0).
# Linha 2: WCS | WCO reversos (signal2 lidera) — o outro lado do eixo de lag,
#          cada painel ainda unilateral. Cada WPS tem um eixo estreito a direita
#          compartilhando o eixo de frequencia, com o mesmo limite de potencia.
fig = plt.figure(figsize=(20, 13.5))
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

# WCS (orientacao normal: signal1 x signal2)
ax_wcs = fig.add_subplot(gs[1, 0:2])
im2 = ax_wcs.pcolormesh(NN, F, np.abs(WCS_alt), cmap='coolwarm', shading='auto', vmax=wcs_vmax)
fig.colorbar(im2, ax=ax_wcs, label='Cross Power', fraction=0.046, pad=0.04)
ax_wcs.set_title(f"Wavelet Cross Spectrum — {NAME1} ({FLOW1}) x {NAME2} ({FLOW2}) from $R_{{xy}}(t)$")
ax_wcs.set_xlabel("Sample")
ax_wcs.set_ylabel("Frequency [Hz]")
apply_freq_axis(ax_wcs)

# WCO (orientacao normal)
ax_wco = fig.add_subplot(gs[1, 2:4])
im3 = ax_wco.pcolormesh(NN, F, WCO_alt, cmap='plasma', shading='auto', vmin=0, vmax=1)
fig.colorbar(im3, ax=ax_wco, label='Coherence', fraction=0.046, pad=0.04)
ax_wco.set_title(f"Wavelet Coherence — {NAME1} ({FLOW1}) x {NAME2} ({FLOW2})")
ax_wco.set_xlabel("Sample")
ax_wco.set_ylabel("Frequency [Hz]")
apply_freq_axis(ax_wco)

# WCS reverso (signal2 x signal1 — o outro lado do eixo de lag)
ax_wcs_r = fig.add_subplot(gs[2, 0:2])
im4 = ax_wcs_r.pcolormesh(NN, F, np.abs(WCS_rev), cmap='coolwarm', shading='auto', vmax=wcs_vmax)
fig.colorbar(im4, ax=ax_wcs_r, label='Cross Power', fraction=0.046, pad=0.04)
ax_wcs_r.set_title(f"Wavelet Cross Spectrum (reversed) — {NAME2} ({FLOW2}) x {NAME1} ({FLOW1}) from $R_{{yx}}(t)$")
ax_wcs_r.set_xlabel("Sample")
ax_wcs_r.set_ylabel("Frequency [Hz]")
apply_freq_axis(ax_wcs_r)

# WCO reverso
ax_wco_r = fig.add_subplot(gs[2, 2:4])
im5 = ax_wco_r.pcolormesh(NN, F, WCO_rev, cmap='plasma', shading='auto', vmin=0, vmax=1)
fig.colorbar(im5, ax=ax_wco_r, label='Coherence', fraction=0.046, pad=0.04)
ax_wco_r.set_title(f"Wavelet Coherence (reversed) — {NAME2} ({FLOW2}) x {NAME1} ({FLOW1})")
ax_wco_r.set_xlabel("Sample")
ax_wco_r.set_ylabel("Frequency [Hz]")
apply_freq_axis(ax_wco_r)


if SHOW:
    plt.show()
else:
    out_dir = os.path.join(HERE, "%s_%s__x__%s_%s" % (NAME1, FLOW1, NAME2, FLOW2))
    os.makedirs(out_dir, exist_ok=True)
    fig_sig.savefig(os.path.join(out_dir, "signals.png"), dpi=120, bbox_inches="tight")
    fig.savefig(os.path.join(out_dir, "wavelet.png"), dpi=120, bbox_inches="tight")
    print("saved: %s" % out_dir)
