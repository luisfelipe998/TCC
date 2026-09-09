import os
import sys

import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, MaxNLocator, EngFormatter
from matplotlib.gridspec import GridSpec
import pandas as pd
import pywt

# --- Signal source (v8: v7 + estimador ENVIESADO da v4) ---
# Le sinais reais de src/data_processed (canal optico Voltage_2 da aba Sheet1,
# fs = 1 kHz). Igual a v7 (sinal 1 = media das autocorrelacoes de TODAS as runs de
# um padrao, hold-out da run usada como sinal 2), mas usando o estimador enviesado
# (divisor 1/N) da v4 no lugar do nao-enviesado (1/(N-m)).
# Construcao das entradas (todas truncadas ao menor N comum do grupo):
#   Rxx = R_bar_xx = (1/K) sum_k autocorr(run_k)           # media da classe (hold-in)
#   Ryy =            autocorr(holdout)                       # run individual crua
#   Rxy = (1/K) sum_k crosscorr(run_k, holdout)             # media dos cross pareados
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.join(HERE, "../", "data_processed")
OS_COL = "Voltage_2"
sampling_rate = 1000  # Hz

# Faixa de frequencia de interesse (Hz) e ticks a exibir.
FREQ_MIN, FREQ_MAX = 1.5, 128.0
FREQ_TICKS = [1.5, 3, 5, 10, 20, 30, 50, 70, 100, 128]


def load_signal(name):
    """Le o canal optico (Voltage_2) de data_processed/<name>/dados.xlsx.

    Remove o offset DC (media grande do OS domina a CWT); nao normaliza amplitude.
    """
    path = os.path.join(DATA_ROOT, name, "dados.xlsx")
    df = pd.read_excel(path, engine="openpyxl")
    if OS_COL not in df.columns:
        raise KeyError("Canal '%s' ausente em %s. Colunas: %s" % (OS_COL, path, list(df.columns)))
    sig = pd.to_numeric(df[OS_COL], errors="coerce").to_numpy(float)
    return sig - np.nanmean(sig)


def folder_flow_type(name):
    """Tipo de escoamento (ou 'synthetic') = nome do .txt marcador da pasta."""
    folder = os.path.join(DATA_ROOT, name)
    for p in sorted(glob.glob(os.path.join(folder, "*.txt"))):
        stem = os.path.splitext(os.path.basename(p))[0]
        if stem != "dados":
            return stem
    return name


def list_runs():
    """Pastas de data_processed que tem dados.xlsx (uma run cada)."""
    out = []
    for d in sorted(os.listdir(DATA_ROOT)):
        if os.path.exists(os.path.join(DATA_ROOT, d, "dados.xlsx")):
            out.append(d)
    return out


def group_members(pattern):
    """Runs de data_processed cujo padrao de escoamento e `pattern`."""
    return [r for r in list_runs() if folder_flow_type(r) == pattern]


# --- Autocorrelation and Cross-correlation (estimador ENVIESADO, divisor 1/N) ---
# Estimador enviesado da v4 (divisor fixo 1/N no lugar do nao-enviesado 1/(N-m)):
#   R_xx(m) = (1/N) * sum_n x[n]*conj(x[n+m]) ; idem R_xy com y no lugar de x.
# Como a soma em m tem so N-m termos, dividir por N (e nao por N-m) equivale a
# aplicar um peso triangular (N-m)/N que vai a zero quando m -> N. Assim a cauda
# de lag alto decai suavemente em vez de estourar em variancia (no estimador
# nao-enviesado, m ~ N deixa 1-2 termos sobre N-m ~ 1, gerando ruido que a CWT
# espalha na banda baixa e satura a WCO). O preco e um vies que subestima a
# magnitude nos lags altos; os lags pequenos/medios (onde vive o sinal) ficam ~iguais.
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


# --- Argumentos ---
# Uso: python wavelet_alternativa_v2.py <PATTERN> <HOLDOUT_RUN>
#   PATTERN     : padrao de escoamento (ex.: churn) — define o conjunto da media.
#   HOLDOUT_RUN : run individual usada como sinal 2 (mantida fora da media).
# Sem args: usa um default para inspecao rapida.
# Flag --show: comportamento interativo (plt.show). Sem ela, salva em subpasta.
SHOW = "--show" in sys.argv
_pos = [a for a in sys.argv[1:] if not a.startswith("--")]
PATTERN = _pos[0] if len(_pos) > 0 else "churn"
HOLDOUT = _pos[1] if len(_pos) > 1 else None

members = group_members(PATTERN)
if not members:
    raise SystemExit("Nenhuma run com padrao '%s' em data_processed." % PATTERN)
if HOLDOUT is None:
    HOLDOUT = members[0]
if HOLDOUT not in members:
    raise SystemExit(
        "Holdout '%s' nao pertence ao padrao '%s'. Membros: %s"
        % (HOLDOUT, PATTERN, members)
    )

# Conjunto da media = todos os membros do padrao EXCETO o holdout (hold-out).
hold_in = [m for m in members if m != HOLDOUT]
if not hold_in:
    raise SystemExit(
        "Padrao '%s' tem apenas 1 run (%s); sem hold-in para a media."
        % (PATTERN, HOLDOUT)
    )

FLOW1 = PATTERN                 # sinal 1 = assinatura media da classe
FLOW2 = folder_flow_type(HOLDOUT)  # sinal 2 = run individual (mesmo padrao)
NAME1 = "%s_mean_K%d" % (PATTERN, len(hold_in))
NAME2 = HOLDOUT

# Carrega hold-in e holdout; trunca tudo ao menor N comum do grupo.
sig_hold_in = [load_signal(m) for m in hold_in]
sig_holdout = load_signal(HOLDOUT)
num_samples = min([len(s) for s in sig_hold_in] + [len(sig_holdout)])
sig_hold_in = [s[:num_samples] for s in sig_hold_in]
sig_holdout = sig_holdout[:num_samples]

n = np.arange(num_samples)  # eixo de amostras

# Sinal 1: media das autocorrelacoes das runs hold-in (assinatura da classe,
# Schuck & Bodmann, eq. ~32.18):  Rbar_xx(m) = (1/K) * sum_{k=1}^{K} R_xx^(k)(m).
Rxx = np.mean([biased_autocorrelation(s) for s in sig_hold_in], axis=0)
# Sinal 2: autocorrelacao da run holdout crua.
Ryy = biased_autocorrelation(sig_holdout)
# Cross: media das cross-correlacoes pareadas (cada run hold-in x o holdout):
#   Rbar_xy(m) = (1/K) * sum_{k=1}^{K} R_xy^(k)(m).
Rxy = np.mean(
    [biased_crosscorrelation(s, sig_holdout) for s in sig_hold_in], axis=0
)

# --- Scales chosen so the CWT covers exactly FREQ_MIN..FREQ_MAX Hz ---
# Morlet complexa 'cmorB-C': a frequencia (em unidades de fs) de uma escala a e
# C/a, logo f[Hz] = C*fs/a  ->  a = C*fs/f. Escalas log-espacadas amostram a
# faixa de frequencia de forma uniforme em log.
wavelet = 'cmor2.0-0.9549'
center_freq = pywt.central_frequency(wavelet)          # C em unidades de fs
freqs = np.logspace(np.log10(FREQ_MIN), np.log10(FREQ_MAX), 200)  # Hz, crescente
scales = center_freq * sampling_rate / freqs           # escala correspondente

# Transformada wavelet continua: WPS de Rxx/Ryy (auto) e WCS de Rxy (cross).
#   W(a,b) = (1/sqrt(a)) * sum_t R(t) * conj(psi((t-b)/a))
WPS1_alt, _ = pywt.cwt(Rxx.real, scales, wavelet, sampling_period=1 / sampling_rate)
WPS2_alt, _ = pywt.cwt(Ryy.real, scales, wavelet, sampling_period=1 / sampling_rate)
WCS_alt, _ = pywt.cwt(Rxy.real, scales, wavelet, sampling_period=1 / sampling_rate)

# --- Coherence (Alternative WCO Definition) ---
# Coerencia wavelet normalizada pela potencia em b0=0 (primeira coluna):
#   WCO(a,b) = |W_xy(a,b)|^2 / ( |W_xx(a,b0)| * |W_yy(a,b0)| )
WCO_alt = np.abs(WCS_alt) ** 2 / (np.abs(WPS1_alt[:, [0]]) * np.abs(WPS2_alt[:, [0]]))

# Potencia media (media ao longo das amostras) de cada WPS, por frequencia.
mean_power1 = np.mean(np.abs(WPS1_alt), axis=1)
mean_power2 = np.mean(np.abs(WPS2_alt), axis=1)

# Limite comum do eixo x dos dois paineis laterais de "Mean power": pega o maior
# pico entre as duas curvas e aplica POWER_MARGIN (folga de 10%) para o pico nao
# encostar na borda direita. Ambos os paineis usam este mesmo limite, deixando as
# potencias das duas WPS comparaveis lado a lado. So afeta a escala do grafico.
POWER_MARGIN = 1.1
power_xlim = max(mean_power1.max(), mean_power2.max()) * POWER_MARGIN

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
plt.plot(n, Rxx.real, label='Class mean autocorr $\\bar{R}_{xx}$')
plt.plot(n, sig_holdout, label='Holdout run (raw)', alpha=0.7)
plt.title(f"Signals — {NAME1} ({FLOW1}) vs {NAME2} ({FLOW2})")
plt.xlabel("Sample")
plt.ylabel("Amplitude")
plt.legend()
plt.grid(True)

# --- Four wavelet charts on one image ---
fig = plt.figure(figsize=(20, 9))
gs = GridSpec(2, 4, figure=fig, width_ratios=[9, 2.5, 9, 2.5],
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
ax_wcs.set_xlabel("Sample")
ax_wcs.set_ylabel("Frequency [Hz]")
apply_freq_axis(ax_wcs)

# WCO
ax_wco = fig.add_subplot(gs[1, 2:4])
im3 = ax_wco.pcolormesh(NN, F, WCO_alt, cmap='plasma', shading='auto', vmin=0, vmax=1)
fig.colorbar(im3, ax=ax_wco, label='Coherence', fraction=0.046, pad=0.04)
ax_wco.set_title(f"Wavelet Coherence — {NAME1} ({FLOW1}) x {NAME2} ({FLOW2})")
ax_wco.set_xlabel("Sample")
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
