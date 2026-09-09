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

# --- Signal source (v12: v10 + both-sides lag, like v5.1) ---
# Le sinais reais de src/data_processed (canal optico Voltage_2 da aba Sheet1,
# fs = 1 kHz). Igual a v10 (media das autocorrelacoes de TODAS as runs de um
# padrao, hold-out da run usada como sinal 2; estimador enviesado 1/N; filtro
# passa-baixas 50 Hz), mas calcula WCS/WCO nos dois sentidos do eixo de lag
# (sinal 1 lidera E sinal 2 lidera), como v5.1.
# Construcao das entradas (todas truncadas ao menor N comum do grupo):
#   Rxx = R_bar_xx = (1/K) sum_k autocorr(run_k)           # media da classe (hold-in)
#   Ryy =            autocorr(holdout)                       # run individual crua
#   Rxy = (1/K) sum_k crosscorr(run_k, holdout)             # media dos cross pareados
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.join(HERE, "../", "data_processed")
OS_COL = "Voltage_2"
sampling_rate = 1000  # Hz

# Frequencia de corte do filtro passa-baixas (Hz) aplicado aos sinais antes da
# analise. Butterworth ordem 4, aplicado com fase zero (filtfilt).
LOWPASS_HZ = 50.0
LOWPASS_ORDER = 4

# Faixa de frequencia de interesse (Hz) e ticks a exibir.
FREQ_MIN, FREQ_MAX = 0.5, 128.0
FREQ_TICKS = [0.5, 1, 1.5, 3, 5, 10, 20, 30, 50, 70, 100, 128]


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

FLOW1 = PATTERN                    # sinal 1 = assinatura media da classe
FLOW2 = folder_flow_type(HOLDOUT)  # sinal 2 = run individual (mesmo padrao)
NAME1 = "%s_mean_K%d" % (PATTERN, len(hold_in))
NAME2 = HOLDOUT

# Carrega hold-in e holdout; trunca tudo ao menor N comum do grupo.
sig_hold_in = [load_signal(m) for m in hold_in]
sig_holdout = load_signal(HOLDOUT)
num_samples = min([len(s) for s in sig_hold_in] + [len(sig_holdout)])
sig_hold_in = [s[:num_samples] for s in sig_hold_in]
sig_holdout = sig_holdout[:num_samples]

# --- Filtro passa-baixas (corte em LOWPASS_HZ) ---
# Butterworth ordem P: |H(f)|^2 = 1 / (1 + (f/fc)^(2P)), fc = LOWPASS_HZ.
# Aplicado com filtfilt (fase zero). Remove conteudo acima de fc de cada run
# (hold-in e holdout) antes de formar as autocorrelacoes e a media de classe.
_b, _a = butter(LOWPASS_ORDER, LOWPASS_HZ / (sampling_rate / 2.0), btype="low")
sig_hold_in = [filtfilt(_b, _a, s) for s in sig_hold_in]
sig_holdout = filtfilt(_b, _a, sig_holdout)

n = np.arange(num_samples)  # eixo de amostras

# Sinal 1: media das autocorrelacoes das runs hold-in (assinatura da classe,
# Schuck & Bodmann, eq. ~32.18):  Rbar_xx(m) = (1/K) * sum_{k=1}^{K} R_xx^(k)(m).
Rxx = np.mean([biased_autocorrelation(s) for s in sig_hold_in], axis=0)
# Sinal 2: autocorrelacao da run holdout crua.
Ryy = biased_autocorrelation(sig_holdout)

# --- Scales chosen so the CWT covers exactly FREQ_MIN..FREQ_MAX Hz ---
wavelet = 'cmor2.0-0.9549'
center_freq = pywt.central_frequency(wavelet)
freqs = np.logspace(np.log10(FREQ_MIN), np.log10(FREQ_MAX), 200)
scales = center_freq * sampling_rate / freqs


def cwt_of(seq):
    """CWT da parte real de uma sequencia (auto/cross-correlacao)."""
    W, _ = pywt.cwt(seq.real, scales, wavelet, sampling_period=1 / sampling_rate)
    return W


def wps_mean(sig_list):
    """WPS da classe: media das magnitudes das CWTs individuais (eq. 32.18).

    (1/K) * sum_k |CWT{R_xx^k}| — cada run tem sua CWT calculada separadamente,
    depois as magnitudes sao mediadas. Difere de |CWT{media das correlacoes}|
    quando as runs tem energia em posicoes de lag distintas.
    """
    return np.mean([np.abs(cwt_of(biased_autocorrelation(s))) for s in sig_list], axis=0)


def wcs_mean(sig_x_list, sig_y):
    """WCS: media das magnitudes das CWTs individuais de cada cross-correlacao (eq. 32.20).

    (1/K) * sum_k |CWT{R̂_xy^k}| — averaging incoerente, conforme o paper.
    Difere de |CWT{media das cross-correlacoes}| quando as runs tem estrutura
    de fase variavel entre si.
    """
    return np.mean([np.abs(cwt_of(biased_crosscorrelation(s, sig_y))) for s in sig_x_list], axis=0)


def cross_terms_mean(sig_hold_in_list, sig_holdout_arr, forward=True):
    """WCS e WCO (definicao alternativa, normalizada em b0=0).

    forward=True : media_classe lidera (lag >= 0) — cada run hold-in x holdout.
    forward=False: holdout lidera (lado reverso) — holdout x cada run hold-in.
    WPS calculado por wps_mean (eq. 32.18); WCS por wcs_mean (eq. 32.20).
    Ambos usam averaging incoerente (media das magnitudes), conforme o paper.
    """
    if forward:
        WCS = wcs_mean(sig_hold_in_list, sig_holdout_arr)
        WPSx = wps_mean(sig_hold_in_list)
        WPSy = np.abs(cwt_of(biased_autocorrelation(sig_holdout_arr)))
    else:
        WCS = np.mean(
            [np.abs(cwt_of(biased_crosscorrelation(sig_holdout_arr, s))) for s in sig_hold_in_list],
            axis=0,
        )
        WPSx = np.abs(cwt_of(biased_autocorrelation(sig_holdout_arr)))
        WPSy = wps_mean(sig_hold_in_list)
    WCO = WCS ** 2 / (WPSx[:, [0]] * WPSy[:, [0]])
    return WCS, WCO


# Auto-espectros (WPS).
# WPS1: media das magnitudes das CWTs individuais das runs hold-in (eq. 32.18).
# WPS2: CWT da autocorrelacao da run holdout (unica realizacao).
WPS1_alt = wps_mean(sig_hold_in)
WPS2_alt = np.abs(cwt_of(Ryy))

# Termos cruzados nos dois sentidos do eixo de lag.
WCS_alt, WCO_alt = cross_terms_mean(sig_hold_in, sig_holdout, forward=True)
WCS_rev, WCO_rev = cross_terms_mean(sig_hold_in, sig_holdout, forward=False)

# Potencia media (media ao longo das amostras) de cada WPS, por frequencia.
mean_power1 = np.mean(np.abs(WPS1_alt), axis=1)
mean_power2 = np.mean(np.abs(WPS2_alt), axis=1)

POWER_MARGIN = 1.1
power_xlim = max(mean_power1.max(), mean_power2.max()) * POWER_MARGIN

# Escala de cor comum para WCS nos dois sentidos.
wcs_vmax = max(np.abs(WCS_alt).max(), np.abs(WCS_rev).max())

# Malha amostra x frequencia para pcolormesh.
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

# --- Six wavelet charts (3 rows x 2 panels each) ---
# Linha 0: WPS1 (+ painel lateral de potencia media) | WPS2 (+ painel lateral).
# Linha 1: WCS | WCO na orientacao normal (classe lidera, lag >= 0).
# Linha 2: WCS | WCO reversos (holdout lidera) — o outro lado do eixo de lag.
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
ax_mp2.set_xlim(0, power_xlim)
ax_mp2.xaxis.set_major_locator(MaxNLocator(3))
ax_mp2.xaxis.set_major_formatter(EngFormatter())
ax_mp2.set_title("Mean power")
ax_mp2.set_xlabel("Mean |WPS|")
apply_freq_axis(ax_mp2)
plt.setp(ax_mp2.get_yticklabels(), visible=False)

# WCS (orientacao normal)
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

# WCS reverso (holdout x classe — outro lado do eixo de lag)
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
