"""
extract_descriptors.py — descritor por RUN para a estatística mínima (A1/A2/A3).

Análise estatística sobre a base v8 (intra-padrão) / v8_1 (inter-padrão). NÃO é uma nova
versão de wavelet: apenas extrai UM descritor numérico por ensaio a partir do mesmo motor
Schuck & Bodmann de `wavelet_v8.py`, que permanece intocado.

Reuso do motor: `wavelet_v8.py` NÃO tem guard `if __name__ == "__main__"` — importá-lo executa
o pipeline de CLI no import e quebra. Por isso as funções necessárias são COPIADAS aqui
(byte-a-byte equivalentes) e um assert de identidade confirma, em runtime, que o WPS extraído
bate com o que o motor produziria.

Descritor (FIXO): perfil de escala do WPS auto na coluna de lag-zero b=0 (índice N-1 do WPS
bilateral), truncado a 7001 px, whitened POR ESCALA.

  Whitening (fórmula fixa): w = log(|WPS_b0| + eps); ajusta-se uma reta de w contra
  log(escala) (proxy do decaimento 1/f); o descritor é o RESÍDUO dessa reta. Remove o piso
  espectral comum e deixa só o desvio de forma por escala.

Saída: descriptors_cropped.npz e descriptors_center.npz com
  Z          (19, 200)  descritor whitened por run
  wps_b0     (19, 200)  |WPS b=0| cru (para a máscara de banda de A3)
  scales     (200,)     escalas CWT (px)
  labels     (19,)      padrão
  q, m       (19,)      fluxo de calor, vazão (do nome da pasta mXXqYY)
  shift_mean (19,)      shift médio de costura (metadata do xlsx)
  n_original (19,)      comprimento antes da truncagem
  runs       (19,)      nome da pasta
"""

import os
import re
import glob
import numpy as np
import pandas as pd
import pywt
from scipy.signal import detrend as scipy_detrend, fftconvolve

HERE      = os.path.dirname(os.path.abspath(__file__))
DATA_ROOT = os.path.join(HERE, "..", "data_processed")
OS_COL    = "Voltage_2"

# eixo de frequência espacial — IDÊNTICO ao motor (wavelet_v8.py L18-19, L165)
FREQ_MIN  = 1 / 3000
FREQ_MAX  = 1 / 50
WAVELET   = "cmor2.0-0.9549"

TRUNC_PX  = 7001                       # menor run de churn (m05q10_001)
PATTERNS  = {"churn", "wavy annular"}  # só os dois com n >= 2 relevantes à pergunta


# ── Funções COPIADAS do motor (wavelet_v8.py) ────────────────────────────────

def _load_and_preprocess(path):
    """Cópia de wavelet_v8._load_and_preprocess (L44-53): z-score único."""
    df = pd.read_excel(path, engine="openpyxl", sheet_name="Sheet")
    if OS_COL not in df.columns:
        raise KeyError("Canal '%s' ausente em %s." % (OS_COL, path))
    sig = pd.to_numeric(df[OS_COL], errors="coerce").to_numpy(float)
    sig = sig[~np.isnan(sig)]
    sig = scipy_detrend(sig, type="linear")
    sig = sig - sig.mean()
    sig = sig / sig.std(ddof=1)
    return sig


def folder_flow_type(name):
    """Cópia de wavelet_v8.folder_flow_type (L68-74)."""
    folder = os.path.join(DATA_ROOT, name)
    for p in sorted(glob.glob(os.path.join(folder, "*.txt"))):
        stem = os.path.splitext(os.path.basename(p))[0].lstrip("_")
        if stem != "dados":
            return stem
    return name


def _biased_autocorr_loop(x, N):
    """Cópia EXATA do laço O(N^2) do motor (wavelet_v8.compute_wavelet L146-147)."""
    return np.array([np.sum(x[:N - m] * x[m:]) / N for m in range(N)])


def _biased_autocorr_fft(x, N):
    """Autocorrelação enviesada via FFT (≡ laço direto, ~4e-16). Acelera O(N^2) → O(N log N)."""
    full = fftconvolve(x, x[::-1], mode="full")   # comprimento 2N-1, pico em índice N-1
    return full[N - 1:] / N                        # lags 0..N-1


def _wps_from_signal(sig, use_fft=True):
    """Reproduz o núcleo autocorr+CWT do motor e devolve (|WPS b=0| (200,), scales (200,)).

    Espelha wavelet_v8.compute_wavelet L143-171, mas para UM sinal (auto), extraindo só a
    coluna de lag-zero (índice N-1) do WPS bilateral.
    """
    N = len(sig)
    Rxx = _biased_autocorr_fft(sig, N) if use_fft else _biased_autocorr_loop(sig, N)
    Rxx_full = np.concatenate([Rxx[N - 1:0:-1], Rxx])
    assert Rxx_full.shape[0] == 2 * N - 1

    center_freq = pywt.central_frequency(WAVELET)
    freqs  = np.logspace(np.log10(FREQ_MIN), np.log10(FREQ_MAX), 200)
    scales = center_freq / freqs                    # sampling_rate = 1.0

    WPS, _ = pywt.cwt(Rxx_full.real, scales, WAVELET, sampling_period=1.0)
    wps_b0 = np.abs(WPS[:, N - 1])                   # coluna de lag-zero = b=0
    return wps_b0, scales


# ── Whitening por escala (fórmula fixa) ──────────────────────────────────────

def whiten_per_scale(wps_b0, scales):
    """w = log(|WPS|+eps); resíduo da reta w ~ log(escala) (remove piso 1/f)."""
    w = np.log(wps_b0 + 1e-12)
    ls = np.log(scales)
    A = np.vstack([ls, np.ones_like(ls)]).T
    coef, *_ = np.linalg.lstsq(A, w, rcond=None)
    return w - A @ coef


# ── Descoberta de runs ───────────────────────────────────────────────────────

def list_runs():
    runs = []
    for d in sorted(os.listdir(DATA_ROOT)):
        path = os.path.join(DATA_ROOT, d)
        if (os.path.isdir(path)
                and os.path.exists(os.path.join(path, "__dados_cropped_v11.xlsx"))
                and os.path.exists(os.path.join(path, "__dados_center_v13.xlsx"))):
            runs.append(d)
    return runs


def read_shift_mean(path):
    """Lê shift_mean_px da aba metadata do xlsx."""
    try:
        meta = pd.read_excel(path, engine="openpyxl", sheet_name="metadata", header=None)
    except Exception:
        return np.nan
    for _, row in meta.iterrows():
        if str(row.iloc[0]) == "shift_mean_px":
            return float(row.iloc[1])
    return np.nan


def parse_m_q(run_name):
    """mXXqYY_ZZZ -> (m, q). Ex.: m05q10_001 -> (5, 10)."""
    mobj = re.match(r"m0?(\d+)q0?(\d+)", run_name)
    if not mobj:
        return np.nan, np.nan
    return int(mobj.group(1)), int(mobj.group(2))


# ── Extração por variante (cropped / center) ─────────────────────────────────

def extract_variant(runs, xlsx_name, identity_checked):
    Z_rows, b0_rows, labels, qs, ms, shifts, n_orig, names = [], [], [], [], [], [], [], []
    scales_ref = None

    for r in runs:
        pattern = folder_flow_type(r)
        if pattern not in PATTERNS:
            continue
        xlsx = os.path.join(DATA_ROOT, r, xlsx_name)
        sig = _load_and_preprocess(xlsx)
        n0 = len(sig)
        sig = sig[:TRUNC_PX]                          # truncagem a 7001 px

        wps_b0, scales = _wps_from_signal(sig, use_fft=True)
        if scales_ref is None:
            scales_ref = scales

        # Assert de identidade: FFT-autocorr ≡ laço direto do motor (uma vez, no 1º run).
        if not identity_checked["done"]:
            wps_loop, _ = _wps_from_signal(sig, use_fft=False)
            err = np.max(np.abs(wps_b0 - wps_loop))
            assert err < 1e-9, "WPS b=0 (FFT) diverge do laço direto do motor: %g" % err
            print("  [assert] identidade FFT vs laço direto do motor OK (erro max=%.2e)" % err)
            identity_checked["done"] = True

        m, q = parse_m_q(r)
        Z_rows.append(whiten_per_scale(wps_b0, scales))
        b0_rows.append(wps_b0)
        labels.append(pattern); qs.append(q); ms.append(m)
        shifts.append(read_shift_mean(xlsx)); n_orig.append(n0); names.append(r)
        print("  %-16s %-14s q=%-3s m=%-3s n0=%-6d cortados=%d"
              % (r, pattern, q, m, n0, max(0, n0 - TRUNC_PX)))

    return dict(
        Z=np.array(Z_rows), wps_b0=np.array(b0_rows), scales=scales_ref,
        labels=np.array(labels), q=np.array(qs, float), m=np.array(ms, float),
        shift_mean=np.array(shifts, float), n_original=np.array(n_orig), runs=np.array(names),
    )


def main():
    runs = list_runs()
    print("runs com ambos xlsx:", len(runs))
    identity_checked = {"done": False}

    for variant, xlsx_name in [("cropped", "__dados_cropped_v11.xlsx"),
                               ("center",  "__dados_center_v13.xlsx")]:
        print("\n=== variante %s ===" % variant)
        d = extract_variant(runs, xlsx_name, identity_checked)
        out = os.path.join(HERE, "descriptors_%s.npz" % variant)
        np.savez(out, **d)
        n_ch = int((d["labels"] == "churn").sum())
        n_an = int((d["labels"] == "wavy annular").sum())
        print("salvo %s | churn=%d wavy annular=%d | Z shape=%s"
              % (out, n_ch, n_an, d["Z"].shape))
        assert n_ch == 10 and n_an == 9, "esperado churn=10, annular=9; obtido %d/%d" % (n_ch, n_an)


if __name__ == "__main__":
    main()
