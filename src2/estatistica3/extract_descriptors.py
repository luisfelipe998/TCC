"""
extract_descriptors.py — descritor por RUN, estendido aos QUATRO padrões (churn, anular
ondulado, plug, slug).

Evolução direta de `../estatistica2/extract_descriptors.py`. O MÉTODO do descritor é
IDÊNTICO (mesmo motor Schuck & Bodmann, mesmo whitening, mesma truncagem); a ÚNICA mudança
é o conjunto de padrões analisado: de {churn, wavy annular} para {churn, wavy annular,
plug, slug}. Isso realiza a Figura 17 (matriz de distância espectral) com os quatro padrões.

Rótulos: SOMENTE rótulos exatos (decisão do dono, 2026-10-06). Pastas com rótulo composto
(plug-slug, slug-churn, smooth annular, dispersed bubbles) e sintéticas são DESCARTADAS.
Resultado: churn=10, wavy annular=9, plug=4, slug=3 → 26 ensaios.

Descritor (FIXO, idêntico à estatistica2): perfil de escala do WPS auto de Schuck na coluna
de lag-zero b=0 (índice N-1 do WPS bilateral), truncado a 7001 px, whitened POR ESCALA
(resíduo da reta log|WPS| ~ log(escala), que remove o decaimento 1/f de fundo).

O motor `../v8_spatial_cropped_and_center/wavelet_v8.py` permanece INTOCADO; as funções do
núcleo são copiadas aqui (byte-a-byte equivalentes) com assert de identidade FFT≡laço em runtime.

Saída: descriptors_{variant}_{floor}.npz (variant ∈ cropped/center; floor ∈ piso50/piso10) com
  Z          (26, 200)  descritor whitened por run
  wps_b0     (26, 200)  |WPS b=0| cru
  scales     (200,)     escalas CWT
  labels     (26,)      padrão
  q, m       (26,)      fluxo de calor, vazão
  shift_mean (26,)      shift médio de costura (NaN no center)
  n_original (26,)      comprimento antes da truncagem
  runs       (26,)      nome da pasta
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

# eixo de frequência espacial — IDÊNTICO ao motor (wavelet_v8.py)
FREQ_MIN      = 1 / 3000
FREQ_MAX_BASE = 1 / 50                  # piso de análise do motor (λ >= 50 px)
WAVELET       = "cmor2.0-0.9549"

TRUNC_PX  = 7001                        # menor run de churn (m05q10_001) — herdado da estatistica2

# QUATRO padrões, rótulo EXATO. Composto/sintético fica de fora.
PATTERNS  = {"churn", "wavy annular", "plug", "slug"}

FLOORS = {"piso50": 1 / 50, "piso10": 1 / 10}   # piso10 = sonda de robustez (banda de costura)

# contagem esperada por padrão (trava de sanidade)
EXPECTED = {"churn": 10, "wavy annular": 9, "plug": 4, "slug": 3}


# ── Funções COPIADAS do motor (wavelet_v8.py) ────────────────────────────────

def _load_and_preprocess(path):
    """Cópia de wavelet_v8._load_and_preprocess: z-score único."""
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
    """Cópia de wavelet_v8.folder_flow_type."""
    folder = os.path.join(DATA_ROOT, name)
    for p in sorted(glob.glob(os.path.join(folder, "*.txt"))):
        stem = os.path.splitext(os.path.basename(p))[0].lstrip("_")
        if stem != "dados":
            return stem
    return name


def _biased_autocorr_loop(x, N):
    """Cópia EXATA do laço O(N^2) do motor (wavelet_v8.compute_wavelet)."""
    return np.array([np.sum(x[:N - m] * x[m:]) / N for m in range(N)])


def _biased_autocorr_fft(x, N):
    """Autocorrelação enviesada via FFT (≡ laço direto, ~4e-16). Acelera O(N^2) → O(N log N)."""
    full = fftconvolve(x, x[::-1], mode="full")   # comprimento 2N-1, pico em índice N-1
    return full[N - 1:] / N                        # lags 0..N-1


def _wps_from_signal(sig, use_fft=True, freq_max=FREQ_MAX_BASE):
    """Reproduz o núcleo autocorr+CWT do motor e devolve (|WPS b=0| (200,), scales (200,))."""
    N = len(sig)
    Rxx = _biased_autocorr_fft(sig, N) if use_fft else _biased_autocorr_loop(sig, N)
    Rxx_full = np.concatenate([Rxx[N - 1:0:-1], Rxx])
    assert Rxx_full.shape[0] == 2 * N - 1

    center_freq = pywt.central_frequency(WAVELET)
    freqs  = np.logspace(np.log10(FREQ_MIN), np.log10(freq_max), 200)
    scales = center_freq / freqs                    # sampling_rate = 1.0

    WPS, _ = pywt.cwt(Rxx_full.real, scales, WAVELET, sampling_period=1.0)
    wps_b0 = np.abs(WPS[:, N - 1])                   # coluna de lag-zero = b=0
    return wps_b0, scales


# ── Whitening por escala (fórmula fixa, idêntica à estatistica2) ─────────────

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


# ── Extração por variante × piso ──────────────────────────────────────────────

def extract_variant(runs, xlsx_name, freq_max, identity_checked):
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

        wps_b0, scales = _wps_from_signal(sig, use_fft=True, freq_max=freq_max)
        if scales_ref is None:
            scales_ref = scales

        # Assert de identidade: FFT-autocorr ≡ laço direto do motor (uma vez por piso).
        if not identity_checked.get(freq_max, False):
            wps_loop, _ = _wps_from_signal(sig, use_fft=False, freq_max=freq_max)
            err = np.max(np.abs(wps_b0 - wps_loop))
            assert err < 1e-9, "WPS b=0 (FFT) diverge do laço direto do motor: %g" % err
            print("  [assert] identidade FFT vs laço direto do motor OK (erro max=%.2e)" % err)
            identity_checked[freq_max] = True

        m, q = parse_m_q(r)
        Z_rows.append(whiten_per_scale(wps_b0, scales))
        b0_rows.append(wps_b0)
        labels.append(pattern); qs.append(q); ms.append(m)
        shifts.append(read_shift_mean(xlsx)); n_orig.append(n0); names.append(r)

    return dict(
        Z=np.array(Z_rows), wps_b0=np.array(b0_rows), scales=scales_ref,
        labels=np.array(labels), q=np.array(qs, float), m=np.array(ms, float),
        shift_mean=np.array(shifts, float), n_original=np.array(n_orig), runs=np.array(names),
    )


def _regression_anchor(d):
    """Confere que o subconjunto {churn, wavy annular} reproduz EXATAMENTE a baseline da
    estatistica2 (piso50/cropped). Garante que estender para 4 padrões não alterou o descritor
    dos 2 padrões originais — o acréscimo de plug/slug é aditivo, não perturba os existentes."""
    base = os.path.join(HERE, "..", "estatistica2", "descriptors_cropped_piso50.npz")
    if not os.path.exists(base):
        print("  [âncora] baseline estatistica2 ausente (%s) — pulado." % base)
        return
    b = np.load(base, allow_pickle=True)
    # alinhar por nome de run (estatistica2 só tem churn+annular)
    name_to_row = {n: i for i, n in enumerate(d["runs"])}
    maxerr = 0.0
    for i, n in enumerate(b["runs"]):
        j = name_to_row[str(n)]
        maxerr = max(maxerr, float(np.max(np.abs(d["Z"][j] - b["Z"][i]))))
    assert maxerr < 1e-9, "subconjunto 2-padrões diverge da baseline estatistica2: max|Δ|=%g" % maxerr
    print("  [âncora] subconjunto {churn, anular} ≡ estatistica2 piso50/cropped (max|Δ|=%.2e)" % maxerr)


def main():
    runs = list_runs()
    print("runs com ambos xlsx:", len(runs))
    identity_checked = {}

    for floor_name, freq_max in FLOORS.items():
        for variant, xlsx_name in [("cropped", "__dados_cropped_v11.xlsx"),
                                   ("center",  "__dados_center_v13.xlsx")]:
            print("\n=== %s / %s (freq_max=1/%g) ===" % (floor_name, variant, 1.0 / freq_max))
            d = extract_variant(runs, xlsx_name, freq_max, identity_checked)
            counts = {p: int((d["labels"] == p).sum()) for p in EXPECTED}
            assert counts == EXPECTED, "contagem inesperada: %s (esperado %s)" % (counts, EXPECTED)
            if floor_name == "piso50" and variant == "cropped":
                _regression_anchor(d)
            out = os.path.join(HERE, "descriptors_%s_%s.npz" % (variant, floor_name))
            np.savez(out, **d)
            print("salvo %s | %s | Z shape=%s"
                  % (os.path.basename(out), counts, d["Z"].shape))


if __name__ == "__main__":
    main()
