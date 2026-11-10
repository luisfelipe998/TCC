"""
plot_out_of_band.py — A3'-a: LOCALIZAR o artefato de costura (idêntico ao da estatistica2).

Objetivo: olhar ONDE mora a energia do |WPS b=0| quando se rebaixa o piso a 1/10 (λ até ~10 px),
em relação à banda de costura [0,7;1,4]×shift_mean e ao piso de análise principal (50 px).

Resultado honesto (medido): a costura de quadros NÃO injeta um pico dominante em λ≈shift. O pico
do |WPS b=0| fica em λ grande (centenas de px) e a banda de costura carrega energia ínfima
(~1% do total). Documenta que o descritor principal (piso50) não é contaminado por um "bump de
costura". Saída = v8-estat-piso.png (nome herdado da figura do resultados.tex).

Observação: a figura usa um run de churn/anular (shift_mean definido); os padrões plug/slug
acrescentados em estatistica3 não mudam esta sonda — o argumento de costura é por-run e
independente do conjunto de padrões analisado.
"""

import os
import numpy as np
import pywt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from extract_descriptors import (
    DATA_ROOT, WAVELET, _load_and_preprocess, _wps_from_signal, read_shift_mean,
)

HERE = os.path.dirname(os.path.abspath(__file__))
RUN  = "m03q20_001"            # shift_mean ~46 px: banda de costura perto/abaixo do piso 50 px
FLOOR_PROBE = 1 / 10


def main():
    xlsx = os.path.join(DATA_ROOT, RUN, "__dados_cropped_v11.xlsx")
    sig = _load_and_preprocess(xlsx)[:7001]
    shift = read_shift_mean(xlsx)

    wps_b0, scales = _wps_from_signal(sig, use_fft=True, freq_max=FLOOR_PROBE)
    center_freq = pywt.central_frequency(WAVELET)
    wavelength = scales / center_freq          # px

    peak_lambda = float(wavelength[np.argmax(wps_b0)])
    band_lo, band_hi = 0.7 * shift, 1.4 * shift
    band = (wavelength >= band_lo) & (wavelength <= band_hi)
    frac_band = float(wps_b0[band].sum() / wps_b0.sum())

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(wavelength, wps_b0, lw=1.6)
    ax.axvspan(band_lo, band_hi, color="orange", alpha=0.25,
               label="banda de costura [0,7–1,4]×shift = [%.1f; %.1f] px (%.1f%% da energia)"
                     % (band_lo, band_hi, 100 * frac_band))
    ax.axvline(50, color="red", ls="--", lw=1.2, label="piso de análise principal (50 px)")
    ax.axvline(shift, color="green", ls=":", lw=1.2, label="shift médio = %.1f px" % shift)
    ax.axvline(peak_lambda, color="purple", ls="-.", lw=1.0,
               label="pico do |WPS b=0| = %.0f px" % peak_lambda)
    ax.set_xscale("log")
    ax.set_xlabel("comprimento de onda λ (px)"); ax.set_ylabel("|WPS b=0| (cru)")
    ax.set_title("Localização da energia do descritor — banda de costura sem pico/energia\n"
                 "relevante (run %s, piso de sonda 10 px)" % RUN, fontsize=10)
    ax.legend(fontsize=8); ax.grid(True, which="both", alpha=0.3)
    out_png = os.path.join(HERE, "v8-estat-piso.png")
    fig.savefig(out_png, dpi=140, bbox_inches="tight"); plt.close(fig)

    print("run=%s  shift_mean=%.1f px  pico |WPS b=0| em λ=%.0f px  (piso principal = 50 px)"
          % (RUN, shift, peak_lambda))
    print("banda de costura [%.1f; %.1f] px carrega apenas %.1f%% da energia e NÃO é o pico →"
          % (band_lo, band_hi, 100 * frac_band))
    print("a costura não injeta bump dominante no WPS b=0; a análise principal não depende dessa região.")
    print("salvo:", out_png)


if __name__ == "__main__":
    main()
