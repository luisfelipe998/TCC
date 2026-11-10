# v5 Spatial Wavelet Analysis — Key Findings v2

Spatial domain (1 px/sample), wavelet `cmor2.0-0.9549` (ω₀≈6). WCO computed from bilateral
autocorrelations `Rxx_full` / `Ryy_full` (length 2N-1) and bilateral cross-correlation `Rxy_full`.

---

## 1. Preprocessing Chain

Steps applied to `Voltage_2` (optical sensor column) before wavelet analysis.

| Step | Operation | Rationale |
|------|-----------|-----------|
| 1. Load | Read `Voltage_2` from `__dados_cropped_v11.xlsx` | Cropped rows 20–59 of panorama (flow core only, excludes tube walls) |
| 2. Invert | `sig = 255 − sig` | Liquid phase → high values, gas phase → low values (matches physical interpretation) |
| 3. Linear detrend | `scipy.signal.detrend(sig, type="linear")` | Removes illumination gradient along panorama (camera/light inhomogeneity); verified negligible residual (<4.5% of signal range) |
| 4. Zero-mean | `sig = sig − mean(sig)` | Ensures DC=0 after detrend (detrend leaves near-zero mean; this is a no-op in practice) |
| 5. Z-score | `sig = sig / std(sig)` | Normalises amplitude across runs; makes correlation coefficients unitless and comparable |
| 6. Biased autocorrelation | `Rxx[m] = (1/N) Σ x[n]x[n+m]`, `m=0..N-1` | Biased estimator (divides by N, not N-m); suppresses noise at large lags; Rxx[0]=1 after z-score |
| 7. Bilateral symmetrisation | `Rxx_full = [Rxx[N-1..1] ∣ Rxx[0..N-1]]` (length 2N-1) | Provides full symmetric context to the CWT; b=0 column at center index N-1 sees equal lag history on both sides — required by Schuck & Bodmann eq. 32.14/32.17 |
| 8. CWT | `pywt.cwt(Rxx_full, scales, "cmor2.0-0.9549")` | Complex Morlet ω₀≈6; scales mapped to 50–3000 px log-spaced; sampling period = 1.0 px |
| 9. WCO | `WCO = \|WCS\|² / (\|WPS1[:,N-1]\| × \|WPS2[:,N-1]\|)` at b=0 | Normalised wavelet coherence (Schuck & Bodmann eq. 32.17); denominator evaluated at lag=0 center column |

---

## 2. Key Findings

### Dominant spatial scales per flow pattern

The 200–500 px vs 500–1000 px power ratio is the strongest cross-pattern discriminator.

| Pattern | Dominant band | Characteristic |
|---------|--------------|----------------|
| Wavy annular | 200–500 px (40%) | Periodic liquid film wave spacing; consistent ridge across all 8 runs |
| Plug | 200–500 px (33%) + 500–1000 px (31%) | Two-band signature: bubble cap diameter + body length |
| Slug / Churn | 500–1000 px (36–39%) | Elongated bubble bodies; slug tighter and more localized than churn |
| Plug-slug | 1000–2000 px (52%) | Very large transitional structures; single run |
| Slug-churn | 500–1000 px (27%), heavy tail at 1k–3k px | Transition regime; less spatially organized than pure slug or churn |
| Dispersed bubbles | Spread 50–2400 px, slight peak 500–1000 px | No dominant scale; many small randomly-spaced bubbles |

### WCO same-pattern coherence

Same-pattern WCO at lag=0 is consistently high (~0.7–1.0) within each pattern's dominant band,
confirming that the dominant spatial scale is a stable fingerprint across independent runs of
the same operating condition. Cross-pattern WCO is lower and less organized.

### Wavy annular

WPS shows a broad, diffuse ridge spanning 100–500 px, consistent across all 8 runs regardless
of mass/heat flux. Cross-run WCS shows recurring hot-spots at 300–500 px spaced ~3000–5000 px
apart in lag — the annular wave pattern is spatially quasi-periodic and repeatable. WCO at lag=0
is moderate (~0.5–0.7), meaning the scale fingerprint is consistent but not perfectly locked
between independent runs.

### Churn

WPS has a concentrated locus near 500–700 px that shifts between runs (different bubble
coalescence regimes across q″=5–80 kW/m²). WCO at lag=0 is high and broad (~0.8–1.0 in the
300–1000 px band) — despite run-to-run peak shifts, churn always concentrates power in this
region, making it a reliable signature.

### Slug

WPS is tighter and more localized than churn, with a sharp peak near 500–700 px. The three
non-outlier runs agree well (WCO ~0.8–1.0, narrow at lag=0). Two outlier runs (m05q05_002,
m06q05_001) show strongly shifted distributions, suggesting slug body length is sensitive to
operating conditions.

### Plug

Two visible WPS ridges: 200–500 px (bubble cap) and 500–1000 px (bubble body). WCO at lag=0
spans both bands (~0.8–1.0) and extends to ±2000 px lag — plug bubbles are large enough that
spatial autocorrelation reaches far. The dual-band signature distinguishes plug from slug
(500–1000 px dominant) and wavy annular (200–500 px dominant).

### Plug-slug

Single run m03q01_001, unusually large structures concentrated above 1000 px. No cross-run
analysis possible.

### Slug-churn

Two runs (m04q05, m07q10) are quite different: m04q05 peaks at ~1000 px while m07q10 has a
broader ridge at 200–500 px. WCO between them is moderate (~0.5–0.6 at lag=0), spread across
a wide band — consistent with a transition regime where spatial structure is less organized.

### Dispersed bubbles

Single run m07q01_001. Power spread from 50 px to the COI limit (~2400 px) with no dominant
ridge — consistent with many small, randomly spaced bubbles without a characteristic spatial
scale. COI limits the 2k–3k px band.

---

## 3. Caveats

| Run | Issue | Treatment |
|-----|-------|-----------|
| `m03q80_001` (smooth annular) | Panorama reconstruction failed — mean shift clamped to 1 px, N=501, λ_max=378 px | Excluded from all analysis |
| `m05q05_002`, `m06q05_001` (slug) | Anomalous power distributions; likely unusual flow conditions or recording artefacts | Excluded from pattern means |
| `dispersed bubbles` | Single run; COI limit at ~2400 px (2k–3k band indicative only) | No intra-pattern coherence analysis possible |
| `plug-slug` | Single run | No intra-pattern coherence analysis possible |
| All patterns | Biased autocorrelation estimator applies implicit Bartlett taper (suppresses large-lag estimates) | Consistent across all runs; no correction needed for relative comparisons |
| COI formula | CWT input is bilateral Rxx_full (length 2N-1); T&C (1998) eq.19 limit is λ_max=(2N-1)/(2√2·cf) | All patterns except dispersed bubbles: λ_max >> 3000 px → full axis reliable |
