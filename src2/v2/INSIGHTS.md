# INSIGHTS — v2 Bilateral Wavelet Cross-Coherence Analysis

**Method:** Two-sided cross-correlation `Rxy_full` (length 2N−1, lags −(N−1)…0…+(N−1)) fed to a
single CWT. WCO normalized by WPS at lag zero (Schuck & Bodmann eq. 32.17). Biased estimator
(÷N). Wavelet: `cmor2.0-0.9549` (ω₀ ≈ 6.0). fs = 1000 Hz. Frequency range: 0.5–60 Hz.

x-axis convention: **lag < 0 → NAME2 leads NAME1 | lag > 0 → NAME1 leads NAME2**.

---

## 1. Synthetic (synthetic1 × synthetic2)

### Signals
Two isolated Morlet-like transient pulses embedded in broadband noise. synthetic1 pulse peaks
near sample 250; synthetic2 pulse peaks near sample 650 — a delay of ~400 samples (0.4 s).

### WPS
Both WPS concentrate power sharply at lag ≈ 0–250 samples, frequency band 10–30 Hz. Mean power
peaks at ~20 Hz. Power drops to noise floor beyond lag ~500 samples. The two WPS are nearly
identical in shape, confirming both signals carry the same transient structure.

### WCS
A single compact bright blob at lag ≈ +400 samples (positive side, right of the dashed line),
frequency band 10–25 Hz, peak cross-power ~0.030. The blob is localized and symmetric in
frequency, with no secondary peaks. Placement on the positive side confirms **synthetic1 leads
synthetic2 by ~400 ms**, consistent with the visual offset in the signals plot.

### WCO
A single high-coherence region (WCO ≈ 0.9–1.0, yellow) at lag ≈ +400 samples, 10–25 Hz. The
rest of the map is near zero. At frequencies below ~2 Hz a diffuse horizontal band of moderate
coherence (~0.2–0.3) spans the full lag range — this is a low-frequency artifact of the biased
estimator and should not be interpreted.

### Interpretation
The method correctly localizes the known delay: positive lag ≈ 400 samples = 0.4 s, frequency
~20 Hz. WCO peak near 1.0 validates the normalization. This pair serves as a positive control —
if the method fails here, it fails everywhere. **Method validated.**

---

## 2. Plug (m05q01_001 × m06q01_001)

### WPS
Both signals show dominant power at 1–3 Hz, with a bright peak near lag 0 decaying slowly over
~2000 samples. m05q01_001 has substantially higher peak power (~400) vs m06q01_001 (~175),
suggesting different signal amplitudes or plug passage rates between the two mass flux conditions.
Secondary structure is visible at 0.5–1 Hz. Mean power profiles peak at ~1.5–2 Hz for both.

### WCS
Cross-power is distributed in elongated horizontal bands at 1–3 Hz, spanning the full lag range
symmetrically (−6000 to +6000 samples). Multiple repetitive lobes appear at ~2000-sample intervals
(~2 s), consistent with the quasi-periodic passage of plug units. No single dominant lag: the
cross-power is roughly symmetric around zero, with slightly stronger lobes on the negative side
(lag ≈ −1000 to −2000 samples), tentatively indicating m06q01_001 leads slightly.

### WCO
WCO in the 0.5–3 Hz band reaches 0.6–0.8 (orange-yellow) across a wide lag range (±4000
samples). The coherence is notably **symmetric** around the zero-lag line, indicating neither run
consistently leads the other. Coherence drops sharply above 5 Hz. The 0.5 Hz band shows
persistent coherence ~0.4 across all lags, reflecting the long periodicity of plug flow.

### Interpretation
Plug flow has a clear quasi-periodic signature at 1–2 Hz (plug passage frequency). The symmetric
bilateral coherence suggests that across two different runs the periodicity is present at similar
rates, but no consistent phase lead exists — the runs are independent realizations of the same
regime, not causally linked. The method captures the characteristic frequency of plug flow well.

---

## 3. Slug (m04q01_001 × m06q05_001)

### WPS
Both signals concentrate power at 0.5–2 Hz, peak near lag 0. m04q01_001 has higher power (~120)
vs m06q05_001 (~100), similar in scale. Mean power peaks sharply at ~1 Hz. Compared to plug,
slug WPS has a sharper, narrower frequency peak — slug units are more regular in size/spacing.

### WCS
A dominant lobe centered near lag **0 to +1000 samples** (positive side), frequency 0.5–1.5 Hz,
peak cross-power ~80. This asymmetry is notable: cross-power is clearly stronger on the positive
side, suggesting **m04q01_001 leads m06q05_001**. Secondary lobes appear at +3000 and +5000
samples, confirming quasi-periodicity at ~2–3 s intervals. Negative-lag side is weaker and more
diffuse.

### WCO
WCO reaches near 1.0 (bright yellow) at lag 0 to +500 samples, 0.5–1.5 Hz — the strongest
localized coherence peak of all real flow patterns. A clear asymmetry: positive lags dominate,
reinforcing that m04q01_001 leads. Coherence remains above 0.4 up to lag +5000 samples at 0.5–1
Hz, indicating very regular slug periodicity. At 2–5 Hz, moderate coherence (~0.3) is scattered.

### Interpretation
Slug flow is the **most coherent and most periodic** of the four real patterns. The consistent
positive-lag dominance suggests the two runs are phase-aligned in a statistically meaningful way,
or that the mass flux difference between m04 and m06 creates a repeatable phase offset. The ~1 Hz
fundamental and ~2–3 s repetition period are characteristic of slug flow at this condition.

---

## 4. Churn (m04q20_001 × m05q20_001)

### WPS
Both signals show power broadly distributed from 0.5 to 10 Hz, with a concentration at 1–3 Hz.
The WPS is visually noisier than plug/slug — many small isolated power spots scattered across the
lag axis rather than a clean decaying ridge. This reflects churn's chaotic, irregular nature.
m04q20_001 has significantly higher peak power (~100) vs m05q20_001 (~20). Mean power profiles
peak at 1.5–2 Hz.

### WCS
Cross-power lobes at 1–3 Hz, roughly symmetric around zero, with peak lobes at lag ≈ ±1000–2000
samples. The bilateral symmetry here is cleaner than plug — lobes appear at approximately equal
strength on both sides, consistent with no dominant lead direction. Peak cross-power ~35 (much
lower than slug/plug in absolute terms). Sporadic elevated cross-power at 3–5 Hz adds to the
disordered appearance.

### WCO
Coherence reaches 0.6–0.8 in the 1–3 Hz band at lags ±1000–2000 samples. Notably, the 0.5 Hz
band shows **continuous high coherence (WCO ≈ 1.0, bright yellow)** across the entire lag range
— this is a very-low-frequency artifact (period >> signal length) that should be disregarded.
Above 3 Hz, coherence is fragmented and low (~0.2–0.3). The coherence map has a diffuse,
patchy appearance consistent with irregular flow.

### Interpretation
Churn flow lacks a clean periodic structure — the WCO map is fragmented where slug/plug show
crisp bands. The moderate coherence at 1–3 Hz indicates some residual periodicity (oscillating
column), but it is much weaker and less regular than slug or plug. The bilateral symmetry
confirms no consistent phase lead between runs. Churn is the most difficult pattern to
characterize with this method due to its chaotic nature.

---

## 5. Wavy Annular (m04q80_001 × m05q80_001)

### WPS
Both signals show a sharp, narrow power ridge at **~3 Hz** extending across a very long lag axis
(N ≈ 18000 samples). The 3 Hz peak is the most frequency-specific of all patterns — the wavy
annular regime has a well-defined wave passage frequency. Mean power peaks sharply at 3 Hz with
secondary structure at 1.5 Hz. The WPS ridge persists to large lags (>10000 samples), indicating
strong stationarity of the 3 Hz oscillation.

### WCS
Cross-power concentrates in a compact region near **lag 0, frequency 2–4 Hz**, with peak ~25.
Unlike plug/slug, the lobes do not repeat at large lags — cross-power decays rapidly beyond ±5000
samples. There is a slight positive-lag bias (m04q80_001 leads marginally), but the peak is close
to zero lag, suggesting the two runs are nearly phase-synchronized at 3 Hz. Moderate cross-power
patches appear at 1–2 Hz.

### WCO
WCO peaks at ~0.8 near lag 0–1000 samples, 2–4 Hz. The coherence map is the most spatially
compact of all patterns — concentrated in a small region rather than spread across many lags.
This reflects the high stationarity and regularity of wavy annular flow. Coherence falls off
faster with increasing lag compared to plug/slug, possibly because the wave amplitude varies more
between independent runs at different mass flux. Below 1.5 Hz, coherence is low and diffuse.

### Interpretation
Wavy annular flow has the most spectrally pure signature: a single dominant frequency at ~3 Hz.
The compact WCO peak near zero lag suggests the two runs are well-synchronized in frequency even
across different conditions (m04 vs m05). This is the pattern most likely to be cleanly
separable from the others by frequency content alone. The 3 Hz peak is a strong candidate as a
discriminating feature for classification.

---

## Cross-Pattern Comparison

| Pattern      | Dominant freq (WPS) | WCO peak (approx) | Lag symmetry       | Periodicity quality |
|-------------|--------------------|--------------------|--------------------|--------------------|
| Synthetic   | ~20 Hz             | ~1.0               | Strongly asymmetric (+400 s) | N/A (test signal) |
| Plug        | 1–2 Hz             | 0.6–0.8            | Symmetric          | Moderate, multi-lobe |
| Slug        | ~1 Hz              | ~1.0               | Asymmetric (+ side) | High, sharp        |
| Churn       | 1–3 Hz             | 0.4–0.6            | Symmetric          | Low, fragmented    |
| Wavy annular| ~3 Hz              | ~0.8               | Near-symmetric     | High, compact      |

### Key findings

1. **Method validation:** The synthetic pair correctly localizes the 400 ms delay at 20 Hz with
   WCO ≈ 1.0, confirming the bilateral two-sided CWT of `Rxy_full` works as intended.

2. **Frequency separation:** Wavy annular (~3 Hz) and slug/plug (~1 Hz) occupy distinct frequency
   bands. Churn overlaps with both. This supports frequency-based classification.

3. **Slug is the most coherent real pattern.** Its WCO near 1.0 at zero-to-moderate positive
   lags indicates high inter-run regularity and a measurable phase lead between the two mass flux
   conditions. This is physically plausible: slug units form periodically and the formation rate
   is condition-dependent.

4. **Churn is the least coherent.** Fragmented WCO, lower peak values, no consistent lag
   structure — consistent with the chaotic, oscillating-column nature of churn flow. This also
   makes it the hardest to classify and the most likely source of confusion with neighboring
   patterns in a ML classifier.

5. **Lag asymmetry as a feature.** Slug shows clear positive-lag dominance; plug is symmetric;
   wavy annular peaks near zero. The lag position of the WCO peak (positive, zero, or symmetric)
   is a potential discriminating feature beyond frequency alone.

6. **Large-lag artifact warning.** All patterns show elevated coherence at the very lowest
   frequencies (~0.5 Hz) across the full lag range. This is a biased-estimator artifact (few
   independent cycles fit in the signal length at 0.5 Hz) and should not be interpreted as
   physical coherence.
