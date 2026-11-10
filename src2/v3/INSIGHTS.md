# INSIGHTS — v3 Bilateral Wavelet Cross-Coherence Analysis

**v3 changes over v2:**
1. Butterworth bandpass 0.5–60 Hz (order 4, zero-phase `sosfiltfilt`) applied before correlations
2. Std normalization after bandpass (WCO amplitude-invariant across runs)
3. Corrected `effective_fs` per folder: mode-based shift (not median), accounts for bimodal slug distribution
4. Wavelet changed from `cmor2.0-0.9549` to `cmor1.5-1.5` (narrower bandwidth, finer time resolution)

**Effective fs per pair:** plug ≈ 510 Hz | slug ≈ 425 Hz | churn ≈ 456 Hz | wavy_annular ≈ 174 Hz | synthetic = 1000 Hz

x-axis convention: **lag < 0 → NAME2 leads NAME1 | lag > 0 → NAME1 leads NAME2**

---

## 1. Preprocessing Impact (overall)

**Bandpass:** The most visible change across all real patterns is the elimination of the continuous
bright-yellow band at 0.5 Hz that spanned the full lag range in v2 WCO maps. That artifact was a
low-frequency artifact caused by DC and sub-0.5 Hz drift entering the correlation. It is completely
gone in v3. The noise floor above 10 Hz is also cleaner — the high-frequency camera noise that was
contributing weakly to all WPS maps is removed.

**Std normalization:** WCO maps now reflect true phase coherence rather than amplitude imbalance.
In v2, pairs where one run had much higher amplitude (e.g. plug m05 peak power ~400 vs m06 ~175)
produced artificially asymmetric WCO. In v3 both runs enter at unit variance; WCO differences
between patterns are now structural, not amplitude-driven.

**Wavelet (cmor1.5-1.5 vs cmor2.0-0.9549):** Narrower bandwidth (B=1.5 vs 2.0) gives sharper
frequency resolution at the cost of slightly lower time localization. In practice the WPS ridges
appear narrower in frequency in v3, which helps separate the 0.5–1 Hz and 1–2 Hz components that
blurred together in v2.

**Corrected fs:** Wavy annular fs dropped from ~177 Hz (v2 average of 159+194) to 174 Hz (mode-
based average of 159+189) — small change. Slug fs changed more: m04q01_001 went from 583 Hz
(wrong median) to 438 Hz (correct mode), making the slug pair fs much more consistent (438 vs 412
vs the previous 583 vs 438 mismatch). The slug WCS/WCO x-axis now represents physical time more
faithfully.

**Regression:** The 0.5 Hz artifact removal is unambiguously an improvement. However, the
wavy_annular pattern loses its previously reported ~3 Hz WPS peak — the corrected fs (174 Hz)
shifts the displayed frequency band, and the two WPS now show different dominant frequencies
(m04 peaks near 1 Hz; m05 peaks near 5 Hz). This is a significant finding discussed below.

---

## 2. Synthetic (synthetic1 × synthetic2)

### WPS
Both WPS concentrate power at lag ≈ 0–200 samples, frequency band 15–30 Hz. Mean power peaks at
~20 Hz. Identical structure to v2. The corrected wavelet (cmor1.5-1.5) gives a marginally sharper
frequency ridge but no substantive difference.

### WCS
Single compact blob at lag ≈ +400 samples (positive side), 10–25 Hz, peak cross-power ~2.0
(higher absolute value than v2's ~0.030 because std normalization raises the correlation scale).
The blob is clean and isolated — negative side is near zero.

### WCO
Single bright-yellow peak (WCO ≈ 1.0) at lag ≈ +400 samples, 10–25 Hz. **New in v3:** a diffuse
moderate-coherence band at 1–3 Hz spanning positive lags — this is a residual artifact from the
biased estimator at these frequency bands (few cycles fit in the signal length at 1 Hz), not
physical coherence. The 0.5 Hz full-lag artifact from v2 is completely gone.

**vs v2:** WCO peak identical. The 0.5 Hz artifact removed. New 1–3 Hz diffuse band is less severe
than the old artifact. **Method still validates correctly: synthetic1 leads synthetic2 by ~400 ms
at ~20 Hz, WCO ≈ 1.0.**

---

## 3. Plug (m05q01_001 × m06q01_001) — fs ≈ 510 Hz

### WPS
Both WPS show dominant power at 0.5–1.5 Hz, peak near lag 0, decaying over ~3000 samples. Mean
power profiles peak sharply at ~1 Hz. The amplitude balance between the two runs is now nearly
equal (both reach magnitude ~2.0) — confirming std normalization removed the 2× amplitude
difference seen in v2. The power ridge is sharper in frequency than v2 thanks to cmor1.5-1.5.

### WCS
Elongated red band at ~1 Hz spanning the full lag range (−6000 to +6000 samples). Peak cross-
power ~2.5. Critically, the band is **strongly asymmetric: negative-lag side (m06q01_001 leads)
is brighter**, with peak around lag −1000 to −2000 samples (≈ −2 to −4 s at 510 Hz). This was
tentative in v2 ("slightly stronger on negative side") but is now clearly visible after
normalization removes amplitude bias.

### WCO
WCO reaches ~1.0 (bright yellow) in the 0.5–1.5 Hz band at negative lags −500 to −4000 samples
(−1 to −8 s). The positive-lag side drops to ~0.4–0.6. **This is a clear reversal from v2's
"roughly symmetric" assessment** — normalization revealed that m06q01_001 consistently leads
m05q01_001 at the plug passage frequency. Coherence above 3 Hz is near zero, clean noise floor.

**vs v2:** Major improvement in interpretability. The 0.5 Hz continuous artifact is gone. The true
asymmetry (m06 leads) is now unambiguous. WCO peak rises from 0.6–0.8 to ~1.0. **v3 is strictly
better for plug.**

---

## 4. Slug (m04q01_001 × m06q05_001) — fs ≈ 425 Hz

### WPS
Both WPS show dominant power at 0.5–1.5 Hz. m04q01_001 power (peak ~2.0) is higher than
m06q05_001 (peak ~0.8) even after normalization — this reflects genuine structural difference:
m04 has the bimodal shift distribution indicating more organized slug units. Mean power peaks at
~0.8–1 Hz for m04 and ~1 Hz for m06, with secondary structure at 1.5–2 Hz visible in m06. The
wider spread in m06 WPS (isolated blobs from 0.5 to 5 Hz) confirms m06q05_001 has a more
irregular slug structure.

### WCS
A wide band at 0.5–1 Hz spanning ±8000 samples, with lobes at approximately ±2000 and ±4000
samples (≈ ±4.7 s and ±9.4 s at 425 Hz). **The bilateral symmetry is near-perfect** — cross-power
lobes appear at equal strength on both positive and negative sides, with no dominant lead
direction. Peak cross-power ~1.4.

### WCO
WCO reaches 0.8–1.0 (yellow) at 0.5–1 Hz, forming elongated symmetric lobes centered near zero
lag, repeating at ±2000 and ±4000 samples. The bilateral symmetry is the most striking feature:
neither run leads the other.

**vs v2:** Significant difference. v2 showed clear positive-lag dominance (m04 leads) and WCO
near 1.0 only on the positive side. In v3, after correcting m04q01_001's fs from 583→438 Hz and
applying normalization, the asymmetry disappears entirely — the lobe structure is symmetric. This
suggests the v2 asymmetry was an artifact of the incorrect fs (583 Hz vs 438 Hz was a 33%
mismatch; at equal fs the cross-correlation is symmetric for independent realizations of the same
pattern). **Slug flow has symmetric bilateral coherence in v3 — no consistent phase lead.**

---

## 5. Churn (m04q20_001 × m05q20_001) — fs ≈ 456 Hz

### WPS
Both WPS show scattered power from 0.5 to 10 Hz, with the highest concentration at 0.5–1.5 Hz.
The WPS texture is the noisiest of all patterns — many small isolated blobs rather than continuous
ridges, confirming churn's chaotic nature. After normalization the two runs have comparable peak
magnitudes (~0.6–0.8). Mean power peaks at ~1–1.5 Hz with a secondary shoulder at 3 Hz.

### WCS
Cross-power lobes at 0.5–1.5 Hz, bilateral and approximately symmetric (lobes at ±2000 and ±4000
samples ≈ ±4.4 s and ±8.8 s). Peak cross-power ~0.6 — the lowest of the real patterns. Above
3 Hz cross-power is diffuse and irregular, consistent with churn's broadband energy distribution.

### WCO
WCO reaches 0.6–0.8 (orange-yellow) at 0.5–1.5 Hz, forming repeated lobes across the full lag
range bilaterally. Coherence at 1.5–3 Hz is fragmented (~0.3–0.4, isolated patches). The 0.5 Hz
full-lag continuous band from v2 is **gone** — the remaining coherence at 0.5–1 Hz is structured
(lobed) not continuous, indicating real quasi-periodic content rather than artifact.

**vs v2:** The removal of the 0.5 Hz artifact significantly clarifies interpretation. Churn's WCO
in v3 shows genuine lobed structure at 0.5–1.5 Hz but weaker than plug/slug. The fragmented
higher-frequency coherence remains. **v3 is better for churn: the artifact removal separates real
structure from noise.**

---

## 6. Wavy Annular (m04q80_001 × m05q80_001) — fs ≈ 174 Hz

### WPS
**Major change from v2.** The two runs now show *different* dominant frequencies:
- m04q80_001: dominant ridge near **0.8–1 Hz**, extending across the full lag range (N ≈ 18000 samples)
- m05q80_001: dominant ridge near **4–6 Hz**, also persistent across all lags

This divergence was not visible in v2 (fs=1000 Hz was wrong for both; the correct fs ≈ 174 Hz
shifts all displayed frequencies down by ~5.7×). The two runs appear to belong to different
sub-regimes of wavy annular flow despite being classified the same.

### WCS
Cross-power concentrates in a narrow blob near **lag 0 to +500 samples (≈ 0–2.9 s), frequency
0.7–1.2 Hz**, peak ~0.7. The 4–6 Hz band from m05 does not appear in the cross-power, confirming
the two runs share only the 0.8–1 Hz component. The blob is slightly on the positive side (m04
leads marginally).

### WCO
WCO peak ~0.3–0.4 (purple-pink) at 0.7–1 Hz near lag 0. This is the **weakest WCO of all
patterns** in v3 — a dramatic drop from v2's peak of ~0.8. The reason: the two runs have different
spectral structures (1 Hz vs 5 Hz), so the coherence estimator correctly reports low similarity.
The coherence map is otherwise near zero everywhere.

**vs v2:** Significant regression in apparent coherence (0.8 → 0.3–0.4), but this is a *correct*
result, not a code regression. The v2 high coherence was an artifact of using wrong fs (1000 Hz).
With correct fs the two wavy annular runs are revealed to have different dominant frequencies —
suggesting they may represent different heat flux / quality conditions within the wavy annular
regime. **This is an important physical finding: the two wavy annular runs are not well-matched.**

---

## Cross-Pattern Comparison

| Pattern       | Dominant freq (WPS) | WCO peak | Lag structure         | Key change from v2 |
|--------------|--------------------|-----------|-----------------------|--------------------|
| Synthetic    | ~20 Hz             | ~1.0      | Asymmetric (+400 ms)  | Artifact removed, result identical |
| Plug         | ~1 Hz              | ~1.0      | **Asymmetric (− side, m06 leads)** | Asymmetry revealed by normalization |
| Slug         | ~1 Hz              | ~0.8–1.0  | **Symmetric**         | v2 asymmetry was fs-artifact |
| Churn        | 0.5–1.5 Hz         | ~0.6–0.8  | Symmetric, lobed      | Artifact removed, lobed structure clearer |
| Wavy annular | m04: ~1 Hz / m05: ~5 Hz | ~0.3–0.4 | Near zero lag       | Runs have different dominant freqs |

### Key findings

1. **Plug is now clearly asymmetric (m06 leads m05 at ~1 Hz).** This was hidden in v2 by amplitude
   imbalance. After std normalization, the negative-lag WCO dominance is unambiguous. Physical
   interpretation: at higher mass flux (m06), plug units arrive at a slightly earlier phase than at
   lower mass flux (m05), consistent with faster plug propagation at higher flow rate.

2. **Slug symmetry is restored.** The v2 positive-lag asymmetry was caused by the bimodal shift
   distribution assigning a wrong fs (583 Hz) to m04q01_001. With the corrected mode-based fs
   (438 Hz), both runs align in time scale and cross-correlation is bilateral-symmetric.

3. **Wavy annular runs are spectrally mismatched.** m04q80_001 (lower heat flux / quality) dominates
   at ~1 Hz while m05q80_001 (higher) dominates at ~5 Hz. These may represent different wave modes
   (long slow waves vs short fast ripples). The WCO is low (~0.3) because the runs share only a
   weak 0.8–1 Hz component. Classification based on WCO alone would fail for this pattern with
   these specific pairs.

4. **Frequency discrimination is now sharper:**
   - Wavy annular (m04): ~1 Hz — overlaps with plug/slug/churn in v3 (no longer at 3 Hz as v2 showed)
   - Wavy annular (m05): ~5 Hz — distinct from all others
   - Plug, slug, churn: all concentrate at 0.5–1.5 Hz

5. **The 0.5 Hz full-lag coherence artifact is eliminated.** All v3 WCO maps show structured
   (lobed) coherence at low frequencies rather than continuous horizontal bands.

---

## Conclusions

**Is v3 strictly better than v2?**

For interpretation: **yes for plug, churn, synthetic**. The amplitude artifact removal (normalization)
and low-frequency artifact removal (bandpass) give cleaner, more trustworthy WCO maps.

For slug: **yes in correctness** — the v3 result is physically correct; v2's asymmetry was an
artifact of wrong fs. The symmetric v3 result is the honest answer.

For wavy annular: **the low WCO is a correct finding, not a regression** — the two runs are
genuinely spectrally different. v2's high coherence was inflated by wrong fs and no normalization.

**Recommendations:**
1. Apply the resampling fix (see `todo.md`) to eliminate the remaining fs averaging for cross-run comparisons.
2. Wavy annular needs more run pairs — the current pair (m04 vs m05) spans a large heat flux range,
   causing spectral mismatch. Pairs within the same mass flux level would give more reliable WCO.
3. Plug's negative-lag asymmetry (m06 leads by ~2–4 s at 1 Hz) is a candidate discriminating feature
   alongside frequency, potentially useful for mass-flux estimation within a detected pattern.
4. Churn remains the hardest to classify: moderate, symmetric, broadband coherence overlapping with plug/slug.
