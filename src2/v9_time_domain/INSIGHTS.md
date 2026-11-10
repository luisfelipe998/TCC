# v9 Time-Domain Wavelet Analysis — Insights

This document summarizes the conclusions of the v9 time-domain wavelet
cross-spectrum / coherence study. The signal is the per-frame flow-core mean
brightness (rows 20–59 of each raw camera frame collapsed to one scalar per
frame), analyzed with the same biased-correlogram → CWT → WPS/WCS/WCO
machinery as the spatial-domain v8 pipeline. All 91 within-pattern pairs across
five flow patterns (churn, plug, slug, slug-churn, wavy annular) were computed
with zero failures.

## Headline

**Across all five flow patterns, the v9 time-domain wavelet coherence shows no
genuine shared cross-run signature.** FFT dominant frequencies scatter with the
operating point (mass/heat flux) rather than clustering by pattern, and the
biased-correlogram WCO sits at its saturating null floor with no frequency
localization. Only slug shows a marginal ~2-sigma elevation, and that rests on
n = 3 and is not convincing.

This is consistent with the entire v5–v9 lineage, where the same
biased-correlogram estimator repeatedly returned a negative verdict.

## What the signal actually measures

The clips are very short: T = 0.07–0.13 s (N = 287–501 samples after resampling
to a common 4000 fps). The record length pins the cone-of-influence (COI) floor
at 10.6 / 13.2 / 18.5 Hz depending on N. As a result, the **true physical
pattern frequencies fall below the resolvable floor**:

| Pattern | Physical frequency | Resolvable? |
|---|---|---|
| Plug | < 1 Hz | No (far below floor) |
| Slug | ~2.8 Hz | No |
| Churn | 15–28 Hz | Mostly no |
| Annular | 2–20 Hz | Mostly no |

So v9 does **not** measure the slow structure periodicity of each flow regime.
It measures the ~20–200 Hz **intra-structure texture** band that survives above
the COI floor. The pattern-discriminating information that lives in the slow
structure is invisible to this short-clip time-domain analysis.

## Dominant frequency tracks flux, not pattern

Per-run FFT dominant frequencies scatter and follow the operating point:

- **Churn** (10 runs): 41.8–153.3 Hz, a loose cluster near ~80 Hz in 3 of 10
  runs but no shared line.
- **Plug** (4 runs): 41.8–97.6 Hz, no common peak.
- **Slug** (3 runs): the lowest and most consistent — 27.9 / 27.9 / 55.7 Hz,
  all top-3 lists in the 27.9–97.6 Hz band, and the highest 10–30 Hz energy
  fraction of the dataset (consistent with slow intermittent slug flow).
- **Wavy annular** (9 runs): dominant frequency rises **monotonically with
  flux**, from ~16 Hz (m03 low-q) to 195 Hz (m07q80). "Wavy annular" as labeled
  spans a flux-driven frequency range, not one characteristic tone.

Because dominant frequency is a flux coordinate, any frequency-based separator
would largely re-learn mass/heat flux rather than the flow pattern.

## Per-pattern coherence

Mean zero-lag WCO in the resolvable band, averaged across within-pattern pairs:

| Pattern | Runs / Pairs | Mean zero-lag WCO | Verdict |
|---|---|---|---|
| Churn | 10 / 45 | 0.294 | no signature |
| Wavy annular | 9 / 36 | 0.217 | no signature |
| Plug | 4 / 6 | 0.432 | no signature (COI-edge artifact) |
| Slug | 3 / 3 | 0.403 | weak signature (not convincing) |
| Slug-churn | 2 / 1 | 0.460 | no signature (n = 1 pair) |

- **Churn** — 0.294 at zero lag, 0.137 over all lags: near the estimator floor,
  no frequency-localized elevated sub-band.
- **Wavy annular** — 0.217, the lowest of all patterns: broadband, no frequency
  lock. The single highest pair (m05q80 × m06q40 = 0.44) is isolated and modest.
- **Plug** — 0.432 looks elevated but is a COI-edge artifact: per-pair peaks all
  glue to ~19–23 Hz (the 18.5 Hz resolvable-band edge), several peak values
  exceed 1.0 (max 1.28) — impossible for true normalized coherence, confirming
  estimator saturation. Peak location does not coincide with any FFT dominant
  frequency (42–98 Hz). All-lag mean is only 0.129.
- **Slug** — 0.403 (per-pair 0.396 / 0.437 / 0.378) vs a same-N white-noise
  floor of 0.269 ± 0.069, i.e. only ~2 sigma above baseline; the all-lag mean
  (0.141) sits at/below the noise floor (0.172 ± 0.051). Backed by partial FFT
  agreement (two runs lock at 27.9 Hz, third at 55.7 Hz with 27.9 Hz as 2nd
  peak) and a shared low-band tendency. Real but modest, and n = 3 is too small
  to trust.
- **Slug-churn** — single pair, 0.460, but the peak lands at 18.7 Hz (COI-floor
  edge) and the two runs disagree on dominant frequency and even in character
  (one slug-like, one churn-like). n = 1 pair gives no statistical basis.

## Estimator caveat (read before quoting any WCO number)

The WCO estimator is degenerate; its absolute magnitude must **not** be read as
physical coherence. In `wavelet_v9.py` line 233:

```python
WCO = |WCS|**2 / (|WPS1[:, N-1]| * |WPS2[:, N-1]|)
```

the cross-spectrum WCS is taken at all lags but normalized only by the
**zero-lag** (N-1) column of the two auto-spectra. This biased-correlogram ratio
saturates toward and beyond 1 at zero lag (peak WCO > 1.0, up to 1.28 for plug —
impossible for true coherence) and decays away from zero lag. That is exactly
why every pattern shows a large zero-lag / all-lag gap (e.g. plug 0.432 vs
0.129; slug-churn 0.460 vs 0.130). Consequently "high" zero-lag WCO is **not**
evidence of shared structure — it is the estimator's null behavior.

The only defensible readings are (a) relative comparison against a same-N
white-noise floor (done for slug), and (b) whether elevation is
frequency-localized away from the COI edge (it never is — peaks glue to the
~19–23 Hz floor edge).

## Resample-to-4000 fps note

The dataset was filmed at 4000 / 5000 / 7000 fps, so cross-pattern pairs mixed
sampling rates on one axis. Every run is anti-alias downsampled (polyphase FIR,
`scipy.signal.resample_poly`) to the common 4000 fps floor. This was a necessary
bookkeeping step so pairs share one frequency axis and one COI-floor logic — but
it did **not** create physical signal. Temporal resolution is still bounded by
each clip's short duration, which is what pins the COI floor and keeps the slow
pattern frequencies unresolvable. Downsampling discards no analyzed signal
(< 1.2% of power sits above 1000 Hz, and the axis caps at 1000 Hz anyway),
keeps every sample real (no interpolated invention), and preserves clip duration
T — so the per-run COI floor, set by T not by fps, is unchanged. Its only
material effect is enabling clean cross-run pairing (91/91 pairs, 0 failures)
and making the COI floor depend on record length (N = 401 → 13.24 Hz; N = 287 →
18.51 Hz) rather than native fps. The verdict did not change after resampling.

## Limitations

- Clips are far too short (T = 0.07–0.13 s, N = 287–501) to resolve the true
  physical pattern frequencies; most sit below the 10.6 / 13.2 / 18.5 Hz COI
  floor.
- The biased-correlogram WCO estimator saturates toward/above 1 and normalizes
  cross-lag WCS by only the zero-lag auto-spectra, so magnitudes are not
  interpretable as true coherence and > 1.0 peaks occur.
- Reported coherence peaks consistently land on the COI-floor edge (~19–23 Hz),
  a boundary/estimator artifact rather than physical structure.
- The signal is only per-frame mean brightness — a single scalar per frame —
  discarding all spatial structure the v8 pipeline captures.
- Small samples for several patterns: slug n = 3 runs / 3 pairs, slug-churn
  n = 2 runs / 1 pair; no statistical basis for a per-pattern signature there.
- Dominant frequency confounds with operating point (mass/heat flux), so
  frequency-based separation would partly re-learn flux rather than flow
  pattern.
- Only one channel (Voltage_2 brightness proxy) is analyzed; detrend + z-score
  then truncation to common min length discards low-frequency content and
  shortens records further.
- Resampling to a common 4000 fps unifies dt but does not add real temporal
  resolution beyond each clip's native acquisition.

## Recommendation

Treat v9 time-domain wavelet coherence as a **complement, not a classifier**.
Honest verdict for the thesis: it is weakly useful at most and not sufficient
alone — four of five patterns sit at the biased-WCO null floor, the lone
weak-signature (slug, ~27.9 Hz low-frequency tendency, ~2 sigma over a
white-noise floor) rests on n = 3 and does not generalize, and "high" coherence
numbers are estimator saturation / COI-edge artifacts.

The genuinely defensible findings are: (1) the texture band lives at ~20–200 Hz,
and (2) dominant frequency tracks mass/heat flux (wavy annular's monotonic
16 → 195 Hz rise is the cleanest example) — both flux coordinates, not pattern
labels. Because the clips are too short to resolve the slow slug / plug /
annular periodicity, the pattern-discriminating information lives in the
**spatial domain (v8)**, which captures the slow structure. Report v9 as
evidence that time-domain texture and coherence do not separate patterns on
their own, and lean on the spatial pipeline for classification. Do not present
the 0.40–0.46 zero-lag WCO values as coherence without the estimator-saturation
and COI-edge caveats.
