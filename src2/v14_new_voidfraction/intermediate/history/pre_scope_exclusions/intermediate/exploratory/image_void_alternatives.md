# Image-only void alternatives: numerical audit

All numerical estimators use raw frames from `src2/data_processed`. Regime labels are attached after computation. No estimator is tuned to a desired regime ranking. 120 uniformly spaced frames and every eighth axial column are sampled per sequence (501 original frames). These correlated observations do not constitute independent experimental replicates.

**Result:** replacing the median with a mean is necessary for a time/space average, but neither statistic calibrates brightness into gas volume. Dark-interface envelopes and closed-contour filling recover additional gas-outline information. They remain optical proxies because detected shadow positions are refracted interface images; a geometric diameter cannot yet be assigned to those shadows.

## Per-run results

BC values are center bright-width / 60 pixels. Envelope values below are zero-filled *conditional circular reconstructions*, E[detected × (apparent interface separation/60)²]. Zero means failed detection or no detected center gas, not a proven all-liquid section. They are not calibrated void fractions.

| Run | Label | BC mean | BC median | Interface detected | Envelope surrogate | Filled dark projected area | Template detected |
|---|---|---:|---:|---:|---:|---:|---:|
| m03q01_001 | plug-slug | 0.383 | 0.383 | 0.680 | 0.335 | 0.612 | 0.754 |
| m03q05_001 | churn | 0.271 | 0.333 | 0.725 | 0.353 | 0.800 | 0.867 |
| m03q10_001 | wavy annular | 0.286 | 0.350 | 0.452 | 0.245 | 0.811 | 0.718 |
| m03q20_001 | wavy annular | 0.285 | 0.367 | 0.232 | 0.134 | 0.740 | 0.566 |
| m03q40_001 | wavy annular | 0.475 | 0.500 | 0.115 | 0.070 | 0.495 | 0.476 |
| m03q80_001 | smooth annular | 0.488 | 0.500 | 0.024 | 0.015 | 0.434 | 0.135 |
| m04q01_001 | slug | 0.311 | 0.317 | 0.735 | 0.333 | 0.656 | 0.885 |
| m04q05_001 | slug-churn | 0.268 | 0.300 | 0.821 | 0.330 | 0.695 | 0.929 |
| m04q10_001 | churn | 0.258 | 0.317 | 0.831 | 0.359 | 0.799 | 0.922 |
| m04q20_001 | churn | 0.255 | 0.317 | 0.835 | 0.371 | 0.779 | 0.939 |
| m04q40_001 | wavy annular | 0.265 | 0.333 | 0.588 | 0.307 | 0.855 | 0.835 |
| m04q80_001 | wavy annular | 0.291 | 0.350 | 0.163 | 0.097 | 0.838 | 0.562 |
| m05q01_001 | plug | 0.399 | 0.333 | 0.600 | 0.277 | 0.644 | 0.736 |
| m05q05_001 | slug | 0.370 | 0.317 | 0.651 | 0.241 | 0.628 | 0.737 |
| m05q05_002 | dispersed bubbles | 0.363 | 0.283 | 0.543 | 0.172 | 0.606 | 0.641 |
| m05q10_001 | churn | 0.231 | 0.283 | 0.840 | 0.331 | 0.800 | 0.971 |
| m05q20_001 | churn | 0.238 | 0.283 | 0.830 | 0.342 | 0.806 | 0.983 |
| m05q40_001 | churn | 0.250 | 0.300 | 0.704 | 0.337 | 0.872 | 0.993 |
| m05q80_001 | wavy annular | 0.280 | 0.350 | 0.269 | 0.152 | 0.878 | 0.991 |
| m06q01_001 | plug | 0.310 | 0.300 | 0.693 | 0.276 | 0.666 | 0.761 |
| m06q05_001 | slug | 0.253 | 0.283 | 0.783 | 0.294 | 0.665 | 0.846 |
| m06q10_001 | churn | 0.210 | 0.267 | 0.814 | 0.320 | 0.720 | 0.943 |
| m06q20_001 | churn | 0.207 | 0.250 | 0.724 | 0.321 | 0.839 | 0.978 |
| m06q40_001 | wavy annular | 0.233 | 0.300 | 0.491 | 0.256 | 0.878 | 0.974 |
| m06q80_001 | wavy annular | 0.243 | 0.333 | 0.286 | 0.158 | 0.885 | 0.955 |
| m07q01_001 | dispersed bubbles | 0.677 | 0.833 | 0.002 | 0.000 | 0.387 | 0.092 |
| m07q05_001 | plug | 0.335 | 0.300 | 0.646 | 0.222 | 0.618 | 0.703 |
| m07q10_001 | slug-churn | 0.249 | 0.283 | 0.725 | 0.278 | 0.673 | 0.808 |
| m07q20_001 | churn | 0.196 | 0.233 | 0.842 | 0.304 | 0.740 | 0.936 |
| m07q40_001 | churn | 0.196 | 0.233 | 0.775 | 0.331 | 0.840 | 0.973 |
| m07q80_001 | wavy annular | 0.225 | 0.300 | 0.459 | 0.241 | 0.878 | 0.968 |

## What each alternative resolves

1. **Mean instead of median.** Physical averaged void is an arithmetic expectation of instantaneous gas area. v11.3 returns the median bright width across frames and columns, which can ignore the less common phase. This audit records both statistics and five consecutive frame-block means. The bright-width numerator still does not equal gas diameter.
2. **Dark-interface envelope.** Smooth radial intensity with σ=1 pixel; locate upper/lower dark extrema within rows 19:38 and 42:70. Accept interior extrema when both are ≥15 gray levels below the center. Repeat at 10 and 20 levels. All gates are identical across runs. This uses the boundary, including a bright interior, rather than declaring every shadow pixel to be gas. An absence of accepted interfaces is an unresolved observation.
The CSV also records an expanded lower search, rows 42:73, rejecting row 72. It exposes substantial search-window dependence when gas shadows overlap the lower wall/reflection band. Better detection coverage does not establish that the optical shadow is the physical gas boundary.
3. **Circular reconstruction.** If an undistorted axial view identifies the physical connected gas-core diameter d_g, and the core is concentric and axisymmetric, α(t,z)=(d_g/D)². Averaging the square matters: E[(d_g/D)²] is generally different from E[d_g/D]². Absolute px→mm scale cancels from this ratio. A second camera is not mathematically mandatory under that assumption. The present bright band and dark shadows have no validated mapping to d_g, however; in churn, off-axis fragments further violate this reconstruction. Apparent diameter candidates 56, 60 and 64 pixels quantify sensitivity, not an independently measured tube diameter.
4. **Filled dark contours.** Otsu-segment rows 18:71, close 3×3 gaps, then fill enclosed holes. This captures interiors of some dark-outlined bubbles. Annular interfaces may connect to image edges and fail to enclose any region, while merging shadows can falsely enclose liquid. The projected image area is not the circular cross-sectional area.
5. **Liquid-like template.** Build a robust normalized radial profile from `m07q01_001`, normalize each column using its upper eight background rows, and detect paired deficits below 65% of the template. This explicitly exposes the reference dependence. `m07q01_001` contains bubbles and is not a validated liquid-only calibration photo; cross-run exposure, gamma and reflections can differ. The resulting template detector is a feasibility probe.

## Negative control and identifiability

The liquid-like `m07q01_001` sequence has visible dispersed bubbles but a nearly featureless bright center. Its bright-core midpoint median can be high because that operation measures brightness shape, not phase occupancy. A per-run adaptive threshold necessarily finds image classes even without a separable gas/liquid signal. Likewise an own-run temporal percentile is not automatically a liquid baseline: a persistent annular core can occupy the same pixel throughout the record.

The 80-row acquisition ROI includes the wall/reflection bands. The windows 20:60 or 10:70 are analysis crops, not demonstrated physical tube diameters. Even an apparent wall spacing is affected by cylindrical refraction. A proper inverse model requires the tube's optical geometry (including glass outer diameter and refractive index), camera/illumination geometry, fluid refractive index, and validation against liquid-only and known-volume states or another independent void measurement.

## Source check

Schepperle et al. (2024), Sensors 24, 3363, DOI 10.3390/s24113363, uses a **rectangular** channel (depth 0.5 mm, width 1.5 mm; §2.1), contour extraction, and histogram comparison against a **separate liquid-only image of the same subsection acquired before boiling** (§2.2, pp. 6–7). Its Eq. (3) counts pixels *after phase labeling*. It does not validate replacing contour/phase classification with uncalibrated bright-pixel counting or an own-run 20th percentile in a circular tube.

**Defensible next step:** preserve these descriptors as optical diagnostics; recover run-specific inlet enthalpy/quality and obtain an independent thermodynamic void prediction. For direct image void, annotate gas outlines using a real liquid calibration image and validate the circular/refraction model. Correcting a rank order is not a calibration.
