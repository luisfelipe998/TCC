# Flip-rate audit: spatial record versus native camera time

All calculations use `src2/data_processed` only. A flip is a change in an optical threshold state. These states have not been independently validated as gas and liquid.

## Exact definitions

The v11.3 quantity is `sum(abs(diff(binary_mask, axis=1))) / [40*(panorama_width-1)]`. It is a spatial, per-pixel adjacent-column change fraction. It does not use a Schmitt trigger, and its values near 0.02–0.04 are not counts per 1000 px. All 26 retained panoramas reproduce the existing v11.3 function exactly.

The native-time alternative uses all ordered acquisition frames, camera column x=508 and, separately, x={254,381,508,635,762}, rows 20:60. A fixed Otsu threshold is calculated from all raw core pixels within a run. `r_t = flips/[40*C*(N-1)]`; `r_t*fps` is optical threshold switches per pixel per second, not bubble-passage frequency. No image registration or panorama stitching enters this estimate.

Camera `_config.xsv` files provide `Rate`, `Frames` and `Step`; filenames are checked against the recorded step and frame count. All runs contain 501 consecutive frames. Rates are 4000, 5000 or 7000 fps, and elapsed first-to-last-frame duration is 0.125, 0.100 or 0.0714286 s. The zero-based time interval has N−1 comparisons.

## Run-level group medians

The current user scope retains 26 acquisitions; the existing >20% panorama-core saturation rule leaves 25 statistical acquisitions: dispersed bubbles 1, plug 3, slug 3, churn 10, wavy annular 8. Mixed labels, smooth annular and the excluded replicate do not enter any per-run computation or grouped statistic. Supplied retained labels are unchanged; a single bubble run does not establish within-pattern variability.

| Pattern | n | spatial rate | native x508 rate | native 5-station rate | native optical switches/pixel/s |
|---|---:|---:|---:|---:|---:|
| dispersed bubbles | 1 | 0.006735 | 0.027300 | 0.026010 | 182.070 |
| plug | 3 | 0.017680 | 0.126450 | 0.131760 | 922.320 |
| slug | 3 | 0.017372 | 0.113650 | 0.123610 | 865.270 |
| churn | 10 | 0.029995 | 0.170975 | 0.175255 | 1226.785 |
| wavy annular | 8 | 0.017681 | 0.149350 | 0.155065 | 939.155 |

## Interpretation and sensitivity

Mann–Whitney p-values below retain the two-sided asymptotic method with continuity correction and tie correction where ties occur. Each comparison saves its method, sample sizes and tie flag. This rank test compares distributions; it is not generally a test of medians alone.

- `spatial_flip_rate`: churn median 0.0299945, annular 0.0176814; AUC(churn>annular)=0.863, exploratory two-sided run-level Mann–Whitney p=0.01133 (n=10 versus 8).
- `raw_x508_flip_rate`: churn median 0.170975, annular 0.14935; AUC(churn>annular)=0.688, exploratory two-sided run-level Mann–Whitney p=0.1976 (n=10 versus 8).
- `raw_x508_flip_hz_per_pixel`: churn median 1196.82, annular 910.05; AUC(churn>annular)=0.688, exploratory two-sided run-level Mann–Whitney p=0.1976 (n=10 versus 8).
- `raw_5probe_flip_rate`: churn median 0.175255, annular 0.155065; AUC(churn>annular)=0.637, exploratory two-sided run-level Mann–Whitney p=0.3508 (n=10 versus 8).
- `raw_5probe_flip_hz_per_pixel`: churn median 1226.79, annular 939.155; AUC(churn>annular)=0.637, exploratory two-sided run-level Mann–Whitney p=0.3508 (n=10 versus 8).
- `raw_spatial_median3_flip_rate`: churn median 0.172775, annular 0.153255; AUC(churn>annular)=0.637, exploratory two-sided run-level Mann–Whitney p=0.3508 (n=10 versus 8).
- `raw_rate_over_shuffle`: churn median 0.665583, annular 0.709394; AUC(churn>annular)=0.287, exploratory two-sided run-level Mann–Whitney p=0.1426 (n=10 versus 8).

The churn>annular median direction in native optical switches/pixel/s survives 21/21 fixed ROI/threshold variants (ROIs 10:70, 20:60, 25:55; Otsu, Otsu±10, fixed 70/90/110, median). This sweep evaluates optical activity, not true phase intermittency. Group p-values are exploratory and not corrected for the many alternatives.

Spatial versus native-time-Hz run-level rank correlation: Spearman rho=0.887, p=3.539e-09. Different flow advection speeds and different frame rates make their numerical values non-interchangeable.

`raw_flip_share_intensity_change_le5` reports the share of switches caused by changes of at most five gray levels. Such switches may include noise, JPEG effects and weak optical fluctuations; the value diagnoses threshold fragility but does not prove that any individual switch is noise. A 3×3 spatial median filter is an explicit nuisance-control variant.

Time shuffling does not make flip-rate disappear. It usually increases changes by destroying temporal correlation. The correct pixel-specific shuffle null is the mean of `2*p_pixel*(1-p_pixel)*N/(N-1)`; the old global-alpha Bernoulli null confounds a static radial profile with temporal occupancy. A complement of the binary mask preserves flip-rate exactly, so flipping gas/bright polarity cannot fix this descriptor.

The observation window is short. Native pixel activity is directly measurable, but long-term stationarity and low-frequency passage statistics are not established. Camera exposure, built-in filtering, intensity clipping and JPEG compression remain measurement limits; no phase-velocity estimate or assumed Taylor conversion is used here.

## Reproduce

`python3 src2/v14_new_voidfraction/scripts/raw_flip_analysis.py` regenerates the per-run values, ROI/threshold sensitivity, group statistics and diagnostic plots. Invariants cover constant and alternating records, complement and time-reversal invariance, N−1 normalization, pixel-specific shuffle null, actual filename/config consistency, and exact all-run reproduction of v11.3.

Outputs: `flip_rate_per_run.csv`, `flip_rate_by_pattern.csv`, `flip_rate_sensitivity.csv`, `flip_rate_group_sensitivity.csv`, `flip_rate_audit.json`, `flip_rate_spatial_vs_native.png`, `flip_rate_native_kymographs.png`.
