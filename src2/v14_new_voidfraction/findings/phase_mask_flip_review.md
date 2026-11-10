# Flips of the gas-silhouette candidate

The current scope contains 26 retained acquisitions and 25 statistical acquisitions passing the saturation criterion across 5 patterns. Mixed labels, smooth annular and the excluded replicate are omitted before computing masks.

These measurements use the same inferred phase mask as the panorama void 
reconstruction. They are spatial changes in a gas-silhouette candidate, 
whose gas labeling remains conditional. The user confirmed rows 18:71 as the full tube interior.

Pixel flip rate = mask changes over all 53 rows / [53 × (W−1)]. 
Fixed-row flip rate = changes in mask row y=40 / (W−1). Each complete isolated 
gas segment along that row contributes two interfaces; entry count and flip count 
are saved separately. Neither quantity is a temporal event frequency.

The legacy center_* column names mean fixed image row y=40; this row has not been validated as the physical tube center.

| Pattern | n | Candidate pixel flip rate | Row y=40 flip rate | Row y=40 gas duty |
|---|---:|---:|---:|---:|
| dispersed bubbles | 1 | 0.002464 | 0.004416 | 0.0136 |
| plug | 3 | 0.009730 | 0.010645 | 0.7052 |
| slug | 3 | 0.010419 | 0.010018 | 0.7985 |
| churn | 10 | 0.011190 | 0.005160 | 0.9603 |
| wavy annular | 8 | 0.001133 | 0.000228 | 0.9995 |

The original v11.3 optical flip rates use raw dark/bright classes and 
40 rows; these candidate-mask values use reconstructed silhouettes and 
53 rows. They measure different masks and must retain separate names.

The same reconstruction that makes the annular core continuous also 
suppresses internal brightness fluctuations. This is an intended phase 
definition, but can suppress real liquid fragments inside churn gas 
envelopes. It is not independent validation of a desired ranking.

Thresholds 0.60–0.80 are all reported in the sensitivity CSV. Group tests 
are exploratory comparisons with no multiple-comparison correction. 
Mann–Whitney uses the exact method when the smaller group has at most eight acquisitions and there are no pooled ties; otherwise it uses the asymptotic method with tie and continuity corrections. Each saved comparison records the actual method, tie flag and both sample sizes. The rank test concerns distributions, not medians alone, and does not validate absolute phase calibration.

Retained supplied labels are unchanged. The single retained dispersed-bubble acquisition does not establish within-pattern variability.

Arithmetic checks reproduce two interfaces for a complete isolated 
segment and zero for a constant gas record, and verify complement/reversal 
invariance for all 26 reconstructed masks.

Reproduce: `python3 src2/v14_new_voidfraction/scripts/phase_mask_flip_analysis.py`.
