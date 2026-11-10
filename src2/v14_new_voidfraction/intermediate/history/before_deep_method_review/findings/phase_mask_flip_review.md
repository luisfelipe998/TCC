# Flips of the gas-silhouette candidate

The current scope contains 26 retained acquisitions and 25 unsaturated statistical acquisitions across 5 patterns. Mixed labels, smooth annular and the excluded replicate are omitted before computing masks.

These measurements use the same inferred phase mask as the panorama void 
reconstruction. They are spatial changes in a gas-silhouette candidate, 
whose gas labeling and radial crop remain conditional.

Pixel flip rate = mask changes over all 53 rows / [53 × (W−1)]. 
Center flip rate = changes in mask row 40 / (W−1). Each complete isolated 
center gas segment contributes two interfaces; entry count and flip count 
are saved separately. Neither quantity is a temporal event frequency.

| Pattern | n | Candidate pixel flip rate | Center flip rate | Center gas duty |
|---|---:|---:|---:|---:|
| dispersed bubbles | 1 | 0.003202 | 0.004416 | 0.0136 |
| plug | 3 | 0.010045 | 0.011825 | 0.7084 |
| slug | 3 | 0.009750 | 0.011302 | 0.8015 |
| churn | 10 | 0.009980 | 0.005295 | 0.9606 |
| wavy annular | 8 | 0.000963 | 0.000228 | 0.9995 |

The original v11.3 optical flip rates use raw dark/bright classes and 
40 rows; these candidate-mask values use reconstructed silhouettes and 
53 rows. They measure different masks and must retain separate names.

The same reconstruction that makes the annular core continuous also 
suppresses internal brightness fluctuations. This is an intended phase 
definition, but can suppress real liquid fragments inside churn gas 
envelopes. It is not independent validation of a desired ranking.

Thresholds 0.60–0.80 are all reported in the sensitivity CSV. Group tests 
are exploratory comparisons, not evidence of absolute phase calibration. 
Retained supplied labels are unchanged. The single retained dispersed-bubble acquisition does not establish within-pattern variability.

Arithmetic checks reproduce two interfaces for a complete isolated 
segment and zero for a constant gas record, and verify complement/reversal 
invariance for all 26 reconstructed masks.

Reproduce: `python3 src2/v14_new_voidfraction/scripts/phase_mask_flip_analysis.py`.
