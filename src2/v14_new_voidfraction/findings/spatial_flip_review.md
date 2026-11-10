# Optical flip-rate sensitivity: saved panorama pixel boundaries

The current scope contains 26 retained acquisitions and 25 statistical acquisitions passing the saturation criterion. Mixed labels, smooth annular and the user-excluded replicate are omitted before processing.

This analysis uses only the saved `src2/data_processed/*/__panorama_v11.png` images. `flip_rate` counts changes in a fixed binary optical threshold mask between neighboring axial pixels, averaged across rows. It is a spatial boundary density and is invariant when the dark/bright labels are swapped. Its value does not certify that those pixels are gas or liquid.

`F = number_of_pixel_switches / [H*(W-1)]`. The quantity `1000*F` is mean pixel-state changes per 1000 adjacent axial comparisons per row; the numerator counts all rows, so it is not a count of bubbles or liquid slugs. The v11.3 baseline uses rows 20:60 and a fixed, per-panorama Otsu threshold. There is no Schmitt trigger and no native-time sampling rate in this calculation.

| Pattern | n | median F | changes/1000 axial comparisons per row |
|---|---:|---:|---:|
| dispersed bubbles | 1 | 0.006735 | 6.735 |
| plug | 3 | 0.017680 | 17.680 |
| slug | 3 | 0.017372 | 17.372 |
| churn | 10 | 0.029995 | 29.995 |
| wavy annular | 8 | 0.017681 | 17.681 |

Mann–Whitney p-values below retain the two-sided asymptotic method with continuity correction and tie correction where ties occur. Saved comparisons record the method, sample sizes and tie flag. These exploratory rank tests compare distributions, not medians alone; no multiple-comparison correction is applied.

Churn versus wavy annular: AUC=0.8625; exploratory two-sided run-level Mann–Whitney p=0.01133, n=10 and 8. The churn>annular median direction holds in 24/24 ROI/threshold/filter variants.

The variants use three radial ROIs (10:70, 20:60, 25:55), Otsu, Otsu±10 gray levels, fixed thresholds 70/90/110, the image median, and Otsu after a 3×3 median filter. None was selected to enforce a regime ordering. Sensitivity values and actual pixel switch counts are recorded for all 26 retained acquisitions. Saturation exclusion is the existing >20% threshold on baseline core pixels and is kept constant across variants.

The saved panorama is stitched from native-frame strips with unequal measured widths and feathering at seams. Spatial density can therefore depend on advection speed, stitching and image noise as well as interface geometry. The supplemental raw-frame audit shows that churn has a higher median native optical flip activity than annular, but substantial run-level overlap; the native five-station comparison has two-sided p=0.3508. A claim that churn is universally more temporally intermittent does not follow from the stronger spatial result.

The retained dispersed-bubble group contains 1 acquisition. Its baseline spatial optical flip rate is 0.006735; one acquisition does not characterize within-pattern variability. The group total is the current unsaturated retained scope, not the earlier unfiltered count.

No absolute gas volume fraction, phase velocity or liquid-reference calibration is inferred by this script. The computed spatial flip values are reproducible optical observables. `flip_share_intensity_change_le5` diagnoses small-amplitude boundary crossings but does not identify which boundaries are noise.

Run `python3 src2/v14_new_voidfraction/scripts/spatial_flip_sensitivity.py` to regenerate `flip_rate_spatial_sensitivity.csv`, `flip_rate_spatial_group_sensitivity.csv` and this review. Native time checks and exact all-run v11.3 reproduction are documented in `flip_rate_review.md`.
