# Panorama shadow-envelope reconstruction

**A usable image solution is possible as a conditional geometric estimate.** This method identifies dark gas outlines against a liquid-like radial template, fills the gas interior, and includes separately enclosed off-center bubbles. It does not equate dark pixels alone, or bright pixels alone, with gas. The result is a physically meaningful reconstruction formula applied to an optical silhouette; its absolute calibration remains conditional.

## Reference inferred from available panoramas

The three plug panoramas `m05q01_001`, `m06q01_001`, `m07q05_001` contain visible liquid intervals. In each, normalize radial profiles by the 75th percentile of rows 30:50; select the 5% of columns with the smallest radial standard deviation in rows 20:68. Median selected profiles give the liquid-like template. Selection is based on profile structure, with no flow-label fitting. `m07q01_001` is excluded from reference construction and is tested independently. The JSON saves exact candidate column coordinates.

The inferred reference has wall/reflection minima at image rows **15 and 75**, an apparent separation of **60 px**. We use D_app=60 px, and repeat with 56 and 64 px. These are apparent optical dimensions: the minima have not been validated as physical lumen edges.

## Segmentation and physical reconstruction

Smooth normalized radial profiles by σ=1 px and divide by the inferred template. Locate paired upper/lower shadow deficits below ratio 0.70, expand each contiguous deficit to its outer shadow edge, and fill between them. Reject paired fills when the top and bottom belong to different closed-shadow components and their accepted axial interval is shorter than 20 px: this prevents bridging liquid between unrelated bubble outlines. Long open interfaces remain eligible, and separate closed gas outlines remain eligible regardless of the paired-fill rule. Add closed dark-outlined components away from radial crop edges, after 3×3 closing and hole filling. Identical parameters apply to every panorama.

The projected gas-silhouette fraction is A_g,projected/(D_app L_record). The **conditional circular void estimate** is E_z[Σ_i(d_i(z)/D_app)²], where d_i is each separate silhouette's radial width. This follows from A_g(z)=Σ_i π d_i²/4 and A_tube=πD²/4. It requires circular/axisymmetric gas cross-sections and negligible overlap in viewing depth. The arithmetic mean of squares is used, not the median width or square of its mean. The px→mm factor cancels in the ratio.

## Per-run estimates and sensitivity

These numbers are **conditional optical reconstructions**, not independently calibrated measurements. Saturated runs are retained and explicitly marked unsupported. The threshold interval is a perturbation range, not a statistical confidence interval.

| Run | Pattern | Conditional α | Missing-width interval | Threshold 0.60–0.80 range | Reference variant range | Center gas duty (spatial) | Status |
|---|---|---:|---:|---:|---:|---:|---|
| m03q01_001 | plug-slug | 0.559 | 0.559–0.727 | 0.545–0.568 | 0.559–0.560 | 0.767 | unsupported_saturation |
| m03q05_001 | churn | 0.755 | 0.755–0.973 | 0.746–0.761 | 0.755–0.755 | 0.999 | conditional_optical_geometry |
| m03q10_001 | wavy annular | 0.777 | 0.777–0.997 | 0.775–0.778 | 0.777–0.777 | 1.000 | conditional_optical_geometry |
| m03q20_001 | wavy annular | 0.780 | 0.780–1.000 | 0.780–0.780 | 0.780–0.780 | 1.000 | conditional_optical_geometry |
| m03q40_001 | wavy annular | 0.765 | 0.765–1.000 | 0.753–0.774 | 0.764–0.766 | 1.000 | unsupported_saturation |
| m03q80_001 | smooth annular | 0.780 | 0.780–1.000 | 0.780–0.780 | 0.780–0.780 | 1.000 | unsupported_saturation |
| m04q01_001 | slug | 0.551 | 0.551–0.672 | 0.528–0.569 | 0.551–0.552 | 0.831 | conditional_optical_geometry |
| m04q05_001 | slug-churn | 0.549 | 0.549–0.684 | 0.521–0.573 | 0.548–0.550 | 0.881 | conditional_optical_geometry |
| m04q10_001 | churn | 0.693 | 0.693–0.888 | 0.677–0.707 | 0.693–0.694 | 0.989 | conditional_optical_geometry |
| m04q20_001 | churn | 0.706 | 0.706–0.926 | 0.690–0.719 | 0.705–0.707 | 0.996 | conditional_optical_geometry |
| m04q40_001 | wavy annular | 0.766 | 0.766–0.989 | 0.762–0.769 | 0.766–0.766 | 1.000 | conditional_optical_geometry |
| m04q80_001 | wavy annular | 0.780 | 0.780–1.000 | 0.779–0.780 | 0.780–0.780 | 1.000 | conditional_optical_geometry |
| m05q01_001 | plug | 0.498 | 0.498–0.622 | 0.483–0.511 | 0.498–0.500 | 0.708 | conditional_optical_geometry |
| m05q05_001 | slug | 0.410 | 0.410–0.493 | 0.392–0.430 | 0.409–0.411 | 0.688 | conditional_optical_geometry |
| m05q05_002 | dispersed bubbles | 0.308 | 0.308–0.344 | 0.290–0.328 | 0.308–0.309 | 0.627 | conditional_optical_geometry |
| m05q10_001 | churn | 0.597 | 0.597–0.711 | 0.570–0.627 | 0.597–0.599 | 0.947 | conditional_optical_geometry |
| m05q20_001 | churn | 0.640 | 0.640–0.772 | 0.604–0.666 | 0.639–0.642 | 0.973 | conditional_optical_geometry |
| m05q40_001 | churn | 0.728 | 0.728–0.907 | 0.700–0.744 | 0.727–0.729 | 0.993 | conditional_optical_geometry |
| m05q80_001 | wavy annular | 0.769 | 0.769–0.976 | 0.748–0.775 | 0.769–0.770 | 0.999 | conditional_optical_geometry |
| m06q01_001 | plug | 0.466 | 0.466–0.572 | 0.443–0.491 | 0.466–0.467 | 0.748 | conditional_optical_geometry |
| m06q05_001 | slug | 0.476 | 0.476–0.568 | 0.447–0.509 | 0.475–0.477 | 0.801 | conditional_optical_geometry |
| m06q10_001 | churn | 0.542 | 0.542–0.673 | 0.508–0.585 | 0.542–0.544 | 0.880 | conditional_optical_geometry |
| m06q20_001 | churn | 0.650 | 0.650–0.841 | 0.608–0.686 | 0.648–0.651 | 0.935 | conditional_optical_geometry |
| m06q40_001 | wavy annular | 0.720 | 0.720–0.921 | 0.672–0.746 | 0.718–0.723 | 0.976 | conditional_optical_geometry |
| m06q80_001 | wavy annular | 0.753 | 0.753–0.968 | 0.717–0.765 | 0.752–0.754 | 0.993 | conditional_optical_geometry |
| m07q01_001 | dispersed bubbles | 0.002 | 0.002–0.002 | 0.001–0.003 | 0.002–0.002 | 0.014 | conditional_optical_geometry |
| m07q05_001 | plug | 0.373 | 0.373–0.417 | 0.351–0.400 | 0.373–0.375 | 0.667 | conditional_optical_geometry |
| m07q10_001 | slug-churn | 0.466 | 0.466–0.573 | 0.439–0.501 | 0.466–0.467 | 0.761 | conditional_optical_geometry |
| m07q20_001 | churn | 0.516 | 0.516–0.601 | 0.477–0.557 | 0.514–0.517 | 0.881 | conditional_optical_geometry |
| m07q40_001 | churn | 0.636 | 0.636–0.796 | 0.588–0.673 | 0.636–0.638 | 0.948 | conditional_optical_geometry |
| m07q80_001 | wavy annular | 0.724 | 0.724–0.927 | 0.671–0.747 | 0.720–0.726 | 0.975 | conditional_optical_geometry |

## Adversarial checks

The sensitivity CSV repeats thresholds 0.60/0.65/0.70/0.75/0.80, leaves out each of the three reference panoramas in turn, changes reference candidate selection to 1% and 10%, narrows the radial crop, and changes the open-interface minimum interval from 20 to 10/30 px. The negative-control dispersed panorama is never a training reference. Five spatial-block means expose within-record variability. These checks measure internal robustness; they cannot substitute for an independent void calibration.

Held-out liquid-like candidate checks are saved separately: reference construction excludes the tested panorama. Candidate projected foreground is 0.017% for m05q01, 0.158% for m06q01, and 0.027% for m07q05. Their conditional area contributions are all below 0.0001. These are inferred liquid candidates rather than certified all-liquid ground truth; the independent dispersed control tests the same detector outside reference construction.

The overlay figure makes false phase assignments reviewable. The independent control should show its small bubbles with little foreground between them. Plug/slug images should show both shadow boundaries and bright interiors, with liquid slugs empty. Annular should keep the connected core, while churn retains uncertain fragmented outlines.

**A censoring limit matters for annular:** filling all 53 rows in the interior crop imposes a maximum visible estimate (53/60)²=0.7803. Several annular runs approach this exact ceiling. Their apparent agreement or low threshold scatter must not be mistaken for precise void determination. The CSV records the fraction of columns touching either radial crop boundary and a missing-width interval: α_visible ≤ α_if_full_width ≤ E[α_visible+1_clipped(1−α_visible)]. This interval assumes correct phase labels, circular reconstruction, and an undistorted diameter; it is not an unconditional physical confidence interval.

## Remaining physical assumptions

A single corrected side view can reconstruct an axisymmetric gas body; a second view and absolute mm scale are not invariably necessary. The present cylindrical tube's refraction can move shadow edges away from physical interface positions. The glass/illumination inverse map remains unknown, and irregular churn or overlapping dispersed bubbles can violate circular reconstruction. Those are the principal limits on calling the numerical result the true void fraction.

A stitched panorama supplies a spatial record. Its mean is spatial occupancy under the panorama reconstruction's frozen-advection assumption; it is not automatically a time average at a fixed point. A physical time average requires native frame timestamps and fixed-column phase masks. No conversion of panorama length to independent sample count or seconds is made here.

Source: Schepperle et al. (2024), Sensors 24, 3363, DOI 10.3390/s24113363, supports contour extraction and comparison against liquid appearance (§2.2); its experiment is rectangular and uses a separate liquid-only photograph. Our inferred reference and circular reconstruction are extensions and must be reported as such.
