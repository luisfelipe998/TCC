# Panorama shadow-envelope reconstruction

The current scope contains 26 acquisitions across 5 retained patterns, with 25 acquisitions passing the saturation criterion in grouped statistics. Mixed labels, smooth annular and the user-excluded replicate are omitted before processing. Supplied retained labels are unchanged.

**A usable image solution is possible as a conditional geometric estimate.** This method identifies dark gas outlines against a liquid-like radial template, fills the gas interior, and includes separately enclosed off-center bubbles. It does not equate dark pixels alone, or bright pixels alone, with gas. The result is a physically meaningful reconstruction formula applied to an optical silhouette; its absolute calibration remains conditional.

## Reference inferred from available panoramas

The three visually reviewed plug panoramas `m05q01_001`, `m06q01_001`, `m07q05_001` provide liquid-looking candidate intervals. In each, normalize radial profiles by the 75th percentile of rows 30:50; select columns at or below the 5th percentile of radial standard deviation in rows 20:68. Median selected profiles give the liquid-like template. Candidate columns are selected from profile structure; the three reference panoramas were chosen by inspection. `m07q01_001` is excluded from reference construction but was inspected during detector development, so it is an internal check rather than independent validation. The JSON saves exact candidate column coordinates.

The inferred reference has wall/reflection minima at image rows **15 and 75**, an apparent separation of **60 px**. We use D_app=60 px, and repeat with 56 and 64 px. These are apparent optical dimensions: the minima have not been validated as physical lumen edges.

## Segmentation and physical reconstruction

Smooth normalized radial profiles by σ=1 px and divide by the inferred template. Locate paired upper/lower shadow deficits below ratio 0.70, expand each contiguous deficit to its outer shadow edge, and fill between them. Reject paired fills from different connected shadow components when their accepted axial interval is shorter than 20 px; this reduces false bridges but does not prove phase identity. Long open interfaces remain eligible. Separately, apply 3×3 morphological closing on a zero-padded crop, restore the crop, reject components touching a radial edge, retain components of at least 4 pixels, and fill enclosed holes. Axial record-edge contact is allowed. Components need not contain a hole: unresolved dark spots remain gas candidates, and genuine bubbles connected to a radial edge can be missed. Identical parameters apply to every panorama.

The projected gas-silhouette fraction is A_g,projected/(D_app L_record). The **conditional circular void estimate** is E_z[Σ_i(d_i(z)/D_app)²], where d_i is each separate silhouette's radial width. This follows from A_g(z)=Σ_i π d_i²/4 and A_tube=πD²/4. It requires circular/axisymmetric gas cross-sections and negligible overlap in viewing depth. The arithmetic mean of squares is used, not the median width or square of its mean. The px→mm factor cancels in the ratio.

## Per-run estimates and sensitivity

These numbers are **conditional optical reconstructions**, not independently calibrated measurements. Saturated runs are retained and explicitly marked unsupported. The threshold interval is a perturbation range, not a statistical confidence interval.

| Run | Pattern | Conditional α | Edge-contact extension scenario | Threshold 0.60–0.80 range | Reference variant range | Fixed-row y=40 gas duty | Status |
|---|---|---:|---:|---:|---:|---:|---|
| m03q05_001 | churn | 0.754 | 0.754–0.972 | 0.745–0.760 | 0.753–0.754 | 0.999 | conditional_optical_geometry |
| m03q10_001 | wavy annular | 0.777 | 0.777–0.996 | 0.775–0.778 | 0.777–0.777 | 1.000 | conditional_optical_geometry |
| m03q20_001 | wavy annular | 0.780 | 0.780–1.000 | 0.780–0.780 | 0.780–0.780 | 1.000 | conditional_optical_geometry |
| m03q40_001 | wavy annular | 0.765 | 0.765–1.000 | 0.753–0.773 | 0.764–0.766 | 1.000 | unsupported_saturation |
| m04q01_001 | slug | 0.550 | 0.550–0.671 | 0.527–0.566 | 0.550–0.551 | 0.829 | conditional_optical_geometry |
| m04q10_001 | churn | 0.691 | 0.691–0.886 | 0.676–0.704 | 0.691–0.691 | 0.989 | conditional_optical_geometry |
| m04q20_001 | churn | 0.704 | 0.704–0.925 | 0.689–0.718 | 0.704–0.705 | 0.996 | conditional_optical_geometry |
| m04q40_001 | wavy annular | 0.766 | 0.766–0.989 | 0.761–0.769 | 0.766–0.766 | 1.000 | conditional_optical_geometry |
| m04q80_001 | wavy annular | 0.780 | 0.780–1.000 | 0.779–0.780 | 0.780–0.780 | 1.000 | conditional_optical_geometry |
| m05q01_001 | plug | 0.497 | 0.497–0.621 | 0.482–0.510 | 0.497–0.499 | 0.705 | conditional_optical_geometry |
| m05q05_001 | slug | 0.407 | 0.407–0.490 | 0.390–0.427 | 0.406–0.408 | 0.680 | conditional_optical_geometry |
| m05q10_001 | churn | 0.595 | 0.595–0.709 | 0.568–0.624 | 0.594–0.596 | 0.945 | conditional_optical_geometry |
| m05q20_001 | churn | 0.637 | 0.637–0.769 | 0.601–0.663 | 0.637–0.639 | 0.972 | conditional_optical_geometry |
| m05q40_001 | churn | 0.727 | 0.727–0.906 | 0.698–0.743 | 0.726–0.728 | 0.993 | conditional_optical_geometry |
| m05q80_001 | wavy annular | 0.769 | 0.769–0.976 | 0.747–0.774 | 0.769–0.770 | 0.999 | conditional_optical_geometry |
| m06q01_001 | plug | 0.463 | 0.463–0.570 | 0.441–0.486 | 0.463–0.464 | 0.738 | conditional_optical_geometry |
| m06q05_001 | slug | 0.473 | 0.473–0.566 | 0.446–0.505 | 0.472–0.474 | 0.799 | conditional_optical_geometry |
| m06q10_001 | churn | 0.538 | 0.538–0.669 | 0.506–0.579 | 0.538–0.539 | 0.872 | conditional_optical_geometry |
| m06q20_001 | churn | 0.647 | 0.647–0.839 | 0.605–0.683 | 0.645–0.648 | 0.935 | conditional_optical_geometry |
| m06q40_001 | wavy annular | 0.719 | 0.719–0.920 | 0.671–0.746 | 0.716–0.721 | 0.976 | conditional_optical_geometry |
| m06q80_001 | wavy annular | 0.752 | 0.752–0.967 | 0.716–0.765 | 0.751–0.753 | 0.993 | conditional_optical_geometry |
| m07q01_001 | dispersed bubbles | 0.001 | 0.001–0.001 | 0.001–0.002 | 0.001–0.001 | 0.014 | conditional_optical_geometry |
| m07q05_001 | plug | 0.370 | 0.370–0.414 | 0.349–0.394 | 0.370–0.372 | 0.659 | conditional_optical_geometry |
| m07q20_001 | churn | 0.510 | 0.510–0.595 | 0.474–0.551 | 0.509–0.512 | 0.878 | conditional_optical_geometry |
| m07q40_001 | churn | 0.634 | 0.634–0.794 | 0.586–0.670 | 0.633–0.635 | 0.948 | conditional_optical_geometry |
| m07q80_001 | wavy annular | 0.723 | 0.723–0.926 | 0.670–0.747 | 0.719–0.725 | 0.975 | conditional_optical_geometry |

## Adversarial checks

The sensitivity CSV repeats thresholds 0.60/0.65/0.70/0.75/0.80, leaves out each of the three reference panoramas in turn, changes reference candidate selection to 1% and 10%, narrows the radial crop, and changes the open-interface minimum interval from 20 to 10/30 px. The dispersed panorama is never a reference. Five spatial-block means expose within-record variability. These checks measure internal robustness; they cannot substitute for independent void calibration.

Liquid-like candidate checks are saved separately, with reference construction excluding the tested panorama. Projected foreground fractions: m05q01_001 0.0168%, m06q01_001 0.0940%, m07q05_001 0.0272%. Conditional area contributions: m05q01_001 0.0000084, m06q01_001 0.0000427, m07q05_001 0.0000106. The candidates were inferred from the same dataset and are not certified all-liquid ground truth.

The overlay figure makes possible false phase assignments reviewable. The dispersed internal check should show small bubbles with little foreground between them. Plug/slug images should show both shadow boundaries and bright interiors, with liquid slugs empty. Annular should keep the connected core, while churn retains uncertain fragmented outlines. Fixed-row duty uses image row y=40, not a validated tube center; the midpoint between apparent wall minima y=15 and y=75 is y=45.

**A censoring limit matters for annular:** filling all 53 rows in the interior crop imposes a computational maximum (53/60)²=0.7803. This follows from the chosen crop and diameter and is not a physical maximum of void fraction. Several annular runs approach it; apparent agreement or low threshold scatter must not be mistaken for precision. The legacy CSV columns `conditional_missing_width_alpha_lower/upper` now explicitly denote an edge-contact extension scenario: retain unflagged columns and let each flagged column extend up to a full tube, E[α_visible+1_edge(1−α_visible)]. This is only an upper bound if no gas is wholly unseen outside the crop and no gas objects are rejected or missed, in addition to correct phase labels, circular geometry and an undistorted diameter. A bubble entirely outside the crop can increase real void while leaving both scenario values unchanged. The interval is not a confidence interval or a general bound on total gas.

## Remaining physical assumptions

A single corrected side view can reconstruct an axisymmetric gas body; a second view and absolute mm scale are not invariably necessary. The present cylindrical tube's refraction can move shadow edges away from physical interface positions. The glass/illumination inverse map remains unknown, and irregular churn or overlapping dispersed bubbles can violate circular reconstruction. Those are the principal limits on calling the numerical result the true void fraction.

A stitched panorama supplies a spatial record. Its mean is spatial occupancy under the panorama reconstruction's frozen-advection assumption; it is not automatically a time average at a fixed point. A physical time average requires native frame timestamps and fixed-column phase masks. No conversion of panorama length to independent sample count or seconds is made here.

Source: Schepperle et al. (2024), Sensors 24, 3363, DOI 10.3390/s24113363, supports contour extraction and comparison against liquid appearance (§2.3); its experiment is rectangular and uses a separate liquid-only photograph. Our inferred reference and circular reconstruction are extensions and must be reported as such.
