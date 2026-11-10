# Panorama shadow-envelope reconstruction

The current scope contains 26 acquisitions across 5 retained patterns, with 25 acquisitions passing the saturation criterion in grouped statistics. Mixed labels, smooth annular and the user-excluded replicate are omitted before processing. Supplied retained labels are unchanged.

**A usable image solution is possible as a conditional geometric estimate.** This method identifies dark gas outlines against a liquid-like radial template, fills the gas interior, and includes separately enclosed off-center bubbles. It does not equate dark pixels alone, or bright pixels alone, with gas. The result is a physically meaningful reconstruction formula applied to an optical silhouette; its absolute calibration remains conditional.

## Reference inferred from available panoramas

The three visually reviewed plug panoramas `m05q01_001`, `m06q01_001`, `m07q05_001` provide liquid-looking candidate intervals. In each, normalize radial profiles by the 75th percentile of rows 30:50; select columns at or below the 5th percentile of radial standard deviation in rows 20:68. Median selected profiles give the liquid-like template. Candidate columns are selected from profile structure; the three reference panoramas were chosen by inspection. `m07q01_001` is excluded from reference construction but was inspected during detector development, so it is an internal check rather than independent validation. The JSON saves exact candidate column coordinates.

The user confirmed the complete tube interior as Python slice 18:71: rows 18–70, 53 pixel cells. We therefore use D=53 px. The profile minima at rows 15/75 are optical extrema, not the lumen boundaries. Their separation no longer defines D. This corrects image normalization but does not calibrate interface refraction.

## Segmentation and physical reconstruction

Smooth normalized radial profiles by σ=1 px and divide by the inferred template. Locate paired upper/lower shadow deficits below ratio 0.70, expand each contiguous deficit to its outer shadow edge, and fill between them. Reject paired fills from different connected shadow components when their accepted axial interval is shorter than 20 px; this reduces false bridges but does not prove phase identity. Long open interfaces remain eligible. Separately, apply 3×3 morphological closing on a zero-padded crop, restore the crop, reject components touching a radial edge, retain components of at least 4 pixels, and fill enclosed holes. Axial record-edge contact is allowed. Components need not contain a hole: unresolved dark spots remain gas candidates, and genuine bubbles connected to a radial edge can be missed. Identical parameters apply to every panorama.

The projected gas-silhouette fraction is A_g,projected/(D L_record). The **conditional circular void estimate** is E_z[Σ_i(d_i(z)/D)²], where d_i is each separate silhouette's radial width. This follows from A_g(z)=Σ_i π d_i²/4 and A_tube=πD²/4. It requires circular/axisymmetric gas cross-sections and negligible overlap in viewing depth. The arithmetic mean of squares is used, not the median width or square of its mean. The px→mm factor cancels in the ratio.

## Per-run estimates and sensitivity

These numbers are **conditional optical reconstructions**, not independently calibrated measurements. Saturated runs are retained and explicitly marked unsupported. The threshold interval is a perturbation range, not a statistical confidence interval.

| Run | Pattern | Conditional α | Threshold 0.60–0.80 range | Reference variant range | Fixed-row y=40 gas duty | Status |
|---|---|---:|---:|---:|---:|---|
| m03q05_001 | churn | 0.966 | 0.955–0.974 | 0.966–0.966 | 0.999 | conditional_optical_geometry |
| m03q10_001 | wavy annular | 0.996 | 0.993–0.997 | 0.996–0.996 | 1.000 | conditional_optical_geometry |
| m03q20_001 | wavy annular | 1.000 | 0.999–1.000 | 1.000–1.000 | 1.000 | conditional_optical_geometry |
| m03q40_001 | wavy annular | 0.980 | 0.965–0.991 | 0.979–0.982 | 1.000 | unsupported_saturation |
| m04q01_001 | slug | 0.705 | 0.675–0.725 | 0.705–0.706 | 0.829 | conditional_optical_geometry |
| m04q10_001 | churn | 0.885 | 0.866–0.902 | 0.885–0.886 | 0.989 | conditional_optical_geometry |
| m04q20_001 | churn | 0.903 | 0.883–0.920 | 0.903–0.904 | 0.996 | conditional_optical_geometry |
| m04q40_001 | wavy annular | 0.981 | 0.976–0.986 | 0.981–0.982 | 1.000 | conditional_optical_geometry |
| m04q80_001 | wavy annular | 0.999 | 0.999–1.000 | 0.999–0.999 | 1.000 | conditional_optical_geometry |
| m05q01_001 | plug | 0.636 | 0.618–0.653 | 0.636–0.640 | 0.705 | conditional_optical_geometry |
| m05q05_001 | slug | 0.522 | 0.500–0.547 | 0.521–0.523 | 0.680 | conditional_optical_geometry |
| m05q10_001 | churn | 0.762 | 0.728–0.799 | 0.762–0.764 | 0.945 | conditional_optical_geometry |
| m05q20_001 | churn | 0.816 | 0.771–0.850 | 0.816–0.819 | 0.972 | conditional_optical_geometry |
| m05q40_001 | churn | 0.931 | 0.895–0.952 | 0.930–0.933 | 0.993 | conditional_optical_geometry |
| m05q80_001 | wavy annular | 0.985 | 0.958–0.993 | 0.985–0.987 | 0.999 | conditional_optical_geometry |
| m06q01_001 | plug | 0.593 | 0.565–0.623 | 0.593–0.595 | 0.738 | conditional_optical_geometry |
| m06q05_001 | slug | 0.606 | 0.572–0.647 | 0.606–0.608 | 0.799 | conditional_optical_geometry |
| m06q10_001 | churn | 0.689 | 0.648–0.742 | 0.689–0.691 | 0.872 | conditional_optical_geometry |
| m06q20_001 | churn | 0.829 | 0.775–0.876 | 0.827–0.831 | 0.935 | conditional_optical_geometry |
| m06q40_001 | wavy annular | 0.921 | 0.860–0.956 | 0.918–0.924 | 0.976 | conditional_optical_geometry |
| m06q80_001 | wavy annular | 0.964 | 0.917–0.980 | 0.962–0.965 | 0.993 | conditional_optical_geometry |
| m07q01_001 | dispersed bubbles | 0.002 | 0.001–0.003 | 0.002–0.002 | 0.014 | conditional_optical_geometry |
| m07q05_001 | plug | 0.475 | 0.447–0.505 | 0.474–0.476 | 0.659 | conditional_optical_geometry |
| m07q20_001 | churn | 0.654 | 0.608–0.706 | 0.652–0.656 | 0.878 | conditional_optical_geometry |
| m07q40_001 | churn | 0.812 | 0.751–0.859 | 0.811–0.814 | 0.948 | conditional_optical_geometry |
| m07q80_001 | wavy annular | 0.926 | 0.858–0.957 | 0.922–0.930 | 0.975 | conditional_optical_geometry |

## Adversarial checks

The sensitivity CSV repeats thresholds 0.60/0.65/0.70/0.75/0.80, leaves out each of the three reference panoramas in turn, changes reference candidate selection to 1% and 10%, narrows the radial crop, and changes the open-interface minimum interval from 20 to 10/30 px. The dispersed panorama is never a reference. Five spatial-block means expose within-record variability. These checks measure internal robustness; they cannot substitute for independent void calibration.

Liquid-like candidate checks are saved separately, with reference construction excluding the tested panorama. Projected foreground fractions: m05q01_001 0.0191%, m06q01_001 0.1064%, m07q05_001 0.0308%. Conditional area contributions: m05q01_001 0.0000108, m06q01_001 0.0000548, m07q05_001 0.0000135. The candidates were inferred from the same dataset and are not certified all-liquid ground truth.

The overlay figure makes possible false phase assignments reviewable. The dispersed internal check should show small bubbles with little foreground between them. Plug/slug images should show both shadow boundaries and bright interiors, with liquid slugs empty. Annular should keep the connected core, while churn retains uncertain fragmented outlines. Fixed-row duty uses image row y=40, not a validated tube center; the midpoint of the confirmed pixel-center rows 18–70 is y=44.

**Complete-interior normalization:** both the maximum mask width and D are 53 pixels. A fully filled lumen gives α=(53/53)²=1. The previous 0.7803 maximum resulted from using D=60 with this 53-row ROI and is superseded. Border contact alone does not imply gas outside the known lumen. The previous missing-width extension scenarios and D56/D64 alternatives are retired. Values near 1 still cannot prove a gas-only cross-section: finite resolution can hide liquid films, and shadow envelopes, refraction and circular reconstruction remain conditional.

## Remaining physical assumptions

A single corrected side view can reconstruct an axisymmetric gas body; a second view and absolute mm scale are not invariably necessary. The present cylindrical tube's refraction can move shadow edges away from physical interface positions. The glass/illumination inverse map remains unknown, and irregular churn or overlapping dispersed bubbles can violate circular reconstruction. Those are the principal limits on calling the numerical result the true void fraction.

A stitched panorama supplies a spatial record. Its mean is spatial occupancy under the panorama reconstruction's frozen-advection assumption; it is not automatically a time average at a fixed point. A physical time average requires native frame timestamps and fixed-column phase masks. No conversion of panorama length to independent sample count or seconds is made here.

Source: Schepperle et al. (2024), Sensors 24, 3363, DOI 10.3390/s24113363, supports contour extraction and comparison against liquid appearance (§2.3); its experiment is rectangular and uses a separate liquid-only photograph. Our inferred reference and circular reconstruction are extensions and must be reported as such.
