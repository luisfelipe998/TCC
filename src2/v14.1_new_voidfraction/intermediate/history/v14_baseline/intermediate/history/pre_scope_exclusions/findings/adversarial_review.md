# Independent adversarial review of the panorama analysis

The liquid-reference silhouette method is a substantial improvement over raw bright/dark pixel counting: it distinguishes visually liquid gaps from enclosed gas bodies in the available panoramas. The geometry values remain **apparent, conditional spatial void estimates**. The principal numerical limitation is now cropped interface width, rather than the earlier brightness-polarity inversion.

All checks below use saved `src2/data_processed/*/__panorama_v11.png` images only. Reproduce them with `python3 src2/v14_new_voidfraction/scripts/adversarial_review.py`; [adversarial_checks.json](../results/validation/adversarial_checks.json) records the exact segmentation source SHA-256. The reviewed source hash is `daa2d99165520f387c578d23d939e638ca60a39e7c633183243a60cb267e66e5`.

The pre-reorganization final rerun of the reproducible checks uses source hash `9e4d7a7c1b5c9ef7bf2d18818cc2b298600f7e41f9006b6774e8b5c89221e835`. It reproduces the reported anchor and crop checks; the final version also exposes the open-interface span as an explicit sensitivity parameter.

## 1. Independent visual phase checks support the new segmentation

I inspected exact interior patches before computing their labels. These are deliberately clear gas or liquid interiors, not difficult boundary pixels or an unbiased segmentation test set. Coordinates are half-open intervals in the original panorama; they are retained for review.

| Panorama | x interval | y interval | Visual interpretation | New gas-mask fraction |
|---|---:|---:|---|---:|
| m07q01_001 | 180:250 | 24:65 | Liquid; not a reference run | 0.000 |
| m07q01_001 | 580:650 | 24:65 | Liquid; not a reference run | 0.000 |
| m06q01_001 | 455:515 | 24:65 | Liquid gap between plugs | 0.000 |
| m06q01_001 | 150:300 | 30:55 | Gas-body interior | 1.000 |
| m05q05_002 | 505:540 | 24:65 | Liquid; not a reference run | 0.000 |
| m05q05_002 | 580:610 | 30:55 | Gas-body interior; not a reference run | 1.000 |

All six results remain the same for ratio gates 0.60/0.65/0.70/0.75/0.80. No gate was fitted to these patches. The m07q01 and m05q05_002 checks are outside the three reference panoramas. This directly addresses the old bright-core failure: an almost uniform bright liquid interval is no longer interpreted as a large gas core. These selected checks support phase labeling, not exact physical gas diameter.

The automatically inferred reference was also replaced independently with the radial profile of each visibly liquid patch: m06q01 x455:515, m05q05_002 x505:540 and m07q01 x180:250. The resulting conditional estimates were:

| Run | Range across automatic and three manual liquid references |
|---|---:|
| m07q01_001 | 0.001836–0.001867 |
| m05q05_002 | 0.30785–0.31233 |
| m06q01_001 | 0.46525–0.47090 |
| m06q05_001 | 0.47561–0.48168 |
| m06q20_001 | 0.64298–0.65634 |
| m06q80_001 | 0.74908–0.75572 |
| m03q20_001 | 0.780092–0.780121 |

The agreement supports using the liquid intervals already present in the images. References introduced for this perturbation are sensitivity checks; their own images are not independent negative controls in those variants. Near-constant m03q20 values are caused primarily by crop clipping, so their apparent reference robustness must not be mistaken for diameter accuracy.

## 2. A concrete false gas bridge was identified and fixed

The initial silhouette method accepted paired upper/lower shadows even when they belonged to different bubbles. In the visibly aerated liquid slug of m06q05, x300:560, it filled liquid between unrelated objects at x312:315, x421:429 and x448:455. For example x313 had its upper extremum at row26 and lower extremum at row54, in different connected dark components, but the original method filled all 29 intervening pixels. This was a false gas bridge visible in the image, not merely an absence of external calibration.

The revised estimator retains pairs with the same connected outline or a sufficiently long open-interface interval (20 columns), preserving continuous annular cores. I checked the corrected source: no different-component paired acceptance remains in that liquid-slug test interval. The reviewable image is [adversarial_slug_gap_mask.png](../figures/validation/adversarial_slug_gap_mask.png). The 20-column rule is an optical continuity heuristic; it should be reported, and can still miss fragmented true outlines or retain a long false bridge.

## 3. Crop clipping limits numerical diameter estimates, especially annular

The default mask occupies rows18:71: at most53 pixels. With D_app=60 the largest connected circular reconstruction is `(53/60)^2 = 0.780278`. A full-height mask consequently gives this same number irrespective of the actual gas core beyond the crop. This affects most annular columns and some churn columns.

I repeated segmentation with the same reference, same threshold and D_app=60, changing only the radial mask crop:

| Run | Rows20:70 | Rows18:71, default | Rows16:73 | Rows15:75 |
|---|---:|---:|---:|---:|
| m07q01_001 | 0.00170 | 0.00184 | 0.00190 | 0.00192 |
| m05q05_002 | 0.29939 | 0.30785 | 0.30995 | 0.31167 |
| m06q01_001 | 0.44259 | 0.46555 | 0.47751 | 0.47761 |
| m06q05_001 | 0.45836 | 0.47561 | 0.48299 | 0.48332 |
| m06q20_001 | 0.61263 | 0.64997 | 0.67934 | 0.68585 |
| m06q80_001 | 0.68308 | 0.75265 | 0.81917 | 0.85824 |
| m03q20_001 | 0.69444 | 0.78009 | 0.89746 | 0.98009 |

For m06q80, crop alone changes alpha by approximately0.175; for m03q20 by0.286. Those changes dwarf the template and ratio-gate effects. A wider crop is not automatically correct: it can include wall/reflection shadows. The table establishes how much the **apparent geometry** depends on where the optical boundary is placed. The recorded missing-width intervals are conditional truncation bounds under the disk model, not physical confidence intervals.

Changing D_app from60 to56 multiplies every diameter-squared estimate by1.14796; changing it to64 multiplies by0.878906. These changes express denominator uncertainty. The reference's wall/reflection minima do not independently establish the actual lumen diameter or the optical gas-interface mapping.

## 4. The reconstructed mask supports a better-defined spatial flip-rate

The new phase-mask flip calculation counts actual mask transitions with an N−1 denominator and retains entries/exits separately. Its complement and reversal invariants are appropriate. Use the saved counts and denominators rather than labeling the old v11.3 number a Schmitt-trigger rate or a count per1000 pixels. The old and new mask rates describe different binary fields.

Center-row flips of the reconstructed silhouette count changes in center gas occupancy. They distinguish passage of gas bodies from fluctuations of brightness inside a persistent core. However, center-row flips miss off-center dispersed bubbles, and filling an envelope can hide liquid droplets or fragments within it. A full-tube volume fraction and phase intermittency require those fragments to be resolved too.

Void and flip-rate are not statistically independent evidence when extracted from the same mask. If there are n_g gas samples and n_l liquid samples in a binary record of N columns, the maximum possible flip count is at most `2 min(n_g,n_l)`; occupancy near1 necessarily limits flips. Annular's near-continuous center gas core therefore produces a small phase-flip rate by construction. That is a legitimate phase-occupancy statement, but does not independently validate its estimated gas diameter or establish absence of waves/chaos.

## 5. Recommended interpretation of the numerical answer

Report the new per-run **apparent circular reconstruction** alongside the projected silhouette and spatial phase-flip count/rate. Keep the actual supplied labels for both bubble acquisitions and show their very different observed morphology. The phase checks above make this more than a ranking chosen to match theory: liquid and gas interiors are correctly separated under the inspected controls, and a demonstrated false bridge was corrected.

Do not call the default annular value around0.77 a calibrated true void fraction: many columns sit near the artificial0.7803 ceiling. Report crop/diameter sensitivity prominently; threshold intervals alone substantially understate uncertainty. Plug/slug estimates are more stable to crop in the inspected examples, but their three-run groups overlap and remain conditional on circular cross-sections. Irregular churn and depth overlap remain less constrained. A correct image-only result is a conditional optical/geometric estimate with explicit limitations, not a universal lookup of exact volumetric alpha by regime.

<!-- v14-reorganization-provenance:start -->
## Current automated checks after folder reorganization

The current segmentation source SHA-256 is `0d1ddf2d1222fe6677d52733626ea5205f23b2a8edf3e3b39be5d750fa82d70c`. The regenerated phase-mask, silhouette-classifier, adversarial and final-verification metadata all match that source. Compared with the preserved pre-move snapshot, CSV tables and all binary masks are byte-identical; numerical JSON contents are unchanged. The recorded changes are source-hash metadata caused by relocating the code.

See the [preservation and navigation check](layout_verification.md), [current adversarial results](../results/validation/adversarial_checks.json) and [final arithmetic verification](../results/validation/final_verification.json). This note records automated rerun evidence. The manual visual review above belongs to the historical source hashes recorded earlier. Physical calibration and the stated crop/geometry limitations remain unresolved.
<!-- v14-reorganization-provenance:end -->
