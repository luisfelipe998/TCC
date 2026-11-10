# Adversarial phase, reference and crop review

Current scope: 26 acquisitions in the five approved single patterns. Only retained acquisitions contribute to the checks below. The [selection manifest](../results/validation/dataset_selection.json) documents exclusions.

Current segmentation source SHA-256: a63db1153e6f73d334ffce09416ab2794d8a360f5474f27aa48c6b9e9c098933.

## Selected visual interior checks

These existing gas/liquid interpretations come from the previous visual review; excluded-run patches are removed. Their current mask fractions are recomputed below.

| Run | x interval | y interval | Visual interpretation | Gas-map fraction |
|---|---|---|---|---:|
| m07q01_001 | [180, 250] | [24, 65] | liquid, held-out run | 0.000000 |
| m07q01_001 | [580, 650] | [24, 65] | liquid, held-out run | 0.000000 |
| m06q01_001 | [455, 515] | [24, 65] | liquid plug gap | 0.000000 |
| m06q01_001 | [150, 300] | [30, 55] | gas body | 1.000000 |

**Limit:** selected clear patches are not an unbiased boundary-validation set. The retained gas-positive patch is from a plug reference panorama. Detector design used these data; all visual checks are internal development checks.

![Retained visual interior patches](../figures/validation/adversarial_phase_anchor_patches.png)

## Liquid-reference alternatives

Automatic reference construction still uses the same three retained plug panoramas. Manual alternatives now use only retained liquid patches.

| Run | Conditional alpha range across retained references |
|---|---:|
| m07q01_001 | 0.001768–0.001794 |
| m06q01_001 | 0.593138–0.599387 |
| m06q05_001 | 0.606295–0.613793 |
| m06q20_001 | 0.828813–0.836740 |
| m06q80_001 | 0.963743–0.967828 |
| m03q20_001 | 0.999762–0.999799 |

## Removal of known-interior rows

The user confirmed rows 18:71 as the complete tube interior: 53 rows and D=53 px. Boundary source: User-confirmed full tube interior on 2026-10-23: Python slice 18:71, rows 18 through 70, 53 pixel cells. Optical profile minima are not used as lumen boundaries. Refraction mapping remains uncalibrated. Every variant below uses that same diameter; smaller slices discard known-interior rows. Segmentation is rerun for each slice, so changes also include contour decisions near the artificial crop edge.

| Run | Rows 20:70 | Rows 19:70 | Rows 19:71 | Rows 18:71 — full lumen |
|---|---:|---:|---:|---:|
| m07q01_001 | 0.001359 | 0.001768 | 0.001768 | 0.001768 |
| m06q01_001 | 0.563705 | 0.573729 | 0.587768 | 0.593138 |
| m06q05_001 | 0.584368 | 0.592043 | 0.602513 | 0.606295 |
| m06q20_001 | 0.781266 | 0.794997 | 0.819604 | 0.828813 |
| m06q80_001 | 0.874763 | 0.903916 | 0.937621 | 0.963743 |
| m03q20_001 | 0.889987 | 0.925943 | 0.962412 | 0.999762 |

**Interpretation:** these variants test removal of rows from the confirmed lumen, not alternative physical diameters. Their analytic maximum is (included rows/53)²; the full-lumen maximum is 1. No outside-lumen expansion contributes to these ranges. Reflections outside the confirmed lumen are not additional gas area.

**Limit:** user-confirmed boundaries do not independently calibrate phase segmentation, refraction or circular reconstruction. A full gas mask does not prove a physically all-gas cross-section. These perturbation ranges are not confidence intervals.

Reference-profile minima near y=15 and y=75 are optical extrema, not the confirmed lumen boundaries. The fixed probe remains y=40; the midpoint of the included row centers is y=44. The probe is not called the physical center.

## Corrected false-bridge case

The liquid slug in m06q05_001, x300:560, has 0 accepted paired fills between different closed-outline components. The original false bridge and its correction remain documented in the [historical review](../intermediate/history/pre_scope_exclusions/findings/adversarial_review.md).

![Current corrected liquid-gap overlay](../figures/validation/adversarial_slug_gap_mask.png)

The connection/interval rule can still miss fragmented true gas outlines or retain a long false bridge. Both void and flip-rate come from the same inferred map; their agreement is not independent physical validation.

Use the [current machine checks](../results/validation/adversarial_checks.json), [final arithmetic verification](../results/validation/final_verification.json) and [post-regeneration adversarial review](post_regeneration_adversarial_review.md). The previous wider-scope evidence is historical, not a current result.
