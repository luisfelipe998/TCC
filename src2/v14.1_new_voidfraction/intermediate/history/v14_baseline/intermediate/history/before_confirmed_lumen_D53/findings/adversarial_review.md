# Adversarial phase, reference and crop review

Current scope: 26 acquisitions in the five approved single patterns. Only retained acquisitions contribute to the checks below. The [selection manifest](../results/validation/dataset_selection.json) documents exclusions.

Current segmentation source SHA-256: b8ea65b2838a3cc372760de6b05d1542d6677c2b334328498cb52c1ea3d98b61.

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
| m07q01_001 | 0.001380–0.001400 |
| m06q01_001 | 0.462812–0.467689 |
| m06q05_001 | 0.473078–0.478929 |
| m06q20_001 | 0.646705–0.652889 |
| m06q80_001 | 0.751987–0.755175 |
| m03q20_001 | 0.780092–0.780121 |

## Crop sensitivity

| Run | Rows 20:70 | Rows 18:71, default | Rows 16:73 | Rows 15:75 |
|---|---:|---:|---:|---:|
| m07q01_001 | 0.001060 | 0.001380 | 0.001712 | 0.001919 |
| m06q01_001 | 0.439846 | 0.462812 | 0.475083 | 0.477211 |
| m06q05_001 | 0.455970 | 0.473078 | 0.481133 | 0.482588 |
| m06q20_001 | 0.609604 | 0.646705 | 0.676224 | 0.683018 |
| m06q80_001 | 0.682559 | 0.751987 | 0.818113 | 0.857515 |
| m03q20_001 | 0.694437 | 0.780092 | 0.897443 | 0.980033 |

**Limit:** the default 53-row crop imposes the conditional maximum (53/60)^2 = 0.7803. Wider crops can include wall/reflection shadows. These ranges describe sensitivity, not physical confidence intervals.

## Corrected false-bridge case

The liquid slug in m06q05_001, x300:560, has 0 accepted paired fills between different closed-outline components. The original false bridge and its correction remain documented in the [historical review](../intermediate/history/pre_scope_exclusions/findings/adversarial_review.md).

![Current corrected liquid-gap overlay](../figures/validation/adversarial_slug_gap_mask.png)

The connection/interval rule can still miss fragmented true gas outlines or retain a long false bridge. Both void and flip-rate come from the same inferred map; their agreement is not independent physical validation.

Use the [current machine checks](../results/validation/adversarial_checks.json), [final arithmetic verification](../results/validation/final_verification.json) and [post-regeneration adversarial review](post_regeneration_adversarial_review.md). The previous wider-scope evidence is historical, not a current result.
