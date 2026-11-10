# Independent review of phase-conditioned motion

## Executive result

**The conditional calculations reproduce, but the stricter phase-interior checks do not identify a usable liquid/gas material velocity ratio.** Four recordings support the declared gas-passage/liquid-looking-passage descriptor with threshold consensus, a liquid-reference guard and lag-2 motion confirmation. All four are plug/slug recordings. None passes the stricter local ratio requirements when both matched patches must have at least 95% phase-consensus occupancy; none passes after additionally restricting to eroded interiors.

The independent verifier checked **1,177,284 numeric values**, with zero discrepancy at the saved precision. This verifies mask representation, phase-state decisions, persistent-episode grouping, motion selection, means and ratios, including the later probe-position and liquid-rich-gate sensitivity tests. It does not establish correct physical phase labeling.

| Checked evidence | Scope | Result |
|---|---:|---|
| Native mask representations and inventory | 39,078 threshold masks from 13,026 frames | Packed/unpacked shapes, lumen constraints, counts, fractions and threshold disagreement reproduce. |
| Fixed-probe states and motion gating | 65,130 frame-state rows and 130,000 comparisons | Every state, short-episode rejection, accepted pair and rate conversion reproduces. |
| Conditional aggregates | 3,016 episode rows, 520 phase summaries, 260 ratio rows | Counts, means, episode-equal means, ratios and withholding reasons reproduce. |
| Phase-labeled local patches | 46,800 tracks, 312 phase summaries, 156 ratios | Source/target mask purity, erosion, assignment, frame-equal means and ratio decisions reproduce. |
| Probe-position / liquid-rich sensitivity | 260 ratio rows and 3,594 episode rows | Probe occupancy, episodes, counts, means, ratios and all eligibility decisions reproduce; the center-strict baseline matches the main consensus analysis. |

The liquid-reference error was independently recomputed from **all 13,026 original images** using an explicit Gaussian convolution. The 7×7 interior erosion was independently reconstructed as a spatial minimum filter, rather than importing the analysis erosion call. The verifier imports neither motion-analysis script.

## What the actual ratios show

The following table uses threshold consensus, the liquid-reference guard and lag-2 confirmation. The first ratio divides the two accepted-frame means; the second gives each accepted episode equal weight before division. Both are conditional image-motion descriptors.

| Run | Pattern | Gas pairs / episodes | Liquid-looking pairs / episodes | Ratio of frame means | Ratio of episode-equal means |
|---|---|---:|---:|---:|---:|
| m05q01_001 | plug | 273 / 9 | 59 / 8 | 0.979927 | 0.983224 |
| m05q05_001 | slug | 253 / 13 | 56 / 9 | 0.885696 | 0.934547 |
| m06q01_001 | plug | 245 / 15 | 31 / 6 | 0.994581 | 0.988144 |
| m07q05_001 | plug | 221 / 18 | 42 / 8 | 0.988544 | 1.008465 |

The averaging convention changes the slug example appreciably and moves one plug descriptor across one. Values below one must not be interpreted as evidence that gas material travels more slowly than liquid material: the matched global template can contain both phases, and the estimates describe different selected image events.

The phase-interior challenge is substantial. For background-subtracted local patches, accepted liquid-looking tracks decrease from **236** at 80% single-threshold purity to **44** at 95% threshold-consensus purity, then **4** with 100% occupancy of eroded consensus interiors. Those four interior tracks all occur in the dispersed-bubble recording. No recording has the required information for both phase means under either stricter local variant. Gas-candidate tracks remain more numerous: **1,725 → 1,314 → 707**.

This loss of liquid-side information is consistent with interface contamination and weak interior texture. It does not uniquely prove the cause of every earlier match, because strict masks also shrink the available area and change which samples survive. It leaves independent liquid material motion unsupported by the tested setup.

## Probe-position and liquid-rich-gate sensitivity

| Variant | Eligible adjacent ratios | Eligible lag-2-confirmed ratios |
|---|---:|---:|
| Center, inferred gas occupancy ≤5% for liquid-looking state | 4 | 4 |
| Left probe, same strict phase gate | 3 | 3 |
| Right probe, same strict phase gate | 1 | 0 |
| Center, inferred gas occupancy ≤20% for liquid-rich state | 6 | 6 |
| Center, inferred gas occupancy ≤40% for liquid-rich state | 6 | 6 |

All three probes are 24 pixels wide: x=496:520, x=208:232 and x=624:648, respectively. The global correlation template occupies x=240:776 in its source frame. In particular, the left probe is outside that template; moving the probe changes which remote global shifts are conditioned, rather than creating a co-located material velocity measurement. Probe sensitivity is therefore an important warning about the observable's definition, not evidence of true phase acceleration between stations.

Relaxing the liquid gate restores two sensitive slug recordings, m04q01_001 and m06q05_001. Their lag-2-confirmed descriptors are approximately **0.398–0.399** and **0.751–0.758**, respectively. Those newly accepted intervals can contain 20% or 40% inferred gas occupancy. They cannot be relabeled as pure-liquid measurements. The added numerical availability does not resolve phase identity or explain these values as physical slip.

The right strict probe yields no lag-2-confirmed reportable ratios. That means the stated gates lack adequate evidence for both phase candidates there; it is not evidence that liquid or gas is physically absent from the experiment.

## Physical meaning of the proposed calculation

Separately averaging accepted image-motion estimates during inferred gas-dominated and liquid-looking intervals is a mathematically valid conditional descriptor. It becomes a gas/liquid material velocity ratio only if each underlying displacement demonstrably tracks the corresponding phase's material. The available segmentation and image correspondence do not establish that identification by themselves.

For a gas mask B and its liquid complement L=1−B, the optical-flow equations are identical: L_t+v·∇L=−(B_t+v·∇B). Complementary masks therefore do not provide two independent velocity constraints. A liquid interval bounded by moving bubbles inherits those gas–liquid boundaries; its displacement does not determine the motion of liquid parcels inside it. The same issue affects correlation of complementary masks, whose mean-subtracted arrays differ only by sign.

Conditioning a broad grayscale-profile velocity on the phase at a narrow fixed probe has an additional interpretation limit: the matched profile can contain both phases and remote interfaces, even while the probe is liquid-looking. Such results should be named motion during inferred gas/liquid passages. Requiring both the source and matched destination patch to lie within the same phase candidate, and excluding interface neighborhoods, tests a stricter identification hypothesis. It may reveal a lack of trackable interior texture rather than supply two velocities.

## Phase identification and pattern-specific handling

- The current v14 geometry records the owner-confirmed lumen as rows 18:71, 53 pixel cells. The older executive-summary text still describes D=60; the current mask code and provenance should control this experiment. No circular-area reconstruction is needed for temporal phase gating.
- The v14 shadow-envelope detector marks inferred gas silhouettes. Its unmarked complement is not independently verified liquid: it can contain missed gas, unresolved bubbles, saturated interiors, or ambiguous optical regions. A three-state label—gas candidate, liquid-looking candidate, unresolved—is preferable to forced binary assignment.
- In dispersed-bubble flow, a whole-section gas-dominance gate can miss every gas event despite visible individual bubbles. Track identified bubble objects separately; do not infer liquid speed from the background or claim that absent gas-dominated intervals mean absent gas.
- Plug/slug sequences offer the clearest alternation. Exclude entry/exit transitions and verify phase consistency over the actual motion-support region. Different front and rear speeds can reflect growth/deformation rather than slip.
- Churn has rapid topology changes, mixed projected occupancy and occlusion. Retain an unresolved class and short event-level tracks. A forced full-record ratio can conceal rejection and selection bias.
- Annular flow contains gas core and liquid film simultaneously. A central temporal probe can remain gas-dominated for the full record while liquid is present at the wall. There may be no liquid-passing interval to average. Region-specific wave/interface motion is a different observable and needs separate interpretation.
- Saturation removes phase-identification information and cannot be repaired by calling the saturated region liquid or assigning zero velocity.

## Averaging, ratios and outliers

Define the estimator explicitly. A mean over accepted frame pairs weights long, trackable passages more heavily; a mean of episode means weights passages equally. Neither automatically gives a phase mass-weighted or volume-weighted velocity. The phase means can come from different times and different spatial content. Report their denominators, episode counts, acceptance fractions and spatial support.

Compute per-recording ratios first. A ratio of means differs from a mean of instantaneous ratios, and pooling frames across recordings makes long/easily tracked records dominate a pattern result. Pattern summaries should retain per-recording values, the number of contributing recordings, and the reason for every missing ratio. One dispersed-bubble acquisition does not establish pattern-level variability.

Preserve signed velocities. A denominator close to zero, opposite directions, or cancellation of changing signs can make a ratio unstable or physically different. Replacing all velocities by absolute values can conceal that problem. Require adequate information in both phase candidates before reporting a ratio. Absence of usable liquid-looking tracks is an unidentified denominator, not zero liquid velocity.

Threshold sensitivity, transition exclusion and interface-erosion sensitivity should be retained as actual alternative analyses. Do not choose a setting because its ratio matches an expected slip range or produces cleaner pattern clusters. Outliers should first be tied to correspondence, boundary crossings, low texture, saturation, missing phase passages or sampling; magnitude alone is not a reason to delete them.

The 501 frames span only 0.07143–0.125 seconds. Adjacent comparisons overlap, and accepted phase passages are temporally correlated. Frame-by-frame resampling cannot justify a narrow independent-sample confidence interval. Episode-level or temporal-block summaries are more informative; very few episodes still do not support stable uncertainty estimates.

Within a recording, a common axial pixel scale cancels from a ratio of two image velocities. This does not resolve phase identification or phase-dependent optical distortions. Across recordings, equal image dimensions do not prove equal magnification.

## Reproduction and mathematical control

Run after the native-mask and conditional-motion analyses:

    python3 -B src2/v15_gas_liquid_velocity/phase_conditioned_motion/scripts/independent_verify.py

The [verification JSON](../results/independent_verification.json) preserves every eligibility-group result and acceptance total. The [native precheck](../results/independent_phase_precheck.json) records projected mask occupancy and threshold disagreement, with no claim of true void fraction.

A separate [complementary-mask control](../results/mask_complement_control.json) translates binary interfaces by +17 px. Direct normalized correlation estimates +17 px for both the gas mask and its liquid complement, giving ratio 1. Their complete correlation curves differ by at most 1.22×10⁻¹⁵. This demonstrates the shared-information problem mathematically; it does not predict the experimental grayscale-conditioned ratios. Reproduce with scripts/mask_complement_control.py.
