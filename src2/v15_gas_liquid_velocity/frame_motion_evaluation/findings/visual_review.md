# Visual review of native frame sequences

## Executive summary

**The selected native sequences visibly preserve moving structures between frames.** Clear bubble and elongated-interface examples justify attempting displacement estimates in px/frame and px/s. The same images do not establish independent liquid and gas material velocities. Churn shows substantial evolving texture, and the annular examples contain broad moving brightness and saturated regions that can confound correspondence.

This is a targeted visual review of **six recordings**, covering all five retained labels and a second annular example with known saturation. It is not a visual certification of all 26 retained recordings or an independent velocity ground truth. The original classifications are retained.

## Reproducible visual evidence

The contact sheets display exact native frames **200, 202, 204, 206, 208 and 210** at their original 1016 × 80 dimensions and stored grayscale levels. Brightness is not normalized or contrast-enhanced. Labels give elapsed time relative to frame 200. The x–time figures use **all 501 frames**, averaging only native rows 20:60, with a fixed grayscale range of 0–255. No panorama, temporal interpolation or motion correction is used.

| Recording and retained label | Visible evidence in the selected sequence | Interpretation limit |
|---|---|---|
| [m07q01_001 — dispersed bubbles](../figures/visual_sequence_m07q01_001.png) | Several separated, recognizable bubble outlines move toward smaller x coordinates. Shapes and transverse positions also change. | Bubble tracking is plausible; centroid or boundary motion is not automatically the gas mass-weighted mean velocity. Small dark features are not identified liquid tracers. |
| [m05q01_001 — plug](../figures/visual_sequence_m05q01_001.png) | The elongated central object and its rounded boundary remain recognizable over ten frame intervals and move left. A liquid interval shares the neighboring gas-object boundaries. | Strong evidence for interface translation. Tracking the liquid interval's boundaries would reuse the same interfaces rather than independently track liquid parcels. |
| [m05q05_001 — slug](../figures/visual_sequence_m05q05_001.png) | Elongated objects and a pair of smaller irregular objects retain visible correspondence while moving left. Rear boundaries and widths deform. | A single translation is an approximation; different interfaces can have different apparent speeds. These camera observations require no assignment of the unmatched LabVIEW sensor blocks. |
| [m05q20_001 — churn](../figures/visual_sequence_m05q20_001.png) | Broad structures move left, but internal brightness, small interfaces and shapes evolve appreciably across the sequence. | Local short-interval correspondence may be possible. Long tracks or a single phase-wide speed require stronger checks than visual resemblance. |
| [m04q80_001 — wavy annular](../figures/visual_sequence_m04q80_001.png) | A bright central band and surrounding fine texture change and move through a mostly continuous core. Individual material parcels are not identifiable in this view. | Correlation can follow a brightness envelope or interfacial texture. A successful local match would remain a structure/wave-motion measurement. |
| [m03q40_001 — wavy annular](../figures/visual_sequence_m03q40_001.png) | Large bright regions have little resolved internal contrast; broad intensity patterns evolve. | Saturation removes information. This recording has **23.844%** of central-band native pixels at grayscale ≥250 over all 501 frames, exceeding the pre-existing 20% warning threshold. |

The [dispersed-bubble x–time image](../figures/visual_xt_m07q01_001.png) has narrow inclined traces, and the [plug x–time image](../figures/visual_xt_m05q01_001.png) has persistent inclined boundaries. These are visible evidence of directed image motion. The [churn x–time image](../figures/visual_xt_m05q20_001.png) contains both broad inclined structure and less persistent fine texture. A line slope could provide a further estimator, but no slopes were fitted in this visual review.

## What this does and does not validate

- These examples support trying displacement measurement directly on consecutive native frames. At 7,000 fps, the displayed 200→210 span is 1.429 ms; it is 2.000 ms at 5,000 fps and 2.500 ms at 4,000 fps. Native timing must be used for each run.
- The leftward motion establishes the sign of these examples in the image coordinate system. It does not by itself define a laboratory upstream/downstream axis.
- No verified liquid-phase tracers were identified in these six examples. The visible small dots cannot be assigned to liquid material without establishing whether they are bubbles, particles, reflections or another optical feature. This is a limitation of identification, not proof that every frame in the dataset lacks a potentially useful liquid feature.
- Gas and liquid masks meet at shared interfaces. Comparing their moving boundaries does not create two independent material velocity measurements.
- This review does not rule out all local annular tracking. It identifies failure mechanisms that local patch tests must check. It also does not turn a high correlation or a visually convincing translation into an accuracy percentage.
- The raw native frames and full-duration x–time figures are the evidence. The contact sheets select a short interval for readability; they are not representative sampling for numerical acceptance rates.

## Reproduction and provenance

Run the following from the repository root:

    python3 -B src2/v15_gas_liquid_velocity/frame_motion_evaluation/scripts/visual_sequences.py

The script checks that each selected recording contains contiguous image indices 0–500 and the expected 80 × 1016 shape, obtains the camera rate from its native MII configuration, and saves sequence-frame hashes, times, dimensions and central-band saturation fractions in [visual_sequences_metadata.json](../results/visual_sequences_metadata.json). Source images and existing v15 outputs are not modified.
