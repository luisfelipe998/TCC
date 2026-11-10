# v15 — Gas/liquid velocity investigation

Start with the [EXECUTIVE SUMMARY](findings/EXECUTIVE_SUMMARY.md). It now presents the latest gas-passage/liquid-gap analysis in **px/s**, with dimensionless ratios, and is synchronized with the [detailed phase-conditioned summary](phase_conditioned_motion/findings/EXECUTIVE_SUMMARY.md). The [current comparison chart](phase_conditioned_motion/figures/conditional_passage_velocities.png) and its [plotted values](phase_conditioned_motion/results/conditional_passage_velocities.csv) use the same four accepted recordings and reference settings.

These are signed image-motion averages during inferred gas-like passages and liquid-looking gaps. They do not establish independent gas/liquid material velocities or slip. The latest evaluation tests all 501 frames in each retained recording, v14 phase candidates, phase-interior patches, smaller adaptive patches and dispersed-bubble tracking. Four recordings support the reference conditional comparison; no universal ratio is established for each pattern.

The earlier v12/v12.1, recorded-flow and thermodynamic investigations are preserved in the [conditional model summary](findings/conditional_model_summary.md) and [conditional model tables](findings/velocity_results.md). Their m/s values are predictions under explicit inlet/property/closure assumptions, not conversions of the current pixel-based measurements.

## Contents

| Folder | Contents |
|---|---|
| [findings](findings/README.md) | Current executive summary, detailed phase evidence and historical model/legacy/literature reviews |
| [phase_conditioned_motion](phase_conditioned_motion/README.md) | Latest phase-passage means, ratios, segmentation checks, adaptive patches and bubble trajectories |
| [frame_motion_evaluation](frame_motion_evaluation/README.md) | Native-frame and local-patch motion evidence used by the phase-conditioned analysis |
| [results](results/README.md) | Index to current pixel-based tables and preserved conditional model inputs/results |
| [figures](figures/README.md) | Current px/s comparison and phase diagnostics; separately labeled historical model plots |
| [scripts](scripts/README.md) | Reproducible computations, reports and independent checks |
| [literature](literature/README.md) | Source extracts, rendered evidence, pinned secondary implementation and arithmetic checks |
| intermediate | Frozen v14 alpha input, run logs and disposable plotting caches |

Only 26 of the 31 experimental camera folders are included: dispersed bubbles 1, plug 3, slug 3, churn 10, wavy annular 9. Excluded: m03q01_001, m03q80_001, m04q05_001, m05q05_002 and m07q10_001. Optical saturation is a validity flag, not an extra thermodynamic-scope exclusion. Synthetic folders are not experimental conditions.

All new work is stored here. Original v12/v12.1 and input data are preserved. The existing repository .gitignore excludes src2; these deliverables are local files and have not been committed.

## Reproduce

To rebuild the current summary and comparison chart from the saved phase results, run from the repository root:

    python3 -B src2/v15_gas_liquid_velocity/phase_conditioned_motion/scripts/build_summary.py

This also republishes the main executive summary with corrected relative links. To synchronize only the summary text from the existing phase report:

    python3 -B src2/v15_gas_liquid_velocity/scripts/publish_current_summary.py

To recompute the later experiments from images, first follow the complete [native-frame reproduction sequence](frame_motion_evaluation/README.md#reproduce), then the [phase-conditioned reproduction sequence](phase_conditioned_motion/README.md#reproduction-order). The phase scripts consume the native-frame global/local motion results.

The original model/optical pipeline remains available:

    python3 -B src2/v15_gas_liquid_velocity/scripts/run_all.py

This rebuilds the original legacy, mass-flow, thermodynamic and optical analyses, generates the separately named conditional model report, republishes the **saved** phase-conditioned summary, and checks the artifacts. It does **not** recompute the later native-frame or phase-conditioned experiments. All reproduction uses local source files; no network access is needed. Source literature extracts and the pinned external implementation are already included. Independent legacy/literature narrative reviews remain audit artifacts.

Verified environment: Python 3.13.2; NumPy 2.3.4; pandas 2.3.3; SciPy 1.18.0; matplotlib 3.10.7; CoolProp 8.0.0; openpyxl 3.1.5; Pillow 12.0.0. Exact package versions are also in [requirements.txt](requirements.txt). No new package installation was required.

The model script freezes the v14 image-alpha input at [intermediate/v14_void_input_snapshot.csv](intermediate/v14_void_input_snapshot.csv). Ordinary reruns reuse that snapshot. Explicitly pass --refresh-inputs to model_comparison.py only when intentionally updating that diagnostic input, then regenerate summary/verification. The primary thermodynamic calculation does not depend on the v14 image alpha.

## Interpreting the output

- Current passage results: signed image motion in px/s and ratios of those conditional means. Phase-looking labels and accepted correspondence do not establish liquid/gas material identity. Missing results remain missing, with reasons.
- Current charts: px/s for passage means; dimensionless ratios for comparisons; selected counts and sensitivity describe available support, not calibrated accuracy.
- Historical Jg/Jl: conditional superficial velocities in m/s. Zero inlet quality, nominal heat input and local-property assumptions remain material.
- Historical ug/ul: conditional phase means in m/s from a named void-fraction closure. RA is that model investigation's reporting reference, not a proven accuracy winner.
- Recorded flow: owner-documented g/s and verified original samples, matched to camera records only by operating-condition name. Sensor/video synchronization and calibration uncertainty are unresolved.
- Earlier optical output: signed structure displacement in px/frame or px/s. The preferred candidate is temporal-background-subtracted NCC without high-pass filtering. Files named optical_primary_pairs and optical_per_run retain the high-pass reference experiment and are identified as such.
- Statistical summaries describe the sampled conditions; model spread, accepted-pair IQR and stress-test ranges are not measurement confidence intervals.

The final verification concerns reproducibility, arithmetic, scope, provenance and links. It does not turn conditional estimates into experimentally accurate measurements.
