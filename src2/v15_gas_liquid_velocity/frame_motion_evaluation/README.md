# Native-frame motion feasibility evaluation

Start with the [EXECUTIVE SUMMARY](findings/EXECUTIVE_SUMMARY.md). This extends v15 using the existing 501-frame recordings to evaluate structure displacement in px/frame and px/s. It does not estimate thermodynamic quality or independently measured liquid/gas material velocities.

The 26-run scope and original data are unchanged. Earlier v15 numerical artifacts are preserved. New computations, figures and findings reside in this directory.

## Reproduce

From the repository root, using the existing Python environment:

    python3 -B src2/v15_gas_liquid_velocity/frame_motion_evaluation/scripts/temporal_motion.py
    python3 -B src2/v15_gas_liquid_velocity/frame_motion_evaluation/scripts/local_patch_motion.py
    python3 -B src2/v15_gas_liquid_velocity/frame_motion_evaluation/scripts/local_patch_report.py
    python3 -B src2/v15_gas_liquid_velocity/frame_motion_evaluation/scripts/visual_sequences.py
    python3 -B src2/v15_gas_liquid_velocity/frame_motion_evaluation/scripts/independent_verify.py
    python3 -B src2/v15_gas_liquid_velocity/frame_motion_evaluation/scripts/build_summary.py

The temporal script imports the existing v15 NCC functions; their hash is recorded. The independent verifier does not import either analysis implementation. The summary script reads completed numerical outputs and the independent ROI/lag audit. No new package installation or network access is needed.

## Outputs

- [Executive summary](findings/EXECUTIVE_SUMMARY.md): feasibility, actual results and recommendation.
- [Temporal evaluation](findings/temporal_review.md): all 26 scenarios and fixed correspondence rules.
- [Local 2-D evaluation](findings/local_patch_review.md): independent spatial correspondence test.
- [Native visual sequences](findings/visual_review.md): six contact sheets and six full-duration space–time figures.
- [Independent verification](findings/independent_review.md): arithmetic/provenance checks and interpretation limits.

Signed positive motion is rightward in the stored image. Medians and interquartile ranges describe accepted image-motion events; they are not measurement uncertainty or phase-weighted velocities. Rejected events remain in the CSVs. The reports separate numerical support, temporal consistency, region sensitivity and physical identification.
