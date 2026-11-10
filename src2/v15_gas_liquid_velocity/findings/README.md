# Findings

Read [EXECUTIVE_SUMMARY.md](EXECUTIVE_SUMMARY.md) first. It presents the latest phase-passage results in px/s and dimensionless ratios, synchronized from the phase-conditioned investigation.

| Report | Purpose |
|---|---|
| [Detailed phase-conditioned summary](../phase_conditioned_motion/findings/EXECUTIVE_SUMMARY.md) | Current passage means/ratios, phase-identification limits, outliers and special handling |
| [Native phase masks](../phase_conditioned_motion/findings/native_phase_masks.md) | Applying v14 phase candidates to all 501 frames per recording |
| [Adaptive-patch review](../phase_conditioned_motion/findings/adaptive_patch_review.md) | Smaller patches within phase-candidate interiors and remaining liquid-identification failures |
| [Dispersed-bubble review](../phase_conditioned_motion/findings/closed_contour_review.md) | Tested detector repair and candidate bubble trajectories |
| [Phase independent review](../phase_conditioned_motion/findings/independent_review.md) | Numerical/provenance checks and physical interpretation limits |
| [Native-frame summary](../frame_motion_evaluation/findings/EXECUTIVE_SUMMARY.md) | Underlying global/local displacement feasibility in px/frame and px/s |
| [Historical conditional model summary](conditional_model_summary.md) | Earlier thermodynamic predictions in m/s under explicit assumptions |
| [Conditional model tables](velocity_results.md) | All 26 acquisitions, both flow baselines, model and input sensitivity; m/s predictions |
| [Model method comparison](method_comparison.md) | Physical definitions, tested models and unresolved experimental inputs |
| [Legacy audit](legacy_audit.md) | v12/v12.1 arithmetic reproduction and interpretation corrections |
| [Mass-flow review](massflow_review.md) | Engineering-unit provenance, original-file verification and camera-mapping limits |
| [Optical investigation](optical_investigation.md) | Four estimators, crop/lag tests, synthetic verification and validity flags |
| [Literature review](literature_review.md) | Direct source evidence, correlation discrepancy and source inconsistencies |
| [Model independent review](independent_review.md) | Separate recomputation of 1,352 model rows and displayed model-summary values |

The phase-conditioned build_summary.py generates the current report and px/s comparison from saved phase results; scripts/publish_current_summary.py synchronizes the main executive summary with corrected links. Root write_summary.py generates the historical conditional model report and then republishes the saved current summary. Model predictions keep their actual m/s units. Detailed reviews distinguish algorithm verification from physical validation.
