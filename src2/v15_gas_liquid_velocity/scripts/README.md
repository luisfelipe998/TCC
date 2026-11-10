# Scripts and reproduction

To regenerate the current px/s comparison and executive summary from saved phase results, run from the repository root:

    python3 -B src2/v15_gas_liquid_velocity/phase_conditioned_motion/scripts/build_summary.py

This builds the current chart/values and invokes publish_current_summary.py to synchronize the main report. Running scripts/publish_current_summary.py alone only synchronizes the existing nested summary and its links. It does not recompute motion.

To rerun later experiments from images, follow these complete sequences in order:

1. [Native-frame motion reproduction](../frame_motion_evaluation/README.md#reproduce): temporal motion, local patches, visual evidence and independent checks.
2. [Phase-conditioned reproduction](../phase_conditioned_motion/README.md#reproduction-order): native v14 candidates, conditional means, sensitivity, dispersed repair/tracking, adaptive patches, verification and current summary.

The original model/optical pipeline remains available from any directory:

    python3 -B /path/to/TCC/src2/v15_gas_liquid_velocity/scripts/run_all.py

1. audit_legacy.py: import and independently reproduce original v12/v12.1 values; verify stored tables and unchanged source hashes.
2. massflow_audit.py: check owner-documented flow against every original LabVIEW source block; generate aggregate and provenance tables.
3. model_comparison.py: compute nominal/recorded flow, four closures, inlet/property/heat scenarios, hybrid diagnostic, plots and model checks.
4. optical_investigation.py: native-frame signed correlation methods, crop/lag alternatives, synthetic tests, saved-result verification and optical report.
5. literature/source_numeric_checks.py: independently calculate the source consistency diagnostics.
6. write_summary.py: produce findings/conditional_model_summary.md and the detailed conditional model tables in m/s, then republish the saved phase-conditioned executive summary.
7. independent_model_review.py: reconstruct all saved thermodynamic rows and historical model summaries, including a comparison to the pinned secondary XF function.
8. final_verification.py: scope, headline values, provenance, source preservation, JSON and relative-link checks; output manifest.

run_all.py does not rerun the later frame_motion_evaluation or phase_conditioned_motion experiments; its current-summary publication uses their saved report. It stops at the first failure, retains each stage log and records command/runtime details. It uses the existing interpreter and does not modify original inputs. The plotting cache is under intermediate. No network or fresh source retrieval is part of reproduction.
