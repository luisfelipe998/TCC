# Phase-conditioned image motion

Start with the [EXECUTIVE SUMMARY](findings/EXECUTIVE_SUMMARY.md). This evaluates separate motion averages during gas-like and liquid-looking passages, using current v14 masks applied to native frames. It preserves unknown states and distinguishes conditional image motion from gas/liquid material velocity.

The original 26-run selection and source data are unchanged. Results are local to this directory. Packed native masks are retained under intermediate/masks; unpacking instructions and source hashes are in results/native_phase_metadata.json.

## Reproduction order

Run the following scripts from this directory using python3 -B, in this order (their paths and input roots are resolved from the script location):

    python3 -B scripts/native_phase_masks.py --workers 3
    python3 -B scripts/conditioned_motion.py
    python3 -B scripts/passage_sensitivity.py
    python3 -B scripts/phase_visual_audit.py
    python3 -B scripts/closed_contour_probe.py
    python3 -B scripts/closed_contour_bubble_tracks.py
    python3 -B scripts/closed_contour_report.py
    python3 -B scripts/adaptive_phase_patches.py
    python3 -B scripts/diagnostic_figures.py
    python3 -B scripts/mask_complement_control.py
    python3 -B scripts/independent_verify.py
    python3 -B scripts/verify_adaptive.py
    python3 -B scripts/build_summary.py

The scripts use the existing NumPy/SciPy/Pillow/matplotlib environment. The original global/local motion evidence from frame_motion_evaluation remains an input and is not regenerated here. The independent verifier imports neither conditional estimator and separately reconstructs relevant masks/decisions.

## Main numerical outputs

- passage_frames.csv / passage_pairs.csv / passage_episodes.csv: every tag, uncertainty state, selected comparison and episode.
- passage_phase_summary.csv / passage_ratios.csv: signed means, counts, episode weighting and explicit ratio-withholding reasons.
- phase_patch_tracks.csv / phase_patch_summary.csv / phase_patch_ratios.csv: actual matched-region phase tests at three purity rules.
- passage_sensitivity_ratios.csv: probe location and relaxed liquid-rich definitions, including failed ratios.
- conditional_ratio_by_pattern.csv: per-pattern summaries with actual contributing run counts; these are not physical slip.
- conditional_passage_velocities.csv: exact signed px/s means, ratios and counts used in the headline chart.
- closed_contour_*: dispersed detector repair, visual evidence, candidate bubble tracks and checks.
- adaptive_*: smaller patches selected within phase-interior candidates.

Undefined means/ratios remain NaN with recorded reasons. A missing phase estimate is not zero velocity. Sign is retained; positive x is rightward. No frame-independent confidence interval or accuracy percentage is asserted.
