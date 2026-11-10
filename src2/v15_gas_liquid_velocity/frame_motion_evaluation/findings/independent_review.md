# Independent verification of frame-motion feasibility

## Executive result

**The saved motion calculations and declared support rules reproduce.** The retained 26 camera sequences each have 501 contiguous native frame identifiers; MII and XSV capture rates agree. All 13,000 archived adjacent-pair pixel-rate conversions and their preferred-method summaries reproduce exactly. The new temporal and local-patch outputs pass independent reaggregation.

This establishes implementation consistency, not experimentally accurate gas or liquid material velocity. The strongest wording supported by these checks is **observed image-structure displacement in px/frame and px/s**, with selection and consistency limitations.

| Checked output | Evidence | Result |
|---|---|---|
| Native scope and timing | Original markers, all frame filenames, native MII/XSV configs | 26 retained / 5 excluded; 501 frames each; 18 at 7,000 fps, 6 at 5,000 fps, 2 at 4,000 fps |
| Archived preferred adjacent estimates | Every stored pair and every per-run median/quartile | 13,000 pairs; 7,528 pass; 17 recordings meet coverage and saturation rules; zero arithmetic discrepancy |
| New temporal diagnostics | Every displacement/rate, adjacent-path mean, joint acceptance, tolerance decision, reversal sum and run summary | 38,896 comparisons across lags 1/2/4; 7 recordings pass the declared lag diagnostics |
| New local 2D patch experiment | Every stored conversion/round trip, patch acceptance, pair aggregation and run summary | 15,600 tracks, 2,600 method/pair summaries, 52 method/run summaries; declared results reproduce |

The independent verifier checks 307,619 numeric values in the three analyses. Maximum floating-point discrepancy is approximately 1.42 × 10⁻¹⁴. It imports neither analysis script. Its reaggregation checks declared decision rules, not the correctness of real-image feature identification.

## Interpretation cautions that must remain in the main findings

1. **A passing temporal diagnostic does not prove a unique representative speed.** Seven recordings pass: m04q01_001, m05q01_001, m05q05_001, m06q01_001, m06q05_001, m07q01_001 and m07q05_001. Two have conspicuous region/selection sensitivity in the earlier independently checked ROI/lag study: m04q01_001 has a 36.37% span across supported variants relative to its baseline magnitude, and m06q05_001 has a 19.41% span. Their lag-4 jointly accepted medians also differ from their full accepted adjacent medians: −7.362 → −5.067 px/frame and −15.986 → −18.996 px/frame, respectively. Those medians describe different accepted subsets; this is selection sensitivity, not an inconsistency in the pairwise arithmetic.
2. **The local-patch gate asks a different question from the global profile gate.** Requiring at least three accepted patches across both axial windows and multiple radial windows tests spatially distributed coherent translation. Only m07q05_001 meets that recording-level gate in either local method. Sparse bubbles, varying interface motion and featureless areas can fail that condition despite trackable features elsewhere. The low count must not be described as proof that all other motion estimates are false.
3. **Annular motion remains unresolved by these tested rules.** No annular run meets the local recording-level gate. Background subtraction yields some accepted annular patches and 10 coherent sampled pairs across the nine annular recordings, but this is not sufficient for a representative speed under the declared criteria. Search windows, patch texture, saturation and real deformation constrain this experiment.
4. **The sequence duration is short.** A 501-frame recording spans 500 intervals: 0.07143, 0.100 or 0.125 s at the three native capture rates. The 500 adjacent comparisons overlap; they are not 500 independent experimental repeats.
5. **Separate material phase velocities have not been identified.** Interfaces are shared by liquid and vapor. A gas/liquid mask split cannot independently establish each material velocity from the same moving boundary. Pixel units avoid the need for an axial length calibration, but do not resolve phase identity. Cross-recording px/s comparisons remain conditional on consistent magnification.

## Existing ROI/lag sensitivity

The archived ncc_dynamic sensitivity table contains nine variants per run (three image-row intervals × three frame lags), with 100 sampled pairs per variant. Among the 17 baseline-supported recordings, only five have at least 50 accepted pairs in all nine variants. The two large supported-variant median spans are reported above. Different variants accept different frame pairs, so this is a descriptive sensitivity comparison, not a controlled estimate of algorithm error or an uncertainty interval.

Details: [per-run existing sensitivity](../results/independent_existing_sensitivity.json).

## Reproduction

Run after the temporal and local-patch scripts:

    python3 -B src2/v15_gas_liquid_velocity/frame_motion_evaluation/scripts/independent_verify.py

Verification data: [independent_verification.json](../results/independent_verification.json). Original frame inventory and archived baseline arithmetic are also retained in [independent_baseline_verification.json](../results/independent_baseline_verification.json).

The previously saved analysis scripts and source recordings were not edited by this reviewer.
