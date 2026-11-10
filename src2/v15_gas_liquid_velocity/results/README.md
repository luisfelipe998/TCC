# Result tables

The [current executive summary](../findings/EXECUTIVE_SUMMARY.md) reports phase-conditioned image motion in **px/s** and dimensionless ratios. Its results are stored under phase_conditioned_motion/results:

| Current file | Interpretation |
|---|---|
| [conditional_passage_velocities.csv](../phase_conditioned_motion/results/conditional_passage_velocities.csv) | Exact values plotted in the current comparison: four eligible reference-guard, lag-2-screened recordings |
| [passage_ratios.csv](../phase_conditioned_motion/results/passage_ratios.csv) | Every configuration's signed means in px/frame, ratio decisions, accepted counts and withholding reasons |
| [passage_phase_summary.csv](../phase_conditioned_motion/results/passage_phase_summary.csv) | Phase-candidate motion summaries and native FPS for px/s conversion |
| [conditional_ratio_by_pattern.csv](../phase_conditioned_motion/results/conditional_ratio_by_pattern.csv) | Selected descriptive ratios with actual contributing counts; no universal physical pattern ratio |
| [phase_patch_ratios.csv](../phase_conditioned_motion/results/phase_patch_ratios.csv) | Fixed-patch phase-interior tests and withheld comparisons |
| [adaptive_ratios.csv](../phase_conditioned_motion/results/adaptive_ratios.csv) | Smaller/adaptive phase-interior tests; no recording meets both-phase reporting criteria |

Additional per-frame observations, probe sensitivity, masks, dispersed-object and adaptive-patch outputs are indexed in the [phase-conditioned README](../phase_conditioned_motion/README.md#main-numerical-outputs). The preceding global/local displacement results are under [frame_motion_evaluation](../frame_motion_evaluation/README.md).

The root tables below preserve the earlier model and optical investigation. **All m/s values are conditional thermodynamic predictions**, documented in the [historical model summary](../findings/conditional_model_summary.md); they are not conversions of the latest image-motion results. Optical velocities remain in px/frame or px/s.

| File | Interpretation |
|---|---|
| scenario_inventory.csv | All 31 folder identities, native frame rates, hashes and five exclusion reasons |
| conditional_velocities_per_run.csv | 26 nominal-G, saturated-inlet conditions × four closures |
| recorded_flow_conditional_velocities_per_run.csv | 26 recorded-mean-flow conditions × four closures; camera pairing unverified |
| conditional_velocities_by_pattern.csv | Nominal min/max/mean/median over the sampled conditions |
| recorded_flow_conditional_velocities_by_pattern.csv | Recorded-flow min/max/mean/median |
| thermodynamic_sensitivity.csv | 1,352 rows: 26 conditions × four closures × thirteen input scenarios |
| sensitivity_relative_to_baseline.csv | All input alternatives compared with the nominal RA baseline |
| closure_spread.csv | Same-condition nominal phase-velocity spread across closures; not uncertainty bounds |
| massflow_per_run.csv | Recorded flow aggregates, nominal comparison, hashes and mapping limitations |
| massflow_source_blocks.csv | Original sensor-block statistics and provenance |
| image_alpha_diagnostics.csv | Unvalidated hybrid diagnostic using v14 spatial alpha; saturated inputs withheld |
| optical_preferred_candidate_per_run.csv | Preferred background-subtracted NCC candidate, native timing and validity flags |
| optical_preferred_candidate_pairs.csv | All 13,000 preferred-candidate pairs, including rejected pairs and quality metrics |
| optical_full_method_summary.csv | All four optical methods over all native adjacent pairs |
| optical_primary_pairs.csv | Individual pairs of the high-pass reference experiment |
| optical_per_run.csv | Per-run high-pass reference experiment, not the preferred candidate |
| optical_sensitivity.csv | 936 run/method/crop/lag combinations, 100 sampled pairs each |
| optical_synthetic_trials.csv | Known translations for estimator verification, not simulated phase-flow truth |
| optical_method_agreement.csv | Differences on jointly accepted pairs; selection-biased diagnostic |
| optical_metadata.json / optical_validation.json | Settings, excluded runs, input fingerprints and verification |
| validation | Legacy/mass/model audits, historical model-summary values and final provenance checks |

Model units: G in kg/(m² s), q in kW/m², x/alpha/slip dimensionless. Undefined or rejected values remain empty/NaN, not zero. quality_status, phase-candidate eligibility and image/optical validity fields explain missing results in their respective tables. Excluded camera acquisitions occur only in inventories/audit metadata, never in compared velocity rows.
