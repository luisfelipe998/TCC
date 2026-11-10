# Figures

The current [executive summary](../findings/EXECUTIVE_SUMMARY.md) uses the phase-passage figures below. Velocity values in its main comparison are px/s; ratios are dimensionless.

- [Current passage means and ratios](../phase_conditioned_motion/figures/conditional_passage_velocities.png) ([SVG](../phase_conditioned_motion/figures/conditional_passage_velocities.svg)): four reference-guard, lag-2-screened recordings; signed mean image motion during gas-like passages and liquid-looking gaps. [Plotted values](../phase_conditioned_motion/results/conditional_passage_velocities.csv).
- [Passage diagnostics](../phase_conditioned_motion/figures/passage_time_diagnostics.png): phase-candidate occupancy and global image motion over time; the displayed motion scale is 1000 px/s. This diagnostic uses the consensus adjacent configuration, not the reference-guard headline configuration.
- [Method-dependent ratios](../phase_conditioned_motion/figures/conditional_ratio_comparison.png): dimensionless comparisons across temporal/local motion support.
- [Fixed-patch support](../phase_conditioned_motion/figures/phase_support_attrition.png) and [adaptive support](../phase_conditioned_motion/figures/adaptive_support.png): accepted correspondence counts under stricter phase rules.
- [Liquid-candidate examples](../phase_conditioned_motion/figures/adaptive_liquid_examples.png): visible limits of inferred phase purity.
- [Candidate bubble trajectories](../phase_conditioned_motion/figures/closed_contour_bubble_tracks.png): dispersed-flow object tracking after the tested detector correction.

The earlier thermodynamic plots are preserved with their real **m/s prediction units** and discussed in the [historical conditional model summary](../findings/conditional_model_summary.md). They are not the current image-motion result or a physical-unit conversion of it:

- [Recorded-flow model comparison](recorded_flow_comparison.png): conditional superficial and RA phase predictions; connected hollow nominal-G and filled recorded-flow points.
- [Superficial model velocity](superficial_velocity.png): nominal saturated-inlet baseline, with all 26 selected conditions.
- [Phase model comparison](phase_model_comparison.png): homogeneous, RA and both XF forms on the same nominal inputs.
- [Model inlet sensitivity](inlet_sensitivity.png): assumed input perturbations; undefined two-phase cases are marked n/a.

Earlier optical diagnostics:

- [Optical displacement quality](optical_displacement_quality.png): high-pass reference experiment, accepted-pair medians/IQR and coverage. This is not the preferred optical candidate.
- [Optical method comparison](optical_method_comparison.png): preprocessing and estimator dependence, without physical ground truth.
- [Native correspondence](optical_native_correspondence.png): visual leftward motion and profile alignment; the illustrated pair fails the stated competitor-margin gate and is not an accuracy benchmark.

Figures show selected observations, sample/model spread and feasibility checks, not calibrated uncertainty or proven flow-pattern boundaries. The optical phase labels do not establish separate material velocities. SVG copies accompany the current comparison and principal thermodynamic plots.
