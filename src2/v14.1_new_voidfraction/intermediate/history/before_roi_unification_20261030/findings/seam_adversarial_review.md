# Adversarial review of the v14.1 correction

**PASS.**

- synthetic_and_mutation_controls: PASS
- method_preservation: PASS
- source_column_order_and_immutable_inputs: PASS
- all_thresholds_and_hulls: PASS
- honest_comparison_and_reports: PASS

The tests reject the old reordered ramp and a deliberate one-pixel mutation, verify every native strip against the corrected output, and restrict real input changes to the defective seam columns. Estimator, classifier, geometry and selection preservation is checked against the frozen v14 source. All 130 threshold cases are recomputed by independent boolean segmentation/interval arithmetic; ten before/after hull comparisons use both a different linear-program formulation and direct polygon intersection without an optimizer. Classifier degradations remain in the comparison. Original files are hash-verified.

The correction does not establish better physical accuracy or remove sensitivity to segmentation and reconstruction assumptions. [Detailed evidence](../results/validation/seam_adversarial_checks.json) · [Effectiveness](fix_effectiveness.md)
