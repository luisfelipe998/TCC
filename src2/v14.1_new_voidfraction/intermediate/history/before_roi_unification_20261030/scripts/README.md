# Reproduction stages

[Back to v14.1](../README.md)

1. [reconstruct_inputs.py](reconstruct_inputs.py) — Reconstruct original panoramas as a control, then preserve source-strip order.
2. [panorama_void_alternatives.py](panorama_void_alternatives.py) — Infer liquid references and reconstruct conditional visible gas area.
3. [phase_mask_flip_analysis.py](phase_mask_flip_analysis.py) — Save gas masks and calculate spatial gas-silhouette flips.
4. [raw_flip_analysis.py](raw_flip_analysis.py) — Reproduce legacy optical flips and supplemental native-frame checks.
5. [spatial_flip_sensitivity.py](spatial_flip_sensitivity.py) — Perturb spatial optical thresholds and radial crops.
6. [inventory_and_legacy_audit.py](inventory_and_legacy_audit.py) — Record input hashes and independently audit legacy values.
7. [pattern_identification.py](pattern_identification.py) — Compare fixed morphology classifiers and create the contact sheet.
8. [silhouette_pattern_identification.py](silhouette_pattern_identification.py) — Compare silhouette classifiers with training-fold references.
9. [adversarial_review.py](adversarial_review.py) — Challenge masks with phase anchors, reference choices and crops.
10. [deep_method_checks.py](deep_method_checks.py) — Independently verify masks and equations; compare v14 and v14.1 at D53.
11. [lumen_normalization_checks.py](lumen_normalization_checks.py) — Verify the unchanged D53 geometry on current masks.
12. [image_void_alternatives.py](image_void_alternatives.py) — Reproduce supplemental native-frame alternatives.
13. [final_verification.py](final_verification.py) — Independently verify saved arithmetic, input hashes and classifier scores.
14. [summary_figures.py](summary_figures.py) — Plot void versus spatial flip-rate and crop sensitivity from saved results.
15. [compare_versions.py](compare_versions.py) — Compare every acquisition, all threshold settings and all classifiers.
16. [write_summary.py](write_summary.py) — Generate the executive summary and folder navigation.
17. [seam_adversarial.py](seam_adversarial.py) — Challenge source-column order, preserved methods, thresholds and change claims.
18. [post_regeneration_adversarial.py](post_regeneration_adversarial.py) — Independently challenge the completed analysis, plots and findings.
19. [verify_layout.py](verify_layout.py) — Verify current exclusions, retained measurements and report navigation.

The ordered runner is [run_all.py](run_all.py). Numerical estimator and classifier functions are checked against the frozen baseline by [seam_adversarial.py](seam_adversarial.py).
