# Numerical results

[Back to v14](../README.md)

Use the per-run tables for individual acquisitions. Group summaries use the supplied labels and their documented saturation filter. Conditional visible area is not calibrated absolute gas volume.

| Folder | Main tables |
|---|---|
| [void_fraction/](void_fraction) | [Per-run void estimates](void_fraction/panorama_void_per_run.csv), [grouped estimates](void_fraction/panorama_void_by_pattern.csv) |
| [Chart data](void_fraction/void_fraction_flip_plot_data.csv) | Exact per-run coordinates, labels and exclusions used in the pattern charts |
| [flip_rate/](flip_rate) | [Current spatial silhouette flips](flip_rate/phase_mask_flip_per_run.csv), [grouped silhouette flips](flip_rate/phase_mask_flip_by_pattern.csv), [legacy optical flips](flip_rate/flip_rate_per_run.csv), [legacy groups](flip_rate/flip_rate_by_pattern.csv) |
| [patterns/](patterns) | [Morphology predictions](patterns/pattern_identification_predictions.csv), [scores](patterns/pattern_identification_scores.csv), [silhouette predictions](patterns/silhouette_pattern_predictions.csv), [scores](patterns/silhouette_pattern_scores.csv) |
| [validation/](validation) | [Selection](validation/dataset_selection.json), [method checks](validation/deep_method_checks.json), [method changes](validation/deep_method_changes.csv), [D53 checks](validation/lumen_normalization_checks.json), [D53 changes](validation/lumen_normalization_changes.csv), [final verification](validation/final_verification.json), [current adversarial checks](validation/adversarial_checks.json), [post-summary adversarial audit](validation/post_regeneration_adversarial_checks.json), [input inventory](validation/data_inventory.csv), [liquid holdout controls](validation/panorama_liquid_holdout_checks.csv), [chart checks](validation/summary_chart_checks.json), [scope and navigation check](validation/layout_verification.json) |

Read the [flip definitions](../findings/phase_mask_flip_review.md) before using flip-rate columns: the gas-silhouette rate uses 53 radial rows and the legacy optical panorama rate uses 40. Sensitivity sweeps live in [intermediate/sensitivity/](../intermediate/sensitivity).
