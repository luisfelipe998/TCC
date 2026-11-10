# Silhouette features: independent pattern-identification comparison

This comparison tests whether the new panorama silhouette improves identification. It does not validate its conditional void values as true gas volume. All 26 quality-filtered pure acquisitions are retained, including both dispersed-bubbles labels.

Seven features are fixed before modeling: conditional visible circular area mean, center candidate-gas duty, candidate-mask pixel flip-rate over rows 18:71, center-row flip-rate, visible radial-width standard deviation /60, mean center gas-interval length /60, and longest connected center gas interval divided by record length. A candidate interval touching either axial record edge is explicitly counted as censored. Interval summaries include those visible censored lengths and are not complete bubble sizes.

## Leakage prevention

For every held-out run or mass-flux group, rebuild the liquid-like template from only the fixed reference candidates remaining among that fold's training acquisitions. Feature scaling, class centroids and nearest neighbors likewise use training acquisitions only. The same fold-specific template segments both training and test images. Every fold verifies that the held-out acquisition is absent from its reference. Four different templates are needed across the folds.

Reference candidates were previously chosen from visually reviewed liquid-containing panoramas. That prior human choice is part of the method, so this is not a fully automatic discovery of calibration images. Removing the held-out acquisition prevents direct test-image leakage; it does not create an independent liquid-only calibration.

Models are the same fixed nearest-centroid, k=1 and k=3 nearest-neighbor methods as the earlier morphology comparison. No hyperparameters or feature subsets are optimized, and every model is reported.

| Model | Validation | Earlier morphology accuracy / balanced | Silhouette accuracy / balanced |
|---|---|---:|---:|
| centroid | leave_one_run_out | 0.692 / 0.623 | 0.615 / 0.532 |
| knn1 | leave_one_run_out | 0.731 / 0.508 | 0.577 / 0.330 |
| knn3 | leave_one_run_out | 0.731 / 0.550 | 0.692 / 0.442 |
| centroid | leave_one_mass_flux_out | 0.654 / 0.598 | 0.577 / 0.512 |
| knn1 | leave_one_mass_flux_out | 0.808 / 0.642 | 0.615 / 0.397 |
| knn3 | leave_one_mass_flux_out | 0.692 / 0.535 | 0.692 / 0.525 |

## Recall and all mistakes

**centroid, leave_one_run_out**: dispersed bubbles: 0.000; plug: 0.667; slug: 0.667; churn: 0.700; wavy annular: 0.625.
Errors: m03q05_001 (churn→wavy annular); m05q05_001 (slug→plug); m05q05_002 (dispersed bubbles→slug); m06q01_001 (plug→slug); m06q10_001 (churn→slug); m06q40_001 (wavy annular→churn); m06q80_001 (wavy annular→churn); m07q01_001 (dispersed bubbles→plug); m07q20_001 (churn→slug); m07q80_001 (wavy annular→churn).

**knn1, leave_one_run_out**: dispersed bubbles: 0.000; plug: 0.000; slug: 0.000; churn: 0.900; wavy annular: 0.750.
Errors: m03q05_001 (churn→wavy annular); m04q01_001 (slug→churn); m05q01_001 (plug→slug); m05q05_001 (slug→plug); m05q05_002 (dispersed bubbles→slug); m05q80_001 (wavy annular→churn); m06q01_001 (plug→slug); m06q05_001 (slug→plug); m06q80_001 (wavy annular→churn); m07q01_001 (dispersed bubbles→plug); m07q05_001 (plug→slug).

**knn3, leave_one_run_out**: dispersed bubbles: 0.000; plug: 0.333; slug: 0.000; churn: 1.000; wavy annular: 0.875.
Errors: m04q01_001 (slug→churn); m05q05_001 (slug→plug); m05q05_002 (dispersed bubbles→plug); m05q80_001 (wavy annular→churn); m06q01_001 (plug→slug); m06q05_001 (slug→churn); m07q01_001 (dispersed bubbles→plug); m07q05_001 (plug→slug).

**centroid, leave_one_mass_flux_out**: dispersed bubbles: 0.000; plug: 0.667; slug: 0.667; churn: 0.600; wavy annular: 0.625.
Errors: m03q05_001 (churn→wavy annular); m04q20_001 (churn→wavy annular); m05q05_001 (slug→plug); m05q05_002 (dispersed bubbles→plug); m06q01_001 (plug→slug); m06q10_001 (churn→slug); m06q40_001 (wavy annular→churn); m06q80_001 (wavy annular→churn); m07q01_001 (dispersed bubbles→plug); m07q20_001 (churn→slug); m07q80_001 (wavy annular→churn).

**knn1, leave_one_mass_flux_out**: dispersed bubbles: 0.000; plug: 0.333; slug: 0.000; churn: 0.900; wavy annular: 0.750.
Errors: m03q05_001 (churn→wavy annular); m04q01_001 (slug→churn); m05q05_001 (slug→plug); m05q05_002 (dispersed bubbles→plug); m05q80_001 (wavy annular→churn); m06q01_001 (plug→slug); m06q05_001 (slug→churn); m06q80_001 (wavy annular→churn); m07q01_001 (dispersed bubbles→slug); m07q05_001 (plug→slug).

**knn3, leave_one_mass_flux_out**: dispersed bubbles: 0.000; plug: 1.000; slug: 0.000; churn: 1.000; wavy annular: 0.625.
Errors: m04q01_001 (slug→churn); m05q05_001 (slug→plug); m05q05_002 (dispersed bubbles→plug); m05q80_001 (wavy annular→churn); m06q05_001 (slug→churn); m06q40_001 (wavy annular→churn); m06q80_001 (wavy annular→churn); m07q01_001 (dispersed bubbles→plug).

## Interpretation

Under these fixed features and models, the silhouette does not improve five-pattern identification: all six balanced-accuracy scores are lower than the earlier morphology scores. The nearest-centroid leave-one-run-out balanced accuracy falls from 0.623 to 0.532, and leave-one-mass-flux-out from 0.598 to 0.512. This result applies to the tested comparison, not to every possible future classifier or segmentation model.

A score improvement would show that silhouette information helps match these supplied labels under this validation scheme. It would not establish correct absolute void fraction. The silhouette crop contains 53 rows but D_app=60: a completely filled visible crop gives (53/60)²=0.7802778. Near-ceiling annular estimates are radially censored; width, duty and continuity features can recognize that geometry without resolving the missing gas area.

The candidate-mask flip denominator is 53*(W−1), whereas the legacy dark-pixel baseline uses 40*(W−1). Center-row flips use W−1. Neither can be compared by treating their raw switch numerators as equal quantities. All features derive from a spatial stitched record, without a conversion to seconds.

The bubble class has only two visibly different acquisitions, and plug/slug have three each. A single bubble error moves recall by one half; a plug or slug error moves it by one third. Disclosed per-run mistakes and confusion matrices matter more than a headline accuracy. No labels are changed to improve a score.

Checks cover exact gas-interval length/censoring examples, the fully filled 53-row geometric ceiling, zero flip-rate for a constant mask, and every fold's test-reference exclusion. Segmentation parameters are inherited unchanged from `panorama_void_alternatives.py`; its source hash is saved.

Reproduce with `python3 src2/v14_new_voidfraction/scripts/silhouette_pattern_identification.py`. Outputs: `silhouette_pattern_features.csv`, `silhouette_pattern_predictions.csv`, `silhouette_pattern_fold_references.csv`, `silhouette_pattern_scores.csv`, `silhouette_pattern_results.json`, and `silhouette_pattern_confusion_matrices.png`.
