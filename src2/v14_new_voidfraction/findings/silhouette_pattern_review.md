# Silhouette features: internal pattern-identification comparison

This comparison tests whether the panorama silhouette improves identification within the current user-selected scope. It does not validate conditional void values as true gas volume. All 25 selected, acquisitions passing the saturation criterion are evaluated: dispersed bubbles: 1; plug: 3; slug: 3; churn: 10; wavy annular: 8. Excluded acquisitions do not enter features, references, training, tests or figures.

Seven features are fixed before modeling: conditional visible circular area mean, center candidate-gas duty, candidate-mask pixel flip-rate over rows 18:71, center-row flip-rate, visible radial-width standard deviation /53, mean center gas-interval length /53, and longest connected center gas interval divided by record length. A candidate interval touching either axial record edge is explicitly counted as censored. Interval summaries include those visible censored lengths and are not complete bubble sizes.

## Leakage prevention

For every held-out run or mass-flux group, rebuild the liquid-like template from only the fixed reference candidates remaining among that fold's training acquisitions. Feature scaling, class centroids and nearest neighbors likewise use training acquisitions only. The same fold-specific template segments both training and test images. Every fold verifies that the held-out acquisition is absent from its reference. 4 different templates are used across the folds.

Reference candidates were previously chosen from visually reviewed liquid-containing panoramas. The detector and its settings were developed using this dataset. Removing a held-out acquisition from template construction and classifier fitting prevents direct fitting leakage but cannot undo that prior inspection. These scores are internal cross-validation conditional on the current detector, not untouched external validation or an independent liquid-only calibration.

Models are the same fixed nearest-centroid, k=1 and k=3 nearest-neighbor methods as the earlier morphology comparison. No hyperparameters or feature subsets are optimized, and every model is reported.

| Model | Validation | Earlier morphology accuracy / balanced | Silhouette accuracy / balanced |
|---|---|---:|---:|
| centroid | leave_one_run_out | 0.720 / 0.623 | 0.520 / 0.425 |
| knn1 | leave_one_run_out | 0.800 / 0.575 | 0.720 / 0.530 |
| knn3 | leave_one_run_out | 0.760 / 0.550 | 0.680 / 0.458 |
| centroid | leave_one_mass_flux_out | 0.680 / 0.598 | 0.520 / 0.425 |
| knn1 | leave_one_mass_flux_out | 0.880 / 0.708 | 0.720 / 0.530 |
| knn3 | leave_one_mass_flux_out | 0.720 / 0.535 | 0.680 / 0.552 |

## Recall and all mistakes

**centroid, leave_one_run_out**: dispersed bubbles: 0.000; plug: 0.667; slug: 0.333; churn: 0.500; wavy annular: 0.625.
Errors: m03q05_001 (churn→wavy annular); m04q01_001 (slug→plug); m04q20_001 (churn→wavy annular); m05q05_001 (slug→plug); m05q40_001 (churn→wavy annular); m06q01_001 (plug→slug); m06q10_001 (churn→slug); m06q40_001 (wavy annular→churn); m06q80_001 (wavy annular→churn); m07q01_001 (dispersed bubbles→plug); m07q20_001 (churn→slug); m07q80_001 (wavy annular→churn).

**knn1, leave_one_run_out**: dispersed bubbles: 0.000; plug: 0.667; slug: 0.333; churn: 0.900; wavy annular: 0.750.
Errors: m03q05_001 (churn→wavy annular); m05q05_001 (slug→plug); m05q80_001 (wavy annular→churn); m06q05_001 (slug→churn); m06q80_001 (wavy annular→churn); m07q01_001 (dispersed bubbles→plug); m07q05_001 (plug→slug).

**knn3, leave_one_run_out**: dispersed bubbles: 0.000; plug: 0.333; slug: 0.333; churn: 1.000; wavy annular: 0.625.
Errors: m05q05_001 (slug→plug); m05q80_001 (wavy annular→churn); m06q01_001 (plug→slug); m06q05_001 (slug→churn); m06q40_001 (wavy annular→churn); m06q80_001 (wavy annular→churn); m07q01_001 (dispersed bubbles→plug); m07q05_001 (plug→slug).

**centroid, leave_one_mass_flux_out**: dispersed bubbles: 0.000; plug: 0.667; slug: 0.333; churn: 0.500; wavy annular: 0.625.
Errors: m03q05_001 (churn→wavy annular); m04q01_001 (slug→plug); m04q10_001 (churn→wavy annular); m04q20_001 (churn→wavy annular); m05q05_001 (slug→plug); m06q01_001 (plug→slug); m06q10_001 (churn→slug); m06q40_001 (wavy annular→churn); m06q80_001 (wavy annular→churn); m07q01_001 (dispersed bubbles→plug); m07q20_001 (churn→slug); m07q80_001 (wavy annular→churn).

**knn1, leave_one_mass_flux_out**: dispersed bubbles: 0.000; plug: 0.667; slug: 0.333; churn: 0.900; wavy annular: 0.750.
Errors: m03q05_001 (churn→wavy annular); m05q05_001 (slug→plug); m05q80_001 (wavy annular→churn); m06q05_001 (slug→churn); m06q80_001 (wavy annular→churn); m07q01_001 (dispersed bubbles→slug); m07q05_001 (plug→slug).

**knn3, leave_one_mass_flux_out**: dispersed bubbles: 0.000; plug: 1.000; slug: 0.333; churn: 0.800; wavy annular: 0.625.
Errors: m04q10_001 (churn→wavy annular); m05q05_001 (slug→plug); m05q80_001 (wavy annular→churn); m06q05_001 (slug→churn); m06q40_001 (wavy annular→churn); m06q80_001 (wavy annular→churn); m07q01_001 (dispersed bubbles→plug); m07q20_001 (churn→slug).

## Interpretation

For the 6 fixed model/validation pairs, silhouette balanced accuracy is higher in 1, equal in 0, and lower in 5 than the morphology comparison on the same selected acquisitions. Nearest-centroid leave-one-run-out balanced accuracy is 0.623 for morphology and 0.425 for silhouettes; leave-one-mass-flux-out is 0.598 and 0.425. These observed results do not establish performance for future classifiers or experiments.

A score improvement would show that silhouette information helps match these supplied labels under this validation scheme. It would not establish correct absolute void fraction. The user confirmed all 53 rows in 18:71 as the complete tube interior; D=53 and a fully filled lumen gives α=1. Near-unity estimates may still include unresolved liquid films or optical shadows. No extra fluid area beyond these confirmed boundaries is inferred.

The candidate-mask flip denominator is 53*(W−1), whereas the legacy dark-pixel baseline uses 40*(W−1). Fixed-row y=40 flips use W−1; legacy center_* field names denote this row, whereas the confirmed lumen midpoint is y=44. Raw switch numerators therefore have different denominators. All features derive from a spatial stitched record, without a conversion to seconds.

Only one dispersed-bubble acquisition remains. Its held-out fold has no bubble training example; absent classes have no fitted centroid and cannot be predicted. Every test remains in the score, including unsupported-class tests. Folds missing at least one training class: leave-one-run-out 1, leave-one-mass-flux-out 5. Per-fold availability is saved. Plug and slug have three acquisitions each, so one error changes their recall by one third. No labels are changed to improve a score.

Checks cover exact gas-interval length/censoring examples, the fully filled 53-row lumen with α=1, zero flip-rate for a constant mask, and every fold's test-reference exclusion. Segmentation parameters are inherited unchanged from `panorama_void_alternatives.py`; its source hash is saved.

Reproduce with `python3 src2/v14_new_voidfraction/scripts/silhouette_pattern_identification.py`. Outputs: `silhouette_pattern_features.csv`, `silhouette_pattern_predictions.csv`, `silhouette_pattern_fold_references.csv`, `silhouette_pattern_scores.csv`, `silhouette_pattern_results.json`, and `silhouette_pattern_confusion_matrices.png`.
