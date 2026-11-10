# Image-only pattern identification: feasibility audit

Morphology can help identify regimes without calibrating absolute void fraction. This audit extracts 12 optical/morphology features from each saved panorama. The extraction function receives no labels, mass flux, heat flux or camera timing. Folder labels are used only to train and validate simple classifiers.

The current user-selected scope contains 26 acquisitions. The existing saturation rule leaves 25 five-class tests: dispersed bubbles: 1; plug: 3; slug: 3; churn: 10; wavy annular: 8. Selection is applied before feature extraction and plotting; only selected runs appear in the feature CSV and contact sheet. Saturated selected images are excluded from classifier fitting and scores. Each acquisition is one replicate; axial crops do not increase sample size.

Features describe pixel flip density, dark projected area, radial dark-fraction variability, bright-band width mean/variability/absence across fixed image row y=40, radial contrast median/variability/flat-column fraction, and enclosed dark-contour hole area/size/density. The legacy center_* feature names refer to that fixed row, whose physical centering is unvalidated. Features are not validated gas areas. All masks use fixed per-panorama Otsu; contrast-flatness threshold 0.25 and 3×3 closing are fixed across runs.

## Validation results

Feature standardization and class centroids/neighbors are fit on training runs only. Models are a nearest class centroid and k-nearest neighbors at fixed k=1 and k=3; this script performs no automated parameter or feature selection. All tested models are reported. Balanced accuracy is mean per-class recall; five-class uniform-guess chance is 0.20. The majority-class ordinary-accuracy baseline is 0.400. Missing training classes have no fitted centroid and cannot be predicted; every held-out test is retained in the score.

These are internal cross-validation results conditional on the current feature definitions. The same dataset was inspected during method development; keeping a run out of a fitted classifier does not undo that prior inspection. An untouched acquisition set would be needed for an independent performance estimate.

| Model | Validation | Variant | Accuracy | Balanced accuracy |
|---|---|---|---:|---:|
| centroid | leave_one_run_out | baseline | 0.720 | 0.623 |
| centroid | leave_one_mass_flux_out | baseline | 0.680 | 0.598 |
| knn1 | leave_one_run_out | baseline | 0.800 | 0.575 |
| knn1 | leave_one_mass_flux_out | baseline | 0.880 | 0.708 |
| knn3 | leave_one_run_out | baseline | 0.760 | 0.550 |
| knn3 | leave_one_mass_flux_out | baseline | 0.720 | 0.535 |
| centroid | leave_one_run_out | otsu_minus10 | 0.680 | 0.598 |
| centroid | leave_one_run_out_consistent_variant | otsu_minus10 | 0.640 | 0.532 |
| knn1 | leave_one_run_out | otsu_minus10 | 0.680 | 0.463 |
| knn1 | leave_one_run_out_consistent_variant | otsu_minus10 | 0.800 | 0.575 |
| knn3 | leave_one_run_out | otsu_minus10 | 0.720 | 0.530 |
| knn3 | leave_one_run_out_consistent_variant | otsu_minus10 | 0.800 | 0.575 |
| centroid | leave_one_run_out | otsu_plus10 | 0.600 | 0.465 |
| centroid | leave_one_run_out_consistent_variant | otsu_plus10 | 0.680 | 0.557 |
| knn1 | leave_one_run_out | otsu_plus10 | 0.680 | 0.422 |
| knn1 | leave_one_run_out_consistent_variant | otsu_plus10 | 0.800 | 0.533 |
| knn3 | leave_one_run_out | otsu_plus10 | 0.640 | 0.443 |
| knn3 | leave_one_run_out_consistent_variant | otsu_plus10 | 0.760 | 0.508 |
| centroid | leave_one_run_out | radial_crop15_65 | 0.520 | 0.290 |
| centroid | leave_one_run_out_consistent_variant | radial_crop15_65 | 0.680 | 0.557 |
| knn1 | leave_one_run_out | radial_crop15_65 | 0.680 | 0.375 |
| knn1 | leave_one_run_out_consistent_variant | radial_crop15_65 | 0.760 | 0.508 |
| knn3 | leave_one_run_out | radial_crop15_65 | 0.640 | 0.402 |
| knn3 | leave_one_run_out_consistent_variant | radial_crop15_65 | 0.720 | 0.525 |
| centroid | leave_one_run_out | radial_crop20_60 | 0.160 | 0.080 |
| centroid | leave_one_run_out_consistent_variant | radial_crop20_60 | 0.720 | 0.623 |
| knn1 | leave_one_run_out | radial_crop20_60 | 0.040 | 0.020 |
| knn1 | leave_one_run_out_consistent_variant | radial_crop20_60 | 0.680 | 0.510 |
| knn3 | leave_one_run_out | radial_crop20_60 | 0.280 | 0.187 |
| knn3 | leave_one_run_out_consistent_variant | radial_crop20_60 | 0.720 | 0.535 |
| centroid | leave_one_run_out | first_axial_half | 0.680 | 0.557 |
| centroid | leave_one_run_out_consistent_variant | first_axial_half | 0.680 | 0.557 |
| knn1 | leave_one_run_out | first_axial_half | 0.760 | 0.555 |
| knn1 | leave_one_run_out_consistent_variant | first_axial_half | 0.720 | 0.488 |
| knn3 | leave_one_run_out | first_axial_half | 0.720 | 0.483 |
| knn3 | leave_one_run_out_consistent_variant | first_axial_half | 0.640 | 0.350 |
| centroid | leave_one_run_out | second_axial_half | 0.520 | 0.332 |
| centroid | leave_one_run_out_consistent_variant | second_axial_half | 0.600 | 0.423 |
| knn1 | leave_one_run_out | second_axial_half | 0.640 | 0.443 |
| knn1 | leave_one_run_out_consistent_variant | second_axial_half | 0.760 | 0.508 |
| knn3 | leave_one_run_out | second_axial_half | 0.600 | 0.377 |
| knn3 | leave_one_run_out_consistent_variant | second_axial_half | 0.720 | 0.525 |

## Per-pattern recall and mistakes

**centroid, leave_one_run_out**: dispersed bubbles 0.000; plug 1.000; slug 0.667; churn 0.700; wavy annular 0.750.
Errors: m03q05_001 (churn→wavy annular); m04q10_001 (churn→wavy annular); m04q20_001 (churn→wavy annular); m05q05_001 (slug→plug); m06q40_001 (wavy annular→churn); m07q01_001 (dispersed bubbles→plug); m07q80_001 (wavy annular→churn).

**centroid, leave_one_mass_flux_out**: dispersed bubbles 0.000; plug 1.000; slug 0.667; churn 0.700; wavy annular 0.625.
Errors: m03q05_001 (churn→wavy annular); m04q10_001 (churn→wavy annular); m04q20_001 (churn→wavy annular); m05q05_001 (slug→plug); m06q40_001 (wavy annular→churn); m06q80_001 (wavy annular→churn); m07q01_001 (dispersed bubbles→plug); m07q80_001 (wavy annular→churn).

**knn1, leave_one_run_out**: dispersed bubbles 0.000; plug 0.667; slug 0.333; churn 1.000; wavy annular 0.875.
Errors: m04q40_001 (wavy annular→churn); m05q01_001 (plug→slug); m05q05_001 (slug→plug); m06q05_001 (slug→plug); m07q01_001 (dispersed bubbles→plug).

**knn1, leave_one_mass_flux_out**: dispersed bubbles 0.000; plug 1.000; slug 0.667; churn 1.000; wavy annular 0.875.
Errors: m04q40_001 (wavy annular→churn); m05q05_001 (slug→plug); m07q01_001 (dispersed bubbles→plug).

**knn3, leave_one_run_out**: dispersed bubbles 0.000; plug 1.000; slug 0.000; churn 1.000; wavy annular 0.750.
Errors: m03q10_001 (wavy annular→churn); m04q01_001 (slug→plug); m04q40_001 (wavy annular→churn); m05q05_001 (slug→plug); m06q05_001 (slug→plug); m07q01_001 (dispersed bubbles→plug).

**knn3, leave_one_mass_flux_out**: dispersed bubbles 0.000; plug 0.667; slug 0.333; churn 0.800; wavy annular 0.875.
Errors: m03q10_001 (wavy annular→churn); m04q01_001 (slug→plug); m04q10_001 (churn→wavy annular); m04q20_001 (churn→wavy annular); m05q05_001 (slug→plug); m07q01_001 (dispersed bubbles→slug); m07q05_001 (plug→slug).

Fixed-centroid model run-label permutation check (499 permutations): observed balanced accuracy 0.623, exploratory p=0.002; null mean 0.156. This tests label association within this dataset, not gas calibration or performance on new experiments. Unrestricted run-label permutations assume exchangeability; they do not remove acquisition-condition confounding or reuse of the dataset during method development.

## Limits of identification

The dispersed-bubble class now has one selected acquisition. Holding it out leaves no bubble training example, so these data cannot validate a transferable bubble classifier. The absent class is never assigned a fabricated centroid and the unsupported test remains in the denominator. Folds with at least one missing training class: leave-one-run-out 1, leave-one-mass-flux-out 5. Exact present/absent classes are saved for every fold and prediction. Plug and slug each have three examples, so one error moves their recall by one third.

Leave-one-mass-flux-out checks whether signatures transfer to an unseen acquisition condition instead of merely recognizing correlated brightness/advection conditions. Radial/threshold tests perturb the unseen test image without refitting its training images. Half-panorama crops test whether a local segment changes the assigned class. These are adversarial diagnostic checks and are not additional independent validation sets.

The consistent-variant rows additionally fit and validate each ROI/threshold/crop definition on all training and held-out acquisitions using the same definition. They separate changed-feature robustness from the harsher test-only perturbation distribution shift. None of these model/variant scores is a selected optimum.

A classifier can assign a label to every run even when its error rate is high. Report the full confusion matrices and per-run errors, and do not use its predicted class to manufacture the expected void ordering. Correct image-only absolute void fraction and reliable universal five-pattern identification are separate questions.

Reproduce: `python3 src2/v14.1_new_voidfraction/scripts/pattern_identification.py`. Outputs include all features, per-run predictions, confusion matrices, sensitivity scores, the panorama contact sheet and confusion plots.
