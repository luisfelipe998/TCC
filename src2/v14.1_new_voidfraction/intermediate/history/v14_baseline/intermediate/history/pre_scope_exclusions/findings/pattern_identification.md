# Image-only pattern identification: feasibility audit

Morphology can help identify regimes without calibrating absolute void fraction. This audit extracts 12 optical/morphology features from each saved panorama. The extraction function receives no labels, mass flux, heat flux or camera timing. Folder labels are used only to train and validate simple classifiers.

The primary existing quality rule leaves 26 labeled runs: two bubbles, three plug, three slug, ten churn and eight wavy annular. Transitions, smooth annular and the >20% saturated wavy-annular run are excluded from these five-class models, but remain in the feature CSV and contact sheet. Each acquisition is one replicate; axial crops do not increase sample size.

Features describe pixel flip density, dark projected area, radial dark-fraction variability, central bright-width mean/variability/absence, radial contrast median/variability/flat-column fraction, and enclosed dark-contour hole area/size/density. They are not validated gas areas. All masks use fixed per-panorama Otsu; contrast-flatness threshold 0.25 and 3×3 closing are fixed across runs.

## Validation results

Feature standardization and class centroids/neighbors are fit on training runs only. Models are a nearest class centroid and k-nearest neighbors at fixed k=1 and k=3; there is no parameter or feature search. All tested models are reported. Balanced accuracy is mean per-class recall; five-class chance is 0.20. Ordinary accuracy has a 10/26=0.385 majority-class baseline.

| Model | Validation | Variant | Accuracy | Balanced accuracy |
|---|---|---|---:|---:|
| centroid | leave_one_run_out | baseline | 0.692 | 0.623 |
| centroid | leave_one_mass_flux_out | baseline | 0.654 | 0.598 |
| knn1 | leave_one_run_out | baseline | 0.731 | 0.508 |
| knn1 | leave_one_mass_flux_out | baseline | 0.808 | 0.642 |
| knn3 | leave_one_run_out | baseline | 0.731 | 0.550 |
| knn3 | leave_one_mass_flux_out | baseline | 0.692 | 0.535 |
| centroid | leave_one_run_out | otsu_minus10 | 0.654 | 0.598 |
| centroid | leave_one_run_out_consistent_variant | otsu_minus10 | 0.654 | 0.598 |
| knn1 | leave_one_run_out | otsu_minus10 | 0.692 | 0.483 |
| knn1 | leave_one_run_out_consistent_variant | otsu_minus10 | 0.769 | 0.575 |
| knn3 | leave_one_run_out | otsu_minus10 | 0.692 | 0.530 |
| knn3 | leave_one_run_out_consistent_variant | otsu_minus10 | 0.808 | 0.642 |
| centroid | leave_one_run_out | otsu_plus10 | 0.538 | 0.398 |
| centroid | leave_one_run_out_consistent_variant | otsu_plus10 | 0.692 | 0.623 |
| knn1 | leave_one_run_out | otsu_plus10 | 0.615 | 0.402 |
| knn1 | leave_one_run_out_consistent_variant | otsu_plus10 | 0.692 | 0.442 |
| knn3 | leave_one_run_out | otsu_plus10 | 0.615 | 0.443 |
| knn3 | leave_one_run_out_consistent_variant | otsu_plus10 | 0.692 | 0.483 |
| centroid | leave_one_run_out | radial_crop15_65 | 0.500 | 0.365 |
| centroid | leave_one_run_out_consistent_variant | radial_crop15_65 | 0.615 | 0.532 |
| knn1 | leave_one_run_out | radial_crop15_65 | 0.615 | 0.430 |
| knn1 | leave_one_run_out_consistent_variant | radial_crop15_65 | 0.731 | 0.550 |
| knn3 | leave_one_run_out | radial_crop15_65 | 0.615 | 0.402 |
| knn3 | leave_one_run_out_consistent_variant | radial_crop15_65 | 0.731 | 0.550 |
| centroid | leave_one_run_out | radial_crop20_60 | 0.192 | 0.180 |
| centroid | leave_one_run_out_consistent_variant | radial_crop20_60 | 0.692 | 0.623 |
| knn1 | leave_one_run_out | radial_crop20_60 | 0.077 | 0.120 |
| knn1 | leave_one_run_out_consistent_variant | radial_crop20_60 | 0.654 | 0.510 |
| knn3 | leave_one_run_out | radial_crop20_60 | 0.231 | 0.167 |
| knn3 | leave_one_run_out_consistent_variant | radial_crop20_60 | 0.692 | 0.535 |
| centroid | leave_one_run_out | first_axial_half | 0.654 | 0.557 |
| centroid | leave_one_run_out_consistent_variant | first_axial_half | 0.654 | 0.557 |
| knn1 | leave_one_run_out | first_axial_half | 0.731 | 0.555 |
| knn1 | leave_one_run_out_consistent_variant | first_axial_half | 0.692 | 0.488 |
| knn3 | leave_one_run_out | first_axial_half | 0.692 | 0.483 |
| knn3 | leave_one_run_out_consistent_variant | first_axial_half | 0.615 | 0.350 |
| centroid | leave_one_run_out | second_axial_half | 0.500 | 0.332 |
| centroid | leave_one_run_out_consistent_variant | second_axial_half | 0.577 | 0.465 |
| knn1 | leave_one_run_out | second_axial_half | 0.577 | 0.418 |
| knn1 | leave_one_run_out_consistent_variant | second_axial_half | 0.731 | 0.508 |
| knn3 | leave_one_run_out | second_axial_half | 0.615 | 0.443 |
| knn3 | leave_one_run_out_consistent_variant | second_axial_half | 0.692 | 0.525 |

## Per-pattern recall and mistakes

**centroid, leave_one_run_out**: dispersed bubbles 0.000; plug 1.000; slug 0.667; churn 0.700; wavy annular 0.750.
Errors: m03q05_001 (churn→wavy annular); m04q10_001 (churn→wavy annular); m04q20_001 (churn→wavy annular); m05q05_001 (slug→plug); m05q05_002 (dispersed bubbles→plug); m06q40_001 (wavy annular→churn); m07q01_001 (dispersed bubbles→plug); m07q80_001 (wavy annular→churn).

**centroid, leave_one_mass_flux_out**: dispersed bubbles 0.000; plug 1.000; slug 0.667; churn 0.700; wavy annular 0.625.
Errors: m03q05_001 (churn→wavy annular); m04q10_001 (churn→wavy annular); m04q20_001 (churn→wavy annular); m05q05_001 (slug→plug); m05q05_002 (dispersed bubbles→plug); m06q40_001 (wavy annular→churn); m06q80_001 (wavy annular→churn); m07q01_001 (dispersed bubbles→plug); m07q80_001 (wavy annular→churn).

**knn1, leave_one_run_out**: dispersed bubbles 0.000; plug 0.667; slug 0.000; churn 1.000; wavy annular 0.875.
Errors: m04q01_001 (slug→plug); m04q40_001 (wavy annular→churn); m05q01_001 (plug→slug); m05q05_001 (slug→plug); m05q05_002 (dispersed bubbles→slug); m06q05_001 (slug→plug); m07q01_001 (dispersed bubbles→plug).

**knn1, leave_one_mass_flux_out**: dispersed bubbles 0.000; plug 1.000; slug 0.333; churn 1.000; wavy annular 0.875.
Errors: m04q01_001 (slug→plug); m04q40_001 (wavy annular→churn); m05q05_001 (slug→plug); m05q05_002 (dispersed bubbles→plug); m07q01_001 (dispersed bubbles→plug).

**knn3, leave_one_run_out**: dispersed bubbles 0.000; plug 1.000; slug 0.000; churn 1.000; wavy annular 0.750.
Errors: m03q10_001 (wavy annular→churn); m04q01_001 (slug→plug); m04q40_001 (wavy annular→churn); m05q05_001 (slug→plug); m05q05_002 (dispersed bubbles→slug); m06q05_001 (slug→plug); m07q01_001 (dispersed bubbles→plug).

**knn3, leave_one_mass_flux_out**: dispersed bubbles 0.000; plug 0.667; slug 0.333; churn 0.800; wavy annular 0.875.
Errors: m03q10_001 (wavy annular→churn); m04q01_001 (slug→plug); m04q10_001 (churn→wavy annular); m04q20_001 (churn→wavy annular); m05q05_001 (slug→plug); m05q05_002 (dispersed bubbles→plug); m07q01_001 (dispersed bubbles→slug); m07q05_001 (plug→slug).

Fixed-centroid model run-label permutation check (499 permutations): observed balanced accuracy 0.623, exploratory p=0.002; null mean 0.187. This tests label association within this dataset, not gas calibration or performance on new experiments.

## Limits of identification

The two bubbles acquisitions have incompatible-looking morphology: m07q01_001 is mostly a nearly flat liquid-like profile with small dispersed shadows; m05q05_002 contains large discrete dark-outlined bodies. Their supplied labels are kept unchanged. With only two labeled acquisitions, leaving one out leaves one bubble training example and cannot establish a stable bubble class. Plug and slug likewise have only three independent examples each, so each error moves recall by one third.

Leave-one-mass-flux-out checks whether signatures transfer to an unseen acquisition condition instead of merely recognizing correlated brightness/advection conditions. Radial/threshold tests perturb the unseen test image without refitting its training images. Half-panorama crops test whether a local segment changes the assigned class. These are adversarial diagnostic checks and are not additional independent validation sets.

The consistent-variant rows additionally fit and validate each ROI/threshold/crop definition on all training and held-out acquisitions using the same definition. They separate changed-feature robustness from the harsher test-only perturbation distribution shift. None of these model/variant scores is a selected optimum.

A classifier can assign a label to every run even when its error rate is high. Report the full confusion matrices and per-run errors, and do not use its predicted class to manufacture the expected void ordering. Correct image-only absolute void fraction and reliable universal five-pattern identification are separate questions.

Reproduce: `python3 src2/v14_new_voidfraction/scripts/pattern_identification.py`. Outputs include all features, per-run predictions, confusion matrices, sensitivity scores, the panorama contact sheet and confusion plots.
