# Independent check under the current dataset selection

The current analysis uses 26 acquisitions across 5 single-pattern classes. Mixed labels, smooth annular and m05q05_002 are excluded before computation; only wavy annular represents annular flow. The original input files and their labels are retained.

Independent Otsu arithmetic reproduces all 26 retained optical flip counts and denominators. Maximum absolute rate difference: 0.

The existing >20% saturated-core exclusion leaves 25 statistical acquisitions. The selected per-pattern counts are {'churn': 10, 'wavy annular': 8, 'slug': 3, 'plug': 3, 'dispersed bubbles': 1}. Only m03q40_001 is saturation-unsupported within the retained scope.

The [selection manifest](../results/validation/dataset_selection.json) records the policy and retained runs. [Excluded acquisitions](../results/validation/dataset_exclusions.csv) are bookkeeping only and do not contribute to estimates, references, classifiers or charts.

Previous wider-scope results are archived under intermediate/history/v14_baseline/intermediate/history/pre_scope_exclusions. Their counts and statistics are historical.

The original file conflates two definitions. v11.1 counts Schmitt-trigger 
switches of a one-dimensional dark occupancy. v11.3 counts individual binary 
pixel changes across neighboring panorama columns, averaged over 40 rows. 
The v11.3 values near 0.02 are fractions of adjacent comparisons, 
not automatically switches per 1000 px or temporal frequencies.

Complementing the optical binary mask preserves its flip count exactly. 
Thus swapping dark and bright phase names cannot correct or validate flip-rate. 
Both complement and spatial-reversal invariance pass for every retained image.

The one-view/scale argument also needs correction: an axisymmetric physical 
gas diameter can yield alpha=(d_g/D)^2 with absolute scale canceling. 
Unvalidated optical shadows and unknown cross-sectional shapes are the actual 
obstacles here. The legacy 60-pixel window is not a tube calibration; the current 
primary estimator uses the user-confirmed full interior 18:71 and D=53 pixels.

The data inventory records SHA-256 for every input panorama to make the source 
images used in this check identifiable. The pipeline does not change them.

Reproduce: `python3 src2/v14.1_new_voidfraction/scripts/inventory_and_legacy_audit.py`.
