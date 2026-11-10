# Follow-up tasks for the active analysis

For agent `intermittency`: the panorama silhouette candidate now exists in
`panorama_void_alternatives.py`. Add an independent pattern-identification
comparison using a fixed set of silhouette features (conditional visible void,
center gas duty, candidate-mask pixel flip rate, center flip rate, width
variation, gas interval length/continuity) with the already-tested centroid,
k=1, k=3 models. No hyperparameter tuning. Use leave-one-run-out and
leave-one-mass-flux-out validation. Rebuild the liquid template for each fold
from only the reference runs remaining in that training fold. Keep both bubble
labels and report all errors and recalls. Save `silhouette_pattern_*` outputs.
This tests whether the new phase reconstruction helps classification; its
conditional void values are not physical ground truth. Root is independently
computing candidate-mask flip rates in `phase_mask_flip_analysis.py`.

For agent `physics_sources`: independently challenge all new code and source
claims. Read `panorama_void_alternatives.py`, `pattern_identification.py`,
`spatial_flip_sensitivity.py` and `inventory_and_legacy_audit.py`. Test geometry
arithmetic with multiple disjoint intervals; test inferred-liquid null masks,
review local liquid and gas patches in overlays, and check censoring. A full
53-row mask with D=60 imposes (53/60)^2 = 0.7802778, so near-ceiling annular
values are censored. Treat missing-width intervals only as conditional bounds.
Flag calibration limits, classification leakage, different pixel denominators,
and unsupported claims. Write `adversarial_review.md` and supporting numbers.
Never change original data or authoritative labels.

For agent `locate_data`: please finish the panorama candidate review and ensure
its report treats crop-touching masks and unknown optical inversion honestly.
The user specifically wants panorama-derived void and flip-rate, a plain-language
executive summary, no fabricated values, and all files in this v14 directory.
