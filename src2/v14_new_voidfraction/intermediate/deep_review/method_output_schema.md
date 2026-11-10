# Method audit implementation and output schema

Additional verification completed: all 26 corrected production masks equal an independently reconstructed mask using Boolean OR/AND neighborhood closing (no SciPy morphology). This did not change current saved masks. Component labeling uses SciPy's default four-neighbor connectivity; the closing element is full 3×3.

Production estimator changed only the closing boundary policy numerically: pad the binary crop by one zero pixel, apply 3x3 closing, then crop back before component labels and radial-contact rejection. Axial contact is retained. Documentation now correctly describes generic >=4-pixel shadow candidates, internal (not independent) controls, row40 probe, and edge-contact extension scenarios. Dynamic liquid-holdout numbers replace hardcoded prose.

`deep_method_checks.py` compiles and its `analytic_controls()` passed. It verifies 10 masks by direct nested-loop dilation/erosion, a separate Boolean-neighborhood implementation, and production closing; also radial rejection, axial preservation, ring filling, three uniform illuminations, bounded noise, inferred liquid profile, a bright gas interior, and interval-square geometry.

Run the full script after current panorama/mask outputs regenerate. It strictly checks all saved masks against independent corrected reconstruction, each old mask against the archived old source, per-run alpha against current tables and flip counts against current tables. The archive used is `intermediate/history/before_deep_method_review`.

Outputs registered by root:
- `deep_method_changes.csv`: run, pattern, eligible, alpha_before, alpha_after, alpha_change, flip_count_before, flip_count_after, flip_comparisons, flips_per1000_before, flips_per1000_after, changed_pixels, added_pixels, removed_pixels, fixed_row40_flips_per1000, apparent_midpoint_row45_flips_per1000, input_panorama_sha256, before_mask_sha256, after_mask_sha256.
- `deep_method_checks.json`: passed, completed_utc, controls, selected_acquisitions, statistical_acquisitions, changed_masks, max_abs_alpha_change, source_before_sha256, source_after_sha256, segmentation_source_sha256, group_comparisons, measurements (same row schema as CSV), failures, interpretation.
- `deep_method_review.md`: simple PT-BR audit report generated from the checks with before/after medians and limitations.

Initial diagnostic numbers (prior to full output regeneration): 24/26 masks change; max absolute alpha change 0.0054171545864458714 on m07q20_001. Groups after correction: dispersed 0.0013796909492273732 and 2.4641390393428964 flips/1000; plug 0.46281243766207864 and 9.730399452135678; slug 0.4730783499564802 and 10.41884285484054; churn 0.6418732294956717 and 11.19010973471768; wavy annular 0.7673458119311427 and 1.132855347038008. These come from the actual diagnostic calculation, not final regenerated tables.

The existing first diagnostic evidence is preserved in `intermediate/method_review_20261023/`; its diagnostic script was run before production was modified. Do not blindly rerun that old scratch script as a before/after source; use the new production audit, which loads the archived before source explicitly.
