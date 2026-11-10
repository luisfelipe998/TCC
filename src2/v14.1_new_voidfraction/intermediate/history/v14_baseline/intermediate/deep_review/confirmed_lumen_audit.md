# Confirmed lumen audit work

Implementation complete. The historical morphology ledger remains normalized at D=60 so its before/after comparison isolates the border correction. Current D=53 is checked separately, with unchanged masks and flip metrics against before_confirmed_lumen_D53.

Root: please register scripts/lumen_normalization_checks.py and these new output artifact names, and add its stage after phase-mask results and before post-regeneration checking:

- lumen_normalization_checks.json → results/validation/lumen_normalization_checks.json
- lumen_normalization_changes.csv → results/validation/lumen_normalization_changes.csv
- lumen_normalization_review.md → findings/lumen_normalization_review.md

No pipeline execution by this agent.

## Completed changes

- deep_method_checks.py: current analytic geometry uses D53 and verifies full mask=1; historical morphology comparison still computes alpha_before, alpha_after and alpha_change at D60. Added morphology_comparison_diameter_px=60, current_alpha_D53 and current_normalization_diameter_px=53 per row, and clear report labels. Historical maximum morphology change remains comparable to the old audit; it is not a normalization delta. Removed the outside-current-lumen hidden-bubble control from current geometry claims. The diagnostic second probe now uses confirmed midpoint row44, while primary row40 stays fixed.
- lumen_normalization_checks.py: new independent comparison with before_confirmed_lumen_D53. Requires byte-identical masks, unchanged phase flip metrics/counts/duty, unchanged original image hashes, alpha rescale by (60/53)^2, projected rescale by60/53, mean/median from independent D53 mask arithmetic, all current fractions within0..1, exact user slice/diameter/source, and removal of four obsolete primary fields. Reads current outputs and writes the three artifacts above.
- final_verification.py: asserts exact user geometry, checks area0..1, and invokes the normalization audit. Records boundary provenance and geometry-source hash.
- post_regeneration_adversarial.py: current area computation usesD53, historical morphology rows are checked atD60, separately checks current_alpha_D53. Added independent normalization audit category (also reconstructs the saved normalization ledger), full_lumen_maximum=1 chart check, and user-provenance requirements for both summaries. Removed obsolete demand that current summaries mention0.7803. Historical discussion remains allowed.
- verify_layout.py: checks historical morphology diameter60, independent D53 normalization, saved audit hashes, and describes the two changes separately in scope reports.

## Checks completed before regeneration

All five edited/new scripts compile. Ten morphology controls and the illumination/noise/profile/geometry controls pass under D53. Direct reconstruction of all 26 selected masks equals the before_confirmed_lumen_D53 masks exactly. Independent D53 mean area agrees with saved old D60 alpha times (60/53)^2; the largest absolute numerical difference was 3.3306690738754696e-16. These checks did not modify current output tables.

The full scripts must run after regenerated outputs exist. No PASS claim is made here about regenerated files that have not yet been produced. The new normalization stage must run before the post-regeneration adversarial audit; other audits also invoke its read-only verification function.
