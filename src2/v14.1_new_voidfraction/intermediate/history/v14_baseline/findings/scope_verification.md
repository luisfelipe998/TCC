# Current scope and artifact verification

**PASS:** current scope contains 26 acquisitions across five single-pattern classes. Only wavy annular represents annular flow.

Compared 5 per-run tables against the [earlier outputs](../intermediate/history/pre_scope_exclusions/README.md). Input observables and labels are preserved. Against the [pre-review snapshot](../intermediate/history/before_deep_method_review/README.md), 24 masks change under the independently verified border correction and 2 remain identical. Derived values, groups and classifiers are recomputed.

The user then confirmed the complete lumen as slice 18:71, D=53. Against the [pre-normalization snapshot](../intermediate/history/before_confirmed_lumen_D53/README.md), all 26 masks and 26 sets of flip metrics remain identical. The alpha rescaling by (60/53)^2 is independently verified; the historical morphology comparison remains at D60 on both sides.

Checked 26 active CSV artifacts for excluded runs, labels and reference inputs. Current source hashes match the segmentation script.

Exact checks: [scope_verification.json](../results/validation/scope_verification.json). Independent review after the summary: [post-regeneration adversarial review](post_regeneration_adversarial_review.md).

These checks verify exclusion consistency, documented corrections and arithmetic; absolute gas-volume calibration remains unresolved.

The D60 snapshot preserves original bytes and relative links. Its file hashes are verified, and its links are evaluated in the original v14 layout. Current report links are checked at their current locations.

## Failures

None.
