# Independent adversarial review after regeneration

**Result: PASS.** This check ran after the analysis and executive summary were written.

Independent reading of original folder labels and pixels finds **26 selected acquisitions and 25 unsaturated statistical acquisitions**, across the five approved single-pattern classes. Excluded data are allowed only in explicit exclusion bookkeeping and historical outputs.

| Check | Result |
|---|---|
| scope and all active csvs | PASS |
| independent mask arithmetic and input hashes | PASS |
| retained measurements against archive | PASS |
| calibration and computational metadata scope | PASS |
| groups charts scaling and provenance | PASS |
| heldout classifier tests and missing classes | PASS |
| written executive summary | PASS |

## Findings

The audit tests selection consistency, saved mask arithmetic, retained default measurements, chart coordinates and scaling, held-out predictions and the numbers written in the executive summary. Excluded runs must stay out of current tables, references and classifier tests.

**Limit:** a single retained bubble acquisition cannot validate a transferable bubble classifier. Its unsupported held-out test remains in the scores; absent training classes receive no manufactured centroid.

**Limit:** passing these checks establishes computational and scope consistency, not correct three-dimensional gas shape or calibrated physical gas volume. The crop and optical-mapping limitations remain.

Detailed observations and failures are recorded in [post_regeneration_adversarial_checks.json](../results/validation/post_regeneration_adversarial_checks.json).

