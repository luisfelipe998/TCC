# Reorganization verification

**Result: PASS.** Compared against the preserved pre-move [snapshot](../../reorganization_baseline.json).

Checked 23 CSV tables, 31 binary masks and 13 JSON artifacts. CSVs, masks and historical JSON require identical bytes. Other JSON may differ only in the three declared source-hash fields.

Checked 137 local Markdown link destinations: 0 missing. Found 0 unexpected root/scripts entries. Link fragments are outside this check.

## Observed metadata changes

| Artifact and field | Previous SHA-256 | Current SHA-256 |
|---|---|---|
| `results/validation/adversarial_checks.json` — `$.estimator_sha256` | `9e4d7a7c1b5c9ef7bf2d18818cc2b298600f7e41f9006b6774e8b5c89221e835` | `0d1ddf2d1222fe6677d52733626ea5205f23b2a8edf3e3b39be5d750fa82d70c` |
| `results/validation/final_verification.json` — `$.segmentation_source_sha256` | `9e4d7a7c1b5c9ef7bf2d18818cc2b298600f7e41f9006b6774e8b5c89221e835` | `0d1ddf2d1222fe6677d52733626ea5205f23b2a8edf3e3b39be5d750fa82d70c` |
| `results/validation/phase_mask_flip_metadata.json` — `$.segmentation_source_sha256` | `9e4d7a7c1b5c9ef7bf2d18818cc2b298600f7e41f9006b6774e8b5c89221e835` | `0d1ddf2d1222fe6677d52733626ea5205f23b2a8edf3e3b39be5d750fa82d70c` |
| `results/patterns/silhouette_pattern_results.json` — `$.segmentation_source_sha256` | `9e4d7a7c1b5c9ef7bf2d18818cc2b298600f7e41f9006b6774e8b5c89221e835` | `0d1ddf2d1222fe6677d52733626ea5205f23b2a8edf3e3b39be5d750fa82d70c` |

## Failures

No failures were observed.

Exact comparisons and link locations are in [layout_verification.json](../results/validation/layout_verification.json). This verifies preservation and navigation; it does not establish calibrated physical void fractions.
