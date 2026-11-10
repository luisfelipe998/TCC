# v14: void fraction and flip-rate

Start with the [executive summary](findings/EXECUTIVE_SUMMARY.md) or [resumo em PT-BR](findings/EXECUTIVE_SUMMARY_PT_BR.md) for findings, next steps and limitations. The current analysis uses 26 retained panoramas from src2/data_processed: five single patterns, wavy annular only, with the user-excluded replicate omitted. The [selection manifest](results/validation/dataset_selection.json) records all exclusions.

**Interpretation:** void values are conditional visible circular reconstructions; gas-silhouette panorama flip-rates are spatial changes between neighboring pixels. Absolute physical void fraction remains uncalibrated, and radial clipping affects several annular runs.

## Find what you need

| Folder | Contents |
|---|---|
| [findings/](findings/README.md) | Executive summary, methods, limitations, adversarial review and literature audit |
| [results/](results/README.md) | Per-run and grouped numerical tables, predictions and current verification |
| [figures/](figures/README.md) | Main gas-mask overlays, reference locations, panorama contact sheet and comparison graphics |
| [scripts/](scripts/README.md) | Reproducible analysis code, ordered stages and shared artifact paths |
| [intermediate/](intermediate/README.md) | Features, sensitivity sweeps, binary masks, exploratory work, history, logs and caches |
| [literature/](literature/README.md) | Extracted reference text and inspected pages; original PDFs remain in assets |

## Main results

| Result | Per acquisition | By supplied pattern |
|---|---|---|
| Conditional visible void estimate | [Retained-run table](results/void_fraction/panorama_void_per_run.csv) | [Group summary](results/void_fraction/panorama_void_by_pattern.csv) |
| Spatial gas-silhouette flip-rate | [Counts, rates and duty](results/flip_rate/phase_mask_flip_per_run.csv) | [Group summary](results/flip_rate/phase_mask_flip_by_pattern.csv) |

Legacy brightness-based flip comparisons are in [results/flip_rate/](results/flip_rate). They use a different mask and radial crop; read the [definitions](findings/phase_mask_flip_review.md) before comparing columns.

## Main graphics

- [Void fraction × spatial flip-rate by pattern](figures/void_fraction/void_fraction_vs_flip_rate.png): acquisitions passing the saturation criterion and observed sample spread.
- [Fixed-row y=40 flip-rate version](figures/void_fraction/void_fraction_vs_center_flips.png): a single-chord comparison that misses off-center bubbles.
- [Void crop sensitivity](figures/void_fraction/void_crop_sensitivity.png): how much the result changes when the radial crop changes.
- [Gas-mask overlays](figures/void_fraction/panorama_phase_masks.png): inspect which visible regions the current method counts as gas.
- [Inferred liquid-reference locations](figures/void_fraction/panorama_liquid_reference_candidates.png): inspect the actual reference intervals.
- [All-panorama contact sheet](figures/patterns/pattern_panorama_contact_sheet.png): compare acquisitions and supplied labels.
- [Morphology classifier confusion matrices](figures/patterns/pattern_confusion_matrices.png) and [silhouette classifier comparison](figures/patterns/silhouette_pattern_confusion_matrices.png).
- [Adversarial liquid-gap example](figures/validation/adversarial_slug_gap_mask.png): inspect the corrected false-bridge case.

## Reproduce

From the repository root:

```sh
python3 -B src2/v14_new_voidfraction/scripts/run_all.py
```

Add --include-exploratory to regenerate the earlier raw-frame alternatives as well. The runner resolves inputs and outputs from its own location, so it also works when invoked by absolute path from another directory. It runs sequentially, stops on failure, and saves stage logs plus a run manifest under [intermediate/logs/](intermediate/logs).

See the [script index](scripts/README.md) for individual stages and dependencies. The [final verification](results/validation/final_verification.json) independently checks saved-mask arithmetic, input hashes, labels and classifier scores. The [current scope verification](findings/layout_verification.md) checks retained measurements, current exclusions and navigation against the archived earlier outputs. Historical diagnostics remain in intermediate/history; the [current adversarial checks](results/validation/adversarial_checks.json) identify the current segmentation source.
