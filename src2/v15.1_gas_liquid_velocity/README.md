# v15.1 — Gas/liquid velocity and optical-motion investigation

Start with the [EXECUTIVE SUMMARY](findings/EXECUTIVE_SUMMARY.md) or the [RESUMO EXECUTIVO EM PT-BR](findings/EXECUTIVE_SUMMARY_PT_BR.md), updated on **24 October 2026**. Both now show the latest all-recording chart: **26 optical points, 25 with reference support, 7 stable across all four sensitivity checks, and 0 independently identified material gas/liquid pairs**.

The [complete adversarial investigation](all_scenarios_adversarial/README.md) contains the new methods, controls, numerical tables and independent checks. The earlier [English](findings/EXECUTIVE_SUMMARY_2026-10-23.md) and [PT-BR](findings/EXECUTIVE_SUMMARY_PT_BR_2026-10-23.md) summaries remain archived for provenance.

The main figures report **px/s**. The optional geometric conversion uses **1 mm / 53 pixel cells** and is explicitly conditional. Camera timing is documented; an axial object-space calibration and independently identifiable liquid tracers are unavailable. The example flow-map layout is used without copying its superficial-velocity axis definitions or drawing pattern boundaries.

## Results and figures

| Current artifact | Purpose |
|---|---|
| [All 26 recordings — PT-BR](all_scenarios_adversarial/figures/all_26_selected_wall_PT_BR.png) | Core versus selected-wall optical propagation; support, sensitivity and saturation remain visible |
| [Publication PDF](all_scenarios_adversarial/figures/all_26_selected_wall_PT_BR.pdf) | Same complete chart in vector format |
| [Current method coverage](all_scenarios_adversarial/figures/method_coverage_PT_BR.png) | Includes the new 20-pair two-wall and 25-pair selected-wall results |
| [Exact current plot data](all_scenarios_adversarial/results/final_all_26_plot_data.csv) | Signed values and support for all 26 recordings |
| [Current figure captions](findings/FIGURE_CAPTIONS_PT_BR.md) | Definitions and interpretation limits for the latest figures |

### Earlier results — 23 October 2026

| Artifact | Purpose |
|---|---|
| [Passage scatter in PT-BR](figures/gas_liquid_passages_px_s_PT_BR.png) | Separate motion means during gas-like passages and liquid-looking gaps |
| [Regional scatter in PT-BR](figures/regional_motion_px_s_PT_BR.png) | Supported core/near-wall propagation descriptors; these are spatial regions, not phases |
| [Outer-band alternative in PT-BR](figures/regional_outer_motion_px_s_PT_BR.png) | Same fixed outer bands across all recordings, including one supported annular optical pair |
| [All-run diagnostic in PT-BR](figures/all_regional_candidates_PT_BR.png) | Every computed regional candidate, with failed criteria visibly distinguished |
| [Annular support in PT-BR](figures/annular_region_support_PT_BR.png) | Supported optical propagation by region, with missing estimates preserved |
| [Conditional m/s illustration](figures/gas_liquid_passages_conditional_m_s_PT_BR.png) | Same passage result multiplied by an assumed geometric scale |
| [Earlier figure captions](findings/FIGURE_CAPTIONS_PT_BR_2026-10-23.md) | Captions for the original passage and regional figures |
| [Calibration review](findings/calibration_review.md) | Fresh camera/JPEG/workbook metadata and original de Oliveira/Dias evidence |
| [Adversarial method review](findings/adversarial_method_review.md) | Method challenges, rejected inferences and independent verification |
| [Annular method review](findings/annular_methods.md) | New regional/gradient alternatives, stronger temporal tests and optical-scale sensitivity |
| [Passage table](results/passage_per_run.csv) | All 26 retained recordings, signed means, episode weighting and withholding reasons |
| [Regional candidate table](results/annular_spatial_pairs.csv) | Original fixed-rule regional estimates and support flags |
| [Mass-balance stress table](results/mass_balance_stress.csv) | 702 rejected assumption-dependent inferences, never accepted phase speeds |

All paper figures have **English and PT-BR versions, PDF/SVG vector exports and 600-dpi PNGs**. CSV plot data retain signs and exact selection. A point absent from a supported-only plot is not a zero velocity. The all-run diagnostic deliberately shows failed algorithm candidates as unfilled markers; they are not endorsed measurements.

## Scope and preservation

The original exclusions remain: m05q05_002, smooth annular and double classifications. The [inventory](results/scenario_inventory.csv) retains all 31 camera folder identities and selects 26. Each retained folder supplies 501 consecutive native images, giving 500 frame intervals. The frame rates are 4,000, 5,000 or 7,000 fps.

New analyses use images and metadata from src2/data_processed. The previous native-frame phase masks and passage tracking are reused with recorded hashes and independently reaggregated; those estimators are not silently refitted. The new regional/annular study reads the native JPEGs again. v12, v12.1, v14, v15 and original data are preserved. This directory is local because src2 is excluded by the existing Git configuration.

## Reproduction

The scripts resolve input and output paths from their own locations. They use the existing Python environment; no source download or package installation is required. The previous v15 native-frame/mask results remain explicit inputs.

Run the complete additional investigation from the repository root:

    python3 -B src2/v15.1_gas_liquid_velocity/all_scenarios_adversarial/scripts/run_all.py

Rebuild the main report entry points from the latest saved reports:

    python3 -B src2/v15.1_gas_liquid_velocity/scripts/publish_latest_summary.py

The earlier analysis can still be reproduced with:

    python3 -B src2/v15.1_gas_liquid_velocity/scripts/run_all.py

To rebuild the earlier figures and archived summaries from saved numerical results:

    python3 -B src2/v15.1_gas_liquid_velocity/scripts/paper_figures.py
    python3 -B src2/v15.1_gas_liquid_velocity/scripts/build_summary.py

The summary builder keeps the earlier findings in the dated archives and republishes the latest investigation at the main summary paths. It does not replace the current report with the earlier four-point result.

Numerical verification establishes arithmetic, provenance and correspondence diagnostics. It does not establish calibrated phase-material accuracy. Temporal blocks and frame pairs are not independent experimental repetitions.
