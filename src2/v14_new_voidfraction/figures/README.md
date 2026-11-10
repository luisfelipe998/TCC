# Main graphics

[Back to v14](../README.md)

| Graphic | What to inspect |
|---|---|
| [Void fraction × spatial flip-rate](void_fraction/void_fraction_vs_flip_rate.png) · [SVG](void_fraction/void_fraction_vs_flip_rate.svg) | Eligible single patterns and observed sample spread; saturated acquisition omitted |
| [Fixed-row y=40 comparison](void_fraction/void_fraction_vs_center_flips.png) · [SVG](void_fraction/void_fraction_vs_center_flips.svg) | Gas/liquid switches at row y=40, distinct from the 53-row average |
| [Crop sensitivity](void_fraction/void_crop_sensitivity.png) · [SVG](void_fraction/void_crop_sensitivity.svg) | Effect of removing known-interior rows, with D=53 fixed |
| [Gas-silhouette overlays](void_fraction/panorama_phase_masks.png) | Current inferred envelopes and small shadow components |
| [Liquid-reference candidates](void_fraction/panorama_liquid_reference_candidates.png) | Exact inferred liquid intervals used by the reconstruction |
| [Panorama contact sheet](patterns/pattern_panorama_contact_sheet.png) | All acquisitions and supplied pattern labels |
| [Morphology confusion matrices](patterns/pattern_confusion_matrices.png) | Held-out pattern classification errors |
| [Silhouette confusion matrices](patterns/silhouette_pattern_confusion_matrices.png) | Comparison using reconstructed gas silhouettes |
| [Adversarial phase-anchor patches](validation/adversarial_phase_anchor_patches.png) | Selected visually interpreted liquid/gas interiors |
| [Corrected slug liquid-gap overlay](validation/adversarial_slug_gap_mask.png) | The demonstrated false-bridge case |
| [Spatial versus native optical flips](flip_rate/flip_rate_spatial_vs_native.png) | Supplemental comparison of different sampling bases |
| [Native optical kymographs](flip_rate/flip_rate_native_kymographs.png) | Supplemental time-record checks |

Graphics are grouped into void_fraction, flip_rate, patterns and validation. Earlier alternatives are in [intermediate/exploratory/figures/](../intermediate/exploratory/figures); preliminary review graphics remain in [intermediate/history/](../intermediate/history).

The pattern charts use [these exact coordinates](../results/void_fraction/void_fraction_flip_plot_data.csv); [chart checks](../results/validation/summary_chart_checks.json) record input hashes and arithmetic verification. Colored hulls show observed sample spread, not classification boundaries. Neither axis is independently calibrated physical truth.
