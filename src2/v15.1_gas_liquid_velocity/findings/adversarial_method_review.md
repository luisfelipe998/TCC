# Adversarial method review — v15.1

This review defines interpretation and screening requirements before the new spatial/annular results are selected. It does not claim that an untested option succeeds. The retained scope is 26 recordings: 1 dispersed-bubble, 3 plug, 3 slug, 10 churn and 9 wavy-annular recordings. The original five exclusions remain in force.

## Principal challenge

**Camera timing permits image displacement to be expressed in px/s. It does not identify a material phase or calibrate an axial length.** A gas-versus-liquid scatter requires two defensible observables in the same recording. A requested point for every recording cannot justify filling missing coordinates, borrowing another recording's velocity or labeling optical regions as material phases.

The [previous phase-identification audit](../../v15_gas_liquid_velocity/phase_conditioned_motion/findings/independent_review.md) found four reportable passage-conditioned descriptors, but no two-phase material ratio under its stricter fixed/adaptive interior rules. The [previous local-patch study](../../v15_gas_liquid_velocity/frame_motion_evaluation/findings/local_patch_review.md) did recover 434 background-subtracted annular matches among 2,700 attempts, although no annular recording passed its distributed-agreement run gate. This leaves room for a better-defined regional or wave observable; it does not supply ground truth for selecting the next estimator.

## Units and the reference chart

For an axial displacement d in pixels between frame indices separated by k at native rate f, the estimated image velocity is d f/k px/s. A physical conversion additionally requires a justified axial scale s_x in metres/pixel: u=d f s_x/k. Use actual native indices and timing, not export rate, resampled spreadsheet rate or a campaign-level FPS value.

The current [v14 geometry](../../v14_new_voidfraction/scripts/v14_geometry.py) records the owner-confirmed lumen rows 18:71, containing 53 pixel cells. Combined with D=1 mm, s=0.001/53=1.88679×10⁻⁵ m/pixel is a **conditional geometric scale**. Applying a radial span as an axial ruler assumes the corresponding optical mapping and magnification are justified. Square sensor pixels alone do not calibrate cylindrical refraction or phase-dependent optical distortion. The 80-pixel image height is not the lumen diameter. Do not switch between 52 centre-to-centre spacings and 53 occupied pixel cells. Any optional m/s chart must name the assumed scale in the figure/caption, retain px/s source values and avoid calibrated-accuracy wording.

The supplied reference chart has superficial-velocity axes. Superficial velocity is phase volume flow divided by the complete tube area, J_g=Q_g/A and J_l=Q_l/A. With consistently defined phase areas and area/time-weighted intrinsic means, J_g=alpha u_g and J_l=(1-alpha)u_l. Neither a projected optical mask nor a tracked bubble/interface supplies all those quantities. Total mass flux only gives G=rho_g J_g+rho_l J_l under the assumed densities: one equation does not determine both superficial velocities. Therefore the chart layout can be reproduced, but its physical axis names cannot be copied onto optical descriptors.

## Options worth testing, and the observable each can support

| Option | Potential benefit | Principal adversarial challenge |
|---|---|---|
| Gas-passage versus liquid-looking-gap conditioning | Already supports four plug/slug comparisons; reproducible starting reference | Broad motion template can contain both phases and remote boundaries. Changing the probe/guard changes the sampled events. |
| Central versus near-wall grayscale strips | Can investigate continuous annular films and spatially different motion without requiring temporal phase alternation | A spatial strip is not phase-pure. Both strips can see the same refracted wave/interface; empty-core texture may carry no independent gas information. Label axes by image region. |
| Separate upper/lower interface position or thickness signals h(x,t) | Annular wave propagation can remain observable when full-image correlation is poor | The detected contour can be a reflection or crop boundary. Thickness must be resolved over time; wave celerity is not mean liquid material speed. |
| Space–time correlation or spectral-ridge orientation | Uses persistent motion over multiple frames, with less dependence on one adjacent match | Static wall texture, illumination, several propagation slopes, aliasing and a short record can create ambiguous peaks. Temporal reversal and shuffled-time controls are essential diagnostics. |
| Reciprocal local 2-D matching or feature trajectories | Tests actual matched locations and transverse displacement; useful for isolated dispersed bubbles | Trackability selection favors interfaces and persistent objects; matching quality does not establish liquid tracer identity. |
| Dense optical flow | Can visualize spatially heterogeneous motion | Brightness constancy and regularization produce a field even in uninformative regions; smooth vectors are not evidence of phase material velocity. Must compare to residuals/texture and reject unsupported regions. |
| Converting annular wave speed to film mean using a flow model | Could provide an explicitly assumption-dependent inference | A transport model and independently justified inputs are required. It cannot repair missing data merely by assigning a universal factor. |

For example, in a scalar film model h_t+[q(h)]_x=0, wave celerity is c=dq/dh, whereas liquid mean speed is q/h. A pure quadratic flux law gives c=2 q/h; a cubic law gives c=3 q/h. Mixed pressure/shear effects change the relation. These examples demonstrate non-uniqueness and are not proposed conversion factors for the present horizontal boiling tube. Gas-core velocity remains a separate unknown even if film-wave celerity is measured accurately.

## Prespecified controls for the new experiments

1. **Scope and provenance:** preserve all 26 runs in output inventories, including unsupported cases; verify 501 sorted native images, index spacing, both available FPS fields, mask geometry and original-input hashes. No discarded outlier solely because its velocity or class location is inconvenient.
2. **Correspondence:** record texture, NCC peak, competitor margin, search-boundary hits, saturation and exact spatial support. A peak on a finite search boundary is a limit indication, not an accepted velocity. Reciprocal tracking must start at the matched destination, not another fixed source crop.
3. **Temporal checks:** compare at least two nonzero lags after division by lag, using jointly accepted observations and reporting their coverage. Reverse temporal order to test sign. Shuffled-time or sufficiently separated-frame matches test whether the apparent transport persists without temporal correspondence; nonzero null output must not be silently hidden.
4. **Static and synthetic controls:** a constant/low-texture image must be rejected. A stationary textured image should produce zero or be explicitly treated as non-transport, not be forced to a nonzero speed. Known positive/negative/subpixel shifts check sign and scale. A moving texture plus a stationary wall/background tests whether background subtraction meaningfully changes the selected peak.
5. **Spatial checks:** compare at least two axial positions and separate upper/lower near-wall regions where feasible. Report disagreement and band sensitivity as differences in selected optical observations, not an automatically calibrated error bar. Do not average mutually incompatible bands into one confident velocity.
6. **Phase checks:** for any claimed phase-conditioned track, evaluate both its source and matched target against candidate masks; preserve unknown regions and interface neighborhoods. Visually inspect accepted examples as well as failures and extreme values. Consensus between thresholds does not correct shared segmentation mistakes.
7. **Run support:** prespecify minimum accepted observations, distinct time coverage and direction consistency. Show all denominators and retain reasons for non-reportability. Coverage is a feasibility diagnostic, not estimator accuracy; choose the observable/algorithm before inspecting whether classes form appealing clusters.
8. **Averaging:** keep signed per-pair values. Distinguish mean signed velocity, mean speed, median speed and equal episode weighting. Missing estimates are not zero. Overlapping frames from a 71–125 ms record do not constitute independent experimental replicates.

Exact numerical gates should be declared by each estimator before examining its aggregate output. The prior gates can serve as a baseline, but a new observable need not meet an inappropriate old all-region agreement rule. Any changed criterion must be documented as a change in the descriptor or its reporting support, with sensitivity retained.

## Figure and reporting requirements

- Prefer px/s for the primary image-derived chart unless the metadata audit establishes a defensible axial physical calibration. A separate geometric-scale illustration is permissible with its assumptions printed prominently.
- Plot the requested gas-passage versus liquid-gap descriptor only for runs with both quantities. Include an all-run coverage/status panel or table so unsupported recordings remain visible.
- A broader central-versus-near-wall scatter must have those regional names on its axes and be separate from a gas/liquid passage scatter. Do not pool bubble translation, liquid-gap-conditioned global shifts and annular wave speeds as interchangeable phase velocities.
- A logarithmic display may use the magnitude of the signed mean only after a documented nonzero/direction-consistency rule. Keep the sign in the source table and state that all plotted directions satisfy the rule. Do not take absolute values of every individual estimate and hide reversals or near-zero cancellation. A signed linear scatter is preferable if directions differ.
- Use identical accepted rows, units, axis definitions and figure data in English and PT-BR. Use separate symbols/colors for patterns, legible per-run labels, vector PDF/SVG and high-resolution PNG. Draw no literature flow-regime boundaries or inferred phase regions.
- Report per-recording findings and exceptions; unequal pattern counts and only one dispersed recording prevent a universal calibrated ratio for each pattern. The final adversarial audit must verify numeric source rows, displayed rounding, selection and figure units independently of the plotting code.

## Initial disposition

The strongest defensible extension is to test **regional optical motion and annular wave/interface propagation as explicitly different observables**, retain the existing passage-conditioned reference, and assess whether either can supply additional reliable coordinates. No tested method is yet accepted as the physical-accuracy winner. The final memo will be updated after the saved results and visual examples are available.

## Challenge to the proposed mass-conservation inference

The rearrangement u_l=(G-rho_g alpha u_g)/(rho_l(1-alpha)) is arithmetically valid only when its inputs represent compatible material fluxes and phase-area averages. Using an optical structure speed for u_g, an unvalidated projected/geometric alpha and condition-level unsynchronized G does not establish those requirements. If averages are taken separately, unresolved covariance between phase occupancy and velocity also matters; a product of arbitrary time averages need not equal the mean phase flux.

This calculation can be retained as a **stress test of assumptions**, but not as a method for filling the unmeasured liquid coordinate of the primary chart. In particular:

- State the conditional physical scale, densities and G provenance on every branch. Substituting a different pixel scale changes the inferred liquid result through the assumed gas flux; it is not a direct unit conversion of a measured liquid speed.
- Preserve nonphysical domain flags, negative results and missing values; do not clip them to zero or discard them to improve the graph. A failure challenges the combined assumptions and cannot isolate the faulty input.
- Show at least alpha/segmentation and scale sensitivity. At fixed G, densities and u_g, the derivative with respect to alpha is (G/rho_l-rho_g u_g/rho_l)/(1-alpha)^2. The inference becomes ill-conditioned as alpha approaches one unless the numerator happens to cancel.
- A positive finite inferred liquid speed and exact conservation residual demonstrate internal algebraic consistency, not accuracy or independent validation.

For the regional experiment, apply one declared core/near-wall definition consistently to all 26 recordings. An annular interface-wave result may be a separate useful result, but substituting it only for difficult annular coordinates changes the observable. Keep that method change visible in its own panel/table, rather than presenting it as a direct continuation of the same gas/liquid axes.

## Executed independent check: passage results and assumption stress test

The [independent verifier](../scripts/independent_verification.py) imports no analysis or plotting module. Its [saved results](../results/independent_verification.json) reproduce the 26-run scope, 13,000 pair displacements and phase/lag acceptance decisions, episode identifiers, separate signed/episode-equal means, quartiles, four reportable passage ratios and conditional geometric conversions. All four plotted passage cases have exclusively negative accepted displacements for both appearances, so magnitude-of-signed-mean plotting does not conceal sign reversal in these four cases.

The verifier independently reconstructed radial runs using paired rising/falling mask edges, rather than the generating script's sequential row counter. It processed all **39,078 native masks**, corresponding to **39,703,248 frame/column geometric calculations**, and matched all 78 geometry summary rows. The maximum geometry/occupancy discrepancy was 1.11×10⁻¹⁶. The independent check also recovered all fixed-probe occupancies from those masks. This verifies the saved mask arithmetic, not whether the masks depict the correct material phase or circular cross-section.

For all **702 mass-balance stress rows**, phase mass fluxes independently reproduced the inferred velocities, assumed densities, condition-level flow inputs, scale branches, residuals and limitation flags. The maximum absolute numerical difference in the checked stress values was 5.78×10⁻¹⁰, consistent with floating-point reconstruction of a sensitive denominator. All rows remain rejected as material phase-velocity measurements. The 1,150.442 m/s maximum at the declared threshold-0.70/mean-flow/reference-scale baseline reproduces as an assumption-sensitive inference; it is not a measured liquid speed.

No numerical defect was found in the passage or stress outputs. Physical limitations remain unchanged. The subsequent regional and final figure checks are documented below.

## Regional-method challenge found before final selection

The first ensemble implementation checked whether lag-2/lag-4 peak *locations* agreed with lag-1, but did not independently require those longer-lag peaks to have adequate strength or shuffled-time separation. This affected its initial annular interface support: for example, the m04q40 upper-edge lag-4 peak was only 0.0285; other initially supported edge examples had lag-4 peaks of 0.0247–0.0764. These values cannot inherit the much stronger lag-1 support by location agreement alone. The initial m04q40 upper-inner lag-4 peak was also below 0.10.

The parent and annular agent were asked to preserve the initial decisions and evaluate longer-lag peak/null support explicitly, alongside native time reversal and axial-template-location sensitivity, before finalizing support. This is a methodological qualification, not evidence that every original displacement is wrong. The current ten initially paired regional candidates all have negative core and near-wall signed estimates; their final acceptance must use the revised diagnostics.

## Executed independent check: final regional diagnostics

The [regional verifier](../scripts/independent_regional_verification.py) and its [results](../results/independent_regional_verification.json) independently reconstruct all 148 region decisions, 444 lag summaries, 740 temporal-block rows, 8,436 shuffled-time rows and 888 additional diagnostic rows without importing the estimator. Direct centred Pearson products check 1,772 points of the ensemble correlation curves across all region/lag combinations, 444 individual-pair peak values and 78 reconstructed shuffled-trial peaks. The saved float32 curve/track archive explains the maximum peak-location difference of 0.000185 pixels; direct correlation-value differences stay below 2.98×10⁻⁸ for the float32 comparisons.

The revised support requirements retain **93/148 regional descriptors**, compared with 99 under the initial rules, and retain **10/26 paired inner-band descriptors**. Thus the additional adversarial requirements affect some regional/edge claims without changing the ten primary paired regional points. Numerical checks passed. The new method estimates the location of an averaged normalized-correlation peak; it is not the arithmetic mean of individual displacements, a volume-weighted phase speed or a superficial velocity.

The additional time-reversal/axial-position tests remain fixed-window sensitivity diagnostics, not matched-destination reciprocal tracking. The shuffled trials test temporal organization; their maxima plus an empirical separation margin are feasibility gates, not a formal significance level or calibrated error bound.

## Executed independent check: paired variants and final paper figures

The [paper verifier](../scripts/independent_paper_verification.py) and its [saved results](../results/independent_paper_verification.json) independently reconstructed **61 paired variants and 671 fields**, including the final region-support gates, upper/lower agreement, direction requirements, paired values and ratios. It confirmed **10/26 inner-band pairs**, **12/26 outer-band pairs** and **0/9 gradient-location pairs**. The outer-band alternative includes one annular recording, **m04q80_001**; its broader coverage does not establish greater physical accuracy or material-phase identification. No annular inner-band pair passes the final gates.

The exported plotting tables contain exactly the accepted four passage cases, ten inner-band cases and twelve outer-band cases. All 26 recordings remain in the separate diagnostic candidate table; unsupported candidates are not accepted measurements. The verifier checked all 26 run-status rows and **32 numeric cells** in the English/PT-BR executive passage tables. All currently plotted passage means and regional candidates have negative signed velocities before the declared magnitude conversion.

The verifier parsed the **ten exported scatter SVGs** in both languages and independently reconstructed marker coordinates from rendered axis ticks, rather than invoking plotting code. These cover passage px/s, passage conditional m/s, inner-band px/s, outer-band px/s and the all-recording diagnostic. The largest coordinate discrepancy was **0.000544 px/s**, or **3.02×10⁻⁹ m/s** in the conditional physical-scale plots, attributable to SVG coordinate rounding. Corresponding PNG files have 600 dpi metadata and the expected high-resolution dimensions; corresponding PDF files have valid PDF headers. These checks verify content/format consistency, not experimental accuracy.

Visual review of the final English/PT-BR inner-band, outer-band and all-recording PNGs confirmed that the previously overlapping labels had been separated. Near-coincident measured coordinates remain near-coincident, with labels and leader lines identifying the runs. The passage and conditional-m/s plots were also visually reviewed. The primary figures use px/s; the physical-scale illustration prominently states the uncalibrated **1 mm / 53 px** assumption. The all-recording figure explicitly distinguishes rejected hollow points and search-limit saturation, and discloses that its diagnostic wall coordinate averages upper/lower peaks even when they disagree. It must not be used to infer accepted phase-velocity trends from rejected coordinates.

**Final disposition:** the archived results and bilingual scatter coordinates pass the independent numerical, selection, unit and presentation checks. The defensible outputs are image-structure propagation descriptors with documented empirical support, together with a separate conditional geometric conversion. Neither independent gas/liquid material velocities, superficial velocities, a calibrated uncertainty, an accuracy percentage nor a universal per-pattern gas/liquid ratio has been established from these images and metadata.
