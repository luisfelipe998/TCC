# Independent adversarial review of the completed comparison

This review covers the earlier thermodynamic investigation. References below to its executive table refer to the preserved [conditional model summary](conditional_model_summary.md). The [current executive summary](EXECUTIVE_SUMMARY.md) reports the later phase-passage image-motion analysis in px/s, with its own linked verification.

This is the independent review of the [historical conditional model summary](conditional_model_summary.md) and its model tables. The [current executive summary](EXECUTIVE_SUMMARY.md) reports the later pixel-based phase-passage analysis, whose separate [independent review](../phase_conditioned_motion/findings/independent_review.md) addresses its own evidence and limitations.

**Result: the tested thermodynamic calculations, aggregations and displayed velocity tables agree with independent recomputation. No unresolved numerical defect was found. This is an implementation/reproducibility conclusion, not experimental validation of phase velocities.**

Reproduction: `python3 src2/v15_gas_liquid_velocity/scripts/independent_model_review.py`. Full residuals, checked source/result hashes and reviewed narrative hashes are in [independent_model_review.json](../results/validation/independent_model_review.json). The review does not rewrite the model script or its result tables.

## Scope and independent method

The reviewer independently enumerated source folders and marker files, recovering the same 26 selected acquisitions and five exclusions. The thermodynamic table contains exactly the Cartesian product **26 acquisitions × 13 input scenarios × 4 closures = 1,352 unique rows**, with no duplicate key or scope leakage. Every retained pattern agrees with its acquisition marker.

The independent reconstruction starts from total heater power and mass flow, computes outlet enthalpy, and converts that enthalpy to equilibrium quality using saturated liquid/vapor enthalpies. It does not call `model_comparison.py` to generate the comparison values. It independently obtains properties from CoolProp; that shared property library and the declared physical assumptions remain common inputs.

- **RA:** calculate drift velocity in m/s, then `ug = C0 (Jg+Jl) + Vgj`, `alpha = Jg/ug` and `ul = Jl/(1-alpha)`. Compare this result with both the new saved rows and the separate legacy reconstruction.
- **Homogeneous and secondary Xu–Fang:** execute only the two inspected function definitions extracted with Python AST from the pinned local `fluids_two_phase_voidage.py` snapshot. The full external module is not imported. Pass **g=9.81 m/s² explicitly**, matching the experiment rather than accepting the library default. This independently verifies the implementation with exponent −0.2.
- **Literal 2021 Xu–Fang:** independently derive `S = 1 + 2 Fr_lo^(-2) alpha_h^3.5`, then `ug = Jg + S Jl`, `ul = Jg/S + Jl`. The saved high-resolution source-equation crop was also viewed directly and visibly prints exponent −2. This remains distinct from the secondary −0.2 implementation; source authority is not inferred from numerical agreement.

## Executed checks

| Check | Result |
|---|---:|
| Per-row numeric cells compared | 19,632 |
| Independently recomputed group-summary cells | 1,200 |
| Displayed numeric cells checked in the historical model-summary/per-run reports | 441 |
| Maximum absolute phase-velocity difference | 2.487e−13 m/s |
| Maximum RA difference from the independent legacy reconstruction | 1.776e−14 across compared quantities |
| Mass flow recovered from phase speeds and occupied areas | Maximum residual 1.626e−18 kg/s |
| Full enthalpy balance, including negative-quality cases | Maximum residual 4.263e−14 W |
| Image-alpha diagnostic rows | 78; 75 numerically evaluated, 3 saturation exclusions |
| Inverted RA alpha maximum residual | 8.039e−12 |

The mass check reconstructs `mdot = A [rho_g alpha ug + rho_l (1-alpha) ul]`; it is not limited to adding superficial flux columns. The energy check verifies `mdot (h_out-h_in) = Q_net` independently of the reduced quality-increment expression.

All four closures give the same superficial velocities for the same state. The homogeneous model gives equal gas/liquid phase speeds. Invalid equilibrium qualities retain undefined two-phase outputs; explicit endpoint tests confirm that the absent phase's speed remains undefined, while the present phase has its correct single-phase velocity. The 1/2/5/10 K subcooling tests yield **4/4/7/12** out-of-domain acquisitions, respectively. Those counts are acquisitions, not the four-times-larger number of closure rows.

Recorded mean, median and trimmed-mean G branches use the corresponding recorded-flow input columns. Their Jg agrees with the nominal-G branch to at most 1.87e−14 m/s under zero inlet quality, and their Jl change agrees with `(G_recorded-G_nominal)/rho_l`. This verifies the cancellation without misrepresenting the recorded-flow branch as a synchronized camera measurement.

Both nominal and recorded-flow group tables were reconstructed using independent arithmetic means/medians over **per-acquisition derived values**, with correct n. Mean slip is the mean of the individual slip ratios. The historical model summary correctly labels its table as acquisition medians and differentiates superficial velocities from conditional RA phase predictions. All 52 displayed per-acquisition rows reproduce at their stated precision. Closure-spread ratios and all relative sensitivity percentages were checked separately.

## Narrative and provenance assessment

The reviewed historical model summary, method comparison and numerical-results report correctly preserve the material limitations: unknown inlet enthalpy/net heat/state, condition-only flow/video mapping, unvalidated volumetric alpha, optical length scale, source-equation conflict and absence of independent phase-velocity truth. The recorded-flow evidence is correctly distinguished from placeholder zeros in the image-derived workbooks. RA is presented as a provisional reporting reference, not an accuracy winner.

The reported maximum nominal liquid-velocity model ratio **4.40×** reproduces: the unrounded maximum is 4.403119 at m04q80_001. It describes tested model spread. The secondary Xu–Fang Froude flags are recomputed with the relevant G, yielding **11/26 nominal-G** and **15/26 recorded-mean-G** acquisitions beyond the cited secondary-source limit. That domain statement does not validate the remaining acquisitions or provide error bars.

Source/result hashes agree with the saved metadata, including mass-flow input, the frozen v14 image-alpha snapshot, and the model script. A broad metadata field that originally called all G values nominal was clarified during review to apply specifically to the baseline; the final checked metadata distinguishes recorded-flow branches. No unresolved numerical or metadata correction is requested.

This review independently checks thermodynamic and hybrid numerical results. The optical algorithm's synthetic and real-image correspondence evidence is documented in [optical investigation](optical_investigation.md); the numerical success of this thermodynamic review does not validate optical material velocities. Neither conservation residuals near machine precision nor visually plausible clusters determine true physical accuracy.
