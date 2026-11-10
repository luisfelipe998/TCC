# What can be calculated, and which method to use

This report documents the earlier conditional model investigation. Its m/s predictions and unresolved physical inputs remain relevant to that model. The [current executive summary](EXECUTIVE_SUMMARY.md) presents the later image-motion analysis in px/s, including v14 phase candidates and the limits of separate gas/liquid identification. See the [historical model summary](conditional_model_summary.md) for the results associated with the equations below.

**The evidence supports conditional gas/liquid velocity calculations and optical structure tracking. It does not establish true phase velocities or an experimentally most accurate model.** We compare input provenance, conservation, failure modes and sensitivity. No model is fitted to the pattern labels or to inconsistent published annotations.

## 1. Three different quantities

| Quantity | Definition | Required information | Current status |
|---|---|---|---|
| Superficial gas/liquid velocity | Gas/liquid volume flow divided by whole tube area | Mass flux G, vapor mass quality x, local densities | Calculable under stated inlet/energy/property assumptions |
| Cross-sectional mean phase velocity | ug = Jg/alpha; ul = Jl/(1-alpha) | The above plus cross-sectional volumetric void fraction alpha | Correlation-based estimates; independently calibrated alpha unavailable |
| Visible-structure celerity | Signed displacement / elapsed time | Corresponding image features and native timing; length calibration for m/s | Tested in px/frame and px/s; cannot identify both phase velocities |

G is total mass flow divided by tube area, in kg/(m² s). The documented Voltage_1 is **mass flow in g/s**, despite its voltage-like column name. Convert to kg/s first, then divide by A = pi D²/4. Vapor quality x is a mass-flow fraction, distinct from projected gas area, brightness duty and volumetric void fraction.

At a consistent local two-phase state:

    Jg = G x / rho_g
    Jl = G (1-x) / rho_l
    rho_g Jg + rho_l Jl = G
    ug = Jg / alpha
    ul = Jl / (1-alpha)
    S = ug/ul = x rho_l (1-alpha) / [(1-x) rho_g alpha]

These identities require neither homogeneous flow nor equal phase speeds. No slip enters only if alpha is assigned the homogeneous value. Declining Jl does not by itself establish liquid deceleration: its occupied area may also decline. The [legacy audit](legacy_audit.md) gives a numerical counterexample from the old results.

## 2. Energy balance and the unresolved inlet

For steady mass flow, negligible kinetic/potential energy changes and consistent local properties:

    x_eq,out = (h_in - h_l,out)/h_fg,out + Q_net/(mdot h_fg,out)
    Q_net = eta q'' pi D L
    x_eq,in = (h_in - h_l)/h_fg
    delta_x = eta 4 q'' L / (G D h_fg)

The last two lines use a common nominal saturation state. The heat-loss factor eta applies if the quoted heat flux is heater input rather than already net heat to the fluid. D = 1 mm, heated L = 265 mm, R600a at 27 °C and nominal q'' come from the matching campaign. The 366 mm total length is not the heated length. Equating heated-section outlet state with camera-station state also neglects intervening heat/pressure changes. Per-run inlet enthalpy, net heat loss and absolute comparison-station pressure remain unavailable. A differential pressure signal cannot substitute for absolute pressure in the equation of state.

The nominal baseline reproduces v12/v12.1 with x_eq,in = 0 and eta = 1: saturated liquid at the heated inlet. It is an **upper bound only conditional on a subcooled inlet and appropriate net-heat assumptions**. A two-phase inlet changes the opposite way. Pattern labels and literature alpha ranges do not identify the missing enthalpy.

For x_eq outside [0,1], the scripts preserve the value and flag the two-phase equilibrium calculation as out of domain; they do not clip it into an apparently valid phase velocity. Negative equilibrium quality does not prove absence of vapor in nonequilibrium subcooled boiling. At exactly x = 0 or 1, the absent phase's mean velocity is undefined.

At zero inlet quality and fixed properties:

    Jg = eta 4 q'' L/(D h_fg rho_g)
    Jl = [G - eta 4 q'' L/(D h_fg)]/rho_l

Changing the mass-flow source changes x and Jl but cancels out of Jg. Identical Jg at equal q'' is imposed by this energy model. With nonzero inlet quality, Jg also includes G x_eq,in / rho_g.

## 3. Recorded flow improves provenance, with a mapping limit

The [mass-flow audit](massflow_review.md) follows the documented g/s channel through src/data_processed to original LabVIEW spreadsheets and verifies all samples. The src2 v9/v11/v13 workbooks have zero placeholders in Voltage_0/1; these are not flow measurements. Camera runs and LabVIEW records share condition names, but the files do not prove contemporaneous pairing.

We report a **nominal-G baseline** for comparison with the historical model and a **recorded-mean-flow alternative** with better mass-flow provenance for condition-level analysis. The latter remains conditional on name-based matching, nominal geometry, unknown inlet state and heater assumptions; it is not a synchronized camera measurement. In particular, the suffix-free m05q05 sensor blocks cannot distinguish retained camera _001 from excluded _002. Their aggregate is a sensitivity input for the shared condition, not evidence specifically assigned to _001. No excluded camera images enter any calculation.

Median and 10% trimmed-mean alternatives assess aggregation sensitivity. Original samples are correlated records, not hundreds of thousands of independent experiments. Their spread and small mean/median differences do not quantify instrument uncertainty or camera matching error.

## 4. Void-fraction alternatives

All four models receive the same G, x and properties in each scenario and cannot alter superficial velocities. The [literature review](literature_review.md) provides sources, page evidence and the Xu–Fang exponent discrepancy.

| Tested model | Reason to test it | Interpretation |
|---|---|---|
| Homogeneous | Transparent no-slip benchmark | Forces ug = ul; not a measured finding |
| Rouhani–Axelsson (RA) | Correctly coded in v12.1 and explicitly printed in Dias 2025 | Provisional reporting reference for continuity; no demonstrated accuracy advantage |
| Xu–Fang, literal de Oliveira 2021 Eq. 5 | Closest campaign prints exponent −2 | Source-reproduction diagnostic; not silently treated as canonical Xu–Fang |
| Xu–Fang, fluids implementation | Retrieved secondary implementation uses exponent −0.2 | Alternative with documented provenance; original-source confirmation unresolved |

Implemented expressions:

    alpha_h = (x/rho_g) / [x/rho_g + (1-x)/rho_l]
    C0 = 1 + 0.12 (1-x)
    Vgj = 1.18 (1-x) [g sigma (rho_l-rho_g)]^0.25 / sqrt(rho_l)
    ug_RA = C0 (Jg+Jl) + Vgj
    alpha_RA = Jg / ug_RA
    Fr_lo = G² / (rho_l² g D)
    alpha_XF = [1 + (1 + 2 Fr_lo^p alpha_h^3.5)
                       ((1-x)/x) (rho_g/rho_l)]^(-1)
    p = -2 (literal 2021), or -0.2 (retrieved fluids implementation)

Using the square root of this Fr_lo with the same exponent would change the model. The secondary implementation states Fr_lo = 0.02–145; nominal m06/m07 exceed this. Results retain Fr_lo and a range flag. The literal −2 form approaches the homogeneous model at our large Fr_lo; their agreement is algebraic, not independent corroboration.

RA is a **provisional reporting reference** because it reproduces v12.1 and its equation is legible in a relevant refrigerant study. It is not chosen for attractive clusters. RA's vertical-flow origin and application to related horizontal experiments do not establish its error here. A single best physical closure cannot be selected from these data.

## 5. Inlet, property and heat-input tests

The experiment table includes zero inlet quality; 1/2/5/10 K inlet subcooling at nominal saturation pressure; positive inlet qualities 0.01/0.05; saturation temperatures 25/29 °C; 90% net heat; and three recorded-flow aggregates. CoolProp evaluates subcooled enthalpy at the specified temperature/pressure. Temperature tests recalculate all saturation properties. Positive inlet-quality cases are hypothetical stress tests, not evidence of two-phase inlets in this campaign.

These are one-at-a-time assumptions, **not confidence intervals, sensor tolerances or an exhaustive physical envelope**. Inlet/temperature/heat tests use nominal G; recorded-flow alternatives use zero inlet quality and 27 °C. Interactions are not exhaustively tested. No assumption was fitted to the pattern labels.

## 6. Why v14 image alpha does not resolve phase velocity

The v14 alpha table is copied once into a hashed snapshot to prevent later regeneration silently changing v15. The diagnostic combines its spatially averaged circular-geometry estimate with the nominal energy baseline, testing 56/60/64-pixel diameter assumptions separately. Saturated m03q40 produces no hybrid velocity; it remains in the independent thermodynamic calculation.

This is a conditional diagnostic. The v14 geometry and crop ceiling remain unvalidated; panorama stitching uses image displacement; a spatial average may differ from a time mean at a fixed station; and a ratio of means is not generally the mean instantaneous velocity. Combining assumed x with this alpha can produce very large gas velocities or S < 1. Those results show sensitivity or incompatible assumptions; they do not identify the erroneous input, prove S < 1 impossible, or validate the plausible-looking cases.

The saved RA-implied x and x_in invert the model to match image alpha. They describe a required additional assumption, **not recovered measurements**, and are not used to tune results.

At fixed superficial flux, small void errors propagate approximately as:

    delta_ug/ug ≈ -delta_alpha/alpha
    delta_ul/ul ≈ +delta_alpha/(1-alpha)

Small gas fractions amplify ug uncertainty; near-unity gas fractions amplify ul uncertainty. Exact division cannot turn an uncalibrated mask into accurate phase velocities.

## 7. Recommended calculation and validation

Use the recorded-flow alternative alongside the nominal baseline for the current condition comparison, label each output with its actual inlet assumption, and retain RA plus other closures as sensitivity. The two baselines assume zero inlet quality; inlet-sensitivity experiments explicitly change that assumption. Use signed native-frame NCC separately for structure celerity, with acceptance coverage and alternatives shown. Do not convert panorama length or clamped stitching shifts into phase velocities.

To establish physical accuracy:

1. Match each camera run to its LabVIEW interval and verify flow calibration and actual tube diameter.
2. Recover inlet enthalpy or inlet temperature/pressure/preheater energy, net heat input and camera-station state; recompute x and superficial velocities.
3. Obtain independent cross-sectional void fraction with a defined averaging convention, or phase-resolved velocimetry. Validate optical geometry/scale if images supply these quantities.
4. Compare errors on matched acquisitions against independent measurements, including uncertainty and held-out conditions. Report failures and coverage by pattern; do not optimize class separation.

The user requested autonomous completion with repository evidence. Missing measurements are recorded as limitations and open questions. This investigation completes a reproducible comparison; physical accuracy remains an experimental question.
