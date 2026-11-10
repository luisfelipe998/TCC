# Literature and physical identifiability review

Reviewed 2026-10-23. Scope: single-label scenarios in `src2/data_processed`, excluding `m05q05_002`, smooth annular and compound classifications. This report distinguishes equations, source consistency, and experimental validation. A correct equation with assumed inputs is a conditional estimate; it does not establish the measured velocity.

## 1. What the velocity comparison means

For local total mass flux G, vapor mass quality x and phase densities rho_g/rho_l, the superficial velocities are

```text
j_g = G x / rho_g
j_l = G (1-x) / rho_l.
```

These are volumetric flow rates divided by the **whole** tube area. They are defined independently of slip and void fraction. Dias et al. (2025), [local PDF](../../../assets/dias2025-knn-flow-pattern.pdf), Section 2.2, PDF p.5, Eqs. (6)–(7), supplies these exact definitions. Although Table 2's caption calls them phase velocities, the equations and the surrounding text identify them as **superficial** velocities. Comparing that table directly to `u_g`/`u_l` from v12.1 would mix different quantities.

With cross-sectional gas void fraction alpha, phase-mean velocities are

```text
u_g = j_g / alpha
u_l = j_l / (1-alpha)
S = u_g/u_l = x rho_l (1-alpha) / [alpha rho_g (1-x)].
```

The last identity is de Oliveira et al. (2021), [local PDF](../../../assets/deoliveira2021-flow-patterns-pressure-drop.pdf), PDF p.14, Eq. (17); the first two follow from phase volumetric conservation. For time-dependent flow, `mean(j_g)/mean(alpha)` is a phase-occupancy-weighted velocity under constant-density assumptions, not generally `mean(j_g/alpha)`. Optical averages and thermodynamic averages must refer to the same station and averaging interval. With phase change or appreciable pressure variation, use local quality/densities consistently.

An interface, bubble nose or annular wave has a **structure propagation velocity**. Tracking one does not identify the liquid mean velocity, and an annular wave is not a gas material tracer. Optical flow satisfies image brightness/texture constraints; it does not on its own separate these physical velocities. A single measured structure speed must not be split into gas/liquid speeds merely by assigning bright/dark pixels to phases.

## 2. Thermodynamic inputs: what is supported and what remains unknown

de Oliveira (2021), PDF pp.3–4 and Table 1, describes R600a, Tsat = 27 °C, a horizontal 1 mm internal-diameter tube, 265 mm heated test length, G = 380/510/640/760/890 kg m^-2 s^-1 and q'' = 1/5/10/20/40/80 kW m^-2. These match the nominal map used by v12. They support using that map as an explicit **nominal campaign assumption**. Folder names are not calibrated measurements of G or heat input. The letter q encodes heat flux, **not vapor quality**.

The source's Eqs. (1)–(4), PDF pp.7–8, require preheater inlet enthalpy and power to obtain test-section inlet quality. For uniform heating, constant properties, negligible losses, and negligible kinetic/potential contributions, their first-law balance reduces to

```text
x_in  = [Q_PH/mdot + h_PH,in - h_l]/h_lv
x_out = x_in + 4 q'' L_heated/(G D h_lv).
```

Here a negative equilibrium quality is an enthalpy coordinate for subcooled liquid, not a negative physical gas mass fraction. A full balance with losses uses net heat input. If saturation pressure varies, calculate outlet enthalpy then quality with **outlet** saturated enthalpies; the simple constant-property difference must not silently absorb pressure effects.

Setting `x_in=0` is a useful, reproducible scenario. Section 2, PDF p.4, says that the preheater established a subcooled inlet to the test section in that campaign. If this condition is verified for each local record, and omitted heat losses are nonnegative, the no-loss saturated-inlet estimate is an upper bound on equilibrium outlet quality under the same property assumptions. Without the per-record inlet state and pressure, it is **not an unconditional upper bound**. An inlet already containing vapor would increase outlet quality; an inlet with subcooling would reduce it. Neither direction can be established from the folder names.

Clipping equilibrium x to [0,1] can be useful for flagging single-phase boundaries, but must not disguise an invalid two-phase state or a failed input assumption. For x=0, gas phase-mean velocity is not identified because that phase is absent; similarly liquid phase velocity is not identified at x=1.

**Correction to the earlier data barrier claim:** `src2/data_processed/CLASSIFICACAO_PADRAOES_MANUAL.md`, lines 2–17, explicitly records the data owner's channel identification: `Voltage_0` is pressure drop in kPa, `Voltage_1` is mass flow in g/s and `Voltage_2` is optical voltage. Therefore the v12.1 statement that the spreadsheet channel is necessarily uncalibrated raw voltage is incorrect. The [completed mass-flow audit](massflow_review.md) verifies preservation of the original LabVIEW values and reports condition-level G alternatives. The similarly named image-derived columns under `src2` are zero placeholders. The documented g/s meaning is resolved; instrument calibration uncertainty and the sensor-to-camera mapping remain unresolved. The earlier **unimplemented** claim of a 4–13% cross-check remains unsupported and does not describe the newly audited 26-condition results. A column name alone does not overturn the owner-supplied units. Differential pressure cannot replace the missing absolute outlet pressure. Recovering G also does not recover x_in, preheater power or inlet enthalpy.

For **m05q05**, the sensor blocks are stored under the condition name without the camera suffix. The five available `dados*.xlsx` blocks cannot be assigned independently to camera `_001` versus `_002` from the available metadata. Camera `m05q05_002` is excluded completely. The sensor-based alternative retained alongside `m05q05_001` is therefore explicitly a **condition aggregate**, not an identified `_001` sensor interval; it cannot establish which blocks, if any, coincide with the excluded camera acquisition. The nominal-G baseline avoids this sensor-pairing assumption but still requires the inlet/heat/property assumptions.

The local literature uses REFPROP v9.1 for properties; v12 uses CoolProp. This is a transparent software/property-model choice, not an experimental validation. Independent calculation here with CoolProp 8.0.0 at 300.15 K gives rho_l = 548.1311606 kg/m³, rho_g = 9.6497318 kg/m³, h_lv = 326707.5282 J/kg and sigma = 0.00978195984 N/m. No claim is made that differences between software packages explain source/table inconsistencies.

## 3. Which void-fraction closure is defensible

### Homogeneous baseline

```text
alpha_h = (x/rho_g)/[x/rho_g + (1-x)/rho_l]
u_g = u_l = j_g+j_l; S=1.
```

This is an equal-phase-velocity assumption. It is a useful baseline, not a velocity measurement or a generally reliable model of plug, slug, churn and annular slip. The use of superficial velocities alone in v12 does **not** impose this homogeneous assumption; the erroneous linkage between v12's superficial definitions and a no-slip model should be removed from interpretations.

### Rouhani–Axelsson form used by Dias (2025)

Dias (2025), PDF p.5, Eq. (5), visually checked in [rendered page](../literature/pages/dias2025knn-p05.png):

```text
C0 = 1 + 0.12 (1-x)
drift_term = 1.18 (1-x) [g sigma (rho_l-rho_g)]^0.25 / (G sqrt(rho_l))
alpha_RA = (x/rho_g) / { C0 [x/rho_g + (1-x)/rho_l] + drift_term }.
```

C0 multiplies **both** mixture specific-volume terms. `drift_term` in this form has units m³/kg; multiplying it by G gives the drift contribution in m/s. v12.1's formula matches this source. The article says that the model was initially designed for vertical flow and cites applications to horizontal hydrocarbon mini/microchannels. This makes it a reasonable named engineering estimate in this campaign, with model-form uncertainty. It does not prove that it is the most accurate closure for these images.

The legacy assertion that de Oliveira (2021) also used this RA correlation is incorrect: that paper explicitly used Xu–Fang. Further, RA is a **void-fraction** closure, not a correction to x or to superficial velocities. Applying RA can neither fix unknown inlet enthalpy nor make the nominal energy balance accurate.

### Xu–Fang: source conflict that must stay visible

de Oliveira (2021), PDF p.8, Eq. (5), prints

```text
alpha_printed = [1 + (1 + 2 Fr_lo^(-2) alpha_h^3.5)
                        ((1-x)/x) (rho_g/rho_l)]^(-1).
```

The exponent is visibly **-2**, not -0.2, in a [400-dpi equation crop](../literature/pages/deoliveira2021-eq05-highres.png). That page calls Fr_lo the liquid-only Froude number, without an explicit definition. Thus copying the displayed equation does not by itself fully specify a reproducer.

A publicly inspectable secondary implementation, Caleb Bell's `fluids`, provides **-0.2** and explicitly defines

```text
Fr_lo = G²/(g D rho_l²).
```

See [`Xu_Fang_voidage`, fixed source revision b6b0a04](https://github.com/CalebBell/fluids/blob/b6b0a04d3ecfeb26894a6890fcea90f0e3436057/fluids/two_phase_voidage.py#L1858), lines 1863–1907, also preserved in [local source bytes](../literature/extracts/fluids_two_phase_voidage.py). It cites Xu and Fang (2014), *Applied Thermal Engineering* 64, 242–251, DOI [10.1016/j.applthermaleng.2013.12.032](https://doi.org/10.1016/j.applthermaleng.2013.12.032). The implementation's documented domain is Fr_lo = 0.02–145, rho_g/rho_l = 0.004–0.153, x = 0–1; it reports an AARD of 5% from that source. These domain/error statements are **secondary-source claims**, not a verified uncertainty for the present study. In particular, no 5% error bar should be attached to our predictions. The direct expression has singular endpoints and needs defined limits if x=0 or 1.

At nominal G and 27 °C, the resulting Fr_lo values are:

| G [kg m^-2 s^-1] | Fr_lo | Relative to documented upper limit 145 |
|---:|---:|---|
| 380 | 48.9925 | Within |
| 510 | 88.2475 | Within |
| 640 | 138.9703 | Within |
| 760 | 195.9698 | Outside |
| 890 | 268.7460 | Outside |

The density ratio is 0.0176048, within the documented density-ratio interval. Recompute the flags when using channel-derived G or other temperatures. The assumed x values can lie within [0,1] without proving their physical correctness. The original Xu–Fang PDF is absent from `/assets`. An [OpenAlex lookup](../literature/extracts/xufang2014_openalex.json) returned no open-access/repository full-text location. The original article was **not** independently read in this review.

Consequently test two distinctly named variants: **de Oliveira printed exponent -2** and **Xu–Fang secondary implementation exponent -0.2**. Their difference is a source/formulation sensitivity, not a fitted parameter. The -2 variant will approach the homogeneous result at these Fr values; that numerical property is not proof that the original Xu–Fang correlation has negligible slip. Do not silently select -0.2 as a verified correction, or call either variant the best solely because it agrees with a code library or a desired pattern ordering. Obtain the original article or author clarification before making an authoritative attribution.

### Horizontal intermittent drift-flux fits

Zeghloul and Al-Sarkhi (2023), [local PDF](../../../assets/driftflux2023-void-horizontal-intermittent.pdf), Abstract and Sections 2–3, describe plug/slug data in 19–152.4 mm pipes, including air–water and viscous oil systems up to 600 cP. The reported approximately 16% mean relative absolute error applies to that database. Those geometry/fluid/regime conditions do not validate a 1 mm boiling-R600a application. It is useful for the relations `u_g=C0(j_g+j_l)+Vgr` and `alpha=j_g/u_g`, but using fixed fitted C0 values for every local regime would add unsupported assumptions. Do not transfer those error percentages to our runs.

## 4. Why matching published regime ranges is not an accuracy test

The source papers calculate alpha with a correlation. Agreement between a correlation and another correlation-derived source range is not independent measurement validation. A regime label is many-to-one in operating conditions and does not supply a unique x, alpha, slip ratio or velocity. The classes in these figures are not universal, nonoverlapping calibration intervals.

There is also a concrete internal source inconsistency. de Oliveira (2021), Section 4.4, PDF p.17, reports churn x = 0.35–0.62, alpha = 0.61–0.85 and S = 1.40–2.34. Its own Eq. (17), using 27 °C R600a densities, implies alpha = **0.92893–0.98512** from the printed x/S ranges. This interval is disjoint from 0.61–0.85. Conversely the printed x/alpha ranges require **S >= 5.39754**, above 2.34. These are bounds over the entire reported ranges, so the contradiction does not depend on knowing which values were paired. A rounding explanation or a corrected decimal cannot be justified from the local source.

The v12.1 report's claim that de Oliveira (2021) has no slug alpha is also false: Section 4.3, PDF p.16, explicitly states alpha approximately 0.65 for the plug–slug transition. This is a transition statement, not a universal slug interval. Figure-specific thesis examples likewise should not become strict acceptance bounds on all slug scenarios.

Dias (2025), Table 1, PDF p.4, describes a **different** campaign: G = 516, q'' = 10/20/30/40/60, preheater power 2–58 W, nominal 371 kPa and 27 °C. Table 2, PDF p.5, consists of selected examples. Independently substituting its printed x into its stated equations at G=516 and 27 °C gives:

| Selected example | Printed x | Printed alpha | RA alpha | Printed j_g [m/s] | Calculated j_g [m/s] |
|---|---:|---:|---:|---:|---:|
| Bubbles, small | 0.012 | 0.36 | 0.33922 | 0.7 | 0.64168 |
| Bubbles, larger | 0.025 | 0.50 | 0.50442 | 1.8 | 1.33682 |
| Plug | 0.056 | 0.68 | 0.67283 | 2.4 | 2.99449 |
| Slug | 0.076 | 0.73 | 0.72497 | 3.5 | 4.06395 |
| Churn | 0.220 | 0.85 | 0.85410 | 10.5 | 11.76406 |
| Wavy annular | 0.300 | 0.88 | 0.88150 | 11.7 | 16.04190 |

This is a **source consistency diagnostic**, with the explicit CoolProp substitution; exact original REFPROP properties and unrounded source inputs were not supplied. Some alpha entries are close; several j_g entries do not reproduce. The source's actual local thermodynamic states and any table error are unresolved. Do not assign these differences to rounding, slip or software without evidence, and do not minimize prediction error to this table as if it were local ground truth. Full independent arithmetic and j_l checks are in [source_numeric_checks.json](../literature/source_numeric_checks.json), generated by [source_numeric_checks.py](../literature/source_numeric_checks.py).

The Dias (2024) dissertation text extraction was attempted and is not readable because of its encoded text layer. No thesis numeric annotation is newly adopted as ground truth in this review. Previously transcribed points remain selected points from another record set, even if visually legible.

## 5. Camera timing, calibration and available optical observables

The native per-run `_config.xsv` is the available timing authority. For example, `m05q01_001/_config.xsv` records Rate=7000, Images=501, StartIdx=0, StopIdx=500, Step=1, ROIW=1016 and ROIH=80. Using its `Output_Rate=5` as capture rate would be wrong. Native runs have varying rates; a campaign-level 5000-fps statement from de Oliveira (2021), or 3000 fps from Dias (2025), must not overwrite acquisition metadata. The endpoint-to-endpoint duration is `(N-1)/fps`. Fixed-station tracking and consecutive-frame displacements must use actual frame indices and native time steps. Resampled 4000-Hz spreadsheet signals and stitched panoramas are not substitutes for native video timing.

`ScaleFactor=1.0` in the configs is attached to voltage-channel scaling, not documented camera mm/pixel calibration. The legacy claim that this field proves unit spatial scale is therefore incorrect. Conversely, the absence of a usable spatial calibration still prevents unconditional conversion of image displacement to m/s. Dias (2025), PDF p.4, reports **13.68 micrometres/pixel** for its own 3000-fps, 1016×88-pixel video dataset (p.5). It cannot be transferred automatically to these 1016×80, 501-frame, mixed-rate recordings.

It is legitimate to report a repeatable optical structure velocity in px/frame or px/s. A conditional conversion based on visible full internal diameter D and isotropic optical magnification can be explored if the true internal walls are identified and optical distortion/refraction is quantified. The 80-pixel ROI height is not automatically a 1 mm tube-diameter ruler. Full-field camera calibration or an in-situ ruler is preferable. Registration, axial strips, central and near-interface crops, forward/reverse consistency, multi-frame lag agreement, and synthetic known shifts can test **algorithmic displacement accuracy**; they do not prove phase-mean physical accuracy.

Even a validated projected gas mask is not volumetric alpha in a circular tube without geometry. Under explicitly tested axisymmetry and corrected interface locations, `alpha(z,t)=(d_g(z,t)/D)^2`; integrate local squared diameters, rather than squaring a global mean width. Churn, occlusion, separated bubbles, wall films, lens/refraction effects and saturation can invalidate this construction. Therefore inserting a panorama occupancy into `j_g/alpha` does not turn v12.1 into a measured phase velocity. A panorama may also overrepresent slow or fast structures depending on its reconstruction; use matching native time averages for comparisons.

## 6. Recommended interpretation and missing information

The best-supported current reporting hierarchy is:

1. **Primary comparison:** reproducible j_g × j_l with a clearly labeled assumed inlet state and nominal or owner-labeled channel G, with sensitivity to the inputs that are actually uncertain. No slip closure is needed for these coordinates.
2. **Conditional phase comparison:** RA-derived u_g × u_l as a literature-grounded reference model; homogeneous and the two clearly distinguished XF forms quantify model dependence. Report results as predicted phase means, not measured or independently validated true velocities.
3. **Independent optical comparison:** observed structure displacements with native timing, quality flags and displacement sensitivity. Keep px/s until scale is established; if a conditional scale is used, state it explicitly. Compare the same observable (e.g. bubble translation) across optical methods. A closer fit to predicted j_g or u_g is not evidence that an optical algorithm is physically correct.

No available source supports declaring one method the most accurate in m/s for both phases without an independent velocity reference. Sensitivity ranges, algorithm disagreement and model spread must not be labeled confidence intervals on the true velocity unless a justified measurement-error model is supplied. Independent experimental validation would require matched calibrated flow/enthalpy measurements and validated volumetric alpha, or suitable phase-resolved velocimetry with uncertainty.

Questions that materially limit absolute accuracy remain:

- Can original test/preheater power, inlet temperature/enthalpy and **absolute** pressure be mapped to each recording? Does the verified subcooled-inlet condition hold for every run?
- What are the calibration uncertainty and camera-pairing rule for the owner-confirmed g/s mass-flow channel? In particular, can the m05q05 sensor blocks be assigned to `_001` versus the excluded camera `_002`?
- Is there an original in-situ axial pixel scale or an image containing a reliable geometric reference?
- Can the original Xu–Fang article or source authors resolve the -2 versus -0.2 conflict?
- Is there a matched independent phase velocity or volumetric void-fraction measurement against which to quantify absolute error?

## 7. Evidence provenance

PDF equations reviewed directly: Dias (2025) Eq. (5) and Table 2 on p.5; de Oliveira (2021) Eq. (5) on p.8, including a 400-dpi crop. Text checks: de Oliveira (2021) pp.3–4, 7–8, 14, 16–17; Dias (2025) pp.4–5; Zeghloul and Al-Sarkhi (2023) abstract and Sections 2–3. Native config and owner channel notes were also read directly. Text extractions in [literature/extracts](../literature/extracts) were produced with `pdftotext -layout`; page images with `pdftoppm`. The public `fluids` source is pinned to commit b6b0a04d3ecfeb26894a6890fcea90f0e3436057; retrieval metadata and SHA-256 hashes are preserved. The external code was inspected as text, not imported to generate this report's numbers.
