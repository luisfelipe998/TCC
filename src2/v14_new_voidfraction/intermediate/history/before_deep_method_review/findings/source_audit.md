# Source audit: void fraction and flip-rate

<!-- current-dataset-scope:start -->
**Current dataset scope:** 26 retained acquisitions in five single patterns; annular means wavy annular. Mixed labels, smooth annular and m05q05_002 are excluded from current calculations. Published literature values below remain source evidence, not current dataset measurements.
<!-- current-dataset-scope:end -->
Reviewed directly from `/assets` PDFs on 2026-10-23. Equations and image-table entries below were checked against rendered PDF pages, not copied from earlier project memos. Dataset evidence comes only from `src2`.

## 1. What the closest literature actually establishes

**de Oliveira et al. (2021), IJMF 144, 103793**, [local PDF](../../../../../../assets/deoliveira2021-flow-patterns-pressure-drop.pdf), DOI [10.1016/j.ijmultiphaseflow.2021.103793](https://doi.org/10.1016/j.ijmultiphaseflow.2021.103793). This is the closest operating-condition reference: Table 1, PDF p.4, lists R600a, 27 °C, G = 380/510/640/760/890 kg m⁻² s⁻¹ and q″ = 1/5/10/20/40/80 kW m⁻². Section 2 reports 1 mm internal diameter and 265 mm heated test length. Its quoted 5000 fps camera rate is not the per-file timing authority for `src2`; the actual acquisition configs differ.

Its void fractions are **calculated with a correlation**, not obtained by thresholding high-speed images. Section 3.1, PDF p.8: “The void fraction was evaluated with the correlation proposed by Xu and Fang (2014)”. Eq. (5), visually checked:

```text
alpha_XF = [1 + (1 + 2 Fr_lo^(-2) alpha_h^3.5)
                  * ((1-x)/x) * (rho_v/rho_l)]^(-1)
alpha_h  = [1 + ((1-x)/x) * (rho_v/rho_l)]^(-1)
```

The paper names Fr_lo as the liquid-only Froude number but does not supply its definition; using `Fr_lo = G²/(rho_l² g D)` is the conventional closure and must be declared if this equation is implemented. The homogeneous expression for alpha_h follows from equal phase velocities; it is not an optical calibration. At the large nominal G and 1 mm diameter here, this closure makes XF very close to the homogeneous result.

Eqs. (1)–(4), PDF pp.7–8, show why G and q″ do not determine the outlet quality alone:

```text
x_in  = (Q_PH/mdot + h_PH,in - h_l)/h_lv
x_out = x_in + 4 q″ L_heated/(G D h_lv)
```

The second equation follows algebraically from the first-law balance and cylindrical heated-area/cross-sectional-area ratio. It assumes negligible losses and properties evaluated consistently. Setting x_in = 0 is a scenario. Section 2, PDF p.4, says the PH establishes a “subcooled condition in the inlet of the TS”; if that description applies to a particular record, x_in <= 0 and the zero-inlet result is an upper bound on equilibrium quality. Without proving that inlet condition for each run, it is not an unconditional upper bound.

**Dias et al. (2025), IJMF 185, 105122**, [local PDF](../../../../../../assets/dias2025-knn-flow-pattern.pdf), DOI [10.1016/j.ijmultiphaseflow.2024.105122](https://doi.org/10.1016/j.ijmultiphaseflow.2024.105122). Table 1, PDF p.4, is a different campaign: G = 516, q″ = 10/20/30/40/60, PH power 2–58 W, 3000 fps, 13.68 µm/pixel, 371 kPa and 27 °C. These pixel/fps numbers cannot be transferred to the `src2` runs. Its Table 2 is a few selected images at different operating points, not a universal lookup of void fraction by regime:

| Selected regime | x | alpha |
|---|---:|---:|
| Bubbles, small | 0.012 | 0.36 |
| Bubbles, larger | 0.025 | 0.50 |
| Plug | 0.056 | 0.68 |
| Slug | 0.076 | 0.73 |
| Churn | 0.22 | 0.85 |
| Wavy annular | 0.30 | 0.88 |

Section 2.2, PDF p.5: “The void fraction was determined utilizing the model established by Rouhani and Axelsson (1970)”. Eq. (5), visually checked:

```text
C0 = 1 + 0.12 (1-x)
drift = 1.18 (1-x) [g sigma (rho_l-rho_v)]^0.25/(G sqrt(rho_l))
alpha_RA = (x/rho_v) / [C0 (x/rho_v + (1-x)/rho_l) + drift]
```

The C0 multiplies **both** mixture-volume terms; the existing `src2/v12.1/velocity_map_v12_1.py:217` implementation agrees with the rendered equation. The article explicitly says the original model was designed for vertical flow and cites use in horizontal hydrocarbon mini/microchannels. It is a reasonable alternative model, with model-form uncertainty, not proof of the true alpha of our records. Its inlet-quality Eq. (2), PDF p.4, has a typographical inconsistency; use the internally consistent first-law balance in de Oliveira (2021), not its literal malformed expression.

## 2. Closest-source numeric contradictions prevent using its image tables as truth

The following annotations are legible directly in de Oliveira (2021) figures. Nominal G and q″ match some `src2` folder names, but equal setpoints do not prove identical acquisition/thermodynamic conditions. More importantly, the printed x/alpha pairs fail the same paper's physical consistency relation.

| Source | G | q″ | Printed x | Printed alpha |
|---|---:|---:|---:|---:|
| Fig.7, PDF p.6, plug | 640 | 1 | 0.094 | 0.042 |
| Fig.7, PDF p.6, plug | 760 | 1 | 0.040 | 0.290 |
| Fig.7, PDF p.6, plug | 890 | 5 | 0.690 | 0.240 |
| Fig.12, PDF p.9, slug | 890 | 10 | 0.260 | 0.650 |
| Fig.12, PDF p.9, slug | 760 | 5 | 0.230 | 0.610 |
| Fig.16, PDF p.12, churn | 380 | 5 | 0.420 | 0.850 |
| Fig.16, PDF p.12, churn | 640 | 20 | 0.410 | 0.820 |
| Fig.16, PDF p.12, churn | 890 | 40 | 0.560 | 0.780 |
| Fig.20, PDF p.14, wavy annular | 380 | 20 | 0.680 | 0.860 |
| Fig.20, PDF p.14, wavy annular | 510 | 40 | 0.750 | 0.930 |
| Fig.20, PDF p.14, wavy annular | 760 | 80 | 0.760 | 0.960 |

Eq. (17), PDF p.14, implies the identity

```text
S = x rho_l (1-alpha)/[alpha rho_v (1-x)]
alpha = [1 + S ((1-x)/x) (rho_v/rho_l)]^(-1).
```

Using R600a at the stated 27 °C (`CoolProp`: rho_l = 548.13116 and rho_v = 9.64973 kg/m³), Section 4.4, PDF p.17, simultaneously reports x = 0.35–0.62, alpha = 0.61–0.85, and S = 1.40–2.34. The x/S ranges imply **alpha = 0.92893–0.98512**, disjoint from 0.61–0.85. Conversely the x/alpha ranges imply **S >= 5.39754**, already above its maximum 2.34. Fig.7's x=0.69, alpha=0.24 implies S=400.37, whereas Section 4.2 gives plug S=1.24–1.34. This is not a subtle choice of correlation or a rounding error.

Do not silently change a decimal, fit to these numbers, or claim they validate our absolute alpha. The inconsistency is in the published local PDF and its cause cannot be identified from our data. The standalone alpha ranges may still be quoted as **reported source values**, with this qualification. Also, the earlier project memo's assertion that this article provides no slug alpha is false: Fig.12 gives 0.61/0.65 and Section 4.3, PDF p.16, states alpha approximately 0.65 at the plug–slug transition.

## 3. Correct optical alternative: segment the enclosed phase, then state geometry

**Schepperle, Junaid and Woias (2024), Sensors 24, 3363**, [local PDF](../../../../../../assets/schepperle2024-cv-void-fraction-microchannel.pdf), DOI [10.3390/s24113363](https://doi.org/10.3390/s24113363). Section 2.1, PDF p.3, describes a **rectangular 0.5 mm deep × 1.5 mm wide** DI-water microchannel, unlike the circular tube here. Section 2.3, PDF pp.5–7, uses adaptive thresholding, blur/gamma correction, **contour filling and comparison to a liquid-only reference image**. It identifies liquid contours by histogram correlation >=0.85, then sets gas masks to black and liquid masks to white. That black/white convention is assigned after segmentation; it does **not** mean naturally dark raw pixels are vapor. Eq. (3), PDF p.7, counts mask vapor pixels divided by subsection pixels. Section 2.4, pp.7–8, manually corrects erroneous masks before U-net training. The reported Dice score quantifies agreement with masks, not independent validation of 3D void volume.

This supports a concrete improvement over bright-core thresholding: locate outer gas-interface contours and include their enclosed interiors; separate illumination changes using matched liquid references; inspect masks across plug, slug, churn and annular examples. It does not justify copying a threshold, a U-net trained on another apparatus, or the rectangular projection-to-volume assumption.

For a circular tube, even a perfect side-view gas mask measures a **projected** area fraction. If a vapor core is demonstrably axisymmetric and its corrected local diameter is d_g(z,t), elementary geometry yields

```text
alpha(z,t) = (d_g(z,t)/D)^2
alpha_volume = mean_z [(d_g/D)^2].
```

The ratio is dimensionless, so missing mm/pixel alone does not make an estimate impossible. The real requirements are identifying the full internal tube diameter, optical-distortion correction and testing axisymmetry. For uniform annular film thickness delta, alpha = (1-2 delta/D)^2. For plug/slug tails, integrate the squared local diameter rather than squaring mean projected width. Without axisymmetry or independent depth information, many 3D shapes give the same side-view mask. Churn's occlusion and irregular geometry make a single-view volumetric alpha particularly underdetermined. The current bright central band is not d_g; squaring or rescaling its 0.25–0.35 values does not repair it. These are geometric derivations, not quoted equations from the papers.

**de Oliveira et al. (2019), IJMF 114, 240–257**, [local PDF](../../../../../../assets/deoliveira2019-flow-boiling-hydrocarbons-nonlinear.pdf), DOI [10.1016/j.ijmultiphaseflow.2019.03.015](https://doi.org/10.1016/j.ijmultiphaseflow.2019.03.015), Section 3.1, PDF p.4 (journal p.243), supports exploring radial geometry: horizontal mini/microchannel boiling has “practically the same radial symmetry in interfacial topology” as vertical flow “due to the dominance of capillarity”. This is apparatus-level motivation, not a check of every churn cross-section. Its Section 4, PDF p.9 (journal p.248), also states that its beta parameter was introduced to avoid the “overlapping of flow patterns when presented as function of vapor quality or void fraction”. That directly contradicts treating a strict nonoverlapping alpha-by-regime ladder as a required physical ground truth. This 2019 campaign uses 17 °C in its conclusions and has different operating conditions, so its numeric regime points cannot be substituted for `src2` alpha.

## 4. Flip-rate is not a universal regime ladder

de Oliveira (2021), Sections 4.2–4.5, PDF pp.16–17, separates whole-slug passage from interface waves and dispersed bubbles. Plug PD frequency is 0.5–1 Hz; its OS varies at 2–16 Hz. For slug, the article assigns about 2.8 Hz to slug passage, and 8/12/23/32 Hz to dispersed bubbles between slugs. Churn OS reaches 15–28 Hz with higher mass flux. Annular OS still has 2–16 Hz and approximately 20 Hz wave fronts. Thus a high threshold-crossing rate can come from annular waves, small bubbles, turbulence, brightness or stitching. The papers do not establish `plug > slug > churn > annular` for an arbitrary Schmitt-trigger signal. Section 4.3 attributes plug-to-slug development to coalescence, which itself can reduce event count while increasing gas occupancy.

For a binary, explicitly defined gas state b_j, the descriptive spatial flip count is

```text
N_flip = sum_j |b_j - b_(j-1)|
r_flip,px = 1000 N_flip/(N-1) [transitions/1000 px].
```

Count **both interfaces**. A complete isolated gas segment normally contributes two flips. Blob count, gas-entry count and total flip count differ; a near-equality of blob and flip counts is not validation. Using N instead of N-1 is a small finite-record convention but must be named. A phase occupancy/duty fraction, mean(b_j), is another observable; it is not interchangeable with flip-rate. Low count does not establish a stationary or nonchaotic phase state.

To claim **temporal** phase flips, form b(t) at a fixed axial station in the original frames and use actual native fps. Use multiple fixed stations for sensitivity, with uncertainty from gate and station choice; do not treat stations or panorama pixels as independent replicates. Native config evidence for the current 26 retained acquisitions: frame counts [501]; acquisition counts by native fps {4000: 2, 5000: 6, 7000: 18}. The actual endpoint-to-endpoint record duration is (N−1)/fps. `__timeseries_v9.xlsx` files are already resampled to 4000 Hz, with T_s=N/fs; those metadata are not native frame rates. `ScaleFactor=1` occurs in disabled voltage-acquisition sections with Units=[V], and is **not** a measured camera spatial calibration.

The finite record gives an approximately 1/T = **8–14 Hz** Fourier-bin spacing and covers only 0.035–0.125 cycles of a 0.5–1 Hz process or 0.20–0.35 cycles at 2.8 Hz. A gate can still count any visible fixed-station crossings within the clip, but the count/T is a short-record rate estimate, not a well-established long-run slug frequency. A panorama provides extra spatial detail under reconstruction assumptions; it does not lengthen elapsed acquisition time.

## 5. Additional alternative and practical conclusion

**Zeghloul and Al-Sarkhi (2023), JAFM 16(7), 1499–1510**, [local PDF](../../../../../../assets/driftflux2023-void-horizontal-intermittent.pdf), DOI [10.47176/jafm.16.07.1693](https://doi.org/10.47176/jafm.16.07.1693). Abstract and Sections 1–2 describe a horizontal plug/slug drift-flux fit covering **19–152 mm** pipes and fluids up to 600 cP, with approximately 16% mean relative absolute error. This is not a validated R600a 1 mm boiling-tube correlation. It supplies a reason to compare model forms, not a missing exact alpha.

For the requested panorama-based analysis, the defensible solution is reproducible **optical** mask/occupancy and spatial flip-rate with named gate, tube extent and finite-record sensitivity. The best improvement is contour detection with enclosed-region filling and a liquid-reference comparison, checked on masks across the individual panoramas. Where symmetry is visually supported, also report the **conditional geometric estimate** `mean[(d_g/D)^2]`, with its assumptions and interface-location sensitivity. A gate chosen because it recreates the regime ladder is circular validation. Actual volumetric alpha needs validated geometry; correlation/reference numbers do not cure an unvalidated optical mask. The available papers contain useful methods, but none supplies trustworthy ground-truth alpha and a universal flip-rate for the retained `src2` records. The thermodynamic equations above audit the literature claims and are not a proposal to replace the user's panorama measurements with model values.
