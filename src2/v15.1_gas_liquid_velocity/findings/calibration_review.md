# Calibration and physical meaning of the velocity axes

**Recommended main units: px/s for the measured image displacement, with dimensionless within-recording ratios.** Native timing is documented. No object-space axial metres-per-pixel calibration was found in the inspected camera files, native JPEG metadata or processed-workbook metadata. A conversion based on the 1 mm source tube diameter and the owner-confirmed 53-pixel lumen is possible as an explicitly conditional geometric calculation. It does not establish measured m/s or turn interface motion into superficial velocities.

## Evidence inspected afresh

The audit followed the existing 26-recording selection: mixed classifications, smooth annular and m05q05_002 remain excluded. It read all 78 original camera sidecars (_config.xsv, _config.mii and _config.mme), metadata from all 13,026 native JPEGs and all 104 processed workbooks. It checked frame-name sequences and independently parsed both native FPS fields. It also extracted the original de Oliveira and Dias PDFs again and visually inspected the four cited rendered pages.

| Evidence | Result for the 26 retained recordings | Consequence |
|---|---|---|
| Native Acquisition/Rate and Film speed 1 | 2 recordings at 4,000 fps; 6 at 5,000 fps; 18 at 7,000 fps; both fields agree | Use the native rate for each run. |
| Frame identities | 501 consecutive files, indices 0–500, Step=1 | 500 native intervals; durations 0.125, 0.100 or 0.0714286 s. |
| Camera identity and image dimensions | One camera ID, 19-0111-0149; all native images 1016 × 80 pixels; ROIY=208 | Common acquisition geometry is documented; constant optical magnification is not independently established. |
| Pixel aspect ratio and binning | Aspect ratio 1; Binning=0 in every run | Square image sampling supports an isotropic sensor grid, not equal radial/axial magnification through the tube. |
| Lens type, focal length, focus, reference system, distortion index and correction file | NOVALUE in all 26 .mii files; Movie images corrected=NO | No stored lens/distortion calibration was available. |
| Native JPEG metadata | All 13,026 have no EXIF entries; JFIF unit=0 and density=(1,1) | The JFIF fields express aspect information with no physical density unit. They supply no object-space ruler. |
| ScaleFactor=1 in .xsv | All 16 DaChannel sections per run are disabled and have Units=[V] | These are voltage-channel scales, not camera m/pixel. |
| Output_Rate | 5 in every run | This separate field is not the native capture rate; it conflicts with both native FPS sources. |
| Processed spreadsheets | 78 spatial workbooks specify 1 px/sample; 26 resampled workbooks specify fs_hz=4000 | The processed sheets add no spatial calibration; their time grid must not replace native camera timing. |

The original representative files are [_config.mii](../../data_processed/m05q01_001/_config.mii), [_config.xsv](../../data_processed/m05q01_001/_config.xsv) and [_config.mme](../../data_processed/m05q01_001/_config.mme). Relevant .mii evidence is at lines 8–11, 13–14, 18 and 25–37; native Rate and Output_Rate are at .xsv lines 39–41, exposure at line 48, ROI at lines 63–66, and voltage-channel scaling begins at lines 313–317. The [per-run inventory](../results/calibration_inventory.csv), [native image metadata](../results/calibration_image_metadata.csv), [workbook metadata](../results/calibration_workbook_metadata.csv) and [full provenance](../results/calibration_audit.json) cover every retained run.

Exposure is 62, 127 or 132 microseconds, corresponding to 43.4%–66.0% of a frame period across these runs. Matching consistency does not remove motion blur or shape evolution during that integration. The camera clock itself was not independently calibrated; FPS is documented acquisition metadata, not a certified uncertainty bound.

## What the two source papers support

**de Oliveira et al. (2021), International Journal of Multiphase Flow 144, 103793**, DOI 10.1016/j.ijmultiphaseflow.2021.103793:

- Section 2, **PDF p.3**, specifies a 1.0 mm internal diameter for the preheater; **p.4** specifies that the test section and visualization section have the same internal diameter. The test section heated length is 265 mm. The 136 mm visualization tube length is the physical section length; the article does not identify it as the camera's full image width.
- **PDF p.4**, Table 1, gives R600a, 27 °C, G=380/510/640/760/890 kg/(m² s) and heat fluxes 1/5/10/20/40/80 kW/m². These conditions closely match the existing folder convention. Equal setpoints and a related apparatus do not by themselves identify a particular camera acquisition or its magnification.
- Section 2 on **p.4** reports a 5,000-fps MotionPro Y4S1 camera and says the preheater establishes a subcooled test-section inlet. The present native camera files contain three different rates. The article supplies no transferable per-recording axial pixel ruler or inlet enthalpy for these files.

[Original PDF](../../../assets/deoliveira2021-flow-patterns-pressure-drop.pdf), [fresh text extraction](../intermediate/calibration_literature/deoliveira2021.txt), [rendered p.3](../intermediate/calibration_literature/deoliveira2021_p3.png), [rendered p.4](../intermediate/calibration_literature/deoliveira2021_p4.png).

**Dias et al. (2025), International Journal of Multiphase Flow 185, 105122**, DOI 10.1016/j.ijmultiphaseflow.2024.105122:

- Section 2.1, **PDF p.4**, states “3000 frames/s with a pixel length of 13.68 μm”. This is in the apparatus prose, not Table 1. Table 1 instead gives G=516 kg/(m² s), heat fluxes 10/20/30/40/60 kW/m² and preheater power 2–58 W.
- Section 2.3, **PDF p.5**, describes exactly 200 frames per video at **1016 × 88 pixels**, versus the present **501 frames at 1016 × 80**, with 4,000–7,000 fps. Neither the same camera family nor the same nominal tube diameter verifies an unchanged optical scale. Copying 13.68 μm/pixel is therefore unsupported.
- Section 2.2, **PDF p.5**, computes void fraction with Rouhani–Axelsson, then defines liquid and vapor superficial velocities by **Eqs. (6)–(7)**: j_l=G(1−x)/ρ_l and j_v=Gx/ρ_v. Table 2 lists results of that thermodynamic reduction. The values are not two velocities tracked independently in the videos. The table caption loosely calls them phase velocities; the equations establish the precise superficial-velocity meaning.

[Original PDF](../../../assets/dias2025-knn-flow-pattern.pdf), [fresh text extraction](../intermediate/calibration_literature/dias2025_knn.txt), [rendered p.4](../intermediate/calibration_literature/dias2025_p4.png), [rendered p.5](../intermediate/calibration_literature/dias2025_p5.png).

## Exactly what can be converted

For a displacement Δx in native pixels over k frame intervals:

    v_image [px/s] = Δx_px × FPS / k
    v_image [m/s]  = Δx_px × FPS / k × s_x

The first expression uses known pixel and time coordinates. The second needs the **axial object-space scale s_x in m/pixel**. FPS alone supplies time, not length.

Current [v14_geometry.py](../../v14_new_voidfraction/scripts/v14_geometry.py) records the owner-confirmed full lumen as the half-open slice **18:71: 53 pixel cells**. This supersedes older v14 reports using 60 pixels. The geometric candidate is:

    s_geometric = 0.001 m / 53 px
                = 1.8867924528301888 × 10⁻⁵ m/px
                = 18.867924528301888 μm/px

Using this number axially assumes that the literature's 1 mm diameter applies to these files, that the visible radial span represents that physical diameter, and that its apparent scale equals the axial scale at the tracked feature locations. Square sensor pixels do not prove the last assumption: cylindrical glass/refraction, viewing depth and lens mapping can affect radial and axial directions differently. The confirmed pixel boundaries resolve a geometric indexing ambiguity; they do not calibrate that optical mapping.

An appendix may multiply accepted px/s values by this factor and label every axis and caption **“m/s under the assumed geometric scale D=1 mm, 53 px”**. It would be an assumption-dependent rescaling of the same optical observable. Its apparently realistic magnitude would not validate the assumed scale or phase assignment. No experimentally justified scale uncertainty was found, so no confidence interval on m/s follows from this conversion. The within-run ratio cancels a common scale factor but retains tracking, segmentation and weighting errors.

## Why the requested chart is not automatically a superficial-velocity map

The example's axes J_g and J_l are volumetric phase flow rates divided by the **total tube area**. For an instantaneous cross-section, J_g=α U_g and J_l=(1−α)U_l only when U_g and U_l are the corresponding **phase-area mean material velocities**. With fluctuations, mean superficial flow involves the mean of the product; in general mean(α U_g) is not mean(α)×mean(U_g).

A moving gas outline, gas/liquid passage-conditioned correlation, local optical patch, or liquid-film wave is a different observable. The present images do not independently establish that an optical structure moves with the phase-area mean liquid or gas material. Consequently, multiplying an image speed by a v14 mask fraction would not establish J_g or J_l. V14's projected/circular reconstruction also carries separate three-dimensional geometry assumptions.

For annular flow, the core and film coexist at the same time. Spatially separating central and near-wall signals can estimate different optical motions, including interface waves. Their identification as U_g and U_l would require additional physical assumptions or material-following texture. A wave's propagation speed cannot be set equal to mean film transport solely because it lies at the liquid boundary. No optical estimator or scale conversion repairs an absent second independent material observable.

The defensible chart can copy the **layout** of the example: one point per supported recording, pattern-coded symbols, x=one explicitly named optical motion and y=the other. It should not reuse J_g/J_l labels for these image quantities. A log chart may show magnitudes of same-direction signed means, with the sign convention documented; missing, near-zero or opposed-direction estimates must remain flagged instead of being made positive silently. Additional pattern-specific optical quantities belong in clearly separate panels or figures when their physical meanings differ.

## Adversarial conclusions and remaining gap

1. **Documented FPS is useful and sufficient for px/s**, provided native indices and comparison lag are used. Using workbook fs_hz=4000, Output_Rate=5 or a blanket literature frame rate would mis-time many recordings.
2. **An approximate m/s appendix is possible under declared geometry.** Neither the provided metadata nor the two PDFs independently validates the required axial scale for these recordings.
3. **A phase-material or superficial-velocity chart remains a stronger claim than an optical-motion chart.** Report observed support and unresolved phase identity; do not choose a method because its points resemble expected flow-pattern regions.
4. **The limiting missing evidence is identifiable material transport and a verified axial length mapping**, rather than a missing FPS value. This audit identifies that limitation without asking for additional unavailable experimental records.

**Recomendação em PT-BR:** usar **px/s** no gráfico principal, com nomes explícitos para os movimentos ópticos dos dois eixos. Apresentar **m/s somente como conversão geométrica condicional**, se útil, com D=1 mm e 53 pixels escritos na legenda. Não chamar esses eixos de velocidades superficiais J_g/J_l nem interpretar a razão como escorregamento material gás/líquido sem a identificação física correspondente.

Reproduce with python3 -B src2/v15.1_gas_liquid_velocity/scripts/calibration_audit.py. The script reads and fingerprints original inputs; it writes only inside v15.1. This is a timing/calibration-evidence audit, not experimental validation of velocity accuracy.
