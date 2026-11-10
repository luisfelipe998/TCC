# Independent audit of v12 and v12.1

Audit date: 2026-10-23. The original versions were read and imported without modifying them. Reproduction: `python3 src2/v15_gas_liquid_velocity/scripts/audit_legacy.py`. Machine-readable evidence, original-file SHA-256 hashes, per-run values and workbook headers: [legacy_audit.json](../results/validation/legacy_audit.json).

## 1. Executive finding

The numerical calculations and stored tables reproduce correctly **for their stated nominal inputs and zero inlet-quality assumption**. The principal defects are interpretation, dataset selection and unsupported validation claims. Neither version measures phase velocity from `src2/data_processed`. A correct equation evaluated using an unknown inlet state is not an accurate measurement.

- v12 computes superficial velocities `Jg = G x / rho_g` and `Jl = G (1-x) / rho_l`. These are phase volumetric flow rates divided by the **whole tube area**. They do not require equal phase velocities or a void-fraction correlation.
- v12.1 retains exactly the same `x`, `Jg` and `Jl`, and adds Rouhani–Axelsson (RA) void fraction, then conditional phase-mean velocities `ug = Jg/alpha`, `ul = Jl/(1-alpha)` and slip `S = ug/ul`. These are **model estimates of phase velocities**, not independently measured velocities.
- Current scope contains **26 acquisitions**: dispersed bubbles 1, plug 3, slug 3, churn 10 and wavy annular 9. The legacy 28-row summaries are outside this scope.
- There is no measured phase-velocity reference in these scripts with which to calculate physical velocity error, select a universally most accurate closure, or give a defensible accuracy percentage. Literature-range overlap does not provide that reference. **However, a documented mass-flow signal does exist in `src/data_processed`: the legacy claim that Voltage_1 is uncalibrated raw voltage overlooks owner-confirmed units.** See [mass-flow audit](massflow_review.md).

## 2. What was independently recomputed

The audit obtains R600a properties with CoolProp 8.0.0 at 300.15 K, reconstructs total heater power `Q = q'' pi D L`, obtains vapor mass rate `Q/h_lv`, and divides phase volume rates by tube area. The RA calculation is checked in a different algebraic form, `ug = C0 (Jg+Jl) + Vgj`, followed by `alpha = Jg/ug`. Neither legacy computation function is used to generate these independent values. The library properties and physical assumptions remain shared inputs; this is not an independent thermophysical-property experiment.

| Check | Observed result | Meaning |
|---|---:|---|
| Original retained rows | 28 in each version | Reproduces the historical scope |
| v12 versus v12.1, `x/Jg/Jl` | Maximum difference exactly 0 | RA does not change superficial velocities |
| Independent versus v12.1, all seven quantities | Maximum absolute difference 1.066e-14 | Floating-point agreement |
| Original per-run displayed cells | 84/84 v12; 196/196 v12.1 agree | Stored table reproducibility |
| Original summary displayed cells | 40/40 v12; 50/50 v12.1 agree | Stored arithmetic summaries reproduce |
| `rho_g Jg + rho_l Jl = G`, current scope | Maximum residual 1.137e-13 kg m⁻² s⁻¹ | Mass conservation |
| `alpha ug = Jg` | Maximum residual 2.221e-16 m/s | Correct gas flux definition |
| `(1-alpha) ul = Jl` | Maximum residual 1.111e-16 m/s | Correct liquid flux definition |
| Source and table hashes before/after import | Unchanged | Original artifacts preserved |

Both legacy sets of internal checks also pass. Their pass status verifies the implemented model, not the experimental truth of the inputs or outputs. No unit-conversion or arithmetic averaging error was found in the tested legacy tables.

## 3. Correct dataset selection and its effect

The current audit discovers actual folders and their unique `__*.txt` pattern markers. It excludes the following acquisitions before recomputing summaries:

| Acquisition | Exclusion |
|---|---|
| `m03q01_001` | plug-slug |
| `m03q80_001` | smooth annular |
| `m04q05_001` | slug-churn |
| `m05q05_002` | Explicit user exclusion |
| `m07q10_001` | slug-churn |

v12 merges every label containing `annular` into wavy annular (`velocity_map_v12.py:125–133`) and explicitly adds the second `m05q05` record (`183–210`). v12.1 repeats these choices (`velocity_map_v12_1.py:132–139,230–266`). Both use hardcoded labels and parameter maps, rather than reading acquisition data. The old explanation that `m05q05_002` has a plug marker is stale relative to the present `src2` marker, which says dispersed bubbles; this acquisition is excluded regardless. No relabeling is needed.

The following are **conditional arithmetic means over the 26 selected operating points**, retaining the legacy thermodynamic assumptions. They are not measured class velocities or a recommendation to use RA as ground truth. Units are m/s, except alpha and S.

| Pattern | n | Mean Jg | Mean Jl | Mean alpha RA | Mean ug RA | Mean ul RA | Mean S RA |
|---|---:|---:|---:|---:|---:|---:|---:|
| Dispersed bubbles | 1 | 0.336 | 1.618 | 0.145 | 2.322 | 1.892 | 1.228 |
| Plug | 3 | 0.785 | 1.379 | 0.264 | 2.554 | 1.980 | 1.280 |
| Slug | 3 | 1.233 | 1.140 | 0.402 | 2.784 | 2.036 | 1.358 |
| Churn | 10 | 6.556 | 1.092 | 0.730 | 8.582 | 4.411 | 1.866 |
| Wavy annular | 9 | 17.558 | 0.747 | 0.863 | 19.862 | 5.996 | 3.174 |

Removing the duplicate nominal slug operating point changes its mean Jg from 1.345 to 1.233 m/s and mean ug from 2.908 to 2.784 m/s. Removing smooth annular changes the wavy-annular mean Jg from 18.492 to 17.558 m/s and its median Jg from 20.174 to 13.449 m/s. Its maximum modeled quality falls from 0.6831 to 0.5089, maximum alpha from 0.9541 to 0.9265 and maximum slip from 5.889 to 4.673. These are scope effects, not improvements in physical accuracy.

Each summary averages **per-run derived values**. Mean slip is `mean(ug/ul)`, not `mean(ug)/mean(ul)`. Different G/q combinations and unequal class sizes confound between-pattern means; one dispersed-bubble acquisition cannot characterize its population variance.

## 4. Findings that materially change interpretation

### 4.1 Superficial velocity is not actual phase speed or slip

`v12/INSIGHTS.md:128–133` describes gas traveling at about 27 m/s while liquid travels below 1 m/s, and interprets decreasing Jl as liquid deceleration. Those are superficial velocities. Neither conclusion about actual phase speed follows without alpha. Its own later v12.1 values give a counterexample: for the two retained annular points `m06q40_001` and `m06q80_001`, Jl decreases **1.1498 → 0.9130 m/s**, while RA phase-liquid velocity increases **6.8462 → 8.5885 m/s**.

The statement that v12 assumes gas and liquid travel at the same mean speed (`INSIGHTS.md:186–190`, also `velocity_map_v12.py:44–46,506–508`) is incorrect. Its formulas need no slip assumption. `Jg/Jl` is not the slip ratio: `S = (Jg/Jl)(1-alpha)/alpha`. The v12.1 description correctly explains invariance to slip, but the older contradictory narrative remains.

### 4.2 The inlet enthalpy assumption controls absolute velocities

Both versions hardcode `x = 4 q'' L/(G D h_lv)` (`v12:204–207`; `v12.1:248–251`). The more general steady energy balance is:

```text
x_out,eq = x_in,eq + (Q_to_fluid)/(mdot h_lv)
x_in,eq  = (h_in - h_l)/h_lv
Q_to_fluid = q'' pi D L  only if the specified heat flux is net heat into the fluid
```

It needs an inlet enthalpy and heat-loss convention, and consistently evaluated local thermodynamic properties. Converting a heated-section outlet estimate into the state at the camera station additionally assumes negligible intervening state change. Equilibrium quality may be negative for subcooled liquid; clipping it is not a measurement of vapor mass fraction in nonequilibrium subcooled boiling.

Calling `x_in=0` an upper bound is only justified after establishing `x_in<=0` and appropriate net-heat assumptions for the records. It is not an unconditional bound. Claiming that zero inlet quality **underestimates** slug void fraction (`v12.1:593,608–612`) while calling it an upper bound is unsupported. Literature and data mismatch alone cannot determine the sign of the missing inlet correction.

For fixed properties, zero inlet quality makes `Jg = 4 q'' L/(D h_lv rho_g)` independent of G. Thus vertical lines in the old velocity map are imposed algebraically, not discovered experimentally. With nonzero inlet quality, `G x_in/rho_g` contributes to Jg.

### 4.3 “Real phase velocities” should mean modeled phase means

The v12.1 conversion (`254–261`) is correct for phase-mean velocities at the specified state and void fraction. Its repeated label “velocidades reais” can be mistaken for direct measurement. Use **RA-estimated phase velocities, conditional on nominal inputs and x_in=0**. The audit cannot infer their experimental error.

Phase averages also require a consistent averaging convention: for fluctuating fields, `mean(Jg/alpha)` generally differs from `mean(Jg)/mean(alpha)`. A spatial panorama occupancy, a correlation evaluated at mean quality, and an instantaneous cross-sectional void fraction are not interchangeable without stated assumptions.

### 4.4 Internal assertions and class overlap are not calibration

v12.1 checks that annular maximum alpha exceeds 0.90 and median exceeds 0.85 (`315–324`). These are regime-based plausibility targets. They are not independent velocity observations. All selected churn RA alpha values fall within the quoted 0.61–0.85 band; only **6/9** retained annular values fall inside 0.86–0.99 and **2/3** plug values inside 0.042–0.29. None of these overlap counts establishes measurement accuracy. The generated statement that annular values are “within the range” (`594–597`) is too broad.

Fitting a closure or inlet quality to place each labeled pattern in the expected region would make class separation part of the fitting target. Evaluating on that same ordering would be circular. Compare methods using matched operating conditions, independently measured quality/void fraction/velocity where available, conservation and sensitivity—not the prettiest class clusters.

### 4.5 Several literature assertions remain incorrect

`v12/VALIDACAO_LITERATURA_v12.md:58–71,91–96,106,135` attributes quality and superficial-velocity differences to slip. Slip closure changes alpha and phase velocities at fixed mass quality; it does not repair or alter that quality's energy balance or Jg/Jl. v12.1 corrects this category error in its own introduction, leaving the v12 memo internally outdated.

`v12.1/velocity_map_v12_1.py:44` attributes RA to de Oliveira (2021), whereas that paper's extracted text explicitly names **Xu and Fang (2014)** (`v14_new_voidfraction/literature/extracts/deoliveira2021.txt:312–313`). RA is the cited Dias (2025) alternative. The v12.1 claim that the 2021 slug section has no alpha information (`599–601`) is contradicted by its text reporting **alpha approximately 0.65** (`deoliveira2021.txt:612`).

The newer [v14 source audit](../../v14_new_voidfraction/findings/source_audit.md) also documents contradictions among published x, alpha and S and confirms rendered source pages. Those reported ranges cannot establish the accuracy of our phase velocities. The v15 literature review should be used for source-specific conclusions; this legacy audit does not invent corrected published decimal values.

### 4.6 Legacy data provenance and the optical objection need correction

The v12.1 limitation cites `src/data_processed/<run>/dados.xlsx` (`616–619`), while the present task is about `src2/data_processed`. The audit inventories **104 workbooks across the 26 retained folders**. Their data sheets use position or time plus Voltage_0/1/2, with metadata sheets. **In these image-derived workbooks, Voltage_0/1 are placeholder zeros**, explicitly written by `src2/pipeline_v11.py:185–187`, `pipeline_v9.py:194–197` and `pipeline_v13.py:49–51`. The scripts themselves do not read these workbooks. Nominal G is a lookup (`v12:99–100`, `v12.1:106–107`), not a measured flow-rate result.

Crucially, the original LabVIEW signals under `src/data_processed` are different. `src2/data_processed/CLASSIFICACAO_PADRAOES_MANUAL.md:8–12` documents the owner's confirmation that **Voltage_1 is mass flow in g/s** and **Voltage_0 is differential pressure in kPa**. The old claims that these are only uncalibrated voltages and that measured mass flow does not exist are wrong in light of this metadata. The [new mass-flow audit](massflow_review.md) reproduces every original source block and computes condition-level G. Its mean deviations from nominal G span **−5.01% to +16.85%** for the selected conditions. Therefore the old unsupported **4–13%** numerical claim should not be reinstated; the newly computed values and their condition-level mapping limitation should replace it. No time synchronization with an individual camera acquisition has been established.

The old optical rejection treats `ScaleFactor=1` as a camera spatial-calibration field (`v12:18–22`; v12.1 adversarial memo). In the actual config it belongs to disabled **voltage acquisition channels**, with `Units=[V]` (`data_processed/m05q05_001/_config.xsv:308–327`). It proves neither 1 mm/pixel nor absence of every possible geometric calibration. Native-frame displacement in px/frame can be assessed for tracking consistency; conversion to m/s requires an independently justified optical length scale and native frame timing. A pipe-diameter calibration is conditional on reliable internal-wall location and optical distortion. Interface or wave displacement still is not automatically the bulk gas speed, and unseeded liquid regions do not directly provide liquid material velocity.

## 5. What carries forward

Retain the superficial-velocity definitions, the mass/energy conservation structure and the correctly coded RA equation. Replace hardcoded scope with the actual 26 selected folders, preserve full acquisition identifiers, and report thermodynamic assumptions alongside every numerical velocity table. Treat RA as one candidate closure and use homogeneous/Xu–Fang or other justified alternatives to quantify model sensitivity.

To establish actual accuracy, obtain the inlet thermodynamic state or enough measurements to reconstruct it, confirm the mapping of the available mass-flow records to the camera acquisitions, recover net heater power/heat loss and pressure/temperature at the comparison station, and obtain a defensible independent void-fraction or phase-velocity reference. Pending those inputs, a reproducible conditional estimate and a clearly separated optical structure-speed result are the strongest defensible outputs. No accuracy percentage follows from the present audit.
