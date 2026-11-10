# Mass-flow provenance and sensitivity audit

**Result:** a mass-flow signal with owner-confirmed g/s units exists in the original LabVIEW workbooks. It can support a condition-level alternative to nominal G. The available evidence does not prove that each source block is synchronized with the corresponding camera acquisition.

Reproduce with `python3 src2/v15_gas_liquid_velocity/scripts/massflow_audit.py`. Detailed outputs: [per condition](../results/massflow_per_run.csv), [per original source block](../results/massflow_source_blocks.csv), [hashes, mapping and validation](../results/validation/massflow_audit.json). Original data are not modified.

## 1. Channel identity and the two different data products

`src2/data_processed/CLASSIFICACAO_PADRAOES_MANUAL.md:2–18` records confirmation by the data owner: Voltage_0 is differential pressure in kPa, Voltage_1 is mass flow in g/s, and Voltage_2 is optical luminosity in V. The same statement is in `src/data_processed/CLASSIFICACAO_PADRAOES_MANUAL.md`. The generic Voltage_1 header does not override these documented engineering units. Calling the original mass-flow values raw voltage and demanding a new voltage-to-flow curve disregards the available metadata. The sensor's calibration uncertainty is not independently established by this audit.

The `src2` workbooks are derived from images. They reuse these column names for layout compatibility: `pipeline_v11.py:185–187`, `pipeline_v9.py:194–197`, and `pipeline_v13.py:49–51` write zeros in Voltage_0 and Voltage_1. The audit confirms both columns are entirely zero in all 26 selected `__timeseries_v9.xlsx` files. These zeros are placeholders, not zero physical flow or pressure.

## 2. Provenance, unit conversion and numerical checks

The audit selected 26 camera folders under the requested exclusions and mapped each `mXXqYY_001` to the same **condition name** `src/data_processed/mXXqYY/dados.xlsx`. All 241,000 processed rows contain finite, positive mass-flow values. Their source_file column identifies 241 original workbooks under `src/data/PROGRAMA_LABVIEW/mXXqYY/`. Every referenced source exists, belongs to the expected condition, and all three numeric channels match the processed block exactly. No repeated numeric blocks were found within an individual condition. Each original source block contains 1000 samples.

Original LabVIEW XLSX style definitions are incompatible with the installed openpyxl reader. The audit reads their actual worksheet XML directly, including inline/shared strings, without altering the files. Original and processed file hashes, marker hashes and camera config hashes are recorded. This verifies numeric preservation; it does not by itself prove sensor calibration or camera synchronization.

The conversion uses documented mass-flow units and the assumed 1 mm internal diameter:

```text
A = pi (0.001 m)^2 / 4
G [kg m^-2 s^-1] = mean(Voltage_1 [g/s]) * 1e-3 / A
```

The primary summary is the arithmetic mean over all available samples for that condition. Median and a 10% trimmed mean are robustness checks, not automatically better physical estimators. No samples are clipped, and the audit reports finite/nonfinite and nonpositive counts. Empirical minimum, maximum, 5th/95th percentiles, standard deviation and source-block mean ranges appear in the CSV. These describe stored signal variation; they are not confidence intervals for true mass flow. The diameter assumption enters quadratically and its uncertainty is not supplied here.

## 3. Results for the selected conditions

Delta is 100 × (G_mean/G_nominal − 1). G values are kg m⁻² s⁻¹; mass flow is g/s. All displayed G values are conditional on documented units and D=1 mm.

| Camera condition mapping | n samples | Mean mass flow | Nominal G | Signal mean G | Signal median G | Delta (%) |
|---|---:|---:|---:|---:|---:|---:|
| m03q05_001 | 13000 | 0.334703 | 380 | 426.158 | 426.400 | +12.147 |
| m03q10_001 | 13000 | 0.327583 | 380 | 417.092 | 417.219 | +9.761 |
| m03q20_001 | 15000 | 0.325672 | 380 | 414.659 | 414.714 | +9.121 |
| m03q40_001 | 15000 | 0.330966 | 380 | 421.399 | 421.670 | +10.894 |
| m04q01_001 | 8000 | 0.416589 | 510 | 530.418 | 530.735 | +4.004 |
| m04q10_001 | 17000 | 0.412087 | 510 | 524.686 | 524.892 | +2.880 |
| m04q20_001 | 7000 | 0.411588 | 510 | 524.050 | 523.779 | +2.755 |
| m04q40_001 | 8000 | 0.417661 | 510 | 531.782 | 532.126 | +4.271 |
| m04q80_001 | 12000 | 0.415982 | 510 | 529.645 | 529.901 | +3.852 |
| m05q01_001 | 2000 | 0.477469 | 640 | 607.932 | 607.804 | -5.011 |
| m05q05_001 | 5000 | 0.548824 | 640 | 698.784 | 699.062 | +9.185 |
| m05q10_001 | 9000 | 0.516975 | 640 | 658.233 | 658.442 | +2.849 |
| m05q20_001 | 12000 | 0.510079 | 640 | 649.453 | 649.538 | +1.477 |
| m05q40_001 | 9000 | 0.532375 | 640 | 677.840 | 678.195 | +5.913 |
| m05q80_001 | 10000 | 0.527884 | 640 | 672.122 | 672.353 | +5.019 |
| m06q01_001 | 9000 | 0.677370 | 760 | 862.454 | 863.217 | +13.481 |
| m06q05_001 | 9000 | 0.691019 | 760 | 879.833 | 879.910 | +15.767 |
| m06q10_001 | 6000 | 0.682827 | 760 | 869.402 | 869.060 | +14.395 |
| m06q20_001 | 7000 | 0.681058 | 760 | 867.150 | 867.389 | +14.099 |
| m06q40_001 | 6000 | 0.674416 | 760 | 858.694 | 858.764 | +12.986 |
| m06q80_001 | 4000 | 0.661334 | 760 | 842.037 | 842.350 | +10.794 |
| m07q01_001 | 12000 | 0.802488 | 890 | 1021.759 | 1022.084 | +14.804 |
| m07q05_001 | 8000 | 0.816767 | 890 | 1039.940 | 1040.169 | +16.847 |
| m07q20_001 | 9000 | 0.808263 | 890 | 1029.112 | 1029.040 | +15.631 |
| m07q40_001 | 9000 | 0.780746 | 890 | 994.077 | 994.261 | +11.694 |
| m07q80_001 | 7000 | 0.765219 | 890 | 974.307 | 974.786 | +9.473 |

Mean deviations from nominal span **-5.011% to +16.847%**. The largest absolute difference between mean and median is **0.0884%**; the largest 10% trimmed-mean change relative to the ordinary mean is **0.0411%**. The nominal discrepancy is therefore not removed by choosing a robust central statistic. The old 4–13% agreement claim does not describe these 26 selected conditions and should not be reinstated.

## 4. Mapping and timing limitations

The mapping is to a shared `mXXqYY` operating-condition name. Original sensor paths do not contain the camera suffix `_001`, camera frame identifiers or a timestamp column. Each processed condition combines multiple `dados*.xlsx` source blocks. All inspected original workbook core created/modified date strings use `2006-09-16T00:00:00Z`, which does not establish acquisition timing. The camera configuration has no demonstrated correspondence between sensor row numbers and camera frames. Spreadsheet metadata and numerical similarity are not enough to assert simultaneous measurement.

Consequently, do not assign the 1000-sample block sequence to the 501 camera frames, infer a common duration from file counts, or treat the condition mean as an exact synchronized measurement for a particular image. An acquisition log or owner-confirmed linking rule is needed to strengthen that association. No sensor sampling frequency is inferred here.

The owner-confirmed mass-flow units resolve the physical meaning of the source values. What remains uncertain is the exact acquisition pairing, applicable sensor calibration uncertainty and the diameter at the station used to define G. Differential pressure is available, but does not by itself supply absolute pressure or inlet enthalpy.

## 5. Consequences for velocity calculations

The signal-derived condition mean G is a defensible **sensitivity alternative** and, once acquisition pairing is established, preferable to treating a nominal setpoint as measured flow. This does not identify inlet quality, net heating, local densities, void fraction or phase velocity. All 26 mass-flow records can be evaluated alongside nominal G without fitting to the expected pattern ordering.

Under the legacy zero-inlet-quality, constant-property balance, x is inversely proportional to G while `Jg = 4 q'' L/(D h_lv rho_g)` remains exactly independent of G. Therefore replacing G alone changes x, Jl, alpha and phase velocities but cannot correct a discrepancy in Jg. The exact superficial-liquid change is `(G_signal − G_nominal)/rho_l`. If inlet quality is nonzero, the additional term `G x_in/rho_g` restores G dependence in Jg.

The mean-versus-median agreement supports using the ordinary condition mean for this sensitivity analysis. It does not select RA, Xu–Fang or another void-fraction closure as most accurate, and it does not validate absolute gas or liquid velocities.
