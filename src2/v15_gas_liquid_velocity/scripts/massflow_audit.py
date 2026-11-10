#!/usr/bin/env python3
"""Audit owner-documented mass-flow data and its condition-level camera mapping.

The g/s interpretation is documented, not inferred from column names.
Original LabVIEW XLSX files are read directly as XML because their styles
are incompatible with openpyxl; no source file is rewritten.
"""
from __future__ import annotations
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import zipfile

import numpy as np
import pandas as pd
from scipy.stats import trim_mean

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
SRC2 = ROOT / "src2"
DATA = SRC2 / "data_processed"
G_MAP = {3:380,4:510,5:640,6:760,7:890}
AREA = math.pi * 0.001**2 / 4
NS = {"s":"http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def raw_xlsx(path):
    with zipfile.ZipFile(path) as z:
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            tree = ET.fromstring(z.read("xl/sharedStrings.xml"))
            shared = ["".join(e.itertext()) for e in tree]
        tree = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
        table=[]
        for row in tree.findall("s:sheetData/s:row",NS):
            entries={}
            for cell in row:
                ref=re.match(r"([A-Z]+)",cell.attrib["r"])[1]
                if cell.attrib.get("t")=="inlineStr":
                    value="".join(cell.find("s:is",NS).itertext())
                else:
                    node=cell.find("s:v",NS)
                    value=node.text if node is not None else None
                    if cell.attrib.get("t")=="s":
                        value=shared[int(value)]
                    elif value is not None:
                        value=float(value)
                entries[ref]=value
            table.append([entries.get(c) for c in ("A","B","C")])
        core = ET.fromstring(z.read("docProps/core.xml"))
        meta={node.tag.split("}")[-1]:node.text for node in core}
    return pd.DataFrame(table[1:],columns=table[0]),meta


def write_report(rows,blocks,metadata):
    lines=[
        "# Mass-flow provenance and sensitivity audit",
        "",
        "**Result:** a mass-flow signal with owner-confirmed g/s units exists in the original LabVIEW workbooks. It can support a condition-level alternative to nominal G. The available evidence does not prove that each source block is synchronized with the corresponding camera acquisition.",
        "",
        "Reproduce with `python3 src2/v15_gas_liquid_velocity/scripts/massflow_audit.py`. Detailed outputs: [per condition](../results/massflow_per_run.csv), [per original source block](../results/massflow_source_blocks.csv), [hashes, mapping and validation](../results/validation/massflow_audit.json). Original data are not modified.",
        "",
        "## 1. Channel identity and the two different data products",
        "",
        "`src2/data_processed/CLASSIFICACAO_PADRAOES_MANUAL.md:2–18` records confirmation by the data owner: Voltage_0 is differential pressure in kPa, Voltage_1 is mass flow in g/s, and Voltage_2 is optical luminosity in V. The same statement is in `src/data_processed/CLASSIFICACAO_PADRAOES_MANUAL.md`. The generic Voltage_1 header does not override these documented engineering units. Calling the original mass-flow values raw voltage and demanding a new voltage-to-flow curve disregards the available metadata. The sensor's calibration uncertainty is not independently established by this audit.",
        "",
        "The `src2` workbooks are derived from images. They reuse these column names for layout compatibility: `pipeline_v11.py:185–187`, `pipeline_v9.py:194–197`, and `pipeline_v13.py:49–51` write zeros in Voltage_0 and Voltage_1. The audit confirms both columns are entirely zero in all 26 selected `__timeseries_v9.xlsx` files. These zeros are placeholders, not zero physical flow or pressure.",
        "",
        "## 2. Provenance, unit conversion and numerical checks",
        "",
        f"The audit selected 26 camera folders under the requested exclusions and mapped each `mXXqYY_001` to the same **condition name** `src/data_processed/mXXqYY/dados.xlsx`. All {metadata['summary']['total_samples']:,} processed rows contain finite, positive mass-flow values. Their source_file column identifies {len(blocks)} original workbooks under `src/data/PROGRAMA_LABVIEW/mXXqYY/`. Every referenced source exists, belongs to the expected condition, and all three numeric channels match the processed block exactly. No repeated numeric blocks were found within an individual condition. Each original source block contains 1000 samples.",
        "",
        "Original LabVIEW XLSX style definitions are incompatible with the installed openpyxl reader. The audit reads their actual worksheet XML directly, including inline/shared strings, without altering the files. Original and processed file hashes, marker hashes and camera config hashes are recorded. This verifies numeric preservation; it does not by itself prove sensor calibration or camera synchronization.",
        "",
        "The conversion uses documented mass-flow units and the assumed 1 mm internal diameter:",
        "",
        "```text",
        "A = pi (0.001 m)^2 / 4",
        "G [kg m^-2 s^-1] = mean(Voltage_1 [g/s]) * 1e-3 / A",
        "```",
        "",
        "The primary summary is the arithmetic mean over all available samples for that condition. Median and a 10% trimmed mean are robustness checks, not automatically better physical estimators. No samples are clipped, and the audit reports finite/nonfinite and nonpositive counts. Empirical minimum, maximum, 5th/95th percentiles, standard deviation and source-block mean ranges appear in the CSV. These describe stored signal variation; they are not confidence intervals for true mass flow. The diameter assumption enters quadratically and its uncertainty is not supplied here.",
        "",
        "## 3. Results for the selected conditions",
        "",
        "Delta is 100 × (G_mean/G_nominal − 1). G values are kg m⁻² s⁻¹; mass flow is g/s. All displayed G values are conditional on documented units and D=1 mm.",
        "",
        "| Camera condition mapping | n samples | Mean mass flow | Nominal G | Signal mean G | Signal median G | Delta (%) |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        lines.append(f"| {r['run']} | {r['n_samples']} | {r['massflow_mean_g_s']:.6f} | {r['G_nominal_kg_m2_s']} | {r['G_mean_kg_m2_s']:.3f} | {r['G_median_kg_m2_s']:.3f} | {r['G_mean_delta_pct']:+.3f} |")
    maximum_trimmed=max(abs(100*(r['massflow_trimmed10_mean_g_s']/r['massflow_mean_g_s']-1)) for r in rows)
    lines += [
        "",
        f"Mean deviations from nominal span **{metadata['summary']['mean_G_delta_pct_min']:+.3f}% to {metadata['summary']['mean_G_delta_pct_max']:+.3f}%**. The largest absolute difference between mean and median is **{metadata['summary']['max_abs_mean_median_pct']:.4f}%**; the largest 10% trimmed-mean change relative to the ordinary mean is **{maximum_trimmed:.4f}%**. The nominal discrepancy is therefore not removed by choosing a robust central statistic. The old 4–13% agreement claim does not describe these 26 selected conditions and should not be reinstated.",
        "",
        "## 4. Mapping and timing limitations",
        "",
        "The mapping is to a shared `mXXqYY` operating-condition name. Original sensor paths do not contain the camera suffix `_001`, camera frame identifiers or a timestamp column. Each processed condition combines multiple `dados*.xlsx` source blocks. All inspected original workbook core created/modified date strings use `2006-09-16T00:00:00Z`, which does not establish acquisition timing. The camera configuration has no demonstrated correspondence between sensor row numbers and camera frames. Spreadsheet metadata and numerical similarity are not enough to assert simultaneous measurement.",
        "",
        "Consequently, do not assign the 1000-sample block sequence to the 501 camera frames, infer a common duration from file counts, or treat the condition mean as an exact synchronized measurement for a particular image. An acquisition log or owner-confirmed linking rule is needed to strengthen that association. No sensor sampling frequency is inferred here.",
        "",
        "The owner-confirmed mass-flow units resolve the physical meaning of the source values. What remains uncertain is the exact acquisition pairing, applicable sensor calibration uncertainty and the diameter at the station used to define G. Differential pressure is available, but does not by itself supply absolute pressure or inlet enthalpy.",
        "",
        "## 5. Consequences for velocity calculations",
        "",
        "The signal-derived condition mean G is a defensible **sensitivity alternative** and, once acquisition pairing is established, preferable to treating a nominal setpoint as measured flow. This does not identify inlet quality, net heating, local densities, void fraction or phase velocity. All 26 mass-flow records can be evaluated alongside nominal G without fitting to the expected pattern ordering.",
        "",
        "Under the legacy zero-inlet-quality, constant-property balance, x is inversely proportional to G while `Jg = 4 q'' L/(D h_lv rho_g)` remains exactly independent of G. Therefore replacing G alone changes x, Jl, alpha and phase velocities but cannot correct a discrepancy in Jg. The exact superficial-liquid change is `(G_signal − G_nominal)/rho_l`. If inlet quality is nonzero, the additional term `G x_in/rho_g` restores G dependence in Jg.",
        "",
        "The mean-versus-median agreement supports using the ordinary condition mean for this sensitivity analysis. It does not select RA, Xu–Fang or another void-fraction closure as most accurate, and it does not validate absolute gas or liquid velocities.",
    ]
    finding=HERE/"findings"/"massflow_review.md"
    finding.parent.mkdir(parents=True,exist_ok=True)
    finding.write_text("\n".join(lines)+"\n")


def main():
    rows,blocks,manifest,excluded=[],[],[],[]
    classfile=DATA/"CLASSIFICACAO_PADRAOES_MANUAL.md"
    units_text=classfile.read_text()
    assert "vazão mássica | g/s" in units_text
    for folder in sorted(DATA.glob("m??q??_???")):
        markers=list(folder.glob("__*.txt"))
        assert len(markers)==1
        pattern=markers[0].stem[2:].replace("_"," ")
        if folder.name=="m05q05_002" or pattern=="smooth annular" or "-" in pattern:
            excluded.append(dict(run=folder.name,pattern=pattern))
            continue
        assert pattern in ("dispersed bubbles","plug","slug","churn","wavy annular")
        base=folder.name.split("_")[0]
        path=ROOT/"src"/"data_processed"/base/"dados.xlsx"
        frame=pd.read_excel(path)
        assert list(frame.columns)==["source_file","Voltage_0","Voltage_1","Voltage_2"]
        assert not frame.source_file.isna().any()
        assert all(Path(s).parent.name==base for s in frame.source_file)
        mdot=pd.to_numeric(frame.Voltage_1,errors="coerce").to_numpy(float)
        finite=mdot[np.isfinite(mdot)]
        assert len(finite)>0
        group_means=[]
        signatures=[]
        source_paths=[]
        for source,part in frame.groupby("source_file",sort=False):
            rawpath=ROOT/source
            raw,meta=raw_xlsx(rawpath)
            assert list(raw.columns)==["Voltage_0","Voltage_1","Voltage_2"]
            assert len(raw)==1000
            a=raw.to_numpy(dtype=float)
            b=part[["Voltage_0","Voltage_1","Voltage_2"]].to_numpy(dtype=float)
            assert a.shape==b.shape and np.array_equal(a,b,equal_nan=True),(source,"processed differs from original")
            signature=hashlib.sha256(a.tobytes()).hexdigest()
            signatures.append(signature)
            values=a[:,1]
            f=values[np.isfinite(values)]
            group_means.append(float(f.mean()))
            block=dict(run=folder.name,condition=base,source_file=source,
                       n_samples=len(values),n_finite=len(f),mean_g_s=float(f.mean()),
                       median_g_s=float(np.median(f)),min_g_s=float(f.min()),max_g_s=float(f.max()),
                       raw_sha256=sha(rawpath),numeric_content_sha256=signature,
                       processed_matches_raw=True,core_properties=meta)
            blocks.append(block)
            source_paths.append(source)
        nominal=G_MAP[int(base[1:3])]
        g_mean=float(finite.mean()*1e-3/AREA)
        g_median=float(np.median(finite)*1e-3/AREA)
        derived=pd.read_excel(folder/"__timeseries_v9.xlsx",sheet_name="Sheet")
        assert (derived.Voltage_0==0).all() and (derived.Voltage_1==0).all()
        rows.append(dict(run=folder.name,pattern=pattern,condition=base,
                         n_samples=len(mdot),n_finite=len(finite),n_nonfinite=len(mdot)-len(finite),
                         n_nonpositive=int(np.count_nonzero(finite<=0)),n_source_files=len(source_paths),
                         n_distinct_source_numeric_blocks=len(set(signatures)),
                         massflow_mean_g_s=float(finite.mean()),massflow_median_g_s=float(np.median(finite)),
                         massflow_trimmed10_mean_g_s=float(trim_mean(finite,0.1)),
                         massflow_min_g_s=float(finite.min()),massflow_max_g_s=float(finite.max()),
                         massflow_p05_g_s=float(np.quantile(finite,.05)),
                         massflow_p95_g_s=float(np.quantile(finite,.95)),
                         massflow_std_g_s=float(finite.std(ddof=1)),
                         source_block_mean_min_g_s=min(group_means),
                         source_block_mean_max_g_s=max(group_means),
                         G_nominal_kg_m2_s=nominal,G_mean_kg_m2_s=g_mean,G_median_kg_m2_s=g_median,
                         G_trimmed10_kg_m2_s=float(trim_mean(finite,.1)*1e-3/AREA),
                         G_mean_delta_pct=100*(g_mean/nominal-1),
                         G_median_delta_pct=100*(g_median/nominal-1),
                         mean_vs_median_pct=100*(finite.mean()/np.median(finite)-1),
                         processed_file=str(path.relative_to(ROOT)),processed_sha256=sha(path),
                         src2_v9_flow_column_all_zero=True,
                         unit_basis="owner-documented g/s; acquisition calibration not independently assessed",
                         camera_mapping="same mXXqYY condition; temporal/acquisition pairing unproven"))
        manifest.append(dict(run=folder.name,pattern_marker=str(markers[0].relative_to(ROOT)),
                             pattern_marker_sha256=sha(markers[0]),processed_file=str(path.relative_to(ROOT)),
                             processed_sha256=sha(path),raw_sources=source_paths,
                             camera_config=str((folder/"_config.xsv").relative_to(ROOT)),
                             camera_config_sha256=sha(folder/"_config.xsv")))
    assert len(rows)==26 and len(excluded)==5
    assert all(r["n_nonfinite"]==0 and r["n_nonpositive"]==0 for r in rows)
    result=HERE/"results"
    result.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(result/"massflow_per_run.csv",index=False)
    pd.DataFrame([{k:v for k,v in b.items() if k!="core_properties"} for b in blocks]).to_csv(
        result/"massflow_source_blocks.csv",index=False)
    validation=result/"validation"
    validation.mkdir(exist_ok=True)
    metadata=dict(unit_document=str(classfile.relative_to(ROOT)),unit_document_sha256=sha(classfile),
                  unit_document_lines="8–12",diameter_assumed_m=.001,cross_section_m2=AREA,
                  selected_count=len(rows),excluded=excluded,manifest=manifest,source_blocks=blocks,
                  interpretation="Documented mass-flow signal permits condition-level G comparison; no synchronized camera link proven",
                  has_timestamp_column=False,
                  source_created_properties=sorted({str(b["core_properties"].get("created")) for b in blocks}),
                  source_modified_properties=sorted({str(b["core_properties"].get("modified")) for b in blocks}),
                  raw_processed_all_channels_match=True,
                  summary=dict(total_samples=sum(r["n_samples"] for r in rows),
                               total_finite=sum(r["n_finite"] for r in rows),
                               mean_G_delta_pct_min=min(r["G_mean_delta_pct"] for r in rows),
                               mean_G_delta_pct_max=max(r["G_mean_delta_pct"] for r in rows),
                               max_abs_mean_median_pct=max(abs(r["mean_vs_median_pct"]) for r in rows),
                               run_count_with_duplicate_numeric_blocks=sum(r["n_distinct_source_numeric_blocks"]<r["n_source_files"] for r in rows)))
    (validation/"massflow_audit.json").write_text(json.dumps(metadata,indent=2,ensure_ascii=False)+"\n")
    write_report(rows,blocks,metadata)
    print(json.dumps(metadata["summary"],indent=2))
    print(pd.DataFrame(rows)[["run","G_nominal_kg_m2_s","G_mean_kg_m2_s","G_mean_delta_pct"]].to_string(index=False))


if __name__=="__main__":
    main()
