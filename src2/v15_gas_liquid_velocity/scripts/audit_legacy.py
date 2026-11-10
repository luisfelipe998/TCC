#!/usr/bin/env python3
"""Recompute v12/v12.1 independently without rewriting their artifacts.

These are algebra/reproducibility checks, not experimental accuracy tests.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import statistics
import sys

HERE = Path(__file__).resolve().parents[1]
SRC2 = HERE.parent
OUT = HERE / "results" / "validation" / "legacy_audit.json"
os.environ.setdefault("MPLCONFIGDIR", str(HERE / "intermediate" / "legacy_audit_mplcache"))
sys.dont_write_bytecode = True

import CoolProp
from CoolProp.CoolProp import PropsSI
from openpyxl import load_workbook


def import_source(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_table(path):
    rows = {}
    for line in path.read_text().splitlines():
        if line.startswith("| m"):
            cells = [v.strip() for v in line.strip("|").split("|")]
            rows[cells[0]] = cells
    return rows


def independent(run, pattern, props):
    # Independently evaluate total heater power and phase mass rates.
    # Deliberately do not call either legacy calculation function here.
    m, q = int(run[1:3]), int(run[4:6])
    g = {3: 380, 4: 510, 5: 640, 6: 760, 7: 890}[m]
    diameter, length = 0.001, 0.265
    area = math.pi * diameter**2 / 4
    power = q * 1000 * math.pi * diameter * length
    mdot = g * area
    vapor_mdot = power / props["h_lv"]
    x = vapor_mdot / mdot
    jg = vapor_mdot / (props["rho_g"] * area)
    jl = (mdot - vapor_mdot) / (props["rho_l"] * area)
    c0 = 1 + 0.12 * (1 - x)
    drift_velocity = (1.18 * (1-x) *
                      (9.81 * props["sigma"] * (props["rho_l"] - props["rho_g"]))**0.25
                      / math.sqrt(props["rho_l"]))
    # Multiplying the published RA denominator by G yields alpha=jg/ug.
    ug = c0 * (jg + jl) + drift_velocity
    alpha = jg / ug
    ul = jl / (1-alpha)
    return dict(run=run, pattern=pattern, G=g, q_kW_m2=q, x=x, Jg=jg, Jl=jl,
                alpha=alpha, uv=ug, ul=ul, S=ug/ul, heating_power_W=power)


def summarize(rows):
    result = {}
    for pattern in ("bubble", "plug", "slug", "churn", "wavy annular"):
        subset = [r for r in rows if r["pattern"] == pattern]
        result[pattern] = {"n": len(subset)}
        for key in ("x", "Jg", "Jl", "alpha", "uv", "ul", "S"):
            vals = [r[key] for r in subset]
            result[pattern][key] = dict(min=min(vals), max=max(vals),
                                        mean=statistics.mean(vals), median=statistics.median(vals))
    return result


def main():
    paths = [SRC2 / "v12" / "velocity_map_v12.py",
             SRC2 / "v12.1" / "velocity_map_v12_1.py",
             SRC2 / "v12" / "velocity_table_v12.md",
             SRC2 / "v12.1" / "velocity_table_v12_1.md"]
    before = {str(p.relative_to(SRC2)): digest(p) for p in paths}
    v12, v121 = [import_source(n, p) for n, p in zip(("legacy_v12", "legacy_v121"), paths[:2])]
    props = dict(rho_l=PropsSI("D", "T", 300.15, "Q", 0, "R600a"),
                 rho_g=PropsSI("D", "T", 300.15, "Q", 1, "R600a"),
                 h_lv=PropsSI("H", "T", 300.15, "Q", 1, "R600a") -
                      PropsSI("H", "T", 300.15, "Q", 0, "R600a"),
                 sigma=PropsSI("I", "T", 300.15, "Q", 0, "R600a"))
    all12 = v12.compute_runs(v12.fluid_properties())
    all121 = v121.compute_runs(v121.fluid_properties())
    legacy_checks = {"v12": v12.sanity_checks(all12, v12.fluid_properties()),
                     "v12.1": v121.sanity_checks(all121, v121.fluid_properties())}
    rows12 = [r for r in all12 if r["pattern"]]
    rows121 = [r for r in all121 if r["pattern"]]
    independent_rows = [independent(r["run"], r["pattern"], props) for r in rows121]
    assert len(rows12) == len(rows121) == len(independent_rows) == 28
    max_error = {}
    for key in ("x", "Jg", "Jl", "alpha", "uv", "ul", "S"):
        err = max(abs(a[key]-b[key]) for a,b in zip(rows121, independent_rows))
        max_error[key] = err
        assert err < 1e-11, (key, err)
    cross_version = {key: max(abs(a[key]-b[key]) for a,b in zip(rows12, rows121))
                     for key in ("x", "Jg", "Jl")}
    assert all(v == 0 for v in cross_version.values())
    comparisons = {}
    for path, keys in ((paths[2], (("x",4),("Jg",3),("Jl",4))),
                       (paths[3], (("x",4),("Jg",3),("Jl",4),("alpha",4),
                                   ("uv",3),("ul",4),("S",3)))):
        table = read_table(path)
        assert set(table) == {r["run"] for r in independent_rows}
        checked = 0
        for row in independent_rows:
            for col, (key, digits) in enumerate(keys, start=4):
                assert table[row["run"]][col] == format(row[key], f".{digits}f"), (path, row, key)
                checked += 1
        comparisons[str(path.relative_to(SRC2))] = {"rows":len(table), "numeric_cells_checked":checked,
                                                  "all_match_displayed_precision": True}

    selected, excluded, workbook_inventory = [], [], []
    for folder in sorted((SRC2 / "data_processed").glob("m??q??_???")):
        markers = list(folder.glob("__*.txt"))
        assert len(markers) == 1, folder
        label = markers[0].stem[2:].replace("_"," ")
        reason = ("explicit acquisition exclusion" if folder.name == "m05q05_002" else
                  "smooth annular" if label == "smooth annular" else
                  "dual classification" if "-" in label else None)
        if reason:
            excluded.append(dict(run=folder.name, marker_label=label, reason=reason))
            continue
        pattern = "bubble" if label == "dispersed bubbles" else label
        assert pattern in ("bubble","plug","slug","churn","wavy annular"), pattern
        selected.append(independent(folder.name,pattern,props))
        for path in sorted(folder.glob("*.xlsx")):
            book = load_workbook(path,read_only=True,data_only=True)
            sheets = []
            for sheet in book.worksheets:
                sheets.append(dict(name=sheet.title,rows=sheet.max_row,columns=sheet.max_column,
                                   first_row=list(next(sheet.values,[]))))
            book.close()
            workbook_inventory.append(dict(file=str(path.relative_to(SRC2)),sheets=sheets))
    assert len(selected) == 26 and len(excluded) == 5
    summary_original, summary_selected = summarize(independent_rows), summarize(selected)
    pt = {"bolha dispersa":"bubble","plug":"plug","slug":"slug",
          "churn":"churn","anular ondulado":"wavy annular"}
    for path, columns in ((paths[2], [(k,s,d) for k,d in (("Jg",3),("Jl",4))
                                     for s in ("min","max","mean","median")]),
                          (paths[3], [(k,s,3) for k in ("alpha","S")
                                      for s in ("min","max","mean","median")] +
                                     [("uv","mean",3),("ul","mean",4)])):
        section = path.read_text().split("## Resumo por padrão",1)[1].split("\n## ",1)[0]
        checked = 0
        for line in section.splitlines():
            cells = [v.strip() for v in line.strip("|").split("|")]
            if cells[0] not in pt:
                continue
            expected = summary_original[pt[cells[0]]]
            assert int(cells[1]) == expected["n"]
            for cell,(key,stat,digits) in zip(cells[2:],columns):
                assert cell == format(expected[key][stat],f".{digits}f"), (path,cells[0],key,stat)
                checked += 1
        assert checked == 5*len(columns)
        comparisons[str(path.relative_to(SRC2))]["summary_numeric_cells_checked"] = checked
    residuals = {
        "mass_flux_kg_m2_s": max(abs(props["rho_g"]*r["Jg"] + props["rho_l"]*r["Jl"] - r["G"]) for r in selected),
        "gas_phase_flux_m_s": max(abs(r["alpha"]*r["uv"] - r["Jg"]) for r in selected),
        "liquid_phase_flux_m_s": max(abs((1-r["alpha"])*r["ul"] - r["Jl"]) for r in selected),
    }
    assert max(residuals.values()) < 1e-10
    overlap = {}
    for pattern,lo,hi in (("plug",.042,.29),("churn",.61,.85),("wavy annular",.86,.99)):
        vals=[r["alpha"] for r in selected if r["pattern"]==pattern]
        overlap[pattern]=dict(n=len(vals),inside_reported_literature_range=sum(lo<=v<=hi for v in vals),
                              literature_range=[lo,hi],note="Range overlap is not validation")
    unchanged = all(digest(p) == before[str(p.relative_to(SRC2))] for p in paths)
    assert unchanged
    result=dict(purpose="Legacy reproducibility and dataset-scope audit; not physical validation",
                coolprop_version=CoolProp.__version__,properties=props,
                original_source_hashes=before,original_sources_unchanged=unchanged,
                original_internal_check_messages=legacy_checks,
                independent_max_absolute_difference=max_error,cross_version_max_absolute_difference=cross_version,
                stored_table_comparisons=comparisons,selected_runs=[r["run"] for r in selected],
                excluded_runs=excluded,legacy_pattern_statistics=summary_original,
                selected_pattern_statistics=summary_selected,selected_runs_recomputed=selected,
                conservation_max_absolute_residuals=residuals,literature_range_overlap=overlap,
                src2_workbook_inventory=workbook_inventory,
                physical_accuracy="Not identifiable without matched quality/phase-velocity references; original mass-flow signals are documented separately by massflow_audit.py")
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"artifact":str(OUT),"original_rows":len(independent_rows),
                      "selected_rows":len(selected),"excluded_rows":len(excluded),
                      "max_independent_difference":max(max_error.values()),
                      "original_sources_unchanged":unchanged,
                      "selected_pattern_counts":{p:s["n"] for p,s in summary_selected.items()}},indent=2))


if __name__ == "__main__":
    main()
