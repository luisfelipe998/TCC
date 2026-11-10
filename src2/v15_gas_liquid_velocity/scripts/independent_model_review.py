#!/usr/bin/env python3
"""Independent verification of saved model results; no model outputs are rewritten."""
from __future__ import annotations
import ast
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import statistics

from CoolProp.CoolProp import PropsSI
import numpy as np
import pandas as pd

BASE=Path(__file__).resolve().parents[1]
ROOT=BASE.parents[1]
RESULT=BASE/"results"
AREA=math.pi*.001**2/4
MODELS=("homogeneous","rouhani_axelsson","xu_fang_literal_2021","xu_fang_fluids")
SCENARIOS=("saturated_inlet","subcool_1K","subcool_2K","subcool_5K","subcool_10K",
           "inlet_quality_0.01","inlet_quality_0.05","Tsat_25C","Tsat_29C","net_heat_90pct",
           "recorded_flow_mean","recorded_flow_median","recorded_flow_trimmed10")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def extract_functions(path,names,namespace):
    tree=ast.parse(path.read_text())
    selected=[node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name in names]
    assert {node.name for node in selected}==set(names)
    for node in selected:
        assert not any(isinstance(child,(ast.Import,ast.ImportFrom)) for child in ast.walk(node))
    exec(compile(ast.Module(body=selected,type_ignores=[]),str(path),"exec"),namespace)
    return namespace


def state(scenario,run,flows):
    g={3:380.,4:510.,5:640.,6:760.,7:890.}[int(run[1:3])]
    q=float(run[4:6])*1000
    t=float(scenario[5:-1]) if scenario.startswith("Tsat_") else 27.
    p=PropsSI("P","T",t+273.15,"Q",0,"R600a")
    hliq=PropsSI("H","T",t+273.15,"Q",0,"R600a")
    hgas=PropsSI("H","T",t+273.15,"Q",1,"R600a")
    rl=PropsSI("D","T",t+273.15,"Q",0,"R600a")
    rg=PropsSI("D","T",t+273.15,"Q",1,"R600a")
    sigma=PropsSI("I","T",t+273.15,"Q",0,"R600a")
    hin=hliq
    if scenario.startswith("subcool_"):
        delta=float(scenario.split("_")[1][:-1])
        hin=PropsSI("H","P",p,"T",t+273.15-delta,"R600a")
    if scenario.startswith("inlet_quality_"):
        hin=hliq+float(scenario.split("_")[-1])*(hgas-hliq)
    if scenario.startswith("recorded_flow_"):
        suffix=scenario.replace("recorded_flow_","")
        g=float(flows.loc[run,"G_"+suffix+"_kg_m2_s"])
    power=q*math.pi*.001*.265*(.9 if scenario=="net_heat_90pct" else 1.)
    mdot=g*AREA
    hout=hin+power/mdot
    return dict(G=g,q=q,t=t,rl=rl,rg=rg,sigma=sigma,hin=hin,hout=hout,hl=hliq,hv=hgas,
                x=(hout-hliq)/(hgas-hliq),xin=(hin-hliq)/(hgas-hliq),
                dx=power/mdot/(hgas-hliq),power=power,mdot=mdot,
                Fr=(g/rl)**2/(9.81*.001))


def independent_phase(s,model,secondary):
    x,rl,rg,g=s["x"],s["rl"],s["rg"],s["G"]
    if not 0<x<1:
        return None
    vapor_mass=s["mdot"]*x
    liquid_mass=s["mdot"]-vapor_mass
    jg=vapor_mass/(rg*AREA)
    jl=liquid_mass/(rl*AREA)
    if model=="homogeneous":
        a=secondary["homogeneous"](x,rl,rg)
        ug=ul=jg+jl
    elif model=="rouhani_axelsson":
        vd=1.18*(1-x)*(9.81*s["sigma"]*(rl-rg))**.25/math.sqrt(rl)
        ug=(1+.12*(1-x))*(jg+jl)+vd
        a=jg/ug
        ul=jl/(1-a)
    elif model=="xu_fang_fluids":
        a=secondary["Xu_Fang_voidage"](x,rl,rg,s["mdot"],.001,g=9.81)
        ug=jg/a
        ul=jl/(1-a)
    elif model=="xu_fang_literal_2021":
        ah=secondary["homogeneous"](x,rl,rg)
        slip=1+2*s["Fr"]**(-2)*ah**3.5
        ug=jg+slip*jl
        ul=jg/slip+jl
        a=jg/ug
    return dict(Jg=jg,Jl=jl,alpha=a,ug=ug,ul=ul,slip=ug/ul)


def main():
    paths=[RESULT/"thermodynamic_sensitivity.csv",RESULT/"conditional_velocities_per_run.csv",
           RESULT/"recorded_flow_conditional_velocities_per_run.csv",
           RESULT/"conditional_velocities_by_pattern.csv",RESULT/"recorded_flow_conditional_velocities_by_pattern.csv",
           RESULT/"sensitivity_relative_to_baseline.csv",RESULT/"closure_spread.csv",
           RESULT/"image_alpha_diagnostics.csv",BASE/"scripts/model_comparison.py"]
    before={str(path.relative_to(BASE)):sha(path) for path in paths}
    source=BASE/"literature/extracts/fluids_two_phase_voidage.py"
    secondary=extract_functions(source,["homogeneous","Xu_Fang_voidage"],{"pi":math.pi,"g":9.81})
    selected={}
    excluded=[]
    for folder in sorted((ROOT/"src2/data_processed").glob("m??q??_???")):
        marker=list(folder.glob("__*.txt"))
        assert len(marker)==1
        label=marker[0].stem[2:].replace("_"," ")
        if folder.name=="m05q05_002" or label=="smooth annular" or "-" in label:
            excluded.append(folder.name)
        else:
            selected[folder.name]=label
    assert len(selected)==26 and len(excluded)==5
    table=pd.read_csv(paths[0])
    assert len(table)==26*13*4 and not table.duplicated(["run","scenario","model"]).any()
    assert set(table.run)==set(selected) and set(table.scenario)==set(SCENARIOS) and set(table.model)==set(MODELS)
    expected_keys={(run,scenario,model) for run in selected for scenario in SCENARIOS for model in MODELS}
    assert set(zip(table.run,table.scenario,table.model))==expected_keys
    flows=pd.read_csv(RESULT/"massflow_per_run.csv").set_index("run")
    assert set(flows.index)==set(selected)
    errors=defaultdict(float)
    nchecked=0
    expected_values={}
    invalid=Counter()
    flux_residual=energy_residual=0.
    for row in table.to_dict("records"):
        assert row["pattern"]==selected[row["run"]]
        s=state(row["scenario"],row["run"],flows)
        expected=dict(G_kg_m2_s=s["G"],q_kw_m2=s["q"]/1000,temperature_c=s["t"],
                      rho_l=s["rl"],rho_g=s["rg"],inlet_equilibrium_quality=s["xin"],
                      quality_increment=s["dx"],x_eq=s["x"],Fr_lo=s["Fr"])
        assert bool(row["XF_secondary_Fr_domain"])==(.02<=s["Fr"]<=145)
        phase=independent_phase(s,row["model"],secondary)
        if phase is None:
            assert row["quality_status"]=="outside_two_phase_equilibrium_domain"
            assert all(math.isnan(row[k]) for k in ("Jg","Jl","alpha","ug","ul","slip"))
            invalid[row["scenario"]]+=1
        else:
            assert row["quality_status"]=="two_phase"
            expected.update(phase)
            # Independently recover phase mass rates from saved phase speeds/areas.
            actual_mdot=AREA*(s["rg"]*row["alpha"]*row["ug"]+s["rl"]*(1-row["alpha"])*row["ul"])
            flux_residual=max(flux_residual,abs(actual_mdot-s["mdot"]))
        # Full steady enthalpy balance, including negative equilibrium quality.
        actual_hout=s["hl"]+row["x_eq"]*(s["hv"]-s["hl"])
        energy_residual=max(energy_residual,abs(s["mdot"]*(actual_hout-s["hin"])-s["power"]))
        for k,v in expected.items():
            diff=abs(row[k]-v)
            errors[k]=max(errors[k],diff)
            assert math.isclose(row[k],v,rel_tol=5e-12,abs_tol=5e-12),(row["run"],row["scenario"],row["model"],k,row[k],v)
            nchecked+=1
        expected_values[row["run"],row["scenario"],row["model"]]=expected
    assert flux_residual<1e-15 and energy_residual<1e-10
    for path,scenario in ((paths[1],"saturated_inlet"),(paths[2],"recorded_flow_mean")):
        saved=pd.read_csv(path).sort_values(["run","model"]).reset_index(drop=True)
        expected=table[table.scenario==scenario].sort_values(["run","model"]).reset_index(drop=True)
        pd.testing.assert_frame_equal(saved,expected)
    summary_cells=0
    for path,scenario in ((paths[3],"saturated_inlet"),(paths[4],"recorded_flow_mean")):
        summary=pd.read_csv(path)
        assert len(summary)==20 and not summary.duplicated(["model","pattern"]).any()
        for row in summary.to_dict("records"):
            runs=[run for run,label in selected.items() if label==row["pattern"]]
            assert row["n"]==len(runs)
            values=[expected_values[run,scenario,row["model"]] for run in runs]
            for key in ("G_kg_m2_s","x_eq","Jg","Jl","alpha","ug","ul","slip"):
                if key+"_mean" not in row:
                    continue
                vals=[v[key] for v in values]
                for stat,fn in (("min",min),("max",max),("mean",statistics.mean),("median",statistics.median)):
                    assert math.isclose(row[key+"_"+stat],fn(vals),rel_tol=5e-12,abs_tol=5e-12)
                    summary_cells+=1
    legacy=json.loads((RESULT/"validation/legacy_audit.json").read_text())
    legacy_error=0.
    for r in legacy["selected_runs_recomputed"]:
        here=expected_values[r["run"],"saturated_inlet","rouhani_axelsson"]
        for k,kold in (("x_eq","x"),("Jg","Jg"),("Jl","Jl"),("alpha","alpha"),("ug","uv"),("ul","ul"),("slip","S")):
            legacy_error=max(legacy_error,abs(here[k]-r[kold]))
    assert legacy_error<1e-11
    closure=pd.read_csv(paths[6])
    assert len(closure)==26 and set(closure.run)==set(selected)
    for r in closure.to_dict("records"):
        vals=[expected_values[r["run"],"saturated_inlet",model] for model in MODELS]
        for key in ("ug","ul","alpha"):
            for stat,fn in (("min",min),("max",max)):
                expected=fn(v[key] for v in vals)
                assert math.isclose(r[key+"_"+stat],expected,rel_tol=5e-12,abs_tol=5e-12)
        for key in ("ug","ul"):
            assert math.isclose(r[key+"_max_to_min"],max(v[key] for v in vals)/min(v[key] for v in vals),rel_tol=5e-12)
    hybrid=pd.read_csv(paths[7])
    snapshot=pd.read_csv(BASE/"intermediate/v14_void_input_snapshot.csv").set_index("run")
    assert len(hybrid)==26*3 and not hybrid.duplicated(["run","variant"]).any()
    hybrid_supported=0
    inverse_residual=0.
    for r in hybrid.to_dict("records"):
        col={"D60":"conditional_circle_alpha_spatial_mean","D56":"conditional_alpha_D56","D64":"conditional_alpha_D64"}[r["variant"]]
        a=float(snapshot.loc[r["run"],col])
        assert math.isclose(a,r["alpha_image"],rel_tol=1e-12)
        supported=snapshot.loc[r["run"],"validity"]!="unsupported_saturation" and 0<a<1
        if not supported:
            assert all(math.isnan(r[k]) for k in ("ug","ul","slip","RA_implied_x","RA_implied_x_in"))
            continue
        hybrid_supported+=1
        ref=expected_values[r["run"],"saturated_inlet","rouhani_axelsson"]
        assert math.isclose(a*r["ug"],ref["Jg"],rel_tol=5e-12)
        assert math.isclose((1-a)*r["ul"],ref["Jl"],rel_tol=5e-12)
        assert math.isclose(r["slip"],r["ug"]/r["ul"],rel_tol=5e-12)
        assert math.isclose(r["RA_implied_x_in"],r["RA_implied_x"]-ref["x_eq"],abs_tol=5e-12)
        st=state("saturated_inlet",r["run"],flows)
        st["x"]=r["RA_implied_x"]
        pred=independent_phase(st,"rouhani_axelsson",secondary)["alpha"]
        inverse_residual=max(inverse_residual,abs(a-pred))
    assert inverse_residual<1e-9
    sensitivity=pd.read_csv(paths[5])
    assert len(sensitivity)==26*13 and (sensitivity.model=="rouhani_axelsson").all()
    for r in sensitivity.to_dict("records"):
        ref=expected_values[r["run"],"saturated_inlet","rouhani_axelsson"]
        here=expected_values[r["run"],r["scenario"],"rouhani_axelsson"]
        for key in ("Jg","Jl","ug","ul"):
            if key not in here:
                assert math.isnan(r[key+"_change_pct"])
            else:
                assert math.isclose(r[key+"_change_pct"],100*(here[key]/ref[key]-1),rel_tol=1e-10,abs_tol=1e-10)
    flow_invariance={}
    for stat in ("mean","median","trimmed10"):
        diffs=[]
        for run in selected:
            a=expected_values[run,"saturated_inlet","rouhani_axelsson"]
            b=expected_values[run,"recorded_flow_"+stat,"rouhani_axelsson"]
            diffs.append(abs(a["Jg"]-b["Jg"]))
            assert math.isclose(b["Jl"]-a["Jl"],(b["G_kg_m2_s"]-a["G_kg_m2_s"])/a["rho_l"],abs_tol=1e-12)
        flow_invariance[stat]=max(diffs)
    root_funcs=extract_functions(BASE/"scripts/model_comparison.py",["alpha_model","velocities"],
                                 {"np":np,"GRAVITY":9.81,"D":.001})
    s=state("saturated_inlet",next(iter(selected)),flows)
    p={"rho_l":s["rl"],"rho_g":s["rg"],"sigma":s["sigma"]}
    for model in MODELS:
        for x in (-.01,1.01):
            v=root_funcs["velocities"](x,640,p,model)
            assert all(math.isnan(v[k]) for k in ("Jg","Jl","alpha","ug","ul","slip"))
        zero=root_funcs["velocities"](0,640,p,model)
        one=root_funcs["velocities"](1,640,p,model)
        assert math.isnan(zero["ug"]) and math.isnan(one["ul"])
        assert math.isclose(zero["ul"],640/p["rho_l"]) and math.isclose(one["ug"],640/p["rho_g"])
    metadata=json.loads((RESULT/"validation/model_metadata.json").read_text())
    assert metadata["massflow_input"]["sha256"]==sha(RESULT/"massflow_per_run.csv")
    assert metadata["image_input"]["snapshot_sha256"]==sha(BASE/metadata["image_input"]["snapshot"])
    assert metadata["script_sha256"]==sha(BASE/"scripts/model_comparison.py")
    narrative_cells=0
    narrative_hashes={}
    executive=BASE/"findings/conditional_model_summary.md"
    velocity_report=BASE/"findings/velocity_results.md"
    if executive.exists():
        narrative_hashes[str(executive.relative_to(BASE))]=sha(executive)
        bypattern_nominal=pd.read_csv(paths[3]).query("model == 'rouhani_axelsson'").set_index("pattern")
        bypattern_recorded=pd.read_csv(paths[4]).query("model == 'rouhani_axelsson'").set_index("pattern")
        for line in executive.read_text().splitlines():
            cells=[c.strip() for c in line.strip("|").split("|")]
            if cells[0] not in set(selected.values()) or len(cells)!=7:
                continue
            a,b=bypattern_nominal.loc[cells[0]],bypattern_recorded.loc[cells[0]]
            assert int(cells[1])==int(a['n'])==int(b['n'])
            for cell,value in zip(cells[2:],[a.Jg_median,a.Jl_median,b.Jl_median,b.ug_median,b.ul_median]):
                assert cell==f"{value:.3f}",("historical model summary table",cells[0],cell,value)
                narrative_cells+=1
        assert narrative_cells==25
    if velocity_report.exists():
        narrative_hashes[str(velocity_report.relative_to(BASE))]=sha(velocity_report)
        sections=velocity_report.read_text().split("## ")
        keys=(("G_kg_m2_s",2),("x_eq",5),("Jg",3),("Jl",3),("alpha",4),("ug",3),("ul",3),("slip",3))
        for scenario,prefix in (("recorded_flow_mean","Recorded-mean-flow alternative"),("saturated_inlet","Nominal-G baseline")):
            section=next(part for part in sections if part.startswith(prefix))
            stored=table[(table.scenario==scenario)&(table.model=="rouhani_axelsson")].set_index("run")
            count=0
            for line in section.splitlines():
                if not line.startswith("| m"):
                    continue
                cells=[c.strip() for c in line.strip("|").split("|")]
                r=stored.loc[cells[0]]
                assert r['pattern']==cells[1]
                for cell,(key,digits) in zip(cells[2:],keys):
                    wanted=f"{r[key]:.0f}" if key=="G_kg_m2_s" and scenario=="saturated_inlet" else f"{r[key]:.{digits}f}"
                    assert cell==wanted,(scenario,cells[0],key,cell,wanted)
                    narrative_cells+=1
                count+=1
            assert count==26
    method_report=BASE/"findings/method_comparison.md"
    if method_report.exists():
        narrative_hashes[str(method_report.relative_to(BASE))]=sha(method_report)
    assert all(sha(BASE/path)==digest for path,digest in before.items()),"Inputs changed during audit; rerun"
    notes=[]
    if "temperature_and_G_status" in metadata:
        notes.append("Metadata temperature_and_G_status applies to the nominal baseline; recorded-flow scenario rows have their own G_basis and massflow provenance. Prefer explicit split labels.")
    result=dict(status="passed",scope_count=len(selected),excluded=excluded,
                scenarios=list(SCENARIOS),models=list(MODELS),experiment_rows=len(table),
                numeric_cells_compared=nchecked,summary_cells_compared=summary_cells,
                narrative_numeric_cells_compared=narrative_cells,
                reviewed_narrative_sha256=narrative_hashes,
                max_absolute_differences=dict(errors),
                mass_flow_from_phase_fields_max_residual_kg_s=flux_residual,
                full_enthalpy_balance_max_residual_W=energy_residual,
                legacy_RA_max_absolute_difference=legacy_error,
                hybrid_rows=len(hybrid),hybrid_supported_rows=hybrid_supported,
                inverse_RA_max_void_fraction_residual=inverse_residual,
                invalid_two_phase_acquisitions_by_scenario={k:v//4 for k,v in invalid.items()},
                recorded_G_Jg_invariance_max_difference=flow_invariance,
                secondary_source=dict(file=str(source.relative_to(BASE)),sha256=sha(source),
                                      functions=["homogeneous","Xu_Fang_voidage"],
                                      gravity_explicitly_overridden_m_s2=9.81,
                                      evaluation="Only the two inspected AST function definitions were executed; full source not imported"),
                reviewed_input_sha256=before,minor_metadata_notes=notes,
                limitation="Arithmetic and model implementation verified; no independent experimental phase-velocity accuracy identified")
    (RESULT/"validation/independent_model_review.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k not in ("reviewed_input_sha256","secondary_source")},indent=2))


if __name__=="__main__":
    main()
