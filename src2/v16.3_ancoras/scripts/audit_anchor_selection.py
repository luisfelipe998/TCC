"""Compare optical selection before/after the audit; never computes wavelets."""
import sys
sys.dont_write_bytecode = True
import copy
import importlib.util
import json
from itertools import combinations
from datetime import datetime, timezone

import select_anchors as current
from anchor_diagnostics import assess_anchor
from common import (CFG, INTER, OLD, ROOT, SCRIPTS, inventory, read_csv,
                    sha, verify_source, write_csv, write_json)

OUT = INTER / "auditoria_selecao"
BEFORE = INTER / "historico/antes_auditoria_selecao"


def make_group(module, source, pano, raw, reviews, previous):
    run, pattern = source["run"], source["pattern"]
    if pattern in ("plug", "slug"):
        rows, best = module.event_candidates(run, pattern, pano, raw, reviews, previous)
    else:
        rows, best = module.texture_candidates(run, pattern, raw)
    group = dict(pattern=pattern, rows=rows, best=best)
    if best is not None:
        module.common_window([group])
    return group


def main():
    previous_config = json.loads((BEFORE / "scripts/config.json").read_text())
    # Keep the historic detector/context experiment comparable while allowing
    # the explicitly introduced endpoint tie-break. Preserve the first audit.
    comparable = []
    for configuration in (previous_config, CFG):
        configuration = copy.deepcopy(configuration)
        configuration.pop("versao", None)
        configuration["eventos"].pop("nitidez_faixa_px", None)
        configuration["eventos"].pop("ranking", None)
        comparable.append(configuration)
    assert comparable[0] == comparable[1], "A configuracao alterou mais que o desempate; use uma nova referencia."
    config_changed = sha(SCRIPTS / "config.json") != sha(BEFORE / "scripts/config.json")
    output = OUT / CFG["versao"] if config_changed else OUT
    output.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location("selection_before_audit", BEFORE / "scripts/select_anchors.py")
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    reviews = {r["event_id"]: r for r in read_csv(OLD / "results/revisao_visual.csv")}
    previous = {r["event_id"]: r for r in read_csv(OLD / "results/catalogo_eventos.csv")}
    legacy_texture = {r["anchor_id"] for r in read_csv(OLD / "results/extensao_churn_anular/catalogo_ancoras.csv")
                      if r["selecionada"] == "True"}
    context_changes, runs, groups = [], [], {}
    eligible_checks = 0
    for source in inventory():
        run, pattern = source["run"], source["pattern"]
        pano, raw = verify_source(source)
        before = make_group(old, source, pano, raw, reviews, previous)
        after = make_group(current, source, pano, raw, reviews, previous)
        before_rows = {r["anchor_id"]: r for r in before["rows"]}
        for row in after["rows"]:
            prior = before_rows[row["anchor_id"]]
            if prior["comprimento_minimo"] != row["comprimento_minimo"]:
                context_changes.append(dict(aquisicao=run, padrao=pattern, anchor_id=row["anchor_id"],
                    minimo_antes=prior["comprimento_minimo"], minimo_depois=row["comprimento_minimo"],
                    elegivel_individual_antes=prior.get("elegivel_individual", prior["elegivel"])))
            if not row.get("elegivel_individual", row["elegivel"]):
                continue
            # Direct inequalities, independent of required_window's rounding formula.
            half = row["comprimento_minimo"] // 2
            margin = CFG["contexto"]["margem_evento" if row["tipo"] == "evento_delimitado" else "margem_vale"]
            if row["tipo"] == "evento_delimitado":
                starts = [row["inicio"]] + [row[f"inicio_limiar{s}"] for s in (55, 65)]
                stops = [row["fim"]] + [row[f"fim_limiar{s}"] for s in (55, 65)]
            else:
                starts = [row["base_esquerda"]] + [row[f"base_esquerda_sigma{s}"] for s in (12, 20)]
                stops = [row["base_direita"]] + [row[f"base_direita_sigma{s}"] for s in (12, 20)]
            assert min(starts) - (row["ancora"] - half) >= margin, row["anchor_id"]
            assert row["ancora"] + half - max(stops) >= margin, row["anchor_id"]
            eligible_checks += 1
        a, b = before["best"], after["best"]
        result = dict(aquisicao=run, padrao=pattern, candidatos=len(after["rows"]),
                      ancora_antes=a["ancora"] if a else None, ancora_depois=b["ancora"] if b else None,
                      janela_antes=a["comprimento"] if a else None, janela_depois=b["comprimento"] if b else None)
        result["selecao_preservada"] = (result["ancora_antes"], result["janela_antes"]) == (result["ancora_depois"], result["janela_depois"])
        if b:
            result.update(assess_anchor(b, pano, raw, after["rows"]))
            if pattern not in ("plug", "slug"):
                result["presente_no_catalogo_visual_v16_2"] = b["anchor_id"] in legacy_texture
                inner = current.match_minimum(b["ancora"], current.minima(pano[30:50].mean(axis=0), 16))
                result["ancora_proxima_faixa30_49"] = inner["ancora"] if inner else None
            groups[run] = after
        runs.append(result)
        print(f"{pattern}: {run} | ancora {result['ancora_antes']} -> {result['ancora_depois']} | "
              f"{result.get('alertas_opticos') or 'sem alerta'}", flush=True)

    pair_changes, checked_pairs = [], 0
    for pattern in CFG["padroes"]:
        names = sorted(run for run, group in groups.items() if group["pattern"] == pattern)
        for run_x, run_y in combinations(names, 2):
            paired = [copy.deepcopy(groups[run_x]), copy.deepcopy(groups[run_y])]
            length = current.common_window(paired)
            checked_pairs += 1
            for run, group in zip((run_x, run_y), paired):
                original, candidate = groups[run]["best"], group["best"]
                if candidate is None or (candidate["ancora"], length) != (original["ancora"], original["comprimento"]):
                    pair_changes.append(dict(run_x=run_x, run_y=run_y, aquisicao=run,
                        ancora_individual=original["ancora"], ancora_no_par=candidate["ancora"] if candidate else None,
                        comprimento=length))
    texture = [r for r in runs if r["padrao"] in ("churn", "wavy annular")]
    report = dict(data_utc=datetime.now(timezone.utc).isoformat(), status="PASS",
        escopo="Auditoria optica das aquisicoes; nenhuma wavelet nova nem expansao do piloto.",
        referencia=str(BEFORE.relative_to(ROOT)), configuracao_atual=CFG["versao"],
        desempate_por_contraste_ativo="nitidez_faixa_px" in CFG["eventos"],
        aquisicoes=len(runs), selecoes_preservadas=sum(r["selecao_preservada"] for r in runs),
        candidatos_com_contexto_conferido=eligible_checks,
        minimos_corrigidos=len(context_changes),
        minimos_corrigidos_antes_elegiveis=sum(r["elegivel_individual_antes"] for r in context_changes),
        alteracoes_contexto=context_changes,
        selecao_textura_repetivel=sum(r["repetivel_selecao_completa"] for r in texture),
        total_textura=len(texture),
        selecao_sensivel_sigma=[r["aquisicao"] for r in texture if not r["repetivel_selecao_completa"]],
        origens_conferidas=sum(r.get("coluna_origem_conferida", False) for r in runs),
        ancora_em_coluna_substituida=[r["aquisicao"] for r in runs if r.get("ancora_na_coluna_substituida_feather")],
        pares_com_contexto_conferido=checked_pairs, mudancas_dependentes_par=pair_changes,
        codigo_sha256={str(p.relative_to(ROOT)): sha(p) for p in
                       [SCRIPTS / "select_anchors.py", SCRIPTS / "anchor_diagnostics.py", SCRIPTS / "audit_anchor_selection.py"]})
    write_csv(output / "comparacao_aquisicoes.csv", runs)
    write_csv(output / "correcoes_contexto.csv", context_changes)
    write_json(output / "auditoria.json", report)
    print(f"PASS: {len(runs)} aquisicoes; {eligible_checks} contextos; "
          f"{report['minimos_corrigidos_antes_elegiveis']} minimos elegiveis corrigidos; "
          f"{report['selecao_textura_repetivel']}/{len(texture)} escolhas repetiveis.")


if __name__ == "__main__":
    main()
