"""Prepare individual data for future intra/inter-pattern comparisons; no pairs are calculated."""
import sys
sys.dont_write_bytecode = True
import argparse
import copy
import json
from datetime import datetime, timezone
from itertools import combinations

from common import (BASE, CFG, INTER, OLD, ROOT, SCRIPTS, context_ok, image_path,
                    inventory, read_csv, selected_manifest, sha, verify_source, write_csv, write_json)
from select_anchors import common_window, event_candidates, texture_candidates, event_rank
from anchor_diagnostics import assess_anchor
from wavelet_cache import (CACHE, context_record, get_auto, load_entry, signature,
                           individual_spec, validate_method, validate_length)

OUT = INTER / "preparacao_comparacoes"


def profile_context(record, profile, length, pano, raw, rows):
    """Keep the native anchor; independently assess the exact crop of each profile."""
    result = context_record(record, length)
    local = copy.deepcopy(rows)
    for row in local:
        row["elegivel"] = (row.get("elegivel_individual", row["elegivel"])
                           and row["comprimento_minimo"] <= length
                           and context_ok(row["ancora"], len(raw), length))
    key = event_rank if record["padrao"] in ("plug", "slug") else lambda r:(-r["prominencia"],r["variacao_maxima"],r["ancora"])
    ranked = sorted((r for r in local if r["elegivel"]), key=key)
    position = next(i for i,r in enumerate(ranked,1) if r["anchor_id"] == record["anchor_id"])
    result.update(assess_anchor(result, pano, raw, local))
    result.update(perfil=profile, ranking_nativo=record["ranking"], ranking=position,
                  primeiro_candidato_no_contexto=ranked[0]["anchor_id"],
                  politica_ancora="centro selecionado no contexto nativo, preservado entre perfis",
                  motivo="ancora_nativa_preservada" if profile=="comparavel" else record["motivo"])
    return result


def prepare_all(scope="piloto"):
    validate_method()
    if scope not in ("piloto", "completo"):
        raise ValueError("Escopo deve ser piloto ou completo.")
    OUT.mkdir(exist_ok=True)
    pilot = {r["aquisicao"]: r for r in selected_manifest()["ancoras"]}
    review_path, previous_path = OLD / "results/revisao_visual.csv", OLD / "results/catalogo_eventos.csv"
    reviews = {r["event_id"]: r for r in read_csv(review_path)}
    previous = {r["event_id"]: r for r in read_csv(previous_path)}
    sources = {str(p.relative_to(ROOT.parent)): sha(p) for p in
               (review_path, previous_path, BASE / "results/input_audit.csv")}
    sources[str((INTER/"selecao.json").relative_to(ROOT.parent))] = sha(INTER/"selecao.json")
    anchors, candidates, missing, contexts = [], [], [], {}
    length = validate_length(CFG["comprimento_referencia"])
    selected_sources = [r for r in inventory() if scope=="completo" or r["run"] in pilot]
    for source in selected_sources:
        run, pattern = source["run"], source["pattern"]
        pano, raw = verify_source(source)
        for path in (image_path(run), BASE / "intermediate" / f"{run}_input.npz"):
            sources[str(path.relative_to(ROOT.parent))] = sha(path)
        rows, best = event_candidates(run, pattern, pano, raw, reviews, previous) if pattern in ("plug", "slug") else texture_candidates(run, pattern, raw)
        if best is None:
            missing.append(dict(aquisicao=run, padrao=pattern, motivo="sem_ancora_valida"))
            candidates.extend(rows)
            continue
        group = dict(pattern=pattern, rows=rows, best=best)
        common_window([group])
        best = group["best"]
        best.update(assess_anchor(best, pano, raw, rows))
        sources[best["arquivo_quadro_origem"]] = best["quadro_origem_sha256"]
        if run in pilot:
            for key in ("anchor_id", "ancora", "comprimento", "segmento_inicio", "segmento_fim"):
                assert best[key] == pilot[run][key], (run, key)
        anchors.append(best)
        candidates.extend(rows)
        contexts[run] = {profile:profile_context(best,profile,n,pano,raw,rows)
                         for profile,n in (("nativo",best["comprimento"]),("comparavel",length))}
    if not anchors:
        raise ValueError("Nenhuma ancora disponivel para preparar comparacoes.")
    # A single domain across all patterns, preserving each already-selected center.
    if not all(r["comprimento_minimo"] <= length and context_ok(r["ancora"], r["tamanho_panorama"], length) for r in anchors):
        raise ValueError("As ancoras nao compartilham um contexto viavel; nao trocar centros nem redimensionar sinais.")
    if scope=="piloto" and set(contexts) != set(pilot):
        raise ValueError("A preparacao do piloto deve conter exatamente suas oito aquisicoes.")
    specs = {run:{profile:individual_spec(record) for profile,record in profiles.items()}
             for run,profiles in contexts.items()}
    frozen = dict(data_utc=datetime.now(timezone.utc).isoformat(), antes_dos_espectros=True, escopo=scope,
                  ancoras=anchors, perfis=contexts, especificacoes=specs,
                  comprimento_comparavel=length, fontes_sha256=sources,
                  codigo_sha256={str(p.relative_to(ROOT)): sha(p) for p in
                       [SCRIPTS / name for name in ("config.json", "common.py", "select_anchors.py", "anchor_diagnostics.py",
                                                    "wavelet_core.py", "wavelet_cache.py", "prepare_comparisons.py")]})
    write_json(OUT / "congelamento_individual.json", frozen)
    write_csv(OUT / "ancoras.csv", anchors)
    write_csv(OUT / "candidatos.csv", candidates)
    write_csv(OUT / "contextos_por_perfil.csv", [r for profiles in contexts.values() for r in profiles.values()])
    write_json(OUT / "ausencias.json", missing)
    entries = {}
    for record in anchors:
        run = record["aquisicao"]
        profiles = {}
        for profile, context in contexts[run].items():
            auto = get_auto(context)
            if auto["meta"]["spec"] != specs[run][profile]:
                raise ValueError("O precalculo difere do perfil congelado antes dos espectros.")
            profiles[profile] = dict(cache_key=auto["key"], grade=auto["meta"]["spec"]["grade"],
                                     contexto=context,
                                     metadata_sha256=sha(CACHE / "individuais" / auto["key"] / "metadata.json"))
        entries[run] = dict(padrao=record["padrao"], anchor_id=record["anchor_id"], perfis=profiles)
        print(f"Individual pronto: {record['padrao']} / {run} / {record['comprimento']} e {length} px", flush=True)
    plans = []
    for a, b in combinations(anchors, 2):
        same = a["padrao"] == b["padrao"]
        plans.append(dict(pair_id=f"{a['aquisicao']}__{b['aquisicao']}",
                          run_x=a["aquisicao"], run_y=b["aquisicao"],
                          pattern_x=a["padrao"], pattern_y=b["padrao"],
                          escopo="intra" if same else "inter", perfil_comparavel=True,
                          perfil_nativo_compativel=a["comprimento"] == b["comprimento"],
                          calculado_nesta_preparacao=False))
    write_csv(OUT / "pares_planejados.csv", plans)
    unique = {p["cache_key"] for entry in entries.values() for p in entry["perfis"].values()}
    keys = {p["grade"] for entry in entries.values() for p in entry["perfis"].values()}
    byte_size = sum(p.stat().st_size for key in unique for p in (CACHE / "individuais" / key).iterdir())
    byte_size += sum(p.stat().st_size for key in keys for p in (CACHE / "grades" / key).iterdir())
    index = dict(data_utc=datetime.now(timezone.utc).isoformat(), estado="somente_individuais_precalculados", escopo=scope,
                 congelamento_sha256=sha(OUT / "congelamento_individual.json"), aquisicoes=entries,
                 perfil_comparavel=dict(comprimento=length, lambda_min=CFG["wavelet"]["comprimento_onda_minimo"],
                    lambda_max=length/2, escalas=CFG["wavelet"]["escalas"], dx_px=1,
                    nota="Usar este mesmo perfil em intra e inter para comparar resultados; nativo preserva as janelas por padrao."),
                 individuais_unicos=len(unique), grades=len(keys), bytes_precalculados=byte_size,
                 pares_intra_planejados=sum(r["escopo"] == "intra" for r in plans),
                 pares_inter_planejados=sum(r["escopo"] == "inter" for r in plans), pares_calculados_nesta_etapa=0)
    write_json(OUT / "indice.json", index)
    print(f"Pronto: {len(anchors)} aquisicoes, {len(unique)} contextos individuais, {len(keys)} grades; nenhuma comparacao nova.")
    return index


def _frozen_catalog(index):
    path = OUT / "congelamento_individual.json"
    if sha(path) != index["congelamento_sha256"]:
        raise ValueError("O catalogo individual foi alterado.")
    return json.loads(path.read_text())


def _validate_binding(catalog, frozen, run, profile):
    try:
        entry = catalog["aquisicoes"][run]
        item = entry["perfis"][profile]
        expected_context = frozen["perfis"][run][profile]
        expected_spec = frozen["especificacoes"][run][profile]
    except KeyError as error:
        raise ValueError("Ensaio/perfil ausente do congelamento.") from error
    if (entry["padrao"] != expected_context["padrao"] or entry["anchor_id"] != expected_context["anchor_id"]
            or expected_context["aquisicao"] != run or item["contexto"] != expected_context):
        raise ValueError("Vinculo incompativel entre ensaio, ancora e contexto.")
    if (individual_spec(expected_context) != expected_spec or item["cache_key"] != signature(expected_spec)
            or item["grade"] != expected_spec["grade"]):
        raise ValueError("Vinculo incompativel entre fonte, parametros e cache.")
    if expected_context["diagnosticos_comprimento"] != expected_context["comprimento"]:
        raise ValueError("Diagnosticos pertencem a outra janela.")
    return item, expected_spec


def load_catalog():
    validate_method()
    index = json.loads((OUT / "indice.json").read_text())
    frozen = _frozen_catalog(index)
    for mapping, root in ((frozen["codigo_sha256"], ROOT), (frozen["fontes_sha256"], ROOT.parent)):
        for relative, expected in mapping.items():
            if sha(root / relative) != expected:
                raise ValueError(f"Fonte/parametro alterado; regenere a preparacao: {relative}")
    if set(index["aquisicoes"]) != set(frozen["perfis"]) or index["escopo"] != frozen["escopo"]:
        raise ValueError("Escopo ou inventario do indice difere do congelamento.")
    for run in frozen["perfis"]:
        if set(index["aquisicoes"][run]["perfis"]) != set(frozen["perfis"][run]):
            raise ValueError("Perfis divergentes do congelamento.")
        for profile in frozen["perfis"][run]:
            _validate_binding(index, frozen, run, profile)
    return index


def load_individual(catalog, run, profile="comparavel"):
    validate_method()
    item, expected_spec = _validate_binding(catalog, _frozen_catalog(catalog), run, profile)
    folder = CACHE / "individuais" / item["cache_key"]
    if sha(folder / "metadata.json") != item["metadata_sha256"]:
        raise ValueError("Metadados individuais alterados.")
    auto = load_entry(str(folder))
    if auto["meta"]["spec"] != expected_spec:
        raise ValueError("O conteudo carregado nao corresponde ao ensaio/perfil solicitado.")
    return copy.deepcopy(item["contexto"]), auto


def comparison_spec(catalog, run_x, run_y, profile="comparavel"):
    """Build consistent metadata and image contexts; does not calculate a comparison."""
    from wavelet_cache import compatible_grid
    if run_x == run_y:
        raise ValueError("Os pares principais devem usar aquisicoes diferentes.")
    a, auto_a = load_individual(catalog, run_x, profile)
    b, auto_b = load_individual(catalog, run_y, profile)
    compatible_grid(auto_a, auto_b)
    pattern = a["padrao"] if a["padrao"] == b["padrao"] else f"{a['padrao']} × {b['padrao']}"
    spec = dict(padrao=pattern, pattern_x=a["padrao"], pattern_y=b["padrao"], run_x=run_x, run_y=run_y,
                anchor_x=a["anchor_id"], anchor_y=b["anchor_id"], comprimento=a["comprimento"],
                perfil=profile, escopo="intra" if a["padrao"] == b["padrao"] else "inter")
    return spec, {a["anchor_id"]: a, b["anchor_id"]: b}, (auto_a, auto_b)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--escopo", choices=["piloto", "completo"], default="piloto",
                        help="Completo somente quando a expansao for solicitada; nunca calcula pares.")
    prepare_all(parser.parse_args().escopo)
