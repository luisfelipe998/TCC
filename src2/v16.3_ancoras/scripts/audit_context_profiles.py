"""Recompute optical diagnostics per cached context without changing any choice."""
import sys
sys.dont_write_bytecode = True
import copy
import json
import numpy as np

from common import INTER, OLD, context_ok, load_image, read_csv, write_csv, write_json
from select_anchors import event_candidates, texture_candidates, common_window, event_rank
from anchor_diagnostics import assess_anchor
from prepare_comparisons import load_catalog, load_individual
from audit_method_schuck import OUT, independent_pair, independent_prepare


def main():
    catalog = load_catalog()
    reviews = {r["event_id"]:r for r in read_csv(OLD/"results/revisao_visual.csv")}
    previous = {r["event_id"]:r for r in read_csv(OLD/"results/catalogo_eventos.csv")}
    changes, ranking = [], []
    for run, entry in catalog["aquisicoes"].items():
        pano,raw = load_image(run)
        pattern = entry["padrao"]
        rows,best = event_candidates(run,pattern,pano,raw,reviews,previous) if pattern in ("plug","slug") else texture_candidates(run,pattern,raw)
        group = dict(pattern=pattern,rows=rows,best=best)
        common_window([group])
        for profile in ("nativo","comparavel"):
            record,_ = load_individual(catalog,run,profile)
            n = record["comprimento"]
            local_rows = copy.deepcopy(rows)
            for row in local_rows:
                row["elegivel"] = row["elegivel_individual"] and row["comprimento_minimo"] <= n and context_ok(row["ancora"],len(raw),n)
            fresh = assess_anchor(record,pano,raw,local_rows)
            for field,actual in fresh.items():
                if actual != record[field]:
                    changes.append(dict(aquisicao=run,perfil=profile,comprimento=n,campo=field,
                                        armazenado=record[field],recalculado=actual))
            key = event_rank if pattern in ("plug","slug") else lambda r:(-r["prominencia"],r["variacao_maxima"],r["ancora"])
            winner = min((r for r in local_rows if r["elegivel"]),key=key)
            ranking.append(dict(aquisicao=run,padrao=pattern,perfil=profile,comprimento=n,
                                ancora_congelada=record["ancora"],melhor_candidato_neste_contexto=winner["ancora"],
                                escolha_mesma=winner["anchor_id"]==record["anchor_id"],
                                repetivel_sigma_recalculado=fresh.get("repetivel_selecao_completa")))
    # Resolve whether the literal-estimator excesses also occur in the central
    # lags, rather than only far outside the summarized domain.
    _,a = load_individual(catalog,"m03q05_001","nativo")
    _,b = load_individual(catalog,"m04q10_001","nativo")
    n = len(a["raw"])
    maps = independent_pair(independent_prepare(a["raw"]),independent_prepare(b["raw"]),np.geomspace(50,n/2,200),unbiased=True)
    central = maps["wco"][:,n-33:n+32]
    write_csv(OUT/"diagnosticos_por_perfil.csv",changes)
    write_csv(OUT/"ranking_por_perfil.csv",ranking)
    report = dict(diagnosticos_divergentes=len(changes),
                  mudancas_alem_contagens=[r for r in changes if r["campo"] not in (
                      "costuras_contexto","colunas_substituidas_feather_contexto","outros_candidatos_estaveis_no_contexto","ids_vizinhos_no_contexto")],
                  melhor_ancora_diferente_no_contexto_reduzido=[r for r in ranking if not r["escolha_mesma"]],
                  churn_estimador_literal_centro=dict(maximo=float(np.nanmax(central)),
                       celulas_acima_1=int(np.sum(central>1+1e-10)),total_celulas=central.size),
                  observacao="Somente auditoria: nenhuma ancora, metadado de producao ou figura alterados.")
    write_json(OUT/"auditoria_contextos.json",report)
    print(json.dumps(report,ensure_ascii=False,indent=2))


if __name__ == "__main__":
    main()
