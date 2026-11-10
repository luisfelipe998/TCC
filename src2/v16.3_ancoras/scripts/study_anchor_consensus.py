"""Evaluate a proposed ranking; never updates the selected anchors or wavelets."""
import sys
sys.dont_write_bytecode = True
from datetime import datetime, timezone

from common import (ANCHORS, CFG, INTER, SCRIPTS, WAVELETS, inventory,
                    selected_manifest, sha, verify_source, write_csv, write_json)
from select_anchors import common_window, match_minimum, minima, texture_candidates


def main():
    pilot = selected_manifest()
    pilot_runs = {r["aquisicao"] for r in pilot["ancoras"]}
    protected = [INTER / "selecao.json", INTER / "congelamento.json",
                 SCRIPTS / "select_anchors.py", SCRIPTS / "anchor_diagnostics.py"]
    protected += list(ANCHORS.glob("*.png")) + list(WAVELETS.glob("*.png"))
    protected += list((INTER / "mapas").glob("*.npz"))
    hashes = {str(p): sha(p) for p in protected}
    sigmas = sorted([CFG["textura"]["sigma"], *CFG["textura"]["sigmas_estabilidade"]])
    candidates, comparisons = [], []
    for source in inventory():
        pattern, run = source["pattern"], source["run"]
        if pattern not in ("churn", "wavy annular"):
            continue
        _, raw = verify_source(source)
        rows, best = texture_candidates(run, pattern, raw)
        if best is None:
            continue
        group = dict(pattern=pattern, rows=rows, best=best)
        length = common_window([group])
        best = group["best"]
        detected = {sigma: minima(raw, sigma) for sigma in sigmas}
        local = []
        for row in rows:
            if not row["elegivel"] or row["comprimento_minimo"] > length:
                continue
            record = dict(aquisicao=run, padrao=pattern, anchor_id=row["anchor_id"],
                          ancora=row["ancora"], variacao_maxima=row["variacao_maxima"],
                          comprimento=length)
            for sigma in sigmas:
                matched = match_minimum(row["ancora"], detected[sigma])
                assert matched is not None
                record[f"prominencia_sigma{sigma}"] = matched["prominencia"]
            local.append(record)
        maxima = {s: max(r[f"prominencia_sigma{s}"] for r in local) for s in sigmas}
        for record in local:
            for sigma in sigmas:
                record[f"destaque_relativo_sigma{sigma}"] = record[f"prominencia_sigma{sigma}"] / maxima[sigma]
            record["pior_destaque_relativo"] = min(record[f"destaque_relativo_sigma{s}"] for s in sigmas)
        ordered = sorted(local, key=lambda r: (-r["pior_destaque_relativo"], r["variacao_maxima"],
                                               -r["prominencia_sigma16"], r["ancora"]))
        proposal = ordered[0]
        original = next(r for r in local if r["anchor_id"] == best["anchor_id"])
        comparison = dict(aquisicao=run, padrao=pattern, no_piloto=run in pilot_runs,
                          ancora_atual=best["ancora"], ancora_proposta=proposal["ancora"],
                          muda_ancora=best["ancora"] != proposal["ancora"],
                          pior_destaque_atual=original["pior_destaque_relativo"],
                          pior_destaque_proposto=proposal["pior_destaque_relativo"],
                          vantagem_sobre_segundo=proposal["pior_destaque_relativo"]-ordered[1]["pior_destaque_relativo"] if len(ordered)>1 else None)
        candidates.extend(local)
        comparisons.append(comparison)
        if comparison["muda_ancora"]:
            print(f"{pattern}: {run} | {best['ancora']} -> {proposal['ancora']} | "
                  f"pior destaque {comparison['pior_destaque_atual']:.4f} -> {comparison['pior_destaque_proposto']:.4f}")
    assert hashes == {str(p): sha(p) for p in protected}
    out = INTER / "estudo_consenso"
    out.mkdir(exist_ok=True)
    write_csv(out / "candidatos.csv", candidates)
    write_csv(out / "comparacao.csv", comparisons)
    write_json(out / "resumo.json", dict(data_utc=datetime.now(timezone.utc).isoformat(),
        estado="proposta_diagnostica_nao_aplicada", aquisicoes=len(comparisons),
        mudancas=sum(r["muda_ancora"] for r in comparisons),
        mudancas_no_piloto=sum(r["muda_ancora"] and r["no_piloto"] for r in comparisons),
        regra="Maximizar min_s(P_i(s)/max_j(P_j(s))); s=12,16,20; mesmo conjunto de candidatos elegiveis. "
              "Desempates: menor variacao de posicao, maior prominencia em sigma16, menor coordenada.",
        limite="Estudo de efeito da regra, sem validacao independente de melhoria e sem uso de WCO. "
               "O destaque relativo nao e probabilidade de acerto; o centro permanece o de sigma16.",
        arquivos_protegidos_preservados=True, script_sha256=sha(SCRIPTS / "study_anchor_consensus.py")))
    print(f"Estudo: {len(comparisons)} aquisicoes; {sum(r['muda_ancora'] for r in comparisons)} mudancas; "
          f"{sum(r['muda_ancora'] and r['no_piloto'] for r in comparisons)} no piloto. Seletor e resultados preservados.")


if __name__ == "__main__":
    main()
