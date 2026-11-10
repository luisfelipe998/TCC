"""Optical checks omitted from the first v16.3 pilot; never uses wavelet results."""
import numpy as np
from PIL import Image

from common import BASE, CFG, DATA, ROOT, context_ok, sha


def context_diagnostic_fields():
    """Fields which must never be inherited by a differently sized crop."""
    fields = {"costuras_contexto", "colunas_substituidas_feather_contexto",
              "deslocamentos_nao_positivos_contexto", "outros_candidatos_estaveis_no_contexto",
              "ids_vizinhos_no_contexto", "repetivel_selecao_completa",
              "vantagem_prominencia_sobre_segundo", "alertas_opticos", "diagnosticos_comprimento"}
    for sigma in CFG["textura"]["sigmas_estabilidade"]:
        fields.update((f"ancora_reselecionada_sigma{sigma}", f"repetivel_selecao_sigma{sigma}"))
    return fields


def source_provenance(record, pano, rows):
    """Recover the actual source column, including the v11 feather substitution."""
    run, center = record["aquisicao"], record["ancora"]
    start, stop = record["segmento_inicio"], record["segmento_fim"]
    with np.load(BASE / "intermediate" / f"{run}_input.npz") as cached:
        widths = cached["shifts"].astype(int)
        raw_shifts = np.r_[cached["raw_shifts"], cached["raw_shifts"][-1]]
        seams = cached["seams"].astype(int)
    bounds = np.r_[0, np.cumsum(widths)]
    assert bounds[-1] == pano.shape[1] and len(raw_shifts) == len(widths)
    np.testing.assert_array_equal(seams, bounds[1:-1])
    assert 0 <= start <= center < stop <= pano.shape[1]

    def frame_at(column):
        return int(np.searchsorted(bounds, column, side="right") - 1)

    nominal_frame = frame_at(center)
    offset = int(center - bounds[nominal_frame])
    # With feather=2, [0, 1] replaces only seam-1 with column 1 of the next strip.
    modified = seams[np.minimum(widths[:-1], widths[1:]) >= 4] - 1
    replaced = bool(np.any(modified == center))
    frame = nominal_frame + 1 if replaced else nominal_frame
    if replaced:
        offset = 1
    source = DATA / run / f"ImgA{frame:06d}.jpg"
    with Image.open(source) as image:
        original = np.asarray(image.convert("L"))
    column = int(original.shape[1] // 2 - widths[frame] // 2 + offset)
    if not np.array_equal(pano[:, center], original[:, column]):
        raise ValueError(f"Coluna da ancora nao confere com o quadro de origem: {run}, {center}")

    first, last = frame_at(start), frame_at(stop - 1)
    neighbours = []
    for row in rows:
        if row["anchor_id"] == record["anchor_id"] or not row.get("elegivel_individual", row["elegivel"]):
            continue
        inside = (row["inicio"] < stop and row["fim"] > start) if record["tipo"] == "evento_delimitado" else (
            start <= row["ancora"] < stop)
        if inside:
            neighbours.append(row["anchor_id"])
    return dict(
        arquivo_quadro_origem=str(source.relative_to(ROOT.parent)),
        quadro_ancora=source.name, coluna_quadro_ancora=column,
        quadro_origem_sha256=sha(source), coluna_origem_conferida=True,
        quadro_faixa_nominal=f"ImgA{nominal_frame:06d}.jpg",
        ancora_na_coluna_substituida_feather=replaced,
        distancia_costura_ancora=int(np.min(abs(seams-center))),
        costuras_contexto=int(((seams > start) & (seams < stop)).sum()),
        colunas_substituidas_feather_contexto=int(((modified >= start) & (modified < stop)).sum()),
        deslocamentos_nao_positivos_contexto=int((raw_shifts[first:last+1] <= 0).sum()),
        outros_candidatos_estaveis_no_contexto=len(neighbours),
        ids_vizinhos_no_contexto=";".join(neighbours))


def complete_selection_repeatability(record, raw, rows):
    """Rerun detection, stability filters and ranking with each reference sigma.

    Keep the final window fixed to isolate detector/selection sensitivity.
    This diagnostic does not replace the winner at the configured sigma=16.
    """
    from select_anchors import texture_candidates
    result = {}
    length, center = record["comprimento"], record["ancora"]
    for sigma in CFG["textura"]["sigmas_estabilidade"]:
        candidates, _ = texture_candidates(record["aquisicao"], record["padrao"], raw,
                                            reference_sigma=sigma)
        eligible = [r for r in candidates if r["elegivel"] and r["comprimento_minimo"] <= length
                    and context_ok(r["ancora"], len(raw), length)]
        best = min(eligible, key=lambda r: (-r["prominencia"], r["variacao_maxima"], r["ancora"])) if eligible else None
        result[f"ancora_reselecionada_sigma{sigma}"] = best["ancora"] if best else None
        result[f"repetivel_selecao_sigma{sigma}"] = bool(best and abs(best["ancora"]-center) <= CFG["textura"]["tolerancia"])
    result["repetivel_selecao_completa"] = all(result[f"repetivel_selecao_sigma{s}"]
                                             for s in CFG["textura"]["sigmas_estabilidade"])
    competitors = [r for r in rows if r["elegivel"] and r["anchor_id"] != record["anchor_id"]
                   and r["comprimento_minimo"] <= length and context_ok(r["ancora"], len(raw), length)]
    result["vantagem_prominencia_sobre_segundo"] = record["prominencia"] - max(r["prominencia"] for r in competitors) if competitors else None
    return result


def assess_anchor(record, pano, raw, rows):
    result = source_provenance(record, pano, rows)
    alerts = []
    if result["ancora_na_coluna_substituida_feather"]:
        alerts.append("coluna_da_ancora_substituida_pela_costura_v11")
    if result["deslocamentos_nao_positivos_contexto"]:
        alerts.append("deslocamentos_nao_positivos_no_contexto")
    if record["tipo"] == "vale_optico":
        result.update(complete_selection_repeatability(record, raw, rows))
        if not result["repetivel_selecao_completa"]:
            alerts.append("escolha_sensivel_ao_sigma")
    result["alertas_opticos"] = ";".join(alerts)
    result["diagnosticos_comprimento"] = record["segmento_fim"] - record["segmento_inicio"]
    return result
