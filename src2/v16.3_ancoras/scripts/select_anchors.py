"""Select from optical evidence only. This module never imports wavelets."""
from datetime import datetime, timezone
import numpy as np
from scipy.ndimage import binary_closing, gaussian_filter1d, median_filter
from scipy.signal import find_peaks

from common import (ANCHORS, BASE, CFG, INTER, OLD, ROOT, SCRIPTS,
                    context_ok, image_path, inventory, read_csv, setup, sha,
                    slug, verify_source, write_csv, write_json, required_window)


def intervals(pano, threshold):
    rule = CFG["eventos"]
    lo, hi = CFG["linhas"]
    fraction = (pano[lo:hi + 1] < threshold).mean(axis=0)
    mask = median_filter(fraction, size=rule["mediana"], mode="nearest") >= rule["fracao_escura"]
    size = rule["fechamento"]
    mask = binary_closing(np.pad(mask, size, mode="edge"),
                          structure=np.ones(size))[size:-size]
    changes = np.diff(np.r_[False, mask, False].astype(int))
    return [(int(a), int(b)) for a, b in zip(np.where(changes == 1)[0],
                                           np.where(changes == -1)[0])
            if b - a >= rule["extensao"][0]]


def match_interval(start, stop, alternatives):
    matches = []
    for left, right in alternatives:
        iou = max(0, min(stop, right) - max(start, left)) / (max(stop, right) - min(start, left))
        if iou >= CFG["eventos"]["iou_correspondencia"]:
            matches.append((left, right))
    # A split/merge with multiple matching envelopes is not silently resolved.
    return matches[0] if len(matches) == 1 else None


def endpoint_contrast(raw, left, right, width=None):
    """Median outside-minus-inside intensity; weakest endpoint breaks ties only."""
    width = CFG["eventos"]["nitidez_faixa_px"] if width is None else int(width)
    if width <= 0:
        raise ValueError("A faixa de contraste deve ser positiva.")
    if left < width or right+width > len(raw) or right-left < 2*width:
        return dict(contraste_esquerdo=None, contraste_direito=None, contraste_minimo=None)
    left_contrast = float(np.median(raw[left-width:left]) - np.median(raw[left:left+width]))
    right_contrast = float(np.median(raw[right:right+width]) - np.median(raw[right-width:right]))
    return dict(contraste_esquerdo=left_contrast, contraste_direito=right_contrast,
                contraste_minimo=min(left_contrast, right_contrast))


def event_rank(record):
    return (float(record["variacao_maxima"]), -float(record["contraste_minimo"]), int(record["ancora"]))


def event_candidates(run, pattern, pano, raw, reviews, previous):
    rule = CFG["eventos"]
    alternatives = {q: intervals(pano, q) for q in rule["limiares_estabilidade"]}
    rows = []
    for index, (left, right) in enumerate(intervals(pano, rule["limiar"]), 1):
        ident = f"{run}_E{index:03d}"
        center = (left + right) // 2
        window = required_window(center, left, right, "evento_delimitado")
        review = reviews.get(ident)
        old = previous.get(ident)
        if old is not None:
            assert (left, right, center) == tuple(int(old[k]) for k in ("inicio", "fim", "ancora")), ident
        reasons = []
        if review is None or old is None:
            reasons.append("inspecao_visual_ausente")
        elif review["aceito"] != "True":
            reasons.append(review["motivo"])
        if left == 0 or right == len(raw):
            reasons.append("evento_toca_borda")
        if right - left > rule["extensao"][1]:
            reasons.append("extensao_maior_1536")
        if not context_ok(center, len(raw), window):
            reasons.append("contexto_do_evento_mais_16_indisponivel")
        if context_ok(center, len(raw), window) and raw[center - window // 2:center + window // 2].std() < 1e-10:
            reasons.append("sinal_degenerado")
        record = dict(aquisicao=run, padrao=pattern, anchor_id=ident,
                      tipo="evento_delimitado", ancora=center, inicio=left, fim=right,
                      largura=right-left, comprimento_minimo=window, comprimento=window,
                      tamanho_panorama=len(raw), segmento_inicio=center-window//2, segmento_fim=center+window//2,
                      revisao_visual_aceita=bool(review and review["aceito"] == "True"),
                      motivo_visual=review["motivo"] if review else "pendente",
                      aprovado_regras_base=not reasons)
        shifts = []
        support_left, support_right = left, right
        for threshold, options in alternatives.items():
            other = match_interval(left, right, options)
            record[f"inicio_limiar{threshold}"] = other[0] if other else None
            record[f"fim_limiar{threshold}"] = other[1] if other else None
            change = max(abs(left-other[0]), abs(right-other[1])) if other else None
            record[f"deslocamento_limiar{threshold}"] = change
            if change is None:
                reasons.append(f"correspondencia_ausente_ou_ambigua_limiar{threshold}")
            else:
                shifts.append(change)
                support_left = min(support_left, other[0])
                support_right = max(support_right, other[1])
        # All alternative boundaries must fit around the final (threshold-60) center.
        window = required_window(center, support_left, support_right, "evento_delimitado")
        if not context_ok(center, len(raw), window):
            reasons.append("contexto_dos_tres_limiares_indisponivel")
        record.update(inicio_suporte=support_left, fim_suporte=support_right,
                      comprimento_minimo=window, comprimento=window,
                      segmento_inicio=center-window//2, segmento_fim=center+window//2)
        record.update(variacao_maxima=max(shifts) if len(shifts) == 2 else None,
                      elegivel=not reasons, selecionada=False, ranking=None,
                      motivo=";".join(reasons) if reasons else "elegivel")
        record.update(endpoint_contrast(raw, left, right))
        rows.append(record)
    ranked = sorted((r for r in rows if r["elegivel"]), key=event_rank)
    for rank, row in enumerate(ranked, 1):
        row["ranking"] = rank
    if ranked:
        ranked[0].update(selecionada=True, motivo="estabilidade_com_desempate_por_contraste")
    return rows, ranked[0] if ranked else None


def minima(raw, sigma):
    rule = CFG["textura"]
    smoothed = gaussian_filter1d(np.asarray(raw, dtype=float), sigma, mode="reflect")
    amplitude = float(np.diff(np.percentile(smoothed, [5, 95]))[0])
    threshold = max(rule["prominencia_minima"], rule["fracao_amplitude_robusta"] * amplitude)
    positions, props = find_peaks(-smoothed, prominence=threshold, width=rule["largura"],
                                 distance=rule["distancia"], wlen=rule["janela_prominencia"])
    rows = []
    for j, center in enumerate(positions):
        left_base, right_base = int(props["left_bases"][j]), int(props["right_bases"][j])
        window = required_window(int(center), left_base, right_base, "vale_optico")
        rows.append(dict(ancora=int(center), prominencia=float(props["prominences"][j]),
                         largura=float(props["widths"][j]), inicio=float(props["left_ips"][j]),
                         fim=float(props["right_ips"][j]), limiar_prominencia=threshold,
                         base_esquerda=left_base, base_direita=right_base,
                         comprimento_minimo=window,
                         contexto_valido=context_ok(int(center), len(raw), window)))
    return rows


def match_minimum(center, candidates):
    matches = [r for r in candidates if r["contexto_valido"]
               and abs(r["ancora"] - center) <= CFG["textura"]["tolerancia"]]
    return matches[0] if len(matches) == 1 else None


def texture_candidates(run, pattern, raw, reference_sigma=None):
    rule = CFG["textura"]
    reference_sigma = rule["sigma"] if reference_sigma is None else reference_sigma
    sigmas = {rule["sigma"], *rule["sigmas_estabilidade"]}
    if reference_sigma not in sigmas:
        raise ValueError("Sigma de referencia fora das variantes predefinidas.")
    alternatives = {s: minima(raw, s) for s in sorted(sigmas - {reference_sigma})}
    rows = []
    for candidate in minima(raw, reference_sigma):
        row = dict(candidate)
        center = row["ancora"]
        reasons = [] if row["contexto_valido"] else ["contexto_do_vale_mais_16_indisponivel"]
        shifts = []
        support_left, support_right = row["base_esquerda"], row["base_direita"]
        for sigma, options in alternatives.items():
            other = match_minimum(center, options)
            row[f"ancora_sigma{sigma}"] = other["ancora"] if other else None
            shift = other["ancora"] - center if other else None
            row[f"deslocamento_sigma{sigma}"] = shift
            if shift is None:
                reasons.append(f"vale_ausente_ou_ambiguo_sigma{sigma}")
            else:
                shifts.append(abs(shift))
                support_left = min(support_left, other["base_esquerda"])
                support_right = max(support_right, other["base_direita"])
                row[f"base_esquerda_sigma{sigma}"] = other["base_esquerda"]
                row[f"base_direita_sigma{sigma}"] = other["base_direita"]
        # Recenter the union of all bases; individual radii have different centers.
        window = required_window(center, support_left, support_right, "vale_optico")
        if not context_ok(center, len(raw), window):
            reasons.append("contexto_comum_dos_tres_sigmas_indisponivel")
        row.update(aquisicao=run, padrao=pattern, anchor_id=f"{run}_MIN{center:05d}",
                   tipo="vale_optico", sigma_referencia=reference_sigma,
                   inicio_suporte=support_left, fim_suporte=support_right,
                   comprimento_minimo=window, comprimento=window, tamanho_panorama=len(raw),
                   segmento_inicio=center-window//2, segmento_fim=center+window//2,
                   variacao_maxima=max(shifts) if len(shifts) == 2 else None,
                   elegivel=not reasons, selecionada=False, ranking=None,
                   motivo=";".join(reasons) if reasons else "elegivel")
        rows.append(row)
    ranked = sorted((r for r in rows if r["elegivel"]),
                    key=lambda r: (-r["prominencia"], r["variacao_maxima"], r["ancora"]))
    for rank, row in enumerate(ranked, 1):
        row["ranking"] = rank
    if ranked:
        ranked[0].update(selecionada=True, motivo="maior_prominencia_estavel")
    return rows, ranked[0] if ranked else None


def common_window(groups):
    """Retain each run's own ranking while enforcing equal supported contexts."""
    window = max(g["best"]["comprimento_minimo"] for g in groups)
    for group in groups:
        for row in group["rows"]:
            row.setdefault("elegivel_individual", row["elegivel"])
    preferred = max(window, CFG.get("janelas_por_padrao", {}).get(groups[0]["pattern"], window))
    # If a preferred larger context removes every candidate in a run, back off.
    for trial in range(preferred, window - 1, -CFG["contexto"]["passo"]):
        if all(any(r["elegivel_individual"] and r["comprimento_minimo"] <= trial
                   and context_ok(r["ancora"], r["tamanho_panorama"], trial)
                   for r in group["rows"]) for group in groups):
            window = trial
            break
    for _ in range(20):
        winners = []
        for group in groups:
            options = []
            for row in group["rows"]:
                row.update(selecionada=False, ranking=None)
                row["elegivel"] = row["elegivel_individual"] and context_ok(
                    row["ancora"], row["tamanho_panorama"], max(window, row["comprimento_minimo"]))
                if row["elegivel"]:
                    options.append(row)
                elif row["elegivel_individual"]:
                    row["motivo"] = "contexto_comum_do_par_indisponivel"
            key = event_rank if group["pattern"] in ("plug", "slug") else (
                lambda r: (-r["prominencia"], r["variacao_maxima"], r["ancora"]))
            options.sort(key=key)
            for rank, row in enumerate(options, 1):
                row["ranking"] = rank
                row["motivo"] = "elegivel"
            group["best"] = options[0] if options else None
            if group["best"]:
                winners.append(group["best"])
        if len(winners) != len(groups):
            return window
        required = max(r["comprimento_minimo"] for r in winners)
        if required <= window:
            for row in winners:
                center = row["ancora"]
                row.update(selecionada=True, comprimento=window,
                           segmento_inicio=center-window//2, segmento_fim=center+window//2,
                           motivo="melhor_ancora_individual_com_contexto_comum")
            return window
        window = required
    raise RuntimeError("A escolha do contexto nao convergiu.")


def select_pilot():
    setup()
    from plots import anchor_figure
    from anchor_diagnostics import assess_anchor
    review_file = OLD / "results/revisao_visual.csv"
    previous_file = OLD / "results/catalogo_eventos.csv"
    reviews = {r["event_id"]: r for r in read_csv(review_file)}
    previous = {r["event_id"]: r for r in read_csv(previous_file)}
    sources = {str(p.relative_to(ROOT.parent)): sha(p) for p in
               (review_file, previous_file, BASE / "results/input_audit.csv")}
    all_rows, selected, trials, pairs, diagnostics = [], [], [], [], []
    inv = inventory()
    for pattern in CFG["padroes"]:
        accepted = []
        for source in (r for r in inv if r["pattern"] == pattern):
            run = source["run"]
            pano, raw = verify_source(source)
            for path in (image_path(run), BASE / "intermediate" / f"{run}_input.npz"):
                sources[str(path.relative_to(ROOT.parent))] = sha(path)
            if pattern in ("plug", "slug"):
                rows, best = event_candidates(run, pattern, pano, raw, reviews, previous)
            else:
                rows, best = texture_candidates(run, pattern, raw)
            all_rows.extend(rows)
            trials.append(dict(aquisicao=run, padrao=pattern, candidatos=len(rows),
                               elegiveis_individuais=sum(r["elegivel"] for r in rows),
                               elegiveis=sum(r["elegivel"] for r in rows),
                               status="ancora_selecionada" if best else "sem_ancora_valida"))
            if best is None:
                missing = ANCHORS / f"{slug(pattern)}_{run}_SEM_ANCORA_VALIDA.txt"
                missing.write_text(f"{pattern} | {run}\nNenhum candidato satisfez os criterios de forma, "
                                   "estabilidade e contexto aplicaveis ao padrao.\n"
                                   "Motivos por candidato: ../intermediario/catalogo_candidatos.csv\n",
                                   encoding="utf-8")
                print(f"{pattern}: {run} SEM ANCORA VALIDA", flush=True)
                continue
            accepted.append(dict(best=best, rows=rows, pano=pano, raw=raw, pattern=pattern))
            if len(accepted) == 2:
                window = common_window(accepted)
                for group in accepted:
                    trial = next(t for t in trials if t["aquisicao"] == group["rows"][0]["aquisicao"])
                    trial.update(elegiveis=sum(r["elegivel"] for r in group["rows"]),
                                 comprimento_contexto=window)
                    if group["best"] is None:
                        run = group["rows"][0]["aquisicao"]
                        (ANCHORS / f"{slug(pattern)}_{run}_SEM_ANCORA_VALIDA.txt").write_text(
                            f"{pattern} | {run}\nSem ancora com suporte para o contexto comum de {window} px.\n",
                            encoding="utf-8")
                        next(t for t in trials if t["aquisicao"] == run)["status"] = "sem_ancora_valida"
                accepted = [g for g in accepted if g["best"] is not None]
                if len(accepted) == 2:
                    break
        best_rows = [g["best"] for g in accepted]
        for group in accepted:
            best = group["best"]
            diagnostic = assess_anchor(best, group["pano"], group["raw"], group["rows"])
            diagnostics.append(dict(anchor_id=best["anchor_id"], aquisicao=best["aquisicao"],
                                    padrao=best["padrao"], **diagnostic))
            best.update(diagnostic)
            source = ROOT.parent / diagnostic["arquivo_quadro_origem"]
            sources[str(source.relative_to(ROOT.parent))] = diagnostic["quadro_origem_sha256"]
            selected.append(best)
            anchor_figure(best, group["pano"], group["raw"])
            print(f"{pattern}: {best['anchor_id']} x={best['ancora']} "
                  f"variacao={best['variacao_maxima']} px; N={best['comprimento']}", flush=True)
        if len(accepted) == 2:
            pairs.append(dict(padrao=pattern, run_x=best_rows[0]["aquisicao"],
                              run_y=best_rows[1]["aquisicao"], anchor_x=best_rows[0]["anchor_id"],
                              anchor_y=best_rows[1]["anchor_id"], comprimento=window))
    write_csv(INTER / "catalogo_candidatos.csv", all_rows)
    write_csv(INTER / "ancoras_selecionadas.csv", selected)
    write_csv(INTER / "ensaios.csv", trials)
    write_csv(INTER / "diagnosticos_ancoras.csv", diagnostics)
    used = {r["aquisicao"] for r in selected if r["padrao"] in ("plug", "slug")}
    write_csv(INTER / "revisao_visual_herdada.csv",
              [r for r in reviews.values() if r["event_id"].rsplit("_", 1)[0] in used])
    write_json(INTER / "selecao.json", dict(ancoras=selected, pares=pairs, ensaios=trials,
                                          regra_pares=CFG["ordem_aquisicoes"]))
    frozen_paths = [INTER / name for name in ("catalogo_candidatos.csv", "ancoras_selecionadas.csv",
                    "ensaios.csv", "revisao_visual_herdada.csv", "selecao.json", "diagnosticos_ancoras.csv")]
    frozen_paths += [SCRIPTS / name for name in ("config.json", "common.py", "select_anchors.py", "anchor_diagnostics.py")]
    for name in ("comparacao_janelas.json", "decisao_janelas.json"):
        evidence = INTER / "janelas" / name
        if evidence.exists():
            frozen_paths.append(evidence)
    write_json(INTER / "congelamento.json", dict(data_utc=datetime.now(timezone.utc).isoformat(),
        antes_das_wavelets=True, arquivos_sha256={str(p.relative_to(ROOT)): sha(p) for p in frozen_paths},
        fontes_sha256=sources, total_pares=len(pairs)))
    assert len(pairs) <= 4, "O piloto nao pode expandir para o lote."
    return selected, pairs
