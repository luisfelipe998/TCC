#!/usr/bin/env python3
"""Screen real tracks for conditional liquid-piston/bubble velocity intervals.

No independent liquid measurement or complete-cell mean is produced. Numerical
endpoints are conditional on the physical hypotheses in the accompanying report.
"""
from pathlib import Path
import csv
from collections import Counter
import gzip
import hashlib
import json
import os
import sys
import numpy as np
from scipy.spatial.distance import cdist

BASE = Path(__file__).resolve().parents[1]
SRC2 = BASE.parent
OUT = BASE / 'results'
THRESHOLDS = ('0.6', '0.7', '0.8')
MASK_DIR = SRC2 / 'v15_gas_liquid_velocity/phase_conditioned_motion/intermediate/masks'
sys.path.insert(0, str(BASE / 'scripts'))
from physical_bounds import slug_ratio_bounds


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_gzip(path):
    with gzip.open(path, 'rt') as f:
        return json.load(f)


def save_csv(name, rows):
    if not rows:
        return
    with (OUT / name).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def same_objects(reference, alternate):
    """Unique reciprocal candidate under class, centroid and area gates."""
    if not reference or not alternate:
        return {}
    xy_r = np.array([[o['x'], o['y']] for o in reference])
    xy_a = np.array([[o['x'], o['y']] for o in alternate])
    distance = cdist(xy_r, xy_a)
    area_ratio = np.array([o['area'] for o in alternate])[None, :] / np.array([o['area'] for o in reference])[:, None]
    large_r = np.array([o['object_class'] == 'large' for o in reference])
    large_a = np.array([o['object_class'] == 'large' for o in alternate])
    eligible = (distance <= 10) & (area_ratio >= .5) & (area_ratio <= 2) & large_r[:, None] & large_a[None, :]
    pairs = {}
    for i in np.flatnonzero(eligible.sum(axis=1) == 1):
        j = int(np.flatnonzero(eligible[i])[0])
        if eligible[:, j].sum() == 1:
            pairs[int(i)] = j
    return pairs


def clear_span(gas_union, node, side):
    """Longest fully mask-clear axial interval in a prespecified side window."""
    width = gas_union.shape[1]
    if side == 'left':
        lo, hi = max(0, node['xmin'] - 100), max(0, node['xmin'] - 5)
    elif side == 'right':
        lo, hi = min(width, node['xmax'] + 5), min(width, node['xmax'] + 100)
    else:
        raise ValueError('invalid side')
    if hi <= lo:
        return 0, None, None
    clear = ~np.any(gas_union[18:71, lo:hi], axis=0)
    edges = np.diff(np.r_[False, clear, False].astype(np.int8))
    starts = np.flatnonzero(edges == 1)
    stops = np.flatnonzero(edges == -1)
    if len(starts) == 0:
        return 0, None, None
    k = int(np.argmax(stops - starts))
    return int(stops[k] - starts[k]), int(lo + starts[k]), int(lo + stops[k])


def controls():
    a = [dict(x=100., y=40., area=2000, object_class='large')]
    b = [dict(x=105., y=40., area=2100, object_class='large')]
    assert same_objects(a, b) == {0: 0}
    assert not same_objects(a, b + [dict(x=107., y=40., area=2100, object_class='large')])
    assert not same_objects(a + [dict(x=103., y=40., area=2000, object_class='large')], b)
    assert not same_objects(a, [dict(x=105., y=40., area=6000, object_class='large')])
    assert not same_objects(a, [dict(x=105., y=40., area=2100, object_class='small')])
    m = np.ones((80, 1016), bool)
    node = dict(xmin=200, xmax=400)
    m[18:71, 160:170] = False
    assert clear_span(m, node, 'left') == (10, 160, 170)
    m[30, 165] = True
    assert clear_span(m, node, 'left')[0] == 5
    assert clear_span(m, node, 'right')[0] == 0
    m[18:71, 450:461] = False
    assert clear_span(m, node, 'right') == (11, 450, 461)
    assert slug_ratio_bounds(.75) == (.75, 1.)
    return dict(unique_match=True, ambiguous_matches_rejected=True, area_and_class_checked=True,
                exact_clear_span=True, one_gas_pixel_breaks_clear_column=True,
                meaning='Algorithm controls only; no experimental phase validation.')


def finite_or_blank(value):
    return value if value is not None and np.isfinite(value) else ''


def analyze_run(inv, flat_index):
    run = inv['run']
    trackpath = BASE / 'intermediate' / f'tracks_{run}_0.7.json.gz'
    objectpath = BASE / 'intermediate' / f'objects_{run}.json.gz'
    tracks = read_gzip(trackpath)
    sequences = read_gzip(objectpath)['sequences']
    pattern_ok = inv['pattern'] in ('plug', 'slug')
    maskpath = MASK_DIR / (run + '.npz')
    masks = None
    if pattern_ok:
        with np.load(maskpath, allow_pickle=False) as z:
            packed_union = z['masks_060'] | z['masks_070'] | z['masks_080']
        masks = np.unpackbits(packed_union, axis=-1, count=1016, bitorder='big').astype(bool)
        assert masks.shape == (501, 80, 1016)
    mappings = {}
    rows, node_rows, examples = [], [], []
    for track in tracks:
        nodes = track['nodes']
        tid = int(track['track_id'])
        flat = flat_index[run, tid]
        assert (flat['accepted'] == 'True') == track['accepted']
        assert abs(float(flat['signed_px_frame']) - track['signed_px_frame']) < 1e-9
        assert len(nodes) == track['intervals'] + 1
        class_ok = track['object_class'] == 'large'
        rigidity_ok = bool(track['strict_deformation'])
        volume_ok = track['volume_max_min_ratio'] <= 1.10
        v = track['signed_px_frame']
        centroid_ok = track['volume_area_speed_difference_px_frame'] <= max(1., .1 * abs(v))
        direction_ok = abs(v) > 1e-12 and v * track['volume_centroid_px_frame'] > 0
        matched_nodes = []
        left_count = right_count = 0
        local_node_rows = []
        if pattern_ok and class_ok:
            for node in nodes:
                frame, index = node['frame'], node['index']
                source = sequences['0.7'][frame][index]
                assert abs(source['x'] - node['x']) < 1e-9
                for threshold in ('0.6', '0.8'):
                    key = threshold, frame
                    if key not in mappings:
                        mappings[key] = same_objects(sequences['0.7'][frame], sequences[threshold][frame])
                idx06 = mappings['0.6', frame].get(index)
                idx08 = mappings['0.8', frame].get(index)
                jointly_matched = idx06 is not None and idx08 is not None
                alpha06 = sequences['0.6'][frame][idx06]['central_alpha_conditional'] if idx06 is not None else None
                alpha07 = node['central_alpha_conditional']
                alpha08 = sequences['0.8'][frame][idx08]['central_alpha_conditional'] if idx08 is not None else None
                if jointly_matched:
                    matched_nodes.append([alpha06, alpha07, alpha08])
                left = clear_span(masks[frame], node, 'left')
                right = clear_span(masks[frame], node, 'right')
                left_count += left[0] >= 10
                right_count += right[0] >= 10
                nr = dict(run=run, track_id=tid, frame=frame, reference_object=index,
                          object_060='' if idx06 is None else idx06, object_080='' if idx08 is None else idx08,
                          jointly_matched=jointly_matched, alpha_060=finite_or_blank(alpha06),
                          alpha_070=alpha07, alpha_080=finite_or_blank(alpha08),
                          left_clear_columns=left[0], left_clear_start=left[1], left_clear_stop=left[2],
                          right_clear_columns=right[0], right_clear_start=right[1], right_clear_stop=right[2])
                node_rows.append(nr)
                local_node_rows.append(nr)
        fraction = len(matched_nodes) / len(nodes)
        left_fraction, right_fraction = left_count / len(nodes), right_count / len(nodes)
        geometry_ok = fraction >= .8
        context_ok = max(left_fraction, right_fraction) >= .8
        selected_side = 'left' if left_fraction >= right_fraction else 'right'
        reasons = []
        for ok, reason in ((pattern_ok, 'pattern_outside_plug_slug'), (class_ok, 'not_large'),
                           (track['accepted'], 'base_track_rejected'), (rigidity_ok, 'strict_deformation_failed'),
                           (volume_ok, 'volume_variation_exceeds_10pct'), (centroid_ok, 'centroid_geometry_disagreement'),
                           (direction_ok, 'undefined_or_inconsistent_direction'),
                           (geometry_ok, 'cross_threshold_identity_below_80pct'),
                           (context_ok, 'same_side_mask_clear_context_below_80pct')):
            if not ok:
                reasons.append(reason)
        passed = not reasons
        alpha = np.array(matched_nodes) if matched_nodes else None
        alpha_fields = {}
        for j, threshold in enumerate(('060', '070', '080')):
            alpha_fields['alpha_' + threshold + '_min'] = float(alpha[:, j].min()) if alpha is not None else ''
            alpha_fields['alpha_' + threshold + '_median'] = float(np.median(alpha[:, j])) if alpha is not None else ''
            alpha_fields['alpha_' + threshold + '_max'] = float(alpha[:, j].max()) if alpha is not None else ''
        alpha_lo = float(alpha.min()) if alpha is not None else None
        alpha_hi = float(alpha.max()) if alpha is not None else None
        c_signed = track['volume_centroid_px_frame'] * float(inv['fps'])
        ratio_lo, ratio_hi = slug_ratio_bounds(alpha_lo) if passed else (None, None)
        row = dict(run=run, pattern=inv['pattern'], track_id=tid, object_class=track['object_class'],
                   first_frame=track['first_frame'], last_frame=track['last_frame'], intervals=track['intervals'],
                   base_accepted=track['accepted'], base_rejection=track['rejection'],
                   pattern_eligible=pattern_ok, class_eligible=class_ok, strict_deformation=rigidity_ok,
                   volume_max_min_ratio=track['volume_max_min_ratio'], volume_gate=volume_ok,
                   centroid_speed_difference_px_frame=track['volume_area_speed_difference_px_frame'],
                   centroid_gate=centroid_ok, matched_nodes=len(matched_nodes), track_nodes=len(nodes),
                   direction_gate=direction_ok,
                   cross_threshold_fraction=fraction, cross_threshold_gate=geometry_ok,
                   left_context_fraction=left_fraction, right_context_fraction=right_fraction,
                   chosen_context_side=selected_side if pattern_ok and class_ok else '', context_gate=context_ok,
                   **alpha_fields, alpha_sensitivity_min=finite_or_blank(alpha_lo),
                   alpha_sensitivity_max=finite_or_blank(alpha_hi),
                   candidate_bubble_signed_px_s=c_signed, candidate_bubble_speed_px_s=abs(c_signed),
                   conditional_piston_ratio_lower=finite_or_blank(ratio_lo),
                   conditional_piston_ratio_upper=finite_or_blank(ratio_hi),
                   conditional_piston_speed_lower_px_s=abs(c_signed)*ratio_lo if passed else '',
                   conditional_piston_speed_upper_px_s=abs(c_signed)*ratio_hi if passed else '',
                   conditional_candidate=passed, rejection=';'.join(reasons),
                   interpretation='Conditional liquid-piston interval; not measured material liquid velocity.')
        rows.append(row)
        if passed:
            eligible_examples = [nr for nr in local_node_rows if nr['jointly_matched'] and nr[selected_side + '_clear_columns'] >= 10]
            assert eligible_examples
            example = eligible_examples[len(eligible_examples)//2]
            examples.append(dict(run=run, track_id=tid, node=next(n for n in nodes if n['frame'] == example['frame']),
                                 context=example, side=selected_side, ratio_lower=ratio_lo,
                                 gas_union=masks[example['frame']].copy()))
    considered = [r for r in rows if r['pattern_eligible'] and r['class_eligible']]
    passed_rows = [r for r in rows if r['conditional_candidate']]
    summary = dict(run=run, pattern=inv['pattern'], pattern_eligible=pattern_ok,
                   reference_tracks=len(rows), large_plug_slug_tracks=len(considered),
                   accepted_large_tracks=sum(r['base_accepted'] for r in considered),
                   rigid_volume_geometry_candidates=sum(r['base_accepted'] and r['strict_deformation'] and r['volume_gate'] and r['centroid_gate'] for r in considered),
                   conditional_candidates=len(passed_rows),
                   unique_candidate_intervals=len({i for r in passed_rows for i in range(r['first_frame'],r['last_frame'])}),
                   ratio_lower_min=min((r['conditional_piston_ratio_lower'] for r in passed_rows),default=''),
                   ratio_lower_max=max((r['conditional_piston_ratio_lower'] for r in passed_rows),default=''),
                   status='conditional_candidates' if passed_rows else 'no_eligible_candidate' if pattern_ok else 'model_not_applicable')
    provenance = dict(run=run, tracks_sha256=sha(trackpath), object_catalog_sha256=sha(objectpath),
                      original_phase_masks_sha256=sha(maskpath) if pattern_ok else None)
    return rows, node_rows, summary, examples, provenance


def figures(examples, passed_rows):
    if not examples:
        return []
    os.environ.setdefault('MPLCONFIGDIR', str(BASE/'intermediate'/'conditional_mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from PIL import Image
    inputs = []
    # Select the first passing trajectory from each run, up to six examples.
    selected = []
    seen = set()
    for ex in examples:
        if ex['run'] not in seen:
            selected.append(ex)
            seen.add(ex['run'])
    selected = selected[:6]
    fig, axes = plt.subplots(len(selected), 1, figsize=(13, 2.0*len(selected)), squeeze=False)
    for ax, ex in zip(axes[:, 0], selected):
        node, context = ex['node'], ex['context']
        frame = node['frame']
        path = SRC2/'data_processed'/ex['run']/f'ImgA{frame:06d}.jpg'
        im = np.asarray(Image.open(path).convert('L'))
        ax.imshow(im, cmap='gray', vmin=0, vmax=255)
        ax.contour(ex['gas_union'].astype(float), levels=[.5], colors=['#0072B2'], linewidths=.45)
        ax.add_patch(Rectangle((node['xmin']-.5,node['ymin']-.5),node['width'],node['height'],
                               fill=False,edgecolor='#D55E00',linewidth=1.4))
        start = context[ex['side']+'_clear_start']
        stop = context[ex['side']+'_clear_stop']
        ax.add_patch(Rectangle((start-.5,17.5),stop-start,53,fill=False,edgecolor='#009E73',linewidth=1.5))
        ax.set_title(f"{ex['run']} | trajetória {ex['track_id']} | quadro {frame} | Uₛ/c ∈ [{ex['ratio_lower']:.3f}; 1], condicionado",fontsize=10)
        ax.set_xlim(max(0,min(node['xmin'],start)-20),min(1016,max(node['xmax'],stop)+20))
        ax.set_ylim(79,0)
        ax.set_ylabel('Linha (px)')
        inputs.append(dict(path=str(path.relative_to(SRC2)),sha256=sha(path),track_id=ex['track_id'],frame=frame))
    axes[-1,0].set_xlabel('Coluna (px)')
    fig.tight_layout()
    fig.savefig(BASE/'figures/conditional_piston_native_examples_PT_BR.png',dpi=300)
    fig.savefig(BASE/'figures/conditional_piston_native_examples_PT_BR.pdf')
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7.3,5.6))
    runs = sorted({r['run'] for r in passed_rows})
    colors = plt.get_cmap('tab10')
    for i, run in enumerate(runs):
        rs = [r for r in passed_rows if r['run']==run]
        c = [r['candidate_bubble_speed_px_s'] for r in rs]
        lo = [r['conditional_piston_speed_lower_px_s'] for r in rs]
        hi = [r['conditional_piston_speed_upper_px_s'] for r in rs]
        ax.vlines(c,lo,hi,color=colors(i%10),alpha=.75,linewidth=1.4,label=run)
        ax.scatter(c,lo,marker='_',s=35,color=colors(i%10))
        ax.scatter(c,hi,marker='_',s=35,color=colors(i%10))
    ax.set_xlabel('Translação da bolha candidata |c| (px/s)')
    ax.set_ylabel('Intervalo condicional do pistão líquido |Uₛ| (px/s)')
    ax.ticklabel_format(style='sci',axis='both',scilimits=(0,0),useMathText=True)
    ax.grid(alpha=.2)
    ax.legend(fontsize=8,frameon=False)
    fig.tight_layout()
    fig.savefig(BASE/'figures/conditional_piston_intervals_PT_BR.png',dpi=300)
    fig.savefig(BASE/'figures/conditional_piston_intervals_PT_BR.pdf')
    plt.close(fig)
    return inputs


def report(rows, summaries):
    considered = [r for r in rows if r['pattern_eligible'] and r['class_eligible']]
    passing = [r for r in rows if r['conditional_candidate']]
    run_count = len({r['run'] for r in passing})
    text = f'''# Triagem experimental de velocidades condicionais do pistão líquido

**Resultado: {len(passing)} trajetórias candidatas em {run_count} gravações.** Foram examinadas {len(considered)} trajetórias de objetos grandes nos seis vídeos plug/slug, mantendo na tabela todas as {len(rows)} trajetórias de referência dos 26 vídeos. Os critérios foram registrados no [protocolo](conditional_piston_protocol_PT_BR.md) antes da execução desta etapa.

Os valores publicados são intervalos `U_s/c`, condicionados à reconstrução de área transversal, à interpretação de uma bolha isolada quase rígida com pistão líquido adjacente, à ausência relevante de mudança de fase e à hipótese **0 ≤ velocidade média do filme ≤ velocidade da bolha**. Não são duas velocidades materiais medidas independentemente. `U_s` corresponde somente ao líquido no pistão; a média líquida de uma célula completa ou da gravação não foi calculada.

## Resultados por gravação plug/slug

| Gravação | Padrão | Trajetórias grandes | Aceitas originalmente | Rígidas/volume/centroides | Candidatas condicionais | Intervalos distintos observados |
|---|---|---:|---:|---:|---:|---:|
'''
    for s in summaries:
        if s['pattern_eligible']:
            text += f"| {s['run']} | {s['pattern']} | {s['large_plug_slug_tracks']} | {s['accepted_large_tracks']} | {s['rigid_volume_geometry_candidates']} | {s['conditional_candidates']} | {s['unique_candidate_intervals']} |\n"
    gates = [('Aceitação original da trajetória','base_accepted'),
             ('Deformação estrita','strict_deformation'),
             ('Variação do volume ≤10%','volume_gate'),
             ('Concordância entre centroides','centroid_gate'),
             ('Correspondência simultânea entre limiares ≥80%','cross_threshold_gate'),
             ('Contexto sem máscara na mesma lateral ≥80%','context_gate')]
    text += '\n## Dificuldades observadas\n\n| Critério isolado | Trajetórias grandes que passaram |\n|---|---:|\n'
    for name,key in gates:
        text += f'| {name} | {sum(r[key] for r in considered)}/{len(considered)} |\n'
    failures = Counter(reason for r in considered for reason in r['base_rejection'].split(';') if reason)
    text += '\nAs contagens são marginais; falhas podem se sobrepor. No rastreador original, as rejeições dos objetos grandes foram: ' + ', '.join(f'`{key}`: {value}' for key,value in sorted(failures.items())) + '.\n'
    if considered:
        text += f"\nA maior fração de nós com a mesma identificação nos três limiares foi {max(r['cross_threshold_fraction'] for r in considered):.1%}, abaixo dos 80% exigidos. O simples achado de um trecho adjacente sem máscara não superou a instabilidade das trajetórias e da identidade entre segmentações. Não se enfraqueceram os filtros para converter esse contexto em velocidade líquida.\n"
    text += '''
## Como interpretar o intervalo

A conservação no referencial da bolha fornece `U_s/c = α_b + (1−α_b) U_f/c`. O valor de `α_b` foi aproximado pelo indicador de área transversal circular do corpo da bolha. Em cada nó correspondido, compararam-se os mesmos objetos nos três limiares. O extremo inferior apresentado é o menor indicador entre esses limiares e nós; o superior é 1, imposto pela hipótese sobre o filme. O extremo inferior combina variação temporal e sensibilidade de segmentação: **não é um limite garantido da geometria real nem um intervalo de confiança**.

A triagem de contexto procurou uma mesma lateral da trajetória com ao menos dez colunas contíguas sem gás em nenhuma das três máscaras, em pelo menos 80% dos seus nós. Isso exige ausência de marcação em toda a altura do lúmen projetado. Ainda assim, máscara vazia pode significar falha de detecção e não exclui gás fora do plano de observação. A triagem não comprova a velocidade do filme.

As trajetórias foram condicionadas ainda à estabilidade de comprimento e de volume reconstruído, concordância dos dois centroides e correspondência única entre segmentações. Os controles de correspondência e contexto passaram, mas verificam apenas o algoritmo. A escolha dos critérios seleciona objetos fáceis de seguir e não garante representatividade do conjunto gasoso. As trajetórias da mesma gravação não são réplicas experimentais independentes.

Os valores em pixels/s usam o módulo da velocidade do centroide volumétrico condicional e o FPS nativo. A razão dentro de uma trajetória cancela uma escala espacial comum. Comparações de velocidade física entre gravações ainda exigem escala axial relativa ou hipótese explícita de magnificação igual.
'''
    if passing:
        text += '''
## Figuras

![Intervalos condicionais do pistão líquido](../figures/conditional_piston_intervals_PT_BR.png)

Cada segmento vertical corresponde a uma trajetória candidata; suas pontas indicam o envelope condicional, sem estimativa pontual líquida ou barra de confiança. O eixo horizontal é a translação de uma bolha candidata. Não se deve reutilizar este gráfico como mapa de velocidades superficiais de gás e líquido.

![Contexto dos candidatos em quadros nativos](../figures/conditional_piston_native_examples_PT_BR.png)

Azul: contorno da união das três máscaras originais; laranja: caixa do objeto de referência; verde: trecho adjacente sem pixels marcados como gás. Cada exemplo é o nó central entre os que satisfazem simultaneamente correspondência e contexto, da primeira trajetória aceita de cada gravação. A imagem não valida gás/liquido fora do plano.
'''
    else:
        text += '''
Nenhuma trajetória satisfez simultaneamente os critérios. Não se produziram intervalos ou gráficos com valores preenchidos por hipótese. A rejeição desta triagem não demonstra inexistência de pistões físicos; indica que a evidência e os critérios aplicados não sustentaram esta redução nas trajetórias disponíveis.
'''
    text += '''
## Arquivos e reprodução

- [Todas as trajetórias, critérios e falhas](../results/conditional_piston_tracks.csv).
- [Contexto e correspondências por nó](../results/conditional_piston_nodes.csv).
- [Sumário dos 26 vídeos](../results/conditional_piston_summary.csv).
- [Controles, critérios e proveniência](../results/conditional_piston_validation.json).
- [Derivação física e limites de interpretação](physical_identifiability_review_PT_BR.md).

Executar `python3 -B src2/v15.2_gas_liquid_relative_velocity/scripts/conditional_bubble_bounds.py` após concluir `object_motion.py`. A etapa lê os catálogos de objetos, as máscaras originais v15 e os quadros nativos usados nos exemplos. Não altera imagens nem resultados anteriores.
'''
    (BASE/'findings/conditional_piston_PT_BR.md').write_text(text)


def run():
    checks = controls()
    validation_path = OUT/'object_validation.json'
    if not validation_path.exists():
        raise RuntimeError('object_motion.py has not completed all 26 recordings')
    validation = json.loads(validation_path.read_text())
    assert validation['recordings'] == 26 and validation['frames_per_run'] == 501
    inventory_path = SRC2/'v15.1_gas_liquid_velocity/results/scenario_inventory.csv'
    inv = [r for r in csv.DictReader(inventory_path.open()) if r['included']=='True']
    assert len(inv)==26
    flat_path = OUT/'object_tracks.csv'
    flat = {(r['run'],int(r['track_id'])):r for r in csv.DictReader(flat_path.open()) if float(r['threshold'])==.7}
    rows, nodes, summaries, examples, provenance = [], [], [], [], []
    for item in inv:
        r,n,s,e,p = analyze_run(item,flat)
        rows.extend(r);nodes.extend(n);summaries.append(s);examples.extend(e);provenance.append(p)
        print(item['run'],s['large_plug_slug_tracks'],s['conditional_candidates'],flush=True)
    assert len(rows)==len(flat)
    save_csv('conditional_piston_tracks.csv',rows)
    save_csv('conditional_piston_nodes.csv',nodes)
    save_csv('conditional_piston_summary.csv',summaries)
    passing=[r for r in rows if r['conditional_candidate']]
    image_inputs=figures(examples,passing)
    record=dict(controls=checks,recordings=26,reference_tracks=len(rows),conditional_tracks=len(passing),
                conditional_runs=len({r['run'] for r in passing}),
                criteria=dict(patterns=['plug','slug'],class_required='large',reference_threshold=.7,
                    base_acceptance=True,strict_deformation=True,volume_max_min_limit=1.10,
                    centroid_difference_limit='max(1 px/frame, 10% abs(projected speed))',
                    nonzero_consistent_centroid_direction=True,
                    cross_threshold_distance_px=10,cross_threshold_area_ratio=[.5,2],
                    cross_threshold_unique_reciprocal=True,joint_matched_nodes_fraction=.8,
                    clear_columns=10,clear_context_fraction=.8,context_same_side=True,
                    context_window_left='[xmin-100:xmin-5]',context_window_right='[xmax+5:xmax+100]',
                    mask_union=['masks_060','masks_070','masks_080'],lumen=[18,71],
                    interval_lower='min central_alpha_conditional over jointly matched nodes and three thresholds',
                    interval_upper=1),
                script_sha256=sha(Path(__file__)),inventory_sha256=sha(inventory_path),
                object_validation_sha256=sha(validation_path),object_tracks_sha256=sha(flat_path),
                protocol_sha256=sha(BASE/'findings/conditional_piston_protocol_PT_BR.md'),
                inputs=provenance,native_examples=image_inputs,
                interpretation='Conditional piston intervals only. No independent material liquid velocity; no whole-cell/video liquid mean.')
    considered=[r for r in rows if r['pattern_eligible'] and r['class_eligible']]
    record['large_plug_slug_tracks']=len(considered)
    record['marginal_gate_pass_counts']={key:sum(r[key] for r in considered) for key in
        ('base_accepted','strict_deformation','volume_gate','centroid_gate','direction_gate','cross_threshold_gate','context_gate')}
    record['base_rejection_counts']=dict(Counter(reason for r in considered for reason in r['base_rejection'].split(';') if reason))
    record['maximum_cross_threshold_fraction']=max((r['cross_threshold_fraction'] for r in considered),default=None)
    (OUT/'conditional_piston_validation.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    report(rows,summaries)
    print(json.dumps({k:record[k] for k in ('reference_tracks','conditional_tracks','conditional_runs')},indent=2))


if __name__=='__main__':
    run()
