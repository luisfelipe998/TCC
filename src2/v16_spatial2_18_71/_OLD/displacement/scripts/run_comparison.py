"""Native-frame displacement comparison; write only inside displacement/."""
from __future__ import annotations
import argparse
import hashlib
import html
import io
import json
import os
from pathlib import Path
import re
import time

import numpy as np
import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT.parent
DATA = BASE.parent / 'data_processed'
os.environ.setdefault('MPLCONFIGDIR', '/tmp/tcc-displacement-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from motion_core import (estimate, image_alignment, longest_segment,
                         reconstruct_segment, sample_columns, ROI, FIT_ROWS,
                         CHECK_ROWS, SEARCH, MIN_CORRELATION, MIN_MARGIN,
                         MIN_TEXTURE, ABS_TOLERANCE, REL_TOLERANCE)

OLD = '#b65c2e'
NEW = '#087f8c'
REJECTED = '#aaaeb3'
MIN_SEGMENT_INTERVALS = 8


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_frames(run):
    files = sorted((DATA / run).glob('ImgA*.jpg'))
    ids = [int(re.search(r'ImgA(\d+)\.jpg$', p.name)[1]) for p in files]
    assert ids == list(range(501)), run
    digest = hashlib.sha256(); arrays = []
    for p in files:
        raw = p.read_bytes()
        digest.update(p.name.encode()); digest.update(hashlib.sha256(raw).digest())
        with Image.open(io.BytesIO(raw)) as im:
            arrays.append(np.asarray(im.convert('L')))
    frames = np.stack(arrays)
    assert frames.shape == (501, 80, 1016), (run, frames.shape)
    return frames, digest.hexdigest()


def finite_mean(x):
    x = np.asarray(x, float)
    return float(np.mean(x[np.isfinite(x)])) if np.isfinite(x).any() else np.nan


def finite_median(x):
    x = np.asarray(x, float)
    return float(np.median(x[np.isfinite(x)])) if np.isfinite(x).any() else np.nan


def save_segment(run, frames, old, rows):
    segment = longest_segment(rows)
    if segment is None or segment[1] - segment[0] < MIN_SEGMENT_INTERVALS:
        return None
    first, last = segment
    shifts = [r['new_dx'] for r in rows[first:last]]
    new, coordinates = reconstruct_segment(frames[first:last + 1], shifts)
    offsets = np.r_[0, np.cumsum(old['shifts'])]
    with Image.open(DATA / run / '__panorama_v11.png') as image:
        before = np.asarray(image)[:, offsets[first]:offsets[last + 1]].copy()
    Image.fromarray(before).save(ROOT / 'panoramas' / f'{run}_before.png')
    Image.fromarray(new).save(ROOT / 'panoramas' / f'{run}_after.png')
    coordinate_path = ROOT / 'intermediate' / f'{run}_segment_coordinates.npz'
    np.savez_compressed(coordinate_path, first_frame=first, last_frame=last,
                        displacement_dx=np.asarray(shifts), **coordinates)
    for label, image in [('before', before), ('after', new)]:
        pd.DataFrame(dict(position_px=np.arange(image.shape[1]),
                          gray_mean_rows_18_71=image[ROI].mean(axis=0))).to_csv(
            ROOT / 'signals' / f'{run}_{label}.csv', index=False)
    return dict(first=first, last=last, intervals=last - first,
                before=before, after=new, old_width=before.shape[1],
                new_width=new.shape[1], total_travel_px=coordinates['total_travel_px'])


def example_index(rows):
    good = [r['frame_i'] for r in rows if r['supported']]
    # Fixed mid-record target, not a search for the largest apparent gain.
    return min(good, key=lambda i: (abs(i - 250), i)) if good else 250


def alignment_picture(frames, background, i, dx):
    a = frames[i, ROI].astype(float) - background[ROI]
    b = frames[i + 1, ROI].astype(float) - background[ROI]
    reference = a[:, SEARCH:-SEARCH]
    aligned = sample_columns(b, np.arange(SEARCH, frames.shape[2] - SEARCH) + dx)
    return np.concatenate([reference, aligned], axis=0)


def plot_run(inv, frames, background, old, rows, segment):
    run, pattern = inv.run, inv.pattern
    df = pd.DataFrame(rows); ok = df.supported.to_numpy()
    old_dx = df.old_dx.to_numpy(); candidate = df.candidate_dx.to_numpy()
    clamped = df.old_clamped.to_numpy()
    fig, axes = plt.subplots(4, 2, figsize=(17, 12), layout='constrained',
                             gridspec_kw={'height_ratios': [2.1, 1.5, 1.2, 1.4]})
    axes[0, 0].plot(df.frame_i, -df.old_raw_s, color=REJECTED, lw=.7, label='Estimativa original com sinal corrigido')
    axes[0, 0].plot(df.frame_i, old_dx, color=OLD, lw=1, label='Deslocamento usado no panorama v11')
    axes[0, 0].scatter(df.frame_i[clamped], old_dx[clamped], color='#b21c35', s=18,
                       label=f'Forçados a −1 px: {int(clamped.sum())}', zorder=4)
    axes[0, 1].scatter(df.frame_i[~ok], candidate[~ok], color=REJECTED, s=5,
                       label='Candidato sem suporte suficiente')
    axes[0, 1].plot(df.frame_i, df.new_dx, color=NEW, lw=.8)
    axes[0, 1].scatter(df.frame_i[ok], candidate[ok], color=NEW, s=6,
                       label=f'Verificados: {int(ok.sum())}/500')
    limits = np.r_[old_dx, -df.old_raw_s, candidate]
    limits = limits[np.isfinite(limits)]
    lo, hi = min(limits.min(), 0), max(limits.max(), 0)
    pad = max(2, .06 * (hi - lo))
    for ax, title in zip(axes[0], ['ANTES — v11', 'AGORA — estimativa com verificação']):
        ax.set(title=title, xlabel='Quadro inicial do par', ylabel='Deslocamento (px/quadro)',
               xlim=(0, 499), ylim=(lo - pad, hi + pad))
        ax.axhline(0, color='#444', lw=.5); ax.grid(alpha=.18); ax.legend(fontsize=7, loc='best')
    i = example_index(rows)
    picture_old = alignment_picture(frames, background, i, old_dx[i])
    finite_candidate = np.isfinite(candidate[i])
    picture_new = alignment_picture(frames, background, i, candidate[i]) if finite_candidate else np.zeros_like(picture_old)
    clim = max(1., float(np.quantile(abs(np.r_[picture_old.ravel(), picture_new.ravel()]), .99)))
    for ax, image, title in zip(axes[1], [picture_old, picture_new], [
        f'Alinhamento anterior: {old_dx[i]:.2f} px; r = {df.old_dynamic_correlation.iloc[i]:.3f}',
        (f'Alinhamento atual: {candidate[i]:.2f} px; r = {df.new_dynamic_correlation.iloc[i]:.3f}'
         if ok[i] else 'SEM ESTIMATIVA VERIFICADA — candidato ilustrativo')]):
        ax.imshow(image, cmap='gray', vmin=-clim, vmax=clim, aspect='auto',
                  extent=[SEARCH, frames.shape[2] - SEARCH, 108, 0], interpolation='nearest')
        ax.axhline(54, color='tab:orange', lw=.8)
        ax.set(title=title, xlabel=f'Coluna nativa de A; quadros {i} → {i + 1}',
               yticks=[27, 81], yticklabels=['A', 'B alinhado'])
    if segment:
        xmax = max(segment['old_width'], segment['new_width'])
        profiles = [segment['before'][ROI].mean(axis=0), segment['after'][ROI].mean(axis=0)]
        for col, (label, image, profile, color) in enumerate(zip(
                ['Recorte do panorama v11 original', 'Reconstrução em posições subpixel'],
                [segment['before'], segment['after']], profiles, [OLD, NEW])):
            axes[2, col].imshow(image, cmap='gray', vmin=0, vmax=255, aspect='auto',
                                extent=[0, image.shape[1], 80, 0], interpolation='nearest')
            for boundary in (18, 72): axes[2, col].axhline(boundary, color='tab:orange', lw=.5)
            axes[2, col].set(title=f'{label} — {image.shape[1]} px', xlim=(0, xmax),
                             xlabel=f'Mesmos quadros {segment["first"]}–{segment["last"]}; origem local (px)',
                             ylabel='Linha')
            axes[3, col].plot(np.arange(len(profile)), profile, color=color, lw=.8)
            axes[3, col].set(xlim=(0, xmax), ylim=(0, 255), ylabel='Intensidade média',
                             xlabel='Posição local no panorama (px)', title='Média das linhas 18–71')
            axes[3, col].grid(alpha=.18)
    else:
        for ax in axes[2:].flat:
            ax.axis('off'); ax.text(.5, .5, 'Sem segmento contínuo verificado\ncom pelo menos 8 intervalos.\nPanorama novo não gerado.',
                                    ha='center', va='center', transform=ax.transAxes)
    scope = '' if inv.included else ' | fora do escopo principal de 26 registros'
    fig.suptitle(f'{run} | {pattern}{scope}\nNegativo = movimento para a esquerda; lacunas permanecem sem estimativa. '
                 'Reconstruções usam o maior segmento verificado.', fontsize=12)
    fig.savefig(ROOT / 'figures' / f'{run}_comparison.png', dpi=140)
    plt.close(fig)
    return i


def overview(summary):
    eligible = summary[summary.included]
    fig, axes = plt.subplots(1, 2, figsize=(16, 11), layout='constrained',
                             gridspec_kw={'width_ratios': [1, 1.5]})
    for pattern, g in eligible.groupby('pattern'):
        axes[0].scatter(g.old_median_on_supported, g.new_supported_median, s=45, label=pattern)
    lim = min(-1, float(np.nanmin(summary[['old_median_on_supported', 'new_supported_median']].values))) - 3
    axes[0].plot([lim, 0], [lim, 0], color='#555', ls='--', lw=1)
    axes[0].set(xlim=(lim, 0), ylim=(lim, 0), aspect='equal',
                xlabel='Antes: mediana v11 nos mesmos pares (px/quadro)',
                ylabel='Agora: mediana dos pares verificados (px/quadro)',
                title='Comparação sobre exatamente os mesmos intervalos')
    axes[0].legend(fontsize=8); axes[0].grid(alpha=.2)
    y = np.arange(len(summary)); colors = [NEW if included else '#777' for included in summary.included]
    axes[1].barh(y, summary.supported_fraction * 100, color=colors)
    axes[1].set(yticks=y, yticklabels=[r.run + ('' if r.included else ' *') for r in summary.itertuples()],
                xlim=(0, 110), xlabel='Intervalos com suporte espacial e temporal (%)',
                title='Cobertura da nova estimativa — 500 intervalos por registro')
    axes[1].invert_yaxis(); axes[1].tick_params(axis='y', labelsize=8)
    for i, count in enumerate(summary.supported_count):
        axes[1].text(float(summary.supported_fraction.iloc[i]) * 100 + 1, i, f'{count}/500', va='center', fontsize=7)
    fig.suptitle('Deslocamento de estruturas nas imagens: antes e agora\n* Registro auditado, excluído das comparações principais. '
                 'Correspondência de imagem não é uma medida independente da velocidade física.', fontsize=12)
    fig.savefig(ROOT / 'figures/overview.png', dpi=150); plt.close(fig)


def make_gallery(summary):
    rows = []
    for r in summary.itertuples():
        note = '' if r.included else ' · fora do escopo principal'
        rows.append(f'<section><h2>{html.escape(r.run)} · {html.escape(r.pattern)}{note}</h2>'
                    f'<p>{r.supported_count}/500 intervalos verificados; '
                    f'{r.old_clamped_count} intervalos forçados a 1 px no v11.</p>'
                    f'<a href="figures/{r.run}_comparison.png"><img loading="lazy" '
                    f'src="figures/{r.run}_comparison.png" alt="Comparação {r.run}"></a></section>')
    text = '''<!doctype html><html lang="pt-BR"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Deslocamento — antes e agora</title><style>body{font:16px/1.5 system-ui,sans-serif;max-width:1500px;margin:2rem auto;padding:0 1.2rem;color:#172b34;background:#f7f9fa}h1{font-size:2rem}h2{font-size:1.25rem}section{background:white;border:1px solid #dce4e7;border-radius:12px;padding:1.2rem;margin:1.5rem 0}img{width:100%;height:auto}a{color:#087f8c}.lead{max-width:1000px}</style>
<h1>Deslocamento: antes e agora</h1><p class="lead">Estimativa subpixel nas imagens originais, com remoção do fundo estacionário, verificação em outras linhas, sentido inverso e intervalo de dois quadros. Intervalos sem suporte permanecem sem estimativa. Os panoramas novos representam somente segmentos contínuos verificados; não substituem a análise wavelet existente.</p>
<p><a href="README.md">Método e limitações</a> · <a href="results/run_summary.csv">Resultados por registro</a> · <a href="results/displacements.csv">Todos os deslocamentos</a></p>
<img src="figures/overview.png" alt="Visão geral">'''
    (ROOT / 'comparison.html').write_text(text + ''.join(rows) + '</html>\n')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--runs', nargs='+'); args = ap.parse_args()
    start = time.perf_counter()
    for sub in ['results', 'figures', 'panoramas', 'signals', 'intermediate']:
        (ROOT / sub).mkdir(exist_ok=True)
    inv = pd.read_csv(BASE / 'results/input_audit.csv').sort_values('run')
    if args.runs: inv = inv[inv.run.isin(args.runs)]
    summaries = []; transitions = []; temporal_rows = []; inputs = []
    for item in inv.itertuples():
        frames, frame_hash = load_frames(item.run)
        assert sha(DATA / item.run / '__panorama_v11.png') == item.panorama_sha256
        with np.load(BASE / 'intermediate' / f'{item.run}_input.npz') as archive:
            old = {key: archive[key] for key in archive.files}
        assert len(old['raw_shifts']) == 500 and len(old['shifts']) == 501
        rows, temporal = estimate(frames)
        background = frames.mean(axis=0, dtype=np.float64)
        for i, row in enumerate(rows):
            old_dx = -float(old['shifts'][i])  # v11's s has the opposite image-motion sign.
            row.update(run=item.run, pattern=item.pattern, included=bool(item.included),
                       old_raw_s=float(old['raw_shifts'][i]), old_dx=old_dx,
                       old_clamped=bool(old['raw_shifts'][i] <= 0))
            before = image_alignment(frames, background, i, old_dx)
            after = image_alignment(frames, background, i, row['candidate_dx'])
            row.update({f'old_{key}': value for key, value in before.items()})
            row.update({f'new_{key}': value for key, value in after.items()})
            row['dynamic_correlation_change'] = after['dynamic_correlation'] - before['dynamic_correlation']
        for row in temporal: row.update(run=item.run, pattern=item.pattern)
        segment = save_segment(item.run, frames, old, rows)
        example = plot_run(item, frames, background, old, rows, segment)
        df = pd.DataFrame(rows); good = df[df.supported]
        delta = good.dynamic_correlation_change
        saturation = float(np.mean(frames[:, ROI] >= 250))
        summary = dict(run=item.run, pattern=item.pattern, included=bool(item.included),
                       intervals=500, supported_count=len(good), supported_fraction=len(good) / 500,
                       image_pair_support_count=int(df.base_supported.sum()),
                       old_clamped_count=int(df.old_clamped.sum()),
                       supported_old_clamps=int((df.supported & df.old_clamped).sum()),
                       pairwise_supported_old_clamps=int((df.base_supported & df.old_clamped).sum()),
                       old_all_median=finite_median(df.old_dx),
                       old_median_on_supported=finite_median(good.old_dx),
                       new_supported_median=finite_median(good.new_dx),
                       old_check_correlation=finite_mean(good.old_dynamic_correlation),
                       new_check_correlation=finite_mean(good.new_dynamic_correlation),
                       mean_check_correlation_change=finite_mean(delta),
                       alignment_improved=int((delta > 1e-6).sum()),
                       alignment_worsened=int((delta < -1e-6).sum()),
                       alignment_tied=int((abs(delta) <= 1e-6).sum()),
                       saturation_fraction=saturation, saturation_warning=saturation > .20,
                       segment_first=segment['first'] if segment else np.nan,
                       segment_last=segment['last'] if segment else np.nan,
                       segment_intervals=segment['intervals'] if segment else 0,
                       old_segment_width_px=segment['old_width'] if segment else np.nan,
                       new_segment_width_px=segment['new_width'] if segment else np.nan,
                       new_segment_travel_px=segment['total_travel_px'] if segment else np.nan,
                       example_frame_i=example)
        summaries.append(summary); transitions.extend(rows); temporal_rows.extend(temporal)
        inputs.append(dict(run=item.run, native_frames_content_sha256=frame_hash,
                           panorama_sha256=item.panorama_sha256,
                           input_cache_sha256=sha(BASE / 'intermediate' / f'{item.run}_input.npz')))
        print(f'{item.run}: {len(good)}/500 verificados; '
              f'variação média de correlação 2-D {summary["mean_check_correlation_change"]:+.4f}; '
              f'clamps verificados {summary["supported_old_clamps"]}/{summary["old_clamped_count"]}', flush=True)
    summary = pd.DataFrame(summaries)
    pd.DataFrame(transitions).to_csv(ROOT / 'results/displacements.csv', index=False)
    pd.DataFrame(temporal_rows).to_csv(ROOT / 'results/temporal_checks.csv', index=False)
    summary.to_csv(ROOT / 'results/run_summary.csv', index=False)
    overview(summary); make_gallery(summary)
    meta = dict(status='COMPUTED', runs=len(inv), included_runs=int(inv.included.sum()),
                transition_count=len(transitions), input_records=inputs,
                rows_inclusive=[18, 71], fit_rows=FIT_ROWS.tolist(), check_rows=CHECK_ROWS.tolist(),
                positive_dx='rightward in native images', old_dx='negative of v11 strip displacement',
                fixed_parameters=dict(search=SEARCH, minimum_ncc=MIN_CORRELATION,
                    minimum_margin=MIN_MARGIN, minimum_texture_std=MIN_TEXTURE,
                    absolute_tolerance_px=ABS_TOLERANCE, relative_tolerance=REL_TOLERANCE,
                    minimum_panorama_segment_intervals=MIN_SEGMENT_INTERVALS),
                background='Per-pixel temporal mean over all 501 native frames; offline analysis.',
                support='NCC quality + alternate-row agreement + reversed-central-template consistency + at least one consistent lag-2 check; no contradictory comparable lag-2 check.',
                alignment_metric='2-D Pearson correlation on dynamic alternate rows unused for fitting; also influenced by support selection and spatial correlation.',
                rejection='NaN, never replaced with 1 px, a median, interpolation, absolute value, or expected velocity.',
                panorama='Longest contiguous supported monotone segment; first tie; cumulative subpixel positions, midpoint cells, 1-px sampling, no endpoint extrapolation.',
                before_panorama='Exact crop of the original v11 PNG for the same inclusive native frame indices. Endpoint and seam conventions differ from the new construction.',
                comparison_limits='New/old differ in estimator, row selection, background handling, subpixel sampling and construction. This is not a single-factor ablation or physical ground-truth validation.',
                example_selection='Supported interval nearest frame 250, first tie; frame 250 if none supported.',
                source_scope='All 31 audited inputs; original labels and 26-run inclusion flags preserved.',
                elapsed_seconds=time.perf_counter() - start)
    (ROOT / 'results/method.json').write_text(json.dumps(meta, ensure_ascii=False, indent=2) + '\n')
    print(f'Concluído: {len(inv)} registros, {len(transitions)} intervalos, {time.perf_counter() - start:.1f} s', flush=True)


if __name__ == '__main__': main()
