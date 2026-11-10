#!/usr/bin/env python3
"""Render four complete v18 acquisitions without recomputing signal processing.

Run from any directory with python3 -B gerar_panoramas.py. Only the four PNGs
and origens.json in this directory are written. Inputs and historical outputs
are read only. The acquisitions illustrate existing within-pattern examples.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import platform
import sys

sys.dont_write_bytecode = True
os.environ.setdefault('MPLCONFIGDIR', '/tmp/tcc-v18-wavelet-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.text import Text
from matplotlib.ticker import FuncFormatter, MaxNLocator
import numpy as np
from PIL import Image

DEST = Path(__file__).resolve().parent
REPO = DEST.parents[2]
BASE = REPO/'src2/v18_montagem_18_71'
INTER = BASE/'02_wavelet/intermediario'
SELECTION = (
    ('plug', 'Plug', 'm05q01_001'),
    ('slug', 'Slug', 'm04q01_001'),
    ('churn', 'Churn', 'm04q20_001'),
    ('wavy annular', 'Anular ondulado', 'm06q40_001'),
)
INK, TEAL = '#263c48', '#087f8c'
DPI = 300


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def source(path):
    return dict(arquivo=str(path.relative_to(REPO)), sha256=sha(path))


def decimal(value, _=None):
    return f'{value:g}'.replace('.', ',')


def render(run, pattern, image, raw, prepared, destination):
    """Same three-panel idea as the historical panorama_figure, without notes."""
    fig = plt.figure(figsize=(6.3, 4.2), facecolor='white')
    grid = fig.add_gridspec(3, 1, left=.14, right=.98, top=.94, bottom=.13,
                           height_ratios=[.65, 1., 1.], hspace=.55)
    ax = fig.add_subplot(grid[0])
    picture = ax.imshow(image, cmap='gray', vmin=0, vmax=255, aspect='auto',
        extent=[-.5, image.shape[1]-.5, 79.5, -.5], interpolation='nearest')
    ax.set_xlim(-.5, image.shape[1]-.5)
    ax.set_yticks([0, 79], ['0', '79'])
    ax.set_ylabel('Altura\n(pixels)')
    ax.tick_params(axis='x', labelbottom=False)
    ax.xaxis.set_major_locator(MaxNLocator(5, integer=True))
    for boundary in (17.5, 71.5):
        ax.axhline(boundary, color=TEAL, lw=1)
    ax.set_title('Panorama', loc='left')
    image_axis = ax
    np.testing.assert_array_equal(np.asarray(picture.get_array()), image)
    for slot, values, title, ylabel, color in (
        (grid[1], raw, 'Intensidade média por coluna', 'Intensidade\nmédia', INK),
        (grid[2], prepared, 'Sinal preparado', 'Amplitude\npadronizada', TEAL),
    ):
        ax = fig.add_subplot(slot, sharex=image_axis)
        x = np.arange(len(raw))
        line, = ax.plot(x, values, color=color, lw=.8)
        ax.set(xlim=(-.5, len(raw)-.5), ylabel=ylabel)
        if title == 'Sinal preparado':
            ax.set_xlabel('Posição no panorama (pixels)', labelpad=2)
        else:
            ax.tick_params(axis='x', labelbottom=False)
        ax.set_title(title, loc='left')
        ax.xaxis.set_major_locator(MaxNLocator(5, integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(3))
        ax.yaxis.set_major_formatter(FuncFormatter(decimal))
        ax.grid(axis='y', color='#e8edef', lw=.6)
        np.testing.assert_array_equal(line.get_xdata(), x)
        np.testing.assert_array_equal(line.get_ydata(), values)
    assert not fig.texts  # Panel titles only, no notes or redundant main title.
    assert all(not ax.texts for ax in fig.axes)
    fig.canvas.draw()
    box = fig.get_tightbbox(fig.canvas.get_renderer()).padded(.06)
    insertion_scale = 16/(box.width*2.54)
    texts = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text()]
    min_font = min(t.get_fontsize() for t in texts)*insertion_scale
    insertion_height = box.height*2.54*insertion_scale
    assert min_font >= 8 and insertion_height <= 24
    fig.savefig(destination, dpi=DPI, bbox_inches='tight', pad_inches=.06)
    plt.close(fig)
    with Image.open(destination) as saved:
        saved.verify()
    with Image.open(destination) as saved:
        dimensions = list(saved.size)
    return dict(dimensions_px=dimensions, dpi=DPI, insertion_width_cm=16,
                insertion_height_cm=insertion_height, minimum_font_at_insertion_pt=min_font,
                no_notes_or_observations=True, full_arrays_plotted_without_decimation=True)


def main():
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.titlesize': 10,
        'axes.labelsize': 10, 'xtick.labelsize': 10, 'ytick.labelsize': 10,
        'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .65,
        'savefig.facecolor': 'white', 'axes.titleweight': 'medium', 'axes.titlepad': 6,
    })
    acquisitions_path = INTER/'aquisicoes.json'
    summary_path = BASE/'02_wavelet/resumo_tcc.json'
    manifest_path = BASE/'01_panoramas/manifesto.json'
    historical_plot = REPO/'src2/v16_spatial2_18_71/scripts/figuras.py'
    numerical_manifest = BASE/'02_wavelet/manifesto.json'
    acquisitions = {r['ensaio']:r for r in json.loads(acquisitions_path.read_text())}
    summary = json.loads(summary_path.read_text())
    assert summary['execucao_pares_wavelet_completa']
    pairs = {p['padrao_x']:p for p in summary['pares_intrapadrao_selecionados']}
    inputs = [acquisitions_path, summary_path, manifest_path, historical_plot, numerical_manifest]
    for _, _, run in SELECTION:
        inputs.extend([REPO/acquisitions[run]['wavelet_path'], INTER/'entradas'/f'{run}.npz',
                       INTER/'individuais'/f'{run}.npz', INTER/'individuais'/f'{run}.json'])
    before = {str(p.relative_to(REPO)):sha(p) for p in inputs}
    figures = []
    for pattern, display, run in SELECTION:
        row = acquisitions[run]
        assert row['padrao'] == pattern and pairs[pattern]['ensaio_x'] == run
        panorama_path = REPO/row['wavelet_path']
        raw_path, prepared_path = INTER/'entradas'/f'{run}.npz', INTER/'individuais'/f'{run}.npz'
        assert sha(panorama_path) == row['hashes']['wavelet']
        with Image.open(panorama_path) as image_file:
            image = np.asarray(image_file.convert('L')).copy()
        with np.load(raw_path) as data:
            raw = data['raw'].copy()
        with np.load(prepared_path) as data:
            prepared = data['prepared'].copy()
        individual = json.loads(prepared_path.with_suffix('.json').read_text())
        assert individual['npz_sha256'] == sha(prepared_path)
        expected_raw = image[18:72].sum(axis=0, dtype=np.uint64)/54
        np.testing.assert_array_equal(raw, expected_raw)
        assert image.shape == (80, row['comprimento']) and raw.shape == prepared.shape == (image.shape[1],)
        assert np.isfinite(raw).all() and np.isfinite(prepared).all()
        destination = DEST/f'{pattern.replace(" ", "_")}_{run}.png'
        layout = render(run, display, image, raw, prepared, destination)
        figures.append(dict(**source(destination), pattern=pattern, display_pattern=display, run=run,
            full_length_px=len(raw), layout=layout,
            inputs=[source(panorama_path), source(raw_path), source(prepared_path)],
            verification=dict(raw_equal_to_rows_18_through_71_mean=True,
                              prepared_loaded_unchanged_from_current_individual_npz=True,
                              plotted_arrays_identical_to_loaded_inputs=True,
                              png_integrity_verified=True)))
    assert before == {str(p.relative_to(REPO)):sha(p) for p in inputs}
    provenance = dict(generator=source(Path(__file__).resolve()),
        reproduction='python3 -B TCC/figuras/v18_panoramas/gerar_panoramas.py',
        purpose='Four complete panorama and signal illustrations for section 4.2.',
        selection='First acquisition of each already selected within-pattern pair for plug, slug, churn and wavy annular. '
                  'Illustration only, with no additional statistical selection criterion.',
        historical_layout_reference=source(historical_plot),
        numerical_data_recomputed=False, existing_inputs_and_manifests_unchanged=True,
        input_hashes=before, figures=figures,
        runtime=dict(python=platform.python_version(), numpy=np.__version__, matplotlib=matplotlib.__version__))
    (DEST/'origens.json').write_text(json.dumps(provenance, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({f['run']:f['layout'] for f in figures}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
