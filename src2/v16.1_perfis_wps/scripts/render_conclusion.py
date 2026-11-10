"""Três figuras essenciais da conclusão, derivadas das tabelas auditadas."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/tmp/tcc-conclusao-mpl')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd
from PIL import Image

from profile_core import ROOT, BASE, wc

OUT = ROOT/'figures/conclusao'
RESULTS = ROOT/'results'
TEAL, ORANGE, PURPLE = '#217b80', '#c27b17', '#77629d'
DECIMAL = FuncFormatter(lambda x, _: f'{x:g}'.replace('.', ','))
PERCENT = FuncFormatter(lambda x, _: f'{100*x:.0f}%')


def percentage(value):
    return f'{100*value:.1f}%'.replace('.', ',')


def style_profile(ax):
    ax.set_xscale('log')
    ax.set_xticks([50, 100, 200, 500, 1000, 2000])
    ax.xaxis.set_major_formatter(DECIMAL)
    ax.yaxis.set_major_formatter(DECIMAL)
    ax.set_xlabel('Comprimento de onda equivalente (px)')
    ax.set_ylabel('Resposta relativa\n(densidade em log da escala)')
    ax.grid(alpha=.18)


def main():
    assert json.loads((RESULTS/'conclusion_audit.json').read_text())['status'] == 'PASS'
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
        'axes.titlesize': 12, 'axes.spines.top': False, 'axes.spines.right': False,
        'figure.facecolor': 'white', 'savefig.facecolor': 'white', 'pdf.fonttype': 42})
    normalized = pd.read_csv(RESULTS/'profiles_normalized.csv')
    primary = normalized[(normalized.variante == 'principal') &
                         (normalized.dominio == 'principal')]
    ranks = pd.read_csv(RESULTS/'ranking_summary.csv').set_index('metodo')
    influence = pd.read_csv(RESULTS/'acquisition_deletion.csv')
    cases = pd.read_csv(RESULTS/'selected_cases.csv')
    manifest = dict(created_utc=datetime.now(timezone.utc).isoformat(), figures=[])

    def save(fig, name, detail):
        paths = []
        for ext in ['png', 'pdf']:
            path = OUT/f'{name}.{ext}'
            fig.savefig(path, dpi=190, bbox_inches='tight')
            paths.append(str(path.relative_to(ROOT)))
        manifest['figures'].append(dict(name=name, files=paths, **detail))
        plt.close(fig)

    # Representante de plug já predeterminado no protocolo, sem nova escolha por resultado.
    run = 'm05q01_001'
    config = json.loads((ROOT/'config.json').read_text())
    assert run in config['representatives']
    panorama_path = ROOT.parent/'data_processed'/run/'__panorama_v11.png'
    with Image.open(panorama_path) as image:
        panorama = np.asarray(image.convert('L'))
    with np.load(BASE/'intermediate'/f'{run}_input.npz') as stored:
        raw = stored['raw'].copy()
    n = len(raw)
    profile = primary[primary.aquisicao == run].sort_values('comprimento_onda_px')
    fig, axes = plt.subplots(3, 1, figsize=(10.8, 8.1), layout='constrained',
                             gridspec_kw={'height_ratios': [.8, 1, 1.15]})
    axes[0].imshow(panorama, cmap='gray', vmin=0, vmax=255, aspect='auto',
                   extent=[-.5, n-.5, panorama.shape[0]-.5, -.5])
    for row in [17.5, 71.5]:
        axes[0].axhline(row, color='#d24b3e', lw=1.2)
    axes[0].set_title('1. Imagem: a faixa entre as linhas vermelhas forneceu o sinal', loc='left')
    axes[0].set_ylabel('Linha da imagem')
    axes[0].set_xlabel('Posição no panorama (px)')
    axes[0].set_xlim(-.5, n-.5)
    axes[1].plot(np.arange(n), wc.prepare(raw), lw=.8, color=TEAL)
    axes[1].set_title('2. Sinal: cada coluna virou um número, depois preparado', loc='left')
    axes[1].set_xlabel('Posição no panorama (px)')
    axes[1].set_ylabel('Intensidade\npadronizada')
    axes[1].set_xlim(-.5, n-.5)
    axes[1].yaxis.set_major_formatter(DECIMAL)
    axes[1].grid(alpha=.18)
    axes[2].plot(profile.comprimento_onda_px, profile.densidade_log, color=ORANGE, lw=2)
    axes[2].set_title('3. Perfil: a WPS resumiu as escalas presentes no sinal', loc='left')
    style_profile(axes[2])
    axes[2].set_ylim(bottom=0)
    fig.suptitle('Da imagem ao perfil comparado\nExemplo de bolhas alongadas (plug)', fontsize=15)
    save(fig, '01_imagem_sinal_perfil', dict(aquisicao=run, rows=[18, 71],
        n=n, variant='principal', domain='principal',
        note='Imagem com proporção ajustada para leitura; o perfil é normalizado em log(lambda).'))

    methods = ['wco_v16_18_71', 'wps_prefixo_par', 'wps_principal']
    labels = ['v16: escore WCO\nTrecho comum do par · 50–3.000 px',
              'v16.1: perfis no mesmo trecho\nTrecho comum do par · 50–3.000 px',
              'v16.1: perfis principais\nRegistro inteiro · 50–aprox. 2.296 px']
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.9), sharey=True, layout='constrained')
    plotted = []
    for ax, key, title in zip(axes, ['A_pares', 'A_balanceada'],
                             ['Cada par com o mesmo peso', 'Padrões/combinações com pesos iguais']):
        values = ranks.loc[methods, key].to_numpy()
        ax.barh(np.arange(3), values, height=.5, color=[PURPLE, '#6faaa6', TEAL])
        for i, method in enumerate(methods):
            sensitivity = influence[influence.metodo == method][key]
            low, high = float(sensitivity.min()), float(sensitivity.max())
            ax.plot([low, high], [i, i], color='#262626', lw=1.5)
            ax.plot([low, high], [i, i], '|', color='#262626', ms=9)
            ax.text(values[i]+.012, i+.32, percentage(values[i]), fontsize=11)
            plotted.append(dict(method=method, weighting=key, value=float(values[i]),
                                deletion_min=low, deletion_max=high))
        ax.set_yticks(np.arange(3), labels)
        ax.set_ylim(2.65, -.55)
        ax.set_xlim(0, 1)
        ax.set_xticks(np.linspace(0, 1, 6))
        ax.xaxis.set_major_formatter(PERCENT)
        ax.axvline(.5, color='#555', ls='--', lw=1)
        ax.set_xlabel('Concordância de ordenação A')
        ax.set_title(title, fontsize=11)
        ax.grid(axis='x', alpha=.18)
    fig.suptitle('Ordenação de semelhança: resultado principal e comparação no mesmo trecho\n'
                 'Traços pretos: variação ao retirar uma aquisição; não são intervalos de confiança',
                 fontsize=13)
    save(fig, '02_comparacao_essencial', dict(values=plotted,
        direction='WCO maior = mais semelhante; Hellinger menor = mais semelhante',
        note='50% significa ausência de preferência de ordenação. A não é acurácia.'))

    titles = ['Controle cruzado mais semelhante', 'Par intrapadrão menos semelhante']
    plain_titles = ['Padrões diferentes, perfis próximos', 'Mesmo padrão, perfis diferentes']
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.8), sharey=True, layout='constrained')
    illustrated = []
    global_top = 0.
    for ax, title, plain in zip(axes, titles, plain_titles):
        case = cases[cases.tipo == title].iloc[0]
        for i, run in enumerate([case.aquisicao_x, case.aquisicao_y], start=1):
            g = primary[primary.aquisicao == run].sort_values('comprimento_onda_px')
            global_top = max(global_top, float(g.densidade_log.max()))
            ax.plot(g.comprimento_onda_px, g.densidade_log, color=[TEAL, ORANGE][i-1],
                    lw=2, label=f'Aquisição {i}: {g.padrao.iloc[0]}')
        value = f'{case.H:.3f}'.replace('.', ',')
        ax.set_title(f'{plain}\nDistância de Hellinger = {value}')
        style_profile(ax)
        ax.legend(fontsize=9, loc='upper right')
        illustrated.append(dict(type=title, run_x=case.aquisicao_x,
                                run_y=case.aquisicao_y, H=float(case.H)))
    for ax in axes:
        ax.set_ylim(0, global_top*1.32)
    fig.suptitle('Semelhança de perfil não identifica o padrão em todos os casos\n'
                 'Exemplos extremos escolhidos após os resultados, apenas para explicação', fontsize=13)
    save(fig, '03_limites_semelhanca', dict(cases=illustrated,
        variant='principal', domain='principal'))

    sources = [RESULTS/name for name in ['profiles_normalized.csv', 'ranking_summary.csv',
               'acquisition_deletion.csv', 'selected_cases.csv', 'conclusion_audit.json']]
    sources += [Path(__file__), ROOT/'config.json', panorama_path,
                BASE/'intermediate'/'m05q01_001_input.npz']
    manifest['source_sha256'] = {
        str(p.relative_to(ROOT.parent)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    manifest['output_sha256'] = {
        path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
        for figure in manifest['figures'] for path in figure['files']}
    (RESULTS/'conclusion_figures.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n')
    print('Três figuras de conclusão geradas em PNG e PDF.', flush=True)


if __name__ == '__main__':
    main()
