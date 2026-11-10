"""Visão geral dos caminhos executados no TCC, sem dados experimentais."""
from pathlib import Path
import hashlib
import json
import os

os.environ.setdefault('MPLCONFIGDIR', '/tmp/tcc-ilustracoes-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

HERE = Path(__file__).resolve().parent
INK = '#243842'
COLORS = {'input': '#eef1f3', 'individual': '#e1f0f2',
          'pair': '#fff0da', 'optical': '#eee7f3', 'output': '#eef1f3'}

NODES = [
    ('quadros', 5.7, 10.35, 4.6, .7, 'Quadros das aquisições', 'input'),
    ('espectral', .2, 8.85, 10.15, 1.0,
     'Panoramas da análise espectral\nMédia de intensidade por coluna', 'input'),
    ('optico', 10.85, 8.85, 4.95, 1.0,
     'Panoramas da análise óptica\nJustaposição direta das tiras', 'optical'),
    ('completo', .2, 7.15, 4.85, 1.15,
     'Sinal completo\nPreparação', 'individual'),
    ('corte', 5.5, 7.15, 4.85, 1.15,
     'Dois sinais\nCorte no menor comprimento\nPreparação', 'pair'),
    ('mascara', 10.85, 7.15, 4.95, 1.15,
     'Segmentação\nMáscara de gás inferido', 'optical'),
    ('auto_cwt', .2, 5.5, 4.85, 1.05,
     'Autocorrelação\nCWT', 'individual'),
    ('correlacoes_cwt', 5.5, 5.5, 4.85, 1.05,
     'Duas autocorrelações\nCorrelação cruzada\nCWT', 'pair'),
    ('optica_medidas', 10.85, 5.5, 4.95, 1.05,
     'Fração de vazio estimada\nDensidade de transições', 'optical'),
    ('perfil', .2, 3.85, 4.85, 1.05,
     'Perfil de potência\nno atraso zero', 'individual'),
    ('mapas', 5.5, 3.85, 4.85, 1.05,
     'Mapas de potência,\nespectro cruzado e coerência', 'pair'),
    ('optica_resumos', 10.85, 3.85, 4.95, 1.05,
     'Resumos por ensaio\ne padrão', 'optical'),
    ('pearson', .2, 2.1, 4.85, 1.15,
     'Coeficiente de Pearson\nResumos dentro\ne entre padrões', 'individual'),
    ('leitura_mapas', 5.5, 2.1, 4.85, 1.15,
     'Leitura dos mapas\nResumos da coerência', 'pair'),
    ('interpretacao', 3.3, .35, 9.4, .8,
     'Interpretação conjunta dos resultados', 'output'),
]
EDGES = [
    ('quadros', 'espectral'), ('quadros', 'optico'),
    ('espectral', 'completo'), ('espectral', 'corte'),
    ('optico', 'mascara'), ('completo', 'auto_cwt'),
    ('corte', 'correlacoes_cwt'), ('mascara', 'optica_medidas'),
    ('auto_cwt', 'perfil'), ('correlacoes_cwt', 'mapas'),
    ('optica_medidas', 'optica_resumos'), ('perfil', 'pearson'),
    ('mapas', 'leitura_mapas'), ('pearson', 'interpretacao'),
    ('leitura_mapas', 'interpretacao'), ('optica_resumos', 'interpretacao'),
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9.2,
                         'pdf.fonttype': 42, 'ps.fonttype': 42})
    fig = plt.figure(figsize=(16/2.54, 11.5/2.54))
    ax = fig.add_axes([0, 0, 1, 1], xlim=(0, 16), ylim=(0, 11.5))
    ax.set_axis_off()
    by_id = {node[0]: node for node in NODES}
    for node in NODES:
        name, x, y, width, height, label, group = node
        ax.add_patch(FancyBboxPatch((x, y), width, height,
                     boxstyle='round,pad=0.025,rounding_size=0.09',
                     facecolor=COLORS[group], edgecolor=INK, linewidth=.65,
                     zorder=2))
        ax.text(x+width/2, y+height/2, label, ha='center', va='center',
                linespacing=1.18, color=INK, zorder=3)
    for source, target in EDGES:
        s, t = by_id[source], by_id[target]
        start = (s[1]+s[3]/2, s[2]-.035)
        end = (t[1]+t[3]/2, t[2]+t[4]+.04)
        if source == 'espectral':
            start = (end[0], s[2]-.035)
        if target == 'interpretacao':
            if source == 'pearson':
                end = (4.15, t[2]+t[4]+.04)
            elif source == 'leitura_mapas':
                end = (8.0, t[2]+t[4]+.04)
            else:
                start = (s[1]+s[3]/2, s[2]-.035)
                end = (12.0, t[2]+t[4]+.04)
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle='-|>',
                     mutation_scale=9, linewidth=.8, color=INK,
                     connectionstyle='arc3,rad=0', zorder=1))
    for ext in ['pdf', 'png']:
        fig.savefig(HERE/f'fluxograma.{ext}', dpi=300, facecolor='white')
    plt.close(fig)
    metadata = {
        'tipo': 'esquema do procedimento executado, sem dados experimentais',
        'dimensoes_cm': [16, 11.5], 'fonte_pt': 9.2,
        'nos': [dict(id=n[0],texto=n[5],grupo=n[6]) for n in NODES],
        'arestas': EDGES,
        'criterios': [
            'Panoramas espectrais e ópticos são versões distintas da mesma montagem.',
            'Nos pares, o corte antecede a preparação.',
            'A CWT é aplicada às correlações.',
            'Pearson usa perfis dos sinais completos.',
            'WCO usa as correlações dos trechos comuns, sem receber o descritor completo de Pearson.',
            'A análise óptica parte da máscara de gás inferido.',
        ],
        'script_sha256': sha(Path(__file__)),
        'arquivos_sha256': {f'fluxograma.{e}':sha(HERE/f'fluxograma.{e}')
                            for e in ['pdf', 'png']},
    }
    (HERE/'fluxograma_origem.json').write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2)+'\n')


if __name__ == '__main__':
    main()
