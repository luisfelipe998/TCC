"""Matriz descritiva de Pearson usando somente os 87 pares principais da v16.

Não executa testes de separação, classificação ou novas comparações entre padrões.
Os perfis são extraídos dos mapas salvos de cada par, no lag zero.
"""
from datetime import datetime, timezone
import hashlib
import importlib.metadata as metadata
from itertools import combinations
import json
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent/'v16_spatial2_18_71'
RESULTS = HERE/'results'
FIGURES = HERE/'figures'
os.environ.setdefault('MPLCONFIGDIR','/tmp/estatistica5-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

PATTERNS = ['churn','wavy annular','plug','slug']
WAVELENGTHS = np.geomspace(50,3000,200)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(8*1024*1024),b''):
            h.update(chunk)
    return h.hexdigest()


def save_json(path,value):
    Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+chr(10))


def residual(profile,wavelengths):
    """Resíduo log-linear; preserva a definição do descritor da estatistica4."""
    profile = np.asarray(profile,dtype=float)
    if profile.ndim != 1 or not np.isfinite(profile).all() or np.any(profile < 0):
        raise ValueError('Perfil não finito ou negativo.')
    log_profile = np.log(profile+1e-12)
    design = np.column_stack([np.ones(len(wavelengths)),np.log(wavelengths)])
    coefficients = np.linalg.lstsq(design,log_profile,rcond=None)[0]
    z = log_profile-design@coefficients
    assert np.max(abs(design.T@z)) < 1e-9
    return z


def pearson_distance(x,y):
    """Pearson assinado e distância normalizada; recusa vetores degenerados."""
    x,y = np.asarray(x,dtype=float),np.asarray(y,dtype=float)
    if x.shape != y.shape or x.ndim != 1 or not np.isfinite([x,y]).all():
        raise ValueError('Vetores incompatíveis.')
    a,b = x-x.mean(),y-y.mean()
    na,nb = np.linalg.norm(a),np.linalg.norm(b)
    if na <= 1e-12*max(1.,np.linalg.norm(x)) or nb <= 1e-12*max(1.,np.linalg.norm(y)):
        raise ValueError('Pearson indefinido para vetor constante.')
    r = float(np.sum(a*b)/(na*nb))
    if not -1-1e-12 <= r <= 1+1e-12:
        raise ArithmeticError('Pearson fora do domínio.')
    # Ajuste apenas de arredondamento; correlações negativas são preservadas.
    r = min(1.,max(-1.,r))
    return r,(1-r)/2


def validate_formula():
    x = np.array([-1.,0.,1.])
    for y,expected_r,expected_d in [(x,1.,0.),(-x,-1.,1.),(np.array([1.,-2.,1.]),0.,.5)]:
        r,d = pearson_distance(x,y)
        assert abs(r-expected_r) < 1e-12 and abs(d-expected_d) < 1e-12
    r,d = pearson_distance(x,7*x+13)
    assert abs(r-1) < 1e-12 and d < 1e-12
    try:
        pearson_distance(x,np.ones(3))
    except ValueError:
        pass
    else:
        raise AssertionError('Vetor constante aceito.')
    return dict(status='PASS',cases=['r=1: d=0','r=0: d=0.5','r=-1: d=1',
        'invariância afim positiva','vetor constante recusado'],
        scope='Controles matemáticos, não resultados experimentais.')


def independent_profile_check(pair_rows):
    """Confere os perfis salvos em um par de cada padrão por soma temporal."""
    sys.path.insert(0,str(SOURCE/'scripts'))
    from wavelet_core import prepare,correlation,FOURIER_FACTOR
    checks = []
    for pattern in PATTERNS:
        row = pair_rows[pair_rows.pattern == pattern].iloc[0]
        path = SOURCE/'intermediate'/f'{row.run_x}__{row.run_y}_maps.npz'
        with np.load(path,allow_pickle=False) as maps:
            zero = int(np.flatnonzero(maps['lags'] == 0)[0])
            for run,field in [(row.run_x,'wxx'),(row.run_y,'wyy')]:
                input_path = SOURCE/'intermediate'/f'{run}_input.npz'
                with np.load(input_path,allow_pickle=False) as data:
                    signal = prepare(data['raw'][:int(row.n)])
                rxx = correlation(signal,signal,'biased')
                lag = np.arange(-len(signal)+1,len(signal))
                expected = []
                for a in WAVELENGTHS/FOURIER_FACTOR:
                    t = lag/a
                    kernel = np.pi**(-.25)/np.sqrt(a)*np.exp(-.5*t*t)*np.cos(6*t)
                    expected.append(abs(np.sum(rxx*kernel)))
                saved = abs(maps[field][:,zero])
                relative = float(np.linalg.norm(saved-expected)/np.linalg.norm(saved))
                assert relative < 1e-10
                checks.append(dict(pattern=pattern,run=run,n=int(row.n),relative_l2=relative))
    return checks


def render_matrix(matrix,index):
    n = len(index)
    fig,ax = plt.subplots(figsize=(12,10.7))
    fig.subplots_adjust(left=.07,right=.84,bottom=.17,top=.80)
    cmap = plt.get_cmap('viridis_r').copy()
    cmap.set_bad('#dddddd')
    im = ax.imshow(np.ma.masked_invalid(matrix),cmap=cmap,vmin=0,vmax=1,
                   interpolation='nearest',origin='upper',aspect='equal')
    ax.set_xticks(range(n),range(1,n+1),fontsize=8)
    ax.set_yticks(range(n),range(1,n+1),fontsize=8)
    ax.set(xlabel='Ensaio',ylabel='Ensaio')
    ax.xaxis.label.set_size(13);ax.yaxis.label.set_size(13)
    start = 0
    for pattern in PATTERNS:
        count = int((index.pattern == pattern).sum())
        ax.text((start+(count-1)/2+.5)/n,1.018,f'{pattern} ({count})',
                transform=ax.transAxes,ha='center',va='bottom',fontsize=12)
        if start:
            ax.axvline(start-.5,color='white',lw=2)
            ax.axhline(start-.5,color='white',lw=2)
        start += count
    cax = fig.add_axes([.87,.19,.035,.61])
    cb = fig.colorbar(im,cax=cax,ticks=np.linspace(0,1,6))
    cb.set_label('Distância de Pearson normalizada (1 − r) / 2',fontsize=12,labelpad=13)
    fig.suptitle('Distância de Pearson entre os perfis WPS no lag zero'+chr(10)+
                 '25 ensaios · 87 pares principais da v16 · resíduo log-linear',
                 fontsize=16,y=.94)
    ax.legend(handles=[Patch(facecolor='#dddddd',edgecolor='#aaaaaa',
              label='Cinza: comparação fora dos 87 pares principais')],
              loc='upper center',bbox_to_anchor=(.5,-.078),frameon=False,fontsize=11)
    fig.text(.49,.074,'0: r = 1  |  0,5: r = 0  |  1: r = −1',ha='center',fontsize=11)
    fig.text(.49,.045,'Trechos preservados por par: N = menor comprimento dos dois ensaios.'+chr(10)+
             'Diagonal = identidade; números dos ensaios em resultados/indice_ensaios.csv.',
             ha='center',fontsize=10)
    # O diretório efetivo é results; a legenda em português aponta para o nome real.
    fig.texts[-1].set_text(fig.texts[-1].get_text().replace('resultados/','results/'))
    fig.savefig(FIGURES/'matriz_pearson_normalizada.png',dpi=200)
    fig.savefig(FIGURES/'matriz_pearson_normalizada.pdf')
    plt.close(fig)


def write_report(summary,index,pairs,excluded):
    rows = ['| Nº | Ensaio | Padrão |','| --- | --- | --- |']
    rows += [f'| {r.number} | {r.run} | {r.pattern} |' for r in index.itertuples()]
    minimum = pairs.loc[pairs.distance_normalized.idxmin()]
    maximum = pairs.loc[pairs.distance_normalized.idxmax()]
    report = f'''# Pearson dos pares principais da v16

## Gráfico

[PNG](figures/matriz_pearson_normalizada.png) · [PDF](figures/matriz_pearson_normalizada.pdf)

![Distância de Pearson normalizada](figures/matriz_pearson_normalizada.png)

Pearson r varia de −1 a 1. A dissimilaridade usual 1−r varia de 0 a 2. Para apresentar a faixa 0–1 solicitada, foi usada **d = (1−r)/2**, com limites de cor fixos em 0 e 1. Não houve normalização pelo mínimo e máximo da amostra nem uso de |r| ou r². Correlações negativas, se presentes, continuam negativas no CSV e correspondem a d > 0,5.

| r de Pearson | d normalizada | Leitura |
| --- | --- | --- |
| 1 | 0 | Mesma forma linear após centralização, salvo escala positiva |
| 0 | 0,5 | Sem associação linear entre as coordenadas dos perfis |
| −1 | 1 | Formas linearmente opostas |

Zero não exige igualdade de amplitude dos perfis; Pearson é invariante a deslocamento e escala positiva. A distância é uma dissimilaridade de forma, não um percentual de diferença física ou uma probabilidade. O termo distância não implica que a desigualdade triangular seja satisfeita.

## Escopo: somente os pares principais

Foram usados exatamente os **87 pares** de [pair_methods.csv](../v16_spatial2_18_71/results/pair_methods.csv), filtrados pelo método v16_biased, correspondentes aos 87 mapas em intermediate: churn 45, wavy annular 36, plug 3 e slug 3. Esses pares envolvem **25 ensaios**: 10 churn, 9 wavy annular, 3 plug e 3 slug.

As células cinza não foram calculadas: indicam pares fora desse conjunto, e não distância zero. O conjunto contém apenas comparações intrapadrão. A diagonal foi fixada em zero pela identidade de cada vetor consigo mesmo e não entra na contagem de 87 pares. As duas metades da matriz repetem os mesmos pares por simetria.

O ensaio elegível m07q01_001, rotulado dispersed bubbles, não tem par intrapadrão nos mapas principais e ficou fora da matriz. Não foi reclassificado como plug. Os cinco ensaios excluídos na v16 também permanecem fora. As comparações adicionais do experimento de pré-processamento não foram incorporadas.

## Perfil e cálculo preservados

Para cada par foram extraídos P_x(a)=|wxx(0,a)| e P_y(a)=|wyy(0,a)| dos coeficientes complexos salvos no mapa correspondente. O índice do lag zero foi localizado no vetor de lags; não se tomou o meio da imagem como índice presumido. Foram mantidos os 200 comprimentos de onda de 50 a 3.000 px.

Como na análise de forma anterior, cada perfil foi transformado em log(P+10^−12) e dele foi retirado o ajuste de uma reta contra log(comprimento de onda). Pearson foi calculado entre os dois vetores residuais. Essa operação é remoção de tendência log-linear, sem afirmar branqueamento estatístico ou independência entre escalas. O gráfico não mostra Pearson dos sinais brutos, média da WPS em lag, nem alternative WCO.

A v16 usa os primeiros N=min(N_x,N_y) pixels de cada par, retirando tendência e normalizando depois do corte. Esses comprimentos foram preservados, variando de {summary['minimum_pair_n_px']} a {summary['maximum_pair_n_px']} px. **O mesmo ensaio pode ter um descritor diferente conforme o parceiro**, porque o trecho muda. A matriz é uma tabela simétrica de comparações por par; não é apresentada como matriz de correlação de um único vetor fixo por ensaio. Não foi imposta a janela global de 7.001 px da estatistica4.

## Resultados observados

- r variou de {summary['pearson_r_min']:.6f} a {summary['pearson_r_max']:.6f}.
- d normalizada variou de {summary['distance_min']:.6f} a {summary['distance_max']:.6f}; mediana {summary['distance_median']:.6f}.
- Menor distância: {minimum.run_x} / {minimum.run_y}, {minimum.pattern}, d={minimum.distance_normalized:.6f}, r={minimum.pearson_r:.6f}.
- Maior distância: {maximum.run_x} / {maximum.run_y}, {maximum.pattern}, d={maximum.distance_normalized:.6f}, r={maximum.pearson_r:.6f}.

Os valores baixos descrevem semelhança de forma dentro dos padrões nos trechos usados. **Esta matriz não avalia separação entre padrões**, porque as comparações entre padrões foram omitidas conforme o escopo solicitado. Não foram executados testes de separação, classificação ou novos diagnósticos além de Pearson. Os coeficientes não vêm acompanhados de valores-p que tratem as escalas como observações independentes.

## Identificação dos ensaios

{chr(10).join(rows)}

## Arquivos e reprodução

- [Coeficientes e distâncias dos 87 pares](results/pearson_pairs.csv): r assinado, 1−r, (1−r)/2, N e origem de cada mapa.
- [Matriz 0–1](results/matriz_distancia_normalizada.csv), [matriz r](results/matriz_pearson_r.csv), [comprimento por par](results/matriz_n_px.csv) e [índice dos ensaios](results/indice_ensaios.csv).
- [Perfis por par](results/pair_descriptors.npz): magnitudes e resíduos, com IDs dos dois ensaios.
- [Validação](results/validation.json), [resumo](results/summary.json) e [manifesto com hashes](results/manifest.json).
- [Explicação do gráfico ECDF anterior](EXPLICACAO_ECDF.md).

Na raiz do repositório:

    PYTHONDONTWRITEBYTECODE=1 python3 src2/estatistica5/run_pearson.py

Os mapas e tabelas de origem são somente lidos. Os campos vazios dos CSVs representam comparações não incluídas; a diagonal zero é identidade. As versões efetivamente usadas estão no manifesto.
'''
    (HERE/'README.md').write_text(report)


def main():
    RESULTS.mkdir(exist_ok=True);FIGURES.mkdir(exist_ok=True)
    formula = validate_formula()
    inventory_path = SOURCE/'results/input_audit.csv'
    pairs_path = SOURCE/'results/pair_methods.csv'
    inventory = pd.read_csv(inventory_path)
    source_pairs = pd.read_csv(pairs_path)
    source_pairs = source_pairs[source_pairs.method == 'v16_biased'].copy()
    assert len(source_pairs) == 87 and not source_pairs.duplicated(['run_x','run_y']).any()
    eligible = inventory[inventory.included].set_index('run')
    expected = {p for _,g in eligible.reset_index().groupby('pattern') for p in combinations(sorted(g.run),2)}
    assert set(zip(source_pairs.run_x,source_pairs.run_y)) == expected
    participating = set(source_pairs.run_x)|set(source_pairs.run_y)
    assert len(participating) == 25
    order = [run for pattern in PATTERNS for run in sorted(participating) if eligible.loc[run,'pattern'] == pattern]
    assert len(order) == 25
    index = pd.DataFrame(dict(number=np.arange(1,26),run=order,
                             pattern=[eligible.loc[run,'pattern'] for run in order]))
    index.to_csv(RESULTS/'indice_ensaios.csv',index=False)
    hashes = {str(p.relative_to(HERE.parent)):sha(p) for p in [inventory_path,pairs_path,
               SOURCE/'scripts/wavelet_core.py',SOURCE/'scripts/run_analysis.py']}
    dmat = np.full((25,25),np.nan);rmat = np.full((25,25),np.nan);nmat = np.full((25,25),np.nan)
    np.fill_diagonal(dmat,0.);np.fill_diagonal(rmat,1.)
    locations = {run:i for i,run in enumerate(order)}
    rows,px,py,zx,zy = [],[],[],[],[]
    errors = []
    for row in source_pairs.itertuples():
        assert row.run_x in eligible.index and row.run_y in eligible.index
        assert eligible.loc[row.run_x,'pattern'] == eligible.loc[row.run_y,'pattern'] == row.pattern
        assert row.n == min(eligible.loc[row.run_x,'n'],eligible.loc[row.run_y,'n'])
        path = SOURCE/'intermediate'/f'{row.run_x}__{row.run_y}_maps.npz'
        hashes[str(path.relative_to(HERE.parent))] = sha(path)
        with np.load(path,allow_pickle=False) as maps:
            wavelengths = maps['wavelengths']
            assert np.array_equal(wavelengths,WAVELENGTHS)
            zeros = np.flatnonzero(maps['lags'] == 0)
            assert len(zeros) == 1
            p1,p2 = [abs(maps[key][:,zeros[0]]) for key in ['wxx','wyy']]
        z1,z2 = residual(p1,wavelengths),residual(p2,wavelengths)
        r,d = pearson_distance(z1,z2)
        error = float(abs(r-pearsonr(z1,z2).statistic))
        errors.append(error)
        assert error < 1e-12
        i,j = locations[row.run_x],locations[row.run_y]
        dmat[i,j]=dmat[j,i]=d;rmat[i,j]=rmat[j,i]=r;nmat[i,j]=nmat[j,i]=row.n
        rows.append(dict(run_x=row.run_x,run_y=row.run_y,pattern=row.pattern,
            index_x=i+1,index_y=j+1,n_px=row.n,pearson_r=r,distance_1_minus_r=1-r,
            distance_normalized=d,source_map=str(path.relative_to(HERE.parent))))
        px.append(p1);py.append(p2);zx.append(z1);zy.append(z2)
    pairs = pd.DataFrame(rows)
    pairs.to_csv(RESULTS/'pearson_pairs.csv',index=False)
    for name,matrix in [('matriz_distancia_normalizada',dmat),('matriz_pearson_r',rmat),('matriz_n_px',nmat)]:
        pd.DataFrame(matrix,index=order,columns=order).to_csv(RESULTS/f'{name}.csv',index_label='run')
    np.savez_compressed(RESULTS/'pair_descriptors.npz',run_x=pairs.run_x.to_numpy(dtype=str),
        run_y=pairs.run_y.to_numpy(dtype=str),wavelengths_px=WAVELENGTHS,n_px=pairs.n_px.to_numpy(),
        profile_x=np.array(px),profile_y=np.array(py),residual_x=np.array(zx),residual_y=np.array(zy))
    print('87 pares calculados; conferindo perfis e representação...',flush=True)
    checks = independent_profile_check(source_pairs)
    for check in checks:
        path = SOURCE/'intermediate'/f"{check['run']}_input.npz"
        hashes[str(path.relative_to(HERE.parent))] = sha(path)
    assert np.allclose(dmat,dmat.T,equal_nan=True)
    assert np.count_nonzero(np.isfinite(dmat[np.triu_indices(25,1)])) == 87
    assert np.count_nonzero(np.isnan(dmat)) == 426
    assert np.nanmin(dmat) >= 0 and np.nanmax(dmat) <= 1
    excluded = eligible.loc[~eligible.index.isin(participating)].reset_index()[['run','pattern']]
    summary = dict(source_method='v16_biased',n_pairs=87,n_runs=25,n_scales=200,
        patterns={p:dict(runs=int((index.pattern==p).sum()),pairs=int((pairs.pattern==p).sum())) for p in PATTERNS},
        distance_definition='(1 - signed Pearson r) / 2',color_limits=[0,1],
        minimum_pair_n_px=int(pairs.n_px.min()),maximum_pair_n_px=int(pairs.n_px.max()),
        pearson_r_min=float(pairs.pearson_r.min()),pearson_r_max=float(pairs.pearson_r.max()),
        distance_min=float(pairs.distance_normalized.min()),distance_max=float(pairs.distance_normalized.max()),
        distance_median=float(pairs.distance_normalized.median()),
        excluded_eligible_without_pair=excluded.to_dict('records'),
        uncomputed_off_diagonal_cells=426,diagonal='identidade; não é um par experimental',
        interpretation='Similaridade intrapadrão; sem comparação entre padrões.')
    save_json(RESULTS/'summary.json',summary)
    render_matrix(dmat,index)
    write_report(summary,index,pairs,excluded)
    from PIL import Image
    with Image.open(FIGURES/'matriz_pearson_normalizada.png') as image:
        image.verify()
    for name,digest in hashes.items():
        assert sha(HERE.parent/name) == digest,name
    save_json(RESULTS/'validation.json',dict(status='PASS',formula_controls=formula,
        pearson_vs_scipy_max_abs_error=max(errors),independent_profile_checks=checks,
        exact_source_pair_set=True,matrix_symmetric=True,negative_correlations_preserved=True,
        color_limits=[0,1],source_files_unchanged=True))
    save_json(RESULTS/'manifest.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),
        python=sys.version,packages={p:metadata.version(p) for p in ['numpy','scipy','pandas','matplotlib','Pillow']},
        source_sha256=hashes,script_sha256=sha(Path(__file__)),
        output_sha256={str(p.relative_to(HERE)):sha(p) for p in HERE.rglob('*')
                       if p.is_file() and '__pycache__' not in p.parts and p.name not in ['manifest.json','run_pearson.py']}))
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__ == '__main__':
    main()
