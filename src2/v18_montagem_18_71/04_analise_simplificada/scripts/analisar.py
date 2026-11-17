"""Pearson direto dos perfis completos, mapas e figuras da metodologia atual.

Lê as etapas 02_wavelet e 03_optica e escreve apenas no OUTPUT_ROOT configurado.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import sys
from itertools import combinations

sys.dont_write_bytecode = True
os.environ.setdefault('MPLCONFIGDIR', '/tmp/tcc-v18-pearson-direto-mpl')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import numpy as np
from scipy.stats import pearsonr
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import FuncFormatter, NullLocator

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from runtime import OUTPUT_ROOT, SCRIPTS_ROOT, artifact_path, reference

BASE = OUTPUT_ROOT
OUT = BASE / '04_analise_simplificada'
WAVELET = BASE / '02_wavelet'
INTER = WAVELET / 'intermediario'
DATA = OUT / 'dados'
MAPS = OUT / 'intermediario/pares'
FIGURES = OUT / 'figuras'
CORE = SCRIPTS_ROOT / 'wavelet_core.py'
PLOT = SCRIPTS_ROOT / 'figuras_wavelet.py'
PATTERNS = ['churn', 'wavy annular', 'plug', 'slug', 'dispersed bubble']
PATTERN_LABELS = {'churn':'churn', 'wavy annular':'wavy annular', 'plug':'plug',
                  'slug':'slug', 'dispersed bubble':'dispersed bubbles'}
INPUTS = {}
INPUT_PATHS = {}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def track(path):
    path = Path(path)
    label = reference(path)
    INPUTS[label] = sha(path)
    INPUT_PATHS[label] = path
    return path

def read_json(path):
    return json.loads(track(path).read_text())

def read_csv(path):
    with track(path).open(newline='') as handle:
        return list(csv.DictReader(handle))

def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n')

def write_csv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def source(path):
    return dict(arquivo=reference(path), sha256=sha(path))

def describe(values):
    a = np.asarray(values)
    return dict(n=len(a), media=float(a.mean()), mediana=float(np.median(a)),
                minimo=float(a.min()), maximo=float(a.max()),
                q25=float(np.quantile(a,.25)), q75=float(np.quantile(a,.75)))

def statistics():
    records = read_json(INTER/'aquisicoes.json')
    field_configuration = read_json(INTER/'configuracao.json')
    pairs = read_json(INTER/'pares.json')
    waves = None
    powers, means = [], []
    for row in records:
        with np.load(track(INTER/'individuais'/f"{row['ensaio']}.npz")) as data:
            if waves is None:
                waves = data['wavelengths'].copy()
            np.testing.assert_array_equal(data['wavelengths'], waves)
            powers.append(data['wps_zero'].copy())
            means.append(data['wps_mean_all'].copy())
    power, mean = np.asarray(powers), np.asarray(means)
    assert power.shape == (26,200) and np.isfinite(power).all()
    centered = power-power.mean(axis=1,keepdims=True)
    norm = np.linalg.norm(centered,axis=1,keepdims=True)
    assert np.all(norm>1e-12*np.maximum(1.,np.linalg.norm(power,axis=1,keepdims=True)))
    unit = centered/norm
    matrix = np.clip(unit@unit.T,-1.,1.)
    np.fill_diagonal(matrix,1.)
    np.testing.assert_allclose(matrix,np.corrcoef(power),atol=1e-13,rtol=1e-13)
    errors = [abs(matrix[i,j]-pearsonr(power[i],power[j]).statistic)
              for i,j in combinations(range(26),2)]
    assert max(errors)<1e-12
    runs = [r['ensaio'] for r in records]
    index = {r:i for i,r in enumerate(runs)}
    np.savez_compressed(OUT/'intermediario/descritores.npz',runs=np.asarray(runs),
                        wavelengths=waves,power=power,power_mean_all=mean,
                        centered=centered,unit=unit,r=matrix)
    rows = []
    for pair in pairs:
        a,b = pair['run_x'],pair['run_y']
        i,j = index[a],index[b]
        rows.append(dict(numero_x=i+1,numero_y=j+1,ensaio_x=a,ensaio_y=b,
                         padrao_x=pair['padrao_x'],padrao_y=pair['padrao_y'],tipo=pair['tipo'],
                         comprimento_comum=pair['comprimento'],pearson_r=float(matrix[i,j]),
                         **{key:pair[key] for key in ['wco_maxima','wco_media_integral','wco_maxima_zero',
                            'wco_mediana_zero','pico_lag_px','pico_lambda_px']}))
    assert len(rows)==325 and sum(r['tipo']=='intra' for r in rows)==87
    write_csv(DATA/'pares_pearson_wco.csv', rows)
    write_csv(DATA/'matriz_pearson_todos.csv', [dict(ensaio=run,**{other:float(matrix[i,j])
        for j,other in enumerate(runs)}) for i,run in enumerate(runs)])
    write_csv(DATA/'perfis_espectrais.csv', [dict(ensaio=run,padrao=records[i]['padrao'],
        lambda_px=float(w),wps_zero=float(power[i,k]),potencia_media_todos_lags=float(mean[i,k]),
        valor_centralizado=float(centered[i,k]),descritor_unitario=float(unit[i,k]))
        for i,run in enumerate(runs) for k,w in enumerate(waves)])
    grouped = {}
    for row in rows:
        key = tuple(sorted((row['padrao_x'],row['padrao_y']),key=PATTERNS.index))
        grouped.setdefault(key,[]).append(row)
    groups = [dict(padrao_x=a,padrao_y=b,tipo='intra' if a==b else 'inter',pares=len(group),
                   pearson_media=float(np.mean([r['pearson_r'] for r in group])),
                   pearson_mediana=float(np.median([r['pearson_r'] for r in group])))
              for (a,b),group in sorted(grouped.items(),key=lambda v:(PATTERNS.index(v[0][0]),PATTERNS.index(v[0][1])))]
    aggregate = [dict(tipo=kind,metrica='pearson_r',**describe([r['pearson_r'] for r in rows if r['tipo']==kind]))
                 for kind in ['intra','inter']]
    averages = []
    labels = np.asarray([r['padrao'] for r in records])
    for i,row in enumerate(records):
        others = np.arange(26)!=i
        same = others&(labels==row['padrao'])
        different = labels!=row['padrao']
        averages.append(dict(numero=i+1,ensaio=runs[i],padrao=row['padrao'],
            pearson_media_todos=float(matrix[i,others].mean()),
            pearson_media_intra=float(matrix[i,same].mean()) if same.any() else None,
            pearson_media_inter=float(matrix[i,different].mean()),pares_intra=int(same.sum()),pares_inter=int(different.sum())))
    selected = [max([r for r in rows if r['tipo']=='intra' and r['padrao_x']==pattern],key=lambda r:r['pearson_r'])
                for pattern in PATTERNS[:-1]]
    inter = max([p for p in pairs if p['tipo']=='inter'],key=lambda p:p['wco_maxima'])
    chosen = []
    lookup = {(p['run_x'],p['run_y']):p for p in pairs}
    for a,b in [(r['ensaio_x'],r['ensaio_y']) for r in selected]+[(inter['run_x'],inter['run_y'])]:
        value = dict(lookup[a,b])
        value['pearson_r'] = float(matrix[index[a],index[b]])
        value['arquivo'] = reference(FIGURES / Path(value['arquivo']).name)
        value['origem_metricas_wavelet'] = reference(INTER/'pares'/f'{a}__{b}.json')
        chosen.append(value)
    optical_path = BASE/'03_optica/results/void_fraction/void_fraction_flip_plot_data.csv'
    optical = read_csv(optical_path)
    eligible = {r['run']:r for r in optical if r['eligible_pure_unsaturated']=='True'}
    assert set(eligible).issubset({r['ensaio'] for r in records})
    optical_rows = [dict(ensaio=run,padrao=r['source_void_label'],
                        alpha=float(r['conditional_circle_alpha_spatial_mean']),
                        F1000=float(r['candidate_pixel_flips_per1000px'])) for run,r in eligible.items()]
    optical_examples = []
    for pair in chosen[:2]:
        a,b=pair['run_x'],pair['run_y']
        assert a in eligible and b in eligible
        optical_examples.append(dict(padrao=pair['padrao'],ensaio_x=a,ensaio_y=b,pearson_r=pair['pearson_r'],
            alpha_x=float(eligible[a]['conditional_circle_alpha_spatial_mean']),
            alpha_y=float(eligible[b]['conditional_circle_alpha_spatial_mean']),
            F1000_x=float(eligible[a]['candidate_pixel_flips_per1000px']),
            F1000_y=float(eligible[b]['candidate_pixel_flips_per1000px'])))
    write_csv(DATA/'resumo_por_padrao.csv',groups)
    write_csv(DATA/'resumo_intra_inter.csv',aggregate)
    write_csv(DATA/'medias_por_ensaio.csv',averages)
    write_csv(DATA/'optica_por_ensaio.csv',optical_rows)
    write_csv(DATA/'comparacao_optica_pares_ilustrados.csv',optical_examples)
    write_json(DATA/'selecao_intrapadrao.json',selected)
    write_json(DATA/'pares_ilustrados.json',chosen)
    summary = dict(analise='Pearson direto dos perfis de potência no atraso zero',
        origem_campos_wavelet=reference(WAVELET),ensaios=records,grupos=groups,resumos=aggregate,
        diferenca_global=aggregate[0]['media']-aggregate[1]['media'],medias_por_ensaio=averages,
        pares_intrapadrao_selecionados=selected,pares_ilustrados=chosen,maximo_wco_interpadrões=chosen[-1],
        extensao_minima=min(r['comprimento'] for r in records),extensao_maxima=max(r['comprimento'] for r in records),
        menor_anular=min([r for r in averages if r['padrao']=='wavy annular'],key=lambda r:r['pearson_media_intra']),
        m03q40=next(r for r in averages if r['ensaio']=='m03q40_001'),
        execucao_pares_wavelet_completa=True,campos_wavelet_obtidos_da_etapa_02=True,
        discussion_examples=optical_examples)
    write_json(OUT/'resumo_tcc.json',summary)
    write_json(OUT/'intermediario/verificacao_pearson.json',dict(status='PASS',ensaios=26,escalas=200,
        pares_scipy=325,maior_erro_scipy=max(errors),menor_autovalor=float(np.linalg.eigvalsh(matrix).min()),
        numpy_corrcoef_conferido=True,autopares_excluidos=True,
        descricao='Pearson diretamente sobre wps_zero. Centralização e divisão pela norma são apenas a identidade algébrica de Pearson.'))
    write_json(OUT/'intermediario/configuracao.json',dict(
        ensaios=26,pares=325,intra=87,inter=238,linhas_inclusivas=[18,71],
        wavelengths_px=[50,3000],escalas=200,omega0=field_configuration['omega0'],amostragem_px=1,
        descritor='potência no atraso zero dos sinais completos, sem transformação logarítmica ou ajuste espectral',
        implementacao_pearson='produto escalar de perfis centralizados e divididos pela norma',
        verificacao_pearson=['numpy.corrcoef','scipy.stats.pearsonr, todos os 325 pares'],
        grade_comprimentos_onda='logarítmica, 200 comprimentos de onda entre 50 e 3000 pixels',
        preparo_sinal='inversão de intensidade, retirada de tendência linear, centralização e padronização preservadas',
        autocorrelacao='enviesada, divisor N',wps='abs(CWT(Rxx))',wcs='CWT(Rxy)',
        wco='abs(WCS)^2/(WPSx(0)*WPSy(0))',
        mapas_pares='trecho inicial de comprimento min(Na,Nb), preparado após o corte',
        medias_laterais='média dos módulos em todos os atrasos dos campos completos',
        escolha_intra='maior coeficiente de Pearson direto por padrão com pelo menos dois ensaios',
        escolha_inter='maior coerência pontual entre os 238 pares interpadrões',
        configuracao_campos=source(INTER/'configuracao.json')))
    print('ESTATISTICA_PRONTA '+json.dumps(dict(intra=aggregate[0]['media'],inter=aggregate[1]['media'],
        diferenca=summary['diferenca_global'],pares=selected,optica=optical_examples),ensure_ascii=False),flush=True)
    return records,waves,power,matrix,rows,summary

def pair_maps(pair, core, waves):
    a, b = pair['run_x'], pair['run_y']
    target = MAPS / f'{a}__{b}.npz'
    raw = []
    for run in [a, b]:
        with np.load(track(INTER / 'entradas' / f'{run}.npz')) as data:
            raw.append(data['raw'].copy())
    n = min(map(len, raw))
    x, y = [core.prepare(values[:n], linear=True) for values in raw]
    maps = core.pair(x, y, waves, estimator="biased", engine="morlet")
    co = maps['wco']
    j0 = n - 1
    assert maps['denominator_valid'].all() and np.isfinite(co).all()
    position = np.unravel_index(np.argmax(co), co.shape)
    metrics = dict(comprimento=n, wco_maxima=float(co.max()), wco_minima=float(co.min()),
                   wco_media_integral=float(co.mean()), wco_maxima_zero=float(co[:, j0].max()),
                   wco_mediana_zero=float(np.median(co[:, j0])), pico_lag_px=int(maps['lags'][position[1]]),
                   pico_lambda_px=float(waves[position[0]]))
    reference_metrics = read_json(INTER / 'pares' / f'{a}__{b}.json')
    for key, value in metrics.items():
        np.testing.assert_allclose(value, reference_metrics[key], atol=1e-12, rtol=1e-12, err_msg=key)
    step = max(1, int(np.ceil(len(maps['lags']) / 1800)))
    indices = np.unique(np.r_[np.arange(0, len(maps['lags']), step), j0, len(maps['lags']) - 1])
    compact = {key: maps[key][:, indices] for key in ['wxx', 'wyy', 'wxy', 'wco']}
    compact.update(lags=maps['lags'][indices], wavelengths=waves, scales=maps['scales'], x=x, y=y,
        mean_x_all=np.abs(maps['wxx']).mean(axis=1), mean_y_all=np.abs(maps['wyy']).mean(axis=1),
        max_wxx=np.abs(maps['wxx']).max(), max_wyy=np.abs(maps['wyy']).max(), max_wxy=np.abs(maps['wxy']).max(),
        wps_zero_x=np.abs(maps['wxx'][:, j0]), wps_zero_y=np.abs(maps['wyy'][:, j0]),
        scale_peak_wco=co.max(axis=1), full_lag_count=len(maps['lags']))
    np.savez_compressed(target, **compact)
    write_json(target.with_suffix('.json'), dict(run_x=a, run_y=b,
        modo='campos completos recalculados com o núcleo local; compactação somente para apresentação',
        pearson_r=pair['pearson_r'], origem_metricas=source(INTER / 'pares' / f'{a}__{b}.json'),
        mapa=source(target)))
    return target


def save(fig,path):
    fig.savefig(path,dpi=300,bbox_inches='tight',pad_inches=.08)
    plt.close(fig)

def pattern_label(pattern):
    return r'$\mathit{' + PATTERN_LABELS[pattern].replace(' ', r'\ ') + '}$'

def draw_pearson_matrix(records,matrix):
    """Preserve the matrix style regardless of previously imported renderers."""
    with plt.rc_context(matplotlib.rcParamsDefault):
        return _draw_pearson_matrix(records,matrix)

def _draw_pearson_matrix(records,matrix):
    """Identify the five groups on both axes without changing Pearson values."""
    labels=[r['ensaio'].removesuffix('_001') for r in records]
    colors=['#d34a33','#2377b4','#328652','#9466ac','#b0851c']
    cmap=LinearSegmentedColormap.from_list('pearson',[(0,'#5e3c99'),(.25,'#b2abd2'),(.5,'#fffdf0'),
        (.625,'#d6eef1'),(.75,'#99d4de'),(.875,'#3b93bd'),(1,'#164b72')])
    fig=plt.figure(figsize=(6.6,8.0));ax=fig.add_axes([.20,.32,.75,.62])
    im=ax.imshow(matrix,vmin=-1,vmax=1,cmap=cmap,interpolation='nearest')
    ax.set_xticks(range(len(records)),labels,rotation=90,fontsize=8.5)
    ax.set_yticks(range(len(records)),labels,fontsize=8.5)
    ax.tick_params(length=0,pad=3);ax.tick_params(axis='y',pad=12)
    ax.set_xlabel('Ensaio');ax.set_ylabel('Ensaio')
    handles=[]
    for letter,pattern,color in zip('ABCDE',PATTERNS,colors):
        indices=[i for i,row in enumerate(records) if row['padrao']==pattern]
        start,stop=indices[0],indices[-1]+1
        assert indices==list(range(start,stop)), 'Groups must be contiguous on both axes'
        center=(start+stop-1)/2
        ax.add_patch(Rectangle((start-.5,-1.15),stop-start,.4,
                              facecolor=color,edgecolor='none',clip_on=False))
        ax.add_patch(Rectangle((-1.15,start-.5),.4,stop-start,
                              facecolor=color,edgecolor='none',clip_on=False))
        ax.text(center,-1.70,letter,ha='center',va='center',fontsize=9.5,
                fontweight='bold',color='#243842',clip_on=False)
        ax.text(len(records)+.25,center,letter,ha='center',va='center',fontsize=9.5,
                fontweight='bold',color='#243842',clip_on=False)
        label=f'{letter} · {pattern_label(pattern)} ({len(indices)} ensaio'+('s)' if len(indices)>1 else ')')
        handles.append(Patch(facecolor=color,edgecolor='none',label=label))
        if start:
            ax.axhline(start-.5,color='white',lw=1.4)
            ax.axvline(start-.5,color='white',lw=1.4)
    # Marginal artists must not expand the heatmap's original 26×26 extent.
    ax.set_xlim(-.5,len(records)-.5);ax.set_ylim(len(records)-.5,-.5)
    bar=fig.add_axes([.25,.19,.65,.018])
    cb=fig.colorbar(im,cax=bar,orientation='horizontal',ticks=[-1,-.5,0,.5,1])
    cb.ax.xaxis.set_major_formatter(FuncFormatter(lambda x,_:f'{x:g}'.replace('.',',')))
    cb.set_label('Coeficiente de Pearson r')
    # Matplotlib fills legend columns; this order displays A,B / C,D / E.
    fig.legend(handles=handles[::2]+handles[1::2],loc='lower center',
               bbox_to_anchor=(.55,.015),ncol=2,frameon=False,fontsize=9,
               title='Grupos de ensaios',title_fontsize=9.5,
               handlelength=1.2,columnspacing=1.6,labelspacing=.6)
    path=FIGURES/'matriz_pearson_todos.png';save(fig,path)
    return path

def draw_spectral_profiles(records,waves,power):
    """Draw saved spectral profiles, identifying the two cases discussed in the text."""
    colors=['#d34a33','#2377b4','#328652','#9466ac','#b0851c']
    highlighted={'m03q40_001':'upper right','m06q05_001':'upper left'}
    fig,axes=plt.subplots(5,1,figsize=(6.3,8.2),sharex=True,sharey=True)
    displayed=[]
    for ax,pattern,color in zip(axes,PATTERNS,colors):
        indices=[i for i,r in enumerate(records) if r['padrao']==pattern]
        curves=power[indices]/power[indices].max(axis=1,keepdims=True)
        legend_location=None
        for index,curve in zip(indices,curves):
            run=records[index]['ensaio']
            if run in highlighted:
                ax.plot(waves,curve,color=color,lw=1.2,ls='--',alpha=.95,
                        label=run.removesuffix('_001'))
                legend_location=highlighted[run]
            else:
                ax.plot(waves,curve,color=color,lw=.85,alpha=.42)
        median=np.median(curves,axis=0)
        ax.plot(waves,median,color=color,lw=1.8)
        if legend_location:
            ax.legend(loc=legend_location,frameon=False,fontsize=8.5,handlelength=2.3)
        ax.set_title(f'{pattern_label(pattern)} · {len(indices)} ensaio'+('s' if len(indices)>1 else ''),loc='left',fontsize=10)
        ax.set(xscale='log',xlim=(50,3000),ylim=(0,1.05),ylabel='Potência relativa')
        ax.set_yticks([0,.5,1]);ax.yaxis.set_major_formatter(FuncFormatter(lambda x,_:f'{x:g}'.replace('.',',')))
        ax.set_xticks([50,100,300,1000,3000],['50','100','300','1000','3000'])
        ax.xaxis.set_minor_locator(NullLocator());ax.grid(axis='y',color='#e6ecef',lw=.5)
        displayed.append(dict(padrao=pattern,ensaios=[records[i]['ensaio'] for i in indices],
                              normalizacao_exibicao='cada perfil dividido pelo próprio máximo',
                              mediana_das_curvas_normalizadas=median.tolist()))
    axes[-1].set_xlabel('Comprimento de onda (pixels)')
    fig.subplots_adjust(left=.14,right=.98,top=.97,bottom=.075,hspace=.48)
    profile=FIGURES/'perfis_espectrais_pearson.png';save(fig,profile)
    return profile,displayed

def draw_statistics(records,waves,power,matrix,rows):
    profile,displayed=draw_spectral_profiles(records,waves,power)
    write_json(DATA/'perfis_exibidos.json',displayed)

    matrix_path=draw_pearson_matrix(records,matrix)

    fig,ax=plt.subplots(figsize=(6.3,4.6))
    for column,kind,color in [(1,'intra','#087f8c'),(2,'inter','#b66b27')]:
        values=np.array([r['pearson_r'] for r in rows if r['tipo']==kind])
        jitter=np.random.default_rng(2718+column).uniform(-.17,.17,len(values))
        ax.scatter(column+jitter,values,s=16,alpha=.45,color=color,edgecolor='none')
        q1,q3=np.quantile(values,[.25,.75])
        ax.plot([column-.23,column+.23],[values.mean()]*2,color='#182d37',lw=2)
        ax.plot([column,column],[q1,q3],color='#182d37',lw=2)
    ax.set_xticks([1,2],['Mesmo padrão\n87 pares','Padrões diferentes\n238 pares'])
    ax.set(xlim=(.45,2.55),ylim=(-1.05,1.05),ylabel='Coeficiente de Pearson r')
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x,_:f'{x:g}'.replace('.',',')))
    ax.grid(axis='y',color='#e6ecef',lw=.6)
    fig.subplots_adjust(left=.14,right=.98,top=.98,bottom=.20)
    fig.text(.55,.025,'Traço: média · barra vertical: intervalo interquartil',ha='center',fontsize=9)
    distribution=FIGURES/'correlacoes_dentro_entre_pearson.png';save(fig,distribution)
    return [profile,matrix_path,distribution]

def figures(records,waves,power,matrix,rows,summary):
    track(CORE)
    track(PLOT)
    import wavelet_core as core
    import figuras_wavelet as plots
    plots.ROOT = OUT
    plots.WAVES = waves
    by_id={r['ensaio']:r for r in records}
    entries=[]
    for pair in summary['pares_ilustrados']:
        a,b=pair['run_x'],pair['run_y']
        path=pair_maps(pair,core,waves)
        with np.load(path) as saved:
            maps={k:saved[k] for k in saved.files}
        images=[]
        for run in [a,b]:
            with Image.open(track(artifact_path(by_id[run]['wavelet_path']))) as im:
                images.append(np.asarray(im.convert('L')))
        destination=artifact_path(pair['arquivo'])
        if pair['padrao'] in PATTERN_LABELS:
            title=PATTERN_LABELS[pair['padrao']]
        else:
            title='slug × wavy annular'
        layout=plots.pair_figure(images,maps['x'],maps['y'],maps,runs=(a,b),pattern=title,destination=destination)
        entries.append(dict(tipo=pair['tipo'],ensaio_x=a,ensaio_y=b,pearson_r=pair['pearson_r'],
                            comprimento_comum=pair['comprimento'],layout=layout,
                            origem=source(destination),entrada=source(path)))
        print('FIGURA_PAR '+a+' '+b,flush=True)
    for path in draw_statistics(records,waves,power,matrix,rows):
        entries.append(dict(tipo=path.stem,origem=source(path),entrada=source(OUT/'intermediario/descritores.npz')))
    assert len(entries)==8
    manifest=dict(analise=reference(OUT),script=source(Path(__file__).resolve()),
                  criterio_intra='maior Pearson direto de cada padrão com mais de um ensaio',
                  criterio_inter='maior coerência pontual entre os 238 pares interpadrões',
                  figuras=entries,ajustes_rotulos=['Lag (px) substituído por Atraso (pixels)',
                  'Sinal (σ) substituído por Amplitude padronizada', 'Sinais normalizados substituído por Sinais preparados'],
                  fonte_dos_campos='núcleo local e entradas da etapa 02; cinco mapas calculados e conferidos contra as métricas integrais.')
    write_json(FIGURES/'origens.json',manifest)
    write_json(OUT/'intermediario/figuras.json',manifest)

def main(statistics_only=False):
    INPUTS.clear()
    INPUT_PATHS.clear()
    track(Path(__file__).resolve())
    for path in [DATA,MAPS,FIGURES]:
        path.mkdir(parents=True,exist_ok=True)
    result=statistics()
    if not statistics_only:
        figures(*result)
    assert all(sha(INPUT_PATHS[label])==value for label,value in INPUTS.items())
    manifest=dict(status='estatistica_pronta' if statistics_only else 'complete',script=source(Path(__file__).resolve()),
                  inputs_preserved=True,input_sha256=INPUTS,
                  output_sha256={reference(p):sha(p) for p in sorted(OUT.rglob('*'))
                     if p.is_file() and p.name!='manifesto.json' and not {'__pycache__', 'historico', 'scripts'}.intersection(p.relative_to(OUT).parts)},
                  method='Pearson diretamente sobre WPS no atraso zero dos 26 sinais completos. Sem logaritmo ou remoção de tendência do perfil espectral.',
                  signal_preparation='Inversão, detrend linear e padronização do sinal de intensidade preservados.',
                  fields_source=reference(WAVELET),
                  rows_inclusive=[18,71],optical_analysis_unchanged=True)
    write_json(OUT/'manifesto.json',manifest)
    print('ANALISE_DIRETA_COMPLETA' if not statistics_only else 'ESTATISTICA_DIRETA_COMPLETA',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--statistics-only',action='store_true')
    main(parser.parse_args().statistics_only)
