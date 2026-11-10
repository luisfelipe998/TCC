"""Pearson assinado de descritores fixos: todos os 26 ensaios elegíveis da v16."""
import ast
from datetime import datetime, timezone
import hashlib
import importlib.metadata as metadata
from itertools import combinations
import json
import os
from pathlib import Path
import re
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = ROOT.parent/'v16_spatial2_18_71'
RESULTS = HERE/'results'
FIGURES = HERE/'figures'
os.environ.setdefault('MPLCONFIGDIR','/tmp/estatistica5-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from PIL import Image
from scipy.signal import correlate
from scipy.stats import pearsonr

sys.path.insert(0,str(SOURCE/'scripts'))
import wavelet_core as core

PATTERNS = ['churn','wavy annular','plug','slug','dispersed bubbles']
COLORS = {'churn':'#d34a33','wavy annular':'#2377b4','plug':'#328652',
          'slug':'#9466ac','dispersed bubbles':'#ac8015'}
WAVELENGTHS = np.geomspace(50,3000,200)
FACTOR = 4*np.pi/(6+np.sqrt(38))
PLAN = ROOT/'METODOLOGIA_TODOS_CONTRA_TODOS.md'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(8*1024*1024),b''):
            h.update(chunk)
    return h.hexdigest()


def save_json(path,value):
    Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+chr(10))


def zero_profile(x):
    """CWT da autocorrelação em zero por soma temporal independente da FFT-CWT."""
    n = len(x)
    autocorrelation = correlate(x,x,mode='full',method='fft')/n
    lag = np.arange(-n+1,n)
    values = []
    for a in WAVELENGTHS/FACTOR:
        t = lag/a
        kernel = np.pi**(-.25)/np.sqrt(a)*np.exp(-.5*t*t)*np.cos(6*t)
        values.append(abs(np.sum(autocorrelation*kernel)))
    return np.asarray(values)


def residuals(p,wavelengths=WAVELENGTHS):
    p = np.asarray(p,dtype=float)
    if p.ndim != 2 or not np.isfinite(p).all() or np.any(p < 0):
        raise ValueError('Perfis inválidos.')
    logp = np.log(p+1e-12)
    design = np.column_stack([np.ones(len(wavelengths)),np.log(wavelengths)])
    coefficients = np.linalg.lstsq(design,logp.T,rcond=None)[0].T
    z = logp-coefficients@design.T
    assert np.max(abs(z@design)) < 1e-9
    return z,coefficients


def pearson_matrix(z):
    z = np.asarray(z,dtype=float)
    if z.ndim != 2 or not np.isfinite(z).all():
        raise ValueError('Descritores inválidos.')
    centered = z-z.mean(axis=1,keepdims=True)
    norms = np.linalg.norm(centered,axis=1,keepdims=True)
    if np.any(norms <= 1e-12*np.maximum(1.,np.linalg.norm(z,axis=1,keepdims=True))):
        raise ValueError('Pearson indefinido para vetor constante ou degenerado.')
    unit = centered/norms
    r = unit@unit.T
    assert r.min() >= -1-1e-12 and r.max() <= 1+1e-12
    # Corrigir exclusivamente arredondamento de produto de vetores unitários.
    r = np.clip(r,-1,1)
    np.fill_diagonal(r,1.)
    return r,unit


def validate_formulas():
    x = np.array([-1.,0.,1.])
    z = np.array([x,-x,[1.,-2.,1.],7*x+13])
    r,_ = pearson_matrix(z)
    assert np.allclose(r[0],[1,-1,0,1],atol=1e-12)
    try:
        pearson_matrix(np.array([x,np.ones(3)]))
    except ValueError:
        pass
    else:
        raise AssertionError('Vetor constante aceito.')
    rng = np.random.default_rng(20261024)
    raw = rng.normal(size=127)+np.linspace(2,20,127)
    design = np.column_stack([np.ones(len(raw)),np.arange(len(raw))])
    inverted = 255-raw
    residual = inverted-design@np.linalg.lstsq(design,inverted,rcond=None)[0]
    expected = (residual-residual.mean())/residual.std(ddof=0)
    error_prepare = float(np.max(abs(core.prepare(raw)-expected)))
    assert error_prepare < 1e-11
    x = core.prepare(raw)
    literal = np.array([sum(x[k+m]*x[k] for k in range(len(x)-m))/len(x)
                        for m in range(len(x))])
    error_correlation = float(np.max(abs(core.correlation(x,x)[len(x)-1:]-literal)))
    assert error_correlation < 1e-12
    assert abs(FACTOR-core.FOURIER_FACTOR) < 1e-15
    profiles = np.exp(rng.normal(size=(6,200)))
    z1,_ = residuals(profiles)
    z2,_ = residuals(profiles,WAVELENGTHS/FACTOR)
    assert np.allclose(z1,z2,atol=1e-12)
    return dict(status='PASS',known_pearson_minus_one_zero_one=True,constant_descriptor_rejected=True,
        positive_affine_invariance=True,prepare_vs_independent_ols_max_error=error_prepare,
        autocorrelation_vs_literal_sum_max_error=error_correlation,scale_conversion_checked=True,
        residual_scale_invariance_checked=True,scope='Controles sintéticos, separados dos resultados experimentais.')


def extract_all():
    inventory_path = SOURCE/'results/input_audit.csv'
    inventory = pd.read_csv(inventory_path)
    eligible = inventory[inventory.included].copy()
    assert len(inventory) == 31 and len(eligible) == 26
    assert set(inventory.loc[~inventory.included,'run']) == {
        'm03q01_001','m03q80_001','m04q05_001','m05q05_002','m07q10_001'}
    eligible['pattern_order'] = eligible.pattern.map({p:i for i,p in enumerate(PATTERNS)})
    assert eligible.pattern_order.notna().all()
    eligible = eligible.sort_values(['pattern_order','run']).reset_index(drop=True)
    eligible['number'] = np.arange(1,27)
    common = int(eligible.n.min())
    assert common == 3171
    eligible[['number','run','pattern','n']].to_csv(RESULTS/'indice_ensaios.csv',index=False)
    inventory[~inventory.included][['run','pattern','n']].to_csv(RESULTS/'ensaios_excluidos.csv',index=False)
    paths = [inventory_path,SOURCE/'results/profiles.csv',SOURCE/'results/equal_length_windows.csv',
             SOURCE/'results/analysis_config.json',SOURCE/'scripts/wavelet_core.py']
    hashes = {str(path.relative_to(ROOT.parent)):sha(path) for path in paths}
    full_saved = pd.read_csv(SOURCE/'results/profiles.csv',float_precision='round_trip')
    windows_saved = pd.read_csv(SOURCE/'results/equal_length_windows.csv',float_precision='round_trip')
    p_common,p_full,checks = [],[],[]
    for row in eligible.itertuples():
        image_path = ROOT.parent/'data_processed'/row.run/'__panorama_v11.png'
        input_path = SOURCE/'intermediate'/f'{row.run}_input.npz'
        for path in [image_path,input_path]:
            hashes[str(path.relative_to(ROOT.parent))] = sha(path)
        assert hashes[str(image_path.relative_to(ROOT.parent))] == row.panorama_sha256
        with Image.open(image_path) as image:
            pixels = np.asarray(image.convert('L'))
        assert pixels.shape == (row.height,row.n)
        expected = pixels[18:72].sum(axis=0,dtype=np.uint64)/54
        with np.load(input_path,allow_pickle=False) as archive:
            raw = archive['raw']
        assert np.array_equal(raw,expected)
        assert hashlib.sha256(raw.astype('<f8').tobytes()).hexdigest() == row.raw_array_sha256
        x = core.prepare(raw[:common])
        p = zero_profile(x)
        w,_ = core.transform_correlation(core.correlation(x,x),WAVELENGTHS)
        fft_profile = abs(w[:,common-1])
        fft_error = float(np.linalg.norm(p-fft_profile)/np.linalg.norm(fft_profile))
        assert fft_error < 1e-10
        del w
        cache = windows_saved[(windows_saved.run==row.run) & (windows_saved.window=='prefix')].sort_values('wavelength_px')
        assert len(cache) == 200 and cache.n.eq(common).all()
        assert np.allclose(cache.wavelength_px,WAVELENGTHS,rtol=1e-14)
        cached_error = float(np.linalg.norm(p-cache.zero_lag_wps)/np.linalg.norm(p))
        assert cached_error < 1e-10
        p_complete = zero_profile(core.prepare(raw))
        full = full_saved[(full_saved.run==row.run) & (full_saved.method=='v16_biased')].sort_values('wavelength_px')
        assert len(full) == 200 and full.included.all() and full.pattern.eq(row.pattern).all()
        assert np.allclose(full.wavelength_px,WAVELENGTHS,rtol=1e-14)
        full_error = float(np.linalg.norm(p_complete-full.zero_lag_wps)/np.linalg.norm(p_complete))
        assert full_error < 1e-10
        p_common.append(p);p_full.append(p_complete)
        checks.append(dict(run=row.run,n_original=int(row.n),n_common=common,
            pixels_max_abs_difference=float(np.max(abs(raw-expected))),
            direct_vs_fft_profile_relative_l2=fft_error,
            direct_vs_saved_common_relative_l2=cached_error,direct_vs_saved_full_relative_l2=full_error))
        if row.number % 7 == 0 or row.number == 26:
            print(f'Entradas e perfis conferidos: {row.number}/26',flush=True)
    pd.DataFrame(checks).to_csv(RESULTS/'verificacao_por_ensaio.csv',index=False)
    return eligible,np.asarray(p_common),np.asarray(p_full),hashes,checks


def summarize(index,r_main,r_full):
    rows = []
    for i,j in combinations(range(len(index)),2):
        a,b = index.iloc[i],index.iloc[j]
        rows.append(dict(index_x=i+1,index_y=j+1,run_x=a.run,run_y=b.run,
            pattern_x=a.pattern,pattern_y=b.pattern,same_pattern=bool(a.pattern==b.pattern),
            n_main_px=3171,n_full_x=int(a.n),n_full_y=int(b.n),
            pearson_r=float(r_main[i,j]),pearson_r_full=float(r_full[i,j]),
            delta_full_minus_main=float(r_full[i,j]-r_main[i,j])))
    pairs = pd.DataFrame(rows)
    assert len(pairs) == 325 and pairs.same_pattern.sum() == 87
    assert (~pairs.same_pattern).sum() == 238
    pairs.to_csv(RESULTS/'pearson_325_pares.csv',index=False)
    for name,r in [('matriz_pearson',r_main),('matriz_pearson_completo',r_full)]:
        pd.DataFrame(r,index=index.run,columns=index.run).to_csv(RESULTS/f'{name}.csv',index_label='run')
    aggregate = []
    for field,condition in [('pearson_r','3171_px_inicio'),('pearson_r_full','comprimento_completo')]:
        for same,label in [(True,'Mesmo padrão'),(False,'Padrões diferentes')]:
            values = pairs.loc[pairs.same_pattern==same,field]
            aggregate.append(dict(condition=condition,comparison=label,n_pairs=len(values),
                minimum=float(values.min()),q25=float(values.quantile(.25)),median=float(values.median()),
                q75=float(values.quantile(.75)),maximum=float(values.max()),negative_pairs=int((values<0).sum())))
    aggregate = pd.DataFrame(aggregate)
    aggregate.to_csv(RESULTS/'resumo_dentro_entre_padroes.csv',index=False)
    group_rows = []
    for i,pattern_x in enumerate(PATTERNS):
        for pattern_y in PATTERNS[i:]:
            selected = pairs[(pairs.pattern_x==pattern_x) & (pairs.pattern_y==pattern_y)]
            values = selected.pearson_r
            group_rows.append(dict(pattern_x=pattern_x,pattern_y=pattern_y,n_pairs=len(values),
                minimum=float(values.min()) if len(values) else np.nan,
                median=float(values.median()) if len(values) else np.nan,
                maximum=float(values.max()) if len(values) else np.nan))
    groups = pd.DataFrame(group_rows)
    groups.to_csv(RESULTS/'resumo_combinacoes_padroes.csv',index=False)
    primary = aggregate[aggregate.condition=='3171_px_inicio'].set_index('comparison')
    highest = pairs.loc[pairs.pearson_r.idxmax()]
    lowest = pairs.loc[pairs.pearson_r.idxmin()]
    summary = dict(n_runs=26,n_pairs=325,n_intra_pairs=87,n_inter_pairs=238,
        n_common_px=3171,n_scales=200,wavelength_limits_px=[50,3000],color_limits=[-1,1],
        patterns={p:int((index.pattern==p).sum()) for p in PATTERNS},
        pearson_min=float(pairs.pearson_r.min()),pearson_max=float(pairs.pearson_r.max()),
        negative_pairs=int((pairs.pearson_r<0).sum()),
        median_intra=float(primary.loc['Mesmo padrão','median']),
        median_inter=float(primary.loc['Padrões diferentes','median']),
        median_delta=float(primary.loc['Mesmo padrão','median']-primary.loc['Padrões diferentes','median']),
        highest_pair=dict(run_x=highest.run_x,run_y=highest.run_y,r=float(highest.pearson_r),same_pattern=bool(highest.same_pattern)),
        lowest_pair=dict(run_x=lowest.run_x,run_y=lowest.run_y,r=float(lowest.pearson_r),same_pattern=bool(lowest.same_pattern)),
        inter_pairs_above_intra_median=int((pairs.loc[~pairs.same_pattern,'pearson_r'] >= primary.loc['Mesmo padrão','median']).sum()),
        median_abs_change_with_full=float(pairs.delta_full_minus_main.abs().median()),
        max_abs_change_with_full=float(pairs.delta_full_minus_main.abs().max()),
        inference='Descrição dos coeficientes, sem testes de separação ou classificação.')
    save_json(RESULTS/'summary.json',summary)
    return pairs,aggregate,groups,summary


def verify_matrices(zs,rs):
    verification = {}
    for name,z,r in zip(['common','full'],zs,rs):
        expected = np.corrcoef(z)
        numpy_error = float(np.max(abs(r-expected)))
        scipy_error = max(abs(r[i,j]-pearsonr(z[i],z[j]).statistic)
                          for i,j in combinations(range(26),2))
        eigen_min = float(np.linalg.eigvalsh(r).min())
        assert numpy_error < 1e-12 and scipy_error < 1e-12
        assert r.shape == (26,26) and np.isfinite(r).all()
        assert np.allclose(r,r.T,atol=1e-12) and np.array_equal(r.diagonal(),np.ones(26))
        assert eigen_min > -1e-10
        verification[name] = dict(numpy_max_abs_error=numpy_error,scipy_max_abs_error=float(scipy_error),
            minimum_eigenvalue=eigen_min,finite_cells=int(np.isfinite(r).sum()),
            symmetric=True,unit_diagonal=True,positive_semidefinite_within_tolerance=True)
    return verification


def save_figure(fig,name):
    fig.savefig(FIGURES/f'{name}.png',dpi=180)
    fig.savefig(FIGURES/f'{name}.pdf')
    plt.close(fig)


def draw_matrix(index,r,filename,subtitle):
    fig,ax = plt.subplots(figsize=(14,10.8))
    fig.subplots_adjust(left=.20,right=.855,bottom=.075,top=.84)
    palette = LinearSegmentedColormap.from_list('purple_cream_blue', [
        (0.,'#5e3c99'),(.25,'#b2abd2'),(.5,'#fffdf0'),(.625,'#d6eef1'),
        (.75,'#99d4de'),(.875,'#3b93bd'),(1.,'#164b72')], N=256)
    image = ax.imshow(r,vmin=-1,vmax=1,cmap=palette,interpolation='nearest')
    ax.set_xticks(range(26),range(1,27),fontsize=8)
    ax.set_yticks(range(26),range(1,27),fontsize=8)
    ax.set(xlabel='Ensaio')
    start = 0
    for pattern in PATTERNS:
        count = int((index.pattern==pattern).sum())
        center = start+(count-1)/2
        ax.text(-.065,center,f'{pattern} ({count})',transform=ax.get_yaxis_transform(),
                ha='right',va='center',fontsize=10 if count==1 else 11,clip_on=False)
        if count > 1:
            ax.text((center+.5)/26,1.016,f'{pattern} ({count})',transform=ax.transAxes,
                    ha='center',va='bottom',fontsize=11)
        else:
            ax.annotate('dispersed bubbles (1)',xy=(center,-.5),xytext=(center,-2.3),
                        ha='right',va='bottom',fontsize=10,annotation_clip=False,
                        arrowprops=dict(arrowstyle='-',color='#555555',lw=.8))
        if start:
            ax.axhline(start-.5,color='white',lw=1.7);ax.axvline(start-.5,color='white',lw=1.7)
        start += count
    colorbar = fig.colorbar(image,ax=ax,fraction=.046,pad=.04,ticks=np.linspace(-1,1,9))
    colorbar.set_label('Coeficiente de correlação de Pearson, r',fontsize=12)
    fig.suptitle('Matriz de correlação de Pearson entre ensaios'+chr(10)+subtitle,fontsize=16,y=.965)
    save_figure(fig,filename)


def make_figures(index,p,z,pairs,r,rfull):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    draw_matrix(index,r,'matriz_pearson_todos','Descritores residuais da WPS no lag zero · segmentos de 3.171 px')
    draw_matrix(index,rfull,'matriz_pearson_completo','Descritores residuais da WPS no lag zero · registros completos')
    fig,axes = plt.subplots(1,2,figsize=(13,5),layout='constrained',sharex=True)
    for k,(same,label,color) in enumerate([(True,'Mesmo padrão','#367d56'),(False,'Padrões diferentes','#8264a6')]):
        values = pairs.loc[pairs.same_pattern==same,'pearson_r'].to_numpy()
        y = k+.13*np.sin(np.arange(len(values))*2.399963)
        axes[0].scatter(values,y,s=15,color=color,alpha=.5)
        axes[0].plot([values.min(),values.max()],[k,k],color='black',lw=1)
        axes[0].scatter([np.median(values)],[k],marker='D',s=55,color='black',zorder=5)
    axes[0].set_yticks([0,1],['Mesmo padrão (87 pares)','Padrões diferentes (238 pares)'])
    axes[0].set_title('Sobreposição das correlações observadas')
    for k,pattern in enumerate(PATTERNS):
        values = pairs.loc[(pairs.pattern_x==pattern) & (pairs.pattern_y==pattern),'pearson_r'].to_numpy()
        if len(values):
            axes[1].scatter(values,k+.12*np.sin(np.arange(len(values))*2.399963),s=23,color=COLORS[pattern],alpha=.65)
            axes[1].plot([values.min(),values.max()],[k,k],color='black',lw=1)
            axes[1].scatter([np.median(values)],[k],marker='D',s=50,color='black',zorder=5)
        else:
            axes[1].text(-.93,k,'1 ensaio; sem par intrapadrão',va='center',fontsize=9)
    axes[1].set_yticks(range(5),PATTERNS)
    axes[1].set_title('Comparações dentro de cada padrão')
    axes[1].set_ylim(-.5,4.5)
    for ax in axes:
        ax.set(xlim=(-1,1),xlabel='Pearson r');ax.axvline(0,color='#777777',lw=.7)
        ax.grid(axis='x',alpha=.2);ax.invert_yaxis()
    fig.suptitle('Pontos = pares; losango = mediana; linha = mínimo–máximo; sem intervalos de confiança',fontsize=11)
    save_figure(fig,'correlacoes_dentro_entre')
    fig,axes = plt.subplots(2,5,figsize=(16,7),layout='constrained',sharex=True,sharey='row')
    normalized = p/p.max(axis=1,keepdims=True)
    centered = z-z.mean(axis=1,keepdims=True)
    centered /= np.linalg.norm(centered,axis=1,keepdims=True)
    for column,pattern in enumerate(PATTERNS):
        selection = index.pattern.to_numpy()==pattern
        for row,values in enumerate([normalized,centered]):
            ax = axes[row,column]
            for vector in values[selection]:
                ax.plot(WAVELENGTHS,vector,color=COLORS[pattern],alpha=.4,lw=.8)
            ax.plot(WAVELENGTHS,np.median(values[selection],axis=0),color=COLORS[pattern],lw=2)
            ax.set_xscale('log');ax.grid(alpha=.18)
            if row == 0:ax.set_title(f'{pattern} (n={selection.sum()})',fontsize=11)
            else:ax.set_xlabel('Comprimento de onda (px)')
    axes[0,0].set_ylabel('WPS / máximo de cada perfil')
    axes[1,0].set_ylabel('Descritor centrado / norma')
    fig.suptitle('Perfis individuais em 3.171 px; traço grosso = mediana'+chr(10)+
                 'A normalização da WPS acima serve ao gráfico; Pearson usa o descritor abaixo.',fontsize=12)
    save_figure(fig,'perfis_comparados')
    fig,ax = plt.subplots(figsize=(7.3,6.4),layout='constrained')
    for same,label,color in [(True,'Mesmo padrão','#367d56'),(False,'Padrões diferentes','#8264a6')]:
        selected = pairs[pairs.same_pattern==same]
        ax.scatter(selected.pearson_r,selected.pearson_r_full,label=label,s=20,c=color,alpha=.55)
    ax.plot([-1,1],[-1,1],color='black',lw=1,label='Coeficientes iguais')
    ax.set(xlim=(-1,1),ylim=(-1,1),xlabel='Pearson r: primeiros 3.171 px',
        ylabel='Pearson r: registro completo',title='Sensibilidade descritiva ao trecho usado')
    ax.legend(fontsize=9);ax.grid(alpha=.2)
    save_figure(fig,'sensibilidade_comprimento')


def main():
    RESULTS.mkdir(exist_ok=True);FIGURES.mkdir(exist_ok=True)
    plan_hash = sha(PLAN)
    formula_checks = validate_formulas()
    print('Fórmulas conferidas; iniciando os 26 ensaios.',flush=True)
    index,p,pfull,hashes,checks = extract_all()
    z,coefficients = residuals(p);zfull,coefficients_full = residuals(pfull)
    r,u = pearson_matrix(z);rfull,ufull = pearson_matrix(zfull)
    matrix_checks = verify_matrices([z,zfull],[r,rfull])
    np.savez_compressed(RESULTS/'descritores_por_ensaio.npz',runs=index.run.to_numpy(dtype=str),
        patterns=index.pattern.to_numpy(dtype=str),wavelengths_px=WAVELENGTHS,n_common_px=3171,
        n_original=index.n.to_numpy(),profile=p,residual=z,unit_descriptor=u,coefficients=coefficients,
        profile_full=pfull,residual_full=zfull,unit_descriptor_full=ufull,coefficients_full=coefficients_full)
    pairs,aggregate,groups,summary = summarize(index,r,rfull)
    validation = dict(status='PASS',formula_controls=formula_checks,matrix_checks=matrix_checks,
        images_checked=26,raw_arrays_exact=True,common_profiles_vs_fft_checked=26,
        source_profiles_checked=52,max_profile_fft_relative_l2=max(c['direct_vs_fft_profile_relative_l2'] for c in checks),
        n_unique_pairs=325,n_intra_pairs=87,n_inter_pairs=238,missing_matrix_cells=0,
        source_files_unchanged=True,panorama_reconstruction_reexecuted=False,
        method_plan_sha256=plan_hash)
    save_json(RESULTS/'validation.json',validation)
    make_figures(index,p,z,pairs,r,rfull)
    from reporting import write_reports
    write_reports(index,pairs,aggregate,groups,summary,validation)
    # Verificações finais, sem alterar fontes ou inferir dados ausentes.
    for path in HERE.glob('*.py'):ast.parse(path.read_text())
    for path in FIGURES.glob('*.png'):
        with Image.open(path) as image:image.verify()
    assert len(list(FIGURES.glob('*.png'))) == 5
    broken=[]
    documents=[ROOT/'README.md',ROOT/'RESUMO_EXECUTIVO.md',PLAN]
    for document in documents:
        for target in re.findall(r'[]][(]([^)]+)[)]',document.read_text()):
            path=document.parent/target.split('#')[0]
            if '://' not in target and not path.exists() and path != RESULTS/'manifest.json':
                broken.append(target)
    assert not broken,broken
    assert sha(PLAN) == plan_hash
    for name,digest in hashes.items():assert sha(ROOT.parent/name) == digest,name
    save_json(RESULTS/'manifest.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),
        python=sys.version,packages={p:metadata.version(p) for p in ['numpy','scipy','pandas','matplotlib','Pillow','PyWavelets']},
        source_sha256=hashes,source_files_unchanged=True,method_plan_sha256=plan_hash,
        script_sha256={path.name:sha(path) for path in HERE.glob('*.py')},
        output_sha256={str(path.relative_to(ROOT)):sha(path) for path in
            list(RESULTS.iterdir())+list(FIGURES.iterdir())+documents+[ROOT/'RESUMO_EXECUTIVO.pdf']
            if path.is_file() and path.name != 'manifest.json'},
        final_checks=dict(status='PASS',figures_verified=5,broken_links=broken)))
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__ == '__main__':main()
