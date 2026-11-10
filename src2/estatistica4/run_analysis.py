"""Análise reproduzível dos dados existentes; escreve somente em estatistica4."""
from datetime import datetime, timezone
import ast
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
SOURCE = HERE.parent/'v16_spatial2_18_71'
RESULTS = HERE/'results'
FIGURES = HERE/'figures'
os.environ.setdefault('MPLCONFIGDIR', '/tmp/estatistica4-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy.stats import spearmanr

from statistics_core import (WAVELENGTHS, anosim, binary_arrangements, distances,
                             exact_upper_p, holm, permanova, residual_descriptor,
                             unit_rows, zero_lag_profile)
from validate_analysis import validate
sys.path.insert(0,str(SOURCE/'scripts'))
from wavelet_core import prepare, transform_correlation, correlation
from verify_final import verify_maps

PRIMARY_PATTERNS = ['churn','wavy annular']
COLORS = {'churn':'#cf482c','wavy annular':'#2475b1','plug':'#2d8c55',
          'slug':'#9162ae','dispersed bubbles':'#a97b13'}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8*1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+'\n')


def table(frame, digits=4):
    def cell(value):
        if isinstance(value, (float,np.floating)):
            return f'{value:.{digits}f}'
        return str(value)
    lines = ['| '+' | '.join(map(str,frame.columns))+' |',
             '| '+' | '.join(['---']*len(frame.columns))+' |']
    lines += ['| '+' | '.join(cell(v) for v in row)+' |' for row in frame.itertuples(index=False,name=None)]
    return '\n'.join(lines)


def audit_inputs():
    inventory = pd.read_csv(SOURCE/'results/input_audit.csv')
    assert len(inventory) == 31 and inventory.included.sum() == 26
    assert not inventory.run.duplicated().any()
    assert set(inventory.loc[~inventory.included,'run']) == {
        'm03q01_001','m03q80_001','m04q05_001','m05q05_002','m07q10_001'}
    profiles = pd.read_csv(SOURCE/'results/profiles.csv',float_precision='round_trip')
    windows = pd.read_csv(SOURCE/'results/equal_length_windows.csv',float_precision='round_trip')
    data, checks, hashes = {}, [], {}
    paths = [SOURCE/'results'/f for f in ['input_audit.csv','profiles.csv','per_run_methods.csv',
             'equal_length_windows.csv','pair_methods.csv','analysis_config.json',
             'null_diagnostic.csv','null_replicates.csv']]
    paths += list((SOURCE/'scripts').glob('*.py'))
    paths += [HERE.parent/'estatistica'/f for f in ['descriptors_cropped.npz',
              'descriptors_center.npz','run_stats.py','extract_descriptors.py','stats_results.txt']]
    for path in paths:
        hashes[str(path.relative_to(HERE.parent))] = sha(path)
    for row in inventory.itertuples():
        image_path = HERE.parent/'data_processed'/row.run/'__panorama_v11.png'
        assert sha(image_path) == row.panorama_sha256
        hashes[str(image_path.relative_to(HERE.parent))] = row.panorama_sha256
        path = SOURCE/'intermediate'/f'{row.run}_input.npz'
        hashes[str(path.relative_to(HERE.parent))] = sha(path)
        with np.load(path,allow_pickle=False) as archive:
            data[row.run] = {k:archive[k] for k in archive.files}
        current = data[row.run]
        with Image.open(image_path) as image:
            pixels = np.asarray(image.convert('L'))
        assert pixels.shape == (row.height,row.n)
        expected = pixels[18:72].sum(axis=0,dtype=np.uint64)/54
        assert np.array_equal(current['raw'],expected)
        assert np.array_equal(current['raw_20_59'],pixels[20:60].mean(axis=0))
        assert np.array_equal(current['center'],pixels[40])
        array_hash = hashlib.sha256(expected.astype('<f8').tobytes()).hexdigest()
        assert array_hash == row.raw_array_sha256
        csv_path = SOURCE/row.signal_csv
        assert sha(csv_path) == row.signal_csv_sha256
        hashes[str(csv_path.relative_to(HERE.parent))] = row.signal_csv_sha256
        signal = pd.read_csv(csv_path,float_precision='round_trip')
        assert np.array_equal(signal.position_px,np.arange(row.n))
        assert np.array_equal(signal.gray_mean_rows_18_71,expected)
        recomputed = zero_lag_profile(prepare(expected))
        saved = profiles[(profiles.run == row.run) & (profiles.method == 'v16_biased')].sort_values('wavelength_px')
        assert len(saved) == 200
        assert np.allclose(saved.wavelength_px,WAVELENGTHS,rtol=1e-14)
        rel = np.linalg.norm(recomputed-saved.zero_lag_wps)/np.linalg.norm(saved.zero_lag_wps)
        assert rel < 1e-10, (row.run,rel)
        largest_window_error = 0.
        if row.included:
            common = 3171
            for position,start in [('prefix',0),('center',(row.n-common)//2),('suffix',row.n-common)]:
                computed = zero_lag_profile(prepare(expected[start:start+common]))
                saved_window = windows[(windows.run == row.run) & (windows.window == position)].sort_values('wavelength_px')
                assert len(saved_window) == 200 and saved_window.n.eq(common).all()
                error = np.linalg.norm(computed-saved_window.zero_lag_wps)/np.linalg.norm(saved_window.zero_lag_wps)
                assert error < 1e-10, (row.run,position,error)
                largest_window_error = max(largest_window_error,error)
        checks.append(dict(run=row.run,included=row.included,n=row.n,
            raw_max_abs_error=float(np.max(abs(current['raw']-expected))),
            full_profile_relative_l2=float(rel),window_profile_max_relative_l2=largest_window_error))
    for path in sorted((SOURCE/'intermediate').glob('*_maps.npz')):
        hashes[str(path.relative_to(HERE.parent))] = sha(path)
    pairs = pd.read_csv(SOURCE/'results/pair_methods.csv')
    map_error = verify_maps(inventory,pairs)
    # Conferir os valores-p condicionais já existentes, sem gerar novos sinais.
    null = pd.read_csv(SOURCE/'results/null_diagnostic.csv')
    replicates = pd.read_csv(SOURCE/'results/null_replicates.csv')
    for row in null.itertuples():
        samples = replicates[(replicates.run_x == row.run_x) & (replicates.run_y == row.run_y)]
        assert len(samples) == 199
        p = (1+np.count_nonzero(samples.maximum_wco >= row.observed_max))/200
        assert np.isclose(p,row.conditional_monte_carlo_p)
    frame = pd.DataFrame(checks)
    frame.to_csv(RESULTS/'source_checks_by_run.csv',index=False)
    summary = dict(status='PASS',images=31,signal_csvs=31,full_profiles_recomputed=31,
        window_profiles_recomputed=78,saved_pair_maps_equation_checked=87,
        full_profile_max_relative_l2=float(frame.full_profile_relative_l2.max()),
        window_profile_max_relative_l2=float(frame.window_profile_max_relative_l2.max()),
        saved_maps_equation_max_abs_error=map_error,
        panorama_reconstruction_reexecuted=False,full_pair_transforms_reexecuted=False,
        note='Perfis em zero recalculados independentemente; mapas salvos conferidos na Eq. 32.17. Não se refez a reconstrução dos panoramas.')
    write_json(RESULTS/'source_validation.json',summary)
    inventory['m_code'] = inventory.run.str.extract(r'm(\d+)')[0].astype(int)
    inventory['q_code'] = inventory.run.str.extract(r'q(\d+)')[0].astype(int)
    inventory.to_csv(RESULTS/'inventory.csv',index=False)
    return inventory,data,profiles,windows,hashes,summary


def profile_matrix(frame, runs, **filters):
    subset = frame
    for key,value in filters.items():
        subset = subset[subset[key] == value]
    pivot = subset.pivot(index='run',columns='wavelength_px',values='zero_lag_wps').reindex(runs)
    assert pivot.shape == (len(runs),200) and np.isfinite(pivot.to_numpy()).all()
    assert np.allclose(pivot.columns,WAVELENGTHS,rtol=1e-14)
    return pivot.to_numpy()


def extract_scenarios(primary,data,profiles,windows):
    runs = primary.run.tolist()
    common = int(primary.n.min())
    assert common == 7001
    definitions = [
        ('principal_7001_inicio','raw','prefix',common),
        ('7001_centro','raw','center',common),
        ('7001_fim','raw','suffix',common),
        ('linha40_7001_inicio','center','prefix',common),
        ('referencia20_59_7001_inicio','raw_20_59','prefix',common)]
    scenarios, coordinates = {}, []
    for name,field,position,n in definitions:
        rows = []
        for run in runs:
            raw = data[run][field]
            start = {'prefix':0,'center':(len(raw)-n)//2,'suffix':len(raw)-n}[position]
            rows.append(zero_lag_profile(prepare(raw[start:start+n])))
            coordinates.append(dict(scenario=name,run=run,field=field,start_px=start,
                                    stop_exclusive_px=start+n,n=n,n_original=len(raw)))
        scenarios[name] = np.asarray(rows)
        print('Extraído:',name,flush=True)
    scenarios['comprimento_completo'] = profile_matrix(profiles,runs,method='v16_biased')
    for position,translation in [('prefix','inicio'),('center','centro'),('suffix','fim')]:
        scenarios[f'3171_{translation}'] = profile_matrix(windows,runs,window=position)
        for row in primary.itertuples():
            start = {'prefix':0,'center':(row.n-3171)//2,'suffix':row.n-3171}[position]
            coordinates.append(dict(scenario=f'3171_{translation}',run=row.run,field='raw',
                start_px=start,stop_exclusive_px=start+3171,n=3171,n_original=row.n))
    for row in primary.itertuples():
        coordinates.append(dict(scenario='comprimento_completo',run=row.run,field='raw',
            start_px=0,stop_exclusive_px=row.n,n=row.n,n_original=row.n))
    pd.DataFrame(coordinates).to_csv(RESULTS/'window_coordinates.csv',index=False)
    for name,p in scenarios.items():
        z,coef,logp = residual_descriptor(p)
        np.savez_compressed(RESULTS/f'descriptors_{name}.npz',runs=np.array(runs),
            labels=primary.pattern.to_numpy(dtype=str),q=primary.q_code.to_numpy(),
            m=primary.m_code.to_numpy(),wavelengths_px=WAVELENGTHS,profiles=p,
            residuals=z,log_profiles=logp,coefficients=coef)
    # Conferir também um trecho experimental novo com todas as escalas da FFT.
    raw = data[runs[0]]['raw'][:common]
    x = prepare(raw)
    w,_ = transform_correlation(correlation(x,x),WAVELENGTHS)
    error = np.max(abs(scenarios['principal_7001_inicio'][0]-abs(w[:,common-1])))
    assert error < 1e-10
    return scenarios,float(error)


def analyze(primary,scenarios):
    runs = primary.run.tolist()
    labels = primary.pattern.to_numpy()
    groups = labels == 'churn'
    q, m = primary.q_code.to_numpy(),primary.m_code.to_numpy()
    p = scenarios['principal_7001_inicio']
    z, coefficients,logp = residual_descriptor(p)
    d, u = distances(z),unit_rows(z)
    observed = anosim(d,labels)
    F,r2 = permanova(u,labels)
    assignments = binary_arrangements(groups,q)
    assert len(assignments) == 200 and len(np.unique(assignments,axis=0)) == 200
    assert np.any(np.all(assignments == groups,axis=1))
    for value in np.unique(q):
        assert np.all(assignments[:,q == value].sum(axis=1) == groups[q == value].sum())
    rs = np.array([anosim(d,g) for g in assignments])
    fs = np.array([permanova(u,g)[0] for g in assignments])
    p_r,p_f = exact_upper_p(rs,observed),exact_upper_p(fs,F)
    p_adjusted = holm([p_r,p_f])
    permutation_frame = pd.DataFrame(dict(arrangement=np.arange(200),ANOSIM_R=rs,PERMANOVA_F=fs,
                                         is_observed=np.all(assignments == groups,axis=1)))
    for i,run in enumerate(runs):
        permutation_frame[run] = assignments[:,i]
    permutation_frame.to_csv(RESULTS/'exact_q_permutations.csv',index=False)
    m_assignments = binary_arrangements(groups,m)
    joint = binary_arrangements(groups,np.array([f'{mi}:{qi}' for mi,qi in zip(m,q)]))
    m_p = exact_upper_p([anosim(d,g) for g in m_assignments],observed)
    assert len(joint) == 1
    summary = dict(n_runs=len(runs),n_churn=int(groups.sum()),n_wavy_annular=int((~groups).sum()),
        n_pairs=len(runs)*(len(runs)-1)//2,n_scales=200,common_length_px=7001,
        ANOSIM_R=observed,ANOSIM_p_exact_q=p_r,ANOSIM_p_holm=float(p_adjusted[0]),
        PERMANOVA_F=F,PERMANOVA_R2_marginal=r2,PERMANOVA_p_exact_q=p_f,
        PERMANOVA_p_holm=float(p_adjusted[1]),q_arrangements=len(assignments),
        m_arrangements=len(m_assignments),joint_m_q_arrangements=len(joint),
        ANOSIM_p_exact_m_sensitivity=m_p,
        profile_values_at_or_below_log_epsilon=int(np.count_nonzero(p <= 1e-12)))
    pd.DataFrame(d,index=runs,columns=runs).to_csv(RESULTS/'pearson_distance_matrix.csv',index_label='run')
    pd.DataFrame(1-d,index=runs,columns=runs).to_csv(RESULTS/'pearson_r_matrix.csv',index_label='run')
    pairs = []
    for i,j in combinations(range(len(runs)),2):
        same = labels[i] == labels[j]
        group = 'intra' if same else ('inter_q_vizinho' if abs(q[i]-q[j]) <= 20 else 'inter_q_distante')
        pairs.append(dict(run_x=runs[i],run_y=runs[j],pattern_x=labels[i],pattern_y=labels[j],
            same_pattern=bool(same),delta_q=int(abs(q[i]-q[j])),delta_m=int(abs(m[i]-m[j])),
            category=group,pearson_r=float(1-d[i,j]),distance=float(d[i,j])))
    pairs = pd.DataFrame(pairs)
    pairs.to_csv(RESULTS/'pair_distances.csv',index=False)
    ecdf = pairs.groupby('category').distance.agg(n_pairs='size',median='median',mean='mean',
            minimum='min',maximum='max').reset_index()
    ecdf.to_csv(RESULTS/'distance_summary.csv',index=False)
    robustness = []
    for name,profile in scenarios.items():
        zr,_,_ = residual_descriptor(profile)
        dr = distances(zr)
        robustness.append(dict(scenario=name,distance='1-Pearson(residual)',
            ANOSIM_R=anosim(dr,labels),median_distance_to_primary=float(np.median(
                1-np.sum(unit_rows(zr)*u,axis=1)))))
    for name,other_d in [('Spearman_residuos',distances(z,spearman=True)),
                         ('Pearson_log_sem_remover_inclinacao',distances(logp))]:
        robustness.append(dict(scenario=name,distance=name,ANOSIM_R=anosim(other_d,labels),
                               median_distance_to_primary=np.nan))
    robustness = pd.DataFrame(robustness)
    robustness.to_csv(RESULTS/'robustness.csv',index=False)
    leave_one = []
    for i,run in enumerate(runs):
        keep = np.arange(len(runs)) != i
        leave_one.append(dict(removed_run=run,pattern=labels[i],ANOSIM_R=anosim(d[np.ix_(keep,keep)],labels[keep])))
    leave_one = pd.DataFrame(leave_one)
    leave_one.to_csv(RESULTS/'leave_one_run_out.csv',index=False)
    features = primary[['run','pattern','m_code','q_code','n','recomputed_shift_mean']].copy()
    features['peak_wavelength_px'] = WAVELENGTHS[np.argmax(p,axis=1)]
    features['log_linear_slope'] = coefficients[:,1]
    features['residual_rms'] = np.sqrt(np.mean(z*z,axis=1))
    features.to_csv(RESULTS/'features_primary.csv',index=False)
    feature_correlations = []
    feature_names = ['peak_wavelength_px','log_linear_slope','residual_rms']
    nuisance_names = ['m_code','q_code','n','recomputed_shift_mean']
    for feature in feature_names:
        for nuisance in nuisance_names:
            feature_correlations.append(dict(feature=feature,covariate=nuisance,
                spearman_rho=float(spearmanr(features[feature],features[nuisance]).statistic),
                n_runs=len(features),interpretation='descritiva; sem teste de significância'))
    pd.DataFrame(feature_correlations).to_csv(RESULTS/'feature_covariate_correlations.csv',index=False)
    features[nuisance_names].corr().to_csv(RESULTS/'covariate_pearson_correlations.csv',index_label='variable')
    pd.crosstab(primary.q_code,primary.pattern).to_csv(RESULTS/'pattern_by_q.csv')
    pd.crosstab(primary.m_code,primary.pattern).to_csv(RESULTS/'pattern_by_m.csv')
    summary['shift_length_pearson_r'] = float(np.corrcoef(features.n,features.recomputed_shift_mean)[0,1])
    summary['leave_one_out_R_min'] = float(leave_one.ANOSIM_R.min())
    summary['leave_one_out_R_max'] = float(leave_one.ANOSIM_R.max())
    return summary,z,u,d,pairs,ecdf,robustness,leave_one,features


def complementary(primary,inventory,windows,u,scenarios):
    dispersion = []
    individual = []
    for label in PRIMARY_PATTERNS:
        selected = primary.pattern.to_numpy() == label
        vectors = u[selected]
        radii = np.linalg.norm(vectors-vectors.mean(axis=0),axis=1)
        adjusted = radii*np.sqrt(len(vectors)/(len(vectors)-1))
        dispersion.append(dict(pattern=label,n=len(vectors),mean_radius=float(radii.mean()),
            mean_adjusted_radius=float(adjusted.mean()),median_adjusted_radius=float(np.median(adjusted))))
        individual.extend(dict(run=run,pattern=label,radius=float(r),adjusted_radius=float(a))
                          for run,r,a in zip(primary.loc[selected,'run'],radii,adjusted))
    dispersion = pd.DataFrame(dispersion)
    dispersion.to_csv(RESULTS/'dispersion_summary.csv',index=False)
    individual = pd.DataFrame(individual)
    individual.to_csv(RESULTS/'dispersion_by_run.csv',index=False)
    centered = u-u.mean(axis=0)
    left,singular,right = np.linalg.svd(centered,full_matrices=False)
    coordinates = left*singular
    fraction = singular**2/np.sum(singular**2)
    score_frame = primary[['run','pattern']].copy()
    for i in range(len(singular)):
        score_frame[f'PC{i+1}'] = coordinates[:,i]
    score_frame.to_csv(RESULTS/'pca_scores.csv',index=False)
    pd.DataFrame(dict(component=np.arange(1,len(fraction)+1),explained_fraction=fraction)).to_csv(RESULTS/'pca_variance.csv',index=False)
    all_runs = inventory[inventory.included].sort_values('run')
    all_p = profile_matrix(windows,all_runs.run.tolist(),window='prefix')
    all_z,coefficients,_ = residual_descriptor(all_p)
    all_d = distances(all_z)
    pd.DataFrame(all_d,index=all_runs.run,columns=all_runs.run).to_csv(RESULTS/'all26_3171_distance_matrix.csv',index_label='run')
    all_features = all_runs[['run','pattern','m_code','q_code','n']].copy()
    all_features['window_n'] = 3171
    all_features['peak_wavelength_px'] = WAVELENGTHS[np.argmax(all_p,axis=1)]
    all_features['log_linear_slope'] = coefficients[:,1]
    all_features.to_csv(RESULTS/'features_all26_3171.csv',index=False)
    all_summary = []
    for label,g in all_features.groupby('pattern'):
        idx = np.flatnonzero(all_features.pattern.to_numpy() == label)
        sub = all_d[np.ix_(idx,idx)]
        ds = sub[np.triu_indices(len(idx),1)]
        all_summary.append(dict(pattern=label,n_runs=len(g),n_intra_pairs=len(ds),
            median_intra_distance=float(np.median(ds)) if len(ds) else np.nan,
            peak_min_px=float(g.peak_wavelength_px.min()),
            peak_median_px=float(g.peak_wavelength_px.median()),peak_max_px=float(g.peak_wavelength_px.max())))
    all_summary = pd.DataFrame(all_summary)
    all_summary.to_csv(RESULTS/'all26_descriptive_summary.csv',index=False)
    # Consistência espacial: distâncias entre janelas do mesmo ensaio.
    window_rows = []
    for a,b in combinations(['principal_7001_inicio','7001_centro','7001_fim'],2):
        za,_,_ = residual_descriptor(scenarios[a])
        zb,_,_ = residual_descriptor(scenarios[b])
        distance = 1-np.sum(unit_rows(za)*unit_rows(zb),axis=1)
        for row,value in zip(primary.itertuples(),distance):
            starts = {'principal_7001_inicio':0,'7001_centro':(row.n-7001)//2,'7001_fim':row.n-7001}
            overlap = max(0,7001-abs(starts[a]-starts[b]))/7001
            window_rows.append(dict(run=row.run,pattern=row.pattern,window_a=a,window_b=b,
                distance=float(max(0,value)),overlap_fraction=overlap))
    window_rows = pd.DataFrame(window_rows)
    window_rows.to_csv(RESULTS/'within_run_window_distances.csv',index=False)
    return dispersion,individual,coordinates,fraction,all_summary,all_p,all_runs,window_rows


def check_previous(primary):
    rows = []
    for variant in ['cropped','center']:
        with np.load(HERE.parent/'estatistica'/f'descriptors_{variant}.npz',allow_pickle=False) as old:
            assert np.array_equal(old['runs'],primary.run.to_numpy())
            assert np.array_equal(old['labels'],primary.pattern.to_numpy())
            d = distances(old['Z'])
            R = anosim(d,old['labels'])
            arr = binary_arrangements(old['labels'] == 'churn',old['q'])
            p = exact_upper_p([anosim(d,g) for g in arr],R)
            expected = .1770 if variant == 'cropped' else .1791
            assert abs(R-expected) < .00005
            if variant == 'cropped':
                assert abs(p-.02) < 1e-12
            rows.append(dict(variant=variant,n_runs=len(old['runs']),ANOSIM_R=R,p_exact_q=p,
                n_q_arrangements=len(arr),scope='Descritores arquivados; extração original não reexecutada'))
    frame = pd.DataFrame(rows)
    frame.to_csv(RESULTS/'previous_analysis_rechecked.csv',index=False)
    return frame


def figures(primary,z,d,pairs,robustness,individual,coordinates,fraction,all_p,all_runs,window_rows):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes = plt.subplots(1,2,figsize=(13,5),layout='constrained')
    order = np.argsort(primary.pattern.to_numpy(),kind='stable')
    image = axes[0].imshow(d[np.ix_(order,order)],vmin=0,vmax=2,cmap='viridis')
    axes[0].set_xticks(range(19),primary.run.to_numpy()[order],rotation=90,fontsize=7)
    axes[0].set_yticks(range(19),primary.run.to_numpy()[order],fontsize=7)
    axes[0].set_title('Dissimilaridade 1 − Pearson dos resíduos')
    fig.colorbar(image,ax=axes[0],label='1 − r')
    names = {'intra':'Mesmo padrão','inter_q_vizinho':'Entre padrões, |Δq| ≤ 20',
             'inter_q_distante':'Entre padrões, |Δq| > 20'}
    for category,g in pairs.groupby('category'):
        values = np.sort(g.distance)
        axes[1].step(np.r_[0,values],np.r_[0,np.arange(1,len(values)+1)/len(values)],where='post',
                     label=f'{names[category]} (pares={len(values)})')
    axes[1].set(xlabel='Dissimilaridade 1 − Pearson',ylabel='Fração acumulada',ylim=(0,1.02),
                title='ECDF descritiva; pares compartilham ensaios')
    axes[1].legend(fontsize=9);axes[1].grid(alpha=.2)
    fig.suptitle('19 ensaios; primeiros 7.001 px; 200 comprimentos de onda')
    fig.savefig(FIGURES/'pearson_ecdf.png',dpi=170);plt.close(fig)
    fig,axes = plt.subplots(1,2,figsize=(12,5),layout='constrained')
    for k,label in enumerate(PRIMARY_PATTERNS):
        selection = primary.pattern.to_numpy() == label
        axes[0].scatter(coordinates[selection,0],coordinates[selection,1],label=label,
                        c=COLORS[label],s=55,edgecolors='white')
        values = individual.loc[individual.pattern == label,'adjusted_radius'].to_numpy()
        axes[1].scatter(np.linspace(k-.12,k+.12,len(values)),values,c=COLORS[label],s=45)
        axes[1].plot([k-.18,k+.18],[values.mean()]*2,color='black',lw=2)
    # Somente os dois pontos mais afastados na projeção recebem rótulos.
    # Todos os identificadores e coordenadas estão no CSV.
    for i in np.argsort(np.linalg.norm(coordinates[:,:2],axis=1))[-2:]:
        axes[0].annotate(primary.run.iloc[i],coordinates[i,:2],fontsize=8,
                         xytext=(4,4),textcoords='offset points')
    axes[0].set(xlabel=f'PC1 ({fraction[0]:.1%})',ylabel=f'PC2 ({fraction[1]:.1%})',
                title='PCA dos descritores centrados e normalizados')
    axes[0].legend();axes[0].grid(alpha=.2)
    axes[1].set_xticks([0,1],PRIMARY_PATTERNS)
    axes[1].set(ylabel='Distância euclidiana ao centroide, ajustada',
                title='Dispersão descritiva; traço = média')
    axes[1].grid(axis='y',alpha=.2)
    fig.savefig(FIGURES/'pca_dispersion.png',dpi=170);plt.close(fig)
    fig,ax = plt.subplots(figsize=(10,6),layout='constrained')
    rows = robustness.iloc[::-1]
    ax.scatter(rows.ANOSIM_R,np.arange(len(rows)),color='#2475b1',s=50)
    ax.set_yticks(np.arange(len(rows)),rows.scenario)
    ax.axvline(0,color='black',lw=.8);ax.grid(axis='x',alpha=.2)
    ax.set(xlabel='ANOSIM R',title='Sensibilidade: mesma lista de 19 ensaios; sem selecionar o maior R')
    fig.savefig(FIGURES/'robustness.png',dpi=170);plt.close(fig)
    fig,axes = plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    for label in PRIMARY_PATTERNS:
        selected = primary.pattern.to_numpy() == label
        for row in z[selected]:
            axes[0].plot(WAVELENGTHS,row,c=COLORS[label],alpha=.16,lw=.8)
        axes[0].plot(WAVELENGTHS,z[selected].mean(axis=0),c=COLORS[label],label=label,lw=2)
        values = np.sort(window_rows.loc[window_rows.pattern == label,'distance'])
        axes[1].step(np.r_[0,values],np.r_[0,np.arange(1,len(values)+1)/len(values)],where='post',label=label)
    axes[0].set(xscale='log',xlabel='Comprimento de onda (px)',ylabel='Resíduo de log alternative WPS',
                title='Perfis individuais e média por padrão')
    axes[1].set(xlabel='1 − Pearson entre janelas do mesmo ensaio',ylabel='Fração acumulada',
                title='Consistência espacial; janelas podem se sobrepor')
    for ax in axes:ax.legend();ax.grid(alpha=.2)
    fig.savefig(FIGURES/'profiles_windows.png',dpi=170);plt.close(fig)
    fig,axes = plt.subplots(1,5,figsize=(17,4),layout='constrained',sharey=True)
    normalized = all_p/all_p.max(axis=1,keepdims=True)
    for ax,(label,group) in zip(axes,all_runs.groupby('pattern')):
        selected = all_runs.pattern.to_numpy() == label
        for row in normalized[selected]:
            ax.plot(WAVELENGTHS,row,c=COLORS[label],alpha=.4,lw=.8)
        ax.plot(WAVELENGTHS,np.median(normalized[selected],axis=0),c=COLORS[label],lw=2)
        ax.set(xscale='log',xlabel='Comprimento de onda (px)',title=f'{label} (n={len(group)})')
        ax.grid(alpha=.2)
    axes[0].set_ylabel('alternative WPS / máximo de cada ensaio')
    fig.suptitle('26 ensaios, primeiros 3.171 px; traço grosso = mediana; sem teste entre cinco classes')
    fig.savefig(FIGURES/'all26_profiles.png',dpi=170);plt.close(fig)


def verify_outputs(summary):
    """Cruza os arquivos finais com fórmulas independentes e verifica os links."""
    for script in HERE.glob('*.py'):
        ast.parse(script.read_text())
    with np.load(RESULTS/'descriptors_principal_7001_inicio.npz') as archive:
        z = archive['residuals']
        groups = archive['labels'] == 'churn'
        expected = 1-np.corrcoef(z)
    saved = pd.read_csv(RESULTS/'pearson_distance_matrix.csv',index_col=0).to_numpy()
    assert np.allclose(saved,expected,atol=1e-12,rtol=1e-12)
    iu = np.triu_indices(len(groups),1)
    same = groups[iu[0]] == groups[iu[1]]
    intra,inter = saved[iu][same],saved[iu][~same]
    auc = np.mean((inter[:,None] > intra[None,:]) + .5*(inter[:,None] == intra[None,:]))
    assert np.isclose(2*auc-1,summary['ANOSIM_R'],atol=1e-12)
    permutations = pd.read_csv(RESULTS/'exact_q_permutations.csv')
    assert len(permutations) == 200 and permutations.is_observed.sum() == 1
    for statistic,key in [('ANOSIM_R','ANOSIM_p_exact_q'),('PERMANOVA_F','PERMANOVA_p_exact_q')]:
        p = (permutations[statistic] >= summary[statistic]-1e-12).mean()
        assert np.isclose(p,summary[key],atol=1e-12)
    broken = []
    for document in HERE.glob('*.md'):
        for target in re.findall(r'[]][(]([^)]+)[)]',document.read_text()):
            if '://' not in target and not target.startswith('#'):
                path = document.parent/target.split('#')[0]
                if not path.exists() and path != RESULTS/'manifest.json':
                    broken.append(str(path))
    assert not broken,broken
    figures = list(FIGURES.glob('*.png'))
    assert len(figures) == 5
    for path in figures:
        with Image.open(path) as image:
            image.verify()
    write_json(RESULTS/'final_verification.json',dict(status='PASS',
        independent_pearson_and_anosim_auc_check=True,exported_exact_p_values_checked=True,
        figures_verified=len(figures),broken_local_links=broken,
        observed_arrangement_count=int(permutations.is_observed.sum()),
        experimental_pair_distances=len(intra)+len(inter)))


def main():
    RESULTS.mkdir(exist_ok=True);FIGURES.mkdir(exist_ok=True)
    formula_checks = validate()
    write_json(RESULTS/'formula_validation.json',formula_checks)
    print('Fórmulas verificadas. Auditando dados e perfis salvos...',flush=True)
    inventory,data,profiles,windows,hashes,source_checks = audit_inputs()
    print('Entradas, 109 perfis e 87 mapas verificados.',flush=True)
    primary = inventory[inventory.included & inventory.pattern.isin(PRIMARY_PATTERNS)].sort_values('run').reset_index(drop=True)
    assert len(primary) == 19
    previous = check_previous(primary)
    scenarios,experimental_error = extract_scenarios(primary,data,profiles,windows)
    summary,z,u,d,pairs,ecdf,robustness,leave_one,features = analyze(primary,scenarios)
    dispersion,individual,coordinates,fraction,all_summary,all_p,all_runs,window_rows = complementary(primary,inventory,windows,u,scenarios)
    summary['PCA_two_components_fraction'] = float(fraction[:2].sum())
    summary['experimental_direct_profile_vs_fft_max_abs_error'] = experimental_error
    summary['timestamp_utc'] = datetime.now(timezone.utc).isoformat()
    write_json(RESULTS/'summary.json',summary)
    figures(primary,z,d,pairs,robustness,individual,coordinates,fraction,all_p,all_runs,window_rows)
    from report_results import generate_report
    generate_report()
    verify_outputs(summary)
    # Nenhuma dependência consultada pode ter sido alterada durante a execução.
    for name,value in hashes.items():
        assert sha(HERE.parent/name) == value, f'Fonte alterada durante a execução: {name}'
    write_json(RESULTS/'manifest.json',dict(timestamp_utc=summary['timestamp_utc'],python=sys.version,
        packages={name:metadata.version(name) for name in ['numpy','scipy','pandas','matplotlib','PyWavelets','Pillow']},
        module_versions=dict(numpy=np.__version__,pandas=pd.__version__),
        seed_formula_controls=20261024,experimental_inference='enumeração exata; sem aleatoriedade',
        sources_unchanged=True,source_sha256=hashes,
        scripts_sha256={p.name:sha(p) for p in HERE.glob('*.py')},
        method_plan_sha256=sha(HERE/'METODOLOGIA.md'),
        generated_sha256={str(p.relative_to(HERE)):sha(p)
                          for p in list(RESULTS.iterdir())+list(FIGURES.iterdir())+list(HERE.glob('*.md'))
                          if p.is_file() and p.name != 'manifest.json'}))
    print(json.dumps(summary,indent=2,ensure_ascii=False))


if __name__ == '__main__':
    main()
