#!/usr/bin/env python3
"""Auditoria independente v14.1. Fontes somente leitura; saídas apenas em v17.

Não treina nem executa classificadores. Reconstrói entradas a partir dos JPEGs,
confere aritmética independente e testa dependência de emendas/ocupação/ROI.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'src2/v17/resultados/voidfraction'
OUT.mkdir(parents=True, exist_ok=True)
os.environ['MPLCONFIGDIR'] = str(OUT / 'cache')
os.environ['MPLBACKEND'] = 'Agg'
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / 'src2/v14.1_new_voidfraction/scripts'))

import numpy as np
import pandas as pd
from PIL import Image
from scipy.signal import correlate
from scipy.stats import mannwhitneyu, spearmanr

BASE = ROOT / 'src2/v14.1_new_voidfraction'
DATA = ROOT / 'src2/data_processed'


def records(path):
    return list(csv.DictReader(path.open()))


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def save(name, rows):
    pd.DataFrame(rows).to_csv(OUT / name, index=False)


def read_image(path):
    with Image.open(path) as im:
        return np.asarray(im.convert('L'))


def disk_area(mask):
    """Conta intervalos por diferenças nas bordas; não usa geometry()."""
    # Resultado por coluna por soma de comprimentos quadrados / D^2.
    edges = np.diff(np.pad(mask.astype(np.int8), ((1, 1), (0, 0))), axis=0)
    start_y, start_x = np.where(edges == 1)
    end_y, end_x = np.where(edges == -1)
    start_order = np.lexsort((start_y, start_x))
    end_order = np.lexsort((end_y, end_x))
    assert np.array_equal(start_x[start_order], end_x[end_order])
    widths = end_y[end_order] - start_y[start_order]
    return np.bincount(start_x[start_order], weights=widths.astype(float)**2,
                       minlength=mask.shape[1]) / 53**2


def flip_statistics(mask, seams):
    m = mask[18:71]
    change = m[:, 1:] != m[:, :-1]
    count = int(change.sum())
    row_p = m.mean(axis=1)
    null = float(np.mean(2*row_p*(1-row_p)*m.shape[1]/(m.shape[1]-1)))
    seam_flags = np.zeros(m.shape[1]-1, dtype=bool)
    seam_flags[seams-1] = True
    seam_count = int(change[:, seam_flags].sum())
    inside_count = count-seam_count
    seam_n = int(change[:, seam_flags].size)
    inside_n = int(change[:, ~seam_flags].size)
    return dict(flip_rate=count/change.size, flip_count=count,
                seam_flip_rate=seam_count/seam_n,
                nonseam_flip_rate=inside_count/inside_n,
                seam_flip_share=seam_count/count if count else 0,
                seam_comparison_share=seam_n/change.size,
                rowwise_shuffle_expected_flip=null,
                flip_over_rowwise_shuffle=(count/change.size)/null if null else np.nan,
                projected_fraction=float(m.mean()),
                full_lumen_column_fraction=float(m.all(axis=0).mean()),
                row40_flip=float(np.mean(mask[40, 1:] != mask[40, :-1])),
                row44_flip=float(np.mean(mask[44, 1:] != mask[44, :-1])))


def measure_lags(profiles):
    """Recalcula Pearson de sobreposição, momentos e produto via FFT."""
    n = profiles.shape[1]
    lag = np.arange(-min(200,n//4),min(200,n//4)+1)
    count = n-np.abs(lag)
    la, lb = np.maximum(lag,0), np.maximum(-lag,0)
    shifts = []
    margins = []
    for a,b in zip(profiles[:-1],profiles[1:]):
        pa, pb = np.r_[0,np.cumsum(a)], np.r_[0,np.cumsum(b)]
        qa, qb = np.r_[0,np.cumsum(a*a)],np.r_[0,np.cumsum(b*b)]
        sa,sb=pa[la+count]-pa[la],pb[lb+count]-pb[lb]
        va=qa[la+count]-qa[la]-sa*sa/count
        vb=qb[lb+count]-qb[lb]-sb*sb/count
        numerator=correlate(a,b,mode='full',method='fft')[n-1+lag]-sa*sb/count
        denominator=np.sqrt(np.maximum(va,0)*np.maximum(vb,0))
        score=np.divide(numerator,denominator,out=np.full(len(lag),-np.inf),where=denominator>0)
        idx=int(np.argmax(score))
        shifts.append(int(lag[idx]))
        margins.append(float(score[idx]-np.partition(score,-2)[-2]))
    return np.array(shifts), np.array(margins)


def exact_direct_lag(a,b):
    scores=[]
    lim=min(200,len(a)//4)
    for k in range(-lim,lim+1):
        av,bv=(a[k:],b[:len(b)-k]) if k>=0 else (a[:len(a)+k],b[-k:])
        scores.append(np.corrcoef(av,bv)[0,1])
    return int(np.nanargmax(scores)-lim)


def audit_raw_runs(shifts_table, perrun, summary):
    rows=[]
    for run in perrun.index:
        sub=shifts_table[shifts_table.run==run]
        files=sorted((DATA/run).glob('ImgA*.jpg'))
        frames=np.stack([read_image(p) for p in files])
        assert list(sub.frame)==[p.name for p in files]
        profiles=frames[:,20:60].mean(axis=1)
        lag,margin=measure_lags(profiles)
        widths=np.maximum(np.r_[lag,lag[-1]],1)
        direct_ok=all(exact_direct_lag(profiles[k],profiles[k+1])==lag[k]
                      for k in [0,len(lag)//2,len(lag)-1])
        same=np.array_equal(lag,sub.raw_shift_px.to_numpy()[:-1])
        strips=[frame[:,frames.shape[2]//2-int(w)//2:frames.shape[2]//2-int(w)//2+int(w)]
                for frame,w in zip(frames,widths)]
        reconstructed=np.concatenate(strips,axis=1)
        saved=read_image(BASE/'intermediate/reconstructed_inputs'/run/'__panorama_v11.png')
        cfg={}
        for line in (DATA/run/'_config.xsv').read_text(encoding='latin-1').splitlines():
            if '=' in line:
                k,v=line.split('=',1);cfg[k.strip()]=v.strip()
        native=perrun.loc[run]
        probes=np.array([frames.shape[2]//4,3*frames.shape[2]//8,frames.shape[2]//2,
                         5*frames.shape[2]//8,3*frames.shape[2]//4])
        # Otsu limiar salvo: aritmética da taxa independe da implementação original.
        state=frames[:,20:60,:][:,:,probes]<=int(native.raw_otsu_threshold)
        flips=int(np.count_nonzero(state[1:]!=state[:-1]))
        elapsed=(len(files)-1)*int(cfg['Step'])/float(cfg['Rate'])
        # Mesma janela de 70 ms e mesma separação de 1 ms entre comparações.
        # Amostragem de estados (sem filtro) mede diretamente transições nessa
        # escala; não estima todas as transições que ocorreram entre quadros.
        sample_rate=float(cfg['Rate'])/int(cfg['Step'])
        indices=np.rint(np.arange(71)*sample_rate/1000).astype(int)
        assert np.max(np.abs(indices/sample_rate-np.arange(71)/1000))<1e-12
        common=state[indices]
        common_count=int(np.count_nonzero(common[1:]!=common[:-1]))
        common_rate=common_count/common[1:].size
        hist=np.bincount(frames[:,20:60,:].ravel(),minlength=256).astype(float)
        mass=np.cumsum(hist);moment=np.cumsum(hist*np.arange(256))
        den=mass*(mass[-1]-mass)
        between=np.divide((moment[-1]*mass-moment*mass[-1])**2,den,
                          out=np.full(256,-np.inf),where=den>0)
        independent_threshold=int(np.argmax(between))
        rows.append(dict(run=run,n_frames=len(files),all_raw_lags_match_saved=bool(same),
                         direct_Pearson_three_lags_match=bool(direct_ok),
                         corrected_pixels_exact=bool(np.array_equal(saved,reconstructed)),
                         nonpositive_lags=int((lag<=0).sum()),
                         lags_at_search_limit=int((np.abs(lag)==200).sum()),
                         best_second_margin_median=float(np.median(margin)),
                         native_5probe_count_matches=bool(flips==native.raw_5probe_flip_count),
                         native_otsu_independent_matches=bool(independent_threshold==native.raw_otsu_threshold),
                         elapsed_s=elapsed,native_elapsed_error=abs(elapsed-native.elapsed_s),
                         native_5probe_rate_error=abs(flips/state[1:].size-native.raw_5probe_flip_rate),
                         pattern=native.pattern,
                         eligible=bool(native.usable_existing_20pct_saturation_rule),
                         native_full_hz=float(native.raw_5probe_flip_hz_per_pixel),
                         common_1ms_flip_rate=common_rate,
                         common_1ms_switches_per_pixel_s=common_rate*1000))
        print('raw',run,'pixels',rows[-1]['corrected_pixels_exact'],'lags',same,flush=True)
    save('raw_reconstruction_checks.csv',rows)
    rawdf=pd.DataFrame(rows)
    comps=[]
    for scope,sample in [('all26',rawdf),('eligible25',rawdf[rawdf.eligible])]:
        for feature in ['native_full_hz','common_1ms_switches_per_pixel_s']:
            c=sample.loc[sample.pattern=='churn',feature].to_numpy()
            a=sample.loc[sample.pattern=='wavy annular',feature].to_numpy()
            u=mannwhitneyu(c,a,method='auto',alternative='two-sided')
            comps.append(dict(scope=scope,feature=feature,n_churn=len(c),n_annular=len(a),
                              churn_median=float(np.median(c)),annular_median=float(np.median(a)),
                              auc_churn_gt_annular=float(u.statistic/(len(c)*len(a))),
                              exploratory_uncorrected_p=float(u.pvalue)))
    save('native_common_time_diagnostics.csv',comps)
    summary['raw_reconstruction']={
        'runs':len(rows),'frames':sum(r['n_frames'] for r in rows),
        'all_corrected_pixels_exact':all(r['corrected_pixels_exact'] for r in rows),
        'all_raw_lags_exact':all(r['all_raw_lags_match_saved'] for r in rows),
        'all_direct_subset_lags_exact':all(r['direct_Pearson_three_lags_match'] for r in rows),
        'all_native_5probe_counts_exact':all(r['native_5probe_count_matches'] for r in rows),
        'all_native_otsu_thresholds_exact':all(r['native_otsu_independent_matches'] for r in rows),
        'nonpositive_lags':sum(r['nonpositive_lags'] for r in rows)}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--skip-raw',action='store_true',help='Só para depuração; auditoria completa não usa esta opção.')
    args=parser.parse_args()
    estimator=importlib.import_module('panorama_void_alternatives')
    void=pd.read_csv(BASE/'results/void_fraction/panorama_void_per_run.csv').set_index('run')
    phase=pd.read_csv(BASE/'results/flip_rate/phase_mask_flip_per_run.csv').set_index('run')
    native=pd.read_csv(BASE/'results/flip_rate/flip_rate_per_run.csv').set_index('run')
    shifts=pd.read_csv(BASE/'results/validation/reconstruction_shifts.csv')
    source_hashes=json.loads((BASE/'results/validation/source_input_hashes.json').read_text())
    mismatch=[]
    for name,expected in source_hashes.items():
        p=ROOT/'src2'/name
        actual=sha(p) if p.exists() else 'MISSING'
        if actual!=expected:mismatch.append(dict(path=str(p.relative_to(ROOT)),expected=expected,actual=actual))
    (OUT/'source_hash_drift.json').write_text(json.dumps(mismatch,indent=2)+'\n')
    summary={'source_hashes_checked':len(source_hashes),'source_hash_mismatches':len(mismatch),
             'software':{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__},
             'scope':'26 aquisições; 25 elegíveis no corte de saturação herdado; sem aprendizado de máquina.'}
    ref,selections=estimator.make_reference()
    saved_ref=np.array(json.loads((BASE/'intermediate/features/panorama_liquid_reference.json').read_text())['radial_reference'])
    summary['reference_max_error']=float(np.max(np.abs(ref-saved_ref)))
    checks=[]
    robustness=[]
    monotone=[]
    for run,v in void.iterrows():
        im=estimator.read(run)
        mask,accepted,_=estimator.silhouette(im,ref)
        saved_mask=read_image(BASE/'intermediate/phase_masks'/(run+'_mask.png'))>0
        reversed_mask,_,_=estimator.silhouette(im[:,::-1],ref)
        area=disk_area(saved_mask)
        sub=shifts[shifts.run==run]
        boundaries=np.cumsum(sub.strip_width_px.to_numpy())[:-1]
        metrics=flip_statistics(saved_mask,boundaries)
        errors=dict(run=run,pattern=phase.loc[run,'pattern'],eligible=bool(phase.loc[run,'eligible']),
                    saved_mask_exact=bool(np.array_equal(mask,saved_mask)),
                    axial_reverse_mask_exact=bool(np.array_equal(mask,reversed_mask[:,::-1])),
                    circle_alpha_error=float(abs(area.mean()-v.conditional_circle_alpha_spatial_mean)),
                    projected_fraction_error=abs(metrics['projected_fraction']-v.gas_silhouette_projected_fraction),
                    flip_rate_error=abs(metrics['flip_rate']-phase.loc[run,'candidate_pixel_flip_rate']),
                    alpha=float(area.mean()),saturation=float(v.core_saturated_fraction),
                    boundary_touch_fraction=float(np.mean(mask[18]|mask[70])),
                    **metrics)
        checks.append(errors)
        # Inversão radial simultânea da imagem e referência no lúmen; diagnóstico
        # de dependência de janelas fixas, não prova de erro físico (óptica assimétrica).
        mirrored=im.copy();mirrored[18:71]=im[18:71][::-1]
        mirrorref=ref.copy();mirrorref[18:71]=ref[18:71][::-1]
        mm,_,_=estimator.silhouette(mirrored,mirrorref)
        mirror_area=float(disk_area(mm).mean())
        robust=dict(run=run,pattern=errors['pattern'],eligible=errors['eligible'],
                    original_alpha=errors['alpha'],mirrored_alpha=mirror_area,
                    alpha_mirror_delta=mirror_area-errors['alpha'],
                    row40_flip=metrics['row40_flip'],row44_flip=metrics['row44_flip'],
                    alpha_threshold_range=v.threshold_060_080_max-v.threshold_060_080_min,
                    alpha_reference_range=v.reference_variants_max-v.reference_variants_min,
                    alpha_five_blocks_range=v.five_spatial_block_means_max-v.five_spatial_block_means_min)
        robustness.append(robust)
        previous=None
        for t in [.60,.65,.70,.75,.80]:
            cur,_,_=estimator.silhouette(im,ref,threshold=t)
            if previous is not None:
                removed=previous & ~cur
                monotone.append(dict(run=run,threshold_to=t,removed_pixels=int(removed.sum()),
                                     added_pixels=int((cur & ~previous).sum()),
                                     alpha_delta=float(disk_area(cur).mean()-disk_area(previous).mean())))
            previous=cur
        print('mask',run,'exact',errors['saved_mask_exact'],flush=True)
    save('arithmetic_and_seam_checks.csv',checks)
    save('robustness_probes.csv',robustness)
    save('threshold_monotonicity.csv',monotone)
    df=pd.DataFrame(checks)
    comps=[]
    for scope,sample in [('all26',df),('eligible25',df[df.eligible])]:
        for feature in ['alpha','flip_rate','nonseam_flip_rate','flip_over_rowwise_shuffle','row40_flip','row44_flip']:
            c=sample.loc[sample.pattern=='churn',feature].dropna().to_numpy()
            a=sample.loc[sample.pattern=='wavy annular',feature].dropna().to_numpy()
            test=mannwhitneyu(c,a,alternative='two-sided',method='auto')
            comps.append(dict(scope=scope,feature=feature,n_churn=len(c),n_annular=len(a),
                              churn_median=float(np.median(c)),annular_median=float(np.median(a)),
                              auc_churn_gt_annular=float(test.statistic/(len(c)*len(a))),
                              exploratory_uncorrected_p=float(test.pvalue)))
    save('non_ml_group_diagnostics.csv',comps)
    summary['mask_reproduction_all_exact']=all(r['saved_mask_exact'] for r in checks)
    summary['mask_axial_reverse_all_exact']=all(r['axial_reverse_mask_exact'] for r in checks)
    summary['arithmetic_max_abs_error']={key:max(r[key] for r in checks)
        for key in ['circle_alpha_error','projected_fraction_error','flip_rate_error']}
    summary['threshold_nonmonotone']={'comparisons':len(monotone),
        'comparisons_removing_pixels':sum(r['removed_pixels']>0 for r in monotone),
        'comparisons_lowering_alpha':sum(r['alpha_delta']<-1e-12 for r in monotone)}
    summary['seam_flip_share_range']=[float(df.seam_flip_share.min()),float(df.seam_flip_share.max())]
    summary['seam_flip_share_median']=float(df.seam_flip_share.median())
    summary['seam_comparison_share_median']=float(df.seam_comparison_share.median())
    summary['alpha_flip_spearman_all26']=float(spearmanr(df.alpha,df.flip_rate).statistic)
    if not args.skip_raw:audit_raw_runs(shifts,native,summary)
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps(summary,indent=2,allow_nan=False))


if __name__=='__main__':main()
