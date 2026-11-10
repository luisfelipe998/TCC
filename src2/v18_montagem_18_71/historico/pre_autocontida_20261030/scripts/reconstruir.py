#!/usr/bin/env python3
"""Rebuild the 26 retained panoramas with alignment profiles from rows 18:72.

Only the alignment band changes. The wavelet image keeps the historical v11
seam operation. The optical image concatenates the same strips without it.
The final manifest is published only after all inputs and outputs are checked.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from io import BytesIO
import json
import os
from pathlib import Path
import platform
import time
import warnings

import numpy as np
from PIL import Image
from scipy.fft import next_fast_len, rfft, irfft

PROJECT = Path(__file__).resolve().parents[3]
VERSION = Path(__file__).resolve().parents[1]
OUTPUT = VERSION / '01_panoramas'
DATA = PROJECT / 'src2/data_processed'
V16 = PROJECT / 'src2/v16_spatial2_18_71'
V14 = PROJECT / 'src2/v14.1_new_voidfraction'
ROWS = slice(18, 72)
NEAR_MAXIMUM_TOLERANCE = 5e-11


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def relative(path):
    return path.relative_to(PROJECT).as_posix()


def save_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def literal_pearson(a, b):
    """Original ascending-lag loop, including undefined values and first winner."""
    limit = min(200, len(a)//4)
    best_lag, best_value = 0, -np.inf
    scores = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        for lag in range(-limit, limit+1):
            aa = a[max(lag, 0):len(a)+min(lag, 0)]
            bb = b[max(-lag, 0):len(b)-max(lag, 0)]
            value = float(np.corrcoef(aa, bb)[0, 1])
            scores.append(value)
            if value > best_value:
                best_lag, best_value = lag, value
    finite = np.asarray(scores)[np.isfinite(scores)]
    tied = int(np.count_nonzero(finite == best_value)) if len(finite) else 0
    return best_lag, best_value, tied


def vectorized_pearson(profiles):
    """The same overlap Pearson objective, with stable centered FFT products.

    A full-row constant is removed only to improve arithmetic. Every overlap
    is still centered and normalized separately. Ambiguous numerical maxima
    fall back to the original literal search, preserving its tie behavior.
    """
    profiles = np.asarray(profiles, dtype=np.float64)
    centered = profiles - profiles.mean(axis=1, keepdims=True)
    a, b = centered[:-1], centered[1:]
    width = profiles.shape[1]
    limit = min(200, width//4)
    lags = np.arange(-limit, limit+1)
    overlap = width - abs(lags)
    transform = rfft(centered, n=next_fast_len(2*width-1), axis=1)
    products = irfft(transform[:-1]*transform[1:].conj(), n=next_fast_len(2*width-1), axis=1)
    products = products[:, lags % products.shape[1]]
    a_start, a_stop = np.maximum(lags, 0), width+np.minimum(lags, 0)
    b_start, b_stop = np.maximum(-lags, 0), width-np.maximum(lags, 0)

    def moments(x, start, stop):
        sums = np.c_[np.zeros(len(x)), np.cumsum(x, axis=1)]
        squares = np.c_[np.zeros(len(x)), np.cumsum(x*x, axis=1)]
        return sums[:,stop]-sums[:,start], squares[:,stop]-squares[:,start]

    sa, qa = moments(a, a_start, a_stop)
    sb, qb = moments(b, b_start, b_stop)
    denominator = np.sqrt(np.maximum(qa-sa*sa/overlap, 0)*np.maximum(qb-sb*sb/overlap, 0))
    scores = np.divide(products-sa*sb/overlap, denominator,
                       out=np.full_like(products, -np.inf), where=denominator>0)
    scores[np.isfinite(scores)] = np.clip(scores[np.isfinite(scores)], -1, 1)
    indices = np.argmax(scores, axis=1)
    chosen = lags[indices].copy()
    best = scores[np.arange(len(indices)), indices].copy()
    sorted_scores = np.sort(scores, axis=1)
    gap = np.full(len(best), np.inf)
    both = np.isfinite(sorted_scores[:,-2:]).all(axis=1)
    gap[both] = sorted_scores[both,-1]-sorted_scores[both,-2]
    fallback = (~np.isfinite(best)) | (gap <= NEAR_MAXIMUM_TOLERANCE)
    details = []
    for index in np.flatnonzero(fallback):
        lag, value, tied = literal_pearson(profiles[index], profiles[index+1])
        chosen[index], best[index] = lag, value
        details.append(dict(pair=int(index), exact_tied_maxima=tied,
                            undefined_all_lags=not np.isfinite(value)))
    return chosen, best, gap, details


def synthetic_checks():
    rng = np.random.default_rng(1871)
    a = rng.normal(size=1016)
    plus, minus = rng.normal(size=(2,1016))
    plus[:-200] = a[200:]
    minus[200:] = a[:-200]
    periodic = np.tile([0., 1.], 508)
    cases = [('positive_limit',a,plus,200), ('negative_limit',a,minus,-200),
             ('exact_periodic_ties',periodic,periodic,None),
             ('undefined_constant',np.ones(1016),np.ones(1016),0)]
    report=[]
    for name,a,b,expected in cases:
        found,_,_,fallback=vectorized_pearson(np.stack([a,b]))
        literal,value,ties=literal_pearson(a,b)
        assert found[0] == literal, name
        if expected is not None:
            assert literal == expected, name
        report.append(dict(case=name,lag=literal,literal_agrees=True,
                           exact_tied_maxima=ties, fallback_pairs=len(fallback),
                           best_score=float(value) if np.isfinite(value) else None))
    return report


def make_images(frames, widths):
    center=frames.shape[2]//2
    strips=[frame[:,center-int(width)//2:center-int(width)//2+int(width)]
            for frame,width in zip(frames,widths)]
    clean=np.concatenate(strips,axis=1).copy()
    historical=clean.astype(np.float32)
    seams=np.cumsum(widths)[:-1]
    for index,seam in enumerate(seams,1):
        feather=min(2,min(int(widths[index-1]),int(widths[index]))//2)
        if feather:
            weight=np.linspace(0.,1.,feather)[None,:]
            historical[:,seam-feather:seam] = (
                strips[index-1][:,-feather:]*(1-weight)+strips[index][:,:feather]*weight)
    return np.clip(historical,0,255).astype(np.uint8),clean


def reverse_collage(frames,widths,wavelet,optical):
    """Read every saved strip back against its original frame, including seams."""
    cursor=0
    wavelet_error=optical_error=0
    center=frames.shape[2]//2
    for index,(frame,width) in enumerate(zip(frames,widths)):
        width=int(width)
        source=frame[:,center-width//2:center-width//2+width]
        optical_error+=int(np.count_nonzero(optical[:,cursor:cursor+width]!=source))
        expected=source.copy()
        # With the historical maximum feather=2, the only changed column is
        # the last column before a seam when both adjacent strips are >=4 px.
        if index+1<len(widths) and min(width,int(widths[index+1]))>=4:
            next_start=center-int(widths[index+1])//2
            expected[:,-1]=frames[index+1,:,next_start+1]
        wavelet_error+=int(np.count_nonzero(wavelet[:,cursor:cursor+width]!=expected))
        cursor+=width
    assert cursor==wavelet.shape[1]==optical.shape[1]
    assert wavelet_error==optical_error==0
    return dict(optical_wrong_pixels=optical_error,wavelet_wrong_pixels=wavelet_error,
                checked_pixels_each=int(wavelet.size),all_strips_checked=len(widths))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--check-only',action='store_true',help='Only execute synthetic arithmetic checks.')
    args=parser.parse_args()
    synthetic=synthetic_checks()
    if args.check_only:
        print(json.dumps(synthetic,indent=2));return
    started=time.monotonic()
    final_manifest=OUTPUT/'manifesto.json'
    if final_manifest.exists():
        raise RuntimeError('Final manifest already exists. Do not overwrite a published reconstruction.')
    OUTPUT.mkdir(parents=True,exist_ok=True)
    selection_path=V16/'intermediario/aquisicoes.json'
    selected=json.loads(selection_path.read_text())
    assert len(selected)==26 and [r['numero'] for r in selected]==list(range(1,27))
    before={}

    def preserve(path):
        digest=sha(path)
        before[relative(path)]=digest
        return digest

    for path in [selection_path,Path(__file__),PROJECT/'src2/pipeline_v11.py',
                 V16/'scripts/auditoria_adversarial.py',V16/'scripts/common.py',
                 V14/'scripts/reconstruct_inputs.py',
                 V14/'results/validation/input_reconstruction.json',
                 V14/'results/validation/reconstruction_shifts.csv',
                 V16/'intermediario/auditoria/manifesto_artefatos.json',
                 PROJECT/'src2/v14.2_54_linhas/results/validation/artifact_manifest.json']:
        preserve(path)
    records, shifts_table, audits=[],[],[]
    frame_total=literal_total=0
    for old in selected:
        run=old['ensaio'];folder=DATA/run
        markers=list(folder.glob('__*.txt'))
        assert len(markers)==1
        pattern=markers[0].stem[2:].replace('_',' ')
        expected='dispersed bubbles' if old['padrao']=='dispersed bubble' else old['padrao']
        assert pattern==expected,run
        preserve(markers[0])
        files=sorted(folder.glob('ImgA*.jpg'))
        assert [p.name for p in files]==[f'ImgA{i:06d}.jpg' for i in range(501)]
        frames=[]
        for path in files:
            payload=path.read_bytes()
            before[relative(path)]=hashlib.sha256(payload).hexdigest()
            with Image.open(BytesIO(payload)) as image:
                frames.append(np.asarray(image.convert('L')))
        frames=np.stack(frames)
        assert frames.shape==(501,80,1016)
        profiles=frames[:,ROWS].mean(axis=1)
        raw,scores,gaps,fallback=vectorized_pearson(profiles)
        sample=set([0,250,499])
        sample.update(np.argsort(gaps)[:3].tolist())
        sample.update(np.flatnonzero((raw<=0)|(abs(raw)==200)).tolist())
        sample.update(item['pair'] for item in fallback)
        literal_checks=[]
        for index in sorted(sample):
            exact,value,ties=literal_pearson(profiles[index],profiles[index+1])
            assert raw[index]==exact,(run,index,int(raw[index]),exact)
            error=abs(scores[index]-value) if np.isfinite(value) else 0.
            assert error<=1e-10,(run,index,error)
            literal_checks.append(dict(pair=index,lag=exact,score_error=float(error),
                                       exact_tied_maxima=ties))
        widths=np.maximum(np.r_[raw,raw[-1]],1)
        old_cache=V16/'_OLD/intermediate'/f'{run}_input.npz'
        old_wavelet=folder/'__panorama_v11.png'
        old_optical=V14/'intermediate/reconstructed_inputs'/run/'__panorama_v11.png'
        prior_hashes={name:preserve(path) for name,path in
                      [('wavelet',old_wavelet),('optica',old_optical),('shift_cache',old_cache)]}
        assert prior_hashes['wavelet']==old['panorama_sha256']
        with np.load(old_cache) as cache:
            old_raw=cache['raw_shifts'].copy();old_widths=cache['shifts'].copy()
        assert old_raw.shape==raw.shape and old_widths.shape==widths.shape
        wavelet,optical=make_images(frames,widths)
        destinations={}
        output_hashes={}
        loaded={}
        for variant,array in [('wavelet',wavelet),('optica',optical)]:
            target=OUTPUT/variant/run/'panorama.png'
            target.parent.mkdir(parents=True,exist_ok=True)
            Image.fromarray(array).save(target)
            with Image.open(target) as image:
                loaded[variant]=np.asarray(image).copy()
            np.testing.assert_array_equal(loaded[variant],array)
            destinations[variant]=relative(target)
            output_hashes[variant]=sha(target)
        reverse=reverse_collage(frames,widths,loaded['wavelet'],loaded['optica'])
        changed=int(np.count_nonzero(raw!=old_raw))
        record=dict(ensaio=run,padrao=old['padrao'],numero=old['numero'],
                    comprimento=int(widths.sum()),old_comprimento=old['comprimento'],
                    wavelet_path=destinations['wavelet'],optica_path=destinations['optica'],
                    hashes=output_hashes,shifts_changed_count=changed,
                    supplied_pattern=pattern,frame_count=len(files),height=80,
                    alignment_rows_inclusive=[18,71],alignment_row_count=54,
                    strip_widths_changed_count=int(np.count_nonzero(widths!=old_widths)),
                    nonpositive_raw_lags=int(np.count_nonzero(raw<=0)),
                    old_nonpositive_raw_lags=int(np.count_nonzero(old_raw<=0)),
                    old_hashes=prior_hashes)
        records.append(record)
        for index,(path,width) in enumerate(zip(files,widths)):
            pair_index=min(index,len(raw)-1)
            shifts_table.append(dict(ensaio=run,quadro=path.name,indice_quadro=index,
                deslocamento_bruto=int(raw[pair_index]),largura_tira=int(width),
                deslocamento_bruto_anterior=int(old_raw[pair_index]),
                largura_tira_anterior=int(old_widths[index]),
                ultimo_deslocamento_repetido=index==len(files)-1,
                pearson_maximo=float(scores[pair_index]) if np.isfinite(scores[pair_index]) else '',
                busca_literal_por_empate=any(d['pair']==pair_index for d in fallback)))
        audits.append(dict(ensaio=run,literal_checks=literal_checks,
                           near_maximum_fallbacks=fallback,reverse_collage=reverse))
        frame_total+=len(files);literal_total+=len(literal_checks)
        print(f'{old["numero"]:02d}/26 {run}: {old["comprimento"]} -> {widths.sum()} px, '
              f'{changed}/500 lags changed, {len(literal_checks)} literal checks',flush=True)
    with (OUTPUT/'deslocamentos.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(shifts_table[0]))
        writer.writeheader();writer.writerows(shifts_table)
    differences=[dict(ensaio=r['ensaio'],numero=r['numero'],padrao=r['padrao'],
                      comprimento=r['comprimento'],old_comprimento=r['old_comprimento'],
                      diferenca_comprimento=r['comprimento']-r['old_comprimento'],
                      shifts_changed_count=r['shifts_changed_count'],
                      strip_widths_changed_count=r['strip_widths_changed_count']) for r in records]
    with (OUTPUT/'diferencas.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(differences[0]))
        writer.writeheader();writer.writerows(differences)
    save_json(OUTPUT/'verificacao.json',dict(status='PASS',synthetic_cases=synthetic,
              frame_count=frame_total,literal_pairs_checked=literal_total,runs=audits,
              near_maximum_tolerance=NEAR_MAXIMUM_TOLERANCE,
              changed_variable_only='Alignment profile band, from rows20:60 to rows18:72.',
              method='Overlap Pearson, ascending lags, FFT batched products with literal fallback. '
                     'Every saved pixel checked by reversing the strip assembly against raw frames.',
              preserved_rules=['same26runs','same501frames','sameprovidedlabels',
                               'lags_minus200_to200','first_maximum','minimum_strip_width1',
                               'repeat_last_shift','center_strips','v11_wavelet_seams','clean_optical_seams']))
    changed_inputs=[path for path,digest in before.items() if sha(PROJECT/path)!=digest]
    assert not changed_inputs,changed_inputs
    save_json(OUTPUT/'fontes_sha256.json',dict(algorithm='sha256',files=before,
              before_after_identical=True,checked_files=len(before)))
    output_paths=sorted(p for p in OUTPUT.rglob('*') if p.is_file())
    save_json(OUTPUT/'artefatos_sha256.json',dict(algorithm='sha256',files={
              p.relative_to(OUTPUT).as_posix():sha(p) for p in output_paths}))
    manifest=dict(version='v18_montagem_18_71',complete=True,records=records,
                  processed_acquisitions=26,raw_frames=frame_total,
                  alignment_rows_inclusive=[18,71],alignment_slice=[18,72],
                  preserved_previous_inputs=True,selection_source=relative(selection_path),
                  wavelet_seams='Exact historical v11 rule, including its column overwrite.',
                  optical_seams='Direct concatenation, no overwrite, same strips and widths.',
                  total_changed_pair_lags=sum(r['shifts_changed_count'] for r in records),
                  runtime_seconds=time.monotonic()-started,
                  runtime=dict(python=platform.python_version(),numpy=np.__version__),
                  reconstruction_script_sha256=sha(Path(__file__)))
    temporary=OUTPUT/'manifesto.tmp.json'
    save_json(temporary,manifest)
    os.replace(temporary,final_manifest)
    print(json.dumps({k:v for k,v in manifest.items() if k!='records'},ensure_ascii=False,indent=2),flush=True)


if __name__=='__main__':
    main()
