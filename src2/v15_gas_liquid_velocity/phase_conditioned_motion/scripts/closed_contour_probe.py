#!/usr/bin/env python3
"""Test disabling long open-envelope acceptance in native v14 silhouettes.

Closed-only is a topology restriction, not verified phase ground truth. Source
images and baseline masks are preserved. Every comparison uses identical frames.
"""
from pathlib import Path
import csv
import hashlib
import json
import os
import sys
import time
import numpy as np
from PIL import Image

BASE=Path(__file__).resolve().parents[1]
SRC2=BASE.parent.parent
V14=SRC2/'v14_new_voidfraction'
DATA=SRC2/'data_processed'
OUT=BASE/'results'
FIG=BASE/'figures'
os.environ.setdefault('MPLCONFIGDIR',str(BASE/'intermediate'/'closed_contour_mplcache'))
sys.path.insert(0,str(V14/'scripts'))
import panorama_void_alternatives as detector
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

KEYS=('masks_060','masks_070','masks_080')
THRESHOLDS=(.6,.7,.8)
OPEN_DISABLED=1000000


def write_csv(name,rows):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def overlay(im,mask):
    rgb=np.repeat((im/255.)[...,None],3,axis=2)
    tint=np.array([.02,.45,1.])
    rgb[mask]=.58*rgb[mask]+.42*tint
    return np.clip(rgb,0,1)


def run():
    started=time.monotonic()
    reference_path=OUT/'native_phase_reference.json'
    reference=np.asarray(json.loads(reference_path.read_text())['radial_reference'])
    inventory=list(csv.DictReader((OUT/'native_phase_inventory.csv').open()))
    summary=[];frame_rows=[];examples={};provenance=[]
    for inv in inventory:
        run=inv['run']
        ids=np.arange(501) if run=='m07q01_001' else np.linspace(0,500,20,dtype=int)
        files=sorted((DATA/run).glob('ImgA*.jpg'))
        assert len(files)==501
        path=BASE/'intermediate'/'masks'/(run+'.npz')
        with np.load(path,allow_pickle=False) as z:
            baseline={k:np.unpackbits(z[k][ids],axis=-1,count=1016,bitorder='big').astype(bool) for k in KEYS}
        closed={k:np.empty((len(ids),80,127),dtype=np.uint8) for k in KEYS}
        totals={k:dict(base=0,closed=0,changed=0,added=0) for k in KEYS}
        disagreement=0;baseline_disagreement=0
        digest=hashlib.sha256()
        for j,i in enumerate(ids):
            raw=files[i].read_bytes();digest.update(files[i].name.encode());digest.update(hashlib.sha256(raw).digest())
            im=np.asarray(Image.open(files[i]).convert('L'),dtype=float)
            assert im.shape==(80,1016)
            current=[]
            row=dict(run=run,pattern=inv['pattern'],frame_id=int(i))
            for k,threshold in zip(KEYS,THRESHOLDS):
                mask,_,_=detector.silhouette(im,reference,threshold=threshold,min_open_interval=OPEN_DISABLED)
                assert mask.dtype==bool
                base=baseline[k][j]
                # Disabling the open acceptance term must only remove pixels.
                assert not np.any(mask & ~base)
                closed[k][j]=np.packbits(mask,axis=-1,bitorder='big')
                totals[k]['base']+=int(base.sum());totals[k]['closed']+=int(mask.sum())
                totals[k]['changed']+=int(np.count_nonzero(base!=mask))
                totals[k]['added']+=int(np.count_nonzero(mask & ~base))
                row[k+'_baseline_fraction']=float(base[18:71].mean())
                row[k+'_closed_fraction']=float(mask[18:71].mean())
                current.append(mask)
            d=np.logical_or.reduce(current)^np.logical_and.reduce(current)
            bd=np.logical_or.reduce([baseline[k][j] for k in KEYS])^np.logical_and.reduce([baseline[k][j] for k in KEYS])
            disagreement+=int(d.sum());baseline_disagreement+=int(bd.sum())
            row['closed_threshold_disagreement']=float(d[18:71].mean())
            row['baseline_threshold_disagreement']=float(bd[18:71].mean())
            frame_rows.append(row)
            if run=='m07q01_001' and i in (200,210):
                examples[int(i)]=(im,[baseline[k][j].copy() for k in KEYS],[m.copy() for m in current])
        den=len(ids)*53*1016
        for k,threshold in zip(KEYS,THRESHOLDS):
            summary.append(dict(run=run,pattern=inv['pattern'],threshold=threshold,n_frames=len(ids),
                baseline_projected_gas_fraction=totals[k]['base']/den,
                closed_projected_gas_fraction=totals[k]['closed']/den,
                removed_fraction_of_lumen=totals[k]['changed']/den,
                removed_fraction_of_baseline=totals[k]['changed']/totals[k]['base'] if totals[k]['base'] else float('nan'),
                added_pixels=totals[k]['added'],
                baseline_threshold_disagreement=baseline_disagreement/den,
                closed_threshold_disagreement=disagreement/den))
        if run=='m07q01_001':
            output=OUT/'closed_contour_m07q01_masks.npz'
            np.savez_compressed(output,**closed,frame_ids=ids,mask_shape=np.array([501,80,1016]),
                original_width=1016,thresholds=np.asarray(THRESHOLDS),lumen_rows=np.array([18,71]),
                bitorder=np.array('big'),min_open_interval=OPEN_DISABLED)
            with np.load(output,allow_pickle=False) as z:
                for k in KEYS:assert np.array_equal(z[k],closed[k])
        provenance.append(dict(run=run,frame_ids=ids.tolist(),selected_frame_content_sha256=digest.hexdigest(),
            baseline_npz_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        print(run,'frames',len(ids),'baseline/closed .70',summary[-2]['baseline_projected_gas_fraction'],summary[-2]['closed_projected_gas_fraction'],flush=True)
    assert len(inventory)==26 and len(frame_rows)==1001 and len(summary)==78
    write_csv('closed_contour_summary.csv',summary)
    write_csv('closed_contour_frames.csv',frame_rows)
    for frame,(im,bases,closed) in examples.items():
        fig,axes=plt.subplots(4,2,figsize=(15,6))
        for ax in axes[0]:ax.imshow(im,cmap='gray',vmin=0,vmax=255);ax.axis('off')
        axes[0,0].set_title('Raw native frame')
        axes[0,1].set_title('Same raw frame')
        for idx,threshold in enumerate(THRESHOLDS):
            axes[idx+1,0].imshow(overlay(im,bases[idx]));axes[idx+1,0].axis('off')
            axes[idx+1,1].imshow(overlay(im,closed[idx]));axes[idx+1,1].axis('off')
            axes[idx+1,0].set_title(f'Baseline threshold {threshold:.2f}; gas candidate {bases[idx][18:71].mean():.1%}',fontsize=9)
            axes[idx+1,1].set_title(f'Open-envelope disabled; gas candidate {closed[idx][18:71].mean():.1%}',fontsize=9)
        fig.suptitle(f'm07q01_001 frame {frame}: targeted topology repair\nBlue is inferred gas; less blue does not itself validate phase identity',fontsize=12)
        fig.tight_layout(rect=(0,0,1,.91))
        fig.savefig(FIG/f'closed_contour_m07q01_frame{frame}.png',dpi=180)
        plt.close(fig)
    metadata=dict(method='Current v14 silhouette with min_open_interval=1000000; all other inputs unchanged.',
        thresholds=list(THRESHOLDS),lumen_rows=[18,71],n_total_frames=1001,
        full_record='m07q01_001',other_runs_frame_sample='20 equally spaced native indices including endpoints',
        detector_sha256=hashlib.sha256((V14/'scripts'/'panorama_void_alternatives.py').read_bytes()).hexdigest(),
        reference_sha256=hashlib.sha256(reference_path.read_bytes()).hexdigest(),
        subset_assertion_passed=True,packing_roundtrip_exact=True,inputs=provenance,
        elapsed_s=time.monotonic()-started,
        limitations=['Closed topology is not certified gas identity.',
            'True open/cropped gas objects can be lost.',
            'Other recordings use sparse sampling, not full-record statistics.',
            'Neither source supplies independent liquid-only or phase ground truth.'])
    (OUT/'closed_contour_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print('Complete',metadata['elapsed_s'],flush=True)


if __name__=='__main__':run()
