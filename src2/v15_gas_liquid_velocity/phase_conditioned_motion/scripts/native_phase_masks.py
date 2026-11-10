#!/usr/bin/env python3
"""Apply the current v14 shadow-envelope detector to native video frames.

This preserves source images and previous v14/v15 artifacts. Mask value 1 means
inferred gas silhouette; value 0 is unmarked, not independently proven liquid.
No geometric void fraction is inferred here.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import re
import sys
import time

import numpy as np
from PIL import Image

BASE = Path(__file__).resolve().parents[1]
V15 = BASE.parent
SRC2 = V15.parent
DATA = SRC2 / 'data_processed'
V14 = SRC2 / 'v14_new_voidfraction'
MASKS = BASE / 'intermediate' / 'masks'
OUT = BASE / 'results'
for sub in (MASKS, OUT, BASE/'findings', BASE/'intermediate'/'phase_mask_mplcache'):
    sub.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(BASE/'intermediate'/'phase_mask_mplcache'))
sys.path.insert(0, str(V14/'scripts'))
import panorama_void_alternatives as v14
from v14_geometry import LUMEN_ROWS, BOUNDARY_PROVENANCE

THRESHOLDS = (0.60, 0.70, 0.80)
KEYS = ('masks_060', 'masks_070', 'masks_080')
PATTERNS = {'dispersed bubbles', 'plug', 'slug', 'churn', 'wavy annular'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def selected_runs():
    included, excluded = [], []
    for folder in sorted(DATA.glob('m??q??_???')):
        markers = list(folder.glob('__*.txt'))
        assert len(markers) == 1, folder
        pattern = markers[0].stem[2:].replace('_', ' ')
        row = dict(run=folder.name, pattern=pattern)
        if folder.name == 'm05q05_002' or pattern not in PATTERNS:
            excluded.append(row)
        else:
            included.append(row)
    assert len(included) == 26 and len(excluded) == 5
    return included, excluded


def synthetic_checks(reference):
    # The neutral reference itself has no inferred gas. Uniform images use a
    # neutral radial reference to test the segmentation rather than lens shape.
    rows = []
    for threshold in THRESHOLDS:
        for kind, image, ref in (
            ('uniform', np.full((80,1016),128.), np.ones(80)),
            ('liquid_reference', np.tile(100*reference[:,None],(1,1016)), reference),
        ):
            mask, _, _ = v14.silhouette(image, ref, threshold=threshold)
            assert not mask.any(), (kind, threshold)
            packed = np.packbits(mask, axis=-1, bitorder='big')
            unpacked = np.unpackbits(packed, axis=-1, count=1016, bitorder='big').astype(bool)
            assert np.array_equal(mask, unpacked)
            rows.append(dict(control=kind, threshold=threshold, gas_pixels=int(mask.sum()),
                             pack_roundtrip_exact=True))
    return rows


def process_run(task):
    row, reference = task
    run, pattern = row['run'], row['pattern']
    folder = DATA/run
    files = sorted(folder.glob('ImgA*.jpg'))
    ids = np.array([int(re.search(r'(\d+)\.jpg$',f.name)[1]) for f in files],dtype=np.int32)
    assert len(files) == 501 and np.array_equal(ids,np.arange(501))
    mii = (folder/'_config.mii').read_text()
    fps = float(re.search(r'Film speed 1\s*:([0-9.]+)',mii)[1])
    xsv = (folder/'_config.xsv').read_text()
    assert fps == float(re.search(r'^Rate=([0-9.]+)',xsv,re.M)[1])
    assert int(re.search(r'^Images=([0-9]+)',xsv,re.M)[1]) == 501
    packed = {key: np.empty((501,80,127),dtype=np.uint8) for key in KEYS}
    totals = {key: 0 for key in KEYS}
    disagreement = 0
    nonmonotonic = 0
    saturation = []
    frame_rows = []
    digest = hashlib.sha256()
    direct_checks = 0
    for i,file in enumerate(files):
        raw = file.read_bytes()
        digest.update(file.name.encode())
        digest.update(hashlib.sha256(raw).digest())
        image = np.asarray(Image.open(io.BytesIO(raw)).convert('L'),dtype=float)
        assert image.shape == (80,1016)
        masks = []
        frow = dict(run=run, pattern=pattern, frame_id=int(ids[i]), fps=fps)
        for key,threshold in zip(KEYS,THRESHOLDS):
            mask, _, _ = v14.silhouette(image,reference,threshold=threshold)
            assert mask.dtype == bool and mask.shape == (80,1016)
            assert not mask[:LUMEN_ROWS[0]].any() and not mask[LUMEN_ROWS[1]:].any()
            packed[key][i] = np.packbits(mask,axis=-1,bitorder='big')
            totals[key] += int(mask.sum())
            frow[key+'_lumen_fraction'] = float(mask[LUMEN_ROWS[0]:LUMEN_ROWS[1]].mean())
            masks.append(mask)
            # Direct result vs the saved representation on fixed indices.
            if i in (0,200,500):
                restored = np.unpackbits(packed[key][i],axis=-1,count=1016,bitorder='big').astype(bool)
                assert np.array_equal(restored,mask)
                direct_checks += 1
        disagreement += int(np.logical_or.reduce(masks).sum()-np.logical_and.reduce(masks).sum())
        nonmonotonic += int(((masks[0] & ~masks[1]) | (masks[1] & ~masks[2])).sum())
        sat = float(np.mean(image[20:60]>=250))
        saturation.append(sat)
        frow['core_saturation_fraction'] = sat
        frow['threshold_disagreement_lumen_fraction'] = float(
            (np.logical_or.reduce(masks) ^ np.logical_and.reduce(masks))[LUMEN_ROWS[0]:LUMEN_ROWS[1]].mean())
        frame_rows.append(frow)
    output = MASKS/(run+'.npz')
    # All array names and the full unpacked shape are explicitly stored.
    np.savez_compressed(output, **packed,
        mask_shape=np.array([501,80,1016],dtype=np.int32), original_width=np.int32(1016),
        frame_ids=ids, fps=np.float64(fps), thresholds=np.array(THRESHOLDS,dtype=np.float64),
        lumen_rows=np.array(LUMEN_ROWS,dtype=np.int32), bitorder=np.array('big'),
        unpacked_dtype=np.array('bool'), mask_meaning=np.array('inferred_gas_silhouette'))
    # Reopen the actual NPZ and compare packed bytes to those made from direct
    # detector outputs; this checks every pixel/threshold/frame representation.
    with np.load(output,allow_pickle=False) as saved:
        for key in KEYS:
            assert np.array_equal(saved[key],packed[key])
        assert tuple(saved['mask_shape']) == (501,80,1016)
    den = 501*(LUMEN_ROWS[1]-LUMEN_ROWS[0])*1016
    summary = dict(run=run, pattern=pattern, frames=501, fps=fps,
        duration_s=500/fps, first_frame_id=0, last_frame_id=500,
        height=80, width=1016, lumen_start=LUMEN_ROWS[0], lumen_stop=LUMEN_ROWS[1],
        core_saturation_fraction=float(np.mean(saturation)),
        threshold_disagreement_lumen_fraction=disagreement/den,
        nonmonotonic_lumen_fraction=nonmonotonic/den,
        raw_frame_sha256=digest.hexdigest(), mii_sha256=sha(folder/'_config.mii'),
        xsv_sha256=sha(folder/'_config.xsv'), packed_file=str(output.relative_to(BASE)),
        packed_sha256=sha(output), packed_bytes=output.stat().st_size,
        selected_direct_unpack_checks=direct_checks,
        all_saved_packed_arrays_exact=True)
    for key in KEYS:
        summary[key+'_gas_pixels'] = totals[key]
        summary[key+'_lumen_fraction'] = totals[key]/den
    return summary, frame_rows


def write_csv(path, rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def write_report(summaries, metadata):
    panorama_path = V14/'results'/'void_fraction'/'panorama_void_per_run.csv'
    with panorama_path.open() as f:
        panorama = {r['run']:r for r in csv.DictReader(f)}
    comparison = []
    for row in summaries:
        a = float(row['masks_070_lumen_fraction'])
        b = float(panorama[row['run']]['gas_silhouette_projected_fraction'])
        comparison.append(dict(run=row['run'],pattern=row['pattern'],
            native_temporal_spatial_gas_fraction=a,panorama_spatial_gas_fraction=b,
            difference_percentage_points=100*(a-b)))
    write_csv(OUT/'native_panorama_comparison.csv',comparison)
    lines = [
        '# Native phase-mask transfer from v14', '',
        '**Completed:** the current v14 detector was applied to all 13,026 native frames in the 26 retained recordings, at shadow thresholds 0.60, 0.70 and 0.80. These 39,078 masks are gas-silhouette candidates. Unmarked pixels are not independently confirmed liquid.', '',
        'The exact three-panorama liquid-like radial reference was recomputed and frozen. Native images, v14 scripts and previous outputs were not changed. The current owner-confirmed lumen is rows 18:71, or 53 pixel cells; no geometric void fraction or physical pixel scale is used in this transfer.', '',
        '## Per-recording projected fractions and sensitivity', '',
        'Fractions below average inferred gas pixels over all 501 native frames and all 53 lumen rows. They are projected optical fractions, not vapor quality or volumetric void. Disagreement is the fraction of lumen pixels whose labels differ among the three thresholds.', '',
        '| Recording | Pattern | Threshold 0.60 (%) | Threshold 0.70 (%) | Threshold 0.80 (%) | Disagreement (%) |',
        '|---|---|---:|---:|---:|---:|',
    ]
    for row in summaries:
        values=[100*float(row[k]) for k in ('masks_060_lumen_fraction','masks_070_lumen_fraction',
                                           'masks_080_lumen_fraction','threshold_disagreement_lumen_fraction')]
        lines.append('| '+row['run']+' | '+row['pattern']+' | '+' | '.join(f'{v:.3f}' for v in values)+' |')
    bubble=next(r for r in summaries if r['run']=='m07q01_001')
    c=next(r for r in comparison if r['run']=='m07q01_001')
    lines += ['', '## Concrete transfer difficulties', '',
        f'**Dispersed-bubble m07q01_001 is a major sensitivity case:** the native projected gas fraction changes from {100*float(bubble["masks_060_lumen_fraction"]):.3f}% at threshold 0.60 to {100*float(bubble["masks_070_lumen_fraction"]):.3f}% at 0.70 and {100*float(bubble["masks_080_lumen_fraction"]):.3f}% at 0.80. The existing current-v14 panorama projected fraction is {100*c["panorama_spatial_gas_fraction"]:.3f}%. Native and panorama sampling differ, so this discrepancy does not identify which is correct; it demonstrates that panorama evidence cannot automatically validate native-frame labels.', '',
        'In contrast, m05q01_001 and m03q05_001 have close native and saved-panorama projected fractions. The all-run comparison is saved with explicit temporal/spatial definitions. Do not compare projected fractions to circular-volume alpha, which is a different quantity.', '',
        'The mask operation is not guaranteed monotonic in threshold: thresholding changes component connectivity, crop-contact rejection and hole filling. Observed nonmonotonic pixels reach approximately 1.061% of lumen pixels in a recording. Conservative consensus must intersect/union all three masks rather than assume one contains the others.', '',
        'Gas-dominated churn and annular masks can leave very little detected liquid at a fixed probe. That is a phase-support limitation for a two-class conditional average. It does not establish that no liquid film or entrained liquid is physically present.', '',
        '## Saved schema and verification', '',
        'Each intermediate/masks/RUN.npz contains masks_060, masks_070 and masks_080 as uint8 arrays with shape (501,80,127), packed along the axial axis with big-endian bit order. It also stores mask_shape=(501,80,1016), original_width=1016, frame_ids, fps, thresholds, lumen_rows, bitorder and unpacked_dtype=bool.', '',
        'Unpack with numpy.unpackbits(npz[key], axis=-1, count=int(npz["original_width"]), bitorder="big").astype(bool). Value 1 is inferred gas silhouette; value 0 is unmarked. Pixels outside rows 18:71 are always unmarked and must not become liquid labels.', '',
        'Uniform images under a neutral reference and exact reference-shaped images produce no inferred gas at every threshold. Packed representations were compared with direct detector outputs, then reopened from disk and checked byte-for-byte. Fixed frames 0, 200 and 500 were also unpacked and compared exactly. Source, reference and packed-file hashes are retained. These tests establish faithful application and storage, not phase-label accuracy.', '',
        'Reproduce with python3 -B scripts/native_phase_masks.py --workers 3 from the repository path for this script.', '',
        'Evidence: [inventory](../results/native_phase_inventory.csv), [per-frame fractions](../results/native_phase_frame_fractions.csv), [reference](../results/native_phase_reference.json), [metadata and checks](../results/native_phase_metadata.json), [native/panorama comparison](../results/native_panorama_comparison.csv).',
    ]
    (BASE/'findings'/'native_phase_masks.md').write_text('\n'.join(lines)+'\n')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--workers',type=int,default=3)
    args=parser.parse_args()
    assert 1<=args.workers<=4
    started=time.monotonic()
    included,excluded=selected_runs()
    reference,selections=v14.make_reference()
    assert reference.shape==(80,) and np.all(np.isfinite(reference))
    refmeta=dict(reference_runs=list(v14.REFERENCE_RUNS),
        radial_reference=reference.tolist(),candidate_selection=selections,
        reference_images=[dict(run=r,path=str((DATA/r/'__panorama_v11.png').relative_to(SRC2)),
                               sha256=sha(DATA/r/'__panorama_v11.png')) for r in v14.REFERENCE_RUNS],
        definition='Exact current-v14 liquid-like three-panorama reference. It is inferred, not independently acquired liquid-only calibration.',
        lumen_rows=list(LUMEN_ROWS),boundary_provenance=BOUNDARY_PROVENANCE)
    (OUT/'native_phase_reference.json').write_text(json.dumps(refmeta,indent=2)+'\n')
    controls=synthetic_checks(reference)
    summaries,frame_rows=[],[]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(process_run,(r,reference)):r['run'] for r in included}
        for future in as_completed(futures):
            summary,frames=future.result()
            summaries.append(summary)
            frame_rows.extend(frames)
            print(summary['run'],'done; baseline projected fraction',
                  round(summary['masks_070_lumen_fraction'],4),'threshold disagreement',
                  round(summary['threshold_disagreement_lumen_fraction'],4),flush=True)
    summaries.sort(key=lambda r:r['run'])
    frame_rows.sort(key=lambda r:(r['run'],r['frame_id']))
    assert len(summaries)==26 and len(frame_rows)==26*501
    write_csv(OUT/'native_phase_inventory.csv',summaries)
    write_csv(OUT/'native_phase_frame_fractions.csv',frame_rows)
    metadata=dict(status='complete',included_runs=included,excluded_runs=excluded,
        thresholds=list(THRESHOLDS),lumen_rows=list(LUMEN_ROWS),boundary_provenance=BOUNDARY_PROVENANCE,
        schema=dict(keys=list(KEYS),packed_shape=[501,80,127],packed_dtype='uint8',
                    original_shape=[501,80,1016],axis=-1,bitorder='big',unpacked_dtype='bool',
                    unpack_instruction="np.unpackbits(npz[key], axis=-1, count=int(npz['original_width']), bitorder='big').astype(bool)"),
        source_files=[dict(path=str((V14/'scripts'/name).relative_to(SRC2)),sha256=sha(V14/'scripts'/name))
                      for name in ('panorama_void_alternatives.py','v14_geometry.py','v14_selection.py','v14_paths.py')],
        reference_json_sha256=sha(OUT/'native_phase_reference.json'),synthetic_checks=controls,
        total_frames=26*501,total_detector_calls=26*501*3,threshold_mask_files=26,
        source_images_and_previous_outputs_modified=False,
        interpretation='Mask 1=inferred visible gas silhouette; 0=unmarked, not independently proven liquid. No circular-volume alpha is calculated.',
        limitations=[
            'The v14 detector was developed on panoramas; application to native frames is a transfer requiring visual checks.',
            'Reference liquid-like intervals are inferred and may contain unresolved gas or optical artifacts.',
            'Gas shadows may be missed or merged; unmarked regions cannot be certified liquid.',
            'Agreement across thresholds is segmentation stability, not ground-truth correctness.',
            'Morphology can make masks nonmonotonic with threshold; consensus must use all masks.',
            'Saturated native regions lose information; those records remain inventoried and flagged.'
        ],elapsed_s=time.monotonic()-started)
    (OUT/'native_phase_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    write_report(summaries, metadata)
    print('Complete:',metadata['total_frames'],'frames;',metadata['total_detector_calls'],
          'masks;',round(metadata['elapsed_s'],1),'seconds',flush=True)


if __name__=='__main__':
    main()
