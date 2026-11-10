#!/usr/bin/env python3
"""Independent morphology, geometry and before/after checks; run after outputs."""
from v14_paths import panorama_path, original_panorama_path, BASELINE
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json

import numpy as np
from PIL import Image
from scipy.ndimage import binary_fill_holes, gaussian_filter1d, label

from v14_paths import ROOT, artifact_path, write_markdown
from v14_selection import iter_runs, is_eligible, PATTERN_ORDER
from v14_geometry import (LUMEN_ROWS, LUMEN_SLICE, LUMEN_HEIGHT_PX, DIAMETER_PX,
                          FIXED_PROBE_ROW, BOUNDARY_PROVENANCE)
import panorama_void_alternatives as estimator


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def independent_closing(crop):
    """3x3 dilation then erosion on an infinite zero-extended binary plane."""
    h, w = crop.shape
    padded = np.pad(np.asarray(crop, dtype=bool), 2)
    dilated = np.logical_or.reduce([padded[dy:dy+h+2, dx:dx+w+2]
                                  for dy in range(3) for dx in range(3)])
    return np.logical_and.reduce([dilated[dy:dy+h, dx:dx+w]
                                  for dy in range(3) for dx in range(3)])


def direct_loop_closing(crop):
    h, w = crop.shape
    padded = np.pad(np.asarray(crop, dtype=bool), 2)
    dilated = np.zeros((h+2, w+2), dtype=bool)
    for y in range(h+2):
        for x in range(w+2):
            dilated[y, x] = any(bool(padded[y+dy, x+dx]) for dy in range(3) for dx in range(3))
    closed = np.zeros((h,w), dtype=bool)
    for y in range(h):
        for x in range(w):
            closed[y, x] = all(bool(dilated[y+dy, x+dx]) for dy in range(3) for dx in range(3))
    return closed


def isolated_shadow_candidates(closed):
    components, count = label(closed)
    radial_labels = np.unique(np.r_[components[0], components[-1]])
    sizes = np.bincount(components.ravel(), minlength=count+1)
    retained = (~np.isin(components, radial_labels)) & (sizes[components] >= 4)
    return binary_fill_holes(retained)


def independent_mask(im, reference, threshold=.70, rows=(18,71), min_open_interval=20):
    normalized = im / np.maximum(np.percentile(im[30:50],75,axis=0),1)[None,:]
    ratio = gaussian_filter1d(normalized,1,axis=0) / np.maximum(reference[:,None],.01)
    cols = np.arange(im.shape[1])
    upper = 19 + np.argmin(ratio[19:38],axis=0)
    lower = 42 + np.argmin(ratio[42:71],axis=0)
    accepted = (ratio[upper,cols]<threshold) & (ratio[lower,cols]<threshold)
    shadow = ratio<threshold
    closed = independent_closing(shadow[rows[0]:rows[1]])
    components,_ = label(closed)
    top = components[np.clip(upper-rows[0],0,closed.shape[0]-1),cols]
    bottom = components[np.clip(lower-rows[0],0,closed.shape[0]-1),cols]
    same = (top>0) & (top==bottom)
    intervals,_ = label(accepted)
    sizes = np.bincount(intervals)
    accepted &= same | (sizes[intervals]>=min_open_interval)
    outer_top, outer_bottom = upper.copy(),lower.copy()
    for _ in range(32):
        go = (outer_top>rows[0]) & shadow[np.maximum(outer_top-1,0),cols]
        outer_top[go]-=1
        go = (outer_bottom<rows[1]-1) & shadow[np.minimum(outer_bottom+1,79),cols]
        outer_bottom[go]+=1
    mask = np.zeros_like(im,dtype=bool)
    for y in range(rows[0],rows[1]):
        mask[y]=accepted & (y>=outer_top) & (y<=outer_bottom)
    mask[rows[0]:rows[1]] |= isolated_shadow_candidates(closed)
    return mask


def area_by_intervals(mask, diameter=DIAMETER_PX):
    """Independent interval endpoints, rather than the estimator accumulator."""
    transitions = np.diff(np.pad(mask.astype(np.int8),((1,1),(0,0))),axis=0)
    result=np.zeros(mask.shape[1],dtype=float)
    for x in range(mask.shape[1]):
        starts=np.flatnonzero(transitions[:,x]==1)
        stops=np.flatnonzero(transitions[:,x]==-1)
        result[x]=np.sum((stops-starts)**2)/diameter**2
    return result


def flips(arr):
    return int(np.count_nonzero(arr[...,1:]!=arr[...,:-1]))


def analytic_controls():
    assert LUMEN_ROWS == (18,71) and LUMEN_HEIGHT_PX == 53 and DIAMETER_PX == 53.0
    assert estimator.DIAMETER_PX == 53.0 and FIXED_PROBE_ROW == 40
    cases=[np.zeros((9,12),bool),np.ones((9,12),bool)]
    top=np.zeros((9,12),bool);top[0:5,4:9]=True;cases.append(top)
    left=np.zeros((9,12),bool);left[3:7,0:4]=True;cases.append(left)
    ring=np.zeros((11,13),bool);ring[2:9,3:10]=True;ring[3:8,4:9]=False;cases.append(ring)
    rng=np.random.default_rng(20261023)
    cases.extend(rng.random((9,12))<p for p in (.05,.2,.5,.8,.95))
    for case in cases:
        expected=direct_loop_closing(case)
        assert np.array_equal(independent_closing(case),expected)
        assert np.array_equal(estimator.close_shadow_crop(case),expected)
        assert np.all(~case | expected), 'Zero-extended closing must be extensive'
    assert not isolated_shadow_candidates(estimator.close_shadow_crop(top)).any()
    assert np.array_equal(isolated_shadow_candidates(estimator.close_shadow_crop(left)),left)
    assert isolated_shadow_candidates(estimator.close_shadow_crop(left))[:,0].any()
    assert not isolated_shadow_candidates(estimator.close_shadow_crop(np.ones((9,12),bool))).any()
    ring_filled=isolated_shadow_candidates(estimator.close_shadow_crop(ring))
    assert ring_filled[5,6] and int(ring_filled.sum())==49
    for level in (50,100,200):
        assert not estimator.silhouette(np.full((80,120),level,dtype=float),np.ones(80))[0].any()
    noise=100+rng.uniform(-1,1,size=(80,120))
    assert not estimator.silhouette(noise,np.ones(80))[0].any()
    reference,_=estimator.make_reference()
    assert not estimator.silhouette(np.tile(100*reference[:,None],(1,120)),reference)[0].any()
    body=np.full((80,140),100.,dtype=float)
    body[24:28,40:100]=25;body[62:66,40:100]=25
    body_mask=estimator.silhouette(body,np.ones(80))[0]
    assert body_mask[40,70] and not body_mask[:,10].any() and flips(body_mask[40])==2
    empty=np.zeros((80,10),bool)
    full=empty.copy();full[18:71]=True
    separate=empty.copy();separate[20:30]=True;separate[50:55]=True
    for mask in (empty,full,separate):
        assert np.allclose(estimator.geometry(mask),area_by_intervals(mask),atol=1e-15,rtol=0)
    assert np.allclose(estimator.geometry(full),1.0)
    assert np.allclose(estimator.geometry(separate),(10**2+5**2)/53**2)
    assert np.allclose(area_by_intervals(full,60),(53/60)**2)
    return dict(morphology_reference_cases=len(cases),direct_loop_and_boolean_reduction_agree=True,
                radial_contact_rejected=True,axial_contact_preserved=True,
                interior_hole_filled=True,flat_intensity_controls=3,
                bounded_noise_control=True,liquid_profile_control=True,bright_gas_interior_control=True,
                separated_interval_geometry_control=True,full_lumen_maximum=1.0,
                confirmed_lumen_rows=list(LUMEN_ROWS),confirmed_diameter_px=DIAMETER_PX,
                historical_D60_full_mask_value=(53/60)**2,
                boundary_source=BOUNDARY_PROVENANCE)


def table(path):
    with path.open(newline='') as fp:return list(csv.DictReader(fp))



def main():
    controls=analytic_controls()
    reference,_=estimator.make_reference()
    void={r['run']:r for r in table(artifact_path('panorama_void_per_run.csv'))}
    phase={r['run']:r for r in table(artifact_path('phase_mask_flip_per_run.csv'))}
    records=[]
    for folder in iter_runs():
        run=folder.name;im=estimator.read(run)
        path=artifact_path('phase_masks')/f'{run}_mask.png'
        saved=np.asarray(Image.open(path))==255
        independent=independent_mask(im,reference)
        assert np.array_equal(saved,independent),run
        assert np.array_equal(estimator.silhouette(im,reference)[0],independent),run
        area=area_by_intervals(independent)
        error=abs(float(area.mean())-float(void[run]['conditional_circle_alpha_spatial_mean']))
        count=flips(independent[18:71])
        assert error<1e-14 and count==int(phase[run]['candidate_pixel_flip_count'])
        prior_path=BASELINE/'intermediate/phase_masks'/path.name
        prior=np.asarray(Image.open(prior_path))==255
        old_alpha=float(area_by_intervals(prior).mean())
        records.append(dict(run=run,pattern=phase[run]['pattern'],
            eligible=phase[run]['eligible']=='True',alpha_before=old_alpha,
            alpha_after=float(area.mean()),alpha_change=float(area.mean())-old_alpha,
            diameter_before_px=53,diameter_after_px=53,flip_count_before=flips(prior[18:71]),
            flip_count_after=count,changed_pixels=int(np.count_nonzero(prior!=saved)),
            before_mask_sha256=sha256(prior_path),after_mask_sha256=sha256(path),
            independent_area_error=error,input_panorama_sha256=sha256(panorama_path(run))))
    with artifact_path('deep_method_changes.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    result=dict(passed=True,controls=controls,selected_acquisitions=len(records),
        statistical_acquisitions=sum(r['eligible'] for r in records),
        changed_masks=sum(r['changed_pixels']>0 for r in records),
        max_independent_area_error=max(r['independent_area_error'] for r in records),
        segmentation_source_sha256=sha256(artifact_path('panorama_void_alternatives.py')),
        failures=[],interpretation='Independent morphology and interval arithmetic on corrected inputs. Both versions use D53; changes compare v14 with v14.1, not the historical morphology or D60 correction.')
    artifact_path('deep_method_checks.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    write_markdown(artifact_path('deep_method_review.md'),
        '# Revisão independente dos métodos da v14.1\n\n'
        '**PASS:** os controles de fechamento morfológico, bordas, perfis líquidos, áreas e transições passaram. '
        f'As {len(records)} máscaras foram reproduzidas por operações booleanas independentes; '
        f'{result["changed_masks"]} diferiram da v14 após corrigir a entrada e reconstruir a referência.\n\n'
        'A segmentação, seus parâmetros e D=53 foram preservados. A comparação usa D53 em ambos os lados. '
        'A conta de α continua condicionada à identidade das fases, à refração e a seções circulares.\n\n'
        '[Verificações](deep_method_checks.json) · [Alterações por aquisição](deep_method_changes.csv)\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
