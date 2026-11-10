#!/usr/bin/env python3
"""Targeted native-frame phase-mask review; images are not phase ground truth."""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage import binary_fill_holes, label

BASE=Path(__file__).resolve().parents[1]
SRC2=BASE.parents[1]
DATA=SRC2/'data_processed'
V14=SRC2/'v14_new_voidfraction'/'scripts'
for name in ('results','figures','findings','intermediate'):
    (BASE/name).mkdir(parents=True,exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR',str(BASE/'intermediate'/'visual_mplcache'))
sys.path.insert(0,str(V14))
import panorama_void_alternatives as estimator

RUNS=('m07q01_001','m05q01_001','m05q05_001','m04q01_001',
      'm06q05_001','m05q20_001','m04q80_001','m03q40_001')
FRAMES=(200,210)
THRESHOLDS=(.60,.70,.80)
PROBE_X=(496,520)
LUMEN_Y=(18,71)


def font(size):
    for path in ('/System/Library/Fonts/Supplemental/Arial.ttf',
                 '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'):
        if Path(path).exists(): return ImageFont.truetype(path,size)
    return ImageFont.load_default()


def tint(a,gas,uncertain=None):
    out=np.repeat(a[...,None],3,axis=2).astype(float)
    out[gas]=.6*out[gas]+.4*np.array([0,165,230])
    if uncertain is not None:
        out[uncertain]=.5*out[uncertain]+.5*np.array([245,160,0])
    return Image.fromarray(np.round(out).astype('uint8'))


def isolated_only(ratio,threshold=.70):
    """Existing v14 isolated-shadow branch with the long paired-fill branch off.

    This is a targeted alternative for the known dispersed-bubble recording,
    not a validated replacement for the full detector across all patterns.
    """
    closed=estimator.close_shadow_crop((ratio<threshold)[18:71])
    components,count=label(closed)
    rejected=np.unique(np.concatenate([components[0],components[-1]]))
    sizes=np.bincount(components.ravel(),minlength=count+1)
    retained=(~np.isin(components,rejected))&(sizes[components]>=4)
    out=np.zeros_like(ratio,dtype=bool)
    out[18:71]=binary_fill_holes(retained)
    return out


def dispersed_alternative(reference):
    run='m07q01_001'
    files=sorted((DATA/run).glob('ImgA*.jpg'))
    saved=BASE/'intermediate'/'masks'/f'{run}.npz'
    archived=None
    if saved.exists():
        with np.load(saved) as archive:
            archived=np.unpackbits(archive['masks_070'],axis=-1,count=1016,bitorder='big').astype(bool)
    examples={}
    rows=[]
    for i,p in enumerate(files):
        a=np.asarray(Image.open(p).convert('L'))
        baseline,_,ratio=estimator.silhouette(a.astype(float),reference,threshold=.70)
        if archived is not None:
            assert np.array_equal(baseline,archived[i]),(run,i)
        alternative=isolated_only(ratio)
        rows.append(dict(frame=i,baseline_lumen_fraction=float(baseline[18:71].mean()),
            isolated_lumen_fraction=float(alternative[18:71].mean()),
            baseline_right_edge_fraction=float(baseline[18:71,860:1016].mean()),
            isolated_right_edge_fraction=float(alternative[18:71,860:1016].mean()),
            baseline_probe_fraction=float(baseline[18:71,496:520].mean()),
            isolated_probe_fraction=float(alternative[18:71,496:520].mean())))
        if i in FRAMES:examples[i]=(a,baseline,alternative)
    canvas=Image.new('RGB',(1220,690),'white');draw=ImageDraw.Draw(canvas)
    draw.text((12,9),'m07q01_001 | Dispersed-bubble detector alternative',fill='black',font=font(19))
    draw.text((12,36),'Blue: inferred gas. Alternative keeps isolated-shadow branch and disables long paired filling.',fill='black',font=font(16))
    draw.text((12,59),'Targeted visual diagnostic; reduced mask area is not proof of improved phase accuracy.',fill='black',font=font(16))
    for block,(i,(a,b,c)) in enumerate(examples.items()):
        for k,(name,panel) in enumerate(zip(('Native image','v14 baseline 0.70','Isolated components'),(Image.fromarray(a).convert('RGB'),tint(a,b),tint(a,c)))):
            top=92+block*294+k*91
            canvas.paste(panel,(194,top))
            draw.text((10,top+15),name,fill='black',font=font(14))
            draw.text((10,top+42),f'Frame {i}',fill='black',font=font(14))
    canvas.save(BASE/'figures'/'phase_visual_dispersed_alternative.png')
    result=dict(run=run,frames=len(rows),threshold=.70,
        modification='Use only the existing v14 isolated-shadow component/hole-filling branch; disable paired-envelope fill.',
        status='Targeted exploratory alternative; no independent phase-mask accuracy labels.',
        baseline_archive_exact_match=archived is not None,
        mean_baseline_lumen_fraction=float(np.mean([r['baseline_lumen_fraction'] for r in rows])),
        mean_isolated_lumen_fraction=float(np.mean([r['isolated_lumen_fraction'] for r in rows])),
        mean_baseline_right_edge_fraction=float(np.mean([r['baseline_right_edge_fraction'] for r in rows])),
        mean_isolated_right_edge_fraction=float(np.mean([r['isolated_right_edge_fraction'] for r in rows])),
        per_frame=rows)
    (BASE/'results'/'phase_visual_dispersed_alternative.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def run():
    reference,selections=estimator.make_reference()
    evidence=[]
    for run in RUNS:
        folder=DATA/run
        files=sorted(folder.glob('ImgA*.jpg'))
        assert len(files)==501
        fps=float(re.search(r'^Rate=([0-9.]+)',(folder/'_config.xsv').read_text(),re.M)[1])
        label=next(folder.glob('__*.txt')).stem[2:].replace('_',' ')
        canvas=Image.new('RGB',(1220,690),'white')
        draw=ImageDraw.Draw(canvas)
        draw.text((12,9),f'{run} | {label} | {fps:g} frames/s',fill='black',font=font(19))
        draw.text((12,36),'Blue: inferred gas; orange: disagreement across thresholds 0.60 / 0.70 / 0.80.',fill='black',font=font(16))
        draw.text((12,59),'Unmarked pixels are NOT certified liquid. Dashed red probe: x=496:520, y=18:71.',fill='black',font=font(16))
        for block,index in enumerate(FRAMES):
            file=files[index]
            assert int(re.search(r'(\d+)\.jpg$',file.name)[1])==index
            raw=file.read_bytes()
            a=np.asarray(Image.open(file).convert('L'))
            assert a.shape==(80,1016)
            masks=np.stack([estimator.silhouette(a.astype(float),reference,threshold=t)[0] for t in THRESHOLDS])
            consensus=masks.all(axis=0)
            uncertain=masks.any(axis=0)&~consensus
            row0=92+block*294
            panels=[Image.fromarray(a).convert('RGB'),tint(a,masks[1]),tint(a,consensus,uncertain)]
            names=['Native image','Mask at 0.70','Consensus / disagreement']
            for k,(panel,name) in enumerate(zip(panels,names)):
                top=row0+k*91
                canvas.paste(panel,(194,top))
                draw.text((10,top+8),name,fill='black',font=font(14))
                draw.text((10,top+31),f'Frame {index}',fill='black',font=font(14))
                draw.text((10,top+53),f'{index/fps*1000:.3f} ms',fill='black',font=font(14))
                for x in PROBE_X:
                    for y in range(LUMEN_Y[0],LUMEN_Y[1],7):
                        draw.line((194+x,top+y,194+x,top+min(y+3,70)),fill=(225,35,35),width=1)
                for y in LUMEN_Y:
                    draw.line((194,top+min(y,70),1209,top+min(y,70)),fill=(155,155,155),width=1)
            px=slice(*PROBE_X); py=slice(*LUMEN_Y)
            evidence.append(dict(run=run,pattern=label,frame=index,fps=fps,
                source_file=str(file),source_sha256=hashlib.sha256(raw).hexdigest(),
                source_shape=list(a.shape),native_core_saturation=float(np.mean(a[20:60]>=250)),
                probe_mask_fraction_by_threshold={str(t):float(m[py,px].mean()) for t,m in zip(THRESHOLDS,masks)},
                full_lumen_mask_fraction_by_threshold={str(t):float(m[py].mean()) for t,m in zip(THRESHOLDS,masks)},
                probe_disagreement_fraction=float(uncertain[py,px].mean()),
                lumen_disagreement_fraction=float(uncertain[py].mean()),
                centerpoint_mask_by_threshold={str(t):bool(m[40,508]) for t,m in zip(THRESHOLDS,masks)},
                mask_sha256_by_threshold={str(t):hashlib.sha256(m.tobytes()).hexdigest() for t,m in zip(THRESHOLDS,masks)}))
        canvas.save(BASE/'figures'/f'phase_visual_{run}.png')
    meta=dict(description='Targeted visual phase-identification audit, not independently annotated ground truth.',
        selected_runs=list(RUNS),selected_frames=list(FRAMES),thresholds=list(THRESHOLDS),
        selection='Same fixed frame indices 200 and 210 for eight recordings, including five patterns and two sensitive slug examples.',
        lumen_rows=list(LUMEN_Y),probe_x=list(PROBE_X),probe_point=[508,40],
        detector_sha256=hashlib.sha256((V14/'panorama_void_alternatives.py').read_bytes()).hexdigest(),
        geometry_sha256=hashlib.sha256((V14/'v14_geometry.py').read_bytes()).hexdigest(),
        reference=reference.tolist(),reference_selections=selections,evidence=evidence)
    (BASE/'results'/'phase_visual_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    alt=dispersed_alternative(reference)
    print('Dispersed baseline/isolated projected fractions:',alt['mean_baseline_lumen_fraction'],alt['mean_isolated_lumen_fraction'])
    print('Saved eight raw/mask comparison sheets for 16 fixed native frames.')


if __name__=='__main__':
    run()
