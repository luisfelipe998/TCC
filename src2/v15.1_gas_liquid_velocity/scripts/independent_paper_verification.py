#!/usr/bin/env python3
"""Independently verify paired variants, bilingual data, and rendered SVG positions."""
from pathlib import Path
import hashlib
import json
import math
import re
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from PIL import Image

BASE=Path(__file__).resolve().parents[1]
OUT=BASE/'results'
NS={'s':'http://www.w3.org/2000/svg'}
INPUTS={}
PT={'dispersed bubbles':'Bolhas dispersas','plug':'Plug','slug':'Slug','churn':'Churn','wavy annular':'Anular ondulado'}


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def remember(path):
    INPUTS[str(path.relative_to(BASE))]=sha(path)
    return path


def table(name): return pd.read_csv(remember(OUT/name))


def numeric(text,pt):
    text=text.strip().replace('−','-')
    return float(text.replace('.','').replace(',','.') if pt else text.replace(',',''))


def svg_points(path,pt):
    doc=ET.parse(remember(path)).getroot()
    axes=doc.find(".//s:g[@id='axes_1']",NS)
    assert axes is not None
    transforms=[]
    for axis in ('x','y'):
        positions=[]
        for group in axes.iter('{'+NS['s']+'}g'):
            if not group.attrib.get('id','').startswith(axis+'tick_'): continue
            label=group.find('.//s:text',NS)
            tick=group.find('.//s:use',NS)
            assert label is not None and tick is not None
            positions.append((float(tick.attrib[axis]),numeric(''.join(label.itertext()),pt)))
        assert len(positions)>=2
        coord,values=np.array(positions).T
        slope=(values[-1]-values[0])/(coord[-1]-coord[0])
        intercept=values[0]-slope*coord[0]
        assert np.allclose(values,coord*slope+intercept,atol=.002,rtol=1e-8)
        transforms.append((slope,intercept))
    points=[]
    for group in axes.findall('s:g',NS):
        if not group.attrib.get('id','').startswith('PathCollection_'): continue
        for marker in group.findall('.//s:use',NS):
            x,y=float(marker.attrib['x']),float(marker.attrib['y'])
            points.append((transforms[0][0]*x+transforms[0][1],transforms[1][0]*y+transforms[1][1]))
        # Matplotlib emits hollow polygons as absolute paths rather than uses.
        for marker in group.findall('s:path',NS):
            coordinates=np.array([float(v) for v in re.findall(r'[+-]?[0-9]+(?:[.][0-9]*)?(?:[eE][+-]?[0-9]+)?',marker.attrib['d'])]).reshape(-1,2)
            x,y=(coordinates.min(axis=0)+coordinates.max(axis=0))/2
            points.append((transforms[0][0]*x+transforms[0][1],transforms[1][0]*y+transforms[1][1]))
    text=' '.join(''.join(node.itertext()) for node in doc.findall('.//s:text',NS))
    return np.array(points),text


def check_scatter(stem,expected,unit,extra_markers=0):
    results=[]
    for pt,language in [(False,'EN'),(True,'PT_BR')]:
        path=BASE/'figures'/f'{stem}_{language}.svg'
        actual,text=svg_points(path,pt)
        assert len(actual)==len(expected)+extra_markers,(stem,len(actual),len(expected))
        assert unit in text
        if unit=='px/s': assert 'm/s' not in text
        if 'conditional' in stem:
            assert ('ESCALA CONDICIONAL' if pt else 'CONDITIONAL SCALE') in text
            assert '53' in text
        remaining=list(range(len(actual)))
        max_error=0.
        for point in expected:
            distances=[float(np.max(np.abs(actual[k]-point))) for k in remaining]
            index=int(np.argmin(distances));error=distances[index]
            assert error<(.003 if unit=='px/s' else 2e-7),(stem,language,point,error)
            max_error=max(max_error,error)
            remaining.pop(index)
        png=remember(BASE/'figures'/f'{stem}_{language}.png')
        pdf=remember(BASE/'figures'/f'{stem}_{language}.pdf')
        with Image.open(png) as im:
            assert min(im.info.get('dpi',(0,0)))>=599
            assert im.width>2500 and im.height>2000
        assert pdf.read_bytes().startswith(b'%PDF')
        results.append(dict(stem=stem,language=language,points=len(expected),extra_markers=extra_markers,
                            max_svg_coordinate_difference=max_error,units=unit))
    return results


def main():
    region=table('annular_region_summary.csv').set_index(['run','region'])
    variants=table('annular_variant_pairs.csv')
    spatial=table('annular_spatial_pairs.csv').set_index('run')
    passage=table('passage_per_run.csv').set_index('run')
    assert len(variants)==61 and not variants.duplicated(['run','variant']).any()
    checked_variant_fields=0
    for r in variants.itertuples():
        c,u,l=(region.loc[(r.run,name)] for name in (r.core_region,r.upper_region,r.lower_region))
        wall=(u.signed_velocity_px_frame+l.signed_velocity_px_frame)/2
        agree=math.isfinite(wall) and np.sign(u.signed_velocity_px_frame)==np.sign(l.signed_velocity_px_frame) and abs(u.signed_velocity_px_frame-l.signed_velocity_px_frame)<=max(2.,.25*abs(wall))
        direction=math.isfinite(wall) and abs(wall)>=1 and np.sign(wall)==np.sign(c.signed_velocity_px_frame)
        baseline=bool(c.ensemble_supported and u.ensemble_supported and l.ensemble_supported and agree and direction)
        final=bool(c.robust_ensemble_supported and u.robust_ensemble_supported and l.robust_ensemble_supported and agree and direction)
        assert bool(r.baseline_supported)==baseline and bool(r.recommended_supported)==final
        assert bool(r.upper_lower_agreement)==agree and bool(r.core_wall_same_direction)==direction
        expected=[c.signed_velocity_px_frame,wall,c.signed_velocity_px_s,wall*r.fps,u.signed_velocity_px_s,l.signed_velocity_px_s,
                  c.signed_velocity_px_frame/wall if abs(wall)>=1 else math.nan]
        actual=[r.core_signed_px_frame,r.nearwall_signed_px_frame,r.core_signed_px_s,r.nearwall_signed_px_s,
                r.upper_signed_px_s,r.lower_signed_px_s,r.core_to_nearwall_celerity_ratio]
        assert np.allclose(actual,expected,atol=1e-8,rtol=1e-12,equal_nan=True)
        checked_variant_fields+=11
    variant_counts={}
    for name,g in variants.groupby('variant'):
        variant_counts[name]=dict(total=len(g),baseline=int(g.baseline_supported.sum()),final=int(g.recommended_supported.sum()),
                                 annular_final=int(((g.pattern=='wavy annular')&g.recommended_supported).sum()))
    assert variant_counts['inner_bands']['final']==10
    assert variant_counts['outer_bands']['final']==12
    assert variant_counts['gradient_locations']['final']==0
    pplot=table('paper_passage_plot_data.csv').set_index('run')
    splot=table('paper_regional_plot_data.csv').set_index('run')
    oplot=table('paper_regional_outer_plot_data.csv').set_index('run')
    allplot=table('paper_all_candidate_plot_data.csv').set_index('run')
    assert set(pplot.index)==set(passage[passage.eligible_passage_comparison].index)
    assert set(splot.index)==set(spatial[spatial.robust_spatial_descriptor].index)
    accepted_outer=variants[(variants.variant=='outer_bands')&variants.recommended_supported].set_index('run')
    assert set(oplot.index)==set(accepted_outer.index) and len(oplot)==12
    assert set(allplot.index)==set(spatial.index) and list(allplot.figure_id)==list(range(1,27))
    for name,column in [('gas','gas_plot_magnitude_px_s'),('liquid','liquid_plot_magnitude_px_s')]:
        np.testing.assert_allclose(pplot[column],abs(passage.loc[pplot.index,name+'_mean_signed_px_s']),rtol=1e-13)
        assert (passage.loc[pplot.index,name+'_mean_signed_px_s']<0).all()
    for name in ('core','nearwall'):
        np.testing.assert_allclose(splot[name+'_magnitude_px_s'],abs(spatial.loc[splot.index,name+'_signed_px_s']),rtol=1e-13)
        np.testing.assert_allclose(oplot[name+'_magnitude_px_s'],abs(accepted_outer.loc[oplot.index,name+'_signed_px_s']),rtol=1e-13)
        assert (accepted_outer.loc[oplot.index,name+'_signed_px_s']<0).all()
        np.testing.assert_allclose(allplot[name+'_signed_px_s'],spatial.loc[allplot.index,name+'_signed_px_s'],rtol=1e-13)
        assert (allplot[name+'_signed_px_s']<0).all()
    assert np.array_equal(allplot.robust_spatial_descriptor,spatial.loc[allplot.index].robust_spatial_descriptor)
    figures=[]
    points=pplot[['gas_plot_magnitude_px_s','liquid_plot_magnitude_px_s']].to_numpy()
    figures+=check_scatter('gas_liquid_passages_px_s',points,'px/s')
    figures+=check_scatter('gas_liquid_passages_conditional_m_s',points*.001/53,'m/s')
    figures+=check_scatter('regional_motion_px_s',splot[['core_magnitude_px_s','nearwall_magnitude_px_s']].to_numpy(),'px/s')
    figures+=check_scatter('regional_outer_motion_px_s',oplot[['core_magnitude_px_s','nearwall_magnitude_px_s']].to_numpy(),'px/s')
    figures+=check_scatter('all_regional_candidates',abs(allplot[['core_signed_px_s','nearwall_signed_px_s']].to_numpy()),'px/s',extra_markers=1)
    status=table('all_run_status.csv').set_index('run')
    assert set(status.index)==set(passage.index) and len(status)==26
    assert np.array_equal(status.passage_comparison,passage.loc[status.index].eligible_passage_comparison)
    assert np.array_equal(status.regional_comparison,spatial.loc[status.index].robust_spatial_descriptor)
    assert not status.gas_liquid_material_pair.any() and not status.superficial_velocity_pair.any()
    assert not status.regional_reasons.fillna('').str.contains('nan').any()
    report_cells=0
    report_documents=[]
    for pt,suffix in [(False,''),(True,'_PT_BR')]:
        archived=BASE/'findings'/f'EXECUTIVE_SUMMARY{suffix}_2026-10-23.md'
        path=remember(archived if archived.exists() else BASE/'findings'/f'EXECUTIVE_SUMMARY{suffix}.md')
        report_documents.append(str(path.relative_to(BASE)))
        document=path.read_text()
        seen=set()
        for line in document.splitlines():
            cells=[x.strip() for x in line.strip('|').split('|')]
            if len(cells)==6 and cells[0] in pplot.index:
                r=passage.loc[cells[0]]
                assert cells[1]==(PT[r.pattern] if pt else r.pattern)
                assert numeric(cells[2],pt)==round(r.gas_mean_signed_px_s)
                assert numeric(cells[3],pt)==round(r.liquid_mean_signed_px_s)
                assert numeric(cells[4],pt)==round(r.ratio_of_signed_means,3)
                assert cells[5]==f'{r.gas_pairs}/{r.liquid_pairs}'
                seen.add(cells[0]);report_cells+=4
        assert seen==set(pplot.index)
        assert 'px/s' in document and '53' in document
    for name,digest in INPUTS.items(): assert sha(BASE/name)==digest,('Inputs changed',name)
    result=dict(status='passed',paired_variant_rows=len(variants),checked_variant_fields=checked_variant_fields,
                variant_counts=variant_counts,plot_data_rows=dict(passage=len(pplot),regional=len(splot),regional_outer=len(oplot),all_candidates=len(allplot)),
                svg_checks=figures,bilingual_report_numeric_cells=report_cells,all_run_status_rows=26,input_sha256=INPUTS,
                report_documents=report_documents,report_scope='Original passage and regional results; newer all-recording results have their own independent verification.',
                verifier_sha256=sha(Path(__file__)),
                independence='No generator imported. Variant gates reconstructed from region rows; exported SVG marker positions mapped back to data using axis ticks in each language.',
                limitations='Checks numerical plotting, scope, support selection and figure units, not experimental material velocity or physical calibration. Annotation readability is also reviewed visually.')
    (OUT/'independent_paper_verification.json').write_text(json.dumps(result,indent=2)+chr(10))
    print(json.dumps({k:v for k,v in result.items() if k!='input_sha256'},indent=2))


if __name__=='__main__': main()
