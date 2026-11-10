#!/usr/bin/env python3
"""Verify final plot coordinates from SVG ticks, coverage, links and provenance."""
from pathlib import Path
import ast,hashlib,json,re,sys
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from PIL import Image
BASE=Path(__file__).resolve().parents[1]
NS={'s':'http://www.w3.org/2000/svg'}


def verify_svg(path,data):
    tree=ET.parse(path)
    root=tree.getroot()
    axes=root.find('.//s:g[@id="axes_1"]',NS)
    fits=[]
    for prefix,coordinate in [('xtick_','x'),('ytick_','y')]:
        coords=[];values=[]
        for element in axes.findall('.//s:g',NS):
            if not element.get('id','').startswith(prefix):continue
            t=element.find('.//s:text',NS)
            u=element.find('.//s:use',NS)
            if t is None or u is None:continue
            try:value=float(''.join(t.itertext()))
            except ValueError:continue
            coords.append(float(u.attrib[coordinate]));values.append(value)
        assert len(coords)>=4
        fits.append(np.polyfit(coords,values,1))
    target=np.abs(data[['plot_core_signed_px_s','plot_wall_signed_px_s']].to_numpy())
    found=[];filled=[];errors=[]
    for g in axes.findall('s:g',NS):
        if not g.get('id','').startswith('PathCollection_'):continue
        u=g.find('.//s:use',NS)
        assert u is not None
        point=np.array([np.polyval(fits[0],float(u.attrib['x'])),np.polyval(fits[1],float(u.attrib['y']))])*1000
        distance=np.max(abs(target-point),axis=1)
        at=int(np.argmin(distance));error=float(distance[at])
        assert error<.003,(path.name,point,error)
        style=u.get('style','')
        # The extra black x is a coverage warning at the first candidate.
        if 'stroke: #000000' in style:
            assert not bool(data.iloc[at].internal_supported)
            continue
        white='fill: #ffffff' in style or 'fill: none' in style
        assert (not white)==bool(data.iloc[at].stable_all_variants),(path.name,at,style)
        found.append(at);filled.append(not white);errors.append(error)
    assert sorted(found)==list(range(26))
    text=' '.join(''.join(t.itertext()) for t in root.findall('.//s:text',NS))
    for run in data.run:assert run in text
    assert '10³ px/s' in text
    return {'points':len(found),'filled':sum(filled),'max_coordinate_error_px_s':max(errors)}


def main():
    sys.path.insert(0,str(BASE.parent/'scripts'))
    from publish_latest_summary import check_publication
    main_reports=check_publication()
    data=pd.read_csv(BASE/'results/final_all_26_plot_data.csv')
    source=pd.read_csv(BASE/'results/split_background_pairs.csv')
    source=source[(source.variant=='training_background')&(source.selection_direction=='forward')].set_index('run')
    assert len(data)==26 and data.run.nunique()==26
    old=pd.read_csv(BASE.parent/'results/all_run_status.csv').set_index('run')
    assert set(data.run)==set(old.index)
    assert all(data.pattern==data.run.map(old.pattern))
    for r in data.itertuples():
        for field in ('core_signed_px_s','wall_signed_px_s'):
            assert np.isclose(getattr(r,field),source.loc[r.run,field],rtol=1e-12)
        assert r.internal_supported==source.loc[r.run,'supported']
        assert r.coordinate_source=='jointly_accepted_assessment_blocks'
        assert not r.material_phase_pair
    meta=json.loads((BASE/'results/final_summary_metrics.json').read_text())
    ledger=pd.read_csv(BASE/'results/final_method_coverage.csv')
    for field,value in meta['counts'].items():assert int(ledger[field].sum())==value
    # Reconstruct sensitivity flags from paired candidate files, not plotting code.
    split=pd.read_csv(BASE/'results/split_background_pairs.csv')
    original=pd.read_csv(BASE/'results/selected_wall_per_run.csv')
    stable={}
    for r in data.itertuples():
        accepted=bool(r.internal_supported)
        for variant,direction,frame in [('training_background','reverse',split),
            ('training_background_hp32','forward',split),('training_background_hp8','forward',split),
            ('raw_blocks10','forward',original)]:
            row=frame[(frame.run==r.run)&(frame.variant==variant)&(frame.selection_direction==direction)].iloc[0]
            accepted &= bool(row.supported)
            for field in ('core_signed_px_s','wall_signed_px_s'):
                ref=getattr(r,field);accepted &= bool(abs(row[field]-ref)<=max(2*r.fps,.2*abs(ref)))
        stable[r.run]=accepted
        assert accepted==r.stable_all_variants
    # Exact 3-block/20%-tolerance ablation reference must reproduce the final gate.
    gates=pd.read_csv(BASE/'results/reporting_gate_sensitivity.csv')
    gate=gates[gates.primary_setting].set_index('run')
    assert all(data.internal_supported==data.run.map(gate.supported))
    figures={}
    for language in ('PT_BR','EN'):
        stem='all_26_selected_wall_'+language
        figures[language]=verify_svg(BASE/'figures'/f'{stem}.svg',data)
        im=Image.open(BASE/'figures'/f'{stem}.png')
        assert min(im.info['dpi'])>599
        assert (BASE/'figures'/f'{stem}.pdf').read_bytes().startswith(b'%PDF-')
    links=[]
    for path in BASE.rglob('*.md'):
        for target in re.findall(r'\]\(([^)]+)\)',path.read_text()):
            if target.startswith(('http:','https:','#','mailto:')):continue
            local=target.split('#')[0]
            assert (path.parent/local).exists(),(str(path),target)
            links.append((str(path.relative_to(BASE)),target))
    for path in (BASE/'scripts').glob('*.py'):ast.parse(path.read_text())
    audits={}
    for filename in ('independent_block_verification.json','phase_local_verification.json','spectral_verification.json'):
        obj=json.loads((BASE/'results'/filename).read_text())
        audits[filename]=obj.get('status')
        assert 'pass' in str(obj.get('status')).lower(),(filename,obj.get('status'))
    manifests={}
    for path in sorted(BASE.rglob('*')):
        if not path.is_file() or 'intermediate' in path.parts or '__pycache__' in path.parts:continue
        if path.name in ('final_verification.json','artifact_manifest.json'):continue
        manifests[str(path.relative_to(BASE))]=hashlib.sha256(path.read_bytes()).hexdigest()
    (BASE/'results/artifact_manifest.json').write_text(json.dumps(manifests,indent=2))
    result={'status':'passed','all_points':26,'internally_supported':int(data.internal_supported.sum()),
        'stable_all_variants':sum(stable.values()),'material_phase_pairs':0,'svg_coordinate_checks':figures,
        'local_links_checked':len(links),'audits':audits,'manifest_files':len(manifests),
        'main_reports_match_latest':main_reports,
        'limitations':'Numerical/provenance/presentation verification, not calibration or experimental validation of material phase velocities.'}
    (BASE/'results/final_verification.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
