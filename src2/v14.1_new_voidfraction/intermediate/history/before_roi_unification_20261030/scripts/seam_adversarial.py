"""Adversarial checks for the isolated seam correction; no physical-accuracy claim."""
import ast
from collections import defaultdict
import csv
import hashlib
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.optimize import linprog
from scipy.spatial import ConvexHull

from v14_paths import ROOT, BASELINE, DATA, artifact_path, panorama_path, original_panorama_path, write_markdown
from v14_selection import selected_runs, PATTERN_ORDER
import panorama_void_alternatives as estimator
from deep_method_checks import independent_mask, area_by_intervals
from reconstruct_inputs import center_strips, stitch_strips, legacy_stitch, measured_shift


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def table(name):
    with artifact_path(name).open() as stream:
        return list(csv.DictReader(stream))


def reject(action):
    try:
        action()
    except (AssertionError, ValueError):
        return True
    raise AssertionError('An intentionally incorrect control was accepted')


def same_pixels(a,b):
    assert a.shape==b.shape and np.array_equal(a,b), 'Different source-column sequence'


def synthetic_controls():
    strips=[np.array([[1,2,3,4]],dtype=np.uint8),np.array([[5,6,7,8]],dtype=np.uint8)]
    correct=stitch_strips(strips);old=legacy_stitch(strips)
    assert correct.tolist()==[[1,2,3,4,5,6,7,8]]
    assert old.tolist()==[[1,2,3,6,5,6,7,8]]
    flips=lambda a:int(np.count_nonzero((a[:,1:]>5)!=(a[:,:-1]>5)))
    assert flips(correct)==1 and flips(old)==3
    assert reject(lambda:same_pixels(old,correct))
    # Explicit scalar oracle for odd/even and one-pixel widths.
    frames=[(np.arange(2*18).reshape(2,18)+i*40).astype(np.uint8) for i in range(6)]
    widths=[1,2,3,4,5,6]
    expected=np.array([[int(frame[row,18//2-width//2+j]) for frame,width in zip(frames,widths)
                        for j in range(width)] for row in range(2)],dtype=np.uint8)
    actual=stitch_strips(center_strips(frames,widths));same_pixels(actual,expected)
    assert actual.shape==(2,sum(widths))
    bad=actual.copy();bad[0,0]^=1
    assert reject(lambda:same_pixels(bad,expected))
    assert reject(lambda:center_strips(frames,[0]*6))
    assert reject(lambda:center_strips(frames,[19]*6))
    assert reject(lambda:center_strips(frames,[1]))
    for value in (0,127,255):
        constant=[np.full((3,w),value,dtype=np.uint8) for w in widths]
        assert np.all(stitch_strips(constant)==value)
    # Independent known flip denominator; catch W instead of W-1 and wrong H.
    mask=np.zeros((53,8),bool);mask[:,4:]=True
    flips_count=np.count_nonzero(mask[:,1:]!=mask[:,:-1])
    assert flips_count==53 and flips_count/(53*7)==1/7
    assert flips_count/(53*8)!=1/7 and flips_count/(60*7)!=1/7
    whole=np.zeros((80,8),bool);whole[18:71]=True
    assert np.all(area_by_intervals(whole)==1)
    assert np.all(area_by_intervals(whole,60)<1)
    return dict(monotone_source_columns=correct.tolist()[0],legacy_columns=old.tolist()[0],
                correct_flips=flips(correct),legacy_flips=flips(old),
                odd_even_one_pixel_widths_passed=True,constant_pixels_preserved=True,
                legacy_and_one_pixel_mutations_rejected=True,invalid_widths_rejected=True,
                wrong_flip_denominators_and_D60_detected=True)


def preservation_checks():
    class AllowedInputAndDocumentationChanges(ast.NodeTransformer):
        def visit_ImportFrom(self,node):
            return None if node.module=='v14_paths' else node
        def visit_Constant(self,node):
            if isinstance(node.value,str):
                node.value=node.value.replace('src2/v14_new_voidfraction','src2/v14.1_new_voidfraction')
                node.value=node.value.replace('src2/data_processed/*/__panorama_v11.png',
                                              'intermediate/reconstructed_inputs/*/__panorama_v11.png')
            return node
        def visit_Call(self,node):
            node=self.generic_visit(node)
            if (isinstance(node.func,ast.Attribute) and isinstance(node.func.value,ast.Name)
                    and node.func.value.id=='Image' and node.func.attr=='open' and node.args):
                path=node.args[0]
                if (isinstance(path,ast.BinOp) and isinstance(path.op,ast.Div)
                        and isinstance(path.right,ast.Constant) and path.right.value=='__panorama_v11.png'):
                    run=path.left
                    if (isinstance(run,ast.BinOp) and isinstance(run.op,ast.Div)
                            and isinstance(run.left,ast.Name) and run.left.id=='DATA'):
                        run=run.right
                    node.args[0]=ast.Call(func=ast.Name(id='panorama_path',ctx=ast.Load()),args=[run],keywords=[])
            return node
    modules=['panorama_void_alternatives.py','phase_mask_flip_analysis.py','raw_flip_analysis.py',
             'image_void_alternatives.py','pattern_identification.py','silhouette_pattern_identification.py']
    for filename in modules:
        normalize=lambda path:ast.dump(AllowedInputAndDocumentationChanges().visit(ast.parse(path.read_text())),include_attributes=False)
        assert normalize(artifact_path(filename))==normalize(BASELINE/'scripts'/filename),f'Unrequested pipeline change: {filename}'
    functions={
        'panorama_void_alternatives.py':['normalize','make_reference','close_shadow_crop','silhouette','geometry'],
        'raw_flip_analysis.py':['otsu','state_statistics','config','two_sided_mann_whitney'],
        'pattern_identification.py':['features','predict','evaluate','training_fold_metadata','metrics'],
        'silhouette_pattern_identification.py':['interval_lengths','extract','checks'],
        'summary_figures.py':['hull_indices','alternation_chart'],
    }
    checked=[]
    for filename,names in functions.items():
        current=ast.parse(artifact_path(filename).read_text())
        previous=ast.parse((BASELINE/'scripts'/filename).read_text())
        lookup=lambda tree:{n.name:ast.dump(n,include_attributes=False) for n in tree.body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef))}
        a,b=lookup(current),lookup(previous)
        for name in names:
            assert a[name]==b[name],f'Unrequested method change: {filename}:{name}'
            checked.append(filename+':'+name)
    for filename in ('v14_geometry.py','v14_selection.py'):
        assert artifact_path(filename).read_bytes()==(BASELINE/'scripts'/filename).read_bytes(),filename
    for filename in functions:
        assignments=lambda path:{n.targets[0].id:ast.dump(n.value,include_attributes=False) for n in ast.parse(path.read_text()).body
                                  if isinstance(n,ast.Assign) and len(n.targets)==1 and isinstance(n.targets[0],ast.Name)
                                  and n.targets[0].id.isupper() and n.targets[0].id not in {'HERE','ROOT','DATA','CACHE'}}
        assert assignments(artifact_path(filename))==assignments(BASELINE/'scripts'/filename),f'Changed constants: {filename}'
    assert estimator.THRESHOLD==.70
    return dict(entire_analysis_modules_unchanged_except_input_routing_and_documented_paths=modules,
                unchanged_function_AST=checked,geometry_and_selection_byte_identical=True,analysis_constants_unchanged=True)


def source_column_checks():
    shifts=defaultdict(list)
    for row in table('reconstruction_shifts.csv'):
        shifts[row['run']].append(row)
    reconstruction=json.loads(artifact_path('input_reconstruction.json').read_text())
    recorded={r['run']:r for r in reconstruction['runs']}
    assert set(shifts)==set(recorded)==set(selected_runs())
    checks=[]
    for run,rows in shifts.items():
        paths=sorted((DATA/run).glob('ImgA*.jpg'))
        assert [p.name for p in paths]==[r['frame'] for r in rows] and len(paths)==501
        with Image.open(panorama_path(run)) as image:current=np.asarray(image.convert('L'))
        with Image.open(original_panorama_path(run)) as image:prior=np.asarray(image.convert('L'))
        assert current.shape==prior.shape
        cursor=0
        for i,(path,row) in enumerate(zip(paths,rows)):
            width=int(row['strip_width_px']);raw=int(row['raw_shift_px'])
            assert width==max(raw,1)
            assert (row['last_shift_reused']=='True')==(i==len(rows)-1)
            with Image.open(path) as image:frame=np.asarray(image.convert('L'))
            start=frame.shape[1]//2-width//2
            # No import/call to center_strips or stitch_strips in this real-data oracle.
            assert np.array_equal(current[:,cursor:cursor+width],frame[:,start:start+width]),(run,i)
            cursor+=width
        assert rows[-1]['raw_shift_px']==rows[-2]['raw_shift_px']
        assert cursor==current.shape[1]==recorded[run]['width_px']
        assert sha(panorama_path(run))==recorded[run]['corrected_panorama_sha256']
        assert sha(original_panorama_path(run))==recorded[run]['original_panorama_sha256']
        # Every changed pixel must fall in the one column actually overwritten by v11.
        permitted=np.zeros(current.shape[1],bool)
        widths=[int(r['strip_width_px']) for r in rows]
        for i,boundary in enumerate(np.cumsum(widths)[:-1],1):
            if min(2,min(widths[i-1],widths[i])//2)==2:permitted[boundary-1]=True
        changed=np.any(current!=prior,axis=0)
        assert not (changed & ~permitted).any(),f'Unexpected input modification: {run}'
        assert int(np.count_nonzero(current!=prior))==recorded[run]['changed_pixels']
        if changed.any():assert reject(lambda:same_pixels(prior,current))
        checks.append(dict(run=run,native_strips=len(rows),width=cursor,
                          changed_columns=int(changed.sum()),changes_confined_to_defective_seam_columns=True))
    sources=json.loads(artifact_path('source_input_hashes.json').read_text())
    assert all(sha(ROOT.parent/p)==h for p,h in sources.items()),'Native data or shared pipeline changed'
    frozen=json.loads((ROOT/'intermediate/history/v14_baseline_manifest.json').read_text())
    original=Path(frozen['original_root'])
    for relative,digest in frozen['files'].items():
        assert sha(BASELINE/relative)==sha(original/relative)==digest,f'Baseline/original changed: {relative}'
    return dict(runs=checks,native_strips_checked=sum(r['native_strips'] for r in checks),
                source_files_unchanged=len(sources),baseline_and_original_files_unchanged=len(frozen['files']))


def independent_intersections(points):
    polygons={}
    for pattern in PATTERN_ORDER:
        values=np.array([[float(r['alpha']),float(r['F1000'])] for r in points if r['pattern']==pattern])
        if len(values)>=3:polygons[pattern]=values[ConvexHull(values).vertices]
    result=[]
    # Convex combinations instead of compare_versions' half-space inequalities.
    for a,b in combinations(polygons,2):
        x,y=polygons[a],polygons[b]
        equality=np.vstack([np.hstack([x.T,-y.T]),np.r_[np.ones(len(x)),np.zeros(len(y))],
                            np.r_[np.zeros(len(x)),np.ones(len(y))]])
        solved=linprog(np.zeros(len(x)+len(y)),A_eq=equality,b_eq=[0,0,1,1],bounds=(0,None),method='highs')
        if solved.success:result.append([a,b])
        else:assert solved.status==2,solved.message
    # Also check without an optimizer, so shared LP feasibility tolerance cannot
    # manufacture a reported intersection between narrowly separated hulls.
    def cross(a,b,c):
        return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
    def inside(point,polygon):
        signs=[cross(a,b,point) for a,b in zip(polygon,np.roll(polygon,-1,axis=0))]
        return min(signs)>=0 or max(signs)<=0
    def crosses(a,b,c,d):
        if any(max(min(a[k],b[k]),min(c[k],d[k]))>min(max(a[k],b[k]),max(c[k],d[k])) for k in (0,1)):
            return False
        return cross(a,b,c)*cross(a,b,d)<=0 and cross(c,d,a)*cross(c,d,b)<=0
    geometric=[]
    for a,b in combinations(polygons,2):
        x,y=polygons[a],polygons[b]
        if (any(inside(p,y) for p in x) or any(inside(p,x) for p in y) or
                any(crosses(p,q,r,s) for p,q in zip(x,np.roll(x,-1,axis=0))
                    for r,s in zip(y,np.roll(y,-1,axis=0)))):
            geometric.append([a,b])
    assert result==geometric,'LP tolerance changed a hull-intersection decision'
    return result


def threshold_checks():
    reference,_=estimator.make_reference()
    coordinates=table('threshold_plot_data.csv')
    lookup={(r['run'],float(r['shadow_threshold'])):r for r in coordinates if r['condition']=='v14.1'}
    alphas={(r['run'],float(r['variant'].removeprefix('shadow_ratio_'))):float(r['alpha_conditional'])
            for r in table('panorama_void_sensitivity.csv') if r['variant'].startswith('shadow_ratio_')}
    flips={(r['run'],float(r['shadow_threshold'])):r for r in table('phase_mask_flip_sensitivity.csv')}
    checked=0
    for run in selected_runs():
        image=estimator.read(run)
        for threshold in (.60,.65,.70,.75,.80):
            mask=independent_mask(image,reference,threshold=threshold)
            alpha=float(area_by_intervals(mask).mean())
            rate=float(np.count_nonzero(mask[18:71,1:]^mask[18:71,:-1])/(53*(mask.shape[1]-1)))
            assert abs(alpha-alphas[run,threshold])<1e-12
            assert abs(rate-float(flips[run,threshold]['candidate_pixel_flip_rate']))<1e-15
            if (run,threshold) in lookup:
                assert abs(alpha-float(lookup[run,threshold]['alpha']))<1e-12
                assert abs(rate*1000-float(lookup[run,threshold]['F1000']))<1e-12
            checked+=1
    intersections=table('threshold_hull_intersections.csv')
    assert len(intersections)==10
    for row in intersections:
        subset=[p for p in coordinates if p['condition']==row['condition'] and float(p['shadow_threshold'])==float(row['threshold'])]
        assert len(subset)==int(row['n_acquisitions'])==25
        assert independent_intersections(subset)==json.loads(row['overlap_pairs'])
    return dict(independent_mask_area_flip_threshold_cases=checked,
                hull_intersections_verified_with_distinct_linear_program=10,
                hull_intersections_also_verified_without_optimizer=10,
                thresholds_preserved=[.60,.65,.70,.75,.80],primary_threshold=.70)


def comparison_checks():
    rows=table('fix_comparison_per_run.csv')
    assert {r['run'] for r in rows}==set(selected_runs())
    saved=json.loads(artifact_path('fix_effectiveness.json').read_text())
    selected=[r for r in rows if r['eligible']=='True']
    for row in rows:
        for prefix,delta in [('alpha','alpha_change'),('F1000','F1000_change')]:
            assert abs(float(row[prefix+'_v141'])-float(row[prefix+'_v14'])-float(row[delta]))<1e-12
        assert abs(100*float(row['F1000_change'])/float(row['F1000_v14'])-float(row['flip_relative_change_percent']))<1e-12
    assert max(abs(float(r['alpha_change'])) for r in selected)==saved['max_absolute_alpha_change']
    assert min(float(r['flip_relative_change_percent']) for r in selected)==saved['min_relative_flip_change_percent']
    assert max(float(r['flip_relative_change_percent']) for r in selected)==saved['max_relative_flip_change_percent']
    classifiers=table('classifier_comparison.csv')
    fixed=[r for r in classifiers if r['variant']=='baseline' and r['validation'] in ('leave_one_run_out','leave_one_mass_flux_out')]
    for row in classifiers:
        assert abs(float(row['balanced_accuracy_v141'])-float(row['balanced_accuracy_v14'])-float(row['delta_balanced_accuracy']))<1e-15
    assert len(fixed)==saved['fixed_classifier_comparisons']
    assert sum(float(r['delta_balanced_accuracy'])>1e-12 for r in fixed)==saved['fixed_classifiers_improved']
    assert sum(float(r['delta_balanced_accuracy'])< -1e-12 for r in fixed)==saved['fixed_classifiers_worsened']
    assert sum(abs(float(r['delta_balanced_accuracy']))<=1e-12 for r in fixed)==saved['fixed_classifiers_unchanged']
    for filename in ('EXECUTIVE_SUMMARY.md','EXECUTIVE_SUMMARY_PT_BR.md'):
        text=artifact_path(filename).read_text()
        for name in ('max_absolute_alpha_change','min_relative_flip_change_percent','max_relative_flip_change_percent'):
            formatted=f'{saved[name]:.6f}' if name=='max_absolute_alpha_change' else f'{saved[name]:.2f}'
            assert formatted in text or formatted.replace('.',',') in text
        assert 'threshold_sensitivity.png' in text and 'fix_comparison_per_run.csv' in text
    return dict(acquisitions=len(rows),classifier_variants=len(classifiers),fixed_classifier_comparisons=len(fixed),
                report_change_numbers_checked=True,all_outcomes_reported=True)


def main():
    jobs=[('synthetic_and_mutation_controls',synthetic_controls),('method_preservation',preservation_checks),
          ('source_column_order_and_immutable_inputs',source_column_checks),('all_thresholds_and_hulls',threshold_checks),
          ('honest_comparison_and_reports',comparison_checks)]
    checks=[]
    for name,job in jobs:
        try:checks.append(dict(check=name,passed=True,details=job()))
        except Exception as error:checks.append(dict(check=name,passed=False,error=f'{type(error).__name__}: {error}'))
    passed=all(r['passed'] for r in checks)
    result=dict(passed=passed,checks=checks,script_sha256=sha(__file__),
                limitation='These controls verify implementation and reporting, not true gas phase, refraction, advection validity or physical void accuracy.')
    artifact_path('seam_adversarial_checks.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    write_markdown(artifact_path('seam_adversarial_review.md'),
        '# Adversarial review of the v14.1 correction\n\n'+f'**{"PASS" if passed else "FAIL"}.**\n\n'+
        '\n'.join(f'- {r["check"]}: {"PASS" if r["passed"] else r["error"]}' for r in checks)+
        '\n\nThe tests reject the old reordered ramp and a deliberate one-pixel mutation, verify every native strip against the corrected output, '
        'and restrict real input changes to the defective seam columns. Estimator, classifier, geometry and selection preservation is checked against the frozen v14 source. '
        'All 130 threshold cases are recomputed by independent boolean segmentation/interval arithmetic; ten before/after hull comparisons use both a different linear-program formulation and direct polygon intersection without an optimizer. '
        'Classifier degradations remain in the comparison. Original files are hash-verified.\n\n'
        'The correction does not establish better physical accuracy or remove sensitivity to segmentation and reconstruction assumptions. '
        '[Detailed evidence](seam_adversarial_checks.json) · [Effectiveness](fix_effectiveness.md)\n')
    print(json.dumps(result,indent=2))
    return 0 if passed else 1


if __name__=='__main__':raise SystemExit(main())
