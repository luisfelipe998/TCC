"""Check the unchanged D53 geometry on the newly reconstructed v14.1 masks."""
import csv
import hashlib
import json
import numpy as np
from PIL import Image
from v14_paths import artifact_path, panorama_path, write_markdown
from v14_selection import selected_runs
from v14_geometry import DIAMETER_PX, LUMEN_ROWS, BOUNDARY_PROVENANCE


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_normalization(expected_runs=None):
    from deep_method_checks import area_by_intervals
    assert DIAMETER_PX == 53 and LUMEN_ROWS == (18,71)
    names=set(selected_runs()) if expected_runs is None else set(expected_runs)
    with artifact_path('panorama_void_per_run.csv').open() as stream:
        values={r['run']:r for r in csv.DictReader(stream)}
    with artifact_path('data_inventory.csv').open() as stream:
        inventory={r['run']:r for r in csv.DictReader(stream)}
    assert set(values)==names
    assert set(inventory)==names
    records=[]
    for run in sorted(names):
        assert sha(panorama_path(run))==inventory[run]['sha256'],run
        mask=np.asarray(Image.open(artifact_path('phase_masks')/f'{run}_mask.png'))==255
        assert not mask[:18].any() and not mask[71:].any()
        area53=area_by_intervals(mask,53)
        area60=area_by_intervals(mask,60)
        assert np.all((area53>=0)&(area53<=1))
        error=float(np.max(np.abs(area53-area60*(60/53)**2)))
        assert error<1e-14
        assert abs(area53.mean()-float(values[run]['conditional_circle_alpha_spatial_mean']))<1e-14
        records.append(dict(run=run,current_alpha_D53=float(area53.mean()),
            counterfactual_same_mask_D60=float(area60.mean()),max_rescaling_error=error,
            input_panorama_sha256=sha(panorama_path(run))))
    result=dict(passed=True,selected_acquisitions=len(records),current_masks_checked=len(records),
        corrected_input_hashes_checked=len(records),confirmed_lumen_rows=[18,71],lumen_diameter_px=53,
        full_lumen_maximum=1,alpha_scale_factor=(60/53)**2,max_rescaling_error=max(r['max_rescaling_error'] for r in records),
        boundary_source=BOUNDARY_PROVENANCE,
        segmentation_source_sha256=sha(artifact_path('panorama_void_alternatives.py')),
        geometry_source_sha256=sha(artifact_path('v14_geometry.py')),failures=[],
        interpretation='D53 is unchanged between v14 and v14.1. D60 is only an algebraic control on each current mask; masks and flip rates may change because the panorama input changed.')
    return result,records


def main():
    result,records=verify_normalization()
    artifact_path('lumen_normalization_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    with artifact_path('lumen_normalization_changes.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    write_markdown(artifact_path('lumen_normalization_review.md'),
        '# Verificação da geometria preservada na v14.1\n\n'
        '**PASS:** fatia 18:71, 53 linhas e D=53 preservados. Uma máscara cheia produz α=1. '
        'Para cada máscara atual, verificou-se também α53=α60(60/53)². '
        'D60 é somente um controle algébrico, não uma nova estimativa.\n\n'
        'As mudanças entre v14 e v14.1 vêm da entrada corrigida e da referência reconstruída. '
        'Não se afirma que as máscaras ou os flip-rates sejam iguais entre versões.\n\n'
        '[Verificação](lumen_normalization_checks.json) · [Controle por aquisição](lumen_normalization_changes.csv)\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
