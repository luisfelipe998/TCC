#!/usr/bin/env python3
"""Audit native timing and spatial-calibration evidence without estimating motion."""
from pathlib import Path
from collections import Counter
import configparser
import csv
import hashlib
import json
import math
import re
import runpy
import subprocess

from PIL import Image, ExifTags
import openpyxl

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[1]
DATA = ROOT / 'src2/data_processed'
PRIOR = ROOT / 'src2/v15_gas_liquid_velocity/results/scenario_inventory.csv'
GEOMETRY = ROOT / 'src2/v14_new_voidfraction/scripts/v14_geometry.py'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(name, rows):
    with (BASE / 'results' / name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def mii_fields(path):
    return {line.split(':', 1)[0].strip(): line.split(':', 1)[1].strip()
            for line in path.read_text().splitlines() if ':' in line}


def main():
    with PRIOR.open() as stream:
        inventory = list(csv.DictReader(stream))
    chosen = [r for r in inventory if r['included'] == 'True']
    assert len(chosen) == 26
    assert {r['run'] for r in inventory if r['included'] == 'False'} == {
        'm03q01_001', 'm03q80_001', 'm04q05_001', 'm05q05_002', 'm07q10_001'}
    geometry = runpy.run_path(str(GEOMETRY))
    lumen = tuple(geometry['LUMEN_ROWS'])
    assert lumen == (18, 71) and geometry['DIAMETER_PX'] == 53
    rows, images, provenance, workbook_rows = [], [], [], []
    for selected in chosen:
        folder = DATA / selected['run']
        parser = configparser.ConfigParser(interpolation=None, strict=False)
        parser.optionxform = str
        parser.read(folder / '_config.xsv')
        acquisition = parser['Acquisition']
        info = mii_fields(folder / '_config.mii')
        for field in ('Lens type 1', 'Lens focal length 1', 'Focus 1', 'Distortion index 1',
                      'Correction parameter file 1', 'Reference system 1'):
            assert info[field] == 'NOVALUE', (selected['run'], field, info[field])
        assert info['Movie images corrected 1'] == 'NO'
        assert float(info['Aspect ratio of pixels 1']) == 1
        fps = float(acquisition['Rate'])
        assert fps == float(info['Film speed 1']) == float(selected['native_fps'])
        assert sha(folder/'_config.xsv') == selected['config_sha256']
        assert sha(folder/'_config.mii') == selected['mii_sha256']
        files = sorted(folder.glob('ImgA*.jpg'))
        ids = [int(re.search(r'([0-9]+)$', p.stem).group(1)) for p in files]
        expected = list(range(int(acquisition['StartIdx']), int(acquisition['StopIdx']) + 1,
                              int(acquisition['Step'])))
        assert ids == expected and len(files) == int(acquisition['Frames']) == 501
        assert len(files) == int(info['Number of images 1'])
        exposure_ns = int(acquisition['ExposureNs'])
        assert math.isclose(exposure_ns / 1e9, float(info['Shutter time 1']), abs_tol=5e-7)
        image_digests = hashlib.sha256()
        shapes, modes, exif_count, jfif_values = set(), set(), 0, set()
        for frame_id, path in zip(ids, files):
            digest = sha(path)
            image_digests.update((path.name + ':' + digest + chr(10)).encode())
            with Image.open(path) as image:
                width, height = image.size
                exif = image.getexif()
                named_exif = {str(ExifTags.TAGS.get(k, k)): str(v) for k, v in exif.items()}
                density = image.info.get('jfif_density')
                unit = image.info.get('jfif_unit')
                exif_count += bool(exif)
                shapes.add((width, height))
                modes.add(image.mode)
                jfif_values.add((unit, tuple(density) if density is not None else None))
                images.append(dict(run=selected['run'], frame_index=frame_id, file=path.name,
                                   sha256=digest, width_px=width, height_px=height, mode=image.mode,
                                   jfif_unit=unit, jfif_density=json.dumps(density),
                                   exif_fields=json.dumps(named_exif, sort_keys=True)))
        assert shapes == {(int(acquisition['ROIW']), int(acquisition['ROIH']))} == {(1016, 80)}
        assert exif_count == 0 and jfif_values == {(0, (1, 1))}
        channels = [{k: parser[s].get(k) for k in ('Enable', 'Units', 'ScaleFactor')}
                    for s in parser.sections() if s.startswith('DaChannel_')]
        assert channels and all(c == dict(Enable='0', Units='[V]', ScaleFactor='1.000000000000000')
                                for c in channels)
        workbook_metadata = []
        for path in sorted(folder.glob('*.xlsx')):
            book = openpyxl.load_workbook(path, read_only=True, data_only=True)
            metadata = {}
            if 'metadata' in book.sheetnames:
                for cells in book['metadata'].iter_rows(values_only=True):
                    if len(cells) >= 2 and cells[0] is not None:
                        metadata[str(cells[0])] = cells[1]
            entry = dict(run=selected['run'], file=path.name, sha256=sha(path),
                         sheets=json.dumps(book.sheetnames), metadata=json.dumps(metadata, sort_keys=True))
            workbook_rows.append(entry)
            workbook_metadata.append(entry)
            book.close()
        row = dict(run=selected['run'], pattern=selected['pattern'], native_fps=fps,
                   frames=len(files), interval_count=len(files)-1, duration_s=(ids[-1]-ids[0])/fps,
                   frame_step=int(acquisition['Step']), exposure_ns=exposure_ns,
                   exposure_fraction_of_frame=exposure_ns/1e9*fps,
                   output_rate=float(acquisition['Output_Rate']), width_px=1016, height_px=80,
                   roi_x=int(acquisition['ROIX']), roi_y=int(acquisition['ROIY']),
                   camera_type=info['Camera type 1'], camera_id=info['Camera ID-number 1'],
                   pixel_aspect_ratio=float(info['Aspect ratio of pixels 1']),
                   binning=int(acquisition['Binning']), lens_type=info['Lens type 1'],
                   lens_focal_length=info['Lens focal length 1'], focus=info['Focus 1'],
                   distortion_index=info['Distortion index 1'],
                   movie_images_corrected=info['Movie images corrected 1'],
                   correction_parameter_file=info['Correction parameter file 1'],
                   reference_system=info['Reference system 1'], exif_nonempty_frames=exif_count,
                   jfif_metadata=json.dumps(sorted(jfif_values)),
                   documented_camera_m_per_px=False,
                   config_sha256=sha(folder/'_config.xsv'), mii_sha256=sha(folder/'_config.mii'),
                   mme_sha256=sha(folder/'_config.mme'), ordered_image_hashes_sha256=image_digests.hexdigest())
        rows.append(row)
        provenance.append(dict(run=selected['run'], metadata_files=[str(p.relative_to(ROOT))
                               for p in sorted(folder.glob('_config.*'))], mii_fields=info,
                               mme_fields=mii_fields(folder/'_config.mme'),
                               acquisition_fields=dict(acquisition), voltage_channel_fields=channels,
                               workbooks=workbook_metadata))
    write_csv('calibration_inventory.csv', rows)
    write_csv('calibration_image_metadata.csv', images)
    write_csv('calibration_workbook_metadata.csv', workbook_rows)
    literature_dir = BASE / 'intermediate/calibration_literature'
    literature_dir.mkdir(exist_ok=True)
    papers = []
    for filename, stem, assertions in [
        ('deoliveira2021-flow-patterns-pressure-drop.pdf', 'deoliveira2021',
         [(3, '1.0 mm internal diameter'), (4, '265 mm'), (4, '5000'), (4, 'subcooled condition')]),
        ('dias2025-knn-flow-pattern.pdf', 'dias2025_knn',
         [(4, '13.68'), (4, '3000'), (5, '1016'), (5, '88 pixels'), (5, 'exactly 200 frames')])]:
        pdf = ROOT / 'assets' / filename
        destination = literature_dir / (stem + '.txt')
        subprocess.run(['pdftotext', '-layout', str(pdf), str(destination)], check=True)
        pages = destination.read_text().split(chr(12))
        evidence = []
        for page_number, phrase in assertions:
            assert phrase in pages[page_number-1], (filename, page_number, phrase)
            evidence.append(dict(pdf_page=page_number, phrase=phrase))
        papers.append(dict(file=str(pdf.relative_to(ROOT)), sha256=sha(pdf),
                           extract=str(destination.relative_to(BASE)), extract_sha256=sha(destination),
                           directly_checked_text=evidence))
        render_pages = (3, 4) if stem == 'deoliveira2021' else (4, 5)
        render_stem = 'deoliveira2021' if stem == 'deoliveira2021' else 'dias2025'
        for number in render_pages:
            subprocess.run(['pdftoppm', '-f', str(number), '-l', str(number), '-scale-to', '1600',
                            '-png', '-singlefile', str(pdf), str(literature_dir/f'{render_stem}_p{number}')], check=True)
    scale = .001 / geometry['DIAMETER_PX']
    report = dict(status='passed', selected_recordings=len(rows), native_frames_checked=len(images),
                  fps_counts=dict(sorted(Counter(str(int(r['native_fps'])) for r in rows).items())),
                  image_dimensions=sorted({f"{r['width_px']}x{r['height_px']}" for r in rows}),
                  camera_ids=sorted({r['camera_id'] for r in rows}),
                  pixel_aspect_ratios=sorted({r['pixel_aspect_ratio'] for r in rows}),
                  exif_nonempty_frames=sum(r['exif_nonempty_frames'] for r in rows),
                  workbook_files_checked=len(workbook_rows),
                  spatial_calibration_found=False,
                  spatial_calibration_conclusion='No object-space axial m/pixel calibration in the inspected camera/image/workbook metadata.',
                  recommended_main_units='px/s for accepted image-structure motion; dimensionless within-recording ratios',
                  conditional_geometry=dict(tube_diameter_m=.001, tube_diameter_source='de Oliveira 2021, PDF pp.3-4; campaign association remains assumed',
                      lumen_rows=list(lumen), lumen_pixel_cells=geometry['DIAMETER_PX'],
                      conditional_m_per_px=scale, conditional_um_per_px=scale*1e6,
                      assumptions=['1 mm source tube diameter applies to these recordings',
                                   'owner-confirmed 53 radial pixel cells map to the physical full diameter',
                                   'radial apparent scale also applies axially at tracked feature locations',
                                   'optical magnification is stable across compared recordings'],
                      status='Geometric sensitivity conversion, not measured axial calibration; no accuracy bound established.'),
                  prior_inventory=dict(file=str(PRIOR.relative_to(ROOT)),sha256=sha(PRIOR)),
                  audit_script=dict(file=str(Path(__file__).resolve().relative_to(ROOT)),sha256=sha(Path(__file__))),
                  current_geometry=dict(file=str(GEOMETRY.relative_to(ROOT)),sha256=sha(GEOMETRY)),
                  source_papers=papers, configurations=provenance,
                  interpretation='Timing and metadata provenance audit; it does not validate phase identity, superficial velocities or physical velocity accuracy.')
    (BASE/'results/calibration_audit.json').write_text(json.dumps(report, indent=2) + chr(10))
    print(json.dumps({k:v for k,v in report.items() if k not in ('configurations','source_papers')}, indent=2))


if __name__ == '__main__':
    main()
