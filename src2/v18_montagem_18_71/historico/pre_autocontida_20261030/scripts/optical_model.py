"""Load the current single-ROI optical functions, with source hash verification.

Numerical operations are isolated from the historical supplemental scripts.
"""
from __future__ import annotations

import ast
import hashlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from PIL import Image
from scipy.ndimage import binary_closing, binary_fill_holes, gaussian_filter1d, label

SOURCE_SHA256 = "423ba8dc0fc8dd84ebe40a46b1851247689e0584b2fc95968e214ac3a53781d8"
FUNCTION_NAMES = ("crop_roi", "roi_saturation", "normalize", "make_reference", "close_shadow_crop", "silhouette", "geometry")
LUMEN_ROWS = (18, 72)
LUMEN_HEIGHT_PX = 54
DIAMETER_PX = 54.0
REFERENCE_RUNS = ("m05q01_001", "m06q01_001", "m07q05_001")
THRESHOLD = 0.70


def load_model(source_root: Path, panorama_paths: dict[str, Path]):
    source = source_root / "scripts/panorama_void_alternatives.py"
    source_bytes = source.read_bytes()
    if hashlib.sha256(source_bytes).hexdigest() != SOURCE_SHA256:
        raise RuntimeError("Current segmentation source changed. Review its provenance before running.")
    tree = ast.parse(source_bytes, filename=str(source))
    definitions = [node for node in tree.body
                   if isinstance(node, ast.FunctionDef) and node.name in FUNCTION_NAMES]
    assert tuple(node.name for node in definitions) == FUNCTION_NAMES

    def read(run):
        path = panorama_paths[run]
        with Image.open(path) as image:
            return np.asarray(image.convert("L"), dtype=float)

    scope = dict(np=np, gaussian_filter1d=gaussian_filter1d,
                 binary_closing=binary_closing, binary_fill_holes=binary_fill_holes,
                 label=label, read=read, LUMEN_ROWS=LUMEN_ROWS, LUMEN_HEIGHT_PX=LUMEN_HEIGHT_PX,
                 DIAMETER_PX=DIAMETER_PX, REFERENCE_RUNS=REFERENCE_RUNS, THRESHOLD=THRESHOLD)
    module = ast.Module(body=definitions, type_ignores=[])
    exec(compile(module, str(source), "exec"), scope)
    assert scope["silhouette"].__defaults__ == (THRESHOLD, LUMEN_ROWS, 20)
    assert scope["geometry"].__defaults__ == (DIAMETER_PX,)
    return SimpleNamespace(read=read, **{name: scope[name] for name in FUNCTION_NAMES})
