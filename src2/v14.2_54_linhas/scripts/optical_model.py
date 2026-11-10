"""Load unchanged v14.1 optical functions into an isolated 54-row namespace.

No legacy module is imported or modified. Only the five listed function
definitions are compiled, after the new geometry is supplied, so their default
arguments also use the new geometry. The historical source is hash-locked.
"""
from __future__ import annotations

import ast
import hashlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from PIL import Image
from scipy.ndimage import binary_closing, binary_fill_holes, gaussian_filter1d, label

SOURCE_SHA256 = "68e7d0b7ea4f37a6a4ea474f68e440510a2887ef960a258dd52aa92faa6658f2"
FUNCTION_NAMES = ("normalize", "make_reference", "close_shadow_crop", "silhouette", "geometry")
LUMEN_ROWS = (18, 72)
LUMEN_HEIGHT_PX = 54
DIAMETER_PX = 54.0
REFERENCE_RUNS = ("m05q01_001", "m06q01_001", "m07q05_001")
THRESHOLD = 0.70


def load_model(source_root: Path):
    source = source_root / "scripts/panorama_void_alternatives.py"
    source_bytes = source.read_bytes()
    if hashlib.sha256(source_bytes).hexdigest() != SOURCE_SHA256:
        raise RuntimeError("Historical segmentation source changed. Review its provenance before running.")
    tree = ast.parse(source_bytes, filename=str(source))
    definitions = [node for node in tree.body
                   if isinstance(node, ast.FunctionDef) and node.name in FUNCTION_NAMES]
    assert tuple(node.name for node in definitions) == FUNCTION_NAMES

    def read(run):
        path = source_root / "intermediate/reconstructed_inputs" / run / "__panorama_v11.png"
        with Image.open(path) as image:
            return np.asarray(image.convert("L"), dtype=float)

    scope = dict(np=np, gaussian_filter1d=gaussian_filter1d,
                 binary_closing=binary_closing, binary_fill_holes=binary_fill_holes,
                 label=label, read=read, LUMEN_ROWS=LUMEN_ROWS,
                 DIAMETER_PX=DIAMETER_PX, REFERENCE_RUNS=REFERENCE_RUNS, THRESHOLD=THRESHOLD)
    module = ast.Module(body=definitions, type_ignores=[])
    exec(compile(module, str(source), "exec"), scope)
    assert scope["silhouette"].__defaults__ == (THRESHOLD, LUMEN_ROWS, 20)
    assert scope["geometry"].__defaults__ == (DIAMETER_PX,)
    return SimpleNamespace(read=read, **{name: scope[name] for name in FUNCTION_NAMES})
