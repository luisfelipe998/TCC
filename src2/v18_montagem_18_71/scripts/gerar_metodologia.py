#!/usr/bin/env python3
"""Generate the four manuscript methodology illustrations inside this package.

The flowchart, wavelet profile example and optical-measure diagram are
deterministic didactic examples. The mask example uses the actual optical
panorama and inferred mask of m06q01_001 produced by stages 01 and 03.
Run this stage after reconstruction and optical segmentation.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

from runtime import OUTPUT_ROOT, SCRIPTS_ROOT, reference


ILLUSTRATIONS = (
    ("fluxograma", "diagram", False),
    ("perfis_wavelet", "synthetic_wavelet_example", True),
    ("exemplo_mascara", "experimental_crop", True),
    ("medidas_opticas", "synthetic_mask_example", True),
)


def main() -> None:
    destination = OUTPUT_ROOT / "05_figuras/metodologia"
    destination.mkdir(parents=True, exist_ok=True)
    figures = []
    for name, kind, has_arrays in ILLUSTRATIONS:
        script = SCRIPTS_ROOT / f"metodologia_{name}.py"
        # Each original illustration sets its own drawing style. A separate
        # process prevents one style from leaking into the next illustration.
        subprocess.run([sys.executable, "-B", str(script)], check=True)
        paths = [destination / f"{name}.{suffix}" for suffix in ("pdf", "png")]
        paths.append(destination / f"{name}_origem.json")
        if has_arrays:
            paths.append(destination / f"{name}_dados.npz")
        for path in paths:
            if not path.is_file():
                raise RuntimeError(f"Missing methodology artifact: {path}")
        figures.append({
            "name": name,
            "kind": kind,
            "experimental_data": kind == "experimental_crop",
            "generator": reference(script),
            "generator_sha256": hashlib.sha256(script.read_bytes()).hexdigest(),
            "artifacts": {
                reference(path): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in paths
            },
        })
    (destination / "manifesto.json").write_text(
        json.dumps({"status": "complete", "figures": figures},
                   indent=2, ensure_ascii=False) + "\n"
    )
    print(json.dumps({"methodology_figures": len(figures),
                      "destination": reference(destination)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
