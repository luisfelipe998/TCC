"""Portable paths shared by the self-contained v18 processing stages.

Only DATA_ROOT may point to project data outside this package. OUTPUT_ROOT can
select a fresh candidate directory so published results are never overwritten
before the migration comparison has passed.
"""
from __future__ import annotations

import os
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_ROOT = Path(os.environ.get("TCC_V18_OUTPUT", PACKAGE_ROOT)).resolve()
DATA_ROOT = Path(os.environ.get("TCC_DATA_ROOT", PACKAGE_ROOT.parent / "data_processed")).resolve()
SCRIPTS_ROOT = PACKAGE_ROOT / "scripts"
CONFIG_ROOT = PACKAGE_ROOT / "config"


def artifact_path(value: str | Path) -> Path:
    """Read portable artifact names and the pre-consolidation v18 prefix."""
    value = str(value).replace("\\", "/")
    legacy = "src2/v18_montagem_18_71/"
    if value.startswith(legacy):
        value = value[len(legacy):]
    path = Path(value)
    if path.is_absolute():
        resolved = path.resolve()
    else:
        resolved = (OUTPUT_ROOT / path).resolve()
    if not resolved.is_relative_to(OUTPUT_ROOT):
        raise ValueError(f"Artifact must remain inside the output directory: {value}")
    return resolved


def reference(path: str | Path) -> str:
    """Portable provenance label; never resolve a label as an unchecked path."""
    path = Path(path).resolve()
    for base, prefix in ((OUTPUT_ROOT, ""), (PACKAGE_ROOT, "package/"), (DATA_ROOT, "data_processed/")):
        if path.is_relative_to(base):
            return prefix + path.relative_to(base).as_posix()
    raise ValueError(f"Undeclared external dependency: {path}")
