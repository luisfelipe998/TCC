"""Shared input, output, and report-link paths for the organized v14 analysis."""
from pathlib import Path
import json
import os
import re

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent
DATA = ROOT.parent / "data_processed"
BASELINE = ROOT / "intermediate" / "history" / "v14_baseline"
PANORAMAS = ROOT / "intermediate" / "reconstructed_inputs"
CACHE = ROOT / "intermediate" / "cache"
LAYOUT = json.loads((SCRIPT_DIR / "artifact_layout.json").read_text())


def panorama_path(run):
    """v14.1 spatial input; native frames and original PNGs remain read-only."""
    path = PANORAMAS / Path(run).name / "__panorama_v11.png"
    if not path.is_file():
        raise FileNotFoundError(f"Reconstruct inputs before analysis: {path}")
    return path


def original_panorama_path(run):
    return DATA / Path(run).name / "__panorama_v11.png"


def artifact_path(name):
    """Resolve a known artifact; new artifacts must be registered in the layout."""
    try:
        path = ROOT / LAYOUT[str(name)]
    except KeyError as exc:
        raise KeyError(f"Register {name!r} in scripts/artifact_layout.json") from exc
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def relocate_markdown(text, document):
    """Keep original flat-layout links usable in relocated and regenerated reports."""
    document = Path(document)

    def replace_link(match):
        target = match.group(2)
        if re.match(r"[a-zA-Z][a-zA-Z0-9+.-]*:", target) or target.startswith("#"):
            return match.group(0)
        href, separator, fragment = target.partition("#")
        if href in LAYOUT:
            resolved = ROOT / LAYOUT[href]
        elif href.split("/", 1)[0] in {"findings", "results", "figures", "scripts", "intermediate", "literature"}:
            resolved = ROOT / href
        elif (document.parent / href).exists():
            resolved = document.parent / href
        elif (ROOT / href).exists():
            resolved = ROOT / href
        else:
            return match.group(0)
        relative = Path(os.path.relpath(resolved, document.parent)).as_posix()
        return f"[{match.group(1)}]({relative}{separator}{fragment})"

    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", replace_link, text)
    return re.sub(r"python3 src2/v14\.1_new_voidfraction/(?!scripts/)",
                  "python3 src2/v14.1_new_voidfraction/scripts/", text)


def write_markdown(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path.write_text(relocate_markdown(text, path))
