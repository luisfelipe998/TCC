"""Paths, source checks and small I/O helpers for the bounded pilot."""
from pathlib import Path
import csv
import hashlib
import json
import os

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
OLD = ROOT.parent / "v16.2_wco_eventos"
BASE = ROOT.parent / "v16_spatial2_18_71"
DATA = ROOT.parent / "data_processed"
ANCHORS = ROOT / "01_ancoras"
WAVELETS = ROOT / "02_wavelets"
CONCLUSIONS = ROOT / "03_conclusoes"
INTER = ROOT / "intermediario"
os.environ.setdefault("MPLCONFIGDIR", str(INTER / "matplotlib"))
os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")

import numpy as np
from PIL import Image

CFG = json.loads((SCRIPTS / "config.json").read_text())


def setup():
    for path in (ANCHORS, WAVELETS, CONCLUSIONS, INTER, INTER / "mapas",
                 INTER / "previas", INTER / "matplotlib"):
        path.mkdir(parents=True, exist_ok=True)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_csv(path):
    with Path(path).open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with Path(path).open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2,
                                    allow_nan=False) + "\n", encoding="utf-8")


def slug(pattern):
    return pattern.replace(" ", "_")


def image_path(run):
    return DATA / run / "__panorama_v11.png"


def load_image(run):
    with Image.open(image_path(run)) as image:
        pano = np.asarray(image.convert("L"))
    lo, hi = CFG["linhas"]
    if (not all(isinstance(v, int) and not isinstance(v, bool) for v in (lo,hi))
            or not 0 <= lo <= hi < pano.shape[0]):
        raise ValueError("Faixa de linhas invalida para o panorama; nao truncar silenciosamente.")
    return pano, pano[lo:hi + 1].mean(axis=0)


def context_ok(center, length, window):
    return window // 2 + CFG["margem"] <= center <= length - window // 2 - CFG["margem"]


def required_window(center, left, right, kind):
    rule = CFG["contexto"]
    margin = rule["margem_evento"] if kind == "evento_delimitado" else rule["margem_vale"]
    diameter = 2 * (max(center - left, right - center) + margin)
    return max(rule["comprimento_minimo"], int(np.ceil(diameter / rule["passo"])) * rule["passo"])


def inventory():
    return sorted((r for r in read_csv(BASE / "results/input_audit.csv")
                   if r["included"] == "True" and r["pattern"] in CFG["padroes"]
                   and not r["run"].endswith("_002")), key=lambda r: r["run"])


def verify_source(record):
    run = record["run"]
    assert sha(image_path(run)) == record["panorama_sha256"], run
    pano, raw = load_image(run)
    assert pano.shape == (int(record["height"]), int(record["n"])), run
    with np.load(BASE / "intermediate" / f"{run}_input.npz") as cache:
        np.testing.assert_array_equal(raw, cache["raw"])
    return pano, raw


def selected_manifest():
    frozen = json.loads((INTER / "congelamento.json").read_text())
    for relative, expected in frozen["arquivos_sha256"].items():
        assert sha(ROOT / relative) == expected, f"Selecao modificada: {relative}"
    for relative, expected in frozen["fontes_sha256"].items():
        assert sha(ROOT.parent / relative) == expected, f"Fonte modificada: {relative}"
    return json.loads((INTER / "selecao.json").read_text())
