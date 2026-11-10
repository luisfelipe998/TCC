"""Configuração e procedência da análise espacial integral v16."""
import csv
import hashlib
import json
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
os.environ.setdefault("MPLCONFIGDIR", "/tmp/pearson-v16-comparacao-mpl")
import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT / "_OLD"
INTER = ROOT / "intermediario"
DATA = ROOT.parent / "data_processed"
STATS = ROOT / "04_estatistica"
TABLES = STATS / "dados"
PATTERNS = ["churn", "wavy annular", "plug", "slug", "dispersed bubble"]
WAVES = np.geomspace(50, 3000, 200)
ROW_START, ROW_STOP, ROW_END, ROW_COUNT = 18, 72, 71, 54


def sha(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def read_csv(path):
    with Path(path).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def pattern_name(value):
    return "dispersed bubble" if value == "dispersed bubbles" else value


def slug(value):
    return value.replace(" ", "_")


def load_image(run):
    with Image.open(DATA/run/"__panorama_v11.png") as image:
        return np.asarray(image.convert("L"))


def load_raw(run):
    with np.load(INTER/"entradas"/f"{run}.npz") as archive:
        return archive["raw"].copy()


def records():
    return json.loads((INTER/"aquisicoes.json").read_text())


def pairs():
    return json.loads((INTER/"pares.json").read_text())


def number(value, digits=4):
    return f"{value:.{digits}f}".replace(".", ",")
