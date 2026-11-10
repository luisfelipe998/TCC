"""Renderiza os resultados ópticos de 54 linhas com a apresentação anterior.

Execução: python3 TCC/figuras/v14_2/gerar_figuras_opticas.py

Lê as máscaras e tabelas de v14.2, sem executar a segmentação. As funções de
desenho são extraídas dos geradores históricos, sem executar seus programas.
Uma reprodução dos dados antigos deve coincidir pixel a pixel com as duas
figuras de v14_1 antes da renderização dos novos dados.
"""

import ast
from collections import Counter
import csv
import hashlib
from io import BytesIO
import json
import math
import os
from pathlib import Path
import tempfile

os.environ["MPLCONFIGDIR"] = "/tmp/tcc-v14-2-figuras-mpl"

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
import numpy as np
from PIL import Image


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
SOURCE = REPO / "src2/v14.2_54_linhas"
PREVIOUS = REPO / "src2/v14.1_new_voidfraction"
OLD_FIGURES = HERE.parent / "v14_1"
PANEL_SCRIPT = OLD_FIGURES / "gerar_painel_mascaras.py"
SCATTER_SCRIPT = PREVIOUS / "scripts/summary_figures.py"
RUNS = ("m07q01_001", "m06q01_001", "m06q05_001", "m06q20_001", "m06q80_001")
PATTERNS = ("dispersed bubbles", "plug", "slug", "churn", "wavy annular")
CSV_NAMES = (
    "results/void_fraction/panorama_void_per_run.csv",
    "results/flip_rate/phase_mask_flip_per_run.csv",
    "results/void_fraction/void_fraction_flip_plot_data.csv",
)
PANEL_NAME = "panorama_phase_masks.png"
SCATTER_NAME = "void_fraction_vs_flip_rate"
TITLE = "Panorama corrigido (primeiros ≤2.500 px)"
DPI = 180


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record(path):
    return {"file": str(path.relative_to(REPO)), "sha256": sha256(path)}


def load_rendering(path, functions, namespace, constants=()):
    """Carrega apenas as funções de apresentação e constantes especificadas."""
    tree = ast.parse(path.read_text())
    selected = []
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in functions:
            selected.append(node)
        elif isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in constants:
                    namespace[target.id] = ast.literal_eval(node.value)
    assert {node.name for node in selected} == set(functions)
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(path), "exec"), namespace)
    return namespace


def table(path):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == len({row["run"] for row in rows}), path
    return rows


def plot_rows(source):
    rows = table(source / CSV_NAMES[2])
    for row in rows:
        for key in ("plotted", "eligible_pure_unsaturated", "saturated"):
            assert row[key] in ("True", "False")
            row[key] = row[key] == "True"
    return rows


def read_mask(path):
    with Image.open(path) as image:
        encoded = np.asarray(image.convert("L"))
    assert set(np.unique(encoded)) <= {0, 255}, path
    assert encoded.shape[0] == 80, path
    return encoded == 255


def panel_inputs(source):
    images, inputs = [], []
    for run in RUNS:
        panorama = PREVIOUS / "intermediate/reconstructed_inputs" / run / "__panorama_v11.png"
        mask_path = source / "intermediate/phase_masks" / f"{run}_mask.png"
        with Image.open(panorama) as image:
            intensity = np.asarray(image.convert("L"), dtype=float)
        mask = read_mask(mask_path)
        assert intensity.shape == mask.shape
        images.append((run, intensity, mask))
        inputs.append({"run": run, "panorama": record(panorama), "mask": record(mask_path)})
    return images, inputs


def pixels(data):
    with Image.open(BytesIO(data)) as image:
        return np.asarray(image.convert("RGBA")).copy()


def check_new_data(rows):
    """Confere tabelas, contagens e reconstrução geométrica nas máscaras salvas."""
    void = {r["run"]: r for r in table(SOURCE / CSV_NAMES[0])}
    phase = {r["run"]: r for r in table(SOURCE / CSV_NAMES[1])}
    assert set(void) == set(phase) == {r["run"] for r in rows}
    saved = json.loads((SOURCE / "results/validation/artifact_manifest.json").read_text())["files"]
    for name in CSV_NAMES:
        assert sha256(SOURCE / name) == saved[name], name
    checks = []
    for row in rows:
        run = row["run"]
        v, p = void[run], phase[run]
        relative = f"intermediate/phase_masks/{run}_mask.png"
        mask_path = SOURCE / relative
        assert sha256(mask_path) == saved[relative], run
        mask = read_mask(mask_path)
        assert not mask[:18].any() and not mask[72:].any(), run
        lumen = mask[18:72]
        width = mask.shape[1]
        assert width == int(v["panorama_width_px"]) == int(p["panorama_width"])
        assert row["supplied_pattern"] == p["pattern"] == v["label"].replace("_", " ")
        assert row["plotted"] == row["eligible_pure_unsaturated"] == (p["eligible"] == "True")
        assert row["saturated"] == (float(v["core_saturated_fraction"]) > 0.20)
        assert row["plotted"] == (not row["saturated"])
        for item in (v, row):
            assert int(item["lumen_start_row"]) == 18
            assert int(item["lumen_stop_row_exclusive"]) == 72
            assert float(item["lumen_diameter_px"]) == 54
        for item in (p, row):
            assert item["candidate_mask_rows"] == "18:72"
            assert int(item["candidate_mask_n_rows"]) == 54
        count = int(np.count_nonzero(lumen[:, 1:] != lumen[:, :-1]))
        comparisons = 54 * (width - 1)
        rate = count / comparisons
        for item in (p, row):
            assert int(item["candidate_pixel_flip_count"]) == count
            assert int(item["candidate_pixel_flip_comparisons"]) == comparisons
            assert abs(float(item["candidate_pixel_flip_rate"]) - rate) < 1e-15
        assert abs(float(row["candidate_pixel_flips_per1000px"]) - 1000 * rate) < 1e-12
        # Cada intervalo vertical conectado fornece um diâmetro em pixels.
        changes = np.diff(np.pad(lumen.T.astype(np.int8), ((0, 0), (1, 1))), axis=1)
        start_columns, starts = np.nonzero(changes == 1)
        stop_columns, stops = np.nonzero(changes == -1)
        assert np.array_equal(start_columns, stop_columns)
        alpha = float(np.sum((stops - starts) ** 2)) / (54 ** 2 * width)
        for item in (v, row):
            assert abs(float(item["conditional_circle_alpha_spatial_mean"]) - alpha) < 1e-12
        checks.append({"run": run, "mask": record(mask_path), "width": width,
                       "alpha_from_mask": alpha, "flip_count_from_mask": count,
                       "neighbor_comparisons": comparisons, "flip_rate_from_mask": rate})
    assert len(rows) == 26 and sum(r["plotted"] for r in rows) == 25
    return checks


def main():
    panel = load_rendering(PANEL_SCRIPT, {"render"},
                           {"plt": plt, "np": np, "math": math,
                            "BytesIO": BytesIO, "DPI": DPI})
    scatter = load_rendering(
        SCATTER_SCRIPT,
        {"hull_indices", "save_figure", "portuguese_ticks", "alternation_chart"},
        {"plt": plt, "np": np, "Line2D": Line2D, "FuncFormatter": FuncFormatter,
         "PATTERNS": PATTERNS},
        {"COLORS", "PATTERN_LABELS_PT", "MARKERS"},
    )
    historical = [OLD_FIGURES / PANEL_NAME, OLD_FIGURES / f"{SCATTER_NAME}.png",
                  PANEL_SCRIPT, SCATTER_SCRIPT]
    old_hashes = {str(path): sha256(path) for path in historical}
    old_images, _ = panel_inputs(PREVIOUS)
    baseline, _ = panel["render"](old_images, TITLE)
    assert np.array_equal(pixels(baseline), pixels((OLD_FIGURES / PANEL_NAME).read_bytes()))
    with tempfile.TemporaryDirectory(prefix="tcc-v14-2-baseline-") as temporary:
        scatter["artifact_path"] = lambda name: Path(temporary) / name
        with plt.rc_context({"font.family": "DejaVu Sans", "svg.fonttype": "none"}):
            scatter["alternation_chart"](plot_rows(PREVIOUS), "candidate_pixel_flips_per1000px",
                                         SCATTER_NAME, show_centroids=True)
        assert np.array_equal(pixels((Path(temporary) / f"{SCATTER_NAME}.png").read_bytes()),
                              pixels((OLD_FIGURES / f"{SCATTER_NAME}.png").read_bytes()))

    rows = plot_rows(SOURCE)
    checks = check_new_data(rows)
    images, inputs = panel_inputs(SOURCE)
    new_panel, _ = panel["render"](images, TITLE)
    assert pixels(new_panel).shape == (1800, 2700, 4)
    (HERE / PANEL_NAME).write_bytes(new_panel)
    scatter["artifact_path"] = lambda name: HERE / name
    with plt.rc_context({"font.family": "DejaVu Sans", "svg.fonttype": "none"}):
        hulls, centroids = scatter["alternation_chart"](
            rows, "candidate_pixel_flips_per1000px", SCATTER_NAME, show_centroids=True)
    expected = json.loads((SOURCE / "results/validation/hull_intersections.json").read_text())
    for pattern, actual in hulls.items():
        reference = expected["primary_chart_hulls"][pattern]
        assert actual["eligible_runs"] == reference["eligible_runs"]
        assert actual["shaded"] == reference["shaded"]
        observed = actual["polygon_vertex_runs"]
        vertices = reference["polygon_vertex_runs"]
        # O vértice inicial pode variar sem alterar a mesma sequência cíclica.
        assert observed == vertices or any(
            observed == vertices[i:] + vertices[:i] for i in range(len(vertices)))
    assert all(sha256(Path(path)) == value for path, value in old_hashes.items())
    manifest = {
        "source_run": str(SOURCE.relative_to(REPO)),
        "panorama_source_run": str(PREVIOUS.relative_to(REPO)),
        "generator": record(Path(__file__).resolve()),
        "source_rendering_scripts": [record(PANEL_SCRIPT), record(SCATTER_SCRIPT)],
        "geometry": {"first_row": 18, "last_row_inclusive": 71,
                     "python_slice": "18:72", "row_count": 54, "diameter_px": 54},
        "rendering": {"matplotlib": matplotlib.__version__, "numpy": np.__version__,
                      "dpi": DPI, "segmentation_performed_by_this_generator": False,
                      "original_labels_colors_layout_and_panel_runs_preserved": True},
        "verification": {"old_panel_reproduced_pixel_for_pixel": True,
                         "old_scatter_reproduced_pixel_for_pixel": True,
                         "historical_figure_and_script_hashes_unchanged": True,
                         "new_input_hashes_match_source_manifest": True,
                         "new_csv_counts_and_alpha_match_saved_masks": True,
                         "hulls_match_source_validation": True,
                         "acquisitions": len(rows),
                         "displayed_acquisitions": sum(r["plotted"] for r in rows),
                         "displayed_pattern_counts": dict(Counter(r["supplied_pattern"]
                                                                  for r in rows if r["plotted"]))},
        "input_tables": [record(SOURCE / name) for name in CSV_NAMES],
        "mask_checks": checks,
        "figures": [
            {**record(HERE / PANEL_NAME), "inputs": inputs,
             "size_pixels": [2700, 1800], "first_columns": 2500, "mask_opacity": 0.42},
            {**record(HERE / f"{SCATTER_NAME}.png"),
             "svg": record(HERE / f"{SCATTER_NAME}.svg"),
             "size_pixels": list(pixels((HERE / f"{SCATTER_NAME}.png").read_bytes()).shape[1::-1]),
             "hulls": hulls, "centroids": centroids,
             "centroid_weighting": "equal weight per eligible acquisition"},
        ],
    }
    (HERE / "origens.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(manifest["verification"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
