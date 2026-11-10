"""Renderiza os resultados ópticos após a remontagem usando 54 linhas.

Execução: python3 -B TCC/figuras/v18_optica/gerar_figuras_opticas.py

Lê as máscaras e tabelas da v18, sem executar a segmentação. O painel exibe
somente as linhas 18 a 71. As funções do gráfico de dispersão são extraídas do
gerador histórico, sem executar seu programa ou modificar os dados antigos.
"""

import ast
from collections import Counter
import csv
import hashlib
from io import BytesIO
import json
import os
from pathlib import Path

os.environ["MPLCONFIGDIR"] = "/tmp/tcc-v18-optica-figuras-mpl"

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
import numpy as np
from PIL import Image


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
SOURCE = REPO / "src2/v18_montagem_18_71/03_optica"
PANORAMA_MANIFEST = REPO / "src2/v18_montagem_18_71/01_panoramas/manifesto.json"
SCATTER_SCRIPT = REPO / "src2/v14.1_new_voidfraction/intermediate/history/before_roi_unification_20261030/scripts/summary_figures.py"
RUNS = ("m07q01_001", "m06q01_001", "m06q05_001", "m06q20_001", "m06q80_001")
PATTERNS = ("dispersed bubbles", "plug", "slug", "churn", "wavy annular")
CSV_NAMES = (
    "results/void_fraction/panorama_void_per_run.csv",
    "results/flip_rate/phase_mask_flip_per_run.csv",
    "results/void_fraction/void_fraction_flip_plot_data.csv",
)
PANEL_NAME = "panorama_phase_masks.png"
SCATTER_NAME = "void_fraction_vs_flip_rate"
TITLE = "Panorama (linhas 18–71; primeiros ≤2.500 px)"
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


def panel_inputs():
    images, inputs = [], []
    for run in RUNS:
        records = json.loads(PANORAMA_MANIFEST.read_text())["records"]
        item = next(record for record in records if record["ensaio"] == run)
        panorama = REPO / item["optica_path"]
        mask_path = SOURCE / "intermediate/phase_masks" / f"{run}_mask.png"
        with Image.open(panorama) as image:
            intensity = np.asarray(image.convert("L"), dtype=float)
        mask = read_mask(mask_path)
        assert intensity.shape == mask.shape
        images.append((run, intensity, mask))
        inputs.append({"run": run, "panorama": record(panorama), "mask": record(mask_path)})
    return images, inputs


def render_panel(images):
    """Desenha a ROI comum, mantendo os índices originais das linhas."""
    fig, axes = plt.subplots(len(images), 2, figsize=(15, 2 * len(images)))
    for pair, (run, intensity, mask) in zip(axes, images):
        stop = min(intensity.shape[1], 2500)
        crop = intensity[18:72, :stop]
        gas = mask[18:72, :stop]
        assert crop.shape == gas.shape == (54, stop)
        extent = (-0.5, stop - 0.5, 71.5, 17.5)
        for axis in pair:
            axis.imshow(crop, cmap="gray", vmin=0, vmax=255, aspect="auto",
                        interpolation="nearest", extent=extent)
            axis.set_yticks([18, 45, 71])
            axis.set_xticks([])
        color = np.zeros((54, stop, 4))
        color[..., :3] = (.1, .8, 1)
        color[..., 3] = gas * .42
        pair[1].imshow(color, aspect="auto", interpolation="nearest", extent=extent)
        pair[0].set_ylabel(run, fontsize=8)
    axes[0, 0].set_title(TITLE)
    axes[0, 1].set_title("Gás inferido em azul (mesma região)")
    fig.suptitle("Segmentação dos panoramas com referência líquida inferida")
    fig.tight_layout()
    stream = BytesIO()
    fig.savefig(stream, dpi=DPI)
    plt.close(fig)
    return stream.getvalue()


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
    panoramas = {r["ensaio"]: r for r in json.loads(PANORAMA_MANIFEST.read_text())["records"]}
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
        with Image.open(REPO / panoramas[run]["optica_path"]) as image:
            intensity = np.asarray(image.convert("L"))
        assert intensity.shape == mask.shape
        saturation = float(np.count_nonzero(intensity[18:72] == 255)) / (54 * width)
        assert abs(float(v["core_saturated_fraction"]) - saturation) < 1e-15
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
                       "alpha_from_mask": alpha, "saturation_rows_18_72": saturation,
                       "flip_count_from_mask": count,
                       "neighbor_comparisons": comparisons, "flip_rate_from_mask": rate})
    metadata = json.loads((SOURCE / "results/validation/metadata.json").read_text())
    assert len(rows) == len(panoramas)
    assert sum(r["plotted"] for r in rows) == metadata["eligible_acquisitions"]
    return checks


def main():
    scatter = load_rendering(
        SCATTER_SCRIPT,
        {"hull_indices", "save_figure", "portuguese_ticks", "alternation_chart"},
        {"plt": plt, "np": np, "Line2D": Line2D, "FuncFormatter": FuncFormatter,
         "PATTERNS": PATTERNS},
        {"COLORS", "PATTERN_LABELS_PT", "MARKERS"},
    )
    historical = [SCATTER_SCRIPT]
    old_hashes = {str(path): sha256(path) for path in historical}
    rows = plot_rows(SOURCE)
    checks = check_new_data(rows)
    images, inputs = panel_inputs()
    new_panel = render_panel(images)
    assert pixels(new_panel).shape == (1800, 2700, 4)
    (HERE / PANEL_NAME).write_bytes(new_panel)
    scatter["artifact_path"] = lambda name: HERE / name
    original_save = scatter["save_figure"]

    def save_current_figure(fig, basename):
        # The denominator counts neighboring pairs, as recorded in the CSV.
        axis = fig.axes[0]
        assert "transições por 1.000 pixels" in axis.get_ylabel()
        axis.set_ylabel("Flip-rate espacial\n(transições por 1.000 pares de pixels)",
                        fontsize=axis.yaxis.label.get_fontsize())
        original_save(fig, basename)

    scatter["save_figure"] = save_current_figure
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
    with (SOURCE / "results/summary/optical_group_summary.csv").open(newline="") as stream:
        expected_centroids = list(csv.DictReader(stream))
    for group in expected_centroids:
        actual = centroids[group["pattern"]]
        assert actual["n"] == int(group["n"])
        assert abs(actual["alpha_mean"] - float(group["alpha_arithmetic_mean"])) < 1e-14
        assert abs(actual["flip_rate_per1000px_mean"] - float(group["F1000_arithmetic_mean"])) < 1e-12
    assert all(sha256(Path(path)) == value for path, value in old_hashes.items())
    manifest = {
        "source_run": str(SOURCE.relative_to(REPO)),
        "panorama_manifest": record(PANORAMA_MANIFEST),
        "generator": record(Path(__file__).resolve()),
        "source_rendering_scripts": [record(Path(__file__).resolve()), record(SCATTER_SCRIPT)],
        "geometry": {"first_row": 18, "last_row_inclusive": 71,
                     "python_slice": "18:72", "row_count": 54, "diameter_px": 54},
        "rendering": {"matplotlib": matplotlib.__version__, "numpy": np.__version__,
                      "dpi": DPI, "segmentation_performed_by_this_generator": False,
                      "panel_rows": "18:72",
                      "original_colors_and_panel_runs_preserved": True,
                      "vertical_axis_unit": "transitions per 1000 neighboring pixel pairs",
                      "axis_label_precision_correction": "pixels replaced by neighboring pixel pairs"},
        "verification": {"historical_rendering_script_unchanged": True,
                         "new_input_hashes_match_source_manifest": True,
                         "new_csv_counts_and_alpha_match_saved_masks": True,
                         "saturation_recomputed_on_rows_18_72": True,
                         "hulls_match_source_validation": True,
                         "centroids_match_source_summary": True,
                         "acquisitions": len(rows),
                         "displayed_acquisitions": sum(r["plotted"] for r in rows),
                         "displayed_pattern_counts": dict(Counter(r["supplied_pattern"]
                                                                  for r in rows if r["plotted"]))},
        "input_tables": [record(SOURCE / name) for name in CSV_NAMES],
        "mask_checks": checks,
        "figures": [
            {**record(HERE / PANEL_NAME), "inputs": inputs,
             "size_pixels": [2700, 1800], "first_columns": 2500,
             "rows_start_stop_exclusive": [18, 72], "mask_opacity": 0.42},
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
