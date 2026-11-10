#!/usr/bin/env python3
"""Compara um recorte real com sua máscara óptica já calculada na versão 18.

Não executa segmentação, não modifica as entradas e não retoca pixels. O recorte
fixo das primeiras 520 colunas de m06q01_001 ilustra uma estrutura alongada e a
região adjacente. A figura serve para explicar o método, sem validar a máscara.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path
import platform
import zipfile
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/tcc-v18-metodologia-mpl")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import numpy as np
from PIL import Image, __version__ as pillow_version


from runtime import OUTPUT_ROOT, reference

OUT = OUTPUT_ROOT / "05_figuras/metodologia"
SOURCE = OUTPUT_ROOT
RUN = "m06q01_001"
ROWS = (18, 72)
COLUMNS = (0, 520)
WIDTH_CM = 16.0
HEIGHT_CM = 6.4
GAS = "#247b9b"
OTHER = "#e4e7e9"
INK = "#25343c"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stable_npz(path: Path, arrays: dict[str, np.ndarray]) -> None:
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_STORED) as archive:
        for name in sorted(arrays):
            buffer = io.BytesIO()
            np.lib.format.write_array(buffer, arrays[name], allow_pickle=False)
            entry = zipfile.ZipInfo(name + ".npy", date_time=(1980, 1, 1, 0, 0, 0))
            entry.external_attr = 0o600 << 16
            archive.writestr(entry, buffer.getvalue())


def draw(panorama: np.ndarray, mask: np.ndarray) -> tuple[Path, Path]:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9.0,
        "text.color": INK,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    fig = plt.figure(figsize=(WIDTH_CM / 2.54, HEIGHT_CM / 2.54), facecolor="white")
    image_width = 15.0
    image_height = image_width * panorama.shape[0] / panorama.shape[1]

    def image_axis(bottom_cm: float):
        ax = fig.add_axes((0.6 / WIDTH_CM, bottom_cm / HEIGHT_CM,
                           image_width / WIDTH_CM, image_height / HEIGHT_CM))
        ax.set_xlim(0, panorama.shape[1])
        ax.set_ylim(panorama.shape[0], 0)
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_color("#adb5ba")
            spine.set_linewidth(0.6)
        return ax

    top = image_axis(4.1)
    top.imshow(panorama, cmap="gray", vmin=0, vmax=255, interpolation="nearest",
               extent=(0, panorama.shape[1], panorama.shape[0], 0), aspect="equal")
    top.set_xticks([])
    bottom = image_axis(1.65)
    bottom.imshow(mask, cmap=ListedColormap([OTHER, GAS]), vmin=0, vmax=1,
                  interpolation="nearest", aspect="equal",
                  extent=(0, mask.shape[1], mask.shape[0], 0))
    bottom.set_xticks(np.arange(0, 501, 100))
    bottom.tick_params(axis="x", labelsize=8, length=2.5, pad=2)
    bottom.set_xlabel("Posição longitudinal no recorte (pixels)", fontsize=9, labelpad=3)

    fig.text(0.6 / WIDTH_CM, 5.95 / HEIGHT_CM,
             "a) Imagem de entrada: trecho do panorama", fontsize=10, va="center")
    fig.text(0.6 / WIDTH_CM, 3.46 / HEIGHT_CM,
             "b) Máscara binária: a mesma região da imagem", fontsize=10, va="center")
    fig.legend(handles=[Patch(facecolor=GAS, label="1: gás inferido"),
                        Patch(facecolor=OTHER, edgecolor="#adb5ba", linewidth=0.6,
                              label="0: demais pixels")],
               loc="center", bbox_to_anchor=(0.5, 0.42 / HEIGHT_CM), ncol=2,
               frameon=False, handlelength=1.4, handleheight=0.9, columnspacing=2.4)

    pdf = OUT / "exemplo_mascara.pdf"
    png = OUT / "exemplo_mascara.png"
    fig.savefig(pdf, metadata={"Title": "Exemplo real de imagem e máscara óptica",
                              "Author": "Elaborado pelo autor",
                              "CreationDate": None, "ModDate": None})
    fig.savefig(png, dpi=300, metadata={"Software": "Matplotlib"})
    plt.close(fig)
    return pdf, png


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    panorama_path = SOURCE / "01_panoramas/optica" / RUN / "panorama.png"
    mask_path = SOURCE / "03_optica/intermediate/phase_masks" / f"{RUN}_mask.png"
    eligibility_path = SOURCE / "03_optica/results/void_fraction/void_fraction_flip_plot_data.csv"
    inputs = [panorama_path, mask_path, eligibility_path]
    before = {reference(path): sha256(path) for path in inputs}
    with eligibility_path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    eligible = next(row for row in rows if row["run"] == RUN)
    assert eligible["eligible_pure_unsaturated"] == "True"
    assert eligible["saturated"] == "False"

    with Image.open(panorama_path) as im:
        assert im.mode == "L"
        panorama = np.array(im)
    with Image.open(mask_path) as im:
        assert im.mode == "L"
        encoded = np.array(im)
    assert panorama.shape == encoded.shape == (80, 7786)
    assert set(np.unique(encoded)) == {0, 255}
    assert not encoded[:ROWS[0]].any() and not encoded[ROWS[1]:].any()
    rr, cc = slice(*ROWS), slice(*COLUMNS)
    cropped_panorama = panorama[rr, cc].copy()
    cropped_mask = (encoded[rr, cc] == 255).astype(np.uint8)
    assert cropped_panorama.shape == cropped_mask.shape == (54, 520)
    assert set(np.unique(cropped_mask)) == {0, 1}
    data_path = OUT / "exemplo_mascara_dados.npz"
    stable_npz(data_path, {
        "panorama_uint8": cropped_panorama,
        "mask_binary": cropped_mask,
        "source_rows_start_stop": np.array(ROWS, dtype=np.int64),
        "source_columns_start_stop": np.array(COLUMNS, dtype=np.int64),
    })
    pdf, png = draw(cropped_panorama, cropped_mask)
    after = {reference(path): sha256(path) for path in inputs}
    assert before == after
    with np.load(data_path, allow_pickle=False) as data:
        assert np.array_equal(data["panorama_uint8"], panorama[rr, cc])
        assert np.array_equal(data["mask_binary"], encoded[rr, cc] == 255)

    metadata = {
        "purpose": "Exemplo real e didático do panorama e de sua máscara óptica já calculada; não valida a segmentação.",
        "source_version": "v18_montagem_18_71",
        "run_id": RUN,
        "pattern": eligible["supplied_pattern"],
        "eligible_pure_unsaturated": True,
        "input_shape_height_width": list(panorama.shape),
        "crop_rows_start_stop_exclusive": list(ROWS),
        "crop_columns_start_stop_exclusive": list(COLUMNS),
        "crop_shape_height_width": list(cropped_panorama.shape),
        "crop_selection": "Primeiras 520 colunas do ensaio m06q01_001, incluindo uma estrutura alongada e região adjacente; recorte ilustrativo, sem otimizar concordância de segmentação.",
        "mask_decoding": "PNG monocromático: 255 = gás inferido (1), 0 = demais pixels (0).",
        "display": {"width_cm": WIDTH_CM, "height_cm": HEIGHT_CM,
                    "grayscale_range": [0, 255], "contrast_stretching": False,
                    "interpolation": "nearest", "aspect": "equal",
                    "gas_color": GAS, "other_color": OTHER},
        "segmentation_executed": False,
        "pixel_editing": False,
        "inputs_sha256_before": before,
        "inputs_sha256_after": after,
        "checks": {"inputs_unchanged": True, "same_crop_for_both_panels": True,
                   "mask_only_0_or_255": True, "mask_empty_outside_rows_18_72": True,
                   "saved_arrays_match_source_pixels": True,
                   "eligible_and_unsaturated": True},
        "outputs_sha256": {p.name: sha256(p) for p in (pdf, png, data_path)},
        "script_sha256": sha256(Path(__file__)),
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "matplotlib": matplotlib.__version__, "pillow": pillow_version},
    }
    (OUT / "exemplo_mascara_origem.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"PASS: {RUN}, rows {ROWS}, columns {COLUMNS}; input pixels preserved.")
    print(reference(pdf))
    print(reference(png))


if __name__ == "__main__":
    main()
