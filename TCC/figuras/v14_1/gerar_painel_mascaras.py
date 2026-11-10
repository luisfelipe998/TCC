"""Regenera o painel óptico com as imagens e máscaras já calculadas.

Execução: python3 TCC/figuras/v14_1/gerar_painel_mascaras.py

Reproduz o desenho de panorama_void_alternatives.py:250–265. Somente o título
da coluna esquerda muda de "Panorama original" para "Panorama corrigido".
As máscaras salvas são lidas diretamente, sem executar a segmentação. A
renderização com o título antigo deve reproduzir todos os pixels do PNG de
origem, e a alteração final deve ficar dentro da área dos dois títulos.
"""

from io import BytesIO
from pathlib import Path
import hashlib
import json
import math
import os
import sys

sys.dont_write_bytecode = True
os.environ["MPLCONFIGDIR"] = "/tmp/tcc-painel-mascaras-mpl"

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
SOURCE = REPO / "src2/v14.1_new_voidfraction"
SOURCE_FIGURE = SOURCE / "figures/void_fraction/panorama_phase_masks.png"
SOURCE_SCRIPT = SOURCE / "scripts/panorama_void_alternatives.py"
OUTPUT = HERE / "panorama_phase_masks.png"
MANIFEST = HERE / "origens.json"
RUNS = ("m07q01_001", "m06q01_001", "m06q05_001", "m06q20_001", "m06q80_001")
OLD_TITLE = "Panorama original (primeiros ≤2.500 px)"
NEW_TITLE = "Panorama corrigido (primeiros ≤2.500 px)"
DPI = 180


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record(path):
    return {"file": str(path.relative_to(REPO)), "sha256": sha256(path)}


def read_inputs():
    saved = json.loads((SOURCE / "results/validation/artifact_manifest.json").read_text())["files"]
    images, records = [], []
    for run in RUNS:
        panorama = SOURCE / "intermediate/reconstructed_inputs" / run / "__panorama_v11.png"
        mask_path = SOURCE / "intermediate/phase_masks" / f"{run}_mask.png"
        for path in (panorama, mask_path):
            assert sha256(path) == saved[str(path.relative_to(SOURCE))], path
        with Image.open(panorama) as image:
            intensity = np.asarray(image.convert("L"), dtype=float)
        with Image.open(mask_path) as image:
            encoded = np.asarray(image.convert("L"))
        assert intensity.shape == encoded.shape and intensity.shape[0] == 80, run
        assert set(np.unique(encoded)) <= {0, 255}, run
        images.append((run, intensity, encoded == 255))
        records.append({"run": run, "panorama": record(panorama), "mask": record(mask_path)})
    return images, records


def render(images, title):
    # Mesmas operações de apresentação, dimensões, recortes e cores do original.
    fig, axes = plt.subplots(len(images), 2, figsize=(15, 2 * len(images)))
    for pair, (run, intensity, mask) in zip(axes, images):
        stop = min(intensity.shape[1], 2500)
        pair[0].imshow(intensity[:, :stop], cmap="gray", vmin=0, vmax=255, aspect="auto")
        pair[1].imshow(intensity[:, :stop], cmap="gray", vmin=0, vmax=255, aspect="auto")
        color = np.zeros((80, stop, 4))
        color[..., 0] = .1
        color[..., 1] = .8
        color[..., 2] = 1
        color[..., 3] = mask[:, :stop] * .42
        pair[1].imshow(color, aspect="auto")
        pair[0].set_ylabel(run, fontsize=8)
        for ax in pair:
            ax.set_yticks([0, 40, 79])
            ax.set_xticks([])
    title_artist = axes[0, 0].set_title(title)
    axes[0, 1].set_title("Gás estimado em azul")
    fig.suptitle("Segmentação dos panoramas com referência líquida inferida")
    fig.tight_layout()

    # Limites do título para conferir a extensão das diferenças, sem editar pixels.
    bounds = title_artist.get_window_extent(fig.canvas.get_renderer())
    scale = DPI / fig.dpi
    height = int(round(fig.get_figheight() * DPI))
    padding = 4
    title_box = [math.floor(bounds.x0 * scale) - padding,
                 math.floor(height - bounds.y1 * scale) - padding,
                 math.ceil(bounds.x1 * scale) + padding,
                 math.ceil(height - bounds.y0 * scale) + padding]
    stream = BytesIO()
    fig.savefig(stream, dpi=DPI)
    plt.close(fig)
    return stream.getvalue(), title_box


def pixels(png):
    with Image.open(BytesIO(png)) as image:
        return np.asarray(image.convert("RGBA")).copy()


def main():
    images, inputs = read_inputs()
    original_bytes, old_box = render(images, OLD_TITLE)
    reference = pixels(SOURCE_FIGURE.read_bytes())
    original = pixels(original_bytes)
    assert np.array_equal(original, reference), "O desenho original não foi reproduzido pixel a pixel."

    corrected_bytes, new_box = render(images, NEW_TITLE)
    corrected = pixels(corrected_bytes)
    assert corrected.shape == reference.shape == (1800, 2700, 4)
    existing = pixels(OUTPUT.read_bytes())
    assert np.array_equal(existing, reference) or np.array_equal(existing, corrected), (
        "A figura do TCC contém alterações diferentes desta correção de título."
    )
    changed = np.any(corrected != reference, axis=2)
    ys, xs = np.nonzero(changed)
    assert len(xs) > 0, "Nenhuma alteração de título foi encontrada."
    change_box = [int(xs.min()), int(ys.min()), int(xs.max()) + 1, int(ys.max()) + 1]
    allowed_box = [min(old_box[0], new_box[0]), min(old_box[1], new_box[1]),
                   max(old_box[2], new_box[2]), max(old_box[3], new_box[3])]
    outside = changed.copy()
    outside[allowed_box[1]:allowed_box[3], allowed_box[0]:allowed_box[2]] = False
    assert not outside.any(), "Houve alteração de pixels fora do título."

    manifest = json.loads(MANIFEST.read_text())
    target = next(item for item in manifest["figures"]
                  if item["file"] == str(OUTPUT.relative_to(REPO)))
    target.update({
        "source": str(SOURCE_FIGURE.relative_to(REPO)),
        "source_sha256": sha256(SOURCE_FIGURE),
        "sha256": hashlib.sha256(corrected_bytes).hexdigest(),
        "unchanged_copy": False,
        "presentation_only": True,
        "change": {"old_title": OLD_TITLE, "new_title": NEW_TITLE,
                   "segmentation_recomputed": False,
                   "data_and_crops_unchanged": True},
        "generator": record(Path(__file__).resolve()),
        "source_rendering_script": {**record(SOURCE_SCRIPT), "lines": [250, 265]},
        "inputs": inputs,
        "rendering": {"matplotlib": matplotlib.__version__, "numpy": np.__version__,
                      "dpi": DPI, "size_pixels": [2700, 1800],
                      "first_columns": 2500, "mask_opacity": 0.42},
        "verification": {"original_title_pixels_identical_to_source": True,
                         "changed_pixels": int(changed.sum()),
                         "changed_bbox_xyxy_exclusive": change_box,
                         "allowed_title_bbox_xyxy_exclusive": allowed_box,
                         "changed_pixels_outside_title": 0},
    })
    # Nenhum outro registro do manifesto, dado ou artefato de src2 é modificado.
    OUTPUT.write_bytes(corrected_bytes)
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"file": target["file"], "sha256": target["sha256"],
                      **target["verification"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
