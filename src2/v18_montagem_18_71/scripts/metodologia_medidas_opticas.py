#!/usr/bin/env python3
"""Esquema vetorial das medidas ópticas sobre uma máscara sintética pronta.

Não executa segmentação nem utiliza dados experimentais. Executar este arquivo
recria somente medidas_opticas.pdf/.png, medidas_opticas_dados.npz e
medidas_opticas_origem.json em OUTPUT_ROOT/05_figuras/metodologia.
"""

from __future__ import annotations

import hashlib
import io
import json
import os
from pathlib import Path
import platform
import zipfile

os.environ.setdefault("MPLCONFIGDIR", "/tmp/tcc-medidas-opticas-mpl")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
import numpy as np


from runtime import OUTPUT_ROOT

OUT = OUTPUT_ROOT / "05_figuras/metodologia"
WIDTH_CM = 16.0
HEIGHT_CM = 12.0
GAS = "#247b9b"
OTHER = "#e4e7e9"
INK = "#25343c"
MUTED = "#52616b"
HIGHLIGHT = "#c34a25"
TUBE_DIAMETER_MM = 1.0


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stable_npz(path: Path, arrays: dict[str, np.ndarray]) -> None:
    """Grava um arquivo NumPy com ordem e timestamps de ZIP determinísticos."""
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_STORED) as archive:
        for name in sorted(arrays):
            buffer = io.BytesIO()
            np.lib.format.write_array(buffer, np.asarray(arrays[name]), allow_pickle=False)
            entry = zipfile.ZipInfo(name + ".npy", date_time=(1980, 1, 1, 0, 0, 0))
            entry.compress_type = zipfile.ZIP_STORED
            entry.external_attr = 0o600 << 16
            archive.writestr(entry, buffer.getvalue())


def build_data() -> dict[str, np.ndarray]:
    height, width = 54, 12
    intervals = np.array([[8, 24], [33, 43]], dtype=np.int64)
    row_pattern = np.array([0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0], dtype=np.uint8)
    mask = np.zeros((height, width), dtype=np.uint8)
    for start, stop in intervals:
        mask[start:stop, :] = row_pattern
    column, row = 2, 18
    diameters = intervals[:, 1] - intervals[:, 0]
    centers = np.column_stack((np.zeros(2), height / 2 - intervals.mean(axis=1)))
    areas = np.pi * (diameters / 2.0) ** 2
    tube_area = np.pi * (height / 2.0) ** 2
    alpha_column = float(np.sum(areas) / tube_area)
    alpha_by_column = alpha_column * row_pattern.astype(float)
    transitions = mask[:, 1:] != mask[:, :-1]
    row_transitions = transitions[row]

    assert mask.shape == (54, 12)
    assert np.array_equal(mask[row], row_pattern)
    assert np.array_equal(np.flatnonzero(row_transitions), [1, 4, 6, 8])
    assert int(row_transitions.sum()) == 4
    assert width - 1 == 11
    assert int(transitions.sum()) == 104
    assert np.isclose(alpha_column, (16 / 54) ** 2 + (10 / 54) ** 2)
    assert not np.isclose(alpha_column, mask[:, column].mean())
    assert np.all(np.linalg.norm(centers, axis=1) + diameters / 2 <= height / 2)
    assert np.linalg.norm(centers[0] - centers[1]) >= diameters.sum() / 2

    return {
        "mask": mask,
        "illustrated_column_index": np.array(column),
        "illustrated_row_index": np.array(row),
        "illustrated_column": mask[:, column],
        "illustrated_row": mask[row],
        "gas_intervals_start_stop": intervals,
        "diameters_px": diameters,
        "tube_diameter_px": np.array(height),
        "gas_circle_centers_px": centers,
        "gas_areas_px2": areas,
        "tube_area_px2": np.array(tube_area),
        "alpha_by_column": alpha_by_column,
        "alpha_mean": np.array(alpha_by_column.mean()),
        "alpha_illustrated_column": np.array(alpha_column),
        "horizontal_transitions": transitions,
        "transition_count": np.array(transitions.sum()),
        "neighbor_pair_count": np.array(height * (width - 1)),
        "flip_density": np.array(transitions.mean()),
        "row_transitions": row_transitions,
        "row_transition_count": np.array(row_transitions.sum()),
        "row_neighbor_pair_count": np.array(width - 1),
    }


def draw_figure(data: dict[str, np.ndarray]) -> tuple[Path, Path]:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    })
    fig = plt.figure(figsize=(WIDTH_CM / 2.54, HEIGHT_CM / 2.54))
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set(xlim=(0, WIDTH_CM), ylim=(0, HEIGHT_CM), aspect="equal")
    ax.axis("off")

    def label(x: float, y: float, text: str, **kwargs) -> None:
        ax.text(x, y, text, va="center", **{"color": INK, **kwargs})

    def arrow(x0: float, y0: float, x1: float, y1: float, **kwargs) -> None:
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops={
            "arrowstyle": "-|>", "color": INK, "lw": 1.0,
            "shrinkA": 0, "shrinkB": 0, **kwargs,
        })

    def draw_mask(x: float, top: float, width: float, height: float,
                  *, selected_column: bool) -> None:
        """Desenha todos os pixels, ampliados anisotropicamente para leitura."""
        mask = data["mask"]
        pixel_w, pixel_h = width / mask.shape[1], height / mask.shape[0]
        ax.add_patch(Rectangle((x, top-height), width, height,
                               facecolor=OTHER, edgecolor="none"))
        # Cada bloco continua exatamente correspondente às células da máscara.
        for j in range(mask.shape[1]):
            starts = np.flatnonzero(np.diff(np.r_[0, mask[:, j], 0].astype(int)) == 1)
            stops = np.flatnonzero(np.diff(np.r_[0, mask[:, j], 0].astype(int)) == -1)
            for start, stop in zip(starts, stops):
                ax.add_patch(Rectangle((x+j*pixel_w, top-stop*pixel_h),
                                       pixel_w, (stop-start)*pixel_h,
                                       facecolor=GAS, edgecolor="none"))
        for j in range(1, mask.shape[1]):
            ax.plot([x+j*pixel_w]*2, [top-height, top], color="white", lw=0.35)
        ax.add_patch(Rectangle((x, top-height), width, height,
                               facecolor="none", edgecolor=INK, lw=0.8))
        if selected_column:
            column = int(data["illustrated_column_index"])
            box = (x+column*pixel_w, top-height, pixel_w, height)
        else:
            row = int(data["illustrated_row_index"])
            box = (x, top-(row+1)*pixel_h, width, pixel_h)
        ax.add_patch(Rectangle(box[:2], box[2], box[3], facecolor="none",
                               edgecolor=HIGHLIGHT, lw=1.7, zorder=8))

    label(0.35, 11.60, "a) Da máscara à estimativa de área", fontsize=10.5, weight="bold")
    ax.add_patch(Rectangle((0.35, 10.95), 0.25, 0.25, facecolor=GAS, edgecolor="none"))
    label(0.73, 11.07, "1: gás inferido", fontsize=8.5)
    ax.add_patch(Rectangle((3.55, 10.95), 0.25, 0.25, facecolor=OTHER, edgecolor=MUTED, lw=0.4))
    label(3.93, 11.07, "0: demais pixels", fontsize=8.5)
    label(15.65, 11.07, "Máscara sintética, sem dados experimentais", ha="right", fontsize=8.5)

    for x, title in ((2.12, "1  Máscara"), (5.22, "2  Coluna extraída"),
                     (8.72, "3  Seções circulares"), (13.38, "4  Razão de áreas")):
        label(x, 10.35, title, ha="center", fontsize=8.8, weight="bold")
    label(2.12, 9.91, "54 linhas × 12 colunas", ha="center", fontsize=8.5)

    height = int(data["tube_diameter_px"])
    scale = 3.1 / height
    top = 9.42
    bottom = top - height * scale
    draw_mask(0.60, top, 3.05, 3.1, selected_column=True)
    label(2.12, 5.98, "Coluna destacada", ha="center", color=HIGHLIGHT, fontsize=8.5)
    arrow(3.78, 7.87, 4.28, 7.87)

    column_left, column_width = 4.47, 0.40
    ax.add_patch(Rectangle((column_left, bottom), column_width, height * scale,
                           facecolor=OTHER, edgecolor=HIGHLIGHT, lw=1.0))
    for index, (start, stop) in enumerate(data["gas_intervals_start_stop"], 1):
        y_bottom = top - stop * scale
        y_top = top - start * scale
        ax.add_patch(Rectangle((column_left, y_bottom), column_width, y_top - y_bottom,
                               facecolor=GAS, edgecolor="none"))
        bracket_x = 5.10
        ax.plot([bracket_x, bracket_x], [y_bottom, y_top], color=INK, lw=0.8)
        ax.plot([bracket_x - 0.08, bracket_x + 0.08], [y_bottom, y_bottom], color=INK, lw=0.8)
        ax.plot([bracket_x - 0.08, bracket_x + 0.08], [y_top, y_top], color=INK, lw=0.8)
        subscript = "₁" if index == 1 else "₂"
        label(5.27, (y_bottom + y_top) / 2,
              f"Largura d{subscript}\n{stop - start} pixels", fontsize=8.5, linespacing=1.15)

    circle_x, circle_y, radius = 8.72, (top + bottom) / 2, 1.43
    circle_scale = radius / (height/2)
    arrow(6.68, circle_y, 7.15, circle_y)
    ax.add_patch(Circle((circle_x, circle_y), radius, facecolor=OTHER, edgecolor=INK, lw=0.9))
    for index, (center, diameter) in enumerate(zip(data["gas_circle_centers_px"], data["diameters_px"]), 1):
        gas_center = (circle_x+center[0]*circle_scale, circle_y+center[1]*circle_scale)
        ax.add_patch(Circle(gas_center, diameter*circle_scale/2,
                            facecolor=GAS, edgecolor="none"))
        label(*gas_center, "d₁" if index == 1 else "d₂", color="white", ha="center", fontsize=8.5)
    label(circle_x, 5.99, "Diâmetro D = 1 mm", ha="center", fontsize=8.5)
    label(circle_x, 5.57, "Representação: 54 px", ha="center", fontsize=8.5)
    arrow(10.31, circle_y, 10.77, circle_y)
    label(13.38, 9.04, "Soma das áreas\ninferidas de gás", ha="center", fontsize=8.8, linespacing=1.25)
    ax.plot([11.25, 15.50], [8.40, 8.40], color=MUTED, lw=0.7)
    label(13.38, 8.03, "Área interna do tubo", ha="center", fontsize=8.8)
    label(13.38, 7.15, r"$(16/54)^2 + (10/54)^2$", ha="center", fontsize=10.0)
    alpha = float(data["alpha_illustrated_column"])
    alpha_text = f"{alpha:.4f}".replace(".", "{,}")
    label(13.38, 6.65, "Fração de vazio estimada", ha="center", fontsize=8.5)
    label(13.38, 6.18, r"$\widehat{\alpha}_{\mathrm{opt},j} \approx " + alpha_text + "$",
          ha="center", fontsize=10.0)
    label(13.38, 5.35, "Depois: média das\n12 colunas", ha="center", fontsize=8.5, linespacing=1.25)
    label(0.35, 4.73, "Desenho esquemático, sem escala.",
          fontsize=8.5)

    ax.plot([0.35, 15.65], [4.37, 4.37], color="#c8ced2", lw=0.65)
    label(0.35, 3.98, "b) Da linha da máscara à contagem de transições", fontsize=10.5, weight="bold")
    draw_mask(0.60, 3.44, 3.05, 2.55, selected_column=False)
    label(2.12, 0.44, "Mesma máscara:\nlinha destacada", ha="center", fontsize=8.5, linespacing=1.15)
    label(4.70, 3.44, "Sequência de pixels da linha destacada", fontsize=8.5)
    x0, y0, cell = 4.70, 2.25, 0.89
    cell_height = 0.67
    for j, value in enumerate(data["illustrated_row"]):
        ax.add_patch(Rectangle((x0 + j * cell, y0), cell, cell_height,
                               facecolor=GAS if value else OTHER, edgecolor="white", lw=1.0))
        label(x0 + (j + 0.5) * cell, y0 + cell_height / 2, str(int(value)),
              ha="center", fontsize=10, color="white" if value else INK)
    for j in np.flatnonzero(data["row_transitions"]):
        boundary = x0 + (j + 1) * cell
        arrow(boundary, 3.20, boundary, y0 + cell_height + 0.03, color=HIGHLIGHT)
    row_center_y = 3.44-(int(data["illustrated_row_index"])+0.5)*2.55/height
    arrow(3.78, row_center_y, 4.50, y0+cell_height/2, color=HIGHLIGHT)
    label(10.04, 1.80, "4 mudanças em 11 pares de pixels vizinhos", ha="center", fontsize=9.0)
    label(10.04, 1.28, "Máscara inteira: 104 mudanças em 54 × 11 pares", ha="center", fontsize=8.5)
    flip_text = f"{float(data['flip_density']):.4f}".replace(".", "{,}")
    label(10.04, 0.82, r"Densidade de transições $F = 104/594 \approx " + flip_text + "$",
          ha="center", fontsize=9.0)
    label(10.04, 0.36, "A contagem foi somada nas 54 linhas.", ha="center", fontsize=8.5)

    # Verificação material dos rótulos, independente do dimensionamento do PNG.
    assert min(text.get_fontsize() for text in ax.texts) >= 8.5
    fig.canvas.draw()
    canvas_bbox = fig.bbox
    for text in ax.texts:
        bbox = text.get_window_extent(fig.canvas.get_renderer())
        assert canvas_bbox.contains(bbox.x0, bbox.y0) and canvas_bbox.contains(bbox.x1, bbox.y1), text.get_text()
    pdf = OUT / "medidas_opticas.pdf"
    png = OUT / "medidas_opticas.png"
    fig.savefig(pdf, metadata={"Creator": "Matplotlib", "CreationDate": None, "ModDate": None})
    fig.savefig(png, dpi=300, metadata={"Software": "Matplotlib"})
    plt.close(fig)
    return pdf, png


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    data = build_data()
    npz = OUT / "medidas_opticas_dados.npz"
    stable_npz(npz, data)
    pdf, png = draw_figure(data)
    metadata = {
        "titulo": "Leitura da máscara: áreas e transições (exemplo sintético)",
        "tipo": "esquema didático com máscara binária sintética pronta",
        "dados_experimentais": False,
        "segmentacao_executada": False,
        "dimensoes_cm": [WIDTH_CM, HEIGHT_CM],
        "fonte_minima_pt": 8.5,
        "mascara_shape": list(data["mask"].shape),
        "coluna_ilustrada_indice_base_zero": int(data["illustrated_column_index"]),
        "linha_ilustrada_indice_base_zero": int(data["illustrated_row_index"]),
        "intervalos_gas_inicio_inclusivo_fim_exclusivo": data["gas_intervals_start_stop"].tolist(),
        "diametros_gas_px": data["diameters_px"].tolist(),
        "diametro_tubo_px": int(data["tube_diameter_px"]),
        "diametro_interno_tubo_arranjo_mm": TUBE_DIAMETER_MM,
        "interpretacao_dimensoes": "O tubo do arranjo possui diâmetro interno de 1 mm. O exemplo sintético adota 54 pixels para sua representação radial. Essa associação esquemática não fornece calibração óptica dos pixels experimentais.",
        "representacao_mascara": "A máscara inteira de 54 por 12 células é mostrada nos dois painéis, com ampliação anisotrópica para leitura. A coluna e a linha destacadas correspondem aos índices 2 e 18, respectivamente, em base zero; esses índices não aparecem nos rótulos da figura. Desenho sem escala física.",
        "painel_a": "Máscara completa com coluna destacada; extração dos dois intervalos; estimativa por seções circulares dentro do tubo circular; razão entre as áreas; indicação da média posterior entre colunas.",
        "painel_b": "A mesma máscara com linha destacada; sequência binária correspondente, quatro fronteiras de mudança e onze pares vizinhos; agregação das contagens nas 54 linhas da máscara.",
        "arredondamento_figura": "A fração de vazio estimada na coluna e a densidade de transições são exibidas com quatro casas decimais e sinal de aproximação. Dados salvos conservam a precisão calculada.",
        "centros_circulos_px_relativos_ao_centro_tubo": data["gas_circle_centers_px"].tolist(),
        "hipotese": "Cada intervalo vertical corresponde a uma seção circular. As duas seções estimadas não se sobrepõem e estão contidas no círculo do tubo.",
        "estimador": "Soma das áreas circulares de gás dividida pela área circular do tubo, equivalente à soma dos quadrados das razões entre diâmetros.",
        "alpha_coluna_ilustrada": float(data["alpha_illustrated_column"]),
        "alpha_media_mascara": float(data["alpha_mean"]),
        "fracao_pixels_gas_coluna_nao_usada_como_alpha": float(data["illustrated_column"].mean()),
        "transicoes_linha_ilustrada": int(data["row_transition_count"]),
        "pares_linha_ilustrada": int(data["row_neighbor_pair_count"]),
        "transicoes_mascara": int(data["transition_count"]),
        "pares_mascara": int(data["neighbor_pair_count"]),
        "densidade_transicoes_mascara": float(data["flip_density"]),
        "versoes": {"python": platform.python_version(), "numpy": np.__version__, "matplotlib": matplotlib.__version__},
        "sha256": {path.name: sha256(path) for path in (Path(__file__), npz, pdf, png)},
    }
    (OUT / "medidas_opticas_origem.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"status": "PASS", "alpha_coluna": metadata["alpha_coluna_ilustrada"],
                      "transicoes_linha": metadata["transicoes_linha_ilustrada"],
                      "pares_linha": metadata["pares_linha_ilustrada"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
