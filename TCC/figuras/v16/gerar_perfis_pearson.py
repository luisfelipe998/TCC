"""Reapresenta os perfis espectrais já calculados na run v16.

Execução: python3 TCC/figuras/v16/gerar_perfis_pearson.py

O desenho foi adaptado de scripts/estatistica.py, função render. Os dados são
lidos de perfis_espectrais.csv e conferidos contra descritores.npz, utilizado
na figura original. Não são recalculados sinais, transformadas ou descritores.
Somente a potência da coluna esquerda é dividida pelo máximo para exibição,
como no original. As medianas também reproduzem a apresentação original.
"""

from pathlib import Path
import csv
import hashlib
import json
import os
import sys

sys.dont_write_bytecode = True
os.environ.setdefault("MPLCONFIGDIR", "/tmp/tcc-perfis-pearson-mpl")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.text import Text
from matplotlib.ticker import FuncFormatter, NullLocator
import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
RUN = REPO / "src2/v16_spatial2_18_71"
CSV = RUN / "04_estatistica/dados/perfis_espectrais.csv"
NPZ = RUN / "intermediario/descritores.npz"
OUTPUT = HERE / "perfis_espectrais_pearson.png"
MANIFEST = HERE / "origens_perfis_pearson.json"
GROUPS = [
    ("churn", "churn", "#d34a33", 10),
    ("wavy annular", "Anular ondulado", "#2377b4", 9),
    ("plug", "plug", "#328652", 3),
    ("slug", "slug", "#9466ac", 3),
    ("dispersed bubble", "Bolhas dispersas", "#b0851c", 1),
]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record(path):
    return {"arquivo": str(path.relative_to(REPO)), "sha256": sha256(path)}


def array_sha256(values):
    return hashlib.sha256(np.asarray(values, dtype="<f8").tobytes()).hexdigest()


def decimal(value, _position=None):
    return f"{value:g}".replace(".", ",")


def main():
    inputs = [CSV, NPZ, RUN / "scripts/estatistica.py", RUN / "scripts/figuras.py",
              HERE / "perfis_comparados.png", HERE / "origens.json"]
    input_records = [record(path) for path in inputs]
    original_entries = json.loads((HERE / "origens.json").read_text(encoding="utf-8"))
    original_entry = next(entry for entry in original_entries
                          if entry["destino"] == "TCC/figuras/v16/perfis_comparados.png")
    assert sha256(HERE / "perfis_comparados.png") == original_entry["sha256"]
    assert sha256(REPO / original_entry["origem"]) == original_entry["sha256"]
    with CSV.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    runs = list(dict.fromkeys(row["ensaio"] for row in rows))
    assert len(rows) == 5200 and len(runs) == 26
    grouped_rows = [[row for row in rows if row["ensaio"] == run] for run in runs]
    assert all(len(group) == 200 for group in grouped_rows)
    assert [row for group in grouped_rows for row in group] == rows
    patterns = [group[0]["padrao"] for group in grouped_rows]
    assert all(len({row["padrao"] for row in group}) == 1 for group in grouped_rows)
    waves_by_run = np.array([[float(row["lambda_px"]) for row in group] for group in grouped_rows])
    waves = waves_by_run[0]
    power = np.array([[float(row["wps_zero"]) for row in group] for group in grouped_rows])
    unit = np.array([[float(row["descritor_unitario"]) for row in group] for group in grouped_rows])
    np.testing.assert_array_equal(waves_by_run, np.tile(waves, (26, 1)))
    # O CSV e os dados utilizados pelo gerador original coincidem exatamente.
    with np.load(NPZ) as saved:
        assert saved["runs"].tolist() == runs
        np.testing.assert_array_equal(waves, saved["wavelengths"])
        np.testing.assert_array_equal(power, saved["power"])
        np.testing.assert_array_equal(unit, saved["unit"])
    assert np.isfinite(power).all() and np.isfinite(unit).all()
    display_power = power / power.max(axis=1, keepdims=True)

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10,
        "axes.labelsize": 9, "xtick.labelsize": 9, "ytick.labelsize": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": .65, "savefig.facecolor": "white",
        "axes.titleweight": "medium", "axes.titlepad": 6,
    })
    fig, axes = plt.subplots(5, 2, figsize=(6.3, 8.5))
    fig.suptitle("Perfis espectrais e comparação de Pearson", y=.99, fontsize=11)
    fig.text(.32, .945, "Potência wavelet (WPS)\nno atraso zero", ha="center", va="center", fontsize=9)
    fig.text(.78, .945, "Perfis transformados\nusados no coeficiente\nde Pearson r",
             ha="center", va="center", fontsize=9)
    plotted_groups = []
    for i, (pattern, label, color, expected_count) in enumerate(GROUPS):
        indices = [j for j, name in enumerate(patterns) if name == pattern]
        assert len(indices) == expected_count
        group_manifest = {"padrao_original": pattern, "rotulo": label, "cor": color,
                          "ensaios_na_ordem": [runs[j] for j in indices], "paineis": []}
        for column, data, ylabel in (
            (0, display_power[indices], "WPS / máximo"),
            (1, unit[indices], "Resíduo\nnormalizado"),
        ):
            ax = axes[i, column]
            for profile in data:
                line, = ax.plot(waves, profile, color=color, alpha=.35, lw=.8)
                np.testing.assert_array_equal(line.get_xdata(), waves)
                np.testing.assert_array_equal(line.get_ydata(), profile)
            median = np.median(data, axis=0)
            line, = ax.plot(waves, median, color=color, lw=1.6)
            np.testing.assert_array_equal(line.get_ydata(), median)
            ax.set_title(label + f" · n = {len(indices)}" if column == 0 else "", loc="left", fontsize=9)
            ax.set(xscale="log", xlim=(50, 3000), ylabel=ylabel)
            ax.set_xticks([50, 300, 3000], ["50", "300", "3000"])
            ax.xaxis.set_minor_locator(NullLocator())
            ax.yaxis.set_major_formatter(FuncFormatter(decimal))
            if i == 4:
                ax.set_xlabel("Comprimento de onda λ (px)")
            ax.grid(color="#edf0f0", lw=.5)
            group_manifest["paineis"].append({
                "coluna": column, "valores_sha256_float64_le": array_sha256(data),
                "mediana_sha256_float64_le": array_sha256(median),
                "limites_x": list(ax.get_xlim()), "limites_y": list(ax.get_ylim()),
            })
        plotted_groups.append(group_manifest)
    fig.subplots_adjust(left=.14, right=.97, top=.885, bottom=.065, wspace=.47, hspace=.85)
    fig.canvas.draw()
    bbox = fig.get_tightbbox(fig.canvas.get_renderer()).padded(.06)
    scale = 16 / (bbox.width * 2.54)
    texts = [item for item in fig.findobj(Text) if item.get_visible() and item.get_text()]
    min_font = min(item.get_fontsize() for item in texts) * scale
    assert min_font >= 8
    assert bbox.height * 2.54 * scale <= 24
    fig.savefig(OUTPUT, dpi=300, bbox_inches="tight", pad_inches=.06)
    plt.close(fig)

    assert [record(path) for path in inputs] == input_records
    manifest = {
        "finalidade": "Esclarecer os rótulos dos perfis usados na comparação pelo coeficiente de Pearson r.",
        "script": record(Path(__file__).resolve()),
        "entradas": input_records,
        "saida": record(OUTPUT),
        "substitui_no_texto": original_entry["destino"],
        "figura_original_preservada": original_entry,
        "algoritmo_de_apresentacao": "Adaptado da função render de scripts/estatistica.py. Estilo de scripts/figuras.py.",
        "reprocessamento_cientifico": False,
        "colunas_utilizadas": ["ensaio", "padrao", "lambda_px", "wps_zero", "descritor_unitario"],
        "coluna_esquerda": "wps_zero dividido pelo máximo de cada ensaio, apenas para exibição",
        "coluna_direita": "descritor_unitario lido diretamente do CSV, sem nova transformação",
        "linhas_finas": "26 ensaios na mesma ordem e com as mesmas cores da figura original",
        "linhas_espessas": "Mediana ponto a ponto de cada padrão, como na figura original",
        "grupos": plotted_groups,
        "wavelengths_sha256_float64_le": array_sha256(waves),
        "largura_insercao_cm": 16,
        "fonte_minima_inserida_pt": min_font,
        "bibliotecas": {"numpy": np.__version__, "matplotlib": matplotlib.__version__},
        "verificacoes": ["5200 linhas de CSV, 26 ensaios e 200 comprimentos de onda por ensaio",
                         "Ordem de ensaios e valores do CSV idênticos aos do NPZ original",
                         "Todos os perfis foram desenhados com os valores originais",
                         "Grupos, cores, medianas, limites e disposição dos painéis preservados",
                         "CSV, NPZ, scripts da run, figura anterior e origens.json mantidos com SHA-256 idêntico"],
    }
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"saida": manifest["saida"], "ensaios": len(runs),
                      "fonte_minima_inserida_pt": min_font}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
