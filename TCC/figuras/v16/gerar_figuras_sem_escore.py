"""Apresentação de Pearson sem o painel de escore, a partir do CSV da run v16.

Execução: python3 TCC/figuras/v16/gerar_figuras_sem_escore.py

O painel de Pearson foi adaptado de scripts/estatistica.py, função render.
Foram preservados os valores, a ordem dos pares, as sementes da dispersão
horizontal, os limites verticais, as cores, as médias e os quartis.
A alteração restringiu-se à apresentação. Não são recalculados sinais,
correlações, transformadas, mapas wavelet ou coeficientes de Pearson.
As figuras originais e o manifesto origens.json são apenas lidos.
"""

from pathlib import Path
import csv
import hashlib
import json
import os
import sys

sys.dont_write_bytecode = True
os.environ.setdefault("MPLCONFIGDIR", "/tmp/tcc-figuras-sem-escore-mpl")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.text import Text
from matplotlib.ticker import FuncFormatter
import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
RUN = REPO / "src2/v16_spatial2_18_71"
CSV = RUN / "04_estatistica/dados/pares_pearson_wco.csv"
OUTPUT = HERE / "correlacoes_dentro_entre_pearson.png"
MANIFEST = HERE / "origens_sem_escore.json"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record(path):
    return {"arquivo": str(path.relative_to(REPO)), "sha256": sha256(path)}


def decimal(value, _position=None):
    return f"{value:g}".replace(".", ",")


def main():
    # Verifica toda a coleção anterior sem sobrescrever nenhum arquivo.
    original_manifest = HERE / "origens.json"
    original_entries = json.loads(original_manifest.read_text(encoding="utf-8"))
    for entry in original_entries:
        assert sha256(REPO / entry["origem"]) == entry["sha256"]
        assert sha256(REPO / entry["destino"]) == entry["sha256"]
    inputs = [CSV, RUN / "scripts/estatistica.py", RUN / "scripts/figuras.py",
              original_manifest]
    input_records = [record(path) for path in inputs]
    with CSV.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 325
    pairs = [tuple(sorted((row["ensaio_x"], row["ensaio_y"]))) for row in rows]
    assert len(set(pairs)) == len(pairs)
    assert len({run for pair in pairs for run in pair}) == 26
    assert {row["tipo"] for row in rows} == {"intra", "inter"}

    # Estilo de figuras.py. O painel de estatistica.py foi separado dos demais.
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10,
        "axes.labelsize": 9, "xtick.labelsize": 9, "ytick.labelsize": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": .65, "savefig.facecolor": "white",
        "axes.titleweight": "medium", "axes.titlepad": 6,
    })
    fig, ax = plt.subplots(figsize=(6.3, 4.0))
    fig.suptitle("Correlação de Pearson entre perfis espectrais", y=.99, fontsize=11)
    summary = {}
    for x, kind, color, expected_count in (
        (1, "intra", "#087f8c", 87), (2, "inter", "#b66b27", 238)
    ):
        values = np.array([float(row["pearson_r"]) for row in rows if row["tipo"] == kind])
        assert len(values) == expected_count
        assert np.isfinite(values).all() and ((values >= -1) & (values <= 1)).all()
        rng = np.random.default_rng(2718 + x)
        positions = x + rng.uniform(-.17, .17, len(values))
        points = ax.scatter(positions, values, s=13, alpha=.45, color=color, edgecolor="none")
        np.testing.assert_array_equal(points.get_offsets()[:, 1], values)
        np.testing.assert_array_equal(points.get_offsets()[:, 0], positions)
        q1, q3 = np.quantile(values, [.25, .75])
        mean = float(values.mean())
        ax.plot([x-.23, x+.23], [mean]*2, color="#182d37", lw=2)
        ax.plot([x, x], [q1, q3], color="#182d37", lw=2)
        summary[kind] = {"pares": len(values), "media": mean,
                         "primeiro_quartil": float(q1), "terceiro_quartil": float(q3),
                         "semente_disposicao_horizontal": 2718 + x}
    ax.set_xticks([1, 2], ["Mesmo padrão\n87 pares", "Padrões diferentes\n238 pares"], fontsize=9)
    ax.set(xlim=(.45, 2.55), ylim=(-1, 1), ylabel="Coeficiente de Pearson r")
    ax.yaxis.set_major_formatter(FuncFormatter(decimal))
    ax.grid(axis="y", color="#e6ecef", lw=.6)
    fig.subplots_adjust(left=.14, right=.97, top=.89, bottom=.23)
    fig.text(.53, .025, "Traço horizontal: média · barra vertical: intervalo interquartil",
             ha="center", fontsize=9)
    fig.canvas.draw()
    bbox = fig.get_tightbbox(fig.canvas.get_renderer()).padded(.06)
    insertion_scale = 16 / (bbox.width * 2.54)
    texts = [item for item in fig.findobj(Text) if item.get_visible() and item.get_text()]
    labels = [item.get_text() for item in texts]
    assert not any(term in label.casefold() for label in labels for term in ("escore", "piso", "wco"))
    min_font = min(item.get_fontsize() for item in texts) * insertion_scale
    assert min_font >= 8
    fig.savefig(OUTPUT, dpi=300, bbox_inches="tight", pad_inches=.06)
    plt.close(fig)

    # Todos os insumos e as figuras anteriores devem continuar idênticos.
    assert [record(path) for path in inputs] == input_records
    for entry in original_entries:
        assert sha256(REPO / entry["origem"]) == entry["sha256"]
        assert sha256(REPO / entry["destino"]) == entry["sha256"]
    preserved_names = {
        "plug_m05q01_001.png", "dispersed_bubble_m07q01_001.png",
        "plug_m05q01_001_m06q01_001.png", "slug_m04q01_001_m05q05_001.png",
        "slug_churn_m04q01_001_m07q20_001.png", "perfis_comparados.png",
        "matriz_pearson_todos.png", "churn_m04q20_001_m07q40_001.png",
        "wavy_annular_m06q40_001_m07q80_001.png",
    }
    manifest = {
        "finalidade": "Apresentação da distribuição do coeficiente de Pearson r sem o painel de escore.",
        "script": record(Path(__file__).resolve()),
        "entradas": input_records,
        "algoritmo_de_apresentacao": "Adaptado da função render em scripts/estatistica.py. Estilo de scripts/figuras.py.",
        "reprocessamento_cientifico": False,
        "colunas_utilizadas": ["ensaio_x", "ensaio_y", "tipo", "pearson_r"],
        "saida": record(OUTPUT),
        "substitui_no_texto": "TCC/figuras/v16/correlacoes_dentro_entre.png",
        "resumo_dos_pontos": summary,
        "largura_insercao_cm": 16,
        "fonte_minima_inserida_pt": min_font,
        "bibliotecas": {"numpy": np.__version__, "matplotlib": matplotlib.__version__},
        "figuras_preservadas_sem_escore": [entry for entry in original_entries
                                                 if Path(entry["destino"]).name in preserved_names],
        "figuras_originais_preservadas": len(original_entries),
        "verificacoes": ["325 pares únicos de 26 ensaios", "87 pares intra e 238 pares inter",
                         "Todos os valores de Pearson do CSV foram desenhados sem alteração",
                         "Sementes de disposição horizontal iguais às do gráfico original",
                         "Nenhum rótulo de escore, piso ou WCO no novo gráfico",
                         "Insumos e 12 figuras originais mantidos com SHA-256 idêntico"],
    }
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"saida": manifest["saida"], "pares": summary,
                      "fonte_minima_inserida_pt": min_font}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
