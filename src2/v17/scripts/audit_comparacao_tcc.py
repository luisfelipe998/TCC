"""Confere tabelas numéricas e procedência de figuras do TCC sem alterá-lo."""
from pathlib import Path
import sys
sys.dont_write_bytecode = True
import hashlib
import json
import re
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
PROJECT = BASE.parents[1]
OUT = BASE / "resultados/comparacao_tcc"
OUT.mkdir(parents=True, exist_ok=True)
V16 = PROJECT / "src2/v16_spatial2_18_71/04_estatistica/dados"
checks = []


def record(name, status, details):
    checks.append(dict(verificacao=name, status="PASS" if status else "FAIL", detalhes=details))


def numeric_rows(file, label, ncols):
    text = file.read_text()
    match = re.search(r"\\label\{" + re.escape(label) + r"\}(.*?)\\end\{table\}", text, flags=re.S)
    assert match, label
    result = []
    for line in match[1].splitlines():
        if "&" not in line or not line.rstrip().endswith("\\\\"):
            continue
        cells = line.rstrip().removesuffix("\\\\").split("&")[-ncols:]
        try:
            vals = [float(v.replace("{,}", ",").replace(",", ".").replace("$", "").strip()) for v in cells]
            if len(vals) == ncols:
                result.append(vals)
        except ValueError:
            continue
    return np.array(result)


def compare_table(file, label, expected):
    expected = np.array(expected, float)
    actual = numeric_rows(file, label, expected.shape[1])
    same_shape = actual.shape == expected.shape
    delta = float(np.max(np.abs(actual - expected))) if same_shape else None
    # TCC apresenta quatro casas decimais; comparar arredondamento, não precisão interna.
    ok = same_shape and np.allclose(actual, expected, atol=0.00005001, rtol=0)
    record(label, bool(ok), {"celulas": int(expected.size), "erro_maximo": delta,
                           "forma_tcc": list(actual.shape), "forma_fonte": list(expected.shape)})


def main():
    pairs = pd.read_csv(V16 / "pares_pearson_wco.csv")
    index = pd.read_csv(V16 / "indice_ensaios.csv")
    res = PROJECT / "TCC/resultados.tex"
    disc = PROJECT / "TCC/discussao.tex"
    table_index = re.findall(r"^\s*(\d+)\s*&\s*(m\d+q\d+)\s*&.*?&\s*(\d+)\\\\", res.read_text(), re.M)
    expected_index = [(str(r.numero), r.ensaio.removesuffix("_001"), str(r.comprimento_px))
                      for r in index.itertuples()]
    record("tab:res-ensaios", table_index == expected_index, {"linhas": len(table_index)})
    example_rows = []
    for a, b in [("m04q20_001", "m07q40_001"), ("m06q40_001", "m07q80_001"),
                 ("m05q01_001", "m06q01_001"), ("m04q01_001", "m05q05_001")]:
        p = pairs[(pairs.ensaio_x.isin([a, b])) & (pairs.ensaio_y.isin([a, b]))].iloc[0]
        example_rows.append([p.comprimento_comum, p.pearson_r, p.escore_piso05])
    compare_table(res, "tab:res-exemplos", example_rows)
    names = ["churn", "wavy annular", "plug", "slug", "dispersed bubble"]
    intra, inter = [], []
    for i, a in enumerate(names):
        for b in names[i:]:
            if a == b:
                d = pairs[(pairs.padrao_x == a) & (pairs.padrao_y == a)]
                if len(d):
                    intra.append([sum(index.padrao == a), len(d), d.pearson_r.mean(), d.escore_piso05.mean()])
            else:
                d = pairs[(pairs.padrao_x.isin([a, b])) & (pairs.padrao_y.isin([a, b])) & (pairs.tipo == "inter")]
                inter.append([len(d), d.pearson_r.mean(), d.escore_piso05.mean()])
    compare_table(res, "tab:res-intra", intra)
    compare_table(res, "tab:res-inter", inter)
    composition = pd.read_csv(BASE / "resultados/exploracao/efeito_composicao_tcc.csv")
    global_rows = []
    for metric in ["pearson_r", "escore_piso05"]:
        for weight in ["media", "media_balanceada"]:
            d = composition[(composition.conjunto == "26_todos") & (composition.metrica == metric)].set_index("tipo")
            global_rows.append([d.loc[k, weight] for k in ["intra", "inter", "intra_menos_inter"]])
    compare_table(res, "tab:res-global", global_rows)
    floors = []
    for floor in ["escore_piso01", "escore_piso05", "escore_piso10"]:
        floors.append([pairs[pairs.tipo == k][floor].mean() for k in ["intra", "inter"]])
    compare_table(disc, "tab:disc-pisos", floors)
    controls = pd.read_csv(V16 / "controles_m03q40.csv")
    compare_table(disc, "tab:disc-controles", controls[["pearson_m03q40_outros", "pearson_entre_outros"]].values)
    exclusions = []
    for metric in ["pearson_r", "escore_piso05"]:
        for group, n in [("26_todos", 26), ("25_sem_m03q40", 25)]:
            d = composition[(composition.conjunto == group) & (composition.metrica == metric)].set_index("tipo")
            exclusions.append([n] + [d.loc[k, "media"] for k in ["intra", "inter", "intra_menos_inter"]])
    compare_table(disc, "tab:disc-exclusao", exclusions)
    figures = []
    for entry in json.loads((PROJECT / "TCC/figuras/v16/origens.json").read_text()):
        sha = lambda p: hashlib.file_digest(p.open("rb"), "sha256").hexdigest()
        source, destination = sha(PROJECT / entry["origem"]), sha(PROJECT / entry["destino"])
        figures.append({**entry, "sha256_atual_origem": source, "sha256_atual_destino": destination,
                        "confere": source == destination == entry["sha256"]})
    record("figuras_v16", all(v["confere"] for v in figures), {"arquivos": len(figures)})
    (OUT / "figuras_sha256.json").write_text(json.dumps(figures, ensure_ascii=False, indent=2) + "\n")
    void = pd.read_csv(PROJECT / "src2/v14.1_new_voidfraction/results/void_fraction/panorama_void_per_run.csv")
    e8 = pd.read_csv(PROJECT / "src2/estatistica8/dados/inventario.csv")
    coverage = []
    for row in index.itertuples():
        v = void[void.run == row.ensaio].iloc[0]
        label = v.label.replace("dispersed bubbles", "dispersed bubble")
        coverage.append(dict(ensaio=row.ensaio, padrao_tcc=row.padrao, v14_selecionado=True,
                             v14_validez=v.validity, v14_rotulo_confere=label == row.padrao,
                             estatistica8_selecionado=row.ensaio in set(e8.ensaio)))
    record("coortes", set(void.run) == set(index.ensaio) and set(e8.ensaio) == set(index.ensaio) - {"m07q01_001"}
           and all(r["v14_rotulo_confere"] for r in coverage),
           {"tcc_v16": len(index), "v14_selecionados": len(void), "estatistica8": len(e8),
            "v14_validade": void.validity.value_counts().to_dict()})
    pd.DataFrame(coverage).to_csv(OUT / "coortes.csv", index=False)
    (OUT / "verificacoes.json").write_text(json.dumps(checks, ensure_ascii=False, indent=2) + "\n")
    for c in checks:
        print(c)
    assert all(c["status"] == "PASS" for c in checks)


if __name__ == "__main__":
    main()
