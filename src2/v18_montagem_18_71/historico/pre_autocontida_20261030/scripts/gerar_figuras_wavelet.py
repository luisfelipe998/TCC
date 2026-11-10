"""Aplica aos resultados v18 as funções de apresentação aprovadas na v16.

Os trechos de desenho dos três gráficos estatísticos são carregados por AST.
Seu estilo, dimensões e rótulos ficam preservados, com substituição das entradas.
Não executa processamento de sinais ou o programa principal dos geradores antigos.
"""
import ast
import json
from pathlib import Path
import shutil
import sys
import types

sys.dont_write_bytecode = True
import recalcular_wavelet as calculation
from recalcular_wavelet import (BASE, REPO, OLD, ROOT, INTER, TABLES, WAVES,
    ROW_START, ROW_END, load_image, load_raw, load_module, read_csv, write_json,
    sha, source_record, np)

DEST = REPO / "TCC/figuras/v18_wavelet"
CURRENT = REPO / "TCC/figuras/v16"


def historical_plot_module():
    source = OLD / "scripts/figuras.py"
    tree = ast.parse(source.read_text())
    tree.body = [node for node in tree.body
                 if not (isinstance(node, ast.ImportFrom) and node.module == "common")]
    module = types.ModuleType("figuras")
    module.__dict__.update(ROOT=ROOT, WAVES=WAVES, ROW_START=ROW_START, ROW_END=ROW_END)
    exec(compile(tree, str(source), "exec"), module.__dict__)
    sys.modules["figuras"] = module
    return module


def current_main_drawing(path, updates, return_expression):
    """Extrai o corpo que lê o CSV e desenha, antes dos registros históricos."""
    module = load_module(path, "v18_" + path.stem)
    module.__dict__.update(updates)
    source = path.read_text()
    lines = source.splitlines()
    start = next(i for i, line in enumerate(lines, 1) if line.startswith("    with CSV.open"))
    end = next(i for i, line in enumerate(lines, 1) if line.startswith("    plt.close(fig)"))
    main = next(node for node in ast.parse(source).body if isinstance(node, ast.FunctionDef) and node.name == "main")
    main.name = "render_current_data"
    main.body = [node for node in main.body if node.lineno >= start and node.end_lineno <= end]
    assert main.body and isinstance(main.body[0], ast.With)
    main.body.append(ast.Return(value=ast.parse(return_expression, mode="eval").body))
    tree = ast.fix_missing_locations(ast.Module(body=[main], type_ignores=[]))
    exec(compile(tree, str(path), "exec"), module.__dict__)
    return module.render_current_data(), {"arquivo": str(path.relative_to(REPO)),
                                          "sha256": sha(path), "linhas_desenho": [start, end]}


def matrix_figure(figuras, records):
    source = OLD / "scripts/estatistica.py"
    function = next(node for node in ast.parse(source.read_text()).body
                    if isinstance(node, ast.FunctionDef) and node.name == "render")
    first_for = next(i for i, node in enumerate(function.body) if isinstance(node, ast.For))
    function.name = "matrix_only"
    function.body = function.body[:first_for + 1]
    loop = function.body[-1]
    assert isinstance(loop.iter, ast.List) and len(loop.iter.elts) == 2
    loop.iter.elts = loop.iter.elts[:1]
    function.body.append(ast.Return(value=ast.Name(id="layouts", ctx=ast.Load())))
    scope = dict(np=np, TABLES=TABLES, INTER=INTER, STATS=ROOT / "04_estatistica",
                 records=lambda: records, read_csv=read_csv,
                 combinations=calculation.combinations)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[function], type_ignores=[])),
                 str(source), "exec"), scope)
    return scope["matrix_only"]()


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    records = json.loads((INTER / "aquisicoes.json").read_text())
    summary = json.loads((ROOT / "resumo_tcc.json").read_text())
    assert summary["execucao_pares_wavelet_completa"]
    by_id = {r["ensaio"]: r for r in records}
    originals = [OLD / "scripts/figuras.py", OLD / "scripts/estatistica.py",
                 OLD / "scripts/wavelet_core.py", CURRENT / "gerar_perfis_pearson.py",
                 CURRENT / "gerar_figuras_sem_escore.py"]
    original_hashes = {str(p): sha(p) for p in originals}
    figuras = historical_plot_module()
    outputs, layouts = [], []
    for run in ("m05q01_001", "m07q01_001"):
        row = by_id[run]
        with np.load(INTER / "individuais" / f"{run}.npz") as archive:
            prepared = archive["prepared"].copy()
        layout = figuras.panorama_figure(run, load_image(row), load_raw(run), prepared,
                                        pattern=row["padrao"])
        path = ROOT / layout["arquivo"]
        outputs.append((path, {"tipo": "panorama", "ensaio": run,
                              "entradas": [source_record(REPO / row["wavelet_path"]),
                                           source_record(INTER / "individuais" / f"{run}.npz")]}))
        layouts.append(layout)
    for pair in summary["pares_ilustrados"]:
        a, b = pair["run_x"], pair["run_y"]
        npz_path = INTER / "pares" / f"{a}__{b}.npz"
        with np.load(npz_path) as saved:
            maps = {key: saved[key] for key in saved.files}
        path = ROOT / pair["arquivo"]
        layout = figuras.pair_figure([load_image(by_id[a]), load_image(by_id[b])],
            maps["x"], maps["y"], maps, runs=(a, b), pattern=pair["padrao"], destination=path)
        outputs.append((path, {"tipo": pair["tipo"], "ensaio_x": a, "ensaio_y": b,
                              "pearson_r": pair["pearson_r"], "wco_maxima": pair["wco_maxima"],
                              "comprimento_comum": pair["comprimento"], "entradas": [source_record(npz_path)]}))
        layouts.append(layout)
    layouts.extend(matrix_figure(figuras, records))
    outputs.append((ROOT / "04_estatistica/matriz_pearson_todos.png",
                    {"tipo": "matriz_pearson", "entradas": [source_record(INTER / "descritores.npz")]}))
    profiles, profile_source = current_main_drawing(CURRENT / "gerar_perfis_pearson.py",
        {"CSV": TABLES / "perfis_espectrais.csv", "NPZ": INTER / "descritores.npz",
         "OUTPUT": ROOT / "04_estatistica/perfis_espectrais_pearson.png"},
        "{'grupos': plotted_groups, 'fonte_minima_inserida_pt': min_font}")
    outputs.append((ROOT / "04_estatistica/perfis_espectrais_pearson.png",
                    {"tipo": "perfis", "apresentacao": profiles,
                     "entradas": [source_record(TABLES / "perfis_espectrais.csv"),
                                  source_record(INTER / "descritores.npz")]}))
    distribution, distribution_source = current_main_drawing(CURRENT / "gerar_figuras_sem_escore.py",
        {"CSV": TABLES / "pares_pearson_wco.csv",
         "OUTPUT": ROOT / "04_estatistica/correlacoes_dentro_entre_pearson.png"},
        "{'resumos': summary, 'fonte_minima_inserida_pt': min_font}")
    outputs.append((ROOT / "04_estatistica/correlacoes_dentro_entre_pearson.png",
                    {"tipo": "distribuicoes", "apresentacao": distribution,
                     "entradas": [source_record(TABLES / "pares_pearson_wco.csv")]}))
    entries = []
    for path, metadata in outputs:
        destination = DEST / path.name
        shutil.copyfile(path, destination)
        assert sha(path) == sha(destination)
        entries.append(dict(origem=source_record(path), destino=source_record(destination), **metadata))
    assert len(entries) == 10
    assert all(sha(Path(path)) == digest for path, digest in original_hashes.items())
    manifest = dict(execucao=str(BASE.relative_to(REPO)), script=source_record(Path(__file__).resolve()),
        fontes_historicas_preservadas=[source_record(p) for p in originals],
        desenho_perfis=profile_source, desenho_distribuicoes=distribution_source,
        mesma_apresentacao=True, dados_recalculados=True,
        criterio_intra="maior Pearson de cada padrão com pelo menos dois ensaios",
        criterio_inter="maior coerência pontual dos 238 pares interpadrões",
        figuras=entries, layouts=layouts,
        arquivos_entrada_nao_modificados=True,
        verificacoes="Coeficientes, potências, correlações e máscaras de validade conferidos na execução. Curvas de perfis e pontos da distribuição conferidos contra os arrays lidos.")
    write_json(DEST / "origens.json", manifest)
    write_json(INTER / "auditoria/figuras.json", manifest)
    print(json.dumps({"figuras": len(entries), "manifesto": str(DEST / "origens.json")}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
