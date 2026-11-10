"""Check saved numerical artifacts, report claims, local links and embedded figures.

Does not import the analysis implementation or write to either source study.
"""
from __future__ import annotations

import base64
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np
import pandas as pd
from PIL import Image
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
MAIN = "ancora_2048__hellinger_principal"
FULL = "panorama_completo__hellinger_principal"


def sha(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def rank_A(same, cross):
    n, m = len(same), len(cross)
    ranks = rankdata(np.concatenate([same, cross]), method="average")
    return (ranks[:n].sum() - n*(n+1)/2)/(n*m)


def independent_summary(frame):
    same, cross = frame[frame.mesmo_padrao], frame[~frame.mesmo_padrao]
    balanced = np.mean([rank_A(a.similaridade.to_numpy(), b.similaridade.to_numpy())
        for _, a in same.groupby("grupo") for _, b in cross.groupby("grupo")])
    return rank_A(same.similaridade.to_numpy(), cross.similaridade.to_numpy()), balanced


def main():
    verification = json.loads((OUT / "verificacao.json").read_text())
    assert verification["status"] == "PASS"
    for name, expected in verification["codigo_sha256"].items():
        assert sha(ROOT / name) == expected, f"Código mudou após a análise: {name}"
    frozen = json.loads((OUT / "entradas_congeladas.json").read_text())
    for relative, expected in frozen["sha256"].items():
        assert sha(ROOT.parent / relative) == expected, f"Fonte alterada: {relative}"
    manifest = json.loads((OUT / "artefatos.json").read_text())
    for relative, record in manifest.items():
        path = ROOT / relative
        assert path.stat().st_size == record["bytes"] and sha(path) == record["sha256"], relative
    actual = {str(p.relative_to(ROOT)) for p in ROOT.rglob("*") if p.is_file()
              and p.name not in ("artefatos.json", "verificacao_final.json") and "__pycache__" not in p.parts}
    assert actual == set(manifest)
    summary = pd.read_csv(OUT / "concordancia.csv").set_index("metodo")
    pairs = pd.read_csv(OUT / "pares.csv")
    deletion = pd.read_csv(OUT / "influencia_aquisicoes.csv")
    by_pattern = pd.read_csv(OUT / "por_padrao.csv")
    index = pd.read_csv(OUT / "indice_ensaios.csv").set_index("ensaio")
    windows = pd.read_csv(OUT / "recortes.csv")
    archive = dict(np.load(OUT / "perfis.npz", allow_pickle=False))
    runs = archive["runs"].tolist()
    lookup = {run: i for i, run in enumerate(runs)}
    assert len(runs) == len(set(runs)) == 25
    assert set(runs) == set(index.index)
    assert len(pairs) == 19*300 and len(deletion) == 19*25
    assert len(windows) == 7*25 and len(summary) == 19
    assert (windows.inicio_px >= 0).all()
    assert np.array_equal(windows.fim_exclusivo_px - windows.inicio_px, windows.comprimento_px)
    for row in windows.itertuples():
        assert row.fim_exclusivo_px <= index.loc[row.ensaio, "comprimento_panorama_px"]
    max_area_error, max_H_error, max_A_error = 0., 0., 0.
    for k, condition in enumerate(archive["conditions"]):
        for domain, mask in (("principal", np.ones(200, dtype=bool)), ("conservador", archive["conservative_mask"])):
            axis = archive["wavelengths"][mask]
            p = archive["powers"][k][:, mask]
            mass = archive[f"{condition}__{domain}"]
            area = np.trapezoid(p, x=np.log(axis), axis=1)
            density = p / area[:, None]
            max_area_error = max(max_area_error, float(np.max(abs(np.trapezoid(density, x=np.log(axis), axis=1)-1))))
            assert np.all(np.isfinite(p)) and np.all(p >= 0)
            np.testing.assert_allclose(mass.sum(axis=1), 1., atol=1e-13)
            method = f"{condition}__hellinger_{domain}"
            part = pairs[pairs.metodo == method]
            assert len(part) == 300
            for row in part.itertuples():
                a, b = lookup[row.aquisicao_x], lookup[row.aquisicao_y]
                # Bhattacharyya coefficient; independent of squared vector difference.
                h = np.sqrt(max(0., 1-float(np.sqrt(mass[a]*mass[b]).sum())))
                max_H_error = max(max_H_error, abs(h-row.valor))
    for method, part in pairs.groupby("metodo"):
        assert len(part) == 300 and not part.duplicated(["aquisicao_x", "aquisicao_y"]).any()
        assert set(part.aquisicao_x) | set(part.aquisicao_y) == set(runs)
        assert int(part.mesmo_padrao.sum()) == 87
        for row in part.itertuples():
            assert row.padrao_x == index.loc[row.aquisicao_x, "padrao"]
            assert row.padrao_y == index.loc[row.aquisicao_y, "padrao"]
            assert row.mesmo_padrao == (row.padrao_x == row.padrao_y)
        expected = independent_summary(part)
        observed = summary.loc[method, ["A_pares", "A_balanceada"]].to_numpy(float)
        max_A_error = max(max_A_error, float(np.max(abs(np.asarray(expected)-observed))))
        for run in runs:
            subset = part[(part.aquisicao_x != run) & (part.aquisicao_y != run)]
            assert len(subset) == 276
            expected = independent_summary(subset)
            saved = deletion[(deletion.metodo == method) & (deletion.aquisicao_retirada == run)]
            assert len(saved) == 1
            row = saved.iloc[0]
            assert row.pares_intra + row.pares_entre == 276
            max_A_error = max(max_A_error, float(np.max(abs(np.asarray(expected)-row[["A_pares", "A_balanceada"]].to_numpy(float)))))
    assert max_area_error < 1e-12 and max_H_error < 1e-12 and max_A_error < 1e-12
    # Validate the prominent report conclusions against all saved controls.
    controls = [c + "__hellinger_principal" for c in ["inicio_2048", "centro_2048", "fim_2048", "panorama_completo"]]
    anchor_del = deletion[deletion.metodo == MAIN].set_index("aquisicao_retirada")
    for control in controls:
        control_del = deletion[deletion.metodo == control].set_index("aquisicao_retirada")
        for column in ["A_pares", "A_balanceada"]:
            assert summary.loc[MAIN, column] < summary.loc[control, column]
            assert (anchor_del[column] < control_del[column]).all()
    for column in ["A_pares", "A_balanceada"]:
        a = by_pattern[by_pattern.metodo == MAIN].set_index("padrao")[column]
        b = by_pattern[by_pattern.metodo == FULL].set_index("padrao")[column]
        assert (a < b).all()
        assert summary.loc["ancora_2048__hellinger_conservador", column] > summary.loc[MAIN, column]
        assert summary.loc["ancora_2048__hellinger_conservador", column] < summary.loc["fim_2048__hellinger_conservador", column]
        assert summary.loc["ancora_2048__hellinger_conservador", column] < summary.loc["panorama_completo__hellinger_conservador", column]
    links = 0
    for name in ["README.md", "RELATORIO.md", "PROTOCOLO.md"]:
        source = (ROOT / name).read_text()
        for target in re.findall(r"\]\(([^)]+)\)", source):
            if target.startswith(("https://", "http://", "#")):
                continue
            if target == "results/verificacao_final.json":
                continue  # Written after all checks pass.
            assert (ROOT / target).is_file(), (name, target)
            links += 1
    page = (ROOT / "RESUMO_EXECUTIVO.html").read_text()
    images = re.findall(r"src='data:image/png;base64,([^']+)'", page)
    assert len(images) == 6 and "<html lang='pt-BR'>" in page
    image_hashes = {hashlib.sha256(base64.b64decode(value)).hexdigest() for value in images}
    presentation = json.loads((OUT / "apresentacao.json").read_text())
    assert presentation["codigo_sha256"] == sha(ROOT / "scripts/render_results.py")
    for name, expected in presentation["tabelas_sha256"].items():
        assert sha(OUT / name) == expected
    for figure in presentation["figuras"]:
        path = ROOT / figure["arquivo"]
        assert figure["sha256"] == sha(path) and sha(path) in image_hashes
        with Image.open(path) as image:
            image.verify()
        with Image.open(path) as image:
            assert min(image.size) >= 1200 and list(image.size) == figure["tamanho_px"]
    report = (ROOT / "RELATORIO.md").read_text()
    assert all(method in report for method in summary.index)
    assert "não é um pré-registro cego" in report
    assert "não sustentou substituir" in report and "Não é acurácia" in report
    result = dict(status="PASS", data_utc=datetime.now(timezone.utc).isoformat(),
        arquivos_artefatos=len(manifest), arquivos_fontes=len(frozen["sha256"]),
        pares_verificados=len(pairs), retiradas_verificadas=len(deletion),
        maior_erro_area=max_area_error, maior_erro_H_por_Bhattacharyya=max_H_error,
        maior_erro_A_e_retiradas=max_A_error, links_locais=links, figuras_incorporadas=len(images),
        conclusoes_conferidas=["desvantagem em ambos os pesos nos quatro controles principais",
            "desvantagem persiste nas 25 retiradas pareadas", "desvantagem perante completo nos quatro padrões",
            "faixa conservadora melhora âncoras mas não supera fim/completo"],
        manifesto_sha256=sha(OUT / "artefatos.json"), codigo_sha256=sha(__file__),
        observacao="Verificação de integridade e consistência; não é validação externa do método.")
    (OUT / "verificacao_final.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
