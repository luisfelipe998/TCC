"""Verifica preservação, integridade da entrega e guarda residual sem escrever fontes."""
from pathlib import Path
import ast
import hashlib
import json
import re
import subprocess
import sys
from datetime import datetime, timezone

sys.dont_write_bytecode = True
BASE = Path(__file__).resolve().parents[1]
PROJECT = BASE.parents[1]
OUT = BASE / "resultados"


def sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def main():
    initial = json.loads((OUT / "fontes_antes.json").read_text())
    changed = []
    for r in initial:
        p = PROJECT / r["path"]
        actual = sha(p) if p.is_file() else None
        if actual != r["sha256"]:
            changed.append({"path": r["path"], "sha256_antes": r["sha256"], "sha256_agora": actual})
    initial_paths = {r["path"] for r in initial}
    new_sources = []
    for folder in ["TCC", "src2/estatistica8", "src2/v14.1_new_voidfraction", "src2/v16_spatial2_18_71"]:
        new_sources += [str(p.relative_to(PROJECT)) for p in (PROJECT / folder).rglob("*")
                        if p.is_file() and str(p.relative_to(PROJECT)) not in initial_paths]
    before = json.loads((OUT / "estado_outros_arquivos_antes.json").read_text())
    stat_changes = []
    for row in before:
        p = PROJECT / row["path"]
        st = p.stat() if p.is_file() else None
        if st is None or (st.st_size, st.st_mtime_ns) != (row["size"], row["mtime_ns"]):
            stat_changes.append(row["path"])
    known = {r["path"] for r in before}
    new_other = []
    for folder in ["src", "src2", "TCC", "TCC_1"]:
        for p in (PROJECT / folder).rglob("*"):
            if p.is_file() and BASE not in p.parents and str(p.relative_to(PROJECT)) not in known:
                new_other.append(str(p.relative_to(PROJECT)))
    broken_links = []
    for p in BASE.rglob("*.md"):
        for target in re.findall(r"\]\(([^)]+)\)", p.read_text()):
            if re.match(r"[a-z]+://", target) or target.startswith("#"):
                continue
            target = target.split("#", 1)[0].strip("<>")
            dest = p.parent / target
            # This file is created by the verification itself.
            if dest.resolve() == (OUT / "validacao_entrega.json").resolve():
                continue
            if not dest.exists():
                broken_links.append({"file": str(p.relative_to(BASE)), "target": target})
    scripts = sorted((BASE / "scripts").glob("*.py"))
    for p in scripts:
        ast.parse(p.read_text(), filename=str(p))
    json_files = sorted(BASE.rglob("*.json"))
    for p in json_files:
        json.loads(p.read_text())
    # Exercise the actual corrected function in isolation, without importing pipelines.
    import numpy as np
    p = BASE / "scripts/exploracao_espectral.py"
    node = next(n for n in ast.parse(p.read_text()).body if isinstance(n, ast.FunctionDef) and n.name == "residual")
    namespace = {"np": np}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(p), "exec"), namespace)
    lam = np.geomspace(50, 512, 200)
    guard = {}
    for name, profile in [("constante", np.ones((1, 200))), ("lei_potencia", lam[None, :] ** 1.25)]:
        try:
            namespace["residual"](profile, lam)
            guard[name] = False
        except ValueError:
            guard[name] = True
    real_profiles = 0
    with np.load(OUT / "exploracao/espectros_janelas.npz") as z:
        for key in z.files:
            if key in ["runs", "labels"]:
                continue
            cap = int(key.split("_", 1)[0])
            out = namespace["residual"](z[key], np.geomspace(50, cap, 200))
            assert np.isfinite(out).all()
            real_profiles += len(out)
    git_current = subprocess.run(["git", "status", "--short"], cwd=PROJECT, check=True,
                                 capture_output=True, text=True).stdout
    git_before = (BASE / "logs/git_status_antes.txt").read_text()
    (BASE / "logs/git_status_depois.txt").write_text(git_current)
    git_added_lines = sorted(set(git_current.splitlines()) - set(git_before.splitlines()))
    git_removed_lines = sorted(set(git_before.splitlines()) - set(git_current.splitlines()))
    analysis_changes = [v for v in changed if v["path"].startswith("src2/")]
    text_snapshot = OUT / "comparacao_tcc/snapshot_texto_final.json"
    final_text_checked = False
    if text_snapshot.exists():
        text_rows = json.loads(text_snapshot.read_text())["fontes"]
        final_text_checked = all((PROJECT / r["path"]).is_file()
                                 and sha(PROJECT / r["path"]) == r["sha256"] for r in text_rows)
        final_text_checked &= (BASE / "findings/reconciliacao_texto_atual.md").is_file()
    unexpected_stat_changes = [p for p in stat_changes if not p.startswith("TCC/")]
    computational_ok = not broken_links and all(guard.values())
    source_ok = not (analysis_changes or new_sources or unexpected_stat_changes or new_other)
    if computational_ok and source_ok and not changed and not stat_changes:
        status = "PASS"
    elif computational_ok and source_ok and final_text_checked:
        status = "PASS_WITH_DOCUMENTED_TCC_CHANGES"
    else:
        status = "REVIEW_REQUIRED"
    result = {
        "status": status,
        "verificacoes_computacionais": "PASS" if computational_ok else "FAIL",
        "preservacao_tres_pastas_analise": "PASS" if not analysis_changes else "FAIL",
        "arquivos_tres_pastas_analise": sum(r["path"].startswith("src2/") for r in initial),
        "tcc_snapshot_inicial_identico": not any(r["path"].startswith("TCC/") for r in changed),
        "tcc_versao_final_reconferida_e_hash_estavel": bool(final_text_checked),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "fontes_auditadas_sha256": len(initial), "fontes_modificadas": changed,
        "novos_arquivos_nas_fontes": new_sources,
        "arquivos_externos_conferidos_por_tamanho_mtime": len(before),
        "arquivos_externos_modificados": stat_changes, "novos_arquivos_nos_diretorios_monitorados": new_other,
        "scripts_sintaticamente_validos": len(scripts), "jsons_validos": len(json_files),
        "links_quebrados": broken_links, "guarda_residual": guard,
        "perfis_reais_aceitos_pela_guarda": real_profiles,
        "git_status_identico": git_current == git_before,
        "git_linhas_adicionadas": git_added_lines, "git_linhas_removidas": git_removed_lines,
        "nota_git": "Estado observado, não atribuição de autoria. A trilha desta tarefa escreve apenas em v17. "
                    "Alterações concorrentes no TCC foram reconferidas e estão documentadas; "
                    "adições em assets não integram a auditoria nem a amostra.",
        "alcance": "Verificação computacional e preservação; não validação física, causal ou preditiva.",
    }
    (OUT / "validacao_entrega.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    manifest = []
    for p in sorted(BASE.rglob("*")):
        if not p.is_file() or p == OUT / "artefatos_sha256.json":
            continue
        # Logs redirected to while this process runs are not stable artifacts.
        if p.name == "verificar_entrega.log":
            continue
        manifest.append({"path": str(p.relative_to(BASE)), "bytes": p.stat().st_size, "sha256": sha(p)})
    (OUT / "artefatos_sha256.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    print(f"Artefatos registrados: {len(manifest)}")
    assert result["status"] in ["PASS", "PASS_WITH_DOCUMENTED_TCC_CHANGES"]


if __name__ == "__main__":
    main()
