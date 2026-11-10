"""Consolida a proveniência local dos sinais e campos calculados na etapa 02."""
import json
from pathlib import Path
import platform
import sys

sys.dont_write_bytecode = True
import numpy as np
import scipy
from runtime import OUTPUT_ROOT, SCRIPTS_ROOT, reference
from recalcular_wavelet import ROOT, INTER, sha, write_json, source_record


def main():
    summary = json.loads((ROOT / "resumo_tcc.json").read_text())
    assert summary["execucao_pares_wavelet_completa"]
    records = json.loads((INTER / "aquisicoes.json").read_text())
    pairs = json.loads((INTER / "pares.json").read_text())
    assert len(records) == 26 and len(pairs) == 325
    assert len({(p["run_x"], p["run_y"]) for p in pairs}) == 325
    assert sum(p["tipo"] == "intra" for p in pairs) == 87
    assert all(p["escalas_validas"] == 200 for p in pairs)
    for row in records:
        path = INTER / "individuais" / f"{row['ensaio']}.npz"
        metadata = json.loads(path.with_suffix(".json").read_text())
        assert sha(path) == metadata["npz_sha256"]
        with np.load(path) as values:
            assert len(values["prepared"]) == row["comprimento"]
            assert values["wps_zero"].shape == (200,)
    manifest = dict(status="PASS", complete=True, stage="02_wavelet",
        inputs=source_record(OUTPUT_ROOT / "01_panoramas/manifesto.json"),
        scripts=[source_record(SCRIPTS_ROOT / name) for name in
                 ("wavelet_core.py", "recalcular_wavelet.py", "finalizar_wavelet.py")],
        runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, fft_threads=1),
        counts=dict(individuals=26, pairs=325, intra=87, inter=238, wavelengths=200),
        method="Sinais completos para perfis individuais; prefixos comuns para campos de pares; todos os atrasos.",
        statistics_stage="04_analise_simplificada",
        files={reference(p): sha(p) for p in sorted(ROOT.rglob("*"))
               if p.is_file() and p.name != "manifesto.json"
               and not {"historico", "__pycache__"}.intersection(p.relative_to(ROOT).parts)})
    write_json(ROOT / "manifesto.json", manifest)
    print(json.dumps({"status": "PASS", "files": len(manifest["files"]),
                      "manifesto": reference(ROOT / "manifesto.json")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
