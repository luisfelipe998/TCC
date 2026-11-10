#!/usr/bin/env python3
"""Reproduce the canonical current optical analysis on the v18 panoramas."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parents[1]
CANONICAL = ROOT.parent / "v18_montagem_18_71"


def register_current():
    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    output = CANONICAL / "03_optica"
    manifest = output / "results/validation/artifact_manifest.json"
    files = json.loads(manifest.read_text())["files"]
    for rel, expected in files.items():
        assert sha(output / rel) == expected, rel
    record = dict(status="current_canonical_optical_execution",
        input_panorama_manifest=str((CANONICAL / "01_panoramas/manifesto.json").relative_to(PROJECT)),
        input_panorama_manifest_sha256=sha(CANONICAL / "01_panoramas/manifesto.json"),
        current_output_directory=str(output.relative_to(PROJECT)),
        current_output_manifest=str(manifest.relative_to(PROJECT)),
        current_output_manifest_sha256=sha(manifest),
        shared_numerical_model=str((ROOT / "scripts/panorama_void_alternatives.py").relative_to(PROJECT)),
        shared_numerical_model_sha256=sha(ROOT / "scripts/panorama_void_alternatives.py"),
        roi=[18,72], included_rows="18 through 71 inclusive", height=54, diameter_px=54,
        scope="Current panorama optical analysis used by the TCC; same canonical execution, not independent replication",
        previous_results="intermediate/history/before_roi_unification_20261030",
        previous_panorama_alignment="The preserved reconstructed_inputs used alignment profiles 20:60; current inputs are the v18 optical panoramas aligned on 18:72",
        wavelet_recomputed=False, classifiers_recomputed=False)
    (ROOT / "execucao_atual.json").write_text(json.dumps(record, indent=2, ensure_ascii=False)+"\n")
    return record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--register-only", action="store_true", help="Verify and register existing canonical outputs without recomputing")
    args=parser.parse_args()
    if not args.register_only:
        subprocess.run([sys.executable,"-B",str(CANONICAL/"scripts/process_optical.py")],cwd=PROJECT,check=True)
    print(json.dumps(register_current(),indent=2,ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
