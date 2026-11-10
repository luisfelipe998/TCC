#!/usr/bin/env python3
"""Regenerate the v14.1 correction, comparisons and adversarial checks."""
import argparse
from datetime import datetime, timezone
import json
import os
import subprocess
import sys
import time

from v14_paths import ROOT, SCRIPT_DIR, CACHE, artifact_path


STAGES = (
    ("reconstruct_inputs", "Reconstruct original panoramas as a control, then preserve source-strip order."),
    ("panorama_void_alternatives", "Infer liquid references and reconstruct conditional visible gas area."),
    ("phase_mask_flip_analysis", "Save gas masks and calculate spatial gas-silhouette flips."),
    ("raw_flip_analysis", "Reproduce legacy optical flips and supplemental native-frame checks."),
    ("spatial_flip_sensitivity", "Perturb spatial optical thresholds and radial crops."),
    ("inventory_and_legacy_audit", "Record input hashes and independently audit legacy values."),
    ("pattern_identification", "Compare fixed morphology classifiers and create the contact sheet."),
    ("silhouette_pattern_identification", "Compare silhouette classifiers with training-fold references."),
    ("adversarial_review", "Challenge masks with phase anchors, reference choices and crops."),
    ("deep_method_checks", "Independently verify masks and equations; compare v14 and v14.1 at D53."),
    ("lumen_normalization_checks", "Verify the unchanged D53 geometry on current masks."),
    ("image_void_alternatives", "Reproduce supplemental native-frame alternatives."),
    ("final_verification", "Independently verify saved arithmetic, input hashes and classifier scores."),
    ("summary_figures", "Plot void versus spatial flip-rate and crop sensitivity from saved results."),
    ("compare_versions", "Compare every acquisition, all threshold settings and all classifiers."),
    ("write_summary", "Generate the executive summary and folder navigation."),
    ("seam_adversarial", "Challenge source-column order, preserved methods, thresholds and change claims."),
    ("post_regeneration_adversarial", "Independently challenge the completed analysis, plots and findings."),
    ("verify_layout", "Verify current exclusions, retained measurements and report navigation."),
)
EXPLORATORY_STAGE = (
    "image_void_alternatives",
    "Regenerate earlier raw-frame alternatives in intermediate/exploratory.",
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--include-exploratory", action="store_true",
                        help="compatibility flag; native-frame alternatives are always regenerated")
    parser.add_argument("--list", action="store_true",
                        help="show the ordered stages without running them")
    args = parser.parse_args()
    stages = list(STAGES)
    if args.list:
        for number, (name, description) in enumerate(stages, 1):
            print(f"{number:02d}. {name}.py — {description}")
        return 0

    log_dir = artifact_path("pipeline_logs")
    log_dir.mkdir(parents=True, exist_ok=True)
    matplotlib_cache = CACHE / "matplotlib"
    matplotlib_cache.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["MPLCONFIGDIR"] = str(matplotlib_cache)
    manifest = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "python_executable": sys.executable,
        "include_exploratory": args.include_exploratory,
        "status": "running",
        "stages": [],
    }
    manifest_path = artifact_path("pipeline_run.json")

    def save_manifest():
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    save_manifest()
    for number, (name, _) in enumerate(stages, 1):
        command = [sys.executable, "-B", str(SCRIPT_DIR / (name + ".py"))]
        log = log_dir / (name + ".log")
        print(f"[{number}/{len(stages)}] {name}.py", flush=True)
        started = time.perf_counter()
        with log.open("w") as output:
            completed = subprocess.run(command, cwd=ROOT, env=env, stdout=output,
                                       stderr=subprocess.STDOUT, check=False)
        manifest["stages"].append({
            "script": f"scripts/{name}.py",
            "log": log.relative_to(ROOT).as_posix(),
            "elapsed_seconds": time.perf_counter() - started,
            "returncode": completed.returncode,
        })
        if completed.returncode:
            manifest["status"] = "failed"
            manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
            save_manifest()
            print(f"Stage failed; inspect {log}", file=sys.stderr)
            return completed.returncode
        save_manifest()
    manifest["status"] = "passed"
    manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
    save_manifest()
    print(f"Completed. Start with {ROOT / 'README.md'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
