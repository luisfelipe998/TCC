#!/usr/bin/env python3
"""Run v18 from a minimal copy while denying undeclared file dependencies.

Only scripts/config/04 scripts and executar.py are copied. The output starts
empty, the original project is inaccessible, and data_processed permits data
reads but no Python source or extension imports. No historical output is copied.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import os
from pathlib import Path
import shutil
import site
import subprocess
import sys
import sysconfig
import tempfile

from runtime import DATA_ROOT, PACKAGE_ROOT


def system_roots() -> list[str]:
    candidates = [sys.prefix, sys.base_prefix, sys.exec_prefix, sys.base_exec_prefix,
                  *sysconfig.get_paths().values(), *site.getsitepackages(),
                  site.getusersitepackages(), "/usr", "/System", "/Library/Fonts",
                  "/Library/Apple", "/etc", "/private/etc", "/dev", "/proc",
                  str(Path.home() / "Library/Fonts")]
    return sorted({str(Path(path).resolve()) for path in candidates if path})


def copy_minimal(destination: Path) -> None:
    for directory in ("scripts", "config", "04_analise_simplificada/scripts"):
        source = PACKAGE_ROOT / directory
        if not source.is_dir():
            raise FileNotFoundError(source)
        shutil.copytree(source, destination / directory,
                        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    shutil.copyfile(PACKAGE_ROOT / "executar.py", destination / "executar.py")


def presentation_coverage(output: Path) -> dict:
    """The four methodological diagrams and four raw-signal examples are required."""
    names = ("fluxograma", "perfis_wavelet", "exemplo_mascara", "medidas_opticas")
    required = [f"05_figuras/metodologia/{name}.{extension}"
                for name in names for extension in ("pdf", "png")]
    required += [f"05_figuras/metodologia/{name}_origem.json" for name in names]
    required += [f"05_figuras/metodologia/{name}_dados.npz" for name in names[1:]]
    required += ["05_figuras/metodologia/manifesto.json", "05_figuras/panoramas/origens.json"]
    required += [f"05_figuras/panoramas/{name}.png" for name in (
        "plug_m05q01_001", "slug_m04q01_001", "churn_m04q20_001", "wavy_annular_m06q40_001")]
    missing = [name for name in required if not (output / name).is_file()]
    kinds = {}
    manifest = output / "05_figuras/metodologia/manifesto.json"
    if manifest.is_file():
        kinds = {entry["name"]: entry["kind"]
                 for entry in json.loads(manifest.read_text())["figures"]}
    expected = dict(zip(names, ("diagram", "synthetic_wavelet_example",
                               "experimental_crop", "synthetic_mask_example")))
    return {"complete": not missing and kinds == expected,
            "required_artifacts": len(required), "missing": missing,
            "methodology_sources": kinds}


def prepare_font_cache(cache: Path) -> None:
    """Initialize OS font discovery without giving native tools to the pipeline.

Matplotlib may run fc-list/system_profiler while creating an empty font cache.
Only this library initialization runs before the hook; no package code or
experimental input is loaded. The scientific subprocess receives the cache.
"""
    env = dict(os.environ, PYTHONPATH="", PYTHONDONTWRITEBYTECODE="1",
               MPLCONFIGDIR=str(cache / "matplotlib"),
               XDG_CACHE_HOME=str(cache), TMPDIR=str(cache))
    env.pop("TCC_V18_PORTABILITY_CONFIG", None)
    subprocess.run([sys.executable, "-B", "-c", "import matplotlib.font_manager"],
                   cwd=cache, env=env, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=DATA_ROOT)
    parser.add_argument("--output", type=Path,
                        help="Empty candidate output; defaults to a temporary directory.")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--probe-only", action="store_true",
                        help="Test the hook without running the scientific pipeline.")
    args = parser.parse_args()
    data = args.data_root.resolve()
    if not data.is_dir():
        parser.error(f"Data directory does not exist: {data}")
    work = Path(tempfile.mkdtemp(prefix="tcc-v18-portability-")).resolve()
    package = work / "package"
    copy_minimal(package)
    output = args.output.resolve() if args.output else work / "output"
    if output.exists() and any(output.iterdir()):
        parser.error(f"Output must be empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    support, cache = work / "support", work / "cache"
    support.mkdir(); cache.mkdir()
    prepare_font_cache(cache)
    logs = work / "audit"
    logs.mkdir()
    configuration = {
        "roots": {"package": [str(package)], "data": [str(data)],
                  "output": [str(output)], "cache": [str(cache)],
                  "support": [str(support)], "system": system_roots()},
        "log_directory": str(logs),
    }
    config = support / "guard.json"
    config.write_text(json.dumps(configuration, indent=2) + "\n")
    # Import this local guard in every Python child, including multiprocessing.
    (support / "sitecustomize.py").write_text(
        "import os, sys\n"
        f"sys.path.insert(0, {str(package / 'scripts')!r})\n"
        "try:\n"
        "    from portability_guard import install\n"
        "    install()\n"
        "except BaseException as error:\n"
        "    print('Portability guard could not start: ' + repr(error), file=sys.stderr)\n"
        "    os._exit(97)\n")
    env = dict(os.environ, TCC_V18_PORTABILITY_CONFIG=str(config),
               PYTHONPATH=str(support), PYTHONDONTWRITEBYTECODE="1",
               TCC_DATA_ROOT=str(data), TCC_V18_OUTPUT=str(output),
               MPLCONFIGDIR=str(cache / "matplotlib"),
               XDG_CACHE_HOME=str(cache), TMPDIR=str(cache))
    if args.probe_only:
        # A denied dependency remains recorded even if application code catches
        # the exception. Here exactly these two denials are intentional probes.
        code = (
            "from pathlib import Path\n"
            f"blocked = {[str(PACKAGE_ROOT / 'scripts/runtime.py'), str(data / 'pipeline.py')]!r}\n"
            "for path in blocked:\n"
            "    try:\n        Path(path).read_bytes()\n"
            "    except PermissionError:\n        pass\n"
            "    else:\n        raise AssertionError(path)\n"
            f"Path({str(output / 'probe.txt')!r}).write_text('portable output')\n"
        )
        command = [sys.executable, "-B", "-c", code]
    else:
        command = [sys.executable, "-B", str(package / "executar.py"),
                   "--data-root", str(data), "--output", str(output),
                   "--workers", str(args.workers)]
    result = subprocess.run(command, cwd=package, env=env, check=False)
    events = [json.loads(line) for path in sorted(logs.glob("audit-*.jsonl"))
              for line in path.read_text().splitlines()]
    installed = [event for event in events if event["event"] == "guard_installed"]
    blocked = [event for event in events if not event["allowed"]]
    valid_denials = len(blocked) == 2 if args.probe_only else not blocked
    coverage = None if args.probe_only else presentation_coverage(output)
    passed = (result.returncode == 0 and bool(installed) and valid_denials
              and (args.probe_only or coverage["complete"]))
    validation = output / "validacao"
    validation.mkdir(exist_ok=True)
    report = {
        "status": "PASS" if passed else "FAIL", "probe_only": args.probe_only,
        "returncode": result.returncode, "minimal_package": str(package),
        "data_root": str(data), "output": str(output),
        "original_package_not_allowed": str(PACKAGE_ROOT),
        "copied": ["scripts", "config", "04_analise_simplificada/scripts", "executar.py"],
        "processes_with_guard": len(installed), "blocked_events": blocked,
        "file_events_by_category": dict(Counter(event.get("category") for event in events
                                                 if "category" in event)),
        "audit_logs": str(logs),
        "os_font_cache_prepared_before_hook": True,
        "methodology_presentation_coverage": coverage,
        "scope": "Python file audit for cooperative pipeline; not an OS security sandbox.",
    }
    (validation / "portabilidade.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
