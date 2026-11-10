"""Compare generated scientific outputs with the frozen, portable reference.

The reference is configuration, never a runtime read of a historical execution.
Floating values use atol=rtol=1e-10; discrete values and image pixels are exact.
Run after each stage with --stage panoramas|wavelet|optica|completo.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import sys

sys.dont_write_bytecode = True
import numpy as np
from PIL import Image
from runtime import CONFIG_ROOT, OUTPUT_ROOT, artifact_path, reference

STAGES = ("panoramas", "wavelet", "optica", "completo")
ATOL = RTOL = 1e-10


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def scalar(value):
    if value == "":
        return None
    if value in ("True", "False"):
        return value == "True"
    if re.fullmatch(r"[-+]?\d+", value):
        return int(value)
    try:
        return float(value)
    except ValueError:
        return value


def read_csv(path):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        return list(reader.fieldnames or []), [{k: scalar(v) for k, v in row.items()} for row in reader]


def project(value, selector):
    """Select explicitly documented scientific fields; no name-based exclusions."""
    if selector is None:
        return value
    if "fields" in selector:
        return {key: project(value[key], sub) for key, sub in selector["fields"].items()}
    if "rows" in selector:
        rows = [project(row, selector["rows"]) for row in value]
        if selector.get("sort_by"):
            rows.sort(key=lambda row: tuple(row[key] for key in selector["sort_by"]))
        return rows
    raise ValueError(f"Invalid reference selector: {selector}")


def pixel_record(path):
    with Image.open(path) as im:
        values = np.asarray(im)
    return dict(shape=list(values.shape), dtype=str(values.dtype),
                pixels_sha256=hashlib.sha256(values.tobytes()).hexdigest())


def compare_values(actual, expected, location, failures):
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or set(actual) != set(expected):
            failures.append(dict(location=location, problem="fields differ",
                expected=sorted(expected), actual=sorted(actual) if isinstance(actual, dict) else str(type(actual))))
            return
        for key in expected:
            compare_values(actual[key], expected[key], f"{location}/{key}", failures)
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(actual) != len(expected):
            failures.append(dict(location=location, problem="sequence length differs",
                expected=len(expected), actual=len(actual) if isinstance(actual, list) else str(type(actual))))
            return
        for index, (a, e) in enumerate(zip(actual, expected)):
            compare_values(a, e, f"{location}/{index}", failures)
    elif isinstance(expected, bool) or expected is None or isinstance(expected, str):
        if actual != expected or isinstance(expected, bool) and not isinstance(actual, bool):
            failures.append(dict(location=location, expected=expected, actual=actual))
    elif isinstance(expected, int):
        if isinstance(actual, bool) or not isinstance(actual, (int, float)) or actual != expected:
            failures.append(dict(location=location, problem="discrete value differs", expected=expected, actual=actual))
    elif isinstance(expected, float):
        if not isinstance(actual, (int, float)) or isinstance(actual, bool) or not math.isclose(actual, expected, rel_tol=RTOL, abs_tol=ATOL):
            failures.append(dict(location=location, problem="float differs beyond tolerance", expected=expected, actual=actual))
    else:
        raise TypeError(type(expected))


def array_difference(actual, expected):
    if actual.shape != expected.shape:
        return dict(problem="array shape differs", expected=list(expected.shape), actual=list(actual.shape))
    if actual.dtype.kind != expected.dtype.kind:
        return dict(problem="array dtype kind differs", expected=str(expected.dtype), actual=str(actual.dtype))
    if expected.dtype.kind in "fc":
        equal = np.isclose(actual, expected, atol=ATOL, rtol=RTOL, equal_nan=False)
    else:
        equal = actual == expected
    if np.all(equal):
        return None
    index = tuple(np.argwhere(~equal)[0]) if equal.ndim else ()
    def json_value(value):
        value = value.item()
        return [value.real, value.imag] if isinstance(value, complex) else value
    result = dict(problem="array values differ", first_index=[int(i) for i in index],
                  expected=json_value(expected[index]), actual=json_value(actual[index]))
    if expected.dtype.kind in "fciub":
        result["max_absolute_error"] = float(np.max(np.abs(actual.astype(complex)-expected.astype(complex))))
    return result


def report_value(value):
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    if isinstance(value, dict):
        return {key: report_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [report_value(item) for item in value]
    return value


def compare(stage="completo"):
    baseline_path = CONFIG_ROOT / "referencia_resultados.json"
    arrays_path = CONFIG_ROOT / "referencia_arrays.npz"
    baseline = json.loads(baseline_path.read_text())
    failures, checked = [], []
    if sha(arrays_path) != baseline["arrays_sha256"]:
        failures.append(dict(location=reference(arrays_path), problem="frozen reference array checksum differs"))
    target = STAGES.index(stage)
    with np.load(arrays_path, allow_pickle=False) as frozen:
        for spec in baseline["checks"]:
            if STAGES.index(spec["stage"]) > target:
                continue
            label = spec["path"]
            path = artifact_path(label)
            before = len(failures)
            try:
                if spec["kind"] == "pixels":
                    compare_values(pixel_record(path), spec["expected"], label, failures)
                elif spec["kind"] == "npz":
                    with np.load(path, allow_pickle=False) as values:
                        if set(values.files) != set(spec["arrays"]):
                            failures.append(dict(location=label, problem="array keys differ",
                                expected=sorted(spec["arrays"]), actual=sorted(values.files)))
                        for key, name in spec["arrays"].items():
                            if key not in values:
                                continue
                            difference = array_difference(values[key], frozen[name])
                            if difference:
                                failures.append(dict(location=f"{label}/{key}", **difference))
                elif spec["kind"] == "csv":
                    columns, rows = read_csv(path)
                    expected_columns = set(spec["columns"])
                    allowed_ignored = set(spec.get("ignored_columns", []))
                    if set(columns)-allowed_ignored != expected_columns:
                        failures.append(dict(location=label, problem="CSV scientific columns differ",
                            expected=sorted(expected_columns), actual=sorted(set(columns)-allowed_ignored)))
                    rows = [{key: row[key] for key in spec["columns"]} for row in rows]
                    if spec.get("sort_by"):
                        rows.sort(key=lambda row: tuple(row[key] for key in spec["sort_by"]))
                    compare_values(rows, spec["expected"], label, failures)
                    for index, (actual, expected) in enumerate(zip(rows, spec["expected"])):
                        for key, digits in spec.get("published_decimals", {}).items():
                            if expected[key] is not None and format(actual[key], f".{digits}f") != format(expected[key], f".{digits}f"):
                                failures.append(dict(location=f"{label}/{index}/{key}", problem="published rounded value differs",
                                    expected=format(expected[key], f".{digits}f"), actual=format(actual[key], f".{digits}f")))
                elif spec["kind"] == "json":
                    value = project(json.loads(path.read_text()), spec.get("projection"))
                    compare_values(value, spec["expected"], label, failures)
                else:
                    raise ValueError(spec["kind"])
            except (OSError, KeyError, ValueError, TypeError) as error:
                failures.append(dict(location=label, problem=type(error).__name__, detail=str(error)))
            checked.append(dict(path=label, kind=spec["kind"], passed=len(failures) == before))
    failures = report_value(failures)
    report = dict(status="PASS" if not failures else "FAIL", stage=stage,
        tolerance=dict(absolute=ATOL, relative=RTOL, pixels="exact", discrete="exact", selection="exact"),
        baseline=dict(file=reference(baseline_path), sha256=sha(baseline_path), arrays=reference(arrays_path), arrays_sha256=sha(arrays_path)),
        checked_files=len(checked), failed_checks=sum(not c["passed"] for c in checked),
        checks=checked, failures=failures,
        excluded=baseline["excluded_scope"])
    destination = OUTPUT_ROOT / "validacao" / f"comparacao_{stage}.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+"\n")
    print(json.dumps(dict(status=report["status"], stage=stage, checked_files=len(checked),
        failures=len(failures), report=reference(destination), first_failures=failures[:8]), ensure_ascii=False), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=STAGES, default="completo")
    result = compare(parser.parse_args().stage)
    raise SystemExit(0 if result["status"] == "PASS" else 1)
