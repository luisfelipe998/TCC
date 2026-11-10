"""File dependency audit hook used by verificar_portabilidade.py.

The hook is a reproducibility check for this cooperative Python pipeline, not
an operating-system security boundary. Python children inherit its site hook.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import sys


def install() -> None:
    config_path = os.environ.get("TCC_V18_PORTABILITY_CONFIG")
    if not config_path:
        return
    config = json.loads(Path(config_path).read_text())
    roots = {key: [Path(path).resolve() for path in paths]
             for key, paths in config["roots"].items()}
    log = Path(config["log_directory"])
    log.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(log / f"audit-{os.getpid()}.jsonl",
                         os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600)
    observed = set()

    def classify(path: Path) -> str | None:
        for category in ("package", "data", "output", "cache", "support", "system"):
            if any(path.is_relative_to(root) for root in roots[category]):
                return category
        return None

    def check(event: str, value, writing: bool = False, importing: bool = False):
        if not isinstance(value, (str, bytes, os.PathLike)):
            return
        value = os.fsdecode(value)
        if value.startswith("<") and value.endswith(">"):
            return
        path = Path(value).resolve()
        category = classify(path)
        data_code = category == "data" and (
            importing or path.suffix.lower() in {".py", ".pyc", ".pyo", ".so", ".dylib"})
        allowed = category is not None and not data_code
        if writing:
            allowed = category in {"output", "cache"} or str(path) == "/dev/null"
        key = (event, str(path), writing, importing, allowed)
        if key not in observed:
            observed.add(key)
            entry = {"event": event, "path": str(path), "category": category,
                     "writing": writing, "importing": importing, "allowed": allowed}
            os.write(descriptor, (json.dumps(entry) + "\n").encode())
        if not allowed:
            raise PermissionError(f"PORTABILITY_BLOCKED: {event}: {path}")

    def hook(event, args):
        if event == "open":
            mode = args[1] or ""
            flags = args[2] or 0
            writing = any(char in mode for char in "wax+") or bool(
                flags & (os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_TRUNC | os.O_APPEND))
            check(event, args[0], writing=writing)
        elif event in {"os.listdir", "os.scandir"}:
            check(event, args[0] if args else os.getcwd())
        elif event == "import" and len(args) > 1 and args[1]:
            check(event, args[1], importing=True)
        elif event == "compile" and len(args) > 1:
            check(event, args[1], importing=True)
        elif event in {"os.mkdir", "os.remove", "os.rmdir"}:
            check(event, args[0], writing=True)
        elif event in {"os.rename", "os.replace"}:
            check(event, args[0], writing=True)
            check(event, args[1], writing=True)
        elif event == "subprocess.Popen":
            # Native tools would not inherit Python's audit hook. The v18
            # pipeline needs Python child processes only.
            executable = Path(os.fsdecode(args[0])).resolve()
            if executable != Path(sys.executable).resolve():
                os.write(descriptor, (json.dumps({
                    "event": event, "path": str(executable),
                    "category": classify(executable), "writing": False,
                    "importing": False, "allowed": False,
                    "reason": "non-Python subprocess does not inherit the hook",
                }) + "\n").encode())
                raise PermissionError(f"PORTABILITY_BLOCKED: non-Python child: {executable}")

    sys.addaudithook(hook)
    os.write(descriptor, (json.dumps({"event": "guard_installed",
                                     "pid": os.getpid(), "allowed": True}) + "\n").encode())


if __name__ == "__main__":
    install()
