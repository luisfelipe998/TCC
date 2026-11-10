#!/usr/bin/env python3
"""Rebuild v15 outputs locally, stopping on any failed analysis or verification."""
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import time

BASE = Path(__file__).resolve().parents[1]


def main():
    logs = BASE / 'intermediate/logs'
    logs.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    env['MPLCONFIGDIR'] = str(BASE / 'intermediate/matplotlib')
    stages = ['scripts/audit_legacy.py', 'scripts/massflow_audit.py',
              'scripts/model_comparison.py', 'scripts/optical_investigation.py',
              'literature/source_numeric_checks.py', 'scripts/write_summary.py',
              'scripts/independent_model_review.py', 'scripts/final_verification.py']
    record = dict(started_utc=datetime.now(timezone.utc).isoformat(),
                  interpreter=sys.executable, stages=[], status='running')
    target = logs / 'pipeline_run.json'
    for stage in stages:
        command = [sys.executable, '-B', str(BASE / stage)]
        logfile = logs / (Path(stage).stem + '.log')
        print(f'Running {stage}', flush=True)
        started = time.monotonic()
        with logfile.open('w') as stream:
            completed = subprocess.run(command, cwd=BASE.parents[1], env=env,
                                       stdout=stream, stderr=subprocess.STDOUT)
        record['stages'].append(dict(script=stage, command=command,
                                     elapsed_s=time.monotonic()-started,
                                     returncode=completed.returncode,
                                     log=str(logfile.relative_to(BASE))))
        if completed.returncode:
            record['status'] = 'failed'
            target.write_text(json.dumps(record, indent=2) + '\n')
            print(logfile.read_text(), file=sys.stderr)
            raise SystemExit(completed.returncode)
        target.write_text(json.dumps(record, indent=2) + '\n')
    record['status'] = 'passed'
    record['finished_utc'] = datetime.now(timezone.utc).isoformat()
    target.write_text(json.dumps(record, indent=2) + '\n')
    print('All eight analysis/verification stages passed.', flush=True)


if __name__ == '__main__':
    main()
