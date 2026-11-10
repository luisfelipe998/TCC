#!/usr/bin/env python3
"""Reproduce v15.1 locally, stopping on a failed stage."""
from pathlib import Path
from datetime import datetime, timezone
import json
import os
import subprocess
import sys
import time

BASE=Path(__file__).resolve().parents[1]
STAGES=['prepare_passage_results.py','calibration_audit.py','mass_balance_stress.py',
        'annular_motion.py','annular_report.py','independent_verification.py',
        'independent_regional_verification.py','paper_figures.py','build_summary.py',
        'independent_paper_verification.py','final_verification.py']


def run():
    logs=BASE/'intermediate/logs'
    logs.mkdir(exist_ok=True)
    record=dict(started_utc=datetime.now(timezone.utc).isoformat(),status='running',stages=[])
    env=os.environ.copy()
    env['PYTHONDONTWRITEBYTECODE']='1'
    for stage in STAGES:
        script=BASE/'scripts'/stage
        assert script.exists(),stage
        start=time.monotonic()
        print('Running '+stage,flush=True)
        with (logs/(script.stem+'.log')).open('w') as stream:
            result=subprocess.run([sys.executable,'-B',str(script)],cwd=BASE.parent.parent,
                                  stdout=stream,stderr=subprocess.STDOUT,env=env)
        record['stages'].append(dict(script=stage,returncode=result.returncode,elapsed_s=time.monotonic()-start))
        if result.returncode:
            record['status']='failed'
        (logs/'pipeline_run.json').write_text(json.dumps(record,indent=2)+chr(10))
        if result.returncode:
            raise SystemExit(f'{stage} failed; see {logs/(script.stem+".log")}')
    record.update(status='passed',finished_utc=datetime.now(timezone.utc).isoformat())
    (logs/'pipeline_run.json').write_text(json.dumps(record,indent=2)+chr(10))
    print('All v15.1 stages passed.')


if __name__=='__main__':
    run()
