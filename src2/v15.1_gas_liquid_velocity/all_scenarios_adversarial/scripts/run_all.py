#!/usr/bin/env python3
"""Reproduce every executed branch, preserving all previous v15.1 artifacts."""
from pathlib import Path
import subprocess,sys

SCRIPTS=Path(__file__).resolve().parent
ORDER=(
    'block_ncc.py','selected_wall.py','split_background_ncc.py',
    'phase_local_experiment.py','phase_local_challenge.py','phase_local_review.py',
    'spectral_transport.py','spectral_finish.py',
    'adversarial_identifiability.py','adversarial_new_methods.py',
    'independent_block_verify.py','gate_sensitivity.py','final_figures.py',
    '../../scripts/publish_latest_summary.py','final_verify.py',
)

if __name__=='__main__':
    for name in ORDER:
        print('Running '+name,flush=True)
        subprocess.run([sys.executable,'-B',str(SCRIPTS/name)],check=True,cwd=SCRIPTS.parents[3])
