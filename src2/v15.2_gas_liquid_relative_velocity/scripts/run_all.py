"""Reproduce the v15.2 investigation without changing earlier result versions."""
from pathlib import Path
import subprocess,sys
HERE=Path(__file__).resolve().parent
STAGES=[
    'object_motion.py','meniscus_motion.py','robust_meniscus.py',
    'object_comparisons.py','paired_object_motion.py','physical_bounds.py','conditional_bubble_bounds.py',
    'residual_motion.py','residual_challenge.py','residual_closed_contour.py',
    'residual_visual_check.py','residual_real_check.py','residual_report.py',
    'independent_object_verify.py','independent_comparison_verify.py','write_object_audit_report.py',
    'build_figures.py','build_report.py','final_verify.py']
for stage in STAGES:
    print('Running',stage,flush=True)
    subprocess.run([sys.executable,'-B',str(HERE/stage)],check=True)
