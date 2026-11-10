"""Executa os cálculos da faixa 18–71; relatórios são gerados após a revisão."""
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
ENV=dict(os.environ,MPLCONFIGDIR=os.environ.get('MPLCONFIGDIR','/tmp/v16-18-71-mpl'))

def stage(name,args):
    start=time.perf_counter()
    print(f'Início: {name}',flush=True)
    with (ROOT/'results'/f'{name}.log').open('w') as log:
        subprocess.run([sys.executable,str(ROOT/'scripts'/args[0]),*args[1:]],
                       env=ENV,stdout=log,stderr=subprocess.STDOUT,check=True)
    row=dict(stage=name,elapsed_seconds=time.perf_counter()-start)
    print(f'Concluído: {name}, {row["elapsed_seconds"]:.1f} s',flush=True)
    return row

def primary():
    return [stage(name,args) for name,args in [
        ('run_analysis',['run_analysis.py','--stage','all']),
        ('verify_legacy_baseline',['verify_legacy_baseline.py']),
        ('additional_diagnostics',['additional_diagnostics.py']),
        ('precision_and_windows',['precision_and_windows.py'])]]

def main():
    start=time.perf_counter()
    validation=stage('validate_wavelets',['validate_wavelets.py'])
    with ThreadPoolExecutor(max_workers=2) as pool:
        main_batch=pool.submit(primary)
        statistics=pool.submit(stage,'recompute_unmasked',['recompute_unmasked.py'])
        rows=[validation,*main_batch.result(),statistics.result()]
    result=dict(status='PASS',extraction_rows_inclusive=[18,71],
                total_elapsed_seconds=time.perf_counter()-start,stages=rows,
                note='Dois lotes independentes em paralelo; relatórios e verificação final executados depois.')
    (ROOT/'results/full_run_metadata.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':main()
