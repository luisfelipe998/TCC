"""Recalcula os resumos que dependiam das máscaras, sem mudar as wavelets."""
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]

def main():
    tasks = [
        ('preprocessing_screen', ['preprocessing_experiment.py', '--scales', '64', '--tag', 'screen']),
        ('preprocessing_screen_evaluation', ['evaluate_preprocessing.py', '--tag', 'screen']),
        ('preprocessing_confirm', ['preprocessing_experiment.py', '--scales', '200', '--variants',
                                   'baseline', 'central_band20', 'core_row_median', '--tag', 'confirm']),
        ('preprocessing_confirm_evaluation', ['evaluate_preprocessing.py', '--tag', 'confirm']),
        ('preprocessing_v5', ['compare_v5.py']),
        ('band_comparison', ['compare_bands.py']),
        ('null_diagnostic', ['null_diagnostic.py', '--surrogates', '199']),
    ]
    env = dict(os.environ, MPLCONFIGDIR=os.environ.get('MPLCONFIGDIR', '/tmp/v16-mpl'))
    rows = []
    for name, args in tasks:
        start = time.perf_counter()
        print(f'Início: {name}', flush=True)
        with (ROOT/'results'/f'{name}.log').open('w') as log:
            subprocess.run([sys.executable, str(ROOT/'scripts'/args[0]), *args[1:]],
                           env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        elapsed = time.perf_counter()-start
        rows.append(dict(stage=name, elapsed_seconds=elapsed))
        print(f'Concluído: {name}, {elapsed:.1f} s', flush=True)
        (ROOT/'results/recalculo_sem_mascaras.json').write_text(
            json.dumps(dict(domain='integral, sem COI e sem máscara conservadora',
                            completed_stages=rows), ensure_ascii=False, indent=2)+'\n')

if __name__ == '__main__':
    main()
