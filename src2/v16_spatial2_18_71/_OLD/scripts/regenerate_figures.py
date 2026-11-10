"""Regenera todas as figuras e os relatórios sem repetir os lotes estatísticos."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from plot_labels import PATTERN_PT

ROOT=Path(__file__).resolve().parents[1]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    numeric=list((ROOT/'results').glob('*.csv'))+list((ROOT/'intermediate').glob('*.npz'))
    numeric.append(ROOT/'scripts/wavelet_core.py')
    before={str(p.relative_to(ROOT)):sha(p) for p in numeric}
    figures=list((ROOT/'figures').rglob('*.png'))
    old_figures={str(p.relative_to(ROOT)):(sha(p),p.stat().st_mtime_ns) for p in figures}
    env=dict(os.environ,MPLCONFIGDIR=os.environ.get('MPLCONFIGDIR','/tmp/v16-mpl'))
    stages=[]
    for script in ['render_saved.py','additional_diagnostics.py','precision_and_windows.py',
                   'report_preprocessing.py','finish_outputs.py']:
        start=time.perf_counter()
        print(f'Regenerando: {script}',flush=True)
        with (ROOT/'results'/f'rerender_{Path(script).stem}.log').open('w') as log:
            subprocess.run([sys.executable,str(ROOT/'scripts'/script)],env=env,
                           stdout=log,stderr=subprocess.STDOUT,check=True)
        stages.append(dict(script=script,elapsed_seconds=time.perf_counter()-start))
    changed_numeric=[name for name,value in before.items() if sha(ROOT/name)!=value]
    assert not changed_numeric,changed_numeric
    assert all((ROOT/name).stat().st_mtime_ns!=value[1] for name,value in old_figures.items())
    changed_figures=[name for name,value in old_figures.items() if sha(ROOT/name)!=value[0]]
    result=dict(status='PASS',pattern_labels={k:PATTERN_PT[k] for k in
        ['dispersed bubbles','plug','slug','churn','wavy annular']},other_text_language='pt-BR',
        coi=False,conservative_mask=False,regenerated_figures=len(figures),
        changed_figures=len(changed_figures),unchanged_numeric_files=len(before),
        numerical_data_and_wavelet_core_unchanged=True,stages=stages,
        numeric_sha256=before,previous_figure_sha256={k:v[0] for k,v in old_figures.items()})
    (ROOT/'results/revisao_nomes_padroes.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(f'Concluído: {len(figures)} figuras regeneradas; {len(changed_figures)} imagens alteradas; dados numéricos preservados.',flush=True)

if __name__=='__main__':main()
