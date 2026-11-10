"""Retoma somente diagnósticos, preservando o lote e seus metadados."""
from datetime import datetime,timezone
from pathlib import Path
import hashlib
import json
import resource
import shutil
import time
from analysis_common import ROOT,load_signals
from texture_common import OUT,verify_freeze
from run_texture_analysis import diagnostics

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def atomic_json(path,value):
    temporary=path.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')
    temporary.replace(path)

def main():
    meta_path=OUT/'execucao.json';meta=json.loads(meta_path.read_text())
    if meta['status']=='PASS':
        print('Lote já concluído; nenhum cálculo repetido.');return
    assert meta['status']=='PAUSADO_PELO_USUARIO'
    checkpoint=json.loads((OUT/'checkpoint_pausa.json').read_text())
    for name,digest in checkpoint['sha256'].items():assert sha(ROOT/name)==digest,name
    verify_freeze()
    assert set(meta['etapas'])=={'anchors','fixed'}
    atomic_json(OUT/'verificacao_checkpoint_retomada.json',dict(status='PASS',
        data_utc=datetime.now(timezone.utc).isoformat(),arquivos_conferidos=len(checkpoint['sha256']),
        checkpoint_sha256=sha(OUT/'checkpoint_pausa.json'),sha256_conferidos=checkpoint['sha256']))
    for source,target in [(meta_path,OUT/'execucao_checkpoint_pausa.json'),
                          (ROOT/'RETOMADA.md',ROOT/'RETOMADA_CHECKPOINT.md')]:
        assert not target.exists();shutil.copyfile(source,target)
    meta.update(status='RETOMANDO_DIAGNOSTICOS',retomada_utc=datetime.now(timezone.utc).isoformat())
    atomic_json(meta_path,meta)
    start=time.perf_counter()
    try:
        _,signals=load_signals()
        meta['etapas']['diagnostics']=diagnostics(signals)
        verify_freeze()
    except BaseException:
        meta.update(status='INTERROMPIDO_NA_RETOMADA',interrupcao_retomada_utc=datetime.now(timezone.utc).isoformat())
        atomic_json(meta_path,meta);raise
    meta.update(status='PASS',fim_utc=datetime.now(timezone.utc).isoformat(),
        segundos_total=sum(s['segundos'] for s in meta['etapas'].values()),
        segundos_processo_retomada=time.perf_counter()-start,
        pico_memoria_mib=None,
        pico_memoria_retomada_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,
        nota_tempo='Soma das etapas concluídas; exclui pausa, pós-processamento e tempo parcial dos diagnósticos interrompidos, não cronometrado separadamente.',
        nota_memoria='O pico do primeiro processo não foi salvo antes da interrupção. O pico da retomada refere-se somente ao processo de diagnósticos.',
        retomada_script_sha256=sha(Path(__file__)))
    atomic_json(meta_path,meta)
    print(json.dumps(meta,indent=2,ensure_ascii=False))

if __name__=='__main__':main()
