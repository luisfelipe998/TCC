#!/usr/bin/env python3
"""Recalcula a metodologia completa em uma saída nova e verifica os resultados."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

PACKAGE = Path(__file__).resolve().parent


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def sources():
    paths = [PACKAGE / 'executar.py']
    for directory in ('scripts', '04_analise_simplificada/scripts', 'config'):
        paths.extend(path for path in (PACKAGE / directory).rglob('*')
                     if path.is_file() and '__pycache__' not in path.parts)
    return {path.relative_to(PACKAGE).as_posix(): digest(path) for path in sorted(paths)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, default=PACKAGE.parent / 'data_processed')
    parser.add_argument('--output', type=Path,
                        help='Diretório vazio; padrão: execucoes/<data e hora> dentro da v18.')
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--sem-resumos', action='store_true',
                        help='Omite apenas a etapa complementar06 de resumos de Pearson e coerência.')
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('--workers deve ser positivo.')
    data = args.data_root.resolve()
    if not data.is_dir():
        parser.error(f'Diretório de dados não encontrado: {data}')
    timestamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    output = (args.output or PACKAGE / 'execucoes' / timestamp).resolve()
    if output == PACKAGE or output.is_relative_to(data) or data.is_relative_to(output):
        parser.error('A saída deve ser nova e separada dos dados e do código ativo.')
    if any(output.is_relative_to(PACKAGE / directory) for directory in
           ('scripts', 'config', '04_analise_simplificada/scripts')):
        parser.error('A saída não pode ficar dentro das pastas de código ou configuração.')
    if output.exists() and any(output.iterdir()):
        parser.error(f'A saída já contém arquivos: {output}')
    required = [PACKAGE / 'config' / name for name in ('ensaios.json', 'referencia_resultados.json')]
    if not all(path.is_file() for path in required):
        parser.error('Configuração ou referência de validação ausente; pacote incompleto.')
    output.mkdir(parents=True, exist_ok=True)
    logs = output / 'logs'
    logs.mkdir()
    env = dict(os.environ, TCC_V18_OUTPUT=str(output), TCC_DATA_ROOT=str(data),
               PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
               MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
    env.setdefault('MPLCONFIGDIR', str(output / 'cache' / 'matplotlib'))
    steps = [
        ('montagem', 'scripts/reconstruir.py', []),
        ('comparar_panoramas', 'scripts/comparar_resultados.py', ['--stage', 'panoramas']),
        ('wavelet', 'scripts/recalcular_wavelet.py', ['--workers', str(args.workers)]),
        ('comparar_wavelet', 'scripts/comparar_resultados.py', ['--stage', 'wavelet']),
        ('optica', 'scripts/process_optical.py', []),
        ('comparar_optica', 'scripts/comparar_resultados.py', ['--stage', 'optica']),
        ('pearson_e_mapas', '04_analise_simplificada/scripts/analisar.py', []),
        ('figuras_opticas', 'scripts/gerar_figuras_opticas.py', []),
        ('finalizar_optica', 'scripts/finalize_optical.py', []),
        ('finalizar_wavelet', 'scripts/finalizar_wavelet.py', []),
        ('figuras_panoramas', 'scripts/gerar_panoramas.py', []),
        ('figuras_metodologia', 'scripts/gerar_metodologia.py', []),
        ('comparar_completo', 'scripts/comparar_resultados.py', ['--stage', 'completo']),
    ]
    if not args.sem_resumos:
        steps.append(('resumos_pearson_coerencia', 'scripts/resumir_pares.py',
                      ['--workers', str(args.workers)]))
    before = sources()
    manifest = dict(status='RUNNING', started_utc=timestamp, stages=[],
                    external_data=str(data), source_sha256=before,
                    method='METODOLOGIA.md', workers=args.workers,
                    environment=dict(python=platform.python_version(), libraries={
                        name: version(name) for name in ('numpy', 'scipy', 'Pillow', 'matplotlib', 'PyWavelets')}),
                    production_outputs_automatically_replaced=False)
    manifest_path = output / 'manifesto_execucao.json'

    def save():
        manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')

    save()
    started = time.monotonic()
    try:
        for name, script, arguments in steps:
            print(f'ETAPA {name}', flush=True)
            start = time.monotonic()
            log = logs / f'{name}.log'
            with log.open('w') as stream:
                result = subprocess.run([sys.executable, '-B', str(PACKAGE / script), *arguments],
                                        cwd=PACKAGE, env=env, stdout=stream,
                                        stderr=subprocess.STDOUT, check=False)
            manifest['stages'].append(dict(name=name, script=script, arguments=arguments,
                                          returncode=result.returncode,
                                          seconds=time.monotonic()-start,
                                          log=log.relative_to(output).as_posix()))
            save()
            if result.returncode:
                print(log.read_text()[-10000:], file=sys.stderr, flush=True)
                raise RuntimeError(f'Etapa {name} interrompida; consulte {log}')
            print(f'PASS {name} ({time.monotonic()-start:.1f}s)', flush=True)
        if sources() != before:
            raise RuntimeError('Código/configuração mudou durante a execução; execute novamente com fontes estáveis.')
        manifest.update(status='PASS', complete=True, sources_preserved=True,
                        seconds=time.monotonic()-started)
        save()
    except Exception as error:
        manifest.update(status='FAIL', complete=False, error=str(error),
                        seconds=time.monotonic()-started)
        save()
        print(f'INTERROMPIDO: {error}\nA saída foi preservada para inspeção: {output}', file=sys.stderr)
        return 1
    print(f'COMPLETO: {output}', flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
