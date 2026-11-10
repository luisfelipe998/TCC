"""Delete only generated v16.3 outputs after an explicit request for a fresh run."""
import sys
sys.dont_write_bytecode = True
import argparse
import shutil
from datetime import datetime, timezone
from pathlib import Path

from common import ROOT, INTER, CONCLUSIONS, setup, write_json


def reset_outputs():
    if ROOT.name != "v16.3_ancoras" or ROOT != Path(__file__).resolve().parents[1]:
        raise ValueError("Raiz inesperada; limpeza recusada.")
    if CONCLUSIONS.exists() and any(CONCLUSIONS.iterdir()):
        raise ValueError("A pasta de conclusoes contem arquivos e nao faz parte desta limpeza.")
    targets = [ROOT/name for name in ("01_ancoras", "02_wavelets", "intermediario")]
    inventory = []
    for folder in targets:
        if folder.is_symlink() or folder.resolve().parent != ROOT:
            raise ValueError("Destino de limpeza fora da analise.")
        files = [p for p in folder.rglob("*") if p.is_file()] if folder.exists() else []
        inventory.append(dict(pasta=folder.name, arquivos=len(files), bytes=sum(p.stat().st_size for p in files)))
    started = datetime.now(timezone.utc).isoformat()
    for folder in targets:
        if folder.exists():
            shutil.rmtree(folder)
    setup()
    write_json(INTER/"reinicio.json", dict(inicio_utc=started,
        fim_limpeza_utc=datetime.now(timezone.utc).isoformat(), escopo="piloto",
        pastas_limpas=inventory, calculos_anteriores_reutilizados=0,
        entradas_preservadas=["scripts", "data_processed", "v16.2_wco_eventos", "v16_spatial2_18_71", "assets"],
        estimador="biased", ancoras_previstas=8, pares_previstos=4))
    print(f"Limpeza concluida: {sum(r['arquivos'] for r in inventory)} arquivos gerados removidos; entradas preservadas.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executar", action="store_true", required=True,
                        help="Usar somente quando a remocao dos resultados tiver sido solicitada.")
    parser.parse_args()
    reset_outputs()
