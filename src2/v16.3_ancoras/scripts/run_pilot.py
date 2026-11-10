"""Run explicit phases; never expand this pilot into all acquisition pairs."""
import argparse
import sys
import time
sys.dont_write_bytecode = True

from common import INTER, setup, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--etapa", required=True, choices=["ancoras", "preparar", "wavelets", "figuras", "verificar"])
    args = parser.parse_args()
    setup()
    started = time.perf_counter()
    if args.etapa == "ancoras":
        from select_anchors import select_pilot
        select_pilot()
    elif args.etapa == "preparar":
        from prepare_comparisons import prepare_all
        prepare_all(scope="piloto")
    elif args.etapa == "wavelets":
        from compute_wavelets import calculate
        calculate()
    elif args.etapa == "figuras":
        from compute_wavelets import render_final
        render_final()
    else:
        from verify_pilot import verify
        verify()
    write_json(INTER / f"execucao_{args.etapa}.json",
               dict(etapa=args.etapa, segundos=time.perf_counter()-started, status="concluida"))


if __name__ == "__main__":
    main()
