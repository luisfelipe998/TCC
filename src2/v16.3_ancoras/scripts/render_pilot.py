"""Replace the 12 pilot PNGs using the frozen anchors and existing wavelet maps."""
import sys
sys.dont_write_bytecode = True
import argparse
import json
from datetime import datetime, timezone

from common import (ANCHORS, INTER, WAVELETS, load_image, selected_manifest,
                    sha, slug, write_json)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limpar", action="store_true", help="Remove os PNGs finais atuais antes de regenerar.")
    args = parser.parse_args()
    selection = selected_manifest()
    assert len(selection["ancoras"]) == 8 and len(selection["pares"]) == 4
    preserved = [INTER / "selecao.json", INTER / "congelamento.json", INTER / "catalogo_candidatos.csv",
                 INTER / "metricas_wavelets.json", INTER / "rotulos_visuais.json"]
    preserved += sorted((INTER / "mapas").glob("*.npz"))
    before = {str(p.relative_to(INTER)): sha(p) for p in preserved}
    removed = []
    if args.limpar:
        for folder in (ANCHORS, WAVELETS):
            for path in sorted(folder.glob("*.png")):
                removed.append(str(path.relative_to(folder.parent)))
                path.unlink()
    from plots import anchor_figure
    for record in selection["ancoras"]:
        pano, raw = load_image(record["aquisicao"])
        anchor_figure(record, pano, raw)
        print(f"Âncora: {record['anchor_id']}", flush=True)
    from compute_wavelets import render_final
    render_final()
    assert before == {str(p.relative_to(INTER)): sha(p) for p in preserved}
    diagnostics = [json.loads((INTER / "apresentacao" / (path.stem + ".json")).read_text())
                   for folder in (ANCHORS, WAVELETS) for path in sorted(folder.glob("*.png"))]
    for record in diagnostics:
        assert "favorável" not in record["titulo"].lower()
        assert "desfavorável" not in record["titulo"].lower()
        assert all("(cinza)" not in label for label in record["rotulos"])
    for comparison in selection["pares"]:
        labels = json.loads((INTER / "rotulos_visuais.json").read_text())
        name = f"{slug(comparison['padrao'])}_{comparison['run_x']}_{comparison['run_y']}_{labels[comparison['padrao']]['rotulo']}"
        info = json.loads((INTER / "apresentacao" / f"{name}.json").read_text())
        assert "\n" not in info["titulo"]
        assert info["titulos_paineis"][:2] == [f"A - {comparison['run_x']}", f"B - {comparison['run_y']}"]
    write_json(INTER / "regeneracao_apresentacao.json", dict(
        data_utc=datetime.now(timezone.utc).isoformat(), pngs_removidos=removed,
        pngs_regenerados=len(diagnostics), dados_preservados_sha256=before,
        altura_maxima_a_16cm=max(r["altura_inserida_cm"] for r in diagnostics),
        fonte_minima_a_16cm=min(r["fonte_minima_inserida_pt"] for r in diagnostics)))
    print("12 PNGs regenerados; âncoras, janelas e mapas numéricos preservados.")


if __name__ == "__main__":
    main()
