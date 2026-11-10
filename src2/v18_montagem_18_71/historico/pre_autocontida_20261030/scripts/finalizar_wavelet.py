"""Consolida rastreabilidade e comparação v16→v18 sem recalcular resultados."""
import json
from pathlib import Path
import platform
import sys

sys.dont_write_bytecode = True
import recalcular_wavelet as computation
from recalcular_wavelet import ROOT, REPO, OLD, BASE, TABLES, INTER, read_csv, sha, write_json, source_record
import scipy
import matplotlib
import pywt


def main():
    summary = json.loads((ROOT / "resumo_tcc.json").read_text())
    assert summary["execucao_pares_wavelet_completa"]
    old_records = {r["ensaio"]: r for r in json.loads((OLD / "intermediario/aquisicoes.json").read_text())}
    old_groups = {(r["padrao_x"], r["padrao_y"]): r
                  for r in read_csv(OLD / "04_estatistica/dados/resumo_por_padrao.csv")}
    old_averages = {r["ensaio"]: r for r in read_csv(OLD / "04_estatistica/dados/medias_por_ensaio.csv")}
    old_pairs = read_csv(OLD / "04_estatistica/dados/pares_pearson_wco.csv")
    new_pairs = read_csv(TABLES / "pares_pearson_wco.csv")
    old_peak = max((r for r in old_pairs if r["tipo"] == "inter"), key=lambda r: float(r["wco_maxima"]))
    new_peak = max((r for r in new_pairs if r["tipo"] == "inter"), key=lambda r: float(r["wco_maxima"]))
    individuals = json.loads((INTER / "individuais.json").read_text())
    comparison = dict(
        escopo="Somente mudança da faixa usada para estimar os deslocamentos da montagem. Núcleo espectral e definições estatísticas preservados.",
        globais=summary["comparacao_antigo_novo"],
        extensoes=[dict(ensaio=r["ensaio"], padrao=r["padrao"],
                       anterior=old_records[r["ensaio"]]["comprimento"], atual=r["comprimento"],
                       diferenca=r["comprimento"] - old_records[r["ensaio"]]["comprimento"])
                   for r in summary["ensaios"]],
        grupos=[dict(padrao_x=r["padrao_x"], padrao_y=r["padrao_y"], pares=r["pares"],
                     media_anterior=float(old_groups[r["padrao_x"], r["padrao_y"]]["pearson_media"]),
                     media_atual=r["pearson_media"])
                for r in summary["grupos"]],
        medias_intra_por_ensaio=[dict(ensaio=r["ensaio"], padrao=r["padrao"],
                                    anterior=float(old_averages[r["ensaio"]]["pearson_media_intra"])
                                    if old_averages[r["ensaio"]]["pearson_media_intra"] else None,
                                    atual=r["pearson_media_intra"])
                                for r in summary["medias_por_ensaio"]],
        inter_maximo_anterior=old_peak, inter_maximo_atual=new_peak,
        pares_intra_atual=summary["pares_intrapadrao_selecionados"],
        picos_individuais=[{k: r[k] for k in ("ensaio", "padrao", "pico_wps_zero_lambda_px", "pico_potencia_media_lambda_px")}
                           for r in sorted(individuals, key=lambda r: r["ensaio"])],
        interpretacoes={
            "media_intra_maior_que_inter": summary["diferenca_global"] > 0,
            "diferenca_das_medias_diminuiu": summary["diferenca_global"] < summary["comparacao_antigo_novo"]["diferenca_global"]["antigo"],
            "maiores_medias_intra": [r["padrao_x"] for r in sorted(
                (r for r in summary["grupos"] if r["tipo"] == "intra"), key=lambda r: r["pearson_media"], reverse=True)][:2],
            "menor_media_intra_anular": summary["menor_anular"],
            "par_maximo_inter_mudou": (old_peak["ensaio_x"], old_peak["ensaio_y"]) != (new_peak["ensaio_x"], new_peak["ensaio_y"]),
            "descricao_atual_inter": "slug m04q01 / anular ondulado m05q80, conforme máximo pontual da coerência nos 238 pares entre padrões",
            "inspecao_visual_perfis": "Plug manteve perfis transformados próximos. Slug manteve uma curva discrepante (m06q05). O perfil anular m03q40 continuou com pico WPS em menor comprimento de onda, 171,83 px, enquanto os outros oito picos ficaram entre 353,05 e 480,69 px.",
            "distribuicoes": "Continuou a ampla sobreposição dos coeficientes intra e interpadrões. Não foi calculada significância, acurácia ou ganho de classificação.",
        },
    )
    write_json(ROOT / "comparacao_tcc.json", comparison)
    figures_manifest = REPO / "TCC/figuras/v18_wavelet/origens.json"
    figures = json.loads(figures_manifest.read_text())
    for entry in figures["figuras"]:
        for kind in ("origem", "destino"):
            record = entry[kind]
            assert sha(REPO / record["arquivo"]) == record["sha256"]
    source_paths = [OLD / "scripts/wavelet_core.py", OLD / "scripts/estatistica.py",
                    OLD / "scripts/figuras.py", REPO / "TCC/figuras/v16/gerar_perfis_pearson.py",
                    REPO / "TCC/figuras/v16/gerar_figuras_sem_escore.py"]
    assert sha(source_paths[0]) == computation.CORE_HASH
    manifest = dict(status="PASS", complete=True, stage="02_wavelet",
        inputs=source_record(BASE / "01_panoramas/manifesto.json"),
        historical_sources=[source_record(p) for p in source_paths],
        scripts=[source_record(BASE / "scripts" / name) for name in
                 ("recalcular_wavelet.py", "gerar_figuras_wavelet.py", "finalizar_wavelet.py")],
        figures_manifest=source_record(figures_manifest),
        runtime=dict(python=platform.python_version(), numpy=computation.np.__version__,
                     scipy=scipy.__version__, matplotlib=matplotlib.__version__, pywavelets=pywt.__version__,
                     workers=2, fft_threads=1),
        counts=dict(individuals=26, pairs=325, intra=87, inter=238,
                    wavelengths=200, compact_pair_maps=5, figures_in_manuscript=10),
        files={str(p.relative_to(ROOT)): sha(p) for p in sorted(ROOT.rglob("*"))
               if p.is_file() and p.name != "manifesto.json"},
    )
    write_json(ROOT / "manifesto.json", manifest)
    print(json.dumps({"status": "PASS", "files": len(manifest["files"]),
                      "manifesto": str(ROOT / "manifesto.json")}, ensure_ascii=False))


if __name__ == "__main__":
    main()
