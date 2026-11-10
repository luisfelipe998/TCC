"""Optical-only preliminary window audit; no wavelet values enter selection."""
import sys
sys.dont_write_bytecode = True
import json
import numpy as np
from common import CFG, INTER, context_ok, load_image, read_csv, slug, write_csv, write_json

WINDOWS = [1024, 1280, 1536, 2048, 3072, 4096, 6144, 8192]
OUT = INTER / "janelas"
BASELINE = INTER / "historico/piloto_inicial/dados"


def main():
    OUT.mkdir(exist_ok=True)
    previous = json.loads((BASELINE / "selecao.json").read_text())
    catalog = read_csv(BASELINE / "catalogo_candidatos.csv")
    initial = {r["aquisicao"]: r for r in previous["ancoras"]}
    records = []
    for old in previous["ancoras"]:
        run, pattern = old["aquisicao"], old["padrao"]
        _, raw = load_image(run)
        rows = [r for r in catalog if r["aquisicao"] == run]
        stable = [r for r in rows if r.get("elegivel_individual", r["elegivel"]) == "True"]
        centers = sorted(int(r["ancora"]) for r in stable)
        spacing = float(np.median(np.diff(centers))) if len(centers) >= 3 else None
        for length in WINDOWS:
            eligible = [r for r in stable if length >= int(r["comprimento_minimo"])
                        and context_ok(int(r["ancora"]), len(raw), length)]
            key = (lambda r: (float(r["variacao_maxima"]), int(r["ancora"]))) if pattern in ("plug", "slug") else (
                lambda r: (-float(r["prominencia"]), float(r["variacao_maxima"]), int(r["ancora"])))
            eligible.sort(key=key)
            record = dict(padrao=pattern, aquisicao=run, comprimento=length,
                          fracao_panorama=length/len(raw), candidatos_validos=len(eligible),
                          distancia_mediana_marcos=spacing,
                          ancora_inicial=old["anchor_id"], janela_inicial=old["comprimento"],
                          ancora_inicial_cabe=context_ok(old["ancora"], len(raw), length))
            if eligible:
                best = eligible[0]
                center = int(best["ancora"])
                start, stop = center-length//2, center+length//2
                complete = [r for r in stable if float(r["inicio"]) >= start and float(r["fim"]) <= stop]
                marks = [int(r["ancora"]) for r in complete]
                thirds = np.array_split(raw[start:stop], 4)
                record.update(anchor_id=best["anchor_id"], ancora=center, segmento_inicio=start,
                              segmento_fim=stop, marcos_completos=len(complete),
                              intervalos_entre_marcos=max(0, len(complete)-1),
                              amplitudes_medias_quartos=float(np.ptp([r.mean() for r in thirds])),
                              desvio_sinal=float(raw[start:stop].std()),
                              proeminencia=float(best["prominencia"]) if best.get("prominencia") else None,
                              variacao_maxima=float(best["variacao_maxima"]),
                              cobertura_marcos=(max(marks)-min(marks))/length if len(marks) > 1 else 0.)
            records.append(record)
    write_csv(OUT / "comparacao_janelas.csv", records)
    write_json(OUT / "comparacao_janelas.json", records)
    render(records, initial)
    for pattern in CFG["padroes"]:
        print(pattern)
        for n in WINDOWS:
            group = [r for r in records if r["padrao"] == pattern and r["comprimento"] == n]
            print(n, [(r["aquisicao"], r["candidatos_validos"], r.get("anchor_id"),
                       r.get("marcos_completos"), round(r["fracao_panorama"]*100, 1)) for r in group])


def render(records, initial):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    for pattern in CFG["padroes"]:
        fig, axes = plt.subplots(2, 1, figsize=(8, 6), constrained_layout=True)
        group = [r for r in records if r["padrao"] == pattern]
        for run in sorted({r["aquisicao"] for r in group}):
            data = [r for r in group if r["aquisicao"] == run]
            xs = [r["comprimento"] for r in data]
            axes[0].plot(xs, [r["candidatos_validos"] for r in data], "o-", label=run)
            axes[1].plot(xs, [r.get("marcos_completos", np.nan) for r in data], "o-", label=run)
        axes[0].set_ylabel("Âncoras válidas")
        axes[1].set_ylabel("Marcos no recorte")
        for ax in axes:
            ax.set_xlabel("Janela (px)")
            ax.grid(alpha=.2)
            ax.legend(frameon=False, fontsize=8)
        fig.suptitle(f"{pattern} · Contexto e disponibilidade")
        fig.savefig(OUT / f"{slug(pattern)}_disponibilidade.png", dpi=180, bbox_inches="tight", pad_inches=.05)
        plt.close(fig)
    # The user's concern: compare the actual wavy-annular profiles at increasing context.
    runs = sorted(run for run, r in initial.items() if r["padrao"] == "wavy annular")
    fig, axes = plt.subplots(4, 2, figsize=(11, 8), constrained_layout=True)
    for i, n in enumerate((1024, 2048, 4096, 6144)):
        for j, run in enumerate(runs):
            r = next(r for r in records if r["aquisicao"] == run and r["comprimento"] == n)
            ax = axes[i, j]
            if not r["candidatos_validos"]:
                ax.text(.5, .5, "Sem âncora válida", ha="center", transform=ax.transAxes)
                continue
            _, raw = load_image(run)
            start, stop, center = r["segmento_inicio"], r["segmento_fim"], r["ancora"]
            ax.plot(np.arange(start, stop)-center, raw[start:stop], color="#246a75", lw=.65)
            ax.axvline(0, color="#d5852c", ls="--", lw=.8)
            ax.set_title(f"{run} · {n} px · {r['marcos_completos']} marcos", fontsize=9)
            ax.set_xlabel("Posição relativa (px)", fontsize=8)
            ax.set_ylabel("Cinza", fontsize=8)
    fig.suptitle("wavy annular · Avaliação preliminar das janelas")
    fig.savefig(OUT / "wavy_annular_perfis.png", dpi=180, bbox_inches="tight", pad_inches=.05)
    plt.close(fig)


if __name__ == "__main__":
    main()
