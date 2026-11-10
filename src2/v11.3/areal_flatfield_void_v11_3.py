"""
v11.3 (B — SUPLEMENTAR) — Vazio AREAL flat-field × intermitência (flip-rate).

Análise IRMÃ suplementar do par v11.3. Serve a dois propósitos:
  (i) PROVAR que a normalização remove o reflexo óptico W/M (figura de prova:
      o perfil radial colapsa ~100–300× após o flat-field);
  (ii) dar a métrica de vazio AREAL 2-D (fração de ÁREA da seção, schepperle2024)
      que o dono pediu — algo que o bright-core (cordal) não é.

A normalização (flat-field por BASELINE DE LÍQUIDO = percentil 20 por pixel) é a
variante validada: divide cada quadro por um fundo dominado por líquido que NÃO
acompanha o núcleo de gás móvel, e classifica gás = excesso claro (Fn > 1,3).

LIMITE FÍSICO (documentado, NÃO é bug): mesmo normalizado, o vazio AREAL por
limiar de área clara **inverte churn/annular** — o annular NÃO fica no topo.
Em wavy annular as bandas de interface ESCURAS (linhas ~20/60) persistem mesmo
quando o núcleo claro enche o centro, então o limiar de área clara SUBCONTA o
annular. É limite da vista única + óptica W/M. A ordenação correta annular>churn
só é recuperável pelo bright-core (script A, principal).

RESSALVAS (obrigatórias, repetidas na tabela):
  1. Proxy óptico RELATIVO não-calibrado (sem px→mm, vista única).
  2. Recupera dispersed < plug < slug < churn, mas **annular ABAIXO de churn**
     (não fica no topo) — limite óptico.
  3. dispersed bubbles instável em n=2 (os dois ensaios discordam muito).
  4. Correlação de brilho/saturação residual ~0,3; q-confound não testado.

Não altera v11.1 nem v11.2. Reusa ``void_common_v11_3.py``.

Saídas (em src2/v11.3/):
  areal_flatfield_void_table_v11_3.md
  areal_flatfield_void_map_v11_3.png   (X = vazio areal, Y = flip-rate)
  flatfield_radial_proof_v11_3.png     (perfil radial bruto × flat-field)

Uso:
  cd src2/v11.3 && python3 areal_flatfield_void_v11_3.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import void_common_v11_3 as C


SHIFT = {"bubble": (0.0, 0.004), "plug": (0.0, 0.004), "slug": (0.03, 0.0),
         "churn": (0.0, 0.004), "annular": (0.0, -0.004)}


def main():
    scenarios = C.list_scenarios()
    print(f"{len(scenarios)} ensaios de padrão puro (dispersed bubbles n=2 "
          "incluído)\n")

    rows = []
    for name in scenarios:
        pano = C.load_panorama(name)
        ch = C.alpha_chordal(pano)
        flip = C.temporal_flip_rate(pano)
        frames = C.load_raw_frames(name)

        if frames is not None:
            areal = C.areal_flatfield_void(frames)
            raw_std, flat_std = C.radial_std_rawvsflat(frames)
        else:
            areal, raw_std, flat_std = np.nan, np.nan, np.nan

        label = C.folder_flow_type(name)
        group = C.GROUP_MAP.get(label, label)
        m, q = C.parse_m_q(name)
        dropped = ch["sat_frac"] > C.SAT_LIMIT
        flag = (f"saturado (sat={ch['sat_frac'] * 100:.0f}%)" if dropped else "")

        rows.append({
            "scenario": name, "pattern": label, "group": group, "m": m, "q": q,
            "G": C.G_MAP[m], "areal": areal, "flip_rate": flip["flip_rate"],
            "sat_frac": ch["sat_frac"], "core_mean": ch["core_mean"],
            "raw_std": raw_std, "flat_std": flat_std,
            "N": int(pano.shape[1]), "dropped": dropped, "flag": flag,
        })
        print(f"  {name:14s} {label:18s} areal={areal:.3f} "
              f"flip={flip['flip_rate']:.4f} radial_std {raw_std:.1f}→{flat_std:.2f}"
              + ("  [SAT]" if dropped else ""))

    write_table(rows)
    C.make_void_flip_map(
        rows, "areal",
        "Fração de vazio AREAL flat-field (área de gás/área seção)  "
        "[–, relativo, não-calibrado]",
        "v11.3 (B — suplementar) — vazio AREAL flat-field × flip-rate\n"
        "reflexo W/M removido · annular CAI abaixo de churn (limite óptico, "
        "ver script A)",
        os.path.join(C.HERE, "areal_flatfield_void_map_v11_3.png"), SHIFT)
    make_radial_proof(rows)
    report_stats(rows)


def report_stats(rows):
    print("\n" + "=" * 72)
    print("VEREDITO — vazio AREAL flat-field (ordenação) + prova do reflexo")
    print("=" * 72)
    order = []
    for g in C.GROUP_ORDER:
        vals = C.grp(rows, g, "areal")
        if vals:
            order.append((g, float(np.median(vals)), len(vals)))
    order_sorted = sorted(order, key=lambda t: t[1])
    print("\nordem (baixo→alto): " +
          " < ".join(f"{C.GROUP_PRETTY[g]} {m:.3f}" for g, m, _ in order_sorted))
    ann = C.grp(rows, "annular", "areal")
    chu = C.grp(rows, "churn", "areal")
    if ann and chu:
        print(f"anular {np.median(ann):.3f} vs churn {np.median(chu):.3f} "
              f"⇒ annular {'ACIMA' if np.median(ann) > np.median(chu) else 'ABAIXO'}"
              " de churn (esperado: ABAIXO — limite óptico)")
    raw = np.array([r["raw_std"] for r in rows if not np.isnan(r["raw_std"])])
    flat = np.array([r["flat_std"] for r in rows if not np.isnan(r["flat_std"])])
    print(f"\nprova do reflexo: std radial interior {raw.mean():.1f} (bruto) → "
          f"{flat.mean():.2f} (flat-field) — remoção média ~{raw.mean()/flat.mean():.0f}×")
    # correlação de brilho residual (ressalva #4)
    av = np.array([r["areal"] for r in rows if not r["dropped"]
                   and not np.isnan(r["areal"])])
    cm = np.array([r["core_mean"] for r in rows if not r["dropped"]
                   and not np.isnan(r["areal"])])
    sf = np.array([r["sat_frac"] for r in rows if not r["dropped"]
                   and not np.isnan(r["areal"])])
    if len(av) > 3:
        print(f"corr(vazio, brilho núcleo) = {np.corrcoef(av, cm)[0,1]:+.2f}; "
              f"corr(vazio, saturação) = {np.corrcoef(av, sf)[0,1]:+.2f} "
              "(ressalva: dependência residual de iluminação)")


def write_table(rows):
    out = os.path.join(C.HERE, "areal_flatfield_void_table_v11_3.md")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("# Vazio AREAL flat-field × flip-rate (v11.3 B — SUPLEMENTAR)\n\n")
        fh.write(
            "Análise **suplementar** do par v11.3 (a principal é o bright-core, "
            "script A). Eixo X = fração de vazio **AREAL 2-D** (área de gás / área "
            "da seção, schepperle2024) sobre o quadro **flat-field por baseline de "
            "líquido** (percentil 20 por pixel; gás = excesso claro Fn>1,3). Eixo Y "
            "= taxa de troca temporal (flip-rate, v11.2 M4).\n\n")
        fh.write(
            "**Dois propósitos:** (i) PROVAR que a normalização remove o reflexo "
            "óptico W/M (o perfil radial colapsa ~100–300×, ver "
            "`flatfield_radial_proof_v11_3.png` e a coluna std radial); (ii) dar a "
            "métrica de ÁREA que o bright-core cordal não é.\n\n")
        fh.write(
            "> **LIMITE FÍSICO (não é bug):** mesmo normalizado, o limiar de área "
            "clara **inverte churn/annular** — o annular NÃO fica no topo. As "
            "bandas de interface ESCURAS do wavy annular (linhas ~20/60) persistem "
            "e subcontam o annular. A ordenação correta anular>churn só vem do "
            "bright-core (script A).\n\n")

        n_in = sum(1 for r in rows if not r["dropped"])
        fh.write(f"Padrões: dispersed bubbles, plug, slug, churn, wavy annular — "
                 f"{len(rows)} ensaios ({n_in} na estatística).\n\n")

        fh.write("## Mediana por padrão (saturados fora)\n\n")
        ar = dict((g, m) for g, _, m in C.per_pattern_median(rows, "areal"))
        fr = dict((g, m) for g, _, m in C.per_pattern_median(rows, "flip_rate"))
        ns = dict((g, n) for g, n, _ in C.per_pattern_median(rows, "areal"))

        def f3(x): return "—" if x is None else f"{x:.3f}"
        def f4(x): return "—" if x is None else f"{x:.4f}"
        C.md_table(fh, ["padrão", "n", "vazio AREAL X", "flip-rate Y (M4)"],
                   [[C.GROUP_PRETTY[g], ns[g], f3(ar[g]), f4(fr[g])]
                    for g in C.GROUP_ORDER])
        ordered = sorted([(g, ar[g]) for g in C.GROUP_ORDER if ar[g] is not None],
                         key=lambda t: t[1])
        fh.write("\n**Ordem recuperada (baixo→alto):** "
                 + " < ".join(f"{C.GROUP_PRETTY[g]} {m:.3f}" for g, m in ordered)
                 + ". Correta até slug; o annular cai abaixo do churn (limite "
                   "óptico documentado acima).\n\n")

        fh.write("## Prova da remoção do reflexo W/M\n\n")
        fh.write("Desvio-padrão radial do perfil de intensidade (col-mean por "
                 "linha, média sobre quadros) na janela interior, BRUTO vs "
                 "FLAT-FIELD. O W/M bruto tem std radial grande; após o flat-field "
                 "colapsa ~100–300×. Ver `flatfield_radial_proof_v11_3.png`.\n\n")
        rr = sorted([r for r in rows if not np.isnan(r["raw_std"])],
                    key=lambda r: (r["group"], r["scenario"]))
        C.md_table(fh, ["scenario", "padrão", "std radial BRUTO", "std radial FLAT",
                        "razão"],
                   [[r["scenario"], r["pattern"], f"{r['raw_std']:.1f}",
                     f"{r['flat_std']:.2f}",
                     f"{r['raw_std']/r['flat_std']:.0f}×"] for r in rr])

        fh.write("\n## Por ensaio (transparência)\n\n")
        rr = sorted(rows, key=lambda r: -(r["areal"]
                    if not np.isnan(r["areal"]) else -1))
        C.md_table(fh, ["scenario", "pattern", "G", "q", "vazio areal", "flip-rate",
                        "sat %", "flag"],
                   [[r["scenario"], r["pattern"], r["G"], r["q"],
                     f"{r['areal']:.3f}", f"{r['flip_rate']:.4f}",
                     f"{r['sat_frac']*100:.0f}", r["flag"]] for r in rr])

        fh.write("\n## Ressalvas honestas\n\n")
        fh.write("- **Proxy óptico RELATIVO não-calibrado** (sem px→mm, vista "
                 "única) — não é a fração de área absoluta da seção.\n")
        fh.write("- **annular ABAIXO de churn** (não fica no topo): limite óptico "
                 "da vista única — as bandas de interface escuras do annular "
                 "derrotam o limiar de área clara. A ordenação anular>churn só é "
                 "recuperável pelo bright-core (script A).\n")
        fh.write("- **dispersed bubbles instável (n=2):** os dois ensaios "
                 "discordam muito; não entra em comparação estatística.\n")
        fh.write("- **Dependência de iluminação residual** (corr ~0,3 com brilho/"
                 "saturação); q-confound não testado por subconjunto q-match.\n")
    print(f"\ntabela: {out}")


def make_radial_proof(rows):
    """Perfil radial (col-mean por linha, média sobre quadros) de um anular e um
    churn representativos, BRUTO vs FLAT-FIELD — prova que a normalização remove
    o W/M."""
    def pick(group):
        cand = [r for r in rows if r["group"] == group and not r["dropped"]]
        cand.sort(key=lambda r: abs(r["q"] - 20))
        return cand[0]["scenario"] if cand else None

    ann, chu = pick("annular"), pick("churn")
    if not (ann and chu):
        return
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.5))
    rows_sl = C.INTERIOR_ROWS
    yidx = np.arange(rows_sl.start, rows_sl.stop)
    for name, color, lab in [(ann, C.GROUP_COLOR["annular"], f"anular ({ann})"),
                             (chu, C.GROUP_COLOR["churn"], f"churn ({chu})")]:
        frames = C.load_raw_frames(name)
        win = frames[:, rows_sl, :].astype(float)
        raw_prof = win.mean(axis=(0, 2))
        base = C.liquid_baseline(frames)
        base = np.where(base < 1.0, 1.0, base)
        flat_prof = (win / base[None, :, :]).mean(axis=(0, 2))
        axL.plot(yidx, raw_prof, color=color, lw=1.8, label=lab)
        axR.plot(yidx, flat_prof, color=color, lw=1.8, label=lab)
    axL.set_title("BRUTO — perfil radial W/M\n(núcleo CLARO, interfaces ESCURAS)")
    axR.set_title("FLAT-FIELD (÷ baseline líquido)\nreflexo removido ~100–300×")
    for ax in (axL, axR):
        ax.axvspan(C.CENTER_ROW - 5, C.CENTER_ROW + 5, color="0.85", alpha=0.5,
                   label="núcleo (centro)")
        ax.set_xlabel("linha (seção transversal do tubo)")
        ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.4)
        ax.legend(loc="best", fontsize=8)
    axL.set_ylabel("intensidade média  [0–255]")
    axR.set_ylabel("intensidade normalizada  [–]")
    out = os.path.join(C.HERE, "flatfield_radial_proof_v11_3.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"figura: {out}")


if __name__ == "__main__":
    main()
