"""
v11.3 (A — PRINCIPAL) — Fração de vazio BRIGHT-CORE × intermitência (flip-rate).

Esta é a análise PRINCIPAL e defensável do par v11.3. O eixo X é o vazio
**bright-core** (maior núcleo CLARO contíguo pelo centro, no quadro bruto) — o
método do v11.2 que foi verificado adversarialmente e RECUPERA a ordenação física
anular > churn (o α cordal "gás = escuro" do v11/v11.1 invertia: churn ≳ anular).

Por que bright-core e não "gás = pixel escuro": a óptica do microtubo de 1 mm
retroiluminado é **W/M** — o núcleo de gás aparece CLARO no centro (linha ~40), as
interfaces/meniscos aparecem ESCUROS fora do centro (linhas ~20/60). Contar o
núcleo CLARO respeita a física; contar pixel escuro conta a sombra de interface.

O eixo Y é a **taxa de troca temporal gás↔líquido** (flip-rate, v11.2 método 4):
eixo ORTOGONAL de ESTABILIDADE — churn (caótico) troca ~2× mais que anular
(núcleo estável). Imune a brilho/q.

RESSALVAS (obrigatórias, repetidas na tabela):
  1. Proxy óptico RELATIVO não-calibrado (sem px→mm, vista única) — a ordenação
     ENTRE padrões é legítima; o valor absoluto ~0,35 NÃO é a fração de área real
     (~0,9 areal da literatura).
  2. dispersed bubbles (n=2) sai alto ESPÚRIO (centro coincidentemente claro) —
     anotado, fora da estatística.
  3. plug empata com anular; sobreposição e n pequeno (anular 8 × churn 10).

Não altera v11.1 nem v11.2. Reusa ``void_common_v11_3.py`` (portado do v11.2).

Saídas (em src2/v11.3/):
  brightcore_void_table_v11_3.md
  brightcore_void_map_v11_3.png   (X = bright-core, Y = flip-rate, regiões/padrão)

Uso:
  cd src2/v11.3 && python3 brightcore_void_v11_3.py
"""

import os
import numpy as np

import void_common_v11_3 as C


SHIFT = {"bubble": (0.0, 0.004), "plug": (0.0, 0.004), "slug": (0.03, 0.0),
         "churn": (0.0, 0.004), "annular": (0.0, -0.004)}


def main():
    scenarios = C.list_scenarios()
    print(f"{len(scenarios)} ensaios de padrão puro (dispersed bubbles n=2 "
          "incluído; transições, synthetic*, smooth annular excluídos)\n")

    rows = []
    for name in scenarios:
        pano = C.load_panorama(name)
        ch = C.alpha_chordal(pano)
        flip = C.temporal_flip_rate(pano)
        frames = C.load_raw_frames(name)

        bc = C.bright_core_void(frames) if frames is not None else np.nan
        bc_mid = (C.bright_core_void(frames, use_midpoint=True)
                  if frames is not None else np.nan)
        axc = C.axial_continuity(frames) if frames is not None else np.nan

        label = C.folder_flow_type(name)
        group = C.GROUP_MAP.get(label, label)
        m, q = C.parse_m_q(name)
        dropped = ch["sat_frac"] > C.SAT_LIMIT
        flag = (f"saturado (sat={ch['sat_frac'] * 100:.0f}%)" if dropped else "")

        rows.append({
            "scenario": name, "pattern": label, "group": group, "m": m, "q": q,
            "G": C.G_MAP[m], "alpha_chordal": ch["alpha_chordal"],
            "core_mean": ch["core_mean"], "sat_frac": ch["sat_frac"],
            "bright_core": bc, "bright_core_mid": bc_mid, "axial_cont": axc,
            "flip_rate": flip["flip_rate"], "N": int(pano.shape[1]),
            "dropped": dropped, "flag": flag,
        })
        print(f"  {name:14s} {label:18s} bright_core={bc:.3f} "
              f"flip={flip['flip_rate']:.4f}" + ("  [SAT]" if dropped else ""))

    write_table(rows)
    C.make_void_flip_map(
        rows, "bright_core",
        "Fração de vazio bright-core (núcleo CLARO)  [–, cordal px, relativo]",
        "v11.3 (A) — vazio bright-core × taxa de troca temporal\n"
        "anular > churn no vazio (ordenação física recuperada) · "
        "churn mais instável (flip-rate)",
        os.path.join(C.HERE, "brightcore_void_map_v11_3.png"), SHIFT)
    report_stats(rows)


def report_stats(rows):
    print("\n" + "=" * 72)
    print("VEREDITO (anular × churn) — sat descartada")
    print("=" * 72)
    ann = C.grp(rows, "annular", "bright_core")
    chu = C.grp(rows, "churn", "bright_core")
    U, p = C.mannwhitney_u_greater(ann, chu)
    auc = C.auc_greater(ann, chu)
    print(f"\nbright-core (X):  annular {np.median(ann):.3f} > churn "
          f"{np.median(chu):.3f}  | MWU p(unilateral)={p:.4f} AUC={auc:.3f} "
          f"(n_ann={len(ann)}, n_chu={len(chu)})")
    # flip-rate: churn > annular (estabilidade)
    annf = C.grp(rows, "annular", "flip_rate")
    chuf = C.grp(rows, "churn", "flip_rate")
    Uf, pf = C.mannwhitney_u_greater(chuf, annf)
    aucf = C.auc_greater(chuf, annf)
    print(f"flip-rate (Y):    churn {np.median(chuf):.4f} > annular "
          f"{np.median(annf):.4f}  | MWU p(unilateral)={pf:.4f} AUC={aucf:.3f}")


def write_table(rows):
    out = os.path.join(C.HERE, "brightcore_void_table_v11_3.md")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("# Vazio bright-core × taxa de troca temporal (v11.3 A — PRINCIPAL)\n\n")
        fh.write(
            "Análise **principal** do par v11.3. Eixo X = fração de vazio "
            "**bright-core** (maior núcleo CLARO contíguo pelo centro, quadro "
            "bruto) — método do v11.2 verificado adversarialmente que RECUPERA "
            "anular > churn (o α cordal 'gás = escuro' do v11/v11.1 invertia por "
            "causa da óptica W/M: núcleo de gás CLARO, interfaces ESCURAS). "
            "Eixo Y = **taxa de troca temporal** gás↔líquido (flip-rate, v11.2 M4), "
            "eixo ORTOGONAL de estabilidade (churn caótico > anular estável).\n\n")

        n_in = sum(1 for r in rows if not r["dropped"])
        fh.write(f"Padrões: dispersed bubbles, plug, slug, churn, wavy annular — "
                 f"{len(rows)} ensaios ({n_in} na estatística).\n\n")

        fh.write("## Mediana por padrão (saturados fora)\n\n")
        bc = dict((g, m) for g, _, m in C.per_pattern_median(rows, "bright_core"))
        fr = dict((g, m) for g, _, m in C.per_pattern_median(rows, "flip_rate"))
        al = dict((g, m) for g, _, m in C.per_pattern_median(rows, "alpha_chordal"))
        ac = dict((g, m) for g, _, m in C.per_pattern_median(rows, "axial_cont"))
        ns = dict((g, n) for g, n, _ in C.per_pattern_median(rows, "bright_core"))

        def f3(x): return "—" if x is None else f"{x:.3f}"
        def f4(x): return "—" if x is None else f"{x:.4f}"
        C.md_table(fh,
            ["padrão", "n", "bright-core X (M1)", "flip-rate Y (M4)",
             "continuidade axial", "α cordal (v11.1, invertido)"],
            [[C.GROUP_PRETTY[g], ns[g], f3(bc[g]), f4(fr[g]), f3(ac[g]), f3(al[g])]
             for g in C.GROUP_ORDER])

        fh.write("\n**Leitura:** no eixo X (bright-core) **anular > churn** "
                 "(ordenação física recuperada; o α cordal invertido fica na "
                 "última coluna só para contraste). No eixo Y (flip-rate) "
                 "**churn > anular** — churn é o mais instável. Anular = vazio "
                 "ALTO + estável (canto inferior-direito); churn = vazio alto + "
                 "caótico (topo).\n\n")

        fh.write("## Controle de brilho — bright-core sem nível absoluto\n\n")
        fh.write("`bright-core-mid` usa limiar = (min+max)/2 por coluna, removendo "
                 "o NÍVEL de cinza. Se anular>churn sobrevive, a separação é de "
                 "FORMA do perfil, não de brilho.\n\n")
        bcm = dict((g, m) for g, _, m in C.per_pattern_median(rows, "bright_core_mid"))
        C.md_table(fh, ["padrão", "bright-core (Otsu)", "bright-core-mid (sem nível)"],
                   [[C.GROUP_PRETTY[g], f3(bc[g]), f3(bcm[g])]
                    for g in C.GROUP_ORDER])

        fh.write("\n## Por ensaio (transparência)\n\n")
        rr = sorted(rows, key=lambda r: -(r["bright_core"]
                    if not np.isnan(r["bright_core"]) else -1))
        C.md_table(fh,
            ["scenario", "pattern", "G", "q", "bright-core", "flip-rate",
             "α cordal", "sat %", "flag"],
            [[r["scenario"], r["pattern"], r["G"], r["q"],
              f"{r['bright_core']:.3f}", f"{r['flip_rate']:.4f}",
              f"{r['alpha_chordal']:.3f}", f"{r['sat_frac']*100:.0f}", r["flag"]]
             for r in rr])

        fh.write("\n## Ressalvas honestas\n\n")
        fh.write("- **Ordenação, não magnitude.** bright-core ~0,35 é cordal em "
                 "pixels (proxy óptico RELATIVO não-calibrado); longe do areal "
                 "~0,9 da literatura. Sem px→mm nem segunda vista, o vazio areal "
                 "verdadeiro NÃO é recuperável destas imagens.\n")
        fh.write("- **dispersed bubbles (n=2):** sai alto ESPÚRIO (centro "
                 "coincidentemente claro — a regra 'centro claro' precisa de porta "
                 "de contraste); anotado no mapa, fora da estatística.\n")
        fh.write("- **Sobreposição e n pequeno** (anular 8 × churn 10); plug "
                 "empata com anular. A separação é de tendência central.\n")
        fh.write("- **Eixo Y (flip-rate)** é ESTABILIDADE, não magnitude: põe "
                 "churn no topo (mais caótico), consistente com o anular embaixo "
                 "(estável).\n")
        fh.write("- **Significância é de DIREÇÃO, não de p forte** (re-verificação "
                 "adversarial independente, `holds=true`): a ordem anular>churn "
                 "NUNCA inverte sob perturbação (9 variantes de janela/limiar/"
                 "polaridade) e sobrevive ao q-match (o q-confound joga CONTRA — "
                 "anular é amostrado em q maior mas bright-core correla "
                 "NEGATIVAMENTE com q em ambos os padrões, r≈−0,4..−0,5). MAS o "
                 "p<0,05 é FRÁGIL nas bordas: p=0,022 na janela 10:70, p≈0,002 em "
                 "12:68–15:65, mas sobe a p=0,11 na janela larga 8:72 (dilui nas "
                 "paredes saturadas) e p=0,24 com Otsu+10. Reportar como efeito "
                 "borderline p≈0,02 COM a varredura de perturbação, não como p "
                 "forte. Sem correção de múltiplas comparações.\n")
    print(f"\ntabela: {out}")


if __name__ == "__main__":
    main()
