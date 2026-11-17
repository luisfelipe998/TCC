"""Fluxograma da seção 3.1, com operações separadas em tópicos."""
from pathlib import Path
import hashlib
import json
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/tcc-ilustracoes-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
EXECUTION = ROOT / "src2/v18_montagem_18_71"
WIDTH, HEIGHT, FONT = 16, 22.2, 9.5
INK = "#243842"
COLORS = {"input": "#eef1f3", "individual": "#e1f0f2",
          "pair": "#fff0da", "optical": "#eee7f3", "output": "#eef1f3"}
NODES = []


def node(code, title, items, group):
    NODES.append(dict(id=code, title=title, items=items, group=group))


node("1", "1. Seleção dos ensaios", [
    "26 aquisições de cinco padrões, com rótulos manuais do banco de imagens.",
], "input")
node("2", "2. Reconstrução dos panoramas", [
    "Comparação dos perfis de intensidade e estimativa dos deslocamentos.",
    "Recorte e junção das tiras na ordem de aquisição.",
], "input")
node("3", "3. Entradas de análise (linhas de índices 18 a 71)", [
    r"Sinais espaciais de intensidade $I_i[n]$: média por coluna para os caminhos A e B.",
    "Imagens dos panoramas para a segmentação no caminho C.",
], "input")

node("A1", "A1. Preparo individual", [
    "Inversão da intensidade.",
    "Remoção da tendência\nlinear.",
    "Padronização.",
], "individual")
node("B1", "B1. Preparo dos pares", [
    "Corte no menor\ncomprimento.",
    "Preparo após o corte,\nconforme A1.",
], "pair")
node("C1", "C1. Brilho e referência", [
    "Normalização do brilho\npor coluna.",
    "Referência líquida dos\ntrês ensaios de plug.",
], "optical")

node("A2", "A2. Correlação e CWT", [
    "Autocorrelação estimada\n" + r"$\widehat R_{xx}[m]$.",
    "CWT de Morlet sobre\ntodos os atrasos.",
], "individual")
node("B2", "B2. CWT de Morlet", [
    "Autocorrelações\n" + r"$\widehat R_{xx}[m]$ e $\widehat R_{yy}[m]$.",
    "Correlação cruzada\n" + r"estimada $\widehat R_{xy}[m]$.",
], "pair")
node("C2", "C2. Detecção de sombras", [
    "Suavização da imagem\ne da referência líquida.",
    "Detecção pela razão\nde intensidades.",
], "optical")

node("A3", "A3. Extração do descritor", [
    "Perfil de potência no\n" + r"atraso zero $P_{ik}$.",
    "200 valores, um por\nescala, após a CWT.",
], "individual")
node("B3", "B3. Mapas e curvas", [
    "Potências alternativas\n" + r"$\widetilde P_x(a,b)$ e $\widetilde P_y(a,b)$.",
    "Espectro cruzado\n" + r"$\widetilde C_{xy}(a,b)$.",
    "Perfis de potência\n" + r"média $M_{ik}$.",
], "pair")
node("C3", "C3. Segmentação", [
    "Fechamento de contornos.",
    "Preenchimento das áreas\npara formar a máscara\nde gás inferido.",
], "optical")

node("A4", "A4. Forma dos perfis", [
    "Coeficiente de Pearson\n" + r"$r_{ij}$ dos perfis completos.",
    "325 pares distintos.",
], "individual")
node("B4", "B4. Resumos de cada par", [
    "Coerência alternativa\n" + r"$\widetilde\gamma_{xy}(a,b)$.",
    "Mapa completo: mínimo,\nmédia, mediana, máximo.",
], "pair")
node("C4", "C4. Medidas por ensaio", [
    "Fração de vazio estimada\n" + r"$\widehat\alpha_{\mathrm{opt}}$: seções circulares.",
    "Densidade de transições\n" + r"$F$: alternância horizontal.",
], "optical")

node("4", "4. Comparação por grupo", [
    "Coeficiente de Pearson\n" + r"$r_{ij}$ e médias dos mapas" +
    "\nde coerência alternativa\n" + r"$\widetilde\gamma_{xy}(a,b)$.",
    "Resumos e distribuições\ndentro e entre padrões.",
], "output")
node("5", "5. Seleção dos exemplos", [
    "Quatro pares intrapadrão:\nmaior coeficiente de\n" +
    r"Pearson $r_{ij}$ em cada" + "\ngrupo com ≥2 ensaios.",
    "Um par entre padrões:\nmaior valor pontual da\ncoerência alternativa\n" +
    r"$\widetilde\gamma_{xy}(a,b)$.",
], "output")
node("C5", "C5. Critério de saturação", [
    "Até 20% de pixels com\nintensidade 255 no\npanorama original.",
    "25 ensaios nos resumos.",
], "optical")
node("C6", "C6. Resumos e gráficos", [
    "Médias, medianas e\ngráficos por padrão.",
    "Contornos (≥3 ensaios)\ne suas sobreposições.",
], "optical")
node("6", "6. Interpretação conjunta", [
    "Comparação descritiva dos resultados espectrais e ópticos.",
    "Discussão das semelhanças, diferenças e limitações do método.",
], "output")

# As saídas dos pares alimentam tanto os resumos quanto a seleção dos exemplos.
EDGES = [("1", "2"), ("2", "3"), ("3", "A1"), ("3", "B1"), ("3", "C1"),
         ("A1", "A2"), ("A2", "A3"), ("A3", "A4"),
         ("B1", "B2"), ("B2", "B3"), ("B3", "B4"),
         ("C1", "C2"), ("C2", "C3"), ("C3", "C4"),
         ("C4", "C5"), ("C5", "C6"),
         ("A4", "J"), ("B4", "J"), ("J", "4"), ("J", "5"),
         ("4", "6"), ("5", "6"), ("C6", "6")]
SOURCES = ["METODOLOGIA.md", "config/ensaios.json", "scripts/reconstruir.py",
           "scripts/wavelet_core.py", "scripts/resumir_pares.py",
           "scripts/optical_model.py", "scripts/process_optical.py",
           "04_analise_simplificada/scripts/analisar.py"]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": FONT,
                         "mathtext.fontset": "dejavusans",
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig = plt.figure(figsize=(WIDTH/2.54, HEIGHT/2.54))
    ax = fig.add_axes([0, 0, 1, 1], xlim=(0, WIDTH), ylim=(0, HEIGHT))
    ax.set_axis_off()
    by_id = {n["id"]: n for n in NODES}
    texts = {}
    for n in NODES:
        assert n["items"] and not any("->" in s or "→" in s for s in n["items"])
        title = ax.text(0, 0, n["title"], ha="center", va="top",
                        fontweight="bold", color=INK, zorder=3)
        items = [ax.text(0, 0, s, ha="left", va="top", linespacing=1.10,
                         color=INK, zorder=3) for s in n["items"]]
        bullets = [ax.text(0, 0, "•", ha="left", va="top", color=INK, zorder=3)
                   for _ in items]
        texts[n["id"]] = (title, items, bullets)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    px_cm = fig.dpi/2.54
    # Cada tópico tem marcador próprio; as continuações compartilham o recuo.
    title_gap, item_gap, padding = 4, 3, 7
    for n in NODES:
        title, items, _ = texts[n["id"]]
        heights = [t.get_window_extent(renderer).height for t in items]
        content_h = title.get_window_extent(renderer).height + title_gap
        content_h += sum(heights) + item_gap*(len(items)-1)
        n["min_height"] = (content_h+padding)/px_cm
        n["content_px"] = content_h

    def place(code, x, top, width, height=None):
        n = by_id[code]
        n.update(x=x, y=top-(height or n["min_height"]), width=width,
                 height=height or n["min_height"])

    gap = .245
    top = HEIGHT-.12
    for code in ("1", "2", "3"):
        place(code, .12, top, 15.76)
        top = by_id[code]["y"]-gap
    heading_y = top-.06
    top -= .43
    for row in (("A1", "B1", "C1"), ("A2", "B2", "C2"),
                ("A3", "B3", "C3"), ("A4", "B4", "C4")):
        h = max(by_id[code]["min_height"] for code in row)
        for x, code in zip((.12, 5.50, 10.88), row):
            place(code, x, top, 5.0, h)
        top -= h+gap
    join_y = top-.02
    outputs_top = top-.31
    optical_top = top
    optical_h = by_id["C5"]["min_height"]+gap+by_id["C6"]["min_height"]
    lower_h = max(optical_h, .31+by_id["4"]["min_height"],
                  .31+by_id["5"]["min_height"])
    bottom = optical_top-lower_h
    for x, code in ((.12, "4"), (5.50, "5")):
        place(code, x, outputs_top, 5.0, outputs_top-bottom)
    place("C5", 10.88, optical_top, 5.0)
    place("C6", 10.88, by_id["C5"]["y"]-gap, 5.0,
          by_id["C5"]["y"]-gap-bottom)
    place("6", .12, bottom-.42, 15.76)
    assert by_id["6"]["y"] >= .10, {
        "altura_necessaria_cm": HEIGHT-by_id["6"]["y"]+.12,
        "alturas": {n["id"]: round(n["height"], 3) for n in NODES}}

    for n in NODES:
        x, y, w, h = (n[k] for k in ("x", "y", "width", "height"))
        ax.add_patch(FancyBboxPatch((x, y), w, h,
                     boxstyle="round,pad=0.018,rounding_size=0.06",
                     facecolor=COLORS[n["group"]], edgecolor=INK,
                     linewidth=.65, zorder=2))
        title, items, bullets = texts[n["id"]]
        title_h = title.get_window_extent(renderer).height/px_cm
        cursor = y+h-padding/(2*px_cm)
        title.set_position((x+w/2, cursor))
        cursor -= title_h+title_gap/px_cm
        for item, bullet in zip(items, bullets):
            item.set_position((x+.43, cursor))
            bullet.set_position((x+.16, cursor))
            cursor -= (item.get_window_extent(renderer).height+item_gap)/px_cm

    headings = []
    for x, title in ((2.62, "A · Perfis completos"),
                     (8.00, "B · Mapas dos pares"),
                     (13.38, "C · Análise óptica")):
        headings.append(ax.text(x, heading_y, title, ha="center", va="center",
                                fontweight="bold", color=INK, zorder=3))

    def arrow(start, end):
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>",
                     mutation_scale=7, linewidth=.8, color=INK,
                     shrinkA=0, shrinkB=0, zorder=1))

    for source, target in EDGES:
        if "J" in (source, target):
            continue
        s, t = by_id[source], by_id[target]
        start = (s["x"]+s["width"]/2, s["y"]-.035)
        end = (t["x"]+t["width"]/2, t["y"]+t["height"]+.045)
        if source == "3":
            x = t["x"]+t["width"]-.18
            start, end = (x, start[1]), (x, end[1])
        if target == "6":
            end = (start[0], end[1])
        arrow(start, end)
    ax.plot([2.62, 8.0], [join_y, join_y], color=INK, linewidth=.8, zorder=1)
    for source, target in (("A4", "4"), ("B4", "5")):
        s, t = by_id[source], by_id[target]
        x = s["x"]+s["width"]/2
        ax.plot([x, x], [s["y"]-.035, join_y], color=INK, linewidth=.8, zorder=1)
        ax.plot(x, join_y, "o", color=INK, markersize=2, zorder=1)
        arrow((x, join_y), (x, t["y"]+t["height"]+.045))

    fig.canvas.draw()
    checked, errors = [], []
    for n in NODES:
        box = ax.transData.transform_bbox(
            matplotlib.transforms.Bbox.from_bounds(n["x"], n["y"], n["width"], n["height"]))
        title, items, bullets = texts[n["id"]]
        for artist in (title, *items, *bullets):
            bounds = artist.get_window_extent(renderer)
            right_margin = 3 if artist is title else 4
            if not (bounds.x0 >= box.x0+1 and bounds.x1 <= box.x1-right_margin
                    and bounds.y0 >= box.y0+1 and bounds.y1 <= box.y1-1):
                errors.append((n["id"], artist.get_text(), "fora da caixa",
                               round(bounds.width, 1), round(box.width, 1)))
        for previous, current in zip((title, *items), items):
            if current.get_window_extent(renderer).y1+2 > previous.get_window_extent(renderer).y0:
                errors.append((n["id"], "sobreposição entre tópicos"))
        checked.append(n["id"])
    assert not errors, errors
    for ext in ("pdf", "png"):
        fig.savefig(HERE/f"fluxograma.{ext}", dpi=300, facecolor="white")
    plt.close(fig)
    metadata = {
        "tipo": "esquema do procedimento executado, sem novos dados experimentais",
        "secao": "3.1 — Etapas do processamento",
        "dimensoes_cm": [WIDTH, HEIGHT], "fonte_pt": FONT,
        "nos": NODES, "arestas": EDGES,
        "juncoes": {"J": "Resultados A4 e B4 fornecidos às etapas 4 e 5."},
        "fontes_execucao_sha256": {
            str((EXECUTION/name).relative_to(ROOT)): sha(EXECUTION/name)
            for name in SOURCES},
        "criterios": [
            "O bloco 3 reúne as entradas de análise; as diferenças de montagem permanecem descritas no texto.",
            "Todas as análises usam as linhas 18 a 71.",
            "Cada operação ou resultado do bloco possui um marcador; continuações têm recuo alinhado.",
            "Os conectores indicam o fluxo; os textos dos blocos não usam setas como separadores.",
            "Pearson usa os perfis das autocorrelações dos sinais completos, após a CWT.",
            "Nos pares, o corte antecede o preparo e a CWT é aplicada às correlações.",
            "As quatro estatísticas por mapa usam todas as células válidas.",
            "Os resumos entre pares usam um Pearson e uma média de coerência por par, com igual peso.",
            "A seleção dos exemplos usa resultados individuais dos pares, não os resumos de grupo.",
            "As 26 máscaras e medidas precedem a seleção dos 25 ensaios dos resumos ópticos.",
            "A interpretação conjunta é descritiva; os pares específicos são detalhados no texto.",
        ],
        "verificacao_texto_contido_nas_caixas": checked,
        "script_sha256": sha(Path(__file__)),
        "arquivos_sha256": {f"fluxograma.{e}": sha(HERE/f"fluxograma.{e}")
                            for e in ("pdf", "png")},
    }
    (HERE/"fluxograma_origem.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2)+"\n")
    print(f"Fluxograma gerado: {len(NODES)} blocos; "
          f"{sum(len(n['items']) for n in NODES)} tópicos; fonte {FONT} pt.")


if __name__ == "__main__":
    main()
