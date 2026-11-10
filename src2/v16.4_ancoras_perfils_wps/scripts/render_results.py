"""Generate scientific figures and local reports from the v16.4 result tables."""
from __future__ import annotations

import base64
import hashlib
import html
import json
import os
import sys
from pathlib import Path

sys.dont_write_bytecode = True
os.environ.setdefault("MPLCONFIGDIR", "/tmp/v16_4_ancoras_wps_mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
OUT, FIG = ROOT / "results", ROOT / "figures"
PATTERNS = ["plug", "slug", "churn", "wavy annular"]
LABELS = {"plug": "Bolhas alongadas (plug)", "slug": "Golfadas (slug)",
          "churn": "Agitado (churn)", "wavy annular": "Anular ondulado"}
COLORS = dict(zip(PATTERNS, ["#29815d", "#8854a3", "#d75c37", "#2676a6"]))
CONDITIONS = ["ancora_2048", "inicio_2048", "centro_2048", "fim_2048", "panorama_completo"]
NAMES = {"ancora_2048": "Âncora · 2.048 px", "inicio_2048": "Início · 2.048 px",
         "centro_2048": "Centro · 2.048 px", "fim_2048": "Fim · 2.048 px",
         "panorama_completo": "Panorama completo", "ancora_menos16": "Âncora −16 px",
         "ancora_mais16": "Âncora +16 px"}
MAIN = "ancora_2048__hellinger_principal"
FULL = "panorama_completo__hellinger_principal"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pct(value):
    return f"{100*value:.1f}%".replace(".", ",")


def num(value, digits=3):
    return f"{value:.{digits}f}".replace(".", ",")


def key(condition, domain="principal"):
    return condition + "__hellinger_" + domain


def table_md(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"]*len(headers)) + " |"] +
                     ["| " + " | ".join(map(str, row)) + " |" for row in rows])


def save(fig, name):
    path = FIG / f"{name}.png"
    fig.savefig(path, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=.12)
    plt.close(fig)
    return path


def make_figures(summary, pairs, by_pattern, index, profiles, archive):
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
        "legend.fontsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.spines.top": False, "axes.spines.right": False, "font.family": "DejaVu Sans"})
    figures = []
    # Fixed illustrative choice: first ID of plug in the catalogue order.
    example = index[index.padrao == "plug"].sort_values("ensaio").iloc[0]
    image_path = ROOT.parent / "data_processed" / example.ensaio / "__panorama_v11.png"
    assert digest(image_path) == example.panorama_sha256
    pano = np.asarray(Image.open(image_path).convert("L"))
    raw = pano[18:72].mean(axis=0)
    lo, hi, anchor = int(example.inicio_px), int(example.fim_exclusivo_px), int(example.ancora_px)
    fig, axes = plt.subplots(3, 1, figsize=(6.3, 7.7), layout="constrained",
                             gridspec_kw={"height_ratios": [.9, 1, 1.2]})
    axes[0].imshow(pano, cmap="gray", vmin=0, vmax=255, aspect="auto", extent=[0, len(raw), pano.shape[0], 0])
    axes[0].axvspan(lo, hi, color="#df8236", alpha=.25)
    axes[0].axvline(anchor, color="#d5422f", lw=1)
    axes[0].axhline(18, color="#289bd0", lw=.8)
    axes[0].axhline(72, color="#289bd0", lw=.8)
    axes[0].set(title=f"1. Panorama e recorte — {example.ensaio} (plug)", xlabel="Coluna do panorama (px)", ylabel="Linha")
    axes[1].plot(np.arange(lo, hi)-anchor, raw[lo:hi], color="#37475c", lw=.8)
    axes[1].axvline(0, color="#d5422f", ls="--", lw=.8)
    axes[1].set(title="2. Média das linhas 18–71: entrada antes do preparo",
                xlabel="Posição relativa à âncora (px)", ylabel="Intensidade média")
    for condition, color in (("ancora_2048", "#bd4e2b"), ("panorama_completo", "#246c99")):
        p = profiles[(profiles.ensaio == example.ensaio) & (profiles.condicao == condition) & (profiles.dominio == "principal")]
        axes[2].semilogx(p.lambda_px, p.densidade_log, label=NAMES[condition], color=color, lw=1.5)
    axes[2].set(title="3. Perfil WPS no lag zero, com área igual a um",
                xlabel="Comprimento de onda λ (px; escala log)", ylabel="Densidade por log λ")
    axes[2].set_xticks([50, 100, 200, 500, 1024], ["50", "100", "200", "500", "1.024"])
    axes[2].legend()
    save(fig, "01_entrada_e_perfil")
    figures.append(("01_entrada_e_perfil", "Da âncora ao perfil", "O sombreado marca os 2.048 px usados. A linha vertical marca a âncora e as horizontais delimitam a faixa de linhas. A WPS descreve todo o sinal recortado; λ não localiza uma estrutura no panorama. O panorama foi ajustado em proporção para leitura. Exemplo fixado pelo primeiro ID de plug, sem seleção pelo resultado."))

    fig, axes = plt.subplots(2, 2, figsize=(6.3, 5.7), sharex=True, sharey=True, layout="constrained")
    anchored = profiles[(profiles.condicao == "ancora_2048") & (profiles.dominio == "principal")]
    for ax, pattern in zip(axes.flat, PATTERNS):
        part = anchored[anchored.padrao == pattern]
        for _, group in part.groupby("ensaio"):
            ax.semilogx(group.lambda_px, group.densidade_log, color=COLORS[pattern], alpha=.4, lw=.8)
        mean = part.groupby("lambda_px").densidade_log.mean()
        ax.semilogx(mean.index, mean.values, color=COLORS[pattern], lw=2, label="Média visual")
        ax.set_title(f"{LABELS[pattern]} · n={part.ensaio.nunique()}")
        ax.set_xticks([50, 100, 250, 500, 1024], ["50", "100", "250", "500", "1.024"])
        ax.grid(alpha=.17)
    axes[0, 0].set_ylabel("Densidade por log λ")
    axes[1, 0].set_ylabel("Densidade por log λ")
    for ax in axes[1]:
        ax.set_xlabel("Comprimento de onda λ (px)")
    save(fig, "02_perfis_ancorados")
    figures.append(("02_perfis_ancorados", "Um perfil por aquisição", "Linhas finas: todos os ensaios. Linha grossa: média das densidades apenas para leitura. Cada distância foi calculada entre perfis individuais; essa média não é um modelo de classe. A escala vertical comum evidencia a variabilidade interna."))

    fig, ax = plt.subplots(figsize=(6.3, 4.8), layout="constrained")
    methods = [key(c) for c in CONDITIONS] + ["v163_pearson_ancora", "v163_wco_central"]
    names = [NAMES[c] for c in CONDITIONS] + ["Pearson · âncora", "WCO central · âncora"]
    y = np.arange(len(methods))
    for offset, column, label, color in ((-.18, "A_pares", "Peso igual por par", "#235f83"),
                                        (.18, "A_balanceada", "Peso igual por grupo", "#d57c32")):
        values = summary.loc[methods, column].to_numpy()*100
        ax.barh(y+offset, values, height=.32, color=color, label=label)
        for yy, value in zip(y+offset, values):
            ax.text(value+.8, yy, num(value, 1), va="center", fontsize=8)
    ax.axvline(50, color="#6d6d6d", ls="--", lw=.9)
    ax.set(yticks=y, yticklabels=names, xlim=(0, 100), xlabel="Concordância de ordenação A (%)")
    ax.invert_yaxis()
    ax.legend(loc="lower right")
    save(fig, "03_comparacao_controles")
    figures.append(("03_comparacao_controles", "Âncoras e controles no mesmo conjunto", "Todos usam os mesmos 25 ensaios e 50–1.024 px. As cinco primeiras linhas usam Hellinger; Pearson e WCO são referências da v16.3. Os valores representam A, tornando comparável a ordenação. A linha de 50% indica ausência de preferência; não é teste estatístico nem acurácia."))

    fig, axes = plt.subplots(2, 1, figsize=(6.3, 5.7), sharex=True, sharey=True, layout="constrained")
    limit = max(pairs[pairs.metodo.isin([MAIN, FULL])].valor.max(), .5)
    bins = np.arange(0, np.ceil(limit/.05)*.05+.051, .05)
    for ax, method, title in zip(axes, [MAIN, FULL], ["Âncoras · 2.048 px", "Panoramas completos"]):
        part = pairs[pairs.metodo == method]
        for same, label, color in ((True, "Mesmo padrão · 87 pares", "#287a9b"), (False, "Padrões diferentes · 213 pares", "#cb6a3e")):
            values = part[part.mesmo_padrao == same].valor
            ax.hist(values, bins=bins, weights=np.full(len(values), 100/len(values)),
                    alpha=.4, color=color, label=label, edgecolor=color)
            ax.axvline(values.median(), color=color, ls="--", lw=1)
        ax.set(title=title, ylabel="Pares do grupo (%)")
        ax.legend()
    axes[-1].set_xlabel("Distância Hellinger H — menor significa mais semelhante")
    save(fig, "04_distribuicoes_distancias")
    figures.append(("04_distribuicoes_distancias", "Sobreposição das distâncias", "Cada histograma soma 100% dentro de seu grupo; as linhas tracejadas marcam medianas. As âncoras produziram maior sobreposição. O domínio matemático de H vai de 0 a 1; a figura mostra a faixa observada. Os pares compartilham aquisições."))

    fig, axes = plt.subplots(1, 2, figsize=(6.3, 4.6), sharey=True, layout="constrained")
    y = np.arange(len(CONDITIONS))
    for ax, column, title in zip(axes, ["A_pares", "A_balanceada"], ["Peso igual por par", "Peso igual por grupo"]):
        for offset, domain, label, color in ((-.16, "principal", "50–1.024 px", "#235f83"),
                                            (.16, "conservador", "50–348,7 px", "#499d86")):
            values = np.asarray([summary.loc[key(c, domain), column] for c in CONDITIONS])*100
            ax.barh(y+offset, values, height=.29, label=label, color=color)
            for yy, value in zip(y+offset, values):
                ax.text(value+.8, yy, num(value, 1), va="center", fontsize=7.5)
        ax.axvline(50, color="#777777", ls="--", lw=.8)
        ax.set(xlim=(0, 100), title=title, xlabel="A (%)", yticks=y)
    axes[0].set_yticklabels([NAMES[c].replace(" · 2.048 px", "") for c in CONDITIONS])
    axes[0].invert_yaxis()
    axes[1].legend(loc="upper center", bbox_to_anchor=(.5, -.17), ncol=2, fontsize=7)
    save(fig, "05_sensibilidade_escalas")
    figures.append(("05_sensibilidade_escalas", "A faixa de escalas importa", "A faixa conservadora segue 3a≤N/2 e contém 129 pontos da grade original; as áreas foram normalizadas novamente. As âncoras melhoraram nessa faixa, mas não superaram o panorama completo nem o fim. As duas faixas foram preservadas, sem escolher a melhor após os resultados."))

    fig, axes = plt.subplots(2, 1, figsize=(6.3, 5.6), sharex=True, sharey=True, layout="constrained")
    x = np.arange(4)
    for ax, column, title in zip(axes, ["A_pares", "A_balanceada"], ["Peso igual por par", "Controles balanceados por padrão"]):
        for offset, method, label, color in ((-.17, MAIN, "Âncoras", "#bd4e2b"), (.17, FULL, "Panoramas completos", "#246c99")):
            values = by_pattern[by_pattern.metodo == method].set_index("padrao").loc[PATTERNS, column].to_numpy()*100
            ax.bar(x+offset, values, width=.3, color=color, label=label)
            for xx, value in zip(x+offset, values):
                ax.text(xx, value+1, num(value, 1), ha="center", fontsize=8)
        ax.axhline(50, color="#777777", ls="--", lw=.8)
        ax.set(ylim=(0, 100), ylabel="A por padrão (%)", title=title)
    axes[0].legend(loc="upper right", fontsize=7.5)
    axes[-1].set_xticks(x, ["Plug\n3 pares intra", "Slug\n3 pares intra", "Churn\n45 pares intra", "Anular ondulado\n36 pares intra"])
    save(fig, "06_resultados_por_padrao")
    figures.append(("06_resultados_por_padrao", "A perda aparece nos quatro padrões", "Cada padrão foi avaliado contra os pares inter que o envolvem. A segunda ponderação dá o mesmo peso às três outras classes. Anular ondulado ficou abaixo de 50% nas âncoras. Plug e slug têm somente três aquisições cada."))
    return figures


def render_documents(summary, pairs, by_pattern, deletion, effects, checks, figures):
    main, full = summary.loc[MAIN], summary.loc[FULL]
    conservative = summary.loc[key("ancora_2048", "conservador")]
    main_del = deletion[deletion.metodo == MAIN]
    title = "Âncoras + perfis WPS: combinação coerente, sem melhora geral neste teste"
    lead = ("Faz sentido combinar as duas etapas: a v16.3 seleciona o trecho e a v16.1 compara a forma de sua WPS. "
            "A análise executada, porém, não sustentou substituir os panoramas completos por uma única janela ancorada de 2.048 px. "
            f"Na faixa de 50–1.024 px, A foi {pct(main.A_pares)} nas âncoras e {pct(full.A_pares)} nos panoramas completos. "
            "Os controles de início, centro e fim com o mesmo comprimento também superaram as âncoras.")
    blocks = [("h1", title), ("p", lead), ("p", "A mede concordância de ordenação entre pares. Não é acurácia, probabilidade de classificação correta ou significância estatística.")]

    def section(name):
        blocks.append(("h2", name))

    def paragraph(text):
        blocks.append(("p", text))

    def table(headers, rows):
        blocks.append(("table", headers, rows))

    def figure(index):
        name, caption, description = figures[index]
        blocks.append(("figure", name, caption, description))

    section("Como a proposta foi aplicada")
    paragraph("O caminho foi: panorama → âncora congelada → recorte de 2.048 px → média por coluna das linhas 18–71 → preparo do sinal → autocorrelação → WPS no lag zero → área unitária em log λ → distância de Hellinger. A âncora define qual conteúdo entra no sinal. A normalização e a distância são as mesmas da v16.1.")
    paragraph("Foram preservados 25 ensaios: 3 plug, 3 slug, 10 churn e 9 wavy annular. Cada um forneceu um perfil individual, independente do parceiro, formando 87 pares intra e 213 inter. A aquisição de bolhas dispersas da v16.1 não pertence ao catálogo de âncoras; por isso, foi excluída de todos os controles diretos. Os centros, regras ópticas e aprovações morfológicas vieram da v16.3, sem nova seleção.")
    paragraph("Foi usada a WPS alternativa |CWT(Rxx)| no lag zero, sem elevar novamente ao quadrado. A curva de potência média sobre lags da v16.3 não é essa representação. Tampouco se aplicou outra wavelet sobre um perfil já calculado. O núcleo Morlet das duas versões é idêntico em bytes.")
    figure(0)
    figure(1)
    section("Comparação principal em condições comuns")
    rows = [[NAMES[c], pct(summary.loc[key(c), "A_pares"]), pct(summary.loc[key(c), "A_balanceada"])] for c in CONDITIONS]
    table(["Entrada do perfil Hellinger", "A: peso igual por par", "A: peso igual por grupo"], rows)
    paragraph("Todas essas linhas usam os mesmos ensaios, extração, preparo e 200 comprimentos de onda entre 50 e 1.024 px. Início, centro e fim foram cortados antes do preparo e têm os mesmos 2.048 px das âncoras. Assim, a vantagem dos controles de mesmo tamanho não pode ser explicada apenas por terem mais amostras. O panorama completo modifica comprimento e conteúdo simultaneamente; ele não isola o efeito do comprimento.")
    delta = full.A_pares-main.A_pares
    paragraph(f"A diferença âncora menos panorama completo foi −{num(100*delta, 1)} pontos percentuais por par e {num(100*(main.A_balanceada-full.A_balanceada), 1)} pontos na ponderação balanceada. H médio foi {num(main.media_intra)} intra e {num(main.media_entre)} inter nas âncoras; nos panoramas completos, {num(full.media_intra)} e {num(full.media_entre)}, respectivamente. Menor H significa maior semelhança; essas médias não são limiares de classificação.")
    figure(2)
    figure(3)
    counts = [f"{int(v):,}".replace(",", ".") for v in [main.ordens_favoraveis, main.ordens_contrarias, main.comparacoes]]
    paragraph(f"A principal veio de {counts[0]} ordenações favoráveis, {counts[1]} contrárias e {int(main.empates)} empates, entre 87×213={counts[2]} comparações. Esses números reutilizam apenas 25 aquisições. Na ponderação balanceada, cada um dos quatro grupos intra e cada uma das seis combinações inter recebe peso igual; não são novas observações.")

    section("O que os controles de estabilidade mostraram")
    rows = [[NAMES[c], pct(summary.loc[key(c, "conservador"), "A_pares"]), pct(summary.loc[key(c, "conservador"), "A_balanceada"])] for c in CONDITIONS]
    table(["Entrada — faixa conservadora", "A por par", "A balanceada"], rows)
    paragraph(f"O suporte conservador 3a≤N/2 reteve {checks['escalas_conservadoras']} pontos, até λ={num(checks['lambda_conservadora_max_px'], 2)} px. A das âncoras subiu para {pct(conservative.A_pares)} / {pct(conservative.A_balanceada)}. Nessa faixa, elas superaram início e centro, mas ficaram abaixo do fim e do panorama completo. Há informação útil e sensibilidade à banda; não há ganho geral da seleção por uma única âncora. A faixa conservadora foi definida por suporte, antes de calcular esses controles, e não substituiu a principal.")
    figure(4)
    table(["Perturbação da âncora", "A por par", "A balanceada"],
          [[NAMES[c], pct(summary.loc[key(c), "A_pares"]), pct(summary.loc[key(c), "A_balanceada"])] for c in ["ancora_menos16", "ancora_2048", "ancora_mais16"]])
    paragraph("Deslocar os centros em ±16 px alterou pouco o resultado na faixa principal. Essa checagem testa um pequeno erro de posição; não representa selecionar outro evento nem repetir a seleção completa. Os recortes perturbados não receberam uma nova aprovação morfológica. Na sensibilidade de métrica, o cosseno nas âncoras produziu " + pct(summary.loc["ancora_2048__cosseno_principal", "A_pares"]) + " por par, também sem vantagem sobre Hellinger.")
    paragraph(f"Retirando uma aquisição por vez, A das âncoras variou de {pct(main_del.A_pares.min())} a {pct(main_del.A_pares.max())}; a balanceada variou de {pct(main_del.A_balanceada.min())} a {pct(main_del.A_balanceada.max())}. Em todas as 25 retiradas pareadas, as âncoras permaneceram abaixo de início, centro, fim e panorama completo na faixa principal, nas duas ponderações. Essas são análises de influência, não intervalos de confiança nem validação preditiva.")

    section("Resultados por padrão")
    bp = by_pattern.set_index(["metodo", "padrao"])
    table(["Padrão", "Âncora: A por par", "Completo: A por par", "Âncora: A balanceada", "Completo: A balanceada"],
          [[LABELS[p], pct(bp.loc[(MAIN, p), "A_pares"]), pct(bp.loc[(FULL, p), "A_pares"]),
            pct(bp.loc[(MAIN, p), "A_balanceada"]), pct(bp.loc[(FULL, p), "A_balanceada"])] for p in PATTERNS])
    paragraph("A desvantagem perante os panoramas completos apareceu nos quatro padrões. Anular ondulado teve o resultado mais desfavorável das âncoras. A maior A de slug não é evidência de vantagem da âncora: o controle completo do mesmo padrão também foi maior. Por padrão, os controles inter envolvem somente aquela classe; esses números não são a decomposição simples da A global.")
    figure(5)

    section("Relação com as análises anteriores")
    table(["Comparação nos mesmos 25 ensaios", "A por par", "A balanceada"],
          [["Hellinger nas âncoras · 50–1.024 px", pct(main.A_pares), pct(main.A_balanceada)],
           ["Pearson residual nas âncoras · 50–1.024 px", pct(summary.loc["v163_pearson_ancora", "A_pares"]), pct(summary.loc["v163_pearson_ancora", "A_balanceada"])],
           ["WCO média central nas âncoras · 50–1.024 px", pct(summary.loc["v163_wco_central", "A_pares"]), pct(summary.loc["v163_wco_central", "A_balanceada"])],
           ["v16.1 original, restrita aos 25 · 50–2.295,94 px", pct(summary.loc["v161_principal_25", "A_pares"]), pct(summary.loc["v161_principal_25", "A_balanceada"])]])
    paragraph("Hellinger ficou ligeiramente acima de Pearson e WCO no agregado das âncoras, mas a direção dessa vantagem não persistiu em todas as retiradas de aquisição. Isso não constitui melhora robusta. A WCO usada aqui é a média central da v16.3; não é o escore de máximos que obteve 21,5% na comparação antiga da v16.1.")
    paragraph("Os 76,5% / 81,6% publicados na v16.1 usavam 26 ensaios, incluindo uma aquisição de bolhas dispersas. Restrita aos mesmos 25, aquela análise produz 74,2% / 71,4%, mas conserva uma faixa mais ampla. O controle direto desta nova análise é o panorama completo recalculado/conferido em 50–1.024 px: 76,3% / 74,6%. Não atribuímos à âncora mudanças de amostra ou de faixa.")

    section("Interpretação e caminho recomendado")
    paragraph("A combinação foi tecnicamente válida e relevante para testar a hipótese. Os resultados atuais desfavorecem uma única âncora de 2.048 px como substituta do panorama completo para esse descritor. O panorama completo continua sendo a referência mais promissora entre as entradas avaliadas na faixa principal. A v16.4 preserva uma análise exploratória com resultado desfavorável, útil para justificar essa decisão.")
    paragraph("Uma explicação plausível é que um evento opticamente bem delimitado ou um vale proeminente não representa toda a diversidade do escoamento. A WPS da autocorrelação resume todo o recorte e perde informações de localização e fase; colocar uma estrutura no centro não torna o perfil específico daquela estrutura. Os controles sustentam a importância do conteúdo escolhido, mas não demonstram essa explicação causalmente.")
    paragraph("Se a investigação continuar, o próximo passo mais informativo é comparar conjuntos de múltiplos eventos/janelas por aquisição, com a mesma cobertura e regra para todos os padrões, e avaliar a estabilidade entre esses eventos. Esse experimento ainda não foi executado. Uma avaliação de classificação futura precisará de seleção aplicável sem conhecer a classe e de separação por aquisição/condição experimental; os pares ou janelas da mesma aquisição não poderão atravessar treino e teste.")
    paragraph("A seleção atual usa conhecimento do padrão: plug/slug dependem de eventos com morfologia aprovada, enquanto churn/anular usam vales. O presente estudo é descritivo, condicionado a esses rótulos. Os panoramas têm comprimentos desiguais, costuras e condições operacionais confundidas com os padrões. K=1, poucas réplicas em plug/slug e um evento por ensaio limitam qualquer generalização. λ está em pixels do panorama, sem calibração para milímetros, hertz ou comprimento físico de bolha.")

    section("Verificação, arquivos e reprodução")
    paragraph(f"A verificação numérica passou: {checks['perfis_verificados']} perfis de {checks['condicoes']} condições, {checks['metodos']} representações/referências e 300 pares por método. Para os 150 recortes, a autocorrelação direta e a Morlet temporal foram confrontadas com o cálculo por FFT em todas as escalas; os 25 perfis completos temporais já auditados foram confrontados com a implementação da v16.1. A maior diferença de WPS foi {checks['maior_erro_wps']:.2e}; a de A por formulações independentes foi {checks['maior_erro_A']:.2e}.")
    paragraph(f"Foram congelados os hashes de {checks['arquivos_fontes_congelados']} arquivos consumidos. Os {checks['arquivos_anteriores_preservados']} arquivos das duas análises anteriores preservaram nomes, tamanhos e mtimes durante a execução, e os hashes das entradas foram reconferidos. Não foram refeitas a revisão morfológica, a montagem dos panoramas ou a auditoria dos 300 mapas WCO; suas verificações anteriores foram reutilizadas. Os mapas cruzados não foram recalculados.")
    table(["Arquivo", "Conteúdo"], [["PROTOCOLO.md / config.json", "Regras e momento de definição dos controles"],
        ["results/indice_ensaios.csv / recortes.csv", "Âncoras, fontes e coordenadas de cada janela"],
        ["results/perfis.npz / perfis_normalizados.csv", "WPS, massas, densidades e áreas"],
        ["results/pares.csv / concordancia.csv", "Todos os pares e todos os resultados"],
        ["results/por_padrao.csv / grupos_pares.csv", "Desfechos por padrão e combinação"],
        ["results/influencia_aquisicoes.csv / efeitos_pareados.csv", "Retiradas de aquisição e diferenças sob a mesma retirada"],
        ["results/verificacao.json / verificacao_perfis.csv", "Controles numéricos e erros observados"],
        ["results/entradas_congeladas.json / artefatos.json", "Procedência das fontes e manifesto da entrega"],
        ["figures/", "Seis figuras científicas PNG a 300 dpi"]])
    blocks.append(("code", "PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/v16_4_ancoras_wps_mpl python3 src2/v16.4_ancoras_perfils_wps/scripts/run_analysis.py\nPYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/v16_4_ancoras_wps_mpl python3 src2/v16.4_ancoras_perfils_wps/scripts/render_results.py\nPYTHONDONTWRITEBYTECODE=1 python3 src2/v16.4_ancoras_perfils_wps/scripts/verify_results.py"))
    paragraph("Executar a partir da raiz do projeto, com as dependências de requirements.txt. As pastas v16.1, v16.3 e data_processed devem permanecer disponíveis como irmãs. A segunda etapa apenas renderiza resultados existentes. As regras foram registradas após a triagem inicial dos perfis disponíveis; este estudo não é um pré-registro cego. Não foram alterados parâmetros para favorecer as âncoras.")
    section("Todas as representações e referências")
    table(["Identificador reproduzível", "A por par", "A balanceada"],
          [[method, pct(row.A_pares), pct(row.A_balanceada)] for method, row in summary.iterrows()])

    md, markup = [], []
    for block in blocks:
        kind, *content = block
        if kind in ("h1", "h2"):
            md.append(("# " if kind == "h1" else "## ") + content[0])
            markup.append(f"<{kind}>{html.escape(content[0])}</{kind}>")
        elif kind == "p":
            md.append(content[0])
            markup.append(f"<p>{html.escape(content[0])}</p>")
        elif kind == "table":
            headers, rows = content
            md.append(table_md(headers, rows))
            markup.append("<div class='table-wrap'><table><thead><tr>" + "".join(f"<th>{html.escape(str(v))}</th>" for v in headers) +
                "</tr></thead><tbody>" + "".join("<tr>" + "".join(f"<td>{html.escape(str(v))}</td>" for v in row) + "</tr>" for row in rows) + "</tbody></table></div>")
        elif kind == "figure":
            name, caption, description = content
            md.append(f"![{caption}](figures/{name}.png)\n\n{description}")
            encoded = base64.b64encode((FIG / f"{name}.png").read_bytes()).decode()
            markup.append(f"<figure><img src='data:image/png;base64,{encoded}' alt='{html.escape(caption)}'><figcaption><strong>{html.escape(caption)}.</strong> {html.escape(description)}</figcaption></figure>")
        elif kind == "code":
            md.append("```sh\n" + content[0] + "\n```")
            markup.append("<pre><code>" + html.escape(content[0]) + "</code></pre>")
    (ROOT / "RELATORIO.md").write_text("\n\n".join(md) + "\n", encoding="utf-8")
    css = """body{margin:0;background:#f2f4f6;color:#172631;font:17px/1.65 system-ui,sans-serif}main{max-width:990px;margin:30px auto;background:white;padding:32px 44px;border-radius:8px}h1{font-size:2em;line-height:1.2}h2{margin-top:2em;line-height:1.3;color:#245a77}p{max-width:86ch}.table-wrap{overflow:auto}table{border-collapse:collapse;width:100%;font-size:.87em;line-height:1.45}th,td{text-align:left;padding:10px;border-bottom:1px solid #d9e1e6}th{background:#eef4f6}figure{margin:32px 0}img{display:block;width:100%;max-width:740px;height:auto;margin:auto}figcaption{font-size:.9em;color:#425561;margin-top:12px}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#eef2f5;padding:16px;font-size:.75em}@media(max-width:680px){main{padding:18px;margin:0}body{font-size:16px}h1{font-size:1.65em}}@media print{body{background:white;font-size:10pt}main{padding:0;margin:0}figure{break-inside:avoid}h2{break-after:avoid}img{max-height:20cm;object-fit:contain}}"""
    (ROOT / "RESUMO_EXECUTIVO.html").write_text("<!doctype html><html lang='pt-BR'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'><title>" + html.escape(title) + "</title><style>" + css + "</style></head><body><main>" + "\n".join(markup) + "</main></body></html>\n", encoding="utf-8")
    readme = f"""# v16.4 — Âncoras como entrada dos perfis WPS

**A combinação faz sentido metodológico, mas uma única âncora de 2.048 px não melhorou os resultados neste conjunto.** A pasta foi criada para uma análise exploratória controlada; não representa uma nova abordagem validada.

Abra o [relatório visual, com imagens incorporadas](RESUMO_EXECUTIVO.html) ou o [relatório em Markdown](RELATORIO.md).

{table_md(['Entrada · Hellinger em 50–1.024 px', 'A por par', 'A balanceada'], [[NAMES[c], pct(summary.loc[key(c), 'A_pares']), pct(summary.loc[key(c), 'A_balanceada'])] for c in CONDITIONS])}

A é concordância de ordenação, **não acurácia**. Todos os controles usam os mesmos 25 ensaios, com 87 pares intra e 213 inter. Os controles de início/centro/fim têm o mesmo comprimento das âncoras. A desvantagem principal persistiu nas 25 retiradas de aquisição, nas duas ponderações. Na faixa conservadora de 50–348,69 px, as âncoras chegaram a {pct(conservative.A_pares)} por par, ainda abaixo do panorama completo e da janela final na mesma faixa.

A proposta foi aplicada à média por coluna do **recorte ancorado**, seguida de WPS no lag zero, normalização de área e Hellinger da v16.1. Não se utilizou a potência média sobre lags. A seleção da âncora depende do padrão conhecido; o estudo não valida classificação automática.

## Reprodução

Na raiz do projeto, com [requirements.txt](requirements.txt):

```sh
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/v16_4_ancoras_wps_mpl python3 src2/v16.4_ancoras_perfils_wps/scripts/run_analysis.py
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/v16_4_ancoras_wps_mpl python3 src2/v16.4_ancoras_perfils_wps/scripts/render_results.py
PYTHONDONTWRITEBYTECODE=1 python3 src2/v16.4_ancoras_perfils_wps/scripts/verify_results.py
```

As pastas irmãs v16.1, v16.3, v16_spatial2_18_71 e data_processed continuam necessárias; o catálogo da v16.3 também verifica suas fontes da v16.2. A primeira execução congela as entradas; mudanças posteriores são recusadas. Nenhum mapa ou resultado das análises anteriores é reescrito. O protocolo declara a triagem que antecedeu os controles.

- [Protocolo e limitações](PROTOCOLO.md) · [configuração](config.json).
- [Resultados de todos os métodos](results/concordancia.csv) · [resultados por padrão](results/por_padrao.csv).
- [Todos os pares](results/pares.csv) · [retiradas de aquisição](results/influencia_aquisicoes.csv).
- [Verificação numérica](results/verificacao.json) · [fontes congeladas](results/entradas_congeladas.json).
- [Manifesto de artefatos](results/artefatos.json) · [verificação final](results/verificacao_final.json).

Recomendação: manter o panorama completo como referência. Múltiplas âncoras/janelas por aquisição são uma hipótese para outro experimento, ainda não executado.
"""
    (ROOT / "README.md").write_text(readme, encoding="utf-8")


def main():
    FIG.mkdir(exist_ok=True)
    checks = json.loads((OUT / "verificacao.json").read_text())
    assert checks["status"] == "PASS"
    summary = pd.read_csv(OUT / "concordancia.csv").set_index("metodo")
    pairs = pd.read_csv(OUT / "pares.csv")
    by_pattern = pd.read_csv(OUT / "por_padrao.csv")
    index = pd.read_csv(OUT / "indice_ensaios.csv")
    profiles = pd.read_csv(OUT / "perfis_normalizados.csv")
    archive = dict(np.load(OUT / "perfis.npz", allow_pickle=False))
    deletion = pd.read_csv(OUT / "influencia_aquisicoes.csv")
    effects = pd.read_csv(OUT / "efeitos_pareados.csv")
    figures = make_figures(summary, pairs, by_pattern, index, profiles, archive)
    render_documents(summary, pairs, by_pattern, deletion, effects, checks, figures)
    presentation = dict(figuras=[dict(arquivo=f"figures/{name}.png", titulo=title, leitura=caption,
        sha256=digest(FIG / f"{name}.png"), tamanho_px=list(Image.open(FIG / f"{name}.png").size))
        for name, title, caption in figures], imagens_incorporadas_html=len(figures),
        codigo_sha256=digest(__file__), tabelas_sha256={p.name: digest(p) for p in sorted(OUT.glob("*.csv"))})
    (OUT / "apresentacao.json").write_text(json.dumps(presentation, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    paths = [p for p in ROOT.rglob("*") if p.is_file() and p.name not in ("artefatos.json", "verificacao_final.json") and "__pycache__" not in p.parts]
    manifest = {str(p.relative_to(ROOT)): dict(sha256=digest(p), bytes=p.stat().st_size) for p in sorted(paths)}
    (OUT / "artefatos.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Gerados README, relatório, HTML independente e {len(figures)} figuras.")


if __name__ == "__main__":
    main()
