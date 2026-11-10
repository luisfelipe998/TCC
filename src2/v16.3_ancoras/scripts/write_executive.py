"""Write the requested concise Portuguese report from verified output tables."""
import sys
sys.dont_write_bytecode=True
import json
from collections import Counter
from common import ROOT, INTER


def number(value,digits=3):
    return f"{value:.{digits}f}".replace(".",",")


def main():
    stats=json.loads((ROOT/"04_estatistica/dados/resultados.json").read_text())
    sv=json.loads((INTER/"estatistica/verificacao.json").read_text())
    wv=json.loads((INTER/"lote_inter/verificacao.json").read_text())
    assert sv["status"]==wv["status"]=="PASS"
    assert (stats["ensaios"],stats["intra"],stats["inter"])==(25,87,213)
    aggregates={(v["tipo"],v["metrica"]):v for v in stats["resumos"]}
    groups={(v["padrao_x"],v["padrao_y"]):v for v in stats["grupos"]}
    labels=json.loads((INTER/"lote_inter/rotulos_visuais.json").read_text())
    counts=Counter(v["rotulo"] for v in labels.values())
    intra_labels=json.loads((INTER/"lote_intra/rotulos_visuais.json").read_text())
    assert sum(v["rotulo"]=="favoravel" for v in intra_labels.values())==36
    ri=aggregates["intra","pearson_r"]; re=aggregates["inter","pearson_r"]
    wi=aggregates["intra","wco_media_central"]; we=aggregates["inter","wco_media_central"]
    table=["| Padrão | Ensaios | Contexto nativo | Pares intra | Pearson mediano | WCO central mediana |",
           "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for pattern,count,native in (("plug",3,"2.048 px"),("slug",3,"2.048 px"),
                                 ("churn",10,"4.096 px"),("wavy annular",9,"4.096 px")):
        row=groups[pattern,pattern]
        table.append(f"| {pattern} | {count} | {native} | {row['pares']} | {number(row['r_mediana'])} | {number(row['wco_mediana'])} |")
    intertable=["| Padrões diferentes | Pares | Pearson mediano | WCO central mediana |",
                "| --- | ---: | ---: | ---: |"]
    for a,b in (("plug","slug"),("plug","churn"),("plug","wavy annular"),
                ("slug","churn"),("slug","wavy annular"),("churn","wavy annular")):
        row=groups.get((a,b),groups.get((b,a)))
        intertable.append(f"| {a} × {b} | {row['pares']} | {number(row['r_mediana'])} | {number(row['wco_mediana'])} |")
    worst=stats["maiores_wco_inter"][0]
    best=stats["menores_wco_inter"][0]
    highest=stats["maiores_r_inter"][0]
    fullintra=aggregates["intra","pearson_r_completo"]["mediana"]
    fullinter=aggregates["inter","pearson_r_completo"]["mediana"]
    text=f"""# Resumo executivo — âncoras, wavelets e correlação

**Foram analisados 25 ensaios, com uma âncora por ensaio, 87 pares intra padrão e 213 pares inter padrão. Houve diferenças locais, mas Pearson e WCO apresentaram sobreposição entre padrões: os resultados não demonstraram uma separação clara das quatro classes.**

## 01 — Âncoras e janelas

As âncoras foram escolhidas antes das wavelets, usando a média e a fração escura das linhas **18–71**. Em plug/slug, a seleção combinou eventos de morfologia aprovada visualmente na v16.2, estabilidade das fronteiras nos limiares 55/60/65 e contraste das extremidades. Em churn/wavy annular, venceu o vale mais proeminente em σ = 16, entre candidatos localizáveis também em σ = 12 e 20, com contexto completo.

Os 25 ensaios tiveram âncora válida. Os contextos nativos permaneceram nas imagens de `01_ancoras`. Para comparar qualquer par em condições iguais, utilizou-se o perfil comum de **2.048 px**, preservando o centro de cada âncora e as margens exigidas, sem reamostragem. Essa mesma janela foi usada nas estatísticas principais.

{chr(10).join(table)}

![churn — panorama, contexto nativo de 4.096 px e média das linhas 18–71](../01_ancoras/churn_m03q05_001_MIN02256.png)

O vale escolhido permaneceu localmente estável nos 19 ensaios de churn/wavy annular, com deslocamento máximo de **4 px**, abaixo da tolerância de 16 px. Entretanto, refazer todo o ranking com outro σ mudou o vencedor em **5 dos 19 contextos nativos**. Quatro centros também deixariam de ser os primeiros se a seleção fosse reiniciada diretamente em 2.048 px; foram preservados os vencedores nativos. Assim, “melhor âncora” depende das regras e da janela, sem significar um ótimo universal. A aprovação de forma em plug/slug ainda depende das anotações visuais.

## 02 — Comparações dentro do mesmo padrão

Os **87 pares** utilizaram Morlet, 200 comprimentos de onda de **50 a 1.024 px** e a correlação enviesada com divisor N. Foram mantidas as definições alternativas de Schuck–Bodmann: WPS como módulo da CWT da autocorrelação; WCS como CWT da correlação cruzada; WCO como seu módulo ao quadrado dividido pelas WPS individuais no lag zero. A curva **potência média** resumiu a WPS nos lags com suporte geométrico.

A WCO variou conforme a escala e o ensaio. A mediana das médias centrais dos pares intra foi **{number(wi['mediana'])}**, com valores de {number(wi['minimo'])} a {number(wi['maximo'])}. Faixas locais de alta coerência ocorreram mesmo quando a média foi baixa. A revisão visual anterior marcou 36 pares como favoráveis e 51 como desfavoráveis; esses rótulos não foram usados como taxa de acerto.

![plug — exemplo de coerência localizada dentro do mesmo padrão](../02_wavelets/plug_m05q01_001_m06q01_001_favoravel.png)

## 03 — Comparações entre padrões diferentes

Foram gerados todos os **213 pares**, sem misturar grades de 2.048 e 4.096 px. As WPS individuais previamente verificadas foram reutilizadas; apenas os termos cruzados de cada par precisaram ser calculados.

{chr(10).join(intertable)}

A mediana inter da WCO média central foi **{number(we['mediana'])}**, próxima da intra ({number(wi['mediana'])}); **{stats['inter_wco_acima_mediana_intra']} dos 213 pares inter** atingiram ou superaram a mediana intra. Churn × wavy annular teve a maior mediana inter de WCO, inclusive acima das medianas internas de ambos os padrões. Isso dificultou distinguir essas classes por esse resumo.

O par churn **m06q20_001** × wavy annular **m04q40_001** apresentou a maior média central inter, **{number(worst['wco_media_central'])}**, ilustrada abaixo. No outro extremo, [slug m05q05_001 × wavy annular m07q80_001](../03_wavelets_inter/slug_wavy_annular_m05q05_001_m07q80_001_favoravel.png) apresentou **{number(best['wco_media_central'])}**. São exemplos extremos do conjunto completo, não representantes escolhidos para resumir cada classe.

![churn × wavy annular — maior WCO média central entre padrões](../03_wavelets_inter/churn_wavy_annular_m06q20_001_m04q40_001_desfavoravel.png)

Nos nomes dos arquivos inter, **favoravel significa WCO predominantemente baixa**, conforme solicitado. A inspeção resultou em {counts['favoravel']} favoráveis e {counts['desfavoravel']} desfavoráveis. A avaliação considerou a localização e a extensão das faixas coerentes; não corresponde a um corte automático na média nem a significância estatística.

## 04 — Matriz de Pearson e sensibilidade

Reproduziu-se o descritor de `estatistica5`: **WPS no lag zero, em log, após retirada da tendência log-linear em função de λ**. O Pearson assinado comparou a forma desses resíduos, com um descritor fixo por ensaio. Não foi calculado diretamente sobre a WCO ou sobre a curva potência média.

![Matriz completa de Pearson dos 25 ensaios](../04_estatistica/matriz_pearson_todos.png)

Os números dos eixos estão identificados no [índice dos ensaios](../04_estatistica/INDICE_ENSAIOS.md). Roxo indica correlação negativa; creme, próxima de zero; azul, positiva. A diagonal e os pares espelhados foram excluídos dos resumos.

A mediana de Pearson foi **{number(ri['mediana'])} intra** e **{number(re['mediana'])} inter**. Houve {ri['negativos']} pares intra e {re['negativos']} inter com r negativo. **{stats['inter_r_acima_mediana_intra']} pares inter** atingiram ou superaram a mediana intra. Wavy annular teve a maior mediana interna ({number(groups['wavy annular','wavy annular']['r_mediana'])}); plug, a menor ({number(groups['plug','plug']['r_mediana'])}), com apenas três ensaios disponíveis nessa classe.

O maior r inter ocorreu em **churn m07q20_001 × plug m06q01_001: {number(highest['pearson_r'])}**. Nesse par, a WCO média central foi {number(highest['wco_media_central'])}. Isso mostra a diferença entre os indicadores: forma espectral semelhante não exige coerência cruzada elevada junto ao alinhamento.

![Distribuições de Pearson e WCO dentro e entre padrões](../04_estatistica/correlacoes_dentro_entre.png)

Na figura, cada ponto é um par; os traços indicam mediana e quartis. A **WCO média central** foi calculada nas 200 escalas e nos lags −32 a +32 px, todos com suporte válido. O resumo é descritivo e pode esconder faixas locais importantes.

Como controle, repetiu-se Pearson com os panoramas completos, mantendo a mesma faixa de λ. As medianas passaram a **{number(fullintra)} intra** e **{number(fullinter)} inter**. A mudança absoluta mediana de r foi **{number(stats['delta_r']['mediana'])}**, chegando a **{number(stats['delta_r']['maximo'])}**. O resultado dependeu do contexto espacial; comprimento e conteúdo mudaram juntos. A [figura de sensibilidade](../04_estatistica/sensibilidade_comprimento.png) e a [matriz de controle](../04_estatistica/matriz_pearson_completo.png) foram preservadas. A janela principal não foi substituída pelo resultado mais favorável.

## Confiabilidade e interpretação

**As verificações numéricas passaram.** Os 213 mapas inter foram confrontados com correlação direta e Morlet no domínio temporal; a maior diferença de WCO foi {wv['maior_erro_wco']:.2e}. Os 300 coeficientes de cada matriz foram conferidos com SciPy e NumPy, com diferença máxima de {max(v['maior_erro_scipy'] for v in sv['matrizes']):.2e}. As fontes, âncoras, precomputações individuais e resultados intra anteriores permaneceram intactos. Os PNGs foram verificados para inserção legível em A4 retrato, sem exportação em PDF.

A auditoria anterior já havia reconstruído os 25 panoramas a partir de 12.525 quadros e conferido os 44 contextos numéricos individuais. A coluna de costura coincidente com uma âncora de plug teve efeito pontual pequeno no teste realizado; isso não equivale a avaliar a remoção de todas as costuras.

O estudo permaneceu **descritivo, com K = 1 e um evento por ensaio**. Os 300 pares compartilham apenas 25 aquisições; não são 300 réplicas independentes. As maiores escalas têm poucas repetições no recorte. Não foram inferidas causalidade, significância ou acurácia de classificação. As âncoras e os mapas ficaram consistentes e reutilizáveis, mas Pearson e WCO, isoladamente nesses recortes, não bastaram para demonstrar distinção clara entre os padrões.

[Método, fórmulas e dados](../04_estatistica/METODOLOGIA.md) · [Auditoria das âncoras e intra](../intermediario/auditoria_lote_completo/relatorio.txt) · [Verificação inter](../intermediario/lote_inter/verificacao.json) · [Verificação estatística](../intermediario/estatistica/verificacao.json)
"""
    dest=ROOT/"05_conclusao"
    dest.mkdir(exist_ok=True)
    (dest/"Resumo Executivo.md").write_text(text,encoding="utf-8")
    print(f"Resumo escrito: {len(text.split())} palavras; 5 figuras incorporadas.")


if __name__=="__main__":
    main()
