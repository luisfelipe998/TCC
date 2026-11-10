"""Documentação e galerias da análise executada, a partir das tabelas verificadas."""
from html import escape
from itertools import combinations
import json
import sys
from urllib.parse import quote
sys.dont_write_bytecode = True
from common import *


def write(path,text):
    Path(path).write_text(text.strip()+"\n",encoding="utf-8")


def gallery(filename,title,items):
    cards=[]
    for i,(label,relative,detail) in enumerate(items):
        target=quote("../../"+relative)
        cards.append(f'<details data-search="{escape(label.lower(),quote=True)}" {"open" if i==0 else ""}>'
            f'<summary>{escape(label)}<small>{escape(detail)}</small></summary>'
            f'<a href="{target}"><img src="{target}" alt="{escape(label,quote=True)}" loading="lazy"></a></details>')
    write(INTER/"galerias"/filename,'''<!doctype html>
<html lang="pt-BR"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>'''+escape(title)+'''</title><style>
body{font:16px system-ui,sans-serif;color:#263c48;background:#f4f6f8;margin:0}main{max-width:1040px;margin:30px auto;padding:0 18px}
h1{font-size:26px}nav{display:flex;gap:18px;flex-wrap:wrap}a{color:#086879}input{padding:12px;font:inherit;width:min(94%,620px);border:1px solid #bac8cf;border-radius:5px}
details{margin:16px 0;padding:16px;background:white;border:1px solid #dce3e8;border-radius:7px}summary{cursor:pointer;font-weight:600;line-height:1.6}
small{display:block;font-size:14px;font-weight:400;color:#536774}img{display:block;width:100%;max-width:900px;height:auto;margin:18px auto 0}p{line-height:1.6}
</style></head><body><main><nav><a href="01_panoramas.html">Panoramas</a><a href="02_wavelets.html">Intrapadrão</a><a href="03_wavelets_inter.html">Entre padrões</a>
<a href="../../05_conclusao/Resumo%20Executivo.md">Resumo executivo</a></nav>
<h1>'''+escape(title)+'''</h1><p>Abra um item para visualizar a figura. Clique na imagem para abrir o PNG original.</p>
<label for="filter">Filtrar por ensaio ou padrão</label><p><input id="filter" type="search" placeholder="Ex.: m05q01_001 ou plug"></p>
<p id="count"></p>'''+"\n".join(cards)+'''
<script>const field=document.getElementById('filter'),items=[...document.querySelectorAll('details')];
function filter(){let total=0;const term=field.value.trim().toLowerCase();for(const item of items){item.hidden=!item.dataset.search.includes(term);if(!item.hidden)total++;}
document.getElementById('count').textContent=total+' de '+items.length+' figuras';}field.addEventListener('input',filter);filter();</script>
</main></body></html>''')


def reference_pages(acquisitions,pair_rows):
    index=["# Índice dos ensaios", "", "A mesma ordem foi utilizada nas matrizes e nos perfis estatísticos.", "",
           "| Nº | Ensaio | Padrão | Extensão (px) |", "| ---: | --- | --- | ---: |"]
    for row in acquisitions:
        index.append(f"| {row['numero']} | {row['ensaio']} | {row['padrao']} | {row['comprimento']} |")
    write(STATS/"INDICE_ENSAIOS.md","\n".join(index))
    averages=read_csv(TABLES/"medias_por_ensaio.csv")
    table=["# Correlação média por ensaio", "",
        "Calculou-se a média aritmética assinada de Pearson contra os outros 25 ensaios. A diagonal foi excluída. As médias intra e inter utilizaram somente os respectivos parceiros. O padrão dispersed bubble não teve segundo ensaio para uma média intrapadrão.", "",
        "| Ensaio | Padrão | Média geral | Média intra | Média inter |", "| --- | --- | ---: | ---: | ---: |"]
    for row in averages:
        intra=number(float(row['pearson_media_intra'])) if row['pearson_media_intra'] else "—"
        table.append(f"| {row['ensaio']} | {row['padrao']} | {number(float(row['pearson_media_todos']))} | {intra} | {number(float(row['pearson_media_inter']))} |")
    write(STATS/"MEDIAS_POR_ENSAIO.md","\n".join(table))
    method=r'''# Metodologia da análise espacial integral

## Aquisições e extração

Foram auditados 31 panoramas v11 e as respectivas planilhas originais. Mantiveram-se os 26 ensaios elegíveis da v16: 10 de churn, 9 de wavy annular, 3 de plug, 3 de slug e 1 de dispersed bubble. Foram calculados os 325 pares únicos, dos quais 87 foram intrapadrão e 238 entre padrões. Os cinco registros anteriormente excluídos permaneceram na auditoria, sem ingresso nas estatísticas.

O sinal foi extraído pela média aritmética das linhas 18–71, inclusive, totalizando 54 linhas. Os pixels foram novamente lidos dos PNGs originais e conferidos com os sinais arquivados. A montagem dos panoramas não foi refeita. As planilhas originais representaram a extração anterior 20–59 e serviram à conferência de procedência, sem substituir o sinal 18–71.

Cada panorama completo originou um perfil espectral individual fixo. Para os mapas de um par, empregaram-se os primeiros N = min(Na, Nb) pixels de ambos os sinais. O corte precedeu a inversão 255 − intensidade, a remoção da tendência linear, a centralização e a divisão pelo desvio padrão populacional. Não foram utilizadas âncoras, alinhamento de eventos ou escolha de trecho em função do resultado.

## Correlações e transformadas

Utilizaram-se a amostragem espacial de 1 px, uma realização, Morlet com ω₀ = 6 e 200 comprimentos de onda igualmente espaçados em escala logarítmica entre 50 e 3.000 px. A conversão foi λ = 4πa/(6 + √38). A faixa nativa foi preservada inclusive nos registros mais curtos; não foi substituída pela regra N/2 da análise de âncoras.

As correlações bilaterais usaram o estimador enviesado:

$$R_{xy}[m]=\frac{1}{N}\sum_{n\;\mathrm{válidos}}x[n+m]\,\overline{y[n]}.$$

O divisor N constituiu a adaptação já adotada pela v16. As definições alternativas utilizadas foram:

$$WPS_x(m,a)=|\mathrm{CWT}\{R_{xx}\}(m,a)|,$$

$$WCS_{xy}(m,a)=\mathrm{CWT}\{R_{xy}\}(m,a),$$

$$WCO_{xy}(m,a)=\frac{|WCS_{xy}(m,a)|^2}{WPS_x(0,a)\,WPS_y(0,a)}.$$

A WPS foi o módulo da transformada da autocorrelação; esse módulo não foi novamente elevado ao quadrado. A WCS foi armazenada complexa e apresentada em módulo. O denominador da WCO foi avaliado no lag zero. O piso numérico de 10⁻¹⁰ dos máximos individuais protegeu o denominador; todas as 200 escalas foram válidas em todos os pares. Os valores de WCO não foram truncados artificialmente em 1.

Foram utilizados todos os lags de −(N−1) a +(N−1), sem COI ou máscara geométrica. A potência média ao lado de cada WPS correspondeu à média aritmética do módulo em todos os lags. Os eixos dos sinais representaram a extensão do panorama; os eixos dos mapas, o lag da correlação. Essas coordenadas não foram tratadas como equivalentes.

## Resumo da coerência

Preservou-se o escore principal da v16. Para cada escala, calculou-se o máximo da WCO em todos os lags. Esses máximos foram ponderados por √(WPSa(0) WPSb(0)), retendo as escalas em que ambas as WPS atingiram pelo menos 5% dos respectivos máximos. As sensibilidades aos pisos de 1% e 10%, o máximo integral e a WCO ponderada no lag zero também foram registrados.

O piso de potência limitou somente as escalas empregadas no escore, sem remover regiões dos mapas. O máximo de cada escala pôde ocorrer em um lag diferente; o escore não representou um único alinhamento entre os sinais. A média integral da WCO foi fornecida como medida adicional e não substituiu o escore principal.

## Pearson dos perfis espectrais

Empregou-se a definição estatística da v16.3 sobre a entrada e a faixa espectral nativas da v16. Para cada ensaio, extraiu-se a WPS no lag zero do panorama completo e ajustou-se:

$$z_{ik}=\log(P_{ik}+10^{-12})-\widehat c_i-\widehat d_i\log(\lambda_k),\qquad P_{ik}=WPS_i(0,a_k).$$

O ajuste de c e d foi feito por mínimos quadrados nas mesmas 200 escalas. Calculou-se Pearson entre os resíduos centralizados, preservando o sinal da correlação. Cada descritor ficou fixo para todos os parceiros do ensaio. Não foi usado Pearson diretamente das intensidades, da WCO ou da curva de potência média.

A matriz utilizou 26 ensaios e λ de 50–3.000 px. Esse domínio diferiu da comparação anterior entre versões, que utilizou 25 ensaios comuns e λ de 50–1.024 px; os valores não foram comparados diretamente como efeito isolado de âncoras. A diagonal foi unitária e foi excluída das médias. Cada média geral por ensaio reuniu seus 25 parceiros.

## Agregações e alcance

As médias por par atribuíram o mesmo peso a cada par. As médias balanceadas atribuíram peso igual aos quatro padrões com pares intrapadrão e às dez combinações distintas de padrões. Dispersed bubble, com um ensaio, entrou nos 25 controles entre padrões que o envolveram, mas não teve repetibilidade intrapadrão estimada.

Os pares compartilharam aquisições e as escalas foram correlacionadas. Os resumos foram descritivos; não foram produzidos valores-p ou intervalos de confiança tratando pares ou escalas como réplicas independentes. Não houve treinamento de classificador nem estimativa de acurácia nesta análise. As unidades em pixels não foram convertidas em milímetros ou frequência temporal.

## Apresentação e verificação

As métricas e curvas médias foram calculadas com todos os lags, antes de qualquer redução para apresentação. Os arquivos de mapas em `intermediario/pares` conservaram aproximadamente 1.800 posições de lag, incluindo zero e os extremos, os coeficientes complexos, os sinais preparados e os resumos integrais. Esses mapas reduzidos serviram à reprodução das figuras; os máximos integrais não devem ser recalculados apenas deles.

As figuras foram dimensionadas para inserção a 16 cm de largura, com fonte mínima de 8 pt e altura máxima de 24 cm. Foram usadas as mesmas escalas de cor para as duas WPS de cada par e o intervalo 0–1 para exibição da WCO. Os nomes de arquivo identificaram os pares, sem classificação automática como favoráveis ou desfavoráveis.

Os 26 perfis, as métricas dos 325 pares e os mapas intrapadrão foram confrontados com os dados arquivados. Também foram verificadas a correlação por somas diretas, a WPS no lag zero por integração temporal independente da Morlet, a definição de Pearson contra a v16.3, todos os 325 coeficientes por SciPy, as matrizes por NumPy e a integridade dos arquivos arquivados.

Referências de implementação: [núcleo matemático preservado](../scripts/wavelet_core.py), [auditoria original das equações](../_OLD/EQUATION_AUDIT.md), [auditoria dos cálculos atuais](../intermediario/auditoria/calculo.json) e [auditoria estatística](../intermediario/auditoria/estatistica.json).
'''
    write(STATS/"METODOLOGIA.md",method)
    gallery("01_panoramas.html","Panoramas completos · 26 ensaios",[
        (f"{r['padrao']} · {r['ensaio']}",f"01_panoramas/{slug(r['padrao'])}_{r['ensaio']}.png",
         f"Extensão: {r['comprimento']} px · média das linhas 18–71") for r in acquisitions])
    for kind,filename,title in (("intra","02_wavelets.html","Comparações intrapadrão · 87 pares"),
                                ("inter","03_wavelets_inter.html","Comparações entre padrões · 238 pares")):
        gallery(filename,title,[(f"{r['padrao']} · {r['run_x']} × {r['run_y']}",r['arquivo'],
            f"Prefixo comum: {r['comprimento']} px · escore ponderado: {number(r['escore_piso05'])}")
            for r in pair_rows if r['tipo']==kind])


def conclusion(acquisitions,pair_rows):
    result=json.loads((TABLES/"resultados.json").read_text())
    summary={(v['tipo'],v['metrica']):v for v in result['resumos']}
    same_r,diff_r=summary['intra','pearson_r'],summary['inter','pearson_r']
    same_s,diff_s=summary['intra','escore_piso05'],summary['inter','escore_piso05']
    group_rows=[v for v in result['grupos'] if v['tipo']=='intra']
    counts=[]
    for pattern in PATTERNS:
        count=sum(r['padrao']==pattern for r in acquisitions)
        counts.append(f"| {pattern} | {count} | {count*(count-1)//2} |")
    per_pattern=[f"| {v['padrao_x']} | {v['pares']} | {number(v['pearson_media'])} | {number(v['escore_media'])} |" for v in group_rows]
    per_pattern.append("| dispersed bubble | 0 | — | — |")
    summaries=[]
    for metric,label in (("pearson_r","Pearson"),("escore_piso05","Escore de coerência")):
        a,b=summary['intra',metric],summary['inter',metric]
        for field,aggregation in (("media","Peso igual por par"),("media_balanceada","Peso igual por grupo")):
            summaries.append(f"| {label} | {aggregation} | {number(a[field])} | {number(b[field])} | {number(a[field]-b[field])} |")
    lookup={(r['run_x'],r['run_y']):r for r in pair_rows}
    previous=read_csv(OLD/'results/preprocessing_pooled_comparison.csv')
    old_base=next(r for r in previous if r['variant']=='baseline')
    assert abs(float(old_base['mean_score'])-same_s['media'])<1e-12
    variants={'v5':'v5 aplicado à extração 18–71','baseline':'v16: média das linhas 18–71',
              'central_band20':'Controle fixo: média das linhas 30–49',
              'core_row_median':'Controle: mediana das linhas 18–71'}
    prior_table=[f"| {variants[v['variant']]} | {number(float(v['mean_score']),6)} | {number(float(v['balanced_mean']),6)} | {number(float(v['delta_vs_v5']),6)} |" for v in previous]
    bands=read_csv(OLD/'results/band_comparison_summary.csv')
    bands={r['band']:r for r in bands if r['variant']=='baseline' and r['metric']=='weighted_scale_peak_power05'}
    nulls=[]
    for row in read_csv(OLD/'results/null_diagnostic.csv'):
        current=lookup[row['run_x'],row['run_y']]
        assert abs(current['wco_maxima']-float(row['observed_max']))<1e-12
        assert int(row['surrogates'])==199
        assert abs(float(row['conditional_monte_carlo_p'])-(1+int(row['exceedances']))/200)<1e-12
        label=row['pattern'] if current['tipo']=='intra' else 'Controle entre padrões'
        nulls.append(f"| {label} | {number(current['wco_maxima'])} | {number(float(row['conditional_monte_carlo_p']),3)} |")
    example_sections=[]
    for pattern in PATTERNS:
        first=next(r for r in acquisitions if r['padrao']==pattern)
        run=first['ensaio']
        example_sections += [f"### {pattern}", "", f"![Panorama completo e sinais de {run}](../01_panoramas/{slug(pattern)}_{run}.png)", ""]
        same=next((r for r in pair_rows if r['tipo']=='intra' and r['padrao_x']==pattern),None)
        if same:
            example_sections += [f"![Análise wavelet do par {same['run_x']} e {same['run_y']}](../{same['arquivo']})", "",
                f"O primeiro par em ordem alfabética, {same['run_x']} × {same['run_y']}, apresentou escore ponderado de **{number(same['escore_piso05'])}** e WCO máxima de **{number(same['wco_maxima'])}**.", ""]
        else:
            example_sections += ["A única aquisição de dispersed bubble participou de 25 comparações entre padrões. Não houve segundo ensaio independente para comparação intrapadrão.", ""]
    inter=sorted([r for r in pair_rows if r['tipo']=='inter'],key=lambda r:r['escore_piso05'])
    extremes=[]
    for label,row in (("Menor escore entre padrões",inter[0]),("Maior escore entre padrões",inter[-1])):
        extremes += [f"### {label}","",f"![{row['padrao']}](../{row['arquivo']})","",
                     f"{row['run_x']} × {row['run_y']} ({row['padrao']}): escore **{number(row['escore_piso05'])}**.",""]
    body=f'''# Resumo executivo da análise v16

A análise completa foi refeita com **26 ensaios e 325 pares**, mantendo a metodologia espacial da v16 e a apresentação aprovada a partir da v16.3. Foram gerados 26 panoramas com seus sinais, 87 figuras intrapadrão e 238 figuras entre padrões. Os perfis, os resumos dos pares e os mapas intrapadrão reproduziram exatamente os resultados arquivados da implementação original.

O Pearson dos perfis espectrais apresentou média maior dentro do mesmo padrão, **{number(same_r['media'])}**, do que entre padrões, **{number(diff_r['media'])}**. O escore principal de coerência apresentou o comportamento oposto: **{number(same_s['media'])}** dentro e **{number(diff_s['media'])}** entre padrões. A diferença decorreu de medidas com definições e domínios distintos; não foi atribuída a uma mudança das equações.

## Organização e acesso

- [Galeria dos 26 panoramas](../intermediario/galerias/01_panoramas.html).
- [Galeria dos 87 pares intrapadrão](../intermediario/galerias/02_wavelets.html).
- [Galeria dos 238 pares entre padrões](../intermediario/galerias/03_wavelets_inter.html).
- [Metodologia e definições dos resumos](../04_estatistica/METODOLOGIA.md).
- [Índice dos ensaios](../04_estatistica/INDICE_ENSAIOS.md) e [médias de Pearson por ensaio](../04_estatistica/MEDIAS_POR_ENSAIO.md).
- [Tabelas dos pares](../04_estatistica/dados/pares_pearson_wco.csv) e [reprodução](../scripts/README.txt).

O conteúdo anterior foi preservado em `_OLD`. O piloto aprovado foi conservado em `_OLD/piloto_aprovado`. A área ativa reuniu os produtos da análise completa, seus dados e scripts de reprodução.

## Conjunto analisado

| Padrão | Ensaios | Pares intrapadrão |
| --- | ---: | ---: |
{chr(10).join(counts)}

Foram conferidos os 31 panoramas de origem, com manutenção dos cinco registros anteriormente excluídos. A seleção experimental permaneceu a mesma da v16. Os 87 pares intrapadrão e os 238 pares entre padrões compartilharam aquisições, constituindo comparações descritivas dependentes.

## Metodologia aplicada

Extraiu-se o sinal pela média das linhas 18–71 dos panoramas originais. As figuras individuais mostraram os panoramas completos. Nos pares, utilizou-se o prefixo comum de comprimento mínimo entre as duas aquisições; a inversão de intensidade, a remoção da tendência linear e a normalização foram executadas após esse corte.

Foram mantidos Morlet com ω₀ = 6, correlação com divisor N e 200 comprimentos de onda de 50 a 3.000 px. A WPS correspondeu ao módulo da CWT da autocorrelação; a WCS, à CWT da correlação cruzada; a WCO normalizou o módulo quadrático da WCS pelas WPS no lag zero. Todos os lags foram incluídos, sem COI nem máscara geométrica. A potência média ao lado da WPS foi calculada sobre todos esses lags.

As figuras adotaram os títulos por extenso seguidos de WPS, WCS e WCO, “Sinais normalizados”, “Extensão do panorama (px)” e “Altura do panorama (px)”. Não foram inseridas âncoras nem executado alinhamento de eventos. A grade espectral nativa da v16 foi preservada.

## Pearson e coerência

![Correlação de Pearson entre os 26 ensaios](../04_estatistica/matriz_pearson_todos.png)

*A matriz utilizou os perfis de WPS no lag zero dos panoramas completos, após transformação logarítmica e remoção da tendência linear em log(λ). As cores foram fixadas entre −1 e 1. A numeração correspondeu ao índice dos ensaios.*

![Escore de coerência entre os 26 ensaios](../04_estatistica/matriz_escore_wco.png)

*A matriz reuniu os escores ponderados dos máximos de WCO por escala. A diagonal cinza representou comparações do ensaio consigo mesmo, não incluídas no conjunto de pares de aquisições distintas. O piso de 5% foi aplicado ao resumo de potência e não aos mapas.*

| Medida | Agregação | Dentro do padrão | Entre padrões | Dentro − entre |
| --- | --- | ---: | ---: | ---: |
{chr(10).join(summaries)}

A média por par atribuiu peso maior aos padrões com mais combinações. A média balanceada atribuiu peso igual aos quatro padrões com pares intrapadrão e às dez combinações distintas de padrões. Por isso, as duas agregações responderam a ponderações diferentes; ambas foram apresentadas.

![Distribuição das correlações dentro e entre padrões](../04_estatistica/correlacoes_dentro_entre.png)

O contraste médio de Pearson foi positivo, mas houve ampla sobreposição: os valores intrapadrão variaram de {number(same_r['minimo'])} a {number(same_r['maximo'])}, enquanto os valores entre padrões variaram de {number(diff_r['minimo'])} a {number(diff_r['maximo'])}. Assim, correlação elevada também ocorreu entre padrões distintos. A diferença das médias não demonstrou, isoladamente, uma regra de classificação eficaz.

O escore de coerência foi maior entre padrões nas duas agregações. Portanto, este escore não apresentou o comportamento de maior semelhança intrapadrão pretendido para discriminação global. O máximo de cada escala pôde ocorrer em um lag diferente, e os comprimentos comuns dependeram do parceiro. Esse resumo não representou a coincidência de um único evento físico entre aquisições.

## Resultados por padrão

| Padrão | Pares | Pearson médio | Escore médio de coerência |
| --- | ---: | ---: | ---: |
{chr(10).join(per_pattern)}

Plug apresentou a maior correlação média dos perfis dentro do padrão, enquanto slug apresentou o maior escore médio de coerência. Wavy annular apresentou correlação espectral elevada e o menor escore intrapadrão de coerência. Essa diferença confirmou que semelhança da forma espectral e coerência entre sinais não foram medidas equivalentes.

![Perfis espectrais e resíduos utilizados no Pearson](../04_estatistica/perfis_comparados.png)

*As linhas finas representaram cada ensaio; a linha espessa, a mediana do padrão. A normalização pelo máximo foi usada apenas para exibição dos perfis de WPS. Os resíduos unitários foram os descritores utilizados em Pearson. O único ensaio de dispersed bubble não permitiu estimar variabilidade intrapadrão.*

## Exemplos dos padrões

Os exemplos individuais e intrapadrão abaixo foram definidos pela ordem alfabética dos identificadores, sem seleção pelos maiores valores obtidos.

{chr(10).join(example_sections)}

## Exemplos entre padrões

Os dois exemplos seguintes foram selecionados pelos extremos do escore somente para ilustrar sua amplitude. Os 238 pares foram mantidos nas estatísticas e na galeria; os exemplos não constituíram critérios de seleção de dados.

{chr(10).join(extremes)}

## Controles anteriores preservados

Os controles de pré-processamento, faixa de extração e fases aleatórias já executados na v16 foram preservados em `_OLD`. Não foram reexecutados nesta reorganização. A referência principal 18–71 foi recalculada em todos os pares e conferida contra as tabelas desses controles.

### Pré-processamento na entrada 18–71

| Opção | Média dos 87 pares | Média balanceada | Diferença da média por par frente ao v5 |
| --- | ---: | ---: | ---: |
{chr(10).join(prior_table)}

A referência corrigida da v16 apresentou aumento de {number(float(old_base['delta_vs_v5']),6)} na média por par frente ao v5 aplicado à mesma entrada. Esse ganho numérico não comprovou separação entre padrões: o contraste do escore permaneceu negativo. As alternativas conservaram as definições de WPS, WCS e WCO; alteraram somente o pré-processamento. Os pesos desses escores foram definidos pela referência 18–71.

Fonte: [comparação preservada dos pré-processamentos](../_OLD/results/preprocessing_pooled_comparison.csv) e [estudo original](../_OLD/PREPROCESSING.md).

### Faixa de extração

No controle com pesos fixos da referência 20–59, a ampliação para 18–71 alterou a média dos 87 pares de **{number(float(bands['20_59']['pooled_same']),6)}** para **{number(float(bands['18_71']['pooled_same']),6)}**. Essa comparação utilizou os mesmos pesos 20–59 nas duas faixas. Os escores das tabelas principais utilizaram pesos nativos 18–71 e foram apresentados separadamente.

A média balanceada passou de {number(float(bands['20_59']['balanced_same']),6)} para {number(float(bands['18_71']['balanced_same']),6)}, indicando que as conclusões também dependeram da ponderação dos padrões. O resultado desse controle foi misto, sem melhora da média conjunta dos 87 pares.

Fonte: [comparação preservada entre faixas](../_OLD/results/band_comparison_summary.csv).

### Diagnóstico de fases aleatórias

| Exemplo | Máximo observado de WCO | Valor-p condicional |
| --- | ---: | ---: |
{chr(10).join(nulls)}

O controle anterior empregou 199 pares simulados por exemplo e o máximo de WCO no plano integral. Os máximos observados foram novamente reproduzidos pela execução atual. Os valores-p permaneceram condicionais ao modelo de fases aleatórias e não constituíram testes de acurácia, identidade de eventos ou superioridade entre métodos. Não foram acrescentados contornos de significância aos mapas.

Fonte: [diagnóstico preservado](../_OLD/results/null_diagnostic.csv).

## Alcance e conclusão

A análise reproduziu a metodologia integral da v16 e apresentou todos os pares em um formato uniforme. O Pearson apresentou maior semelhança média dentro dos padrões, acompanhada de forte sobreposição com as comparações entre padrões. O escore principal de coerência apresentou médias maiores entre padrões, o que limitou sua utilização direta como medida global de discriminação neste conjunto.

Os resultados não constituíram uma avaliação de classificador. Os pares compartilharam aquisições, os grupos foram desbalanceados e dispersed bubble dispôs de uma única aquisição. As extensões dos panoramas e dos prefixos comuns variaram; o uso de todos os lags conservou regiões de borda e lags com menor quantidade de produtos na correlação. As unidades espaciais do panorama não foram convertidas em tamanho físico de bolhas ou frequência temporal calibrada.

O Pearson atual empregou **26 ensaios e λ de 50–3.000 px**. A comparação anterior entre v16 e âncoras empregou **25 ensaios comuns e λ de 50–1.024 px**. Esses domínios diferentes impediram atribuir a diferença numérica entre tabelas apenas ao uso de âncoras. A comparação controlada anterior permaneceu separada da estatística nativa desta análise.

Verificação: [cálculos e equivalência com a v16](../intermediario/auditoria/calculo.json), [Pearson e agregações](../intermediario/auditoria/estatistica.json). O núcleo matemático permaneceu idêntico ao original. Os mapas foram reduzidos somente para apresentação; todos os resumos foram calculados na resolução integral.
'''
    write(ROOT/'05_conclusao/Resumo Executivo.md',body)


def reproduction():
    write(ROOT/"scripts/README.txt",'''ANÁLISE COMPLETA V16 · LINHAS 18–71

Executar a partir da raiz do projeto, nesta ordem:

python3 src2/v16_spatial2_18_71/scripts/analise.py --workers 2
python3 src2/v16_spatial2_18_71/scripts/estatistica.py
python3 src2/v16_spatial2_18_71/scripts/relatorio.py
python3 src2/v16_spatial2_18_71/scripts/verificar.py

analise.py conferiu 31 entradas, recalculou os 26 perfis individuais e os 325 pares,
e gerou 26 panoramas, 87 figuras intrapadrão e 238 entre padrões. Pares já concluídos
podem ser reutilizados se as assinaturas de código e dados coincidirem. --force
força o recálculo dos pares. O núcleo wavelet_core.py é idêntico ao original.

estatistica.py calculou Pearson, resumos e médias por ensaio, gerando quatro figuras.
relatorio.py escreveu a metodologia, o índice, as galerias e o resumo executivo.
verificar.py conferiu arquivos, fontes, cobertura, equações nos mapas, formatação,
matrizes, links e o conteúdo de _OLD. Os resultados da verificação ficam em
intermediario/auditoria/verificacao_final.json.

Os cálculos usaram todos os lags. Os mapas guardados foram reduzidos somente no
eixo dos lags para apresentação. As curvas médias e máximos salvos vieram das
matrizes integrais e não devem ser recalculados apenas das matrizes reduzidas.

01_panoramas: 26 figuras dos panoramas completos e dos sinais.
02_wavelets: 87 pares intrapadrão no formato aprovado.
03_wavelets_inter: 238 pares entre padrões no mesmo formato.
04_estatistica: quatro figuras, metodologia, índice, médias e tabelas numéricas.
05_conclusao: Resumo Executivo.md, resultados e interpretação.
intermediario: entradas, perfis, mapas, galerias e auditorias necessárias à reprodução.
_OLD: 583 arquivos originais e o piloto aprovado, preservados.

Os controles históricos de pré-processamento, extração e fases aleatórias foram
mantidos em _OLD e citados como resultados anteriores. Os cálculos principais da
versão atual foram reexecutados e confrontados com esses resultados de referência.

Dependências: numpy, scipy, matplotlib e Pillow; PyWavelets é importado pelo núcleo
original para manter os controles da implementação histórica disponíveis.
''')


if __name__ == "__main__":
    reference_pages(records(),pairs())
    conclusion(records(),pairs())
    reproduction()
