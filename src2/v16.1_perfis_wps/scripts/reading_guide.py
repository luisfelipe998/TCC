"""Apresentação didática derivada dos resultados existentes, sem refazer o método."""
from pathlib import Path
import hashlib
import json
import os

os.environ.setdefault('MPLCONFIGDIR', '/tmp/tcc-wps-leitura-mpl')

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'results'
OUT = ROOT / 'figures' / 'leitura'
PATTERNS = ['Bolhas dispersas', 'Bolhas alongadas (plug)', 'Golfadas (slug)',
            'Agitado (churn)', 'Anular ondulado']
SOURCES = ['all_pairs.csv', 'profiles_normalized.csv', 'ranking_summary.csv',
           'ranking_by_pattern.csv', 'selected_cases.csv', 'acquisition_deletion.csv']
TEAL, ORANGE, PURPLE = '#147d80', '#bf701b', '#79629d'


def number(value, digits=3):
    return f'{value:.{digits}f}'.replace('.', ',')


def integer(value):
    return f'{int(value):,}'.replace(',', '.')


def percent(value):
    return number(100 * value, 1) + '%'


def table(headers, rows):
    return ('| ' + ' | '.join(headers) + ' |\n| ' + ' | '.join(['---'] * len(headers))
            + ' |\n' + '\n'.join('| ' + ' | '.join(map(str, row)) + ' |' for row in rows))


def evidence():
    pairs = pd.read_csv(RESULTS / 'all_pairs.csv')
    profiles = pd.read_csv(RESULTS / 'profiles_normalized.csv')
    primary = profiles[(profiles.variante == 'principal') & (profiles.dominio == 'principal')]
    ranks = pd.read_csv(RESULTS / 'ranking_summary.csv').set_index('metodo')
    main = pairs[pairs.metodo == 'wps_principal']
    same = main[main.mesmo_padrao].valor.to_numpy()
    cross = main[~main.mesmo_padrao].valor.to_numpy()
    outcomes = dict(favorable=int((same[:, None] < cross).sum()),
                    ties=int((same[:, None] == cross).sum()),
                    unfavorable=int((same[:, None] > cross).sum()),
                    total=len(same) * len(cross))
    assert sum(outcomes[k] for k in ['favorable', 'ties', 'unfavorable']) == outcomes['total']
    assert np.isclose((outcomes['favorable'] + .5 * outcomes['ties']) / outcomes['total'],
                      ranks.loc['wps_principal', 'A_pares'], atol=1e-12, rtol=0)
    cases = pd.read_csv(RESULTS / 'selected_cases.csv')
    # Mesmos três casos já escolhidos na análise original; nenhuma nova seleção otimizada.
    cases = cases.set_index('tipo').loc[
        ['Perfis próximos com escore WCO baixo', 'Controle cruzado mais semelhante',
         'Par intrapadrão menos semelhante']].reset_index()
    for case in cases.itertuples():
        pair = main[(main.aquisicao_x == case.aquisicao_x) &
                    (main.aquisicao_y == case.aquisicao_y)]
        assert len(pair) == 1 and np.isclose(pair.iloc[0].valor, case.H, atol=1e-12, rtol=0)
    counts = primary.groupby('padrao').aquisicao.nunique().reindex(PATTERNS)
    assert counts.sum() == 26 and len(same) == 87 and len(cross) == 238
    return primary, main, ranks, cases, counts, outcomes


def build_summary():
    primary, pairs, ranks, cases, counts, outcomes = evidence()
    a = ranks.loc['wps_principal']
    wco = ranks.loc['wco_v16_18_71']
    prefix = ranks.loc['wps_prefixo_par']
    pats = pd.read_csv(RESULTS / 'ranking_by_pattern.csv').set_index(['metodo', 'padrao'])
    influence = pd.read_csv(RESULTS / 'acquisition_deletion.csv')
    deletion = influence[influence.metodo == 'wps_principal']
    bubble = deletion[deletion.aquisicao_retirada == 'm07q01_001'].iloc[0]
    count_table = table(['Padrão já atribuído à aquisição', 'Aquisições', 'Pares do mesmo padrão'],
                        [(p, int(counts[p]), int(counts[p] * (counts[p] - 1) // 2)) for p in PATTERNS])
    method_table = table(['Forma de comparação', 'Cada par com o mesmo peso',
                          'Padrões e combinações com pesos iguais'], [
        ('WCO anterior: trecho inicial comum a cada par', percent(wco.A_pares), percent(wco.A_balanceada)),
        ('Perfis: o mesmo trecho e a mesma faixa da WCO', percent(prefix.A_pares), percent(prefix.A_balanceada)),
        ('Perfis principais: registros inteiros', percent(a.A_pares), percent(a.A_balanceada))])
    case_table = table(['Exemplo ilustrado', 'Aquisições', 'Distância H'], [
        (label, f'`{row.aquisicao_x}` e `{row.aquisicao_y}`', number(row.H))
        for label, row in zip(['Mesmo padrão, curvas próximas', 'Padrões diferentes, curvas próximas',
                               'Mesmo padrão, curvas mais diferentes'], cases.itertuples())])
    slug = pats.loc[('wps_principal', 'Golfadas (slug)')]
    slug_wco = pats.loc[('wco_v16_18_71', 'Golfadas (slug)')]
    inside_dominant = sum(int(counts[p] * (counts[p] - 1) // 2)
                          for p in ['Agitado (churn)', 'Anular ondulado'])
    return f'''# Leitura guiada: das imagens à comparação entre padrões

**Os registros do mesmo padrão tenderam a apresentar curvas mais parecidas, mas houve exceções importantes.** A análise encontrou uma forma útil de descrever semelhanças neste conjunto. Ainda não demonstrou identificação automática dos padrões em novos experimentos.

O percurso de leitura é: **imagem → sinal → curva de escalas → distância entre duas curvas → comparação de todos os pares**. As siglas e os detalhes de reprodução aparecem depois dessa explicação.

## 1. Objetivo e dados da comparação

O estudo investigou se imagens de escoamentos com o mesmo padrão também se pareciam quando transformadas em curvas numéricas. Os padrões já estavam identificados na base. Esses rótulos foram usados para avaliar a semelhança calculada; o método não os descobriu automaticamente.

Uma **aquisição** é um registro experimental. Neste estudo, cada aquisição forneceu um panorama, uma imagem comprida montada com partes de imagens capturadas em momentos diferentes. Foram verificadas 31 entradas e mantidas as 26 aquisições previstas no protocolo, com as cinco exclusões anteriores.

{count_table}

Um **par** contém duas aquisições. Com três aquisições de bolhas alongadas, por exemplo, foram possíveis três pares: primeira–segunda, primeira–terceira e segunda–terceira. Ao combinar as 26 aquisições, obtiveram-se **325 pares: 87 do mesmo padrão e 238 de padrões diferentes**. Uma aquisição participou de vários pares; por isso, 325 pares não equivalem a 325 experimentos independentes.

## 2. Da imagem à curva de escalas

Cada coluna do panorama foi resumida pela média dos tons de cinza das linhas 18–71, totalizando 54 linhas. A sequência dessas médias formou o sinal. Depois, inverteram-se as intensidades, retirou-se uma tendência linear e padronizou-se o sinal para reduzir diferenças globais de nível e amplitude.

Para resumir as variações, calculou-se a autocorrelação, que compara o sinal com cópias deslocadas dele próprio. Aplicou-se então a transformada wavelet, que examina oscilações de diferentes escalas. Dela foi extraída uma curva, chamada aqui **perfil de WPS alternativa**, no deslocamento zero. Esse deslocamento pertence à autocorrelação e não significa a primeira coluna da imagem: o cálculo utilizou o registro analisado.

![Figura 1. Um panorama de bolhas alongadas, o sinal extraído e a curva usada na comparação.](figures/conclusao/01_imagem_sinal_perfil.png)

**Leitura da figura 1, de cima para baixo:**

- **Imagem:** as linhas vermelhas delimitam a faixa utilizada. A posição horizontal indica onde está cada coluna do panorama.
- **Sinal:** a curva mostra as intensidades depois do preparo. O eixo horizontal ainda indica posição no panorama.
- **Perfil:** o eixo horizontal mudou de significado. Agora indica a escala espacial da variação, em pixels. À esquerda ficam variações mais curtas; à direita, mais longas.
- **Picos do perfil:** indicam escalas com maior resposta relativa. Um pico próximo de 500 px não localiza uma bolha na coluna 500 e não mede diretamente seu comprimento.

Cada perfil teve sua área ajustada para um, para comparar a distribuição da resposta entre escalas. A altura é uma **densidade relativa**, não uma porcentagem de energia em cada ponto. O eixo é logarítmico: distâncias iguais no desenho correspondem a razões iguais, como dobrar de 100 para 200 px ou de 500 para 1.000 px. A área foi calculada nesse eixo.

A faixa principal foi de **50 a aproximadamente 2.296 px**, definida pelo suporte geométrico do registro mais curto. Os dados não foram convertidos em milímetros ou hertz. O panorama foi exibido com proporção ajustada para leitura.

## 3. Distância entre duas curvas

Para transformar a diferença de forma em um número, aplicou-se a **distância de Hellinger, H**. Ela considera toda a curva normalizada na faixa escolhida, e não apenas o maior pico.

**H = 0 significa curvas idênticas; valores maiores indicam maior diferença.** O limite matemático é 1. Um H de 0,084 não significa 8,4% de erro nem 91,6% de chance de serem do mesmo padrão. Não foi definido um valor de corte para classificar aquisições.

![Figura 2. Três pares reais mostram a relação entre as curvas, a distância e os rótulos já conhecidos.](figures/leitura/01_pares_explicados.png)

**Leitura da figura 2:** cada painel contém duas curvas, uma por aquisição. No primeiro, duas aquisições anulares apresentaram curvas próximas. No segundo, bolhas alongadas e anular ondulado também apresentaram curvas próximas, apesar dos padrões diferentes. No terceiro, duas aquisições anulares apresentaram curvas mais diferentes. As cores e os traços distinguem as aquisições dentro de cada painel.

{case_table}

O contraste é concreto: o par de padrões diferentes teve H = **{number(cases.iloc[1].H)}**, menor que **{number(cases.iloc[2].H)}** do último par do mesmo padrão. Portanto, semelhança da curva não identifica o padrão em todos os casos. Esses três exemplos já haviam sido selecionados após os resultados para explicar seu comportamento; não representam a frequência de cada situação. Essa frequência exige olhar todos os pares.

## 4. Distribuição de todos os pares

![Figura 3. Distribuição das distâncias dos 87 pares do mesmo padrão e dos 238 pares de padrões diferentes.](figures/leitura/02_distribuicao_distancias.png)

**Leitura da figura 3:** o eixo horizontal mostra H; mais à esquerda significa curvas mais parecidas. Cada faixa agrupa distâncias próximas. O eixo vertical mostra a porcentagem dos pares de cada grupo naquela faixa. Cada grupo soma 100%, permitindo comparar grupos de tamanhos diferentes.

Os pares do mesmo padrão tiveram distância média **{number(a.media_intra)}**, contra **{number(a.media_entre)}** entre padrões diferentes. A tendência desejada apareceu, mas as duas distribuições se sobrepuseram: nas mesmas faixas de distância ocorreram os dois tipos de par. As médias assinaladas não são limites de classificação. O gráfico apresenta o intervalo observado; o intervalo matemático de H continua sendo de 0 a 1.

## 5. Significado dos {percent(a.A_pares)}

A análise usou uma estatística chamada **concordância de ordenação, A**. Para entendê-la, considere uma comparação entre dois pares: um do mesmo padrão e outro de padrões diferentes. A ordem esperada ocorreu quando o primeiro teve a menor distância.

![Figura 4. Contagem das comparações que produziram a concordância de ordenação dos perfis principais.](figures/leitura/03_calculo_concordancia.png)

**Leitura da figura 4:** cada um dos 87 pares do mesmo padrão foi comparado a cada um dos 238 pares de padrões diferentes. Isso produziu **87 × 238 = {integer(outcomes['total'])} comparações entre pares**. Em **{integer(outcomes['favorable'])}**, o par do mesmo padrão teve a menor distância; em **{integer(outcomes['unfavorable'])}**, ocorreu a ordem contrária. Houve **{integer(outcomes['ties'])} empates**; a regra atribui meio crédito a cada empate.

**A = ({integer(outcomes['favorable'])} + 0,5 × {integer(outcomes['ties'])}) / {integer(outcomes['total'])} = {percent(a.A_pares)}.**

Portanto, {percent(a.A_pares)} descreveu a frequência da ordenação esperada nessas comparações. **Não é a porcentagem de aquisições classificadas corretamente.** Nenhum classificador foi treinado nessa etapa. As {integer(outcomes['total'])} comparações reutilizaram os mesmos pares e aquisições, sem criar novas observações independentes.

Na escala de A, **50% significa ausência de preferência de ordenação**; acima disso, a ordem esperada predominou; abaixo disso, predominou a ordem contrária. Essa referência não é um teste de significância estatística.

## 6. Comparação com a análise anterior

A análise anterior, denominada **v16**, utilizou um escore de **coerência wavelet alternativa, WCO**. Ele resumiu a correspondência entre dois sinais em escalas e deslocamentos. A **v16.1** passou a comparar diretamente o formato dos perfis. A ideia é semelhante a comparar duas músicas pelo desenho no equalizador: distribuições parecidas de graves e agudos podem existir mesmo quando a sequência de notas é diferente. Nos dados estudados, as variações foram espaciais.

Os dois números têm sentidos diferentes: **na distância H, menor indica mais semelhança; no escore WCO, maior indica mais semelhança**. Por isso, seus valores brutos não foram comparados diretamente. A avaliou a ordenação produzida por cada método.

![Figura 5. Concordância de ordenação dos métodos, com cada par recebendo o mesmo peso.](figures/leitura/04_comparacao_metodos.png)

**Leitura da figura 5:** as barras representam A, e não H ou a WCO bruta. A barra intermediária é o controle de comparação: usou os mesmos trechos iniciais de cada par e a mesma faixa de 50–3.000 px da WCO. Seu resultado de **{percent(prefix.A_pares)}** mostrou que a vantagem de ordenação não dependeu apenas de utilizar registros inteiros ou de reduzir a faixa de escalas.

{method_table}

**As duas colunas respondem à mesma pergunta com pesos diferentes.** Na primeira, cada par recebeu peso igual. Isso deu muita influência a agitado e anular ondulado, que juntos forneceram {inside_dominant} dos 87 pares internos. Na segunda, os quatro padrões com pares internos receberam o mesmo peso entre si; as dez combinações de padrões diferentes também receberam o mesmo peso entre si. Dentro de cada grupo, os pares dividiram esse peso igualmente. Os **{percent(a.A_balanceada)}** vieram dessa reponderação, sem treinamento ou inclusão de dados.

A WCO abaixo de 50% indicou uma ordenação contrária à expectativa de semelhança na direção prevista pelo protocolo. Isso não demonstrou ausência de informação nem inferioridade em qualquer possível classificador. Inverter a direção depois de observar o resultado mudaria a regra de avaliação.

## 7. Alcance e limites da evidência

**Houve uma tendência interna útil, com sobreposição e instabilidade.** Os pontos que delimitam essa conclusão são:

- **O resultado dependeu do trecho.** Em segmentos de 3.171 px, A foi {percent(ranks.loc['inicio3171', 'A_pares'])} no início, {percent(ranks.loc['centro3171', 'A_pares'])} no centro e {percent(ranks.loc['fim3171', 'A_pares'])} no fim. Esse controle avaliou os perfis; não recalculou a WCO nesses segmentos fixos.
- **A vantagem não apareceu em todos os padrões.** Em golfadas, A por par foi {percent(slug.A_pares)} para os perfis e {percent(slug_wco.A_pares)} para a WCO. Esse resultado local compara os pares internos de golfadas apenas aos pares entre padrões que envolvem golfadas.
- **A quantidade de aquisições foi desigual.** Bolhas dispersas tiveram somente uma aquisição, sem par interno para verificar repetibilidade. Ao retirar essa aquisição e recalcular os pesos, A balanceada caiu de {percent(a.A_balanceada)} para {percent(bubble.A_balanceada)}. As combinações que envolviam bolhas respondiam por 40% do peso dos controles balanceados.
- **O resultado dependeu das aquisições presentes.** Retirando uma aquisição por vez, A por par variou de {percent(deletion.A_pares.min())} a {percent(deletion.A_pares.max())}. Essa faixa mede influência dos dados; não é intervalo de confiança nem teste com dados novos.
- **Padrão e condições de coleta não foram separados.** Todos os panoramas anulares foram mais longos que os dos demais padrões. Comprimento, morfologia e condições operacionais podem contribuir para as diferenças observadas. As emendas dos panoramas e os deslocamentos estimados também permaneceram limitações.

**Conclusão para comunicar o resultado:** neste conjunto, a comparação do formato dos perfis colocou pares do mesmo padrão à frente de pares de padrões diferentes com maior frequência que o escore WCO, na direção definida pelo protocolo. A vantagem permaneceu ao utilizar os mesmos trechos de cada par. Os perfis foram úteis como descritores de semelhança, mas a sobreposição e a variação entre aquisições impediram concluir identificação automática ou desempenho em novos experimentos.

O ganho também não deve ser atribuído à correção numérica da v16: ao comparar perfis com a implementação anterior v5 e as mesmas escalas, A por par foi {percent(ranks.loc['v5_mesmas_escalas', 'A_pares'])}, praticamente o mesmo resultado. A principal contribuição desta etapa foi a mudança da forma de comparar as aquisições.

## Complemento 1. Vocabulário e método

| Termo usado nos outros documentos | Significado nesta análise |
| --- | --- |
| Intrapadrão / entre padrões | Duas aquisições do mesmo padrão / de padrões diferentes. |
| Escala ou comprimento de onda equivalente, λ | Tamanho espacial de uma oscilação analisada, em pixels; não mede diretamente o comprimento de uma bolha. |
| Autocorrelação | Comparação do sinal com cópias deslocadas dele próprio. |
| Lag | Deslocamento entre sinais na correlação, em pixels; não é a posição dentro do panorama. |
| CWT | Transformada wavelet contínua: examina a resposta em diferentes escalas e deslocamentos. |
| WPS alternativa | Nesta formulação, módulo da CWT da autocorrelação. O módulo não foi elevado ao quadrado. O perfil principal utilizou sua coluna de lag zero. |
| WCS alternativa | CWT da correlação cruzada entre os dois sinais; conservou valores complexos. |
| WCO alternativa | Módulo da WCS ao quadrado, dividido pelo produto das duas WPS no lag zero. O escore anterior resumiu máximos por escala com pesos de WPS. |
| Normalização por área | Reajuste da curva para que a área em log(λ) seja um. Permitiu comparar proporções de resposta, sem tratar essas proporções como energia física. |
| Prefixo comum do par | Trecho inicial de cada sinal, com ambos recortados ao comprimento do menor registro do par. |
| COI / suporte geométrico | Regra de distância às bordas da correlação. Limitou a faixa principal, sem garantir estacionariedade ou ausência de efeitos de borda. |
| Sensibilidade | Repetição da avaliação com uma mudança específica, como posição do trecho, para observar a estabilidade do resultado. |

A implementação utilizou a formulação alternativa documentada no [relatório técnico](REPORT.md) e no [protocolo](PROTOCOL.md). Aplicou-se Morlet com parâmetro 6 e correlação dividida por N. Essa divisão foi uma adaptação declarada em relação ao estimador do artigo. Esses detalhes pertencem à reprodução do método; o raciocínio da leitura permanece imagem, perfil, distância e ordenação.

## Complemento 2. Leitura dos demais gráficos

| Gráfico original | Elemento a observar | Interpretação e cuidado |
| --- | --- | --- |
| [01: perfis dos cinco padrões](figures/01_perfis_cinco_padroes.png) | Linhas finas são aquisições; linha grossa é a mediana em cada escala; faixa cobre mínimo e máximo. | Examina a variação dentro de cada padrão. A mediana é um resumo, não o perfil de uma aquisição real. A faixa não é intervalo de confiança. |
| [02: comparação das referências](figures/02_comparacao_referencias.png) | Barras mostram A; traços pretos mostram a variação ao retirar uma aquisição. | Amplia a figura 5 com versões históricas e balanceamento. Os traços não são incerteza estatística estimada. |
| [03: matriz de distâncias](figures/03_matriz_distancias.png) | Cruze a linha de uma aquisição com a coluna de outra e consulte a escala de cores. | Cada célula contém H; cores claras indicam menor distância nessa paleta. A diagonal é zero porque compara cada aquisição consigo mesma; foi excluída dos 325 pares. |
| [04: sensibilidades](figures/04_sensibilidades.png) | Cada linha muda uma escolha da análise; círculo e quadrado usam pesos diferentes. | Avalia estabilidade. Variantes com faixas diferentes não são a mesma medida; nenhuma substituiu o método principal por ter resultado maior. |
| [05: resultado por padrão](figures/05_resultados_por_padrao.png) | Compare barras dentro do mesmo padrão. | Usa A local com controles balanceados. A ausência de barra para bolhas significa falta de pares internos, e não resultado zero. |
| [06: comprimento e distância](figures/06_comprimento_e_confundimento.png) | À esquerda, cada ponto é um par; à direita, cada ponto é uma aquisição. | Mostra associação com comprimento; não estabelece causa. |
| [07: panorama, sinal e perfis](figures/07_panorama_sinal_perfis.png) | Leia cada linha da imagem para o sinal e depois para o perfil. | Apresenta um exemplo predeterminado de cada padrão. A curva cinza é uma média em lags, diferente do perfil principal. |
| [08: casos desfavoráveis](figures/08_casos_desfavoraveis.png) | Compare as duas curvas de cada painel. | Reúne os mesmos casos explicados na figura 2 e informa também a WCO. H e escore WCO têm escalas e sentidos diferentes. |
| [09: controles sintéticos](figures/09_controles_sinteticos.png) | Compare sinais construídos com escalas ou fases diferentes. | Verifica o comportamento matemático; não acrescenta experimentos reais. |
| [Mapas WCS/WCO: exemplo de bolhas alongadas](figures/par_m05q01_001__m06q01_001.png) | Eixo horizontal: deslocamento; vertical: escala; cor: valor indicado na barra de cores. | Uma região clara indica maior resposta naquela escala e deslocamento, não a localização de uma bolha. A linha tracejada marca o limite geométrico. O máximo do mapa não é o escore ponderado nem A. |

## Complemento 3. Fontes e reprodução

Os valores e as novas figuras foram extraídos das tabelas já executadas: [perfis](results/profiles_normalized.csv), [todos os pares](results/all_pairs.csv), [ordenação](results/ranking_summary.csv), [resultados por padrão](results/ranking_by_pattern.csv) e [casos ilustrados](results/selected_cases.csv). As {integer(outcomes['total'])} comparações apenas tornam explícita a conta de A já utilizada, sem criar uma avaliação adicional.

O [registro das figuras didáticas](results/reading_guide.json) contém as contagens, as faixas do histograma e os hashes das fontes. Há versões PNG e PDF das novas figuras em [figures/leitura](figures/leitura). A versão [HTML](RESUMO_EXECUTIVO.html) incorpora as cinco figuras principais para leitura sem conexão; links para documentos complementares requerem a pasta original.

Para aprofundar: [comparação entre v16 e v16.1](conclusão.md), [relatório técnico](REPORT.md), [histórico da execução](PASSO_A_PASSO.md), [decisões](DECISIONS.md) e [comandos de reprodução](README.md). As conferências computacionais estão em [final_verification.json](results/final_verification.json); elas não substituem validação experimental independente.
'''


def render_figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, PercentFormatter

    primary, pairs, ranks, cases, counts, outcomes = evidence()
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'axes.titlesize': 12, 'axes.spines.top': False,
                         'axes.spines.right': False, 'figure.facecolor': 'white',
                         'savefig.facecolor': 'white', 'pdf.fonttype': 42})
    decimal = FuncFormatter(lambda value, _: f'{value:g}'.replace('.', ','))
    manifest = dict(scope='Apresentação dos resultados existentes; protocolo e cálculos preservados.',
                    figures=[], order_comparisons=outcomes,
                    acquisitions_by_pattern={p: int(counts[p]) for p in PATTERNS})

    def save(fig, name, **detail):
        paths = []
        for extension in ['png', 'pdf']:
            path = OUT / f'{name}.{extension}'
            fig.savefig(path, dpi=180, bbox_inches='tight')
            paths.append(str(path.relative_to(ROOT)))
        manifest['figures'].append(dict(name=name, files=paths, **detail))
        plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(10, 10.5), sharex=True, sharey=True,
                             layout='constrained')
    titles = ['Mesmo padrão, curvas próximas', 'Padrões diferentes, curvas próximas',
              'Mesmo padrão, curvas mais diferentes']
    for ax, case, title in zip(axes, cases.itertuples(), titles):
        for index, run in enumerate([case.aquisicao_x, case.aquisicao_y]):
            g = primary[primary.aquisicao == run].sort_values('comprimento_onda_px')
            ax.plot(g.comprimento_onda_px, g.densidade_log, color=[TEAL, ORANGE][index],
                    ls=['-', '--'][index], lw=2.3,
                    label=f'{g.padrao.iloc[0]} · {run}')
        ax.set_title(f'{title}     |     H = {number(case.H)}', loc='left', pad=10)
        ax.set_xscale('log')
        ax.set_xticks([50, 100, 200, 500, 1000, 2000])
        ax.xaxis.set_major_formatter(decimal)
        ax.yaxis.set_major_formatter(decimal)
        ax.set_ylim(0, 1.3)
        ax.set_ylabel('Resposta relativa\n(densidade em log da escala)')
        ax.legend(loc='upper right', fontsize=9, frameon=False)
        ax.grid(alpha=.15)
    axes[-1].set_xlabel('Escala espacial da variação (px) · eixo logarítmico')
    fig.suptitle('A distância compara as curvas; os rótulos podem coincidir ou diferir\n'
                 'H menor = curvas mais parecidas · exemplos escolhidos após os resultados', fontsize=13)
    save(fig, '01_pares_explicados', cases=cases.to_dict(orient='records'),
         variant='principal', domain='principal')

    # Mesmas bordas para ambos os grupos, sem suavização e sem aparar observações.
    bin_width = .025
    maximum = float(pairs.valor.max())
    edges = np.arange(0, (int(np.ceil(maximum / bin_width)) + 1) * bin_width, bin_width)
    if edges[-1] < maximum:
        edges = np.append(edges, edges[-1] + bin_width)
    fig, ax = plt.subplots(figsize=(10, 5.2), layout='constrained')
    histogram = []
    for same, color, label, linestyle in [
        (True, TEAL, 'Mesmo padrão', '-'), (False, ORANGE, 'Padrões diferentes', '--')]:
        values = pairs[pairs.mesmo_padrao == same].valor.to_numpy()
        observed, _ = np.histogram(values, bins=edges)
        assert observed.sum() == len(values)
        proportions = observed / len(values)
        ax.stairs(proportions, edges, fill=True, color=color, alpha=.13)
        ax.stairs(proportions, edges, color=color, lw=2.2, ls=linestyle,
                  label=f'{label}: {len(values)} pares · média H = {number(values.mean())}')
        ax.axvline(values.mean(), color=color, lw=1.5, ls=':')
        histogram.append(dict(same_pattern=same, n=len(values), counts=observed.tolist(),
                              proportions=proportions.tolist(), mean=float(values.mean())))
    ax.set_xlim(edges[0], edges[-1])
    ax.set_ylim(0, max(max(group['proportions']) for group in histogram) * 1.45)
    ax.xaxis.set_major_formatter(decimal)
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
    ax.set_xlabel('Distância H · menor = curvas mais parecidas\nIntervalo observado; H pode variar de 0 a 1')
    ax.set_ylabel('Porcentagem dos pares de cada grupo\nem cada faixa de distância')
    ax.legend(frameon=False, loc='upper right', fontsize=10)
    ax.grid(axis='y', alpha=.15)
    ax.set_title('Os grupos se sobrepõem: uma distância isolada não determina o padrão\n'
                 'Cada grupo soma 100% · linhas pontilhadas indicam médias', pad=16)
    save(fig, '02_distribuicao_distancias', bin_edges=edges.tolist(), groups=histogram)

    fig = plt.figure(figsize=(10, 5.7), layout='constrained')
    grid = fig.add_gridspec(2, 1, height_ratios=[1.2, 1])
    top = fig.add_subplot(grid[0])
    top.set_axis_off()
    top.text(.22, .65, '87 pares\ndo mesmo padrão', ha='center', va='center', fontsize=14,
             bbox=dict(boxstyle='round,pad=.8', facecolor='#e1f0f0', edgecolor=TEAL))
    top.text(.5, .65, '×', ha='center', va='center', fontsize=26)
    top.text(.78, .65, '238 pares\nde padrões diferentes', ha='center', va='center', fontsize=14,
             bbox=dict(boxstyle='round,pad=.8', facecolor='#faeedf', edgecolor=ORANGE))
    top.text(.5, .10, f'= {integer(outcomes["total"])} comparações entre pares',
             ha='center', va='center', fontsize=15, fontweight='bold')
    bottom = fig.add_subplot(grid[1])
    fractions = [outcomes[key] / outcomes['total'] for key in ['favorable', 'unfavorable', 'ties']]
    left = 0.
    for fraction, color in zip(fractions, [TEAL, '#bfc6cb', ORANGE]):
        if fraction:
            bottom.barh(0, fraction, left=left, color=color, height=.45)
        left += fraction
    bottom.text(fractions[0] / 2, 0,
                f'{integer(outcomes["favorable"])}: ordem esperada\n{percent(fractions[0])}',
                ha='center', va='center', color='white', fontsize=13, fontweight='bold')
    bottom.text(fractions[0] + fractions[1] / 2, 0,
                f'{integer(outcomes["unfavorable"])}: ordem contrária\n{percent(fractions[1])}',
                ha='center', va='center', fontsize=10)
    bottom.set_xlim(0, 1)
    bottom.set_ylim(-.45, .4)
    bottom.set_yticks([])
    bottom.set_xticks(np.linspace(0, 1, 6))
    bottom.xaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
    bottom.spines['left'].set_visible(False)
    bottom.set_xlabel('Ordem esperada: o par do mesmo padrão teve a menor distância H\n'
                      f'{outcomes["ties"]} empates · comparações reutilizam as mesmas aquisições')
    fig.suptitle(f'De onde vieram os {percent(ranks.loc["wps_principal", "A_pares"])}\n'
                 'Contagem de ordenações entre pares, não de aquisições classificadas', fontsize=14)
    save(fig, '03_calculo_concordancia', counts=outcomes)

    methods = ['wco_v16_18_71', 'wps_prefixo_par', 'wps_principal']
    labels = ['WCO anterior\nTrecho inicial comum · 50–3.000 px',
              'Perfis no mesmo trecho\nTrecho inicial comum · 50–3.000 px',
              'Perfis principais\nRegistro inteiro · 50–aprox. 2.296 px']
    fig, ax = plt.subplots(figsize=(10, 4.5), layout='constrained')
    values = ranks.loc[methods, 'A_pares'].to_numpy()
    ax.barh(range(3), values, height=.48, color=[PURPLE, '#71aaa7', TEAL])
    for i, value in enumerate(values):
        ax.text(value + .015, i, percent(value), va='center', fontsize=13, fontweight='bold')
    ax.set_yticks(range(3), labels, fontsize=10)
    ax.set_ylim(2.6, -.7)
    ax.set_xlim(0, 1)
    ax.set_xticks(np.linspace(0, 1, 6))
    ax.xaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=0))
    ax.axvline(.5, color='#59656b', ls='--', lw=1)
    ax.grid(axis='x', alpha=.15)
    ax.set_xlabel('Concordância A · cada par com o mesmo peso\n'
                  '50%: sem preferência de ordenação · A não é acurácia')
    ax.set_title('Os perfis favoreceram a ordem esperada, inclusive no mesmo trecho\n'
                 'A barra intermediária controla o trecho e a faixa de escalas', pad=16)
    save(fig, '04_comparacao_metodos', values=[dict(method=m, A_pares=float(v))
                                             for m, v in zip(methods, values)])

    source_paths = [RESULTS / name for name in SOURCES] + [ROOT / 'config.json', Path(__file__)]
    manifest['source_sha256'] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in source_paths}
    manifest['output_sha256'] = {
        path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
        for item in manifest['figures'] for path in item['files']}
    (RESULTS / 'reading_guide.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
    print('Quatro figuras didáticas em PNG/PDF; contagens e fontes registradas.', flush=True)


if __name__ == '__main__':
    render_figures()
