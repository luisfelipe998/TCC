"""Resumo executivo gerado exclusivamente a partir dos resultados calculados."""
from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FIGURES = HERE/'figures'


def table(frame,digits=4):
    def cell(value):
        if pd.isna(value):return 'não estimável'
        if isinstance(value,(float,np.floating)):return f'{value:.{digits}f}'
        return str(value)
    lines=['| '+' | '.join(map(str,frame.columns))+' |','| '+' | '.join(['---']*len(frame.columns))+' |']
    lines += ['| '+' | '.join(cell(v) for v in row)+' |' for row in frame.itertuples(index=False,name=None)]
    return chr(10).join(lines)


def executive_pdf(index,pairs,aggregate,s,v):
    """PDF independente com resumo, gráfico, fórmulas e verificação."""
    def page(title,number):
        fig=plt.figure(figsize=(8.27,11.69),facecolor='white')
        fig.text(.075,.955,title,fontsize=19,weight='bold',va='top')
        fig.text(.075,.927,'Análise de correlação entre descritores espectrais',fontsize=9,color='#555555')
        fig.text(.925,.025,str(number),ha='right',fontsize=9,color='#555555')
        return fig
    def paragraph(fig,text,y,width=97,size=10.4):
        wrapped=chr(10).join(textwrap.wrap(text,width=width))
        fig.text(.075,y,wrapped,fontsize=size,va='top',linespacing=1.45)
        return y-.020*len(wrapped.splitlines())-.025
    pdf_path=ROOT/'RESUMO_EXECUTIVO.pdf'
    with PdfPages(pdf_path) as pdf:
        fig=page('Resumo executivo',1)
        y=paragraph(fig,'A análise compreendeu 26 ensaios e 325 pares distintos, dos quais 87 '
            'intrapadrão e 238 entre padrões. A similaridade dos descritores espectrais foi '
            'quantificada pelo coeficiente de correlação de Pearson, no intervalo de −1 a 1.',.885)
        y=paragraph(fig,f'A mediana de r foi {s["median_intra"]:.4f} dentro dos padrões e '
            f'{s["median_inter"]:.4f} entre padrões. Dos 238 pares entre padrões, '
            f'{s["inter_pairs_above_intra_median"]} atingiram ou superaram a mediana intrapadrão. '
            'Há sobreposição: similaridade alta de perfil não identifica, por si só, o padrão de escoamento.',y)
        y=paragraph(fig,f'O r observado variou de {s["pearson_min"]:.4f} a {s["pearson_max"]:.4f}; '
            f'{s["negative_pairs"]} pares tiveram r negativo. Os coeficientes descrevem a forma de '
            'perfis WPS residuais, não correlação entre intensidades brutas, acoplamento físico ou acurácia.',y)
        y=paragraph(fig,'Cada ensaio foi representado por um descritor obtido dos primeiros 3.171 px, '
            'correspondentes ao menor comprimento disponível. Esse descritor foi mantido fixo em '
            'todas as comparações, assegurando a uniformidade do comprimento analisado.',y)
        ax=fig.add_axes([.06,.18,.89,.27]);ax.imshow(plt.imread(FIGURES/'correlacoes_dentro_entre.png'));ax.axis('off')
        paragraph(fig,'Gráfico descritivo: os pares compartilham ensaios. Não foram calculados '
            'testes de separação ou classificadores. A classe dispersed bubbles tem somente um registro.',.13,size=9.7)
        pdf.savefig(fig);plt.close(fig)
        fig=page('Matriz de correlação de Pearson',2)
        ax=fig.add_axes([.015,.13,.97,.77]);ax.imshow(plt.imread(FIGURES/'matriz_pearson_todos.png'));ax.axis('off')
        paragraph(fig,'Coeficientes de correlação entre os descritores espectrais dos 26 ensaios, '
            'calculados em segmentos de 3.171 px. A identificação dos registros é apresentada '
            'ao final do relatório.',.10,size=9.4)
        pdf.savefig(fig);plt.close(fig)
        fig=page('Passo a passo e fórmulas',3)
        y=.88
        steps=[
            ('1. Extração e trecho fixo','Média das 54 linhas 18–71 dos PNGs; primeiros N=3.171 pixels de cada ensaio.',
             r'$s_i[n]=\frac{1}{54}\sum_{h=18}^{71} I_i[h,n]$'),
            ('2. Remoção de tendência e padronização','Inverter intensidade, retirar uma reta por mínimos quadrados e padronizar dentro do trecho.',
             r'$v_i=255-s_i,\quad e_i=v_i-(\hat\alpha_i+\hat\beta_i n),\quad x_i=\frac{e_i-\bar e_i}{\mathrm{std}_N(e_i)}$'),
            ('3. Autocorrelação enviesada','Divisor N, como na versão principal v16; extensão bilateral simétrica para sinais reais.',
             r'$\hat R_{ii}[m]=\frac{1}{N}\sum_{n=0}^{N-1-m}x_i[n+m]x_i[n]$'),
            ('4. Alternative WPS no lag zero','Morlet com frequência angular 6; 200 comprimentos de onda de 50 a 3.000 px.',
             r'$P_i(a)=|\mathrm{CWT}\{\hat R_{ii}\}(0,a)|,\quad \lambda=\frac{4\pi a}{6+\sqrt{38}}$'),
            ('5. Descritor de forma','Retirar a reta ajustada em log(P) versus log(lambda). Isso não garante branqueamento estatístico.',
             r'$z_{ik}=\log(P_i(a_k)+10^{-12})-\hat c_i-\hat d_i\log\lambda_k$'),
            ('6. Pearson direto','Usar o mesmo vetor de cada ensaio em todas as comparações, preservando sinal e escala de r.',
             r'$r_{ij}=\frac{\sum_k(z_{ik}-\bar z_i)(z_{jk}-\bar z_j)}{\sqrt{\sum_k(z_{ik}-\bar z_i)^2}\sqrt{\sum_k(z_{jk}-\bar z_j)^2}}$')]
        for title,description,formula in steps:
            fig.text(.075,y,title,fontsize=11.3,weight='bold',va='top')
            fig.text(.075,y-.026,chr(10).join(textwrap.wrap(description,95)),fontsize=9.2,va='top',linespacing=1.35)
            fig.text(.095,y-.079,formula,fontsize=13.5,va='top')
            y-=.132
        pdf.savefig(fig);plt.close(fig)
        fig=page('Checagem e limites da interpretação',4)
        scipy_error=max(c['scipy_max_abs_error'] for c in v['matrix_checks'].values())
        y=paragraph(fig,f'Validação: PASS. Conferidos 26 PNGs, 26 vetores de entrada, 26 perfis '
            'contra a FFT e 52 perfis contra tabelas anteriores. Todas as 325 comparações de cada '
            f'condição coincidiram com SciPy; maior diferença absoluta {scipy_error:.3e}.',.88)
        y=paragraph(fig,'As matrizes são simétricas, possuem diagonal unitária, 676 células finitas '
            'e autovalores não negativos dentro da tolerância numérica. Foram conferidos casos '
            'conhecidos r=−1,0,1, vetores constantes, pré-processamento, autocorrelação e resíduos.',y)
        y=paragraph(fig,f'No controle com registros completos, a alteração absoluta mediana dos '
            f'coeficientes foi {s["median_abs_change_with_full"]:.4f}, e a maior foi '
            f'{s["max_abs_change_with_full"]:.4f}. O controle mantém um descritor por ensaio, '
            'mas usa comprimentos diferentes. Não foi escolhido o resultado mais favorável.',y)
        ax=fig.add_axes([.21,.205,.60,.36]);ax.imshow(plt.imread(FIGURES/'sensibilidade_comprimento.png'));ax.axis('off')
        paragraph(fig,'As maiores escalas se aproximam do comprimento principal de 3.171 px. '
            'Registro finito, reconstrução espacial e condições operacionais limitam atribuições '
            'físicas. Correção matemática não demonstra generalização, causalidade ou classificação.',.155,size=9.7)
        pdf.savefig(fig);plt.close(fig)
        fig=page('Identificação dos 26 ensaios',5)
        ax=fig.add_axes([.075,.12,.85,.76]);ax.axis('off')
        grid=ax.table(cellText=index[['number','run','pattern','n']].values,
            colLabels=['Nº','Ensaio','Padrão','N original (px)'],loc='center',cellLoc='left',
            colWidths=[.08,.32,.38,.22])
        grid.auto_set_font_size(False);grid.set_fontsize(9);grid.scale(1,1.68)
        for (row,col),cell in grid.get_celld().items():
            cell.set_edgecolor('#dddddd')
            if row==0:cell.set_facecolor('#edf2f6');cell.get_text().set_weight('bold')
        paragraph(fig,'Dados completos, scripts, fórmulas e hashes: '
            'src2/estatistica5/todos_contra_todos/. O método detalhado está em '
            'METODOLOGIA_TODOS_CONTRA_TODOS.md.',.10,size=9.7)
        pdf.savefig(fig);plt.close(fig)


def write_reports(index,pairs,aggregate,groups,s,v):
    main=aggregate[aggregate.condition=='3171_px_inicio'].copy()
    columns={'comparison':'Comparação','n_pairs':'Pares','minimum':'r mínimo',
             'median':'r mediano','maximum':'r máximo','negative_pairs':'r negativos'}
    main=main[list(columns)].rename(columns=columns)
    within=groups[groups.pattern_x==groups.pattern_y][['pattern_x','n_pairs','minimum','median','maximum']]
    within=within.rename(columns={'pattern_x':'Padrão','n_pairs':'Pares','minimum':'r mínimo',
                                'median':'r mediano','maximum':'r máximo'})
    strongest_inter=pairs[~pairs.same_pattern].nlargest(5,'pearson_r')[['run_x','pattern_x','run_y','pattern_y','pearson_r']]
    negative=pairs[pairs.pearson_r<0]
    negative_touching_single=int(((negative.pattern_x=='dispersed bubbles') | (negative.pattern_y=='dispersed bubbles')).sum())
    max_error=max(c['scipy_max_abs_error'] for c in v['matrix_checks'].values())
    report=fr'''# Resumo executivo — correlação entre descritores espectrais

**A análise contemplou 26 ensaios e 325 pares distintos, cuja similaridade espectral foi quantificada pelo coeficiente de correlação de Pearson.** A mediana foi {s['median_intra']:.4f} dentro dos padrões e {s['median_inter']:.4f} entre padrões. Observou-se sobreposição entre as distribuições: {s['inter_pairs_above_intra_median']} dos 238 pares entre padrões atingiram ou superaram a mediana intrapadrão. Esses resultados indicam que a correlação elevada entre perfis, isoladamente, não caracteriza o padrão de escoamento.

[Resumo em PDF](RESUMO_EXECUTIVO.pdf) · [Fórmulas e método detalhado](METODOLOGIA_TODOS_CONTRA_TODOS.md) · [Matriz em PDF](todos_contra_todos/figures/matriz_pearson_todos.pdf)

## Delimitação da amostra e padronização dos descritores

A amostra compreendeu 87 pares intrapadrão e 238 pares entre padrões. A diagonal da matriz representa a correlação de cada descritor consigo mesmo e foi excluída dos resumos estatísticos. Devido à simetria da matriz, cada par foi contabilizado uma única vez.

O ensaio m07q01_001 foi mantido na classe dispersed bubbles e comparado aos demais registros. Como essa classe contém uma única aquisição, sua similaridade intraclasse não pôde ser estimada. Foram preservados os critérios de exclusão dos cinco registros não elegíveis da base original.

Cada ensaio foi representado por um descritor obtido dos primeiros **3.171 px**, correspondentes ao menor comprimento disponível. O mesmo vetor foi utilizado em todas as comparações envolvendo o respectivo ensaio.

## Matriz de correlação de Pearson

![Matriz de correlação de Pearson entre descritores espectrais](todos_contra_todos/figures/matriz_pearson_todos.png)

<!-- Nota de apresentação: a matriz contém 325 pares distintos, sendo 87 intrapadrão e 238 entre padrões. As 676 células foram calculadas; a diagonal r=1 expressa identidade. Na paleta utilizada, roxo indica correlação negativa, creme indica valor próximo de zero e azul indica correlação positiva. Essas informações foram retiradas do rodapé da figura para simplificar sua apresentação. Na versão anterior, as regiões cinzas correspondiam a comparações omitidas e os trechos dependiam do par; a análise atual utiliza um descritor fixo por ensaio. -->

- **r próximo de 1:** forma semelhante, admitindo mudança de amplitude e deslocamento vertical do descritor.
- **r próximo de 0:** pouca associação linear entre as coordenadas espectrais.
- **r próximo de −1:** formas linearmente opostas.

A escala cromática foi fixada no intervalo de −1 a 1. Os coeficientes representam a associação linear entre os **resíduos log-lineares da alternative WPS no lag zero**; sua interpretação se restringe à forma desses descritores espectrais.

## Resultados observados

{table(main)}

A amplitude observada foi de {s['pearson_min']:.6f} a {s['pearson_max']:.6f}; {s['negative_pairs']} pares apresentaram r negativo. Desses, {negative_touching_single} envolveram o único ensaio dispersed bubbles. Isso descreve esse registro, sem demonstrar o comportamento geral dessa classe.

![Correlações dentro e entre padrões](todos_contra_todos/figures/correlacoes_dentro_entre.png)

Os pares compartilham ensaios; as 325 correlações não são 325 aquisições independentes. Pontos, medianas e intervalos mínimo–máximo são descritivos e não incluem inferência de significância.

Dentro de cada padrão:

{table(within)}

Dispersed bubbles tem apenas um registro: sua similaridade intraclasse não pode ser estimada. Plug e slug possuem três registros cada. As contagens e a incerteza de representatividade impedem generalizações fortes a partir das diferenças visuais.

Os cinco maiores coeficientes entre padrões diferentes foram:

{table(strongest_inter)}

Esses pares ilustram por que uma matriz com somente blocos intrapadrão não bastava para avaliar a especificidade da semelhança. O destaque é descritivo; os pares foram ordenados após calcular o conjunto completo e nenhum foi removido.

## O que foi executado, passo a passo

1. Conferiu-se a média das linhas 18–71 diretamente nos 26 PNGs e sua igualdade com os sinais salvos.
2. Selecionaram-se os primeiros 3.171 px de cada ensaio, sem reamostragem.
3. Inverteu-se a intensidade, retirou-se a tendência linear e padronizou-se cada trecho, como na v16.
4. Calculou-se a autocorrelação bilateral enviesada, com divisor N.
5. Aplicou-se a Morlet, com frequência angular 6, em 200 comprimentos de onda de 50 a 3.000 px. Retirou-se a magnitude da CWT da autocorrelação no lag zero.
6. Ajustou-se log(P+10⁻¹²) em função de log(comprimento de onda) e usaram-se os resíduos como descritor de forma.
7. Calculou-se Pearson assinado para todos os pares, com um único descritor por ensaio.
8. Repetiu-se somente Pearson com os registros completos como controle de sensibilidade; os dois resultados são apresentados sem selecionar o mais favorável.

As definições centrais são:

$$P_i(a)=|\mathrm{{CWT}}\{{\widehat R_{{ii}}\}}(0,a)|,\qquad z_{{ik}}=\log(P_i(a_k)+10^{{-12}})-\widehat c_i-\widehat d_i\log\lambda_k,$$

$$r_{{ij}}=\frac{{\sum_k(z_{{ik}}-\overline z_i)(z_{{jk}}-\overline z_j)}}{{\sqrt{{\sum_k(z_{{ik}}-\overline z_i)^2}}\sqrt{{\sum_k(z_{{jk}}-\overline z_j)^2}}}}.$$

O [método detalhado](METODOLOGIA_TODOS_CONTRA_TODOS.md) apresenta também as fórmulas da extração, tendência, normalização, autocorrelação e Morlet. O resíduo log-linear não constitui prova de branqueamento estatístico nem de independência entre escalas.

![Perfis e descritores usados](todos_contra_todos/figures/perfis_comparados.png)

## Sensibilidade ao comprimento

{table(aggregate[['condition','comparison','n_pairs','median','minimum','maximum']])}

Ao usar os registros completos, a alteração absoluta mediana de r foi {s['median_abs_change_with_full']:.4f}; a maior alteração foi {s['max_abs_change_with_full']:.4f}. O controle também usa um descritor fixo por ensaio, porém os comprimentos originais diferem. Mudanças refletem comprimento e conteúdo espacial incluído, sem permitir atribuir toda variação a uma única causa.

![Comparação entre trecho comum e registro completo](todos_contra_todos/figures/sensibilidade_comprimento.png)

A [matriz com registros completos](todos_contra_todos/figures/matriz_pearson_completo.png) foi conservada como controle. A janela principal foi definida pela disponibilidade dos dados, não pelo maior contraste entre padrões. Como 3.000 px se aproxima do comprimento de 3.171 px, as maiores escalas exigem cautela quanto ao registro finito.

## Checagem do método e dos cálculos

**Status: PASS.** Todos os 26 sinais coincidiram com os pixels de origem. As somas temporais da wavelet foram confrontadas com a implementação FFT em todas as 200 escalas de cada ensaio; a maior diferença relativa L2 foi {v['max_profile_fft_relative_l2']:.3e}. Os 52 perfis das duas condições também reproduziram as tabelas da v16.

Os 325 coeficientes de cada condição foram conferidos contra SciPy, com diferença absoluta máxima de {max_error:.3e}, além de NumPy. Foram verificados casos conhecidos r=−1,0,1, recusa de vetores constantes, ordem do pré-processamento, autocorrelação por soma literal, conversão de escala e ortogonalidade dos resíduos. As matrizes são simétricas, têm diagonal unitária, 676 células finitas e autovalores não negativos dentro da tolerância de arredondamento, como exige R=UUᵀ.

[Validação completa](todos_contra_todos/results/validation.json) · [Verificação por ensaio](todos_contra_todos/results/verificacao_por_ensaio.csv) · [Hashes e versões](todos_contra_todos/results/manifest.json)

As fontes permaneceram inalteradas. A reconstrução dos panoramas não foi refeita. Correção matemática não demonstra causalidade, representatividade da amostra ou desempenho de classificação. Não foram calculados testes de separação ou classificadores; Pearson não foi usado para gerar valores-p tomando as escalas como aquisições independentes.

## Identificação dos ensaios

{table(index[['number','run','pattern','n']].rename(columns={'number':'Nº','run':'Ensaio','pattern':'Padrão','n':'N original (px)'}))}
'''
    (ROOT/'RESUMO_EXECUTIVO.md').write_text(report)
    executive_pdf(index,pairs,aggregate,s,v)
    readme=f'''# estatistica5 — Pearson direto, todos contra todos

**Versão atual:** 26 ensaios, 325 pares únicos, escala de Pearson de −1 a 1, sem células não calculadas. Mediana r intrapadrão: {s['median_intra']:.4f}; entre padrões: {s['median_inter']:.4f}. As distribuições apresentam sobreposição; correlação alta de perfil não identifica isoladamente o padrão.

- [Resumo executivo com gráficos](RESUMO_EXECUTIVO.md) · [PDF para compartilhar](RESUMO_EXECUTIVO.pdf)
- [Passo a passo e fórmulas](METODOLOGIA_TODOS_CONTRA_TODOS.md)
- [Matriz principal PNG](todos_contra_todos/figures/matriz_pearson_todos.png) · [PDF vetorial](todos_contra_todos/figures/matriz_pearson_todos.pdf)
- [325 coeficientes](todos_contra_todos/results/pearson_325_pares.csv) · [Matriz CSV](todos_contra_todos/results/matriz_pearson.csv)
- [Índice dos ensaios](todos_contra_todos/results/indice_ensaios.csv) · [Descritores](todos_contra_todos/results/descritores_por_ensaio.npz)
- [Validação](todos_contra_todos/results/validation.json) · [Manifesto](todos_contra_todos/results/manifest.json)

## Reproduzir a versão atual

Na raiz do repositório:

    PYTHONDONTWRITEBYTECODE=1 python3 src2/estatistica5/todos_contra_todos/run_all_pairs.py

A execução confere as fontes e gera todos os resultados, cinco gráficos em PNG/PDF, o resumo executivo em Markdown/PDF e o manifesto. A condição principal usa os primeiros 3.171 px de cada registro; um descritor fixo por ensaio permite a comparação integral. O controle com registros completos é apresentado separadamente.

## Versão anterior

A [análise dos 87 pares](README_87_PARES.md) foi preservada no [arquivo histórico completo](historico_87_pares.zip), incluindo o README original e seu manifesto conferido. Os arquivos antigos em results/ e figures/ e o script run_pearson.py também permanecem disponíveis. Nela havia cinza porque as comparações entre padrões não eram incluídas. Ela usa distância normalizada e trechos dependentes do par; não é a versão atual. O script histórico regrava o README antigo se executado; use o comando acima para reproduzir a comparação completa.

Nenhum ensaio foi reclassificado. O único dispersed bubbles agora participa das comparações com as outras classes, mantendo seu rótulo original. A [explicação da ECDF anterior](EXPLICACAO_ECDF.md) permanece histórica e não descreve os novos coeficientes.
'''
    (ROOT/'README.md').write_text(readme)
