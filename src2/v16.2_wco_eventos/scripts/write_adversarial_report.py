"""Relata a revisão adversarial e integra seus achados à entrega existente."""
from pathlib import Path
from datetime import datetime,timezone
import json
import os
import subprocess
import textwrap
os.environ.setdefault('MPLCONFIGDIR','/tmp/wco-eventos-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/auditoria_adversarial'
FIG=ROOT/'figures/auditoria_adversarial'
R='results/auditoria_adversarial'
F='figures/auditoria_adversarial'
LABEL={'plug':'Pistonado (plug)','slug':'Golfadas (slug)','churn':'Agitado (churn)','wavy annular':'Anular ondulado'}
COLOR={'plug':'#217f9a','slug':'#ce733b','churn':'#7743a7','wavy annular':'#137d69'}
START='<!-- AUDITORIA_ADVERSARIAL_INICIO -->'
END='<!-- AUDITORIA_ADVERSARIAL_FIM -->'


def f(x,d=4):return f'{x:.{d}f}'.replace('.',',')
def sci(x):return f'{x:.2e}'.replace('.',',')
def read(name):return json.loads((OUT/name).read_text())
def table(head,rows):return '\n'.join(['| '+' | '.join(head)+' |','| '+' | '.join(['---']*len(head))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rows])


def append(name,body,early=False):
    p=ROOT/name;text=p.read_text()
    if START in text:
        a=text.index(START);b=text.index(END,a)+len(END);text=text[:a]+text[b:]
    block=START+'\n'+body.strip()+'\n'+END+'\n'
    if early:
        # A primeira seção fica imediatamente depois do parágrafo de abertura.
        first=text.index('\n\n',text.index('\n\n')+2)
        text=text[:first]+'\n\n'+block+text[first:]
    else:text=text.rstrip()+'\n\n'+block
    p.write_text(text.rstrip()+'\n')


def save(fig,name):
    for ax in fig.axes:
        if ax.get_yscale()=='linear':ax.yaxis.set_major_formatter(FuncFormatter(lambda x,_:f'{x:g}'.replace('.',',')))
    fig.savefig(FIG/(name+'.png'),dpi=160,bbox_inches='tight');plt.close(fig)


def plots():
    FIG.mkdir(exist_ok=True)
    df=pd.read_csv(OUT/'repetibilidade_selecao.csv').query("tipo=='minimo'")
    fig,axes=plt.subplots(1,2,figsize=(12,4.4),layout='constrained')
    variants=['sigma12','sigma20','faixa30_49'];x=np.arange(3)
    for ax,pat in zip(axes,['churn','wavy annular']):
        g=df[df.padrao==pat].set_index('alternativa').loc[variants]
        ax.bar(x-.17,g.candidatos_em_64/g.total,width=.34,color='#b1b7bd',label='Candidato próximo do original')
        ax.bar(x+.17,g.selecionados_em_64/g.total,width=.34,color=COLOR[pat],label='Selecionado novamente por estrato')
        for j,r in enumerate(g.itertuples()):
            ax.text(j-.17,r.candidatos_em_64/r.total+.02,f'{r.candidatos_em_64}/{r.total}',ha='center',fontsize=10)
            ax.text(j+.17,r.selecionados_em_64/r.total+.02,f'{r.selecionados_em_64}/{r.total}',ha='center',fontsize=10)
        ax.set_xticks(x,['σ=12','σ=20','Linhas 30–49']);ax.set_ylim(0,1.2);ax.set_ylabel('Fração dos mínimos originais')
        ax.set_title(LABEL[pat]);ax.legend(fontsize=8,loc='lower left')
    fig.suptitle('Localizar um candidato próximo não garante selecionar a mesma âncora\nCorrespondência de mesmo tipo em até 64 px; seleção completa refeita',fontsize=13)
    save(fig,'repetibilidade_selecao')
    df=pd.read_csv(OUT/'sensibilidade_semente.csv');patterns=['plug','slug','churn','wavy annular']
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
    for ax,comp in zip(axes,['fixos','distintos']):
        samples=[];original=[];allpairs=[]
        for pat in patterns:
            control='fixos' if comp=='fixos' else 'intra_distintos' if pat in ['plug','slug'] else 'minimos_distintos'
            g=df[(df.padrao==pat)&(df.comparador==control)]
            samples.append(g[g.variante=='semente_alternativa'].diferenca)
            original.append(g[g.variante=='semente_original'].diferenca.iloc[0]);allpairs.append(g[g.variante=='todos_admissiveis'].diferenca.iloc[0])
        ax.boxplot(samples,positions=np.arange(4),widths=.45,whis=(0,100),showfliers=False,
            medianprops={'color':'#333'},boxprops={'color':'#888'})
        ax.scatter(range(4),original,marker='D',s=45,c=[COLOR[p] for p in patterns],label='Semente original',zorder=4)
        ax.scatter(np.arange(4)+.12,allpairs,marker='*',s=95,color='#c23e28',label='Todos os pares admissíveis',zorder=4)
        ax.axhline(0,color='#999',lw=1);ax.axhline(.05,color='#999',ls=':',lw=1)
        ax.set_xticks(range(4),['Plug','Slug\n1 par equivalente','Churn','Anular']);ax.set_ylabel('Âncoras − controle: WCO média')
        ax.set_title('Controle: '+('janelas fixas' if comp=='fixos' else 'larguras/extensões distintas'));ax.legend(fontsize=8)
    fig.suptitle('Escolha dos parceiros: 100 sementes alternativas, sem mudar as âncoras\nCaixa = quartis; hastes = mínimo/máximo. Não são intervalos de confiança. Linha pontilhada = 0,05.',fontsize=12)
    save(fig,'sensibilidade_parceiros')
    df=pd.read_csv(OUT/'nulo_selecao_completa.csv')
    fig,axes=plt.subplots(1,2,figsize=(11.5,4.4),layout='constrained')
    for r in df.itertuples():axes[0].plot([0,1],[r.fixos,r.ancoras],color='#bbbbbb',lw=.7,alpha=.7)
    axes[0].scatter(np.zeros(len(df)),df.fixos,color='#777',s=14);axes[0].scatter(np.ones(len(df)),df.ancoras,color='#7743a7',s=14)
    axes[0].plot([0,1],[df.fixos.mean(),df.ancoras.mean()],color='#c23e28',marker='D',lw=2,label='Média das 40 repetições')
    axes[0].set_xticks([0,1],['Janelas fixas','Mínimos ancorados']);axes[0].set_ylabel('WCO média por repetição');axes[0].legend(fontsize=8)
    axes[1].hist(df.diferenca,bins=12,color='#7b6b9a',edgecolor='white');axes[1].axvline(0,color='#333',ls=':')
    axes[1].axvline(df.diferenca.mean(),color='#c23e28',label='Média = '+f(df.diferenca.mean()))
    axes[1].set_xlabel('Âncoras − janelas fixas');axes[1].set_ylabel('Repetições');axes[1].legend(fontsize=8)
    axes[1].xaxis.set_major_formatter(FuncFormatter(lambda x,_:f'{x:g}'.replace('.',',')))
    fig.suptitle('Sinais AR(1) independentes: detecção e seleção completas antes da WCO\n40 repetições; φ=0,95. Sem evento físico comum; diagnóstico posterior, sem valor-p.',fontsize=12)
    save(fig,'nulo_selecao')


def page(pdf,title,blocks):
    fig=plt.figure(figsize=(8.27,11.69));fig.text(.08,.95,title,fontsize=22,weight='bold',va='top',color='#173e58')
    y=.83
    for head,body in blocks:
        fig.text(.08,y,head,fontsize=11.7,weight='bold',va='top',color='#137d69');y-=.03
        lines=textwrap.wrap(body,width=87);fig.text(.08,y,'\n'.join(lines),fontsize=10.1,linespacing=1.4,va='top')
        y-=.0205*len(lines)+.025
    assert y>.04,(title,y)
    fig.text(.08,.025,'v16.2_wco_eventos • Revisão adversarial • Detalhes e limites: AUDITORIA_ADVERSARIAL.md',fontsize=7.8,color='#666')
    pdf.savefig(fig);plt.close(fig)


def integrate():
    data=read('dados.json');math=read('formulas.json');diag=read('diagnosticos.json');null=diag['etapas']['synthetic']
    assert all(x['status']=='PASS' for x in [data,math,diag])
    sel=pd.read_csv(OUT/'repetibilidade_selecao.csv').query("tipo=='minimo'")
    seeds=pd.read_csv(OUT/'sensibilidade_semente_resumo.csv')
    allseeds=pd.read_csv(OUT/'sensibilidade_semente.csv')
    seams=pd.read_csv(OUT/'contrastes_sem_feather.csv')
    maxmath=max(r['erro_CWT_soma_direta'] for r in math['comparacoes_independentes'])
    counts={v:int(sel[sel.alternativa==v].selecionados_em_64.sum()) for v in ['sigma12','sigma20','faixa30_49']}
    seedrow=lambda pat,co:seeds[(seeds.padrao==pat)&(seeds.comparador==co)].iloc[0]
    ann=seedrow('wavy annular','fixos');churn=seedrow('churn','fixos')
    select_table=table(['Padrão','Correspondência local σ12 / σ20','Seleção completa σ12 / σ20','Seleção completa: linhas 30–49'],[
        [LABEL[p],f'{int(g.iloc[0].total)}/{int(g.iloc[0].total)} em ambos',
         f'{int(g.set_index("alternativa").loc["sigma12","selecionados_em_64"])}/{int(g.iloc[0].total)}; {int(g.set_index("alternativa").loc["sigma20","selecionados_em_64"])}/{int(g.iloc[0].total)}',
         f'{int(g.set_index("alternativa").loc["faixa30_49","selecionados_em_64"])}/{int(g.iloc[0].total)}'] for p,g in sel.groupby('padrao')])
    seed_table=table(['Padrão','Controle','Mínimo / máximo em 100 sementes','Sementes com diferença positiva'],[
        [LABEL[r.padrao],'Fixos' if r.comparador=='fixos' else 'Larguras/extensões distintas',f'{f(r.minimo)} / {f(r.maximo)}',f'{int(round(100*r.fracao_positiva))}/100'] for r in seeds.itertuples()])
    seam_rows=[]
    for p in ['plug','slug','churn','wavy annular']:
        g=seams[(seams.padrao==p)&(seams.comparador=='fixos')].set_index('variante')
        a,b=g.loc['original','diferenca'],g.loc['sem_feather','diferenca']
        seam_rows.append([LABEL[p],f(a),f(b),f(b-a,6)])
    seam_table=table(['Padrão','Âncoras − fixos original','Sem feather','Mudança'],seam_rows)
    plots()
    verdict=('A revisão confirmou a integridade das entradas e a implementação das fórmulas principais nos testes independentes realizados. '
        'Encontrou, porém, fragilidades de seleção, controles e interpretação, além de um defeito herdado na costura e imprecisões de documentação. '
        'Os resultados principais foram preservados; a conclusão permanece exploratória, sem melhora robusta nem identificação física demonstradas.')
    short=f'''## Revisão adversarial: o que muda na leitura

**A estabilidade é local, condicionada às âncoras originais.** Ao repetir também a seleção por estrato, apenas **{counts['sigma12']}/76 mínimos em σ=12** e **{counts['sigma20']}/76 em σ=20** voltam a ser selecionados em até 64 px. Na faixa 30–49, são **{counts['faixa30_49']}/76**. As perdas significam escolha de outro marco; não deslocamento do mesmo marco por milhares de pixels.

![Correspondência local versus seleção completa]({F}/repetibilidade_selecao.png)

**O pequeno ganho depende de escolhas operacionais.** Em 100 sementes alternativas, o contraste anular − fixos variou de **{f(ann.minimo)} a {f(ann.maximo)}**, sendo positivo em {int(round(100*ann.fracao_positiva))}/100. No churn, variou de **{f(churn.minimo)} a {f(churn.maximo)}**, positivo em 100/100, sempre abaixo de 0,05. Slug não varia com a semente porque só existe um par equivalente; isso não constitui robustez. Os testes não redefiniram a semente original.

![Sensibilidade aos parceiros escolhidos]({F}/sensibilidade_parceiros.png)

**Um ganho pequeno também apareceu sem evento físico comum.** Em 40 pares de sinais AR(1) independentes, com detector e seleção completos, âncoras tiveram média **{f(null['ancoras_media'])}** e fixos **{f(null['fixos_media'])}**: diferença **{f(null['diferenca_media'])}**. É um diagnóstico posterior, sem calibração para os dados reais; não fornece valor-p nem deve ser subtraído do resultado experimental.

![Controle sintético com a seleção completa]({F}/nulo_selecao.png)

**Dados e fórmulas:** 31 panoramas reproduzidos pixel a pixel, 15.531 quadros e 15.500 deslocamentos conferidos; nenhum quadro inteiro duplicado por hash. A soma direta independente da CWT diferiu em no máximo {sci(maxmath)}. A operação de costura de 2 px usa pesos [0,1] e duplica uma coluna do quadro seguinte, em vez de interpolar: sem ela, churn − fixos passou de 0,0283 para 0,0280 e anular de 0,0044 para 0,0046. Essa mudança pequena não valida fisicamente a montagem. O comentário de direção do deslocamento no v11 está invertido, e a fase complexa do legado refletido não deve ser comparada diretamente à principal; as magnitudes e os escores usados permanecem preservados.

Veja a [auditoria completa, com evidências, tratamento e pendências](AUDITORIA_ADVERSARIAL.md). Os protocolos congelados e os resumos históricos continuam disponíveis; a revisão acrescenta diagnósticos e corrige a interpretação, sem substituir os resultados principais.
'''
    full=f'''# Auditoria adversarial — dados, fórmulas, método e aplicação

{verdict}

## Escopo, preservação e critério de leitura

A revisão abrangeu a etapa 1 (plug/slug), a extensão (churn/anular), a cadeia de montagem v11 e a extração 18–71. Os 264 hashes do manifesto anterior foram conferidos antes das alterações. Relatórios e verificadores anteriores foram copiados para [estado_anterior]({R}/estado_anterior); os dados e códigos congelados permaneceram intactos. O [plano dos diagnósticos]({R}/plano_diagnosticos.json) foi registrado antes das novas WCO desta revisão, com os resultados originais já conhecidos. Portanto, estes testes são exploratórios posteriores, não pré-registro ou validação externa.

O revisor foi o mesmo assistente de IA. “Independente” abaixo designa uma implementação matemática separada, não outro pesquisador. PASS indica integridade/consistência computacional; não certifica identidade de bolhas/ondas, rótulos físicos, estacionariedade ou generalização.

## Achados e tratamento

| ID | Gravidade / natureza | Achado | Tratamento e situação |
| --- | --- | --- | --- |
| A1 | Alta / interpretação | Correspondência próxima foi apresentada como estabilidade sem distinguir a seleção completa. | Corrigidos texto e gráfico; teste completo acrescentado. A repetibilidade continua limitada. |
| A2 | Alta / validade do método | Seleção por mínimos e largura também pode produzir ganho em sinais independentes. | Controle com procedimento completo acrescentado; ganho de WCO não é evidência suficiente de identidade física. |
| A3 | Moderada / amostragem | Poucos parceiros por par de aquisições tornam pequenas diferenças sensíveis à semente. | Todos os 1.095 pares intrapadrão admissíveis avaliados; 100 sementes alternativas por etapa. Original preservado. |
| A4 | Moderada / implementação herdada | Feather de dois pixels usa [0,1], duplicando uma coluna; o comentário sobre o sentido do deslocamento está invertido. | Reprodução e contrafactual quantificados; fonte v11 preservada. Correção física da montagem exige nova análise. |
| A5 | Moderada / relatório | A referência de cobertura de 2/3 migrou de aquisições para pares sem explicitar a mudança do denominador. | Texto corrigido e ambas as coberturas tabuladas; nenhum critério foi redefinido depois dos resultados. |
| A6 | Moderada / uso da fase | Refletir o lag da WCS legada sem conjugar não alinha sua fase à convenção principal. | Restrição documentada: comparação v5 por magnitudes/escores. Não houve uso experimental dessa fase nas conclusões. |

Não foi encontrado erro nas fórmulas principais, na faixa de extração ou nas agregações verificadas. A ausência de erro computacional não remove os achados acima nem garante ausência de falhas não detectadas.

## 1. Dados e cadeia de transformação

Foram relidos **31 conjuntos**, com **15.531 quadros**, seus rótulos locais e **31 planilhas**. A ordem e a sequência dos 501 nomes por aquisição foram verificadas. Todos os **15.500 deslocamentos** foram reestimados por Pearson em todos os lags previstos, sem chamar o pipeline v11; coincidiram com os caches. Os 31 panoramas foram reproduzidos pixel a pixel em memória. A média 18–71, calculada por soma inteira/divisor 54, coincidiu exatamente com os sinais da análise; as planilhas foram conferidas contra sua faixa original 20–59, sem confundi-las com a entrada atual. Não houve quadros inteiros repetidos por SHA-256. Foram registrados hashes de {data['fontes_com_hash']} fontes.

Isso confirma a procedência computacional. Os rótulos foram conferidos contra os arquivos locais, sem uma nova classificação física independente. Os cinco registros excluídos foram auditados, mas suas exclusões não foram revertidas. JPEG, iluminação, paredes incluídas na média e reconstrução de uma dinâmica por faixas permanecem limitações. O eixo horizontal é **coluna do panorama montado**; interpretá-lo como distância percorrida exige hipóteses de advecção/deformação não validadas aqui.

### A4: costura e sinal do deslocamento

Para feather=2, `linspace(0,1,2)` gera [0,1]. A última coluna à esquerda da costura recebe a segunda coluna da faixa seguinte, que também permanece à direita. Não há média gradual entre os quadros nesses dois pontos. Com feather=1, o peso é zero e nada muda. A reprodução independente confirmou exatamente esse comportamento, inclusive os deslocamentos não positivos substituídos por 1 px.

O teste contrafactual retirou **somente o feather**, mantendo deslocamentos, âncoras, recortes, parceiros e a agregação. Foram recalculados 458 pares ancorados intrapadrão e 261 fixos:

{seam_table}

O efeito foi pequeno nesses contrastes. Isso não demonstra que os deslocamentos, a evolução das interfaces ou toda a reconstrução sejam fisicamente corretos. Tampouco inclui redetecção das âncoras no panorama sem feather. Não se substituiu a entrada principal depois de conhecer os resultados.

Um deslocamento sintético real de **+7 px à direita** retorna **−7** na função de Pearson do v11, ao contrário de seu comentário. Isso é uma inconsistência de convenção documental: não prova que o escoamento real se desloca no sentido errado e não justifica inverter os dados. A regra de limitar deslocamentos não positivos a 1 px continua uma hipótese operacional, particularmente relevante em m03q40_001 (60/500).

Evidências: [dados por aquisição]({R}/dados_por_aquisicao.csv), [auditoria dos dados]({R}/dados.json), [fontes]({R}/fontes_sha256.json), [contrastes sem feather]({R}/contrastes_sem_feather.csv), código de origem [pipeline_v11.py](../pipeline_v11.py).

## 2. Fórmulas e implementação

O [capítulo original de Schuck Jr. e Bodmann](../../assets/schuck-bardo-wavelet-coherence.pdf), páginas impressas 411, 413 e 414, foi conferido visualmente. A CWT, a Morlet e a conversão de escala correspondem às Eqs. 32.4, 32.6 e 32.8. As WPS usam módulo e a WCO usa o denominador no lag zero, como 32.14–32.17. O divisor 1/N é uma **adaptação enviesada explícita**; o artigo usa 1/(N−|m|) em 32.21–32.23. Não se afirma reprodução literal desse estimador. A remoção de média/tendência também é estimada no próprio recorte.

Foram selecionados o primeiro e o último ID intrapadrão em cada um dos quatro padrões. Nos oito pares, o preparo foi refeito por regressão independente e a correlação por soma direta. A CWT foi conferida por produto escalar explícito com a Morlet temporal em cinco escalas e sete lags, incluindo extremos. A WCO foi conferida nos 65 lags do indicador, nessas cinco escalas. Também foram verificados troca x/y, preenchimento de zeros e escore completo salvo. Maior diferença da CWT: **{sci(maxmath)}**; nenhuma discrepância material. Essa amostragem independente não é recálculo por soma direta de todas as células dos milhares de mapas.

O limite unitário da Morlet principal tem fundamento adicional: com correlação 1/N, o multiplicador espectral H é real e não negativo. A WCS é um produto interno ponderado de X e Y deslocado; Cauchy–Schwarz limita seu módulo ao quadrado pelo produto das energias ponderadas no lag zero. Essa justificativa não se transfere automaticamente ao estimador 1/(N−|m|) nem ao kernel discreto legado. As três ultrapassagens já documentadas do v5 (máximo 1,006157, todas fora de ±32) permanecem preservadas.

### A6: lag alinhado não é fase alinhada

Para sinais reais, trocar Rxy por Ryx reflete o lag e conjuga a CWT. O comparador v5 reflete o mapa depois da CWT para orientar o lag, preservando as magnitudes utilizadas. Essa reflexão sem conjugação **não alinha a fase complexa** à principal. No exemplo testado, o erro relativo frente à CWT legada da Rxy canônica foi {f(math['fase_legado']['erro_relativo_sem_conjugar'])}; ao conjugar caiu para {f(math['fase_legado']['erro_relativo_conjugando'])}, com diferença residual da discretização. Não se deve interpretar essa fase refletida como atraso/fase na convenção principal. A fórmula da Morlet principal não apresenta esse problema.

A ambiguidade do artigo entre a média complexa (32.16) e a média de módulos impressa (32.20) já estava identificada na base e foi confirmada. Com K=1, não muda as magnitudes/WCO; com K>1, pode mudar radicalmente. Duas realizações de fases opostas se cancelam na média complexa e não na média de módulos. As médias descritivas entre pares desta entrega **não são WCO de ensemble**.

Evidências: [conferência das fórmulas]({R}/formulas.json), [erros por par]({R}/formulas_comparacoes.csv), [auditoria anterior das equações](../v16_spatial2_18_71/EQUATION_AUDIT.md).

## 3. Âncoras, seleção e parceiros

### A1: duas noções diferentes de estabilidade

O teste original procura o candidato do mesmo tipo mais próximo do centro congelado, até 64 px. Ele não repete a escolha de maior proeminência por estrato. Refazer o procedimento completo produz:

{select_table}

![Repetibilidade da seleção]({F}/repetibilidade_selecao.png)

As coordenadas dos candidatos correspondentes podem variar até 11 px e, ainda assim, outro candidato vencer a seleção. Quando isso ocorre, não é correto chamar a distância ao novo selecionado de deslocamento físico da mesma estrutura. O novo teste confirma repetibilidade parcial; não invalida a existência das âncoras originais nem altera seu catálogo. Não foram recalculadas todas as WCO de um catálogo totalmente reselecionado: essa é uma análise futura, com novos parceiros e nova interpretação de cobertura.

### A3: sensibilidade ao limite de parceiros

O hash determinístico evita escolher pares pela maior WCO, mas não elimina a variabilidade causada por usar no máximo 12 pares na etapa 1 e 2 na extensão. Foram avaliados **1.095 pares intrapadrão admissíveis** de compatíveis/equivalentes e distintos, sendo 799 novos; os demais escores foram reutilizados. Mantiveram-se as âncoras. Foram aplicadas a semente original e 100 sementes consecutivas adicionais, preservando pesos por par de aquisições e as interseções de comparação:

{seed_table}

![Sensibilidade aos parceiros]({F}/sensibilidade_parceiros.png)

Usar todos os parceiros admissíveis resultou em ganhos frente a fixos de **0,0239 em plug**, **0,2192 em slug**, **0,0298 em churn** e **0,0086 em anular**. Slug continua sustentado por um único par equivalente. A faixa entre sementes mede decisões de amostragem dentro do mesmo catálogo, não incerteza populacional. Não se escolheu uma nova semente ou limite que favorecesse a conclusão.

### A5: denominadores de cobertura

No contraste com mínimos distintos, churn tem **10/10 aquisições**, mas **29/45 pares possíveis**; anular tem **8/9 aquisições**, mas **13/36 pares**. A etapa 1 definiu sua referência de 2/3 por aquisição. O texto da extensão não deve transformá-la em um critério por pares sem declarar a mudança. O resumo foi corrigido; os dois denominadores aparecem em [cobertura_denominadores.csv]({R}/cobertura_denominadores.csv). Ter todas as aquisições representadas não elimina a escassez de parceiros comparáveis.

## 4. A2: controle de todo o procedimento

Os controles sintéticos anteriores testavam formas e o núcleo, mas não isolavam o efeito de escolher extremos em sinais independentes. Agora foram gerados **40 pares independentes de AR(1), φ=0,95**, com 8.192 pontos após aquecimento de 1.024, intensidade média 128 e desvio 20. Foram aplicados o detector, a seleção por estrato, a compatibilidade de largura, o limite de dois parceiros e o contraste com três janelas fixas, antes de consultar a WCO. Todos os 40 casos foram calculáveis; não houve seleção dos casos de maior efeito.

A média ancorada foi **{f(null['ancoras_media'])}**, contra **{f(null['fixos_media'])}** dos fixos; diferença **{f(null['diferenca_media'])}**, positiva em **{int(round(null['calculaveis']*null['fracao_positiva']))}/40** repetições. Assim, mínimos escolhidos em séries independentes de estrutura estatística semelhante também podem gerar ganho. O resultado não prova que todo ganho experimental seja artefato, nem constitui um nulo calibrado para esses dados. Não foram calculados valores-p, nem subtraído esse ganho sintético dos contrastes reais.

![Controle sintético completo]({F}/nulo_selecao.png)

## 5. Aplicação, dependência e limites que permanecem

As agregações principais foram conferidas por par de aquisições e interseção dos controles. Não houve agrupamento silencioso de plug/slug com churn/anular, nem conversão das 15 ausências de σ20 em zero. As variantes do detector e o painel comum de 815 pares preservam os parceiros originais; não representam a variabilidade de um novo catálogo completo.

Pares compartilham aquisições, recortes se sobrepõem e a CWT usa correlações de todo o contexto. A WCO não localiza automaticamente o evento central. A compatibilidade de largura não controla completamente iluminação, contraste, periodicidade, outras estruturas ou costuras. Mínimo × máximo não é um negativo obrigatório, pois a WCO é invariável à inversão global. As janelas fixas tampouco são necessariamente livres de eventos.

Os 200 comprimentos de onda e 65 lags não são réplicas independentes. O piso de denominador é numérico; o suporte geométrico e o critério de dois comprimentos de onda por segmento não certificam potência, estacionariedade ou significância. O limiar de ganho 0,05 é operacional e descritivo. Não há validação externa, revisão física independente, calibração espaço/tempo ou desempenho de classificador. Esses limites já existiam; os novos testes reforçam sua relevância.

## 6. Alterações e próximas decisões científicas

Foram atualizados resumo executivo, PDF integrado, passo a passo, resultados, nota de âncoras/autocorrelação, README, registro de decisões e estado de retomada. A figura de estabilidade original passou a explicitar que trata de correspondência local. As tabelas principais, mapas, protocolos congelados e bases anteriores foram preservados. Os relatórios históricos permanecem como registros da etapa correspondente, com o estado atual indicado nos documentos principais.

Uma análise posterior deverá definir claramente o estimando: repetibilidade de marco óptico, semelhança do contexto, identidade física ou classificação são objetivos diferentes. Para afirmar identidade física, será necessária anotação independente e rastreamento/validação da montagem. Para avaliar estabilidade operacional completa, será necessário refazer seleção e parceiros em cada configuração. Uma nova montagem deve corrigir e testar a costura e a convenção de deslocamento antes de produzir outra versão dos dados. Mais combinações dos mesmos registros não substituem novas aquisições.

Os diagnósticos WCO desta revisão consumiram **{f(sum(v['segundos'] for v in diag['etapas'].values())/60,1)} minutos**; leitura de dados, conferência matemática e redação são etapas separadas. As tabelas contrafactuais exportam somente os escores recalculados, para não misturar outras métricas herdadas; a [revisão de exportação]({R}/revisao_exportacao.json) registra que os escores permaneceram exatamente iguais.

Confira a [verificação adversarial]({R}/verificacao_final.json), a [verificação integrada](results/verificacao_final.json) e o [manifesto atualizado](results/manifesto.json). Os [scripts de auditoria dos dados](scripts/adversarial_data_audit.py), [âncoras](scripts/adversarial_anchor_audit.py), [fórmulas](scripts/adversarial_math_audit.py) e [diagnósticos](scripts/adversarial_diagnostics.py) registram as operações executadas.
'''
    if (ROOT/'results/extensao_churn_anular/verificacao_exemplos.json').exists():
        full+='''
## 7. Conferência posterior dos exemplos de maior e menor escore

Após questionamento sobre as regiões escuras, foram conferidos os seis exemplos contra os IDs e mapas salvos, refeitos seus preparos por regressão e suas correlações por soma direta, e reaplicada a CWT principal. Os quatro extremos estão corretos **na amostra original** de 90 pares compatíveis de churn e 66 de anular: churn 0,0818–0,5007; anular 0,0945–0,4219. Ambos os extremos usam mínimos ópticos por definição do grupo; o escore mede a WCO dos contextos, não a intensidade da âncora.

A legenda anterior não delimitava suficientemente esse universo. Entre todos os admissíveis já calculados na auditoria, são 484 pares de churn (0,0793–0,5007) e 348 de anular (0,0945–0,4476). As imagens preservaram os exemplos originais, com títulos explícitos sobre o número de pares e a diferença entre mínimo óptico e WCO média. Resumo, resultados, nota metodológica e PDF foram atualizados. A conferência passou a integrar o verificador, sem alterar tabelas ou mapas experimentais. [Evidências dos exemplos](results/extensao_churn_anular/verificacao_exemplos.json); [script](scripts/audit_example_scores.py).
'''
    (ROOT/'AUDITORIA_ADVERSARIAL.md').write_text(full)
    append('EXECUTIVE_SUMMARY.md',short,early=True)
    append('RESULTS.md','## Revisão adversarial dos resultados\n\n'+verdict+'\n\n'+select_table+'\n\n'+seed_table+'\n\n'+seam_table+f'\n\nControle com seleção completa: {f(null["diferenca_media"])} de ganho médio em AR(1) independentes, sem valor-p. Consulte [a auditoria completa](AUDITORIA_ADVERSARIAL.md), [sementes]({R}/sensibilidade_semente.csv), [pares admissíveis]({R}/todos_pares_admissiveis.csv) e [o controle sintético]({R}/nulo_selecao_completa.csv).')
    append('PASSO_A_PASSO.md',f'''## Continuação: revisão adversarial solicitada

13. Conferidos os 264 hashes anteriores e preservados os documentos antes das mudanças. Reestimados 15.500 deslocamentos e reproduzidos 31 panoramas a partir de 15.531 quadros. **Limite:** concordância com v11 não valida a reconstrução física.
14. Conferidas as equações no PDF e por somas temporais independentes em oito pares. Explicitadas a adaptação 1/N, a ambiguidade de K>1 e a limitação da fase legada refletida. **Limite:** teste numérico amostral, sem certificação estatística.
15. Repetidas detecção e seleção por estrato: {counts['sigma12']}/76 e {counts['sigma20']}/76 mínimos recuperados em σ12/20. **Correção:** a medida original é correspondência local, não repetibilidade de todo o procedimento.
16. Avaliados todos os 1.095 pares admissíveis intrapadrão e 100 sementes alternativas por etapa. **Limite:** novas combinações dos mesmos registros não são novas aquisições nem intervalos de confiança.
17. Quantificado o efeito de retirar somente o feather e executados 40 controles AR(1) com seleção completa. **Limite:** diagnósticos posteriores, sem nulo calibrado para os dados e sem substituir as entradas originais.
18. Atualizados documentos, gráficos e PDF; explicitados os denominadores de cobertura; conferida a preservação dos resultados principais. A [auditoria adversarial](AUDITORIA_ADVERSARIAL.md) distingue falhas de relato corrigidas de limitações científicas que permanecem.
19. Após questionamento sobre os pontos escuros, conferidos os seis exemplos e os extremos dentro da amostra original. Preparo e correlação foram refeitos; títulos agora distinguem WCO média de intensidade óptica e delimitam 90 pares selecionados de churn e 66 de anular. **Limite:** os extremos de todos os admissíveis da auditoria são diferentes em dois casos. [Conferência dos exemplos](results/extensao_churn_anular/verificacao_exemplos.json).
''')
    append('ANCORAS_E_AUTOCORRELACAO.md',f'''## Ressalvas da revisão adversarial

A correspondência local em σ12/20 não equivale à seleção completa: esta recuperou {counts['sigma12']}/76 e {counts['sigma20']}/76 mínimos originais em até 64 px. Um marco pode continuar localizável e perder a posição de maior proeminência no estrato.

A correlação 1/N e a Morlet principal com multiplicador espectral real não negativo permitem justificar WCO≤1 por Cauchy–Schwarz. Esse limite não valida fisicamente a janela nem se transfere automaticamente ao estimador não enviesado ou ao kernel legado. As Eqs. 32.16 e 32.20 do artigo divergem para K>1; aqui K=1 e as médias entre pares são descritivas.

**Fase legada:** refletir somente o eixo de lag da WCS v5 conserva as magnitudes, mas não a fase na convenção principal. A comparação legada desta entrega é de magnitudes e escores; não se devem extrair atrasos físicos dessa fase refletida.

**Eixo do sinal:** c é uma coluna da montagem, sem calibração de deslocamento físico. O v11 duplica uma coluna na operação de feather de 2 px e seu comentário sobre o sentido do deslocamento está invertido. Os panoramas foram reproduzidos exatamente e o efeito de retirar o feather foi pequeno nos contrastes testados, mas a validade física da reconstrução permanece em aberto.

O procedimento completo também produziu ganho médio de {f(null['diferenca_media'])} em sinais AR(1) independentes. Portanto, localizar extremos e obter WCO relativamente maior não demonstra que se reconheceu o mesmo evento. Veja [evidências e limites da revisão](AUDITORIA_ADVERSARIAL.md).
''')
    append('README.md',f'''## Revisão adversarial incorporada

Leia [AUDITORIA_ADVERSARIAL.md](AUDITORIA_ADVERSARIAL.md) antes de interpretar “estabilidade” ou um pequeno ganho de WCO. A revisão corrigiu a distinção entre correspondência local e seleção completa, os denominadores de cobertura e a descrição da fase legada. Acrescentou auditoria integral dos dados, sensibilidade a parceiros e costuras e controle sintético com o procedimento completo. Os dados principais permanecem congelados.

Resultados e evidências estão em [{R}]({R}); os três gráficos novos, em [{F}]({F}). O gerador write_texture_reports.py integra automaticamente a revisão quando diagnosticos.json está completo. Não execute diretamente a função integrate sobre um PDF já integrado: o ponto de entrada de reprodução é write_texture_reports.py, que recompõe o PDF antes de acrescentar a revisão.

Os scripts adversarial_data_audit.py, adversarial_anchor_audit.py, adversarial_math_audit.py e adversarial_diagnostics.py permitem refazer os diagnósticos em uma cópia separada. O último aceita --stage admissible, seams ou synthetic. O código exato executado e a restrição posterior das colunas exportadas estão registrados em [revisao_exportacao.json]({R}/revisao_exportacao.json). A auditoria integrada chama também verify_adversarial.py, incluindo a preservação das tabelas/matrizes anteriores e os vínculos dos novos documentos.
''')
    summary=dict(data_utc=datetime.now(timezone.utc).isoformat(),conclusao='ANALISE_EXPLORATORIA_COM_RESSALVAS',
        achados=['A1: correspondência local versus seleção completa','A2: ganho em nulo com seleção','A3: sensibilidade a parceiros','A4: feather e direção herdados','A5: denominadores de cobertura','A6: fase legada'],
        selecao_completa_minimos=counts,erro_CWT_independente_maximo=maxmath,
        nulo_ganho=null['diferenca_media'],sementes_anular_positivas=int(round(100*ann.fracao_positiva)),
        arquivos_principais_atualizados=['AUDITORIA_ADVERSARIAL.md','EXECUTIVE_SUMMARY.md','RESUMO_EXECUTIVO.pdf','PASSO_A_PASSO.md','RESULTS.md','ANCORAS_E_AUTOCORRELACAO.md','README.md','RETOMADA.md','DECISIONS.md'])
    (OUT/'conclusao.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
    # Cinco páginas de revisão precedem as vinte páginas da entrega integrada.
    supplement=OUT/'resumo_auditoria.pdf'
    with PdfPages(supplement,metadata={'Title':'Revisão adversarial — âncoras e WCO','Author':'Projeto TCC'}) as pdf:
        page(pdf,'Revisão adversarial\nÂncoras e wavelets',[
            ('Conclusão',verdict),
            ('Correção da estabilidade',f'A correspondência local preserva 76/76 mínimos em σ12/20. Refazer também a seleção por estrato recupera {counts["sigma12"]}/76 e {counts["sigma20"]}/76. Localizar um candidato próximo não garante escolhê-lo novamente.'),
            ('Ganhos pequenos exigem cautela',f'Anular menos fixos mudou de sinal em 20/100 sementes alternativas. O ganho médio em 40 pares AR(1) independentes foi {f(null["diferenca_media"])} após detecção e seleção completas. Não são testes de significância dos dados experimentais.'),
            ('O que foi conferido','31 panoramas reproduzidos exatamente; 15.531 quadros, 15.500 deslocamentos e 31 planilhas conferidos. CWT verificada por somas independentes em oito pares de quatro padrões. PASS computacional não certifica validade física.'),
            ('Preservação','Os resultados principais e protocolos congelados foram preservados. As páginas seguintes distinguem os novos diagnósticos da análise original. Evidências e pendências: AUDITORIA_ADVERSARIAL.md.')])
        for filename,title in [('repetibilidade_selecao','Correspondência local e seleção completa'),('sensibilidade_parceiros','Escolha dos parceiros e fragilidade do ganho'),('nulo_selecao','Controle sintético do procedimento completo')]:
            fig=plt.figure(figsize=(11.69,8.27));ax=fig.add_axes([.025,.09,.95,.80]);ax.imshow(plt.imread(FIG/(filename+'.png')));ax.axis('off')
            fig.suptitle(title,fontsize=16,weight='bold',y=.97);fig.text(.04,.035,'Diagnóstico posterior aos resultados; não substitui o experimento congelado nem fornece validação física.',fontsize=8)
            pdf.savefig(fig);plt.close(fig)
        page(pdf,'Falhas documentadas\ne limites que permanecem',[
            ('Montagem herdada','O feather de 2 px usa pesos [0,1] e duplica uma coluna; o comentário sobre o sentido do deslocamento está invertido. Sem feather, churn menos fixos passa de 0,0283 para 0,0280 e anular de 0,0044 para 0,0046. O efeito pequeno não valida fisicamente o panorama.'),
            ('Cobertura e fase','Aquisições e pares são denominadores distintos. No contraste com larguras distintas: churn 10/10 aquisições e 29/45 pares; anular 8/9 e 13/36. A fase complexa da v5 refletida não está na convenção principal; sua comparação deve ficar nas magnitudes e escores utilizados.'),
            ('Fórmulas e inferência','O estimador 1/N é adaptação explícita do artigo. K=1; médias entre pares não são ensemble. Pares compartilham aquisições e contextos. Não há teste externo, anotação física independente, valor-p calibrado, calibração em mm/tempo ou classificador validado.'),
            ('Estado da entrega','Resumo, passo a passo, resultados e nota do método foram atualizados. Os documentos históricos permanecem identificados. Novas montagens ou critérios devem gerar outra análise e preservar o congelamento atual.')])
    target=OUT/'resumo_integrado_com_auditoria.pdf'
    subprocess.run(['pdfunite',str(supplement),str(ROOT/'RESUMO_EXECUTIVO.pdf'),str(target)],check=True)
    (ROOT/'RESUMO_EXECUTIVO.pdf').write_bytes(target.read_bytes())


if __name__=='__main__':
    from write_texture_reports import main
    main()
