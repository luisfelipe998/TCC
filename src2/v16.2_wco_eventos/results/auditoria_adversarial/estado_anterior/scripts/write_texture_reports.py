"""Integra resultados da extensão aos relatórios, preservando a etapa 1."""
from pathlib import Path
import json
import os
import shutil
import subprocess
import textwrap
os.environ.setdefault('MPLCONFIGDIR','/tmp/wco-eventos-mpl')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from analysis_common import ROOT
from texture_common import OUT,CFG,LABELS,GROUPS,VARIANT_LABELS

F='figures/extensao_churn_anular'
R='results/extensao_churn_anular'
COMP={'fixos':'Janelas fixas','minimos_distintos':'Mínimos de larguras distintas','polaridades_opostas':'Mínimo × máximo'}

def f(x,d=4):
    return 'não estimável' if pd.isna(x) else f'{x:.{d}f}'.replace('.',',')

def table(head,rows):
    return '\n'.join(['| '+' | '.join(head)+' |','| '+' | '.join(['---']*len(head))+' |']+['| '+' | '.join(map(str,row))+' |' for row in rows])

def write(name,text):(ROOT/name).write_text(text.rstrip()+'\n')

def cost_note(timing):
    if timing.get('pico_memoria_mib') is None:
        return (f'As etapas concluídas somaram {f(timing["segundos_total"]/60,1)} minutos de cálculo. '
                'Esse total exclui a pausa e o tempo parcial dos diagnósticos interrompidos, que não foi cronometrado separadamente. '
                'O pico de memória do primeiro processo não foi salvo antes da pausa; '
                f'o processo de retomada dos diagnósticos teve pico de {f(timing["pico_memoria_retomada_mib"],0)} MiB.')
    return f'O lote levou {f(timing["segundos_total"]/60,1)} minutos, com pico de {f(timing["pico_memoria_mib"],0)} MiB.'

def archive():
    for stem in ['EXECUTIVE_SUMMARY','PASSO_A_PASSO','RESULTS','README','ANCORAS_E_AUTOCORRELACAO']:
        target=ROOT/(stem+'_ETAPA1.md')
        if not target.exists():shutil.copyfile(ROOT/(stem+'.md'),target)
    p=ROOT/'RESUMO_EXECUTIVO_ETAPA1.pdf'
    if not p.exists():shutil.copyfile(ROOT/'RESUMO_EXECUTIVO.pdf',p)

def main():
    assert json.loads((OUT/'execucao.json').read_text())['status']=='PASS'
    archive()
    c=json.loads((OUT/'resumo.json').read_text());timing=json.loads((OUT/'execucao.json').read_text())
    cs=pd.read_csv(OUT/'contrastes_por_padrao.csv');cp=pd.read_csv(OUT/'contrastes_pareados.csv')
    stability=pd.read_csv(OUT/'estabilidade_ancoras.csv');coverage=pd.read_csv(OUT/'cobertura_grupos.csv')
    common=pd.read_csv(OUT/'sensibilidade_conjunto_comum.csv');leave=pd.read_csv(OUT/'retirada_aquisicao.csv')
    cat=pd.read_csv(OUT/'catalogo_ancoras.csv').query('selecionada');examples=pd.read_csv(OUT/'exemplos.csv')
    synthetic=pd.read_csv(OUT/'controles_sinteticos.csv').set_index('caso')
    get=lambda p,co='fixos',v='principal':cs[(cs.classe==p)&(cs.comparador==co)&(cs.variante==v)].iloc[0]
    ch,an=get('churn'),get('wavy annular')
    lead=(f'Foi possível determinar âncoras ópticas e aplicar WPS, WCS e WCO alternativas em **todas as 19 aquisições elegíveis** de churn e anular ondulado. '
          f'Foram selecionados **76 mínimos principais e 76 máximos de controle**, antes de calcular a nova WCO. '
          f'Os 76 mínimos mantiveram correspondência nas duas suavizações alternativas, com deslocamento máximo de 11 px. '
          f'**O objetivo operacional foi atendido; não se demonstrou ganho robusto de coerência nem separação entre padrões.** '
          f'Na análise principal, mínimos compatíveis superaram janelas fixas em {f(ch.diferenca)} no churn e {f(an.diferenca)} no anular, abaixo da referência descritiva de 0,05.')
    patterns=table(['Padrão global da aquisição','Aquisições elegíveis','Âncora usada','Situação nesta pasta'],[
        ['Bolhas dispersas',1,'Sem comparação ancorada entre aquisições','Só atlas/referência fixa; sem repetição'],
        ['Pistonado / bolhas alongadas (plug)',3,'Centro do intervalo óptico alongado','Etapa 1: 20 eventos aceitos'],
        ['Golfadas (slug)',3,'Centro do intervalo óptico alongado','Etapa 1: 6 eventos aceitos'],
        ['Agitado (churn)',10,'Mínimo local de intensidade suavizada','Extensão: 40 mínimos + 40 máximos'],
        ['Anular ondulado (wavy annular)',9,'Mínimo local de intensidade suavizada','Extensão: 36 mínimos + 36 máximos']])
    main_table=table(['Padrão','Controle','WCO das âncoras','WCO do controle','Diferença','Pares de aquisições','Fração positiva'],[
        [LABELS[p],COMP[co],f(get(p,co).ancoras),f(get(p,co).controle),f(get(p,co).diferenca),int(get(p,co).pares_aquisicoes),f(100*get(p,co).fracao_positiva,1)+'%']
        for p in CFG['padroes'] for co in COMP])
    stabtable=table(['Padrão / âncoras principais','σ=12: correspondências / mediana / máximo','σ=20: correspondências / mediana / máximo','Faixa 30–49: mediana / máximo'],[
        [LABELS[p]]+[f'{int((s:=stability[(stability.padrao==p)&(stability.tipo=="minimo")&(stability.alternativa==alt)].iloc[0]).disponiveis)}/{int(s.total)}; {f(s.mediana_abs_px,1)} / {f(s.max_abs_px,0)} px' for alt in ['sigma12','sigma20']]+
        [f'{f((s:=stability[(stability.padrao==p)&(stability.tipo=="minimo")&(stability.alternativa=="faixa30_49")].iloc[0]).mediana_abs_px,1)} / {f(s.max_abs_px,0)} px'] for p in CFG['padroes']])
    limits=f'''- **Âncora óptica, identidade física não validada.** O mínimo marca uma depressão do perfil de intensidade; não mede espessura do filme nem identifica uma crista, bolha ou a mesma onda em aquisições diferentes. Os rótulos churn/anular são globais, herdados do inventário.
- **Seleção exploratória.** Regra desenvolvida em quatro aquisições piloto e aplicada no mesmo conjunto, já conhecido; revisão visual assistida por IA, sem especialista independente, concordância entre anotadores ou conjunto externo.
- **Montagem e contexto.** {c['ancoras_em_mistura']} das 152 âncoras caem nas zonas de mistura de 2 px das costuras; as outras 142 tiveram a coluna de origem conferida pixel a pixel. Cada contexto cruza {c['costuras_contexto_min']}–{c['costuras_contexto_max']} costuras. {c['ancoras_com_sobreposicao']} recortes sobrepõem outro recorte selecionado e {c['ancoras_com_outros_marcos']} contêm outro marco selecionado. Considerando só os 76 mínimos principais: {c['minimos_com_sobreposicao']} sobrepõem outro segmento, {c['minimos_com_outros_marcos']} contêm outro marco e {c['minimos_em_mistura']} estão em mistura. Isso impede atribuir toda a WCO ao ponto central. As contagens de proximidade de costuras não são prova de artefato, pois as faixas de montagem são densas.
- **Aquisição com ressalva.** m03q40_001 conserva 60 deslocamentos não positivos herdados; {c['contextos_com_deslocamentos_nao_positivos']} dos 152 contextos selecionados incluem deslocamentos não positivos. A análise principal mantém essa aquisição, e sua retirada é mostrada separadamente; não houve reconstrução dos panoramas.
- **Comparabilidade dos controles.** Mínimos compatíveis têm razão de larguras ≤1,50; distintos, ≥2,00. As larguras medem meia proeminência óptica. Mínimo × máximo pode diferir também em largura e forma, e a WCO é invariável à inversão global de sinal. Janelas fixas podem conter estruturas e não são negativos físicos garantidos.
- **Cobertura e dependência.** Os controles de larguras distintas existem em menos pares de aquisições, especialmente no anular. Cada contraste usa somente sua interseção de aquisições. Pares compartilham registros; não são réplicas independentes. Não foram calculados valores-p nem intervalos de confiança artificiais.
- **Domínio e potência.** N=2.048, K=1, correlação enviesada 1/N, λ de 50 a 1.024 px e lags principais ±32 px. A WPS alternativa usa módulo; WCO não recebe clipping. Validade numérica do denominador não equivale a potência confiável. Não há conversão para mm, Hz ou velocidade, nem validação de classificador.
- **Núcleo legado v5.** Dois mapas ancorados e um fixo apresentaram pequenas ultrapassagens de WCO=1, chegando a {f(c['maximo_wco_v5'],6)}. Todas ocorreram fora de ±32 px. A soma direta das correlações e a CWT por convolução temporal reproduziram os valores, com discrepância máxima de 3,04×10⁻¹³. Foram preservados sem clipping como limitação da variante discreta; a Morlet principal permaneceu ≤1. [Auditoria específica](results/extensao_churn_anular/auditoria_ultrapassagens_v5.json).
'''
    sensitivity_table=table(['Variante','Churn: âncoras − fixos','Anular: âncoras − fixos'],[
        [VARIANT_LABELS[v],f(get('churn',v=v).diferenca),f(get('wavy annular',v=v).diferenca)] for v in CFG['variantes']])
    influence=leave[(leave.retirada=='m03q40_001')&(leave.classe=='wavy annular')&(leave.comparador=='fixos')].iloc[0]
    extext=''
    for p,slug in [('churn','churn'),('wavy annular','anular')]:
        extext+=f'\n### {LABELS[p]} × {LABELS[p]}\n\n'
        for kind,title in [('predeterminado','Exemplo predefinido pelo primeiro ID'),('maior_escore','Maior escore; ilustração posterior'),('menor_escore','Menor escore; ilustração posterior')]:
            r=examples[(examples.classe==p)&(examples.tipo==kind)].iloc[0]
            extext+=f'**{title}: {f(r.escore)}.** x e y pertencem ambos a **{LABELS[p]}**.\n\nPar: `{r.pair_id}`.\n\n![{title}; ambos os sinais de {LABELS[p]}]({F}/par_{slug}_{kind}.png)\n\n'
    executive=f'''# Resumo executivo — âncoras e wavelets alternativas

{lead}

Esta entrega estendeu **v16.2_wco_eventos** sobre a base **v16_spatial2_18_71**, mantendo os cálculos originais de plug/slug. O [resumo da etapa 1](EXECUTIVE_SUMMARY_ETAPA1.md) e seu [PDF original](RESUMO_EXECUTIVO_ETAPA1.pdf) foram preservados. O [PDF integrado](RESUMO_EXECUTIVO.pdf) apresenta a extensão primeiro e a etapa 1 como apêndice histórico.

## Qual padrão é qual

{patterns}

**Churn é o padrão agitado; wavy annular é o anular ondulado.** As cores roxa e verde identificam esses padrões nos novos gráficos. Nas figuras de pares, **x usa linha contínua; y, tracejada**. Cada painel informa o padrão, a aquisição, a âncora e os limites do sinal. As cores dos mapas representam seus valores numéricos, conforme a barra lateral.

## Como conseguimos determinar as âncoras

Da média das linhas **18–71 inclusive**, foram encontrados mínimos e máximos após suavização gaussiana **σ=16 px**, usada somente na detecção. Em cada um de quatro estratos espaciais, selecionou-se o extremo de maior proeminência, respeitando suporte e separação. Proeminência mínima: max(3 níveis de cinza; 25% de P95−P5); largura à meia proeminência: 16–512 px. Dos **1.365 candidatos**, foram mantidas **152 âncoras**, sem consultar WPS/WCS/WCO.

Cada sinal usa **[c−1024, c+1024)**: 2.048 amostras da média original, sem o filtro gaussiano de detecção. A âncora está no índice local 1.024; ela se torna posição relativa zero. Os mínimos de intensidade tendem a aparecer como elevações após a inversão 255−I, sujeitas à remoção de tendência. A CWT alternativa recebe Rxx, Ryy e Rxy desses recortes, e não a sequência de coordenadas das âncoras.

![Âncoras, sinais e WPS alternativa de churn e anular ondulado]({F}/ancoras_dos_padroes.png)

*Exemplos predefinidos: primeiro mínimo selecionado da primeira aquisição de cada padrão. A sombra no detalhe óptico é a largura à meia proeminência, não um contorno físico. O eixo de posição do sinal difere do lag dos mapas.*

## Cobertura e estabilidade

{stabtable}

Todos os **76 mínimos principais** permaneceram localizáveis em σ=12 e σ=20. Nos máximos, houve três ausências em σ=20. Ao trocar apenas a faixa de detecção para 30–49, os 76 mínimos continuaram localizáveis, mas o churn mostrou mudanças de até 47 px; esse diagnóstico não reposicionou a variante wavelet 30–49, que conservou os centros principais.

![Cobertura e deslocamento das âncoras]({F}/cobertura_estabilidade.png)

## O que mudou na WCO

Foram congelados **830 pares de âncoras** entre aquisições: mínimos de larguras compatíveis, mínimos de larguras distintas e mínimo × máximo. A escolha de até dois pares por combinação de aquisições/grupo usou hash determinístico, sem procurar a maior WCO. Primeiro foi calculada a média por par de aquisições; depois, a média do padrão. **A coluna WCO das âncoras muda entre controles porque a interseção de aquisições pode mudar.**

{main_table}

**Churn:** diferença pequena frente a fixos ({f(ch.diferenca)}) e frente a mínimos distintos ({f(get('churn','minimos_distintos').diferenca)}). **Anular ondulado:** diferença quase nula frente a fixos ({f(an.diferenca)}) e negativa frente a mínimos distintos ({f(get('wavy annular','minimos_distintos').diferenca)}). As médias principais não atingiram o efeito descritivo de 0,05. Esse limiar foi definido no protocolo; não é um teste estatístico nem um limiar físico.

A cobertura do contraste frente a fixos foi **45/45 pares possíveis no churn** e **34/36 no anular**. Frente a larguras distintas, caiu para **29/45** e **13/36**, abaixo da referência de cobertura de 2/3. A presença de âncoras nas 19 aquisições não elimina essa limitação de parceiros comparáveis.

![Contrastes pareados por padrão e controle]({F}/contrastes.png)

Entre **churn × anular**, o ganho sobre fixos foi {f(get('entre').diferenca)} em {int(get('entre').pares_aquisicoes)} pares de aquisições. Isso reforça que obter coerência em recortes ancorados não demonstra distinção entre classes; nenhuma acurácia de classificação foi estimada.

![Comparações dentro de cada padrão e entre churn e anular]({F}/entre_padroes.png)

## Sensibilidades e influência das aquisições

{sensitivity_table}

No churn, trocar a média 18–71 pela faixa 30–49 reduziu o contraste frente a fixos de **{f(ch.diferenca)} para {f(get('churn',v='faixa30_49').diferenca)}**, mantendo os centros principais. A posição estável sob mudanças de σ não garante estabilidade da análise quando muda a representação vertical da imagem.

Os parceiros permaneceram congelados. Foram planejadas **6.640 avaliações ancoradas**: **{c['comparacoes_calculadas']} calculáveis** e **{c['comparacoes_nao_estimaveis']} não estimáveis**, todas as ausências por falta de correspondência óptica na variante σ=20. Elas não viraram zeros nem centros originais. A [tabela no conjunto comum de pares](results/extensao_churn_anular/sensibilidade_conjunto_comum.csv) distingue alteração de valor de alteração de cobertura.

![Oito variantes, dois padrões e três controles]({F}/sensibilidades.png)

*O painel usa a interseção dos mesmos **815 pares locais** nas oito variantes, incluindo na linha principal. Apenas 15 pares de controle de polaridade foram retirados desse painel; por isso seu contraste principal de mínimo × máximo pode diferir da tabela com todos os pares. As comparações com fixos e larguras distintas mantêm os mesmos casos. [Dados do painel comum](results/extensao_churn_anular/sensibilidade_conjunto_global.csv).*

Ao retirar m03q40_001 do anular, a diferença frente a fixos passou de {f(an.diferenca)} para **{f(influence.diferenca)}**, em {int(influence.pares_restantes)} pares restantes. Cada retirada eliminou todos os pares incidentes na aquisição; as variações não são intervalos de confiança.

![Influência da retirada de aquisições]({F}/retirada_aquisicao.png)

## Exemplos dos dois sinais e das três medidas wavelet

Os exemplos predefinidos usam o primeiro ID de par compatível. Os extremos abaixo foram escolhidos **somente depois dos resultados, para ilustração**, sem mudar o catálogo, o pareamento ou as conclusões médias. Cada figura mostra os dois recortes, os sinais preparados, as correlações, WPS de x/y, módulo da WCS e WCO. A média principal usa ±32 px de lag; o mapa integral conserva os lags ±2.047 px. Linhas brancas indicam suporte geométrico, sem significado de significância estatística.

Para manter a renderização leve, os mapas nos PNGs exibem um lag a cada 4 px fora da região ±64 px e todos os lags dentro dela. **Os cálculos e os arquivos NPZ preservam as 200 × 4.095 posições completas**, sem essa redução de apresentação.
{extext}
## Controles e limites de interpretação

![Controles sintéticos de forma, fundo comum e polaridade]({F}/controles_sinteticos.png)

Os testes confirmaram WCO invariável ao trocar y por −y; a WCS muda de sinal. Portanto, mínimo × máximo não é um controle que obrigatoriamente deve apresentar WCO baixa. Fundo comum e periodicidade compartilhada também podem gerar coerência, sem equivalência física do marco central.

No controle sintético, mínimos de larguras distintas tiveram WCO média **{f(synthetic.loc['minimos_larguras_distintas','escore'])}**. Ao acrescentar o mesmo fundo periódico aos dois, a média passou a **{f(synthetic.loc['fundo_comum_larguras_distintas','escore'])}**. A contribuição do contexto pode, portanto, dominar a interpretação do marco central.

{limits}

## O que permanece da etapa 1

Plug apresentou diferenças de **0,0275** frente a eventos distintos e **0,0289** frente a fixos, em três pares de aquisições. Slug apresentou **0,1520** e **0,2192**, sustentados por um único par equivalente; ao usar linhas 30–49, o contraste frente a distintos passou a **−0,0427**. Esses resultados continuam exploratórios. **Não foram combinados com churn/anular numa única estimativa**, pois a definição de âncora e de compatibilidade é diferente.

## Execução e rastreabilidade

{cost_note(timing)} A estimativa inicial de aproximadamente 13–30 minutos ficou abaixo do tempo observado. O lote foi pausado por solicitação do usuário após salvar os pares ancorados e todas as referências fixas; a retomada refez somente os diagnósticos interrompidos. Foram reutilizadas 513 comparações fixas principais e calculadas 1.539 sensibilidades fixas, 152 autocomparações e 19 pares da mesma aquisição. As inspeções, figuras, redação e verificação final são etapas separadas. O núcleo matemático permaneceu idêntico ao da base e o cache reproduziu seus coeficientes exatamente nos testes.

Consulte o [passo a passo com limitações](PASSO_A_PASSO.md), a [explicação de âncoras/autocorrelação e limites de x/y](ANCORAS_E_AUTOCORRELACAO.md), o [protocolo da extensão](PROTOCOL_CHURN_ANULAR.md), os [resultados detalhados](RESULTS.md), o [catálogo completo](results/extensao_churn_anular/catalogo_ancoras.csv), a [verificação da extensão](results/extensao_churn_anular/verificacao_final.json) e as [instruções de reprodução](README.md).
'''
    write('EXECUTIVE_SUMMARY.md',executive)
    write_supporting(lead,limits,main_table,sensitivity_table,c,timing,cs,coverage,common)
    write_pdf(lead,limits,c,timing)
    print('Relatórios integrados e PDF gerados; etapa 1 preservada.')

def write_supporting(lead,limits,main_table,sensitivity_table,c,timing,cs,coverage,common):
    steps=f'''# Passo a passo — extensão para churn e anular ondulado

{lead}

O [passo a passo da etapa 1](PASSO_A_PASSO_ETAPA1.md) preserva a execução original de plug/slug e a renomeação para v16.2_wco_eventos. A continuação abaixo foi autorizada pelo usuário na mesma pasta. Não foi alterado o projeto v16.1_perfis_wps.

1. **Preservação do estado anterior.** Os hashes do manifesto anterior foram conferidos e registrados em [estado_anterior.json](results/extensao_churn_anular/estado_anterior.json). Resumo e PDF originais receberam cópias de etapa 1. Tabelas, mapas, figuras originais, catálogo congelado e núcleo matemático foram preservados. **Limitação:** entrega local; src2 permanece ignorado pelo Git.

2. **Escopo e piloto óptico.** Foram mantidas as dez aquisições elegíveis de churn e nove de anular. O piloto usou m03q05_001, m05q40_001, m03q10_001 e m05q80_001; foram examinados panoramas, perfis e quadros de origem. **Limitação:** piloto no conjunto já conhecido, sem dados externos. [Registro do piloto](results/extensao_churn_anular/piloto_optico.csv).

3. **Definição da âncora.** Calculou-se a média das 54 linhas 18–71 inclusive. A suavização gaussiana com σ=16 px permitiu localizar extremos; proeminência mínima max(3; 0,25×(P95−P5)), janela de proeminência 1.025 px, largura à meia proeminência 16–512 px e separação de candidatos de 128 px. A suavização não foi aplicada à entrada wavelet. **Limitação:** o mínimo representa intensidade óptica, não uma interface física. [Protocolo](PROTOCOL_CHURN_ANULAR.md) e [configuração](config_ancoras_textura.json).

4. **Seleção sem métricas wavelet.** A região com suporte foi dividida em quatro estratos. Em cada tipo/estrato, escolheu-se a maior proeminência, com desempate pela menor coluna e separação de 512 px dos selecionados anteriores do mesmo tipo. Foram guardados todos os 1.365 candidatos e os motivos de não seleção. Resultado: quatro mínimos e quatro máximos em cada aquisição, 152 âncoras. **Limitação:** o limite de quatro por tipo é amostragem operacional; não inventário de eventos físicos.

5. **Rastreabilidade e revisão.** Foram produzidas e examinadas as 19 pranchas com oito marcos cada, panorama, perfil e quadro de origem. As 142 âncoras fora de mistura tiveram correspondência exata de pixels com a coluna de origem; dez estavam em zonas de mistura. Coordenadas alternativas σ12/20 e faixa 30–49 foram registradas. **Limitação:** revisão assistida por IA, sem validação independente ou rastreamento de estruturas entre quadros. [Revisão](results/extensao_churn_anular/revisao_visual.csv) e [catálogo](results/extensao_churn_anular/catalogo_ancoras.csv).

6. **Pareamento e congelamento.** Dos 8.208 pares candidatos, foram selecionados 830 por regras de largura/polaridade e hash determinístico. Mínimos compatíveis exigiram razão de larguras ≤1,50; distintos, ≥2,00; a faixa intermediária ficou catalogada. Mínimo × máximo teve um par por orientação. Regras, revisão, coordenadas, parceiros e fontes foram congelados antes da primeira nova WCO. **Limitação:** compatibilidade de largura não implica identidade física. [Congelamento](results/extensao_churn_anular/congelamento.json).

7. **Recorte e autocorrelação.** Cada âncora c definiu [c−1024,c+1024), mantendo 16 px de margem para perturbações. Só depois do recorte foram aplicados 255−I, remoção linear/média e padronização. Rxx/Ryy/Rxy usaram correlação linear completa e divisor 2.048, com 4.095 lags. **Limitação:** a correlação incorpora todo o contexto, inclusive estruturas vizinhas e costuras; zero do lag não é a coordenada global da âncora.

8. **Wavelets alternativas.** O núcleo original foi mantido: Morlet ω0=6, K=1, WPS=|CWT(Rxx)|, WCS=CWT(Rxy), WCO=|WCS|²/(WPSx(0)WPSy(0)), sem clipping. O resumo fez média dos 65 lags ±32 e depois das 200 escalas de λ=50–1.024 px. **Limitação:** o suporte geométrico e a média conservadora não são significância estatística nem calibração física.

9. **Validação e orçamento antes do lote.** Foram conferidos detector constante, centro sintético, translação, polaridade, soma direta da correlação, preparo por regressão independente, atraso conhecido, self e equivalência exata do cache com o núcleo nas oito variantes. Os controles sintéticos incluíram largura, fundo comum e inversão de sinal. O benchmark indicou cerca de 0,11 s/par principal e 0,18 s/par v5; a estimativa de 13–30 minutos e aproximadamente 0,5 GiB foi comunicada antes da execução completa. **Limitação:** testes validam operações, não a interpretação física. [Validação](results/extensao_churn_anular/validacao_numerica.json) e [benchmark](results/extensao_churn_anular/benchmark.json).

10. **Lote experimental, pausa e retomada.** Foram executadas oito variantes: principal, linhas 30–49, só remoção da média, centros opostos ±16 px nas duas orientações, detectores σ12/20 e núcleo literal v5. Das 6.640 avaliações ancoradas planejadas, {c['comparacoes_calculadas']} foram calculadas e {c['comparacoes_nao_estimaveis']} ficaram não estimáveis por ausência óptica em σ20. Foram reutilizados 513 controles fixos principais, calculadas 1.539 sensibilidades fixas e 171 controles de autocomparação/mesma aquisição. O usuário solicitou pausa após anchors/fixed; o checkpoint foi conferido antes de retomar somente diagnostics, sem repetir o lote principal. {cost_note(timing)} **Limitação:** as ausências não receberam centros substitutos; as variantes não foram escolhidas pelo maior escore. [Execução](results/extensao_churn_anular/execucao.json) e [checkpoint conferido na retomada](results/extensao_churn_anular/verificacao_checkpoint_retomada.json).

11. **Agregação e influência.** Primeiro foi feita a média por par de aquisições/grupo, depois por padrão; os contrastes usaram interseções de aquisições. As variantes também foram comparadas nos mesmos IDs de pares locais. A retirada de uma aquisição eliminou todos os seus pares, inclusive em uma análise explícita sem m03q40_001. **Limitação:** sobreposição e aquisições compartilhadas impedem interpretar os pares como réplicas independentes. Não foram produzidos valores-p.

12. **Relatórios e verificação final.** Foram integrados resumo executivo, gráficos, PDF, explicação matemática e arquivos de reprodução. Exemplos predefinidos foram separados dos extremos ilustrativos posteriores. A auditoria reproduziu o catálogo e parceiros, conferiu hashes, mapas, recálculos amostrais, agregações, vínculos e PDF. As figuras e páginas renderizadas foram inspecionadas. **Limitação:** consistência computacional não valida bolhas/ondas, rótulos locais nem generalização. [Verificação da extensão](results/extensao_churn_anular/verificacao_final.json) e [verificação integrada](results/verificacao_final.json).

## Limitações reunidas

{limits}
'''
    write('PASSO_A_PASSO.md',steps)
    coverage_table=table(['Padrão','Grupo','Pares locais','Pares de aquisições','Possíveis'],[
        [LABELS[r.classe],GROUPS[r.grupo],int(r.pares_ancoras),int(r.pares_aquisicoes),int(r.pares_aquisicoes_possiveis)] for r in coverage.itertuples()])
    write('RESULTS.md',f'''# Resultados detalhados — churn e anular ondulado

{lead}

Os [resultados originais de plug/slug](RESULTS_ETAPA1.md) permanecem separados. As definições de âncora e compatibilidade diferem; não foi feita uma média única entre as etapas.

## Cobertura do pareamento

{coverage_table}

## Contrastes principais

{main_table}

Somente os mesmos pares de aquisições entram nos dois lados de cada contraste. Por isso a média das âncoras pode mudar entre comparadores. A tabela numérica também preserva mínimos/máximos por par, a fração de diferenças positivas e o contraste no domínio conservador.

## Sensibilidades

{sensitivity_table}

As variantes de detector e de centros usam os controles fixos principais; faixa, remoção da média e v5 usam referências correspondentes recalculadas. A proteção de denominadores não descartou escalas no lote experimental. Nos pares ancorados, o máximo integral da Morlet foi {f(c['maximo_wco_morlet'],6)} e o do legado v5 foi {f(c['maximo_wco_v5'],6)}, sem clipping; o maior erro da equação foi {c['erro_equacao_maximo']:.2e}. As três pequenas ultrapassagens unitárias do legado (dois mapas ancorados e um fixo) ocorreram fora dos lags do resumo e foram reproduzidas por cálculo independente na [auditoria do v5](results/extensao_churn_anular/auditoria_ultrapassagens_v5.json).

## Tabelas auditáveis

- [Resultados de todos os pares e variantes](results/extensao_churn_anular/resultados_ancoras.csv), incluindo ausências explícitas.
- [Controles fixos calculados e reutilizados](results/extensao_churn_anular/resultados_fixos.csv).
- [Agregados por par de aquisições](results/extensao_churn_anular/agregados_ancoras.csv).
- [Contrastes por par](results/extensao_churn_anular/contrastes_pareados.csv) e [por padrão](results/extensao_churn_anular/contrastes_por_padrao.csv).
- [Sensibilidades no conjunto comum de pares](results/extensao_churn_anular/sensibilidade_conjunto_comum.csv), com principal recalculada no mesmo subconjunto, variante e mudança.
- [Painel das oito variantes na interseção dos mesmos 815 pares locais](results/extensao_churn_anular/sensibilidade_conjunto_global.csv), incluindo a linha principal no mesmo subconjunto.
- [Retirada de uma aquisição](results/extensao_churn_anular/retirada_aquisicao.csv), removendo todos os pares incidentes.
- [Média balanceada entre os dois padrões](results/extensao_churn_anular/contrastes_balanceados.csv); não substitui as conclusões por padrão.
- [Estabilidade óptica](results/extensao_churn_anular/estabilidade_ancoras.csv), [controles sintéticos](results/extensao_churn_anular/controles_sinteticos.csv) e [self/mesma aquisição](results/extensao_churn_anular/controles_mesma_aquisicao.csv).
- [Perfis por escala dos 830 pares principais](results/extensao_churn_anular/perfis_ancoras.npz) e [mapas completos dos exemplos](maps/extensao_churn_anular).

## Limitações

{limits}
''')
    write_method_note(c)
    write_readme(timing)

def write_method_note(c):
    old=(ROOT/'ANCORAS_E_AUTOCORRELACAO_ETAPA1.md').read_text()
    old=old.replace('A seção sobre churn e anular descreve uma possibilidade de extensão **ainda não executada**; não altera o catálogo nem os resultados existentes.',
        'A extensão para churn e anular ondulado **foi executada** com um segundo catálogo, descrito abaixo e no [protocolo da extensão](PROTOCOL_CHURN_ANULAR.md). O catálogo original de plug/slug permaneceu intacto.')
    old=old.replace('## Evento, âncora e limites dos sinais','## Evento, âncora e limites dos sinais\n\nA definição de evento por intervalo [L,R) abaixo pertence à **etapa 1: plug/slug**. Na extensão **churn/anular**, c é um mínimo local do perfil óptico suavizado; não se exige um evento físico delimitado. Os limites [c−1024,c+1024) e as equações de correlação/wavelet são os mesmos nas duas etapas.')
    old=old.replace('Imagens: detectar e revisar eventos sem consultar WCO','Imagens: localizar e revisar marcos ópticos sem consultar WCO')
    old=old.replace('Definir âncoras e parceiros pelo critério visual','Definir âncoras e parceiros pelo critério óptico')
    start=old.index('## Possível extensão a churn e anular — ainda não executada')
    end=old.index('## Fontes dentro da entrega',start)
    extension=f'''## Extensão executada: churn e anular ondulado

Foram localizados **76 mínimos principais e 76 máximos de controle nas 19 aquisições elegíveis**. O marco óptico é um extremo local da média das linhas 18–71, suavizada com gaussiana σ=16 px somente para localização. A proeminência e a largura foram limitadas previamente, e quatro estratos espaciais distribuíram os marcos ao longo do panorama. A seleção não consultou WCO. As 19 pranchas do catálogo foram revisadas visualmente; os 76 mínimos mantiveram correspondência em σ=12/20, com deslocamento máximo de 11 px.

**A suavização não alimenta a CWT.** Depois de achar c, o recorte vem da média original das linhas 18–71. A largura à meia proeminência serve para comparar marcos ópticos; não define os limites de x/y e não equivale a comprimento físico de bolha/onda.

| Padrão e âncora de exemplo | Centro c | Recorte [s,e) | Posição relativa |
| --- | --- | --- | --- |
| Agitado (churn), m03q05_001_MIN02256 | 2.256 | [1.232, 3.280) | −1.024…+1.023 px |
| Anular ondulado, m03q10_001_MIN03119 | 3.119 | [2.095, 4.143) | −1.024…+1.023 px |

Esses são exemplos individuais predefinidos de localização. O [catálogo de pares](results/extensao_churn_anular/pares_ancoras.csv) determina quais âncoras foram efetivamente comparadas. Dois marcos podem ter centros globais diferentes: cada um é centralizado na própria janela, sem deslocamento adicional para maximizar WCO.

### Autocorrelação de um contexto ancorado

Rxx continua comparando as 2.048 amostras de x com uma cópia deslocada de x. Ryy faz o mesmo com y; Rxy cruza os dois recortes. **A âncora não substitui o sinal por um impulso e não entra na soma como peso.** Ela só determina o conteúdo da janela. Rebatizar a origem de uma janela idêntica deixa sua autocorrelação igual; mover a janela troca amostras, tendência estimada e contexto, podendo alterar Rxx/Ryy. A correlação é linear, sem emenda circular.

O marco central é apenas uma parte do conteúdo: {c['ancoras_com_outros_marcos']} segmentos contêm outro marco selecionado e {c['ancoras_com_sobreposicao']} sobrepõem outro segmento selecionado. Portanto, a WCO caracteriza a região ancorada inteira. O ponto de maior coerência no mapa não localiza automaticamente a onda/bolha no panorama.

### Polaridade e limites dos resultados

Trocar y por −y inverte Rxy e WCS, mas deixa Ryy e |WCS|² iguais; a WCO permanece igual. Assim, um controle mínimo × máximo pode apresentar coerência alta. Além disso, ele pode alterar largura e forma, não somente polaridade. Esses controles são diagnósticos, sem obrigação matemática de coerência baixa.

Os 830 pares e oito variantes produziram {c['comparacoes_calculadas']} comparações calculáveis e {c['comparacoes_nao_estimaveis']} ausências ópticas explícitas em σ=20. A ancoragem mostrou viabilidade operacional, mas o ganho médio principal sobre janelas fixas foi pequeno em churn e quase nulo no anular. Não foi comprovada identidade física, superioridade geral ou capacidade de classificar os padrões.

Veja o [resumo com gráficos](EXECUTIVE_SUMMARY.md), o [protocolo específico](PROTOCOL_CHURN_ANULAR.md), a [estabilidade das âncoras](results/extensao_churn_anular/estabilidade_ancoras.csv) e as [coordenadas completas](results/extensao_churn_anular/catalogo_ancoras.csv). A [versão anterior desta nota](ANCORAS_E_AUTOCORRELACAO_ETAPA1.md) permanece como registro histórico da proposta.

'''
    old=old[:start]+extension+old[end:]
    old=old.replace('[Observabilidade por aquisição](results/observabilidade.csv): inclusões e limitações de churn/anular.',
        '[Observabilidade da etapa 1](results/observabilidade.csv): inclusões e limitações antes da extensão; para a situação atual, veja [cobertura da extensão](results/extensao_churn_anular/cobertura_ancoras.csv).')
    write('ANCORAS_E_AUTOCORRELACAO.md',old)

def write_readme(timing):
    write('README.md',f'''# Âncoras e wavelets alternativas — v16.2_wco_eventos

A pasta contém a etapa 1 de plug/slug e a extensão executada para **agitado (churn)** e **anular ondulado (wavy annular)**, ambas baseadas em **v16_spatial2_18_71**. Os arquivos numéricos originais foram preservados; os novos ficam em subpastas extensao_churn_anular.

Comece pelo [resumo executivo com gráficos](EXECUTIVE_SUMMARY.md), [PDF integrado](RESUMO_EXECUTIVO.pdf), [passo a passo com limitações](PASSO_A_PASSO.md) e [explicação de âncoras, autocorrelação e limites de x/y](ANCORAS_E_AUTOCORRELACAO.md). O [resumo original](EXECUTIVE_SUMMARY_ETAPA1.md) e a [reprodução da etapa 1](README_ETAPA1.md) continuam disponíveis.

## Arquivos principais

- [Protocolo congelado da extensão](PROTOCOL_CHURN_ANULAR.md) e [parâmetros](config_ancoras_textura.json).
- [Catálogo completo](results/extensao_churn_anular/catalogo_ancoras.csv): candidatos, seleção, centro, largura óptica, proeminência, estrato, padrão, recorte, costuras, procedência e coordenadas alternativas.
- [Pares selecionados](results/extensao_churn_anular/pares_ancoras.csv) e [todos os pares candidatos](results/extensao_churn_anular/pares_candidatos.csv).
- [Resultados detalhados](RESULTS.md), [resultados ancorados](results/extensao_churn_anular/resultados_ancoras.csv), [controles fixos](results/extensao_churn_anular/resultados_fixos.csv) e [contrastes](results/extensao_churn_anular/contrastes_por_padrao.csv).
- [Pranchas das 19 aquisições](figures/extensao_churn_anular/ancoras), [figuras principais](figures/extensao_churn_anular) e [coeficientes integrais dos exemplos](maps/extensao_churn_anular).
- [Congelamento](results/extensao_churn_anular/congelamento.json), [validação numérica](results/extensao_churn_anular/validacao_numerica.json), [verificação da extensão](results/extensao_churn_anular/verificacao_final.json), [verificação integrada](results/verificacao_final.json) e [manifesto de hashes](results/manifesto.json).

## Reprodução dos relatórios e conferência da entrega

Use Python 3.13.2 em `/Users/i521862/.pyenv/versions/3.13.2/bin/python3`, com as dependências de requirements.txt. O manifesto registra versões, caminhos e hashes. A distribuição PyWavelets declara 1.9.0, enquanto o módulo declara 1.8.0; essa discrepância foi preservada. A Morlet principal usa NumPy/SciPy. Dados e base anteriores precisam estar nos caminhos relativos originais.

Para refazer as agregações, figuras e relatórios a partir dos coeficientes já calculados, use a sequência abaixo na raiz do projeto TCC. Ela conserva os lotes numéricos e o histórico da pausa:

```sh
export MPLCONFIGDIR=/tmp/wco-eventos-mpl
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/summarize_texture_analysis.py
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/render_texture_results.py
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/write_texture_reports.py
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 src2/v16.2_wco_eventos/scripts/verify_final.py
```

O verificador integrado chama a verificação da extensão. Para somente conferir a entrega, basta o último comando: ele inclui reprodução temporária do catálogo e recálculos amostrais nas oito variantes. O manifesto inclui os arquivos de log; ao redirecionar a saída desse comando, use um arquivo fora desta pasta para que o log em gravação não invalide seu próprio hash. A geração do PDF integrado usa o utilitário local pdfunite; a verificação textual usa pdftotext, ambos disponíveis nesta execução.

O lote numérico integral está implementado em [catalog_texture_anchors.py](scripts/catalog_texture_anchors.py), [validate_texture_analysis.py](scripts/validate_texture_analysis.py) e [run_texture_analysis.py](scripts/run_texture_analysis.py), nessa ordem, seguido de [audit_texture_legacy.py](scripts/audit_texture_legacy.py). Uma execução integral deve usar uma cópia separada da pasta, mantendo os caminhos relativos da base e dos dados. Ela gera metadados e tempos próprios; não é a retomada histórica. Os verificadores finais desta entrega também exigem os registros e hashes daquela pausa/retomada, portanto uma nova execução precisa de registro e auditoria próprios, sem substituir nem flexibilizar a conferência histórica. A sequência acima é a reprodução direta suportada para a entrega atual.

{cost_note(timing)} [benchmark.json](results/extensao_churn_anular/benchmark.json) preserva o orçamento anterior. O histórico da pausa e retomada está em [RETOMADA.md](RETOMADA.md); o script resume_texture_analysis.py concluiu somente os diagnósticos, mantendo os metadados das etapas anteriores. O cache é limitado a 16 autocorrelações transformadas, conferidas contra o núcleo original. A revisão visual e o congelamento são insumos históricos: reproduzir o cálculo não cria nova validação independente. Não regenere o congelamento para contornar uma divergência; uma alteração metodológica exige outra versão e declaração de que os resultados atuais já são conhecidos.

## Conservação e escopo

Os relatórios de etapa 1 são arquivos históricos, com o escopo anterior à extensão. Os scripts atuais de relatórios mantêm a integração quando os resultados da extensão existem. O núcleo matemático, as entradas e os resultados originais permanecem preservados. Nenhum lote da pasta de origem é executado por estes scripts.

Todos os valores são descritivos de intensidade óptica, lag e escala em pixels. Não há identidade física comprovada, calibração para mm/tempo, validação externa nem classificação validada. Pares compartilham aquisições e contextos. src2 continua ignorado pelo Git; a entrega é local, sem commit ou sincronização presumidos.
''')

def text_page(pdf,title,blocks,filename=None):
    fig=plt.figure(figsize=(8.27,11.69));fig.patch.set_facecolor('white')
    fig.text(.08,.95,title,fontsize=23,weight='bold',va='top',color='#153c58')
    y=.84
    for heading,body in blocks:
        fig.text(.08,y,heading,fontsize=12,weight='bold',va='top',color='#147d64');y-=.03
        lines=textwrap.wrap(body.replace('**',''),width=87)
        fig.text(.08,y,'\n'.join(lines),fontsize=10.3,va='top',linespacing=1.4)
        y-=.021*len(lines)+.03
    assert y>=.035,(title,y)
    fig.text(.08,.025,'v16.2_wco_eventos • Dados, método e limitações: EXECUTIVE_SUMMARY.md',fontsize=8,color='#666')
    pdf.savefig(fig)
    if filename:fig.savefig(ROOT/F/filename,dpi=130)
    plt.close(fig)

def write_pdf(lead,limits,c,timing):
    extension_pdf=OUT/'resumo_extensao_para_integracao.pdf'
    with PdfPages(extension_pdf,metadata={'Title':'Âncoras em churn e anular ondulado — resumo executivo','Author':'Análise local do projeto TCC','Subject':'Método, resultados e limitações em PT-BR'}) as pdf:
        text_page(pdf,'Âncoras em churn\ne anular ondulado',[
            ('Resultado principal',lead),
            ('O que foi executado',f'19 aquisições; 152 âncoras (76 mínimos e 76 máximos); 830 pares locais; oito variantes. Foram calculadas {c["comparacoes_calculadas"]} comparações ancoradas; {c["comparacoes_nao_estimaveis"]} ficaram ausentes por falta de correspondência óptica em σ=20. Os parceiros foram congelados antes da nova WCO.'),
            ('Como interpretar','Âncora é uma coordenada de intensidade óptica. A análise usa a região inteira de 2.048 amostras ao redor dela. Localizar um mínimo de forma estável não comprova identidade de onda/bolha nem que o contexto apresentará maior coerência.'),
            ('Preservação','Esta extensão acrescenta churn e wavy annular à mesma pasta. Os resultados anteriores de plug/slug foram preservados. As sete páginas originais aparecem após um separador histórico no final deste PDF.')],filename='capa_resumo.png')
        text_page(pdf,'Qual padrão é qual\ne o que a CWT recebe',[
            ('Agitado (churn)','10 aquisições; roxo nos novos gráficos. Quatro mínimos de intensidade e quatro máximos por aquisição: 80 âncoras. x/y são sinais de aquisições diferentes nas comparações principais.'),
            ('Anular ondulado (wavy annular)','9 aquisições; verde nos novos gráficos. Quatro mínimos e quatro máximos por aquisição: 72 âncoras. Os nomes dos padrões são rótulos globais do inventário, não identificação física de cada marco.'),
            ('Etapa 1 preservada','Pistonado / bolhas alongadas (plug): 20 eventos aceitos em três aquisições. Golfadas (slug): seis eventos em três aquisições. Nessas classes a âncora foi o centro de um intervalo óptico alongado. Bolhas dispersas teve um único registro e não permite comparação ancorada entre aquisições.'),
            ('Do panorama aos mapas','Média 18–71 → localizar mínimo com σ=16 px → recortar o perfil original em [c−1024,c+1024) → inverter, remover tendência e padronizar → Rxx/Ryy/Rxy → CWT → WPS/WCS/WCO alternativas.'),
            ('Limites e leitura','Sinais: posição relativa −1.024…+1.023 px. Correlações: lags −2.047…+2.047 px. CWT: 200 comprimentos de onda de 50–1.024 px. Resumo: lags −32…+32. Nas curvas x é contínua e y tracejada; nos mapas a cor indica o valor da medida, conforme a barra lateral.')])
        figures=[('cobertura_estabilidade','Cobertura e estabilidade das âncoras'),('ancoras_dos_padroes','Imagens, localização e WPS dos dois padrões'),
            ('contrastes','Contrastes principais nas mesmas aquisições'),('entre_padroes','Coerência dentro e entre padrões'),('sensibilidades','Sensibilidades e cobertura dos controles'),
            ('retirada_aquisicao','Influência da retirada de uma aquisição'),('controles_sinteticos','Controles sintéticos e limites da WCO'),
            ('par_churn_predeterminado','Churn × churn: exemplo predefinido'),('par_anular_predeterminado','Anular × anular: exemplo predefinido')]
        for name,title in figures:
            img=plt.imread(ROOT/F/(name+'.png'));portrait=img.shape[0]>img.shape[1]
            fig=plt.figure(figsize=(8.27,11.69) if portrait else (11.69,8.27));ax=fig.add_axes([.025,.07,.95,.84]);ax.imshow(img);ax.axis('off')
            fig.suptitle(title,fontsize=16,weight='bold',y=.97,color='#153c58')
            fig.text(.035,.025,'Resultados descritivos. Pares compartilham aquisições; suporte geométrico não é significância estatística.',fontsize=8)
            pdf.savefig(fig);plt.close(fig)
        text_page(pdf,'Limitações e verificações',[
            ('Representação óptica','A âncora não mede espessura de filme nem certifica a mesma onda ou bolha. O detector foi desenvolvido no conjunto já conhecido, sem especialista independente ou teste externo. A troca para linhas 30–49 deslocou mínimos de churn em até 47 px.'),
            ('Contexto e montagem',f'{c["ancoras_com_sobreposicao"]} dos 152 recortes sobrepõem outro recorte selecionado; {c["ancoras_com_outros_marcos"]} contêm outro marco selecionado. Dez âncoras estão em mistura de costuras. m03q40_001 mantém 60 deslocamentos não positivos herdados. Nenhum panorama foi reconstruído.'),
            ('Controles e dependência','Larguras compatíveis: razão ≤1,50; distintas: ≥2,00. Mínimo × máximo também pode diferir em forma/largura. A WCO é invariável à inversão global de sinal. Três mapas do legado v5 excederam ligeiramente 1 (até 1,006157), fora de ±32 px, confirmados por convolução independente. A Morlet principal permaneceu ≤1. Não há valores-p ou classificação validada.'),
            ('Operações conferidas','Cache idêntico ao núcleo original nos testes; correlação por soma direta; preparo por regressão independente; centro e translação do detector; atraso conhecido de +100 px em y recuperado em −100 px; self em lag zero igual a 1. Esses testes validam operações, não identidade física.'),
            ('Custo observado',cost_note(timing)+' Orçamento prévio: aproximadamente 13–30 minutos. Figuras, redação e verificação final foram etapas separadas.')])
        text_page(pdf,'Apêndice histórico\nEtapa 1: plug e slug',[
            ('Escopo das próximas sete páginas','Reprodução do resumo anterior à extensão. Quando o texto original afirma que churn/anular não tiveram eventos ancorados, descreve exclusivamente aquela etapa. A extensão executada e suas limitações estão nas páginas anteriores.'),
            ('Resultados preservados','Plug: diferenças de 0,0275 frente a eventos distintos e 0,0289 frente a fixos. Slug: 0,1520 e 0,2192, sustentados por somente um par equivalente; o contraste frente a distintos caiu a −0,0427 nas linhas 30–49.'),
            ('Definições diferentes','Plug/slug usaram o centro de intervalos ópticos alongados. Churn/anular usaram mínimos locais de intensidade. As duas etapas não foram combinadas numa única estimativa de desempenho.')])
    subprocess.run(['pdfunite',str(extension_pdf),str(ROOT/'RESUMO_EXECUTIVO_ETAPA1.pdf'),str(ROOT/'RESUMO_EXECUTIVO.pdf')],check=True)

if __name__=='__main__':main()
