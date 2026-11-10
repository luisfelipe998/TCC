"""Redige entregáveis em PT-BR com números lidos dos resultados executados."""
import json
from pathlib import Path
import pandas as pd
from profile_core import ROOT, BASE
from reading_guide import build_summary


def f(x, digits=3):
    return f'{x:.{digits}f}'.replace('.',',')


def table(headers, rows):
    return '| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in rows)+'\n'


def main():
    out=ROOT/'results'
    r=pd.read_csv(out/'ranking_summary.csv').set_index('metodo')
    pats=pd.read_csv(out/'ranking_by_pattern.csv')
    influence=pd.read_csv(out/'acquisition_deletion.csv')
    sens=pd.read_csv(out/'within_acquisition_sensitivity.csv')
    meta=json.loads((out/'run_metadata.json').read_text())
    pre=json.loads((out/'preflight.json').read_text())
    ops=json.loads((out/'operations_validation.json').read_text())
    cfg=json.loads((ROOT/'config.json').read_text())
    a,b=r.loc['wps_principal'],r.loc['wco_v16_18_71']
    deletion=influence[influence.metodo=='wps_principal']
    bubble=deletion[deletion.aquisicao_retirada=='m07q01_001'].iloc[0]
    effects=[]
    for test,ref in [('wps_principal','wco_v16_18_71'),('wps_principal','wco_v5_20_59'),
                     ('wps_prefixo_par','wco_v16_18_71'),('v5_prefixo_par','wco_v5_18_71'),
                     ('wps_principal','v5_mesmas_escalas')]:
        effects.append(dict(metodo=test,referencia=ref,
            delta_A_pares=float(r.loc[test,'A_pares']-r.loc[ref,'A_pares']),
            delta_A_balanceada=float(r.loc[test,'A_balanceada']-r.loc[ref,'A_balanceada']),
            variacao_relativa_A_pares_pct=float(100*(r.loc[test,'A_pares']/r.loc[ref,'A_pares']-1)),
            variacao_relativa_A_balanceada_pct=float(100*(r.loc[test,'A_balanceada']/r.loc[ref,'A_balanceada']-1))))
    pd.DataFrame(effects).to_csv(out/'ranking_effects.csv',index=False)
    summary=dict(conclusao='A comparação de formas de WPS alternativa foi mais útil para ordenar pares neste conjunto que o escore WCO herdado. O ganho é interno e descritivo; não demonstramos melhora robusta e geral de classificação entre condições experimentais.',
        base='v16_spatial2_18_71',linhas=[18,71],aquisicoes=26,pares_intra=87,pares_entre=238,
        principal=a.to_dict(),referencia_wco=b.to_dict(),efeitos=effects,
        bolhas_sem_replica=True,exclusao_descritiva_bolhas=bubble.to_dict(),
        limitacoes=['Pares dependentes','Avaliação interna previamente explorada','Sensibilidade ao comprimento e à posição',
                    'Uma aquisição de bolhas dispersas','Três aquisições de plug e de slug',
                    'Condições operacionais e comprimentos associados aos rótulos','Reconstrução espacial sem validação física independente'],
        metrica='A é concordância de ordenação, não acurácia nem WCO')
    (out/'executive_metrics.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
    methods=[('wps_principal','Perfil principal, registro inteiro'),('wps_prefixo_par','Perfil v16, prefixo de cada par'),
             ('v5_prefixo_par','Perfil v5, prefixo de cada par'),('wco_v16_18_71','Escore WCO v16, linhas 18–71'),
             ('wco_v5_18_71','Escore WCO v5, linhas 18–71'),('wco_v5_20_59','Escore WCO v5 original, linhas 20–59')]
    pattern_rows=[]
    for pattern in ['Bolhas dispersas','Bolhas alongadas (plug)','Golfadas (slug)','Agitado (churn)','Anular ondulado']:
        x=pats[(pats.metodo=='wps_principal')&(pats.padrao==pattern)].iloc[0]
        y=pats[(pats.metodo=='wco_v16_18_71')&(pats.padrao==pattern)].iloc[0]
        pattern_rows.append((pattern,int(x.pares_intra),'Não estimável' if x.pares_intra==0 else f(x.A_pares),
                             'Não estimável' if x.pares_intra==0 else f(x.A_balanceada),
                             'Não estimável' if x.pares_intra==0 else f(y.A_pares)))
    pattern_table=table(['Padrão','Pares intra','Perfil: A por par','Perfil: A balanceada','WCO: A por par'],pattern_rows)
    exec_text = build_summary()
    (ROOT/'EXECUTIVE_SUMMARY.md').write_text(exec_text)
    names={m:name for m,name in methods}
    names.update({k:v['label'] for k,v in meta['method_specs'].items()})
    names.update(wco_v16_20_59='Escore WCO v16, linhas 20–59')
    detailed_table=table(['Método','A por par','A balanceada','Média intra','Média entre'],
        [(names.get(m,m),f(x.A_pares,6),f(x.A_balanceada,6),f(x.media_intra,6),f(x.media_entre,6)) for m,x in r.iterrows()])
    sensitivity_table=table(['Comparação na mesma aquisição','H mediano','H máximo'],
        [(names.get(v,v),f(g.H_mesma_aquisicao.median()),f(g.H_mesma_aquisicao.max())) for v,g in sens.groupby('variante')])
    effects_table=table(['Comparação','ΔA por par','ΔA balanceada','Variação relativa A por par'],
        [(f"{names.get(x['metodo'],x['metodo'])} versus {names.get(x['referencia'],x['referencia'])}",
          f(x['delta_A_pares'],6),f(x['delta_A_balanceada'],6),f(x['variacao_relativa_A_pares_pct'],2)+'%') for x in effects])
    report=f'''# Relatório de execução e resultados

## Contrato matemático e fonte primária

Foram conferidas visualmente as páginas impressas 411, 413 e 414 de Schuck Jr. e Bodmann (2019), capítulo 32, DOI 10.1007/978-3-030-16077-7_32, no [PDF local](../../assets/schuck-bardo-wavelet-coherence.pdf). As imagens preservadas constam da [auditoria da base](../v16_spatial2_18_71/EQUATION_AUDIT.md). A CWT adotou normalização 1/√a, Morlet π^(−1/4) exp(−t²/2) exp(i6t), e λ=4πa/(6+√38).

Para K=1, Pₓ(b,a)=|CWT(Rₓₓ)(b,a)|; Cₓᵧ(b,a)=CWT(Rₓᵧ)(b,a); Γₓᵧ(b,a)=|Cₓᵧ(b,a)|²/[Pₓ(0,a)Pᵧ(0,a)]. P é a WPS alternativa, C a WCS alternativa complexa e Γ a WCO alternativa. A correlação cruzada somou x[n+m] conjugado(y[n]); um atraso de y produziu lag negativo. O divisor N foi a adaptação enviesada já autorizada, distinta de N−|m| das Eqs. 32.21–32.23. Não houve ensemble entre aquisições, suavização convencional de Grinsted, alteração do denominador ou clipping.

## Fontes verificadas e reaproveitamento

Foram verificados {pre['manifest_entries_checked']} registros do manifesto da base e congelados hashes de {pre['source_files']} arquivos. Os 31 PNGs e 31 arquivos XLSX coincidiram com os hashes auditados. A extração 18–71 foi recalculada a partir de cada PNG e comparada exatamente ao sinal em cache; 20–59 foi conferida apenas como referência da planilha original. A reconstrução pixel a pixel e a estimação dos deslocamentos são evidências herdadas, não repetidas nesta execução. A reprodução computacional não comprova validade física da costura.

Reutilizaram-se perfis por aquisição e segmentos de N=3.171 da base indicada. A coluna de lag zero e a média integral de todas as 26 aquisições foram recalculadas para verificar o reaproveitamento: erros máximos {meta['max_reused_zero_error']:.3e} e {meta['max_reused_mean_error']:.3e}, respectivamente. As médias com suporte geométrico ausentes na versão atual foram calculadas novamente. Os perfis de 350 combinações únicas de aquisição/comprimento atenderam aos 325 prefixos comuns por par. Não se confundiram perfis de registros inteiros com os prefixos desses pares.

## Representação, quadratura e domínio

Na faixa principal, uⱼ=log(λⱼ). Os pesos trapezoidais foram w₀=(u₁−u₀)/2, wⱼ=(uⱼ₊₁−uⱼ₋₁)/2 para pontos internos e wúltimo=(uúltimo−uantepenúltimo)/2, entendendo o segundo termo como o ponto imediatamente anterior ao último. De forma inequívoca: para J pontos, w[J−1]=(u[J−1]−u[J−2])/2. A área foi ΣwⱼPⱼ e as massas pⱼ=wⱼPⱼ/ΣwₖPₖ. A densidade desenhada foi Pⱼ/área. As massas, densidades, valores brutos, pesos e áreas estão em profiles_normalized.csv e profiles_raw.csv.

Hellinger foi √[Σ(√pⱼ−√qⱼ)²/2]. O cosseno foi 1−ΣwⱼPⱼQⱼ/√[(ΣwⱼPⱼ²)(ΣwⱼQⱼ²)]. A normalização adicional removeu a amplitude da WPS, embora o sinal já tivesse variância unitária. Nenhuma massa foi chamada fração de energia. Valores não finitos, negativos ou de área nula foram rejeitados.

A condição N−1≥√2a em lag zero do menor registro definiu a faixa principal de 187 pontos, 50–2.295,940454 px. A faixa conservadora teve 117 pontos, 50–543,851438 px, segundo 3a≤N/2. A faixa integral reteve os 200 pontos, 50–3.000 px. A borda superior principal foi geometricamente definida e não escolhida pelos rótulos. Há no máximo cerca de 1,38 comprimentos de onda da maior escala no menor registro; suporte geométrico não resolve a escassez de estruturas longas.

Para as médias em lag, as máscaras foram N−1−|b|≥√2a e |b|+3a≤N/2. Cada média dividiu a soma pelo número real de lags válidos; zero lags resultou em valor ausente, não zero artificial. As contagens completas estão em [lag_support_counts.csv](results/lag_support_counts.csv). A média conservadora foi comparada apenas em sua faixa comum menor. Esses controles foram novos e explícitos; os escores WCO herdados permaneceram sem restrição geométrica.

## Desfecho, pesos e referências

Para valores de semelhança S, A foi a média ponderada de 1[S_intra>S_entre]+0,5·1[S_intra=S_entre]. Usou-se S=−H ou −distância cosseno para perfis e S=escore WCO para a referência. A é equivalente à interpretação de ordenação de uma área ROC, sem classificador treinado ou limiar de decisão.

O peso por par deu o mesmo peso aos 87 pares internos e, separadamente, aos 238 controles. O peso balanceado deu 1/4 a cada padrão interno e 1/10 a cada combinação cruzada, dividindo igualmente entre pares dentro de cada grupo. Por padrão, usaram-se apenas seus pares internos e controles que o envolviam, com 1/4 para cada outro padrão. Esse resumo local não é a decomposição exata de A global, pois muda o conjunto de controles.

O escore WCO herdado foi Σa wₐ max_b Γ(b,a)/Σa wₐ, com wₐ=√[Pₓ,ref(0,a)Pᵧ,ref(0,a)] quando ambas as WPS atingiram 5% de seus próprios picos e os denominadores passaram o piso relativo 10⁻¹⁰; fora disso, peso zero. Os máximos vieram de todos os lags. Não é média simples de WCO nem nova equação de Schuck. A coluna exata foi weighted_scale_peak_power05. Pesos nativos 18–71 foram usados nas comparações entre implementações da base nova; pesos nativos 20–59 acompanharam a referência histórica. A comparação entre referências de entradas distintas não isolou somente a extração.

{detailed_table}
As médias intra/entre são distâncias para os perfis e escores para WCO: **não são comparáveis numericamente entre essas duas famílias**. A é o desfecho comum.

{effects_table}
Variações relativas de A são razões aritméticas da mesma estatística, não porcentagens de melhora física ou de acurácia. A principal supera a referência WCO atual, mas o perfil v5 com o mesmo eixo atingiu desempenho praticamente idêntico. A razão de ganho não deve ser atribuída à correção numérica da Morlet.

## Variação por padrão, comprimento e posição

{pattern_table}
{sensitivity_table}
As distâncias entre segmentos e registro inteiro chegaram à ordem das distâncias entre aquisições: a posição afetou a representação. Todas as janelas foram recortadas antes do preparo, preservando 1 px/amostra, sem deformação dos eixos ou procura do trecho mais favorável. A aquisição curta de bolhas gerou três segmentos idênticos, que permaneceram uma única aquisição.

![Comprimentos e distâncias](figures/06_comprimento_e_confundimento.png)

Todos os anulares apresentaram N≥14.224, enquanto os demais chegaram a N≤10.349. A associação descritiva H versus log da razão dos comprimentos foi {f(meta['H_vs_log_length_ratio_spearman_descriptive'])}. Não houve inferência causal nem teste de significância dessa correlação entre pares dependentes. O controle de comprimento igual reduziu parte da separação; não demonstrou ausência de outros fatores de confusão.

## Influência por aquisição e ausência de inferência indevida

Cada uma das 26 retiradas removeu a aquisição e todos os seus pares, deixando 300 pares. Os pesos foram recalculados. A principal por par variou de {f(deletion.A_pares.min(),6)} a {f(deletion.A_pares.max(),6)}; a balanceada, de {f(deletion.A_balanceada.min(),6)} a {f(deletion.A_balanceada.max(),6)}. A retirada de bolhas reduziu de dez para seis as combinações cruzadas e produziu o menor A balanceado. Não se eliminaram bolhas do método principal.

As retiradas mediram influência interna, não erro de generalização. Não houve treinamento, ajuste aos rótulos, intervalo de confiança com pares pseudorreplicados ou permutação de rótulos sem justificativa de permutabilidade. Condições operacionais não foram usadas como estratos estatísticos por falta de um desenho de réplicas independentes validado para esse fim. A seleção de exemplos desfavoráveis ocorreu após os resultados e foi identificada como tal.

## Controles e validação

O controle com componentes de 100 e 500 px produziu H={f(ops['H_componentes_100_500'],6)}. A mudança global de amplitude produziu erro de H={ops['H_amplitude_sinal']:.3e}. No controle de fases distintas, H={f(ops['fases_distintas']['H'],6)} e escore WCO={f(ops['fases_distintas']['escore_WCO'],6)}. Os espectros circulares eram iguais antes do preparo; a autocorrelação linear com bordas e a remoção de tendência não preservam necessariamente identidade exata dos perfis.

![Controles sintéticos](figures/09_controles_sinteticos.png)

Os novos testes cobriram quadratura contra integral analítica, amplitude, identidade, suportes disjuntos, degenerescência, contagens de suporte, ordenação e atraso conhecido. A soma inversa de Fourier usada para calcular somente o lag zero foi comparada à CWT integral em N=3.171, 9.586 e 23.079; erro máximo inferior a 10⁻¹⁵. A otimização alterou o custo, não a equação. O lote principal durou {f(meta['elapsed_seconds'],1)} s, antes das figuras e da auditoria final; a estimativa prévia foi de 120–240 s.

## Limitações e conclusão

Persistiram as limitações do resumo executivo e da reconstrução v11. A auditoria herdada identificou 60/500 deslocamentos não positivos em m03q40_001 e 1/500 em m04q01_001; ambas permaneceram incluídas. Não foi provada a equivalência do panorama a uma fotografia espacial instantânea. Intensidade não foi convertida em fração de vazio. Pixels não foram convertidos em milímetros, frequência temporal ou velocidade.

O ganho de ordenação dos perfis foi claro dentro deste conjunto e persistiu no controle de prefixos por par. Houve perdas locais, dependência da posição e influência substancial do único registro de bolhas no balanceamento. Portanto, não demonstramos melhora robusta e geral entre condições experimentais, nem desempenho de classificação externa. As novas distâncias complementaram a WCS/WCO alternativa sem redefini-las.

As tabelas completas e scripts estão em [README.md](README.md). A execução, os resultados herdados e os novos cálculos foram separados em [PASSO_A_PASSO.md](PASSO_A_PASSO.md). A verificação final e os hashes constam de [final_verification.json](results/final_verification.json) e [artifact_manifest.json](results/artifact_manifest.json).
'''
    # Formulação inequívoca do último peso, sem ambiguidade de nomes de índices.
    report=report.replace('wúltimo=(uúltimo−uantepenúltimo)/2, entendendo o segundo termo como o ponto imediatamente anterior ao último. De forma inequívoca: para J pontos, w[J−1]=(u[J−1]−u[J−2])/2.',
                          'para J pontos, w[J−1]=(u[J−1]−u[J−2])/2 na extremidade final.')
    (ROOT/'REPORT.md').write_text(report)
    print('Resumo executivo, relatório e métricas estruturadas gravados.')


if __name__=='__main__':
    main()
