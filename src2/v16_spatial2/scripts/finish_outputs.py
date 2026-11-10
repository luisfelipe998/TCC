"""Build readable, data-driven summaries and representative panorama figures."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from wavelet_core import *
from plot_labels import PATTERN_PT, METHOD_PT

ROOT=Path(__file__).resolve().parents[1]
LAM=np.geomspace(50,3000,200)
SENSITIVITY_PT={
    'center_vs_core_mean_zero_lag':'Linha central versus média das linhas 20–59',
    'demean_vs_linear_detrend_zero_lag':'Remoção da média versus tendência linear',
    'lower_vs_core_mean_zero_lag':'Faixa inferior versus média das linhas 20–59',
    'upper_vs_core_mean_zero_lag':'Faixa superior versus média das linhas 20–59',
    'pywt_precision_12_same_scales':'PyWavelets, precisão 12, mesmas escalas',
    'pywt_precision_16_same_scales':'PyWavelets, precisão 16, mesmas escalas',
}

def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+
                     ['| '+' | '.join(map(str,row))+' |' for row in rows])

def representative_figures(inv):
    folder=ROOT/'figures/representatives';folder.mkdir(exist_ok=True)
    links=[]
    for pattern,g in inv[inv.included].groupby('pattern'):
        row=g.sort_values('run').iloc[0];name=row.run
        image=np.asarray(Image.open(ROOT.parent/'data_processed'/name/'__panorama_v11.png').convert('L'))
        raw=np.load(ROOT/'intermediate'/f'{name}_input.npz')['raw'];x=prepare(raw);n=len(x)
        w=morlet_cwt(correlation(x,x),LAM/FOURIER_FACTOR);mag=abs(w)
        fig=plt.figure(figsize=(14,9),layout='constrained')
        gs=fig.add_gridspec(3,2,height_ratios=[1,1,3],width_ratios=[5,1.4])
        ax=fig.add_subplot(gs[0,:]);step=max(1,int(np.ceil(n/2000)))
        ax.imshow(image[:,::step],cmap='gray',aspect='auto',extent=[0,n,80,0],vmin=0,vmax=255)
        ax.axhline(20,color='tab:orange',lw=.8);ax.axhline(60,color='tab:orange',lw=.8)
        ax.set(xlabel='Posição no panorama (px)',ylabel='Linha da imagem',title='Panorama original; laranja: linhas de extração 20–59')
        ax=fig.add_subplot(gs[1,:]);ax.plot(np.arange(n),x,lw=.6)
        ax.set(xlabel='Posição no panorama (px)',ylabel='Intensidade padronizada',title='Sinal após remoção da tendência linear e normalização')
        ax=fig.add_subplot(gs[2,0]);lag=np.arange(-n+1,n);step=max(1,int(np.ceil(len(lag)/1800)))
        idx=np.unique(np.r_[np.arange(0,len(lag),step),n-1,len(lag)-1])
        im=ax.pcolormesh(lag[idx],LAM,mag[:,idx],shading='auto',cmap='viridis')
        ax.set(xlabel='Defasagem da correlação (px)',ylabel='Comprimento de onda (px)',yscale='log',ylim=(3000,50),title='alternative WPS = módulo da CWT da autocorrelação enviesada')
        ax.set_yticks([50,100,200,500,1000,2000,3000],labels=['50','100','200','500','1000','2000','3000'])
        fig.colorbar(im,ax=ax,label='Módulo')
        ax=fig.add_subplot(gs[2,1]);ax.plot(mag.mean(axis=1),LAM,label='Todas as defasagens',color='tab:blue')
        ax.set(xlabel='Módulo médio',yscale='log',ylim=(3000,50),title='Média da\nalternative WPS');ax.legend(fontsize=8)
        ax.set_yticks([])
        suffix='\nUm único registro elegível: sem par intrapadrão independente de alternative WCS/alternative WCO' if len(g)==1 else ''
        fig.suptitle(f'{name} | {PATTERN_PT[pattern]} | N={n}{suffix}',fontsize=12)
        target=folder/f'{name}.png';fig.savefig(target,dpi=130);plt.close(fig)
        links.append((pattern,name,f'figures/representatives/{name}.png'))
    return links

def main():
    inv=pd.read_csv(ROOT/'results/input_audit.csv')
    per=pd.read_csv(ROOT/'results/per_run_methods.csv');r=per[per.included]
    pairs=pd.read_csv(ROOT/'results/pair_methods.csv')
    null=pd.read_csv(ROOT/'results/null_diagnostic.csv')
    assert len(null)==5 and (null.surrogates==199).all(), 'O diagnóstico nulo completo precisa terminar antes.'
    assert 'domain' in json.loads((ROOT/'results/null_metadata.json').read_text()), 'Recalcule o diagnóstico nulo sem máscaras antes do relatório.'
    assert len(pairs)==261 and len(inv)==31 and inv.included.sum()==26
    sensitivity=pd.read_csv(ROOT/'results/sensitivity.csv')
    sensitivity=sensitivity[sensitivity.run.isin(inv[inv.included].run)]
    seam=pd.read_csv(ROOT/'results/seam_sensitivity.csv');seam=seam[seam.included]
    profiles=pd.read_csv(ROOT/'results/profiles.csv')
    assert 'max_wco_supported' not in pairs and 'lag_mean_coi' not in profiles, 'Atualize as tabelas sem máscaras antes dos relatórios.'
    p=profiles[(profiles.included)&(profiles.method=='v16_biased')]
    fig,axes=plt.subplots(1,2,figsize=(12,4),layout='constrained')
    for pattern,g in p.groupby('pattern'):
        for ax,field in zip(axes,['zero_lag_wps','lag_mean_all']):
            pivot=g.pivot(index='run',columns='wavelength_px',values=field)
            z=pivot.div(pivot.max(axis=1),axis=0)
            line=ax.plot(pivot.columns,z.median(axis=0),label=f'{PATTERN_PT[pattern]} (n={len(z)})')[0]
            if len(z)>1:ax.fill_between(pivot.columns,z.min(axis=0),z.max(axis=0),color=line.get_color(),alpha=.08)
    for ax,title in zip(axes,['alternative WPS na defasagem zero','Média da alternative WPS em todas as defasagens']):
        ax.set(xlabel='Comprimento de onda (px)',ylabel='Módulo normalizado pelo pico de cada registro',xscale='log',title=title)
        ax.legend(fontsize=7);ax.grid(alpha=.2)
    fig.suptitle('Perfis por padrão: mediana e intervalo observado, registros completos')
    fig.savefig(ROOT/'figures/pattern_profiles.png',dpi=150);plt.close(fig)
    links=representative_figures(inv)
    selfrows=[]
    for method,g in r.groupby('method'):
        selfrows.append([METHOD_PT[method],f'{g.self_wco_max.median():.4f}',f'{g.self_wco_max.max():.4f}'])
    pairrows=[]
    for method,g in pairs.groupby('method'):
        pairrows.append([METHOD_PT[method],f'{g.max_wco.median():.4f}',f'{g.max_wco.max():.4f}',
                         str(int((g.max_wco>1+1e-8).sum()))])
    peaks=[]
    for pattern,g in r[r.method=='v16_biased'].groupby('pattern'):
        peaks.append([PATTERN_PT[pattern],len(g),f'{g.peak_zero_lag_px.min():.1f}–{g.peak_zero_lag_px.max():.1f}',
                      f'{g.peak_zero_lag_px.median():.1f}',f'{g.peak_lag_mean_px.median():.1f}'])
    sensitivity_rows=[]
    for name,g in sensitivity.groupby('comparison'):
        sensitivity_rows.append([SENSITIVITY_PT.get(name,name),f'{100*g.relative_l2.median():.3f}%',f'{100*g.relative_l2.max():.3f}%'])
    nullrows=[[PATTERN_PT.get(x.pattern,'Controle entre padrões'),f'{x.run_x} / {x.run_y}',f'{x.observed_max:.4f}',f'{x.null_q95:.4f}',
               f'{x.conditional_monte_carlo_p:.3f}'] for x in null.itertuples()]
    body='''# Resultados calculados — v16

Esta síntese foi gerada diretamente das tabelas atualizadas. As figuras e as estatísticas atuais usam o plano completo, sem COI nem máscara conservadora. As médias da alternative WPS incluem todos os lags; o denominador da alternative WCO permanece no lag zero. As decisões e demonstrações estão em [REPORT.md](REPORT.md). A versão anterior está no [arquivo histórico com máscaras](results/historico_com_mascaras.zip).

## Cobertura e reprodução

- 31 panoramas reconstruídos e comparados pixel a pixel; 31 planilhas conferidas contra a média das linhas 20–59.
- 26 registros no escopo confirmado; 87 pares intrapadrão; três versões por par.
- Os cinco registros excluídos continuam na auditoria individual, mas não nas tabelas por padrão abaixo.
- 200 comprimentos de onda logaritmicamente espaçados, de 50 a 3000 px; cálculos feitos em todos os lags, sem redução da resolução do sinal.

## Autocoerência: 26 registros elegíveis

O valor em lag zero foi 1 em todos os métodos, por identidade algébrica. As violações abaixo são máximos buscados sobre outros lags e escalas; não são erros dessa identidade.

'''
    body+=table(['Método','Mediana do máximo no plano completo','Maior máximo no plano completo'],selfrows)
    body+='''

A Morlet analítica com correlação enviesada respeitou o limite no plano completo, sem limitar artificialmente os resultados a 1. O PyWavelets do v5 e o estimador não enviesado apresentaram violações. Retirar as máscaras não corrigiu esses extremos nem modificou as equações; apenas deixou o domínio integral explícito nas figuras e estatísticas.

## alternative WCO dos 87 pares

'''
    body+=table(['Método','Mediana do máximo no plano completo','Maior máximo','Pares com máximo > 1'],pairrows)
    body+='''

As medianas de alternative WCO no lag zero ao longo das 200 escalas são baixas; não se confirmou a descrição geral de alternative WCO 0,7–1,0 entre todas as aquisições do mesmo padrão nas notas antigas. Isso não exclui valores altos em bandas específicas de alguns pares. A comparação detalhada encontra-se em [pair_methods.csv](results/pair_methods.csv). O comparador v5 conserva o pré-processamento antes do corte; seu lag é reportado na convenção do artigo após a transformação.

## Escalas dos perfis

'''
    body+=table(['Padrão','Registros','Faixa do pico alternative WPS(0,a), px','Mediana desse pico, px','Mediana do pico da média em lag, px'],peaks)
    body+='''

Os picos de alternative WPS em zero e da média em lag podem ser diferentes no mesmo sinal. Não são duas estimativas intercambiáveis do comprimento de uma bolha. Não houve ajuste das bandas para separar classes, e [pattern_profiles.png](figures/pattern_profiles.png) mostra sobreposição e variabilidade. Os picos são descrições do intervalo investigado, sem teste de significância de cada máximo.

## Sensibilidade da entrada e da CWT

Os erros relativos abaixo são normas L2 entre perfis de alternative WPS no lag zero; não são erros percentuais de cada pixel nem estimativas de incerteza física. A comparação de precisão usa escalas idênticas e corrige a constante de normalização da cmor.

'''
    body+=table(['Comparação','Mediana do erro relativo','Maior erro relativo'],sensitivity_rows)
    body+=f'''

No teste de interpolação de vizinhanças das emendas, {int(seam.interpolation_test_performed.sum())} dos 26 registros permitiram o teste com menos de metade dos pixels alterada. O maior erro relativo do perfil foi {100*seam.seam_interpolation_zero_lag_relative_l2.max():.3f}%. Isso não limita o erro causado por um deslocamento incorreto: interpolar uma emenda não recompõe uma estrutura perdida ou duplicada.

O controle [window_sensitivity.png](figures/diagnostics/window_sensitivity.png) compara começo, centro e fim com N=3171 px, preservando a amostragem. As três janelas de um registro não são réplicas independentes.

## Diagnóstico nulo condicional

Foram recalculados 199 pares de fases aleatórias por exemplo. A estatística foi o máximo em todos os lags e nas 200 escalas de 50–3000 px, aplicando o mesmo domínio ao par observado e a cada realização aleatória. Os valores anteriores, calculados com máscara, foram substituídos.

'''
    body+=table(['Exemplo','Par','Máximo observado','Quantil 95% do modelo','Valor-p condicional'],nullrows)
    body+='''

Esses valores dependem do modelo estacionário de fases aleatórias. Não demonstram sincronização entre aquisições nem identificam causalidade. Não foram feitas decisões de rejeição corrigidas entre os cinco exemplos. A resolução Monte Carlo é 0,005 e nenhum contorno desse diagnóstico foi transferido para os mapas principais. Um máximo incomum sob esse modelo, quando observado, não estabelece ligação física entre registros independentes.

## Verificações e conclusão

Os testes numéricos, a comparação com funções extraídas do próprio v5, os hashes e a cobertura dos artefatos estão nos JSONs de validação em results. O controle de três realizações recuperou −500 amostras de lag e demonstrou que média de módulo e módulo da média não são intercambiáveis.

**Recomendação:** usar a correlação enviesada com a Morlet e a conversão de escala de Schuck, conservar a versão não enviesada como controle e interpretar alternative WPS/alternative WCS/alternative WCO como descritores de lag–escala. A melhoria demonstrada é matemática e numérica; não foi demonstrado ganho de classificação ou reconstrução física do escoamento.
'''
    (ROOT/'RESULTS.md').write_text(body)
    readme='''# v16_spatial2 — análise espacial de Schuck

**Estado atual:** todas as figuras com nomes dos padrões em inglês e demais textos em PT-BR, sem COI nem máscara conservadora. As estatísticas dependentes dessas restrições foram recalculadas no domínio integral; as equações alternativas foram mantidas.

**Entrada:** panorama v11, média das linhas 20–59. **Escopo confirmado:** 26 registros para comparações; auditoria dos 31. **Método principal:** CWT Morlet analítica da correlação enviesada, alternative WPS, alternative WCS, alternative WCO e média da alternative WPS em todos os lags.

## Comece aqui

Veja primeiro a [conclusão executiva honesta](EXECUTIVE_SUMMARY.md), a [verificação das equações com o PDF](EQUATION_AUDIT.md) e o [experimento de pré-processamento](PREPROCESSING.md).

1. [Relatório: equações, fontes, justificativas e lacunas](REPORT.md)
2. [Resultados calculados e comparação com v5](RESULTS.md)
3. [Galeria dos 87 pares](figures/index.html)
4. [Perfis por padrão](figures/pattern_profiles.png)
5. [Diagnóstico de autocorrelação](figures/diagnostics/autocorrelation_estimators.png)
6. [Sensibilidade ao comprimento/posição da janela](figures/diagnostics/window_sensitivity.png)
7. [Controle das definições para múltiplas realizações](figures/diagnostics/ensemble_definition_control.png)

## Panorama, sinal e alternative WPS — exemplos predeterminados

'''
    readme+='\n'.join(f'- [{PATTERN_PT[pattern]}: {name}]({path})' for pattern,name,path in links)
    readme+='''

O padrão dispersed bubble tem somente um registro elegível. Sua alternative WPS é apresentada; não há alternative WCS/alternative WCO intrapadrão independente para esse padrão.

## Saídas atuais

- [EXECUTIVE_SUMMARY.md](EXECUTIVE_SUMMARY.md): conclusão atual, cinco exemplos de padrões e quatro comparações de alternative WCS/alternative WCO.
- [RESULTS.md](RESULTS.md): máximos integrais, perfis e diagnóstico nulo recalculado.
- [PREPROCESSING.md](PREPROCESSING.md): triagem, comparação com v5 e avaliação retirando uma aquisição, sem máscaras.
- [Perfis numéricos](results/profiles.csv), [resultados por registro](results/per_run_methods.csv) e [resultados por par](results/pair_methods.csv): campos integrais e no lag zero; as antigas colunas de COI/sobreposição foram retiradas.
- [Comparação numérica com v5](results/preprocessing_comparison_pairs.csv), [resumos dos candidatos](results/preprocessing_comparison_rankings.csv), [médias dos 87 pares](results/preprocessing_pooled_comparison.csv) e [conclusão estruturada](results/preprocessing_conclusion.json): versão atual com busca em todos os lags.
- [Diagnóstico nulo](results/null_diagnostic.csv) e [metadados](results/null_metadata.json): mesmo domínio completo para observações e sinais aleatórios.
- [Galeria](figures/index.html) e figuras em figures/: nomes dos padrões em inglês e demais textos em PT-BR, sem máscaras.
- [Verificação final](results/final_verification.json) e [manifesto](results/artifact_manifest.json): cobertura, equações e procedência.
- [Tempos do recálculo](results/recalculo_sem_mascaras.json), [migração das tabelas](results/migracao_tabelas_integrais.json) e [preservação dos coeficientes armazenados](results/migracao_cache_sem_mascaras.json): registro desta revisão.

## Organização e reprodução

- scripts/: núcleo matemático, auditoria, execução, controles e verificação final.
- results/: tabelas, parâmetros, testes, logs e manifestos.
- figures/: 87 mapas de pares, 31 comparações individuais e diagnósticos.
- [results/historico_com_mascaras.zip](results/historico_com_mascaras.zip): documentos e saídas anteriores à retirada das máscaras; somente para rastreabilidade.
- literature/: extrações dos PDFs, páginas escaneadas de Schuck e artigo adicional.
- intermediate/: sinais conferidos e matrizes reduzidas apenas para renderização. As estatísticas são calculadas na resolução integral; os mapas armazenados têm aproximadamente 1800 posições de lag. Não use esses mapas reduzidos para recalcular máximos integrais.

Para regenerar somente as figuras e os relatórios a partir dos resultados existentes, sem repetir a triagem ou o diagnóstico nulo:

~~~sh
MPLCONFIGDIR=/tmp/v16-mpl python3 src2/v16_spatial2/scripts/regenerate_figures.py
MPLCONFIGDIR=/tmp/v16-mpl python3 src2/v16_spatial2/scripts/verify_final.py
~~~

Os nomes exibidos dos padrões são dispersed bubble, plug, slug, churn e wavy annular; os demais textos estão em PT-BR. Os identificadores internos das tabelas foram preservados. A [verificação dos rótulos](results/revisao_nomes_padroes.json) registra a regeneração e a preservação dos dados numéricos.

Dependências testadas em requirements.txt. Python 3.13.2 foi usado. Para reproduzir também os cálculos completos, execute a partir da raiz do projeto:

~~~sh
export MPLCONFIGDIR=/tmp/v16-mpl
python3 src2/v16_spatial2/scripts/audit_inputs.py
python3 src2/v16_spatial2/scripts/validate_wavelets.py
python3 src2/v16_spatial2/scripts/run_analysis.py --stage all
python3 src2/v16_spatial2/scripts/verify_legacy_baseline.py
python3 src2/v16_spatial2/scripts/additional_diagnostics.py
python3 src2/v16_spatial2/scripts/precision_and_windows.py
python3 src2/v16_spatial2/scripts/null_diagnostic.py --surrogates 199
python3 src2/v16_spatial2/scripts/preprocessing_experiment.py --scales 64 --tag screen
python3 src2/v16_spatial2/scripts/evaluate_preprocessing.py --tag screen
python3 src2/v16_spatial2/scripts/preprocessing_experiment.py --scales 200 --variants baseline central_band20 core_row_median --tag confirm
python3 src2/v16_spatial2/scripts/evaluate_preprocessing.py --tag confirm
python3 src2/v16_spatial2/scripts/compare_v5.py
python3 src2/v16_spatial2/scripts/report_preprocessing.py
python3 src2/v16_spatial2/scripts/finish_outputs.py
python3 src2/v16_spatial2/scripts/verify_final.py
~~~

Os tempos efetivos desta atualização estão nos logs e metadados em results. O diagnóstico nulo tem duração registrada em null_metadata.json; as etapas de pré-processamento, em preprocessing_*_metadata.json. Tempos históricos não são estimativas de uma nova execução e variam com carga e hardware. A extração/OCR dos PDFs é material de pesquisa preservado, não uma dependência para recalcular as wavelets.

O .gitignore existente ignora src2; os arquivos estão disponíveis no workspace, mas não foram adicionados ao Git. Nenhum commit foi feito.
'''
    (ROOT/'README.md').write_text(readme)
    blocks=['<!doctype html><html lang="pt-BR"><meta charset="utf-8"><title>v16 — pares Schuck</title>',
            '<style>body{font:16px system-ui;max-width:1200px;margin:32px auto;padding:0 16px}img{width:100%;height:auto}article{margin:32px 0}a{color:#125a85}</style>',
            '<h1>v16 — alternative WPS, alternative WCS, alternative WCO e média em lag</h1><p>Morlet de Schuck; correlação enviesada. Plano completo, sem COI nem máscara conservadora. Médias da alternative WPS em todos os lags. Aquisições independentes.</p>',
            '<p><a href="../REPORT.md">Relatório</a> · <a href="../RESULTS.md">Resultados</a></p>']
    for pattern,g in pairs[pairs.method=='v16_biased'].groupby('pattern'):
        blocks.append(f'<h2>{PATTERN_PT[pattern]}</h2>')
        for row in g.itertuples():
            file=f'pairs/{row.run_x}__{row.run_y}.png';label=f'{row.run_x} × {row.run_y}'
            blocks.append(f'<article><h3>{label}</h3><a href="{file}"><img loading="lazy" src="{file}" alt="alternative WPS alternative WCS alternative WCO: {label}"></a></article>')
    blocks.append('</html>');(ROOT/'figures/index.html').write_text('\n'.join(blocks))
    print('Resultados, README, galeria e figuras representativas atualizados: padrões em inglês, demais textos em PT-BR, sem máscaras.')

if __name__=='__main__':main()
