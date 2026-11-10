#!/usr/bin/env python3
"""Audit saved residual decisions and publish PT-BR figures and findings."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parents[1]/'intermediate/residual_mplcache'))
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import residual_motion as study

BASE=study.BASE
METHOD_NAMES={'raw':'Imagem original','fixed_background':'Fundo fixo','interface_distance':'Distância à interface'}
PATTERN_NAMES={'plug':'Pistonado (plug)','slug':'Golfadas (slug)','churn':'Agitado (churn)',
    'wavy annular':'Anular ondulado','dispersed bubbles':'Bolhas dispersas'}
COLORS={'plug':'#0072B2','slug':'#009E73','churn':'#D55E00','wavy annular':'#CC79A7','dispersed bubbles':'#202020'}

def savefig(fig,name):
    for suffix in ('png','pdf','svg'):
        fig.savefig(BASE/'figures'/f'{name}.{suffix}',dpi=300,bbox_inches='tight')
    plt.close(fig)

def run():
    metadata=json.loads((BASE/'results/residual_metadata.json').read_text())
    assert metadata['status']=='complete' and metadata['runs']==26
    rows=pd.read_csv(BASE/'results/residual_summary.csv')
    blocks=pd.read_csv(BASE/'results/residual_blocks.csv')
    lags=pd.read_csv(BASE/'results/residual_lags.csv')
    nulls=pd.read_csv(BASE/'results/residual_nulls.csv')
    pairs=pd.read_csv(BASE/'results/residual_pairs.csv')
    inputs=pd.read_csv(BASE/'results/residual_inputs.csv')
    preds=pd.read_csv(BASE/'results/residual_predictor.csv')
    closed=pd.read_csv(BASE/'results/residual_closed_m07q01_summary.csv')
    closed_pairs=pd.read_csv(BASE/'results/residual_closed_m07q01_pairs.csv')
    assert len(closed)==12 and len(closed_pairs)==6 and not closed_pairs.pair_optical_supported.any()
    closed_fixed=closed[(closed.phase=='liquid_candidate')&(closed.method=='fixed_background')].set_index('erosion_radius')
    closed_dynamic=closed[(closed.phase=='liquid_candidate')&(closed.method=='interface_distance')].set_index('erosion_radius')
    assert (closed[closed.phase=='gas_candidate'].accepted_blocks==0).all()
    assert (closed_fixed.accepted_blocks==6).all() and (closed_dynamic.accepted_blocks==2).all()
    assert len(rows)==312 and len(blocks)==1872 and len(lags)==len(nulls)==5616 and len(pairs)==156
    assert not rows.material_phase_identified.any() and not pairs.material_pair_identified.any()
    assert metadata['script_sha256']==hashlib.sha256(Path(study.__file__).read_bytes()).hexdigest()
    keys=['run','erosion_radius','method','phase','block']
    lag_index=lags.set_index(keys+['lag']).sort_index()
    null_index=nulls.groupby(keys).peak.max()
    decision_errors=[]
    for _,b in blocks.iterrows():
        key=tuple(b[k] for k in keys);lagdata=lag_index.loc[key]
        v=lagdata.loc[1,'shift'];tol=max(2,.2*abs(v))
        okay=all(np.isfinite(p.peak) and p.peak>=.2 and p['std']>=1 and not p.boundary and abs(p['shift']/lag-v)<=tol
            for lag,p in lagdata.iterrows())
        peak=lagdata.loc[1];null=null_index.loc[key]
        okay=bool(okay and b.frames_with_100_pixels>=20 and b.raw_phase_std>=1 and
            np.isfinite(v) and np.isfinite(peak.margin) and peak.margin>=.02 and
            np.isfinite(null) and peak.peak>=null+.05 and b.saturation_fraction<=.2)
        if okay!=b.optical_supported:decision_errors.append(key)
    assert not decision_errors,decision_errors
    pair_errors=[]
    for _,p in pairs.iterrows():
        bs=blocks[(blocks.run==p.run)&(blocks.erosion_radius==p.erosion_radius)&(blocks.method==p.method)]
        accepted=bs.groupby('block').optical_supported.all()
        common=accepted[accepted].index.tolist()
        assert p.common_blocks==len(common) and p.pair_optical_supported==(len(common)>=5)
        assert p.matched_intervals==sum(bs[(bs.phase=='gas_candidate')&(bs.block.isin(common))].frame_intervals)
        for phase,column in [('gas_candidate','gas_candidate_px_frame'),('liquid_candidate','liquid_candidate_px_frame')]:
            bb=bs[(bs.phase==phase)&bs.block.isin(common)]
            expected=np.average(bb.signed_px_frame,weights=bb.frame_intervals) if len(bb) else np.nan
            assert np.isclose(p[column],expected,atol=1e-10,equal_nan=True)
        if p.pair_optical_supported and abs(p.liquid_candidate_px_frame)>2:
            assert np.isclose(p.signed_ratio_gas_candidate_over_liquid_candidate,p.gas_candidate_px_frame/p.liquid_candidate_px_frame)
    mask_metadata=[]
    for _,i in inputs.iterrows():
        file=study.MASKS/(i.run+'.npz')
        assert i.mask_sha256==hashlib.sha256(file.read_bytes()).hexdigest()
        with np.load(file,allow_pickle=False) as z:
            assert np.array_equal(z['frame_ids'],np.arange(501))
            assert np.array_equal(z['mask_shape'],[501,80,1016])
        mask_metadata.append(i.run)
    counts=pairs.groupby(['method','erosion_radius'],sort=False).agg(
        supported_pairs=('pair_optical_supported','sum'),
        records=('run','size'),median_common_blocks=('common_blocks','median')).reset_index()
    counts.to_csv(BASE/'results/residual_counts.csv',index=False)
    sensitivity=[]
    for (method,radius),group in pairs.groupby(['method','erosion_radius']):
        for threshold in range(1,7):
            accepted=group[group.common_blocks>=threshold]
            sensitivity.append(dict(method=method,erosion_radius=radius,minimum_common_blocks=threshold,
                records=len(accepted),minimum_retained_intervals=int(accepted.matched_intervals.min()) if len(accepted) else None,
                promoted_primary_result=(threshold==5)))
    pd.DataFrame(sensitivity).to_csv(BASE/'results/residual_support_sensitivity.csv',index=False)
    stable=[]
    for run,bs in blocks.groupby('run'):
        common=bs.groupby('block').optical_supported.all()
        ids=common[common].index.tolist()
        estimates=[]
        for (method,radius,phase),b in bs[bs.block.isin(ids)].groupby(['method','erosion_radius','phase']):
            estimates.append(dict(method=method,erosion_radius=radius,phase=phase,
                velocity=float(np.average(b.signed_px_frame,weights=b.frame_intervals))))
        agree=True if len(ids)>=5 else False
        for phase in study.PHASES:
            vv=[e['velocity'] for e in estimates if e['phase']==phase]
            if len(vv)!=6:agree=False
            elif max(vv)-min(vv)>max(2,.2*abs(np.mean(vv))):agree=False
        active=bs[bs.method!='raw']
        active_common=active.groupby('block').optical_supported.all()
        active_ids=active_common[active_common].index.tolist()
        active_agree=len(active_ids)>=5
        for phase in study.PHASES:
            vv=[]
            for _,b in active[(active.phase==phase)&active.block.isin(active_ids)].groupby(['method','erosion_radius']):
                vv.append(float(np.average(b.signed_px_frame,weights=b.frame_intervals)))
            if len(vv)!=4 or (len(vv)==4 and max(vv)-min(vv)>max(2,.2*abs(np.mean(vv)))):active_agree=False
        stable.append(dict(run=run,common_blocks_all_six_methods_radii=len(ids),
            jointly_supported_all_six=len(ids)>=5,velocities_agree_all_six=agree,
            common_blocks_four_background_options=len(active_ids),
            jointly_supported_four_background_options=len(active_ids)>=5,velocities_agree_four_background_options=active_agree))
    pd.DataFrame(stable).to_csv(BASE/'results/residual_all_options_stability.csv',index=False)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
        'savefig.facecolor':'white','font.family':'DejaVu Sans'})
    ordering=pairs[['run','pattern']].drop_duplicates().sort_values(['pattern','run'])
    columns=[(m,r) for m in study.METHODS for r in (3,1)]
    matrix=np.array([[int(pairs[(pairs.run==run)&(pairs.method==m)&(pairs.erosion_radius==r)].common_blocks.iloc[0])
        for m,r in columns] for run in ordering.run])
    fig,ax=plt.subplots(figsize=(8.5,9.2))
    im=ax.imshow(matrix,cmap='YlGnBu',vmin=0,vmax=6,aspect='auto')
    for i in range(len(ordering)):
        for j in range(len(columns)):
            ax.text(j,i,str(matrix[i,j]),ha='center',va='center',color='white' if matrix[i,j]>=4 else 'black',fontsize=9)
    short={'raw':'Imagem\noriginal','fixed_background':'Fundo\nfixo','interface_distance':'Distância à\ninterface'}
    ax.set_xticks(range(len(columns)),[f'{short[m]}\nr = {r} px' for m,r in columns],fontsize=9)
    ax.set_yticks(range(len(ordering)),ordering.run)
    for label,pattern in zip(ax.get_yticklabels(),ordering.pattern):label.set_color(COLORS[pattern])
    ax.set_title('Blocos com suporte simultâneo nas duas máscaras candidatas',pad=14)
    cbar=fig.colorbar(im,ax=ax,pad=.02,shrink=.6,ticks=range(7));cbar.set_label('Blocos aprovados (de 6)')
    handles=[Line2D([0],[0],marker='o',color=c,ls='',label=PATTERN_NAMES[p]) for p,c in COLORS.items()]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.014),ncol=3,frameon=False,fontsize=8)
    fig.text(.5,.076,'Suporte: ≥5/6 blocos. Nenhum par é uma velocidade material identificada.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.09,1,1));savefig(fig,'residual_support_all_runs_PT_BR')
    fig,axes=plt.subplots(1,3,figsize=(13.5,4.8),sharex=True,sharey=True)
    supported=pairs[pairs.pair_optical_supported]
    values=np.r_[supported.gas_candidate_px_s.to_numpy(),supported.liquid_candidate_px_s.to_numpy()]/1000
    lower=min(-1,float(values.min())*1.08) if len(values) else -1
    upper=max(1,float(values.max())*1.08) if len(values) else 1
    for ax,method in zip(axes,study.METHODS):
        for _,p in supported[supported.method==method].iterrows():
            ax.scatter(p.gas_candidate_px_s/1000,p.liquid_candidate_px_s/1000,s=50,
                edgecolor=COLORS[p.pattern],facecolor=COLORS[p.pattern] if p.erosion_radius==3 else 'none',
                marker='o' if p.erosion_radius==3 else 's')
            offset=(4,6) if p.erosion_radius==1 else (4,-11)
            align='left'
            if p.run in ('m05q10_001','m07q20_001'):
                offset=(-6,7) if p.erosion_radius==1 else (-6,-11)
                align='right'
            ax.annotate(p.run.replace('_001',''),(p.gas_candidate_px_s/1000,p.liquid_candidate_px_s/1000),
                xytext=offset,textcoords='offset points',fontsize=7,ha=align)
        ax.plot([lower,upper],[lower,upper],color='.7',lw=.7,ls='--')
        ax.axhline(0,color='.85',lw=.5);ax.axvline(0,color='.85',lw=.5)
        if not len(supported[supported.method==method]):
            ax.text(.5,.5,'Nenhum par com ≥5/6 blocos',transform=ax.transAxes,ha='center',fontsize=9)
        ax.set_title(METHOD_NAMES[method]);ax.set_xlabel('Máscara candidata de gás\nMovimento óptico (10³ px/s)')
        ax.set_xlim(lower,upper);ax.set_ylim(lower,upper);ax.grid(alpha=.15)
    axes[0].set_ylabel('Máscara candidata de líquido\nMovimento óptico (10³ px/s)')
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,-.005),ncol=5,frameon=False,fontsize=8)
    fig.text(.5,.055,'Círculo preenchido: erosão de 3 px; quadrado vazio: 1 px. Apenas pares com ≥5/6 blocos comuns.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.08,1,1));savefig(fig,'residual_supported_pairs_PT_BR')
    synth=pd.read_csv(BASE/'results/residual_synthetic.csv')
    challenge=json.loads((BASE/'results/residual_challenge.json').read_text())
    realcheck=json.loads((BASE/'results/residual_real_numerical_check.json').read_text())
    assert realcheck['status']=='passed' and challenge['status']=='passed'
    grouped=counts.set_index(['method','erosion_radius']).supported_pairs
    count_text='\n'.join(f'| {METHOD_NAMES[m]} | {int(grouped[m,3])}/26 | {int(grouped[m,1])}/26 |' for m in study.METHODS)
    pair_text='\n'.join(f'| {p.run} | {METHOD_NAMES[p.method]} | {int(p.erosion_radius)} | {p.gas_candidate_px_s:.1f} | {p.liquid_candidate_px_s:.1f} | {p.signed_ratio_gas_candidate_over_liquid_candidate:.3f} | {int(p.matched_intervals)}/500 |' for p in supported.itertuples())
    if not pair_text:pair_text='| Nenhum par aprovado | — | — | — | — | — | — |'
    pure=synth[synth.case.isin(['interface_halo_without_material_texture','deforming_interface_without_material_texture'])]
    pure_counts=pure.groupby('method').supported.sum().to_dict()
    dynamic_success=int(grouped['interface_distance',3]);raw_success=int(grouped['raw',3]);fixed_success=int(grouped['fixed_background',3])
    review=f'''# Resultado da subtração condicionada à interface

Foram processados os **26 registros, 13.026 imagens nativas**, com 312 combinações de registro, fase candidata, tratamento e erosão. O resultado é uma avaliação de movimento óptico condicionado a regiões, não uma identificação de velocidades materiais. **Nenhum par material líquido/gás foi identificado por esta análise.**

## Resultado observado

| Tratamento | Pares com erosão de 3 px | Pares com erosão de 1 px |
|---|---:|---:|
{count_text}

Cada par precisou ter suporte simultâneo nas duas máscaras em pelo menos cinco dos seis blocos. Foram avaliados 1.872 blocos de fase, 5.616 estimativas por atraso e 5.616 controles embaralhados. A subtração por distância à interface produziu {dynamic_success} pares na referência; a imagem original produziu {raw_success}, e o fundo fixo treinado na outra metade produziu {fixed_success}. Esses números não medem acurácia física.

![Suporte em todos os registros](../figures/residual_support_all_runs_PT_BR.png)

O mapa mostra todos os registros, incluindo ausências. Uma célula com poucos blocos pode indicar textura escassa, ambiguidade, saturação, mudança de forma, inconsistência entre atrasos ou falta de separação do controle embaralhado. A ausência não foi substituída por zero nem por uma velocidade obtida de outra região.

![Pares ópticos com suporte](../figures/residual_supported_pairs_PT_BR.png)

| Registro | Tratamento | Erosão (px) | Gás candidato (px/s) | Líquido candidato (px/s) | Razão com sinal | Intervalos usados |
|---|---|---:|---:|---:|---:|---:|
{pair_text}

Os valores com sinal correspondem ao eixo horizontal da imagem: negativo significa deslocamento para colunas menores. A razão divide os dois movimentos ópticos candidatos, calculados nos mesmos blocos, e só foi emitida quando o denominador teve magnitude superior a 2 px/quadro. Não é razão de velocidades materiais nem escorregamento das fases. Os arquivos CSV preservam também resultados sem suporte; esses valores não foram promovidos ao gráfico.

A linha tracejada do gráfico indica apenas igualdade numérica dos dois movimentos ópticos. Não é uma relação física ajustada aos dados.

## O que os controles demonstraram

As texturas sintéticas independentes foram recuperadas pelos três tratamentos. No controle adicional, a interface deslocou-se a +12 px/quadro enquanto as texturas interiores prescritas se deslocaram a −17,3 e −7,4 px/quadro. A subtração dinâmica preservou as duas velocidades: o maior erro absoluto, incluindo a inversão temporal, foi **{challenge['maximum_dynamic_interface_texture_error_px_frame']:.6f} px/quadro**. Esse erro pertence somente ao sinal sintético prescrito.

Nos quatro casos fase/halo sem textura material — interface móvel e interface deformável, gás e líquido candidatos — a imagem original aprovou {pure_counts.get('raw',0)}/4, o fundo fixo {pure_counts.get('fixed_background',0)}/4 e o resíduo por distância {pure_counts.get('interface_distance',0)}/4. O modelo conseguiu remover esse mecanismo dentro da família sintética testada. Esses halos foram deliberadamente construídos com uma lei de intensidade em função da distância compatível com o preditor; não representam toda óptica possível do tubo.

**A iluminação móvel compartilhada passou pelos critérios nas duas regiões e nos três tratamentos**, com velocidade aparente de aproximadamente +12 px/quadro, apesar de não representar transporte material. A subtração dinâmica não resolveu esse contraexemplo. Portanto, separar a interface e obter duas correlações consistentes ainda não certifica duas velocidades de fase.

Interiores constantes foram rejeitados. Alterar somente as intensidades da metade de avaliação não modificou suas previsões além de arredondamento de ponto flutuante; o maior desvio observado foi {max(challenge['prediction_invariance_error'].values()):.3g} unidades de intensidade. O cálculo FFT foi comparado a Pearson calculado diretamente em cinco deslocamentos; diferença máxima {challenge['direct_pearson_max_error']:.3g}. A auditoria reconstruiu todas as 1.872 decisões de blocos e os 156 pares finais, verificou os pesos temporais, a razão e os metadados das 26 máscaras.

Doze combinações adicionais de imagem estacionária ou quadros aleatórios independentes com máscara móvel foram rejeitadas. Além dos controles sintéticos, {realcheck['checks']} picos de correlação salvos foram recalculados diretamente a partir das imagens nativas de três registros, com os três tratamentos, duas máscaras candidatas quando definidas e os três atrasos. A maior diferença de Pearson foi {realcheck['maximum_pearson_absolute_error']:.3g}, e a maior diferença de desvio-padrão foi {realcheck['maximum_texture_std_absolute_error']:.3g} unidades de intensidade. Esses resultados verificam a implementação e o alinhamento dos dados; não fornecem verdade experimental para velocidade.

## Limitações e decisão

O tratamento por distância é útil como teste adversarial de contaminação por interface: ele removeu halos especificados e preservou texturas independentes no controle sintético. O ganho observado não autoriza chamar o resíduo de líquido ou gás material. Refração, iluminação e deformação podem permanecer; o modelo também pode remover textura real correlacionada com a geometria. As classes pouco povoadas exigiram interpolação; a maior fração de pixels de avaliação nessa condição foi {preds.extrapolated_pixel_fraction.max():.2%}, registrada por metade no inventário do preditor.

A inspeção de quadros reais mostrou uma limitação adicional concreta. No [quadro 200 de m07q01_001](../figures/residual_example_m07q01_001.png), a máscara herdada contém recortes retangulares sem bolha fechada correspondente no fundo visualmente uniforme. O preditor atribui intensidade a esses recortes; a subtração introduz suas bordas no resíduo. Portanto, um erro de segmentação pode gerar textura residual artificial. Este resultado não foi ocultado nem interpretado como traçador. No [quadro 200 de m04q80_001](../figures/residual_example_m04q80_001.png), a máscara líquida de raio 3 está ausente, embora haja bastante variação óptica: o processamento não cria uma segunda fase observável a partir dessa variação. No [quadro 200 de m05q10_001](../figures/residual_example_m05q10_001.png), boa parte da máscara candidata de líquido consiste em regiões pequenas junto às interfaces. São três inspeções ilustrativas, não uma anotação de verdade de referência dos 13.026 quadros.

O treinamento em outra metade foi aplicado somente aos preditores de intensidade desta investigação. As máscaras e sua referência radial foram herdadas da análise anterior; a separação não representa validação externa da segmentação ou de toda a cadeia de processamento.

Como sensibilidade direcionada ao problema de topologia, m07q01_001 foi processado novamente com as máscaras de contornos fechados já calculadas na investigação anterior, mantendo todos os critérios. Essa execução adicional de 12 combinações não foi misturada aos resultados principais dos 26 registros. Nenhuma combinação forneceu um par aprovado. A máscara candidata de gás teve zero blocos aprovados nas seis opções; a de líquido teve seis blocos com fundo fixo e somente dois após subtração por distância, nos dois raios de erosão. As velocidades ópticas com fundo fixo foram {closed_fixed.loc[3,'signed_px_s']:.1f} px/s (raio 3) e {closed_fixed.loc[1,'signed_px_s']:.1f} px/s (raio 1), mas não há um traçador líquido identificado que autorize essa interpretação material. O ajuste topológico evita determinadas máscaras abertas; ele não transforma a textura restante em transporte líquido. Os resultados completos estão na [sensibilidade de contornos fechados](../results/residual_closed_m07q01_summary.csv).

Os três tratamentos não foram combinados para escolher a alternativa que completasse mais pontos. A tabela de [sensibilidade do suporte](../results/residual_support_sensitivity.csv) mostra os efeitos de exigir de um a seis blocos, mas a regra principal permaneceu cinco. A [estabilidade conjunta](../results/residual_all_options_stability.csv) exige os mesmos blocos em todos os tratamentos e erosões, para evitar comparar janelas diferentes. Ela também apresenta separadamente as quatro opções de subtração de fundo, sem exigir que a imagem original seja aprovada: {sum(s['velocities_agree_four_background_options'] for s in stable)} registros tiveram cinco blocos comuns e concordância dos dois movimentos nesse conjunto de quatro opções. A amplitude entre estimativas precisou ficar abaixo de max(2 px/quadro, 20% da magnitude média), para cada fase candidata.

A comparação numérica com a v15.1 não é direta: aqui houve treinamento na outra metade, seis blocos, decisões por bloco e média ponderada de velocidades de blocos comuns. Na v15.1, outras regras e agregações foram usadas. Uma redução de suporte não prova pior rastreamento, assim como um aumento não prova maior acurácia.

Permanece válida a distinção entre **movimento de textura candidata**, **velocidade de um objeto identificável** e **velocidade média material da fase**. Esta análise fornece o primeiro; seu resultado isolado não sustenta o terceiro, nem mesmo em forma relativa intra ou inter padrão. Razões dentro de um registro cancelam somente uma escala espacial comum. Entre registros, ainda faltam confirmação de magnificação comum, mesma população de traçadores e representatividade temporal. Os 501 quadros cobrem somente 71–125 ms.

## Reprodução e arquivos

Executar, nesta ordem: `python3 scripts/residual_motion.py`, `python3 scripts/residual_challenge.py`, `python3 scripts/residual_closed_contour.py`, `python3 scripts/residual_visual_check.py`, `python3 scripts/residual_real_check.py` e `python3 scripts/residual_report.py`, a partir da pasta desta investigação. Os scripts também funcionam chamados por seus caminhos completos.

- [Protocolo](residual_protocol_PT_BR.md).
- [Pares por registro](../results/residual_pairs.csv), [blocos](../results/residual_blocks.csv), [atrasos](../results/residual_lags.csv), [controles](../results/residual_nulls.csv).
- [Resultados sintéticos](../results/residual_synthetic.csv) e [desafio independente](../results/residual_challenge.csv).
- [Recomputação direta com imagens reais](../results/residual_real_numerical_check.csv).
- [Fontes e hashes](../results/residual_inputs.csv), [metadados](../results/residual_metadata.json), [auditoria](../results/residual_audit.json).
- Figuras em PNG, PDF e SVG, com texto em PT-BR.
'''
    (BASE/'findings/residual_review_PT_BR.md').write_text(review)
    study.dump('residual_audit.json',dict(status='passed',records=26,phase_method_radius_combinations=312,
        reconstructed_block_decisions=1872,reconstructed_pairs=156,mask_metadata_checked=len(mask_metadata),
        supported_pairs=counts.to_dict('records'),all_six_option_stable_records=sum(s['velocities_agree_all_six'] for s in stable),
        four_background_option_stable_records=sum(s['velocities_agree_four_background_options'] for s in stable),
        report_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        material_pairs_identified=0))
    print(counts.to_string(index=False));print('audit passed')

if __name__=='__main__':run()
