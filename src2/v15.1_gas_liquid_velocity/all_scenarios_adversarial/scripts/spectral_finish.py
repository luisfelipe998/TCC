#!/usr/bin/env python3
"""Independent saved-result checks and spectral-study figures/report."""
from pathlib import Path
import importlib.util
import hashlib
import json
import os
import numpy as np
import pandas as pd

BASE=Path(__file__).resolve().parents[1]
OLD=BASE.parent
os.environ.setdefault('MPLCONFIGDIR','/tmp/v151_spectral_mplcache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

def main():
    p=pd.read_csv(BASE/'results/spectral_pairs.csv')
    r=pd.read_csv(BASE/'results/spectral_regions.csv')
    b=pd.read_csv(BASE/'results/spectral_blocks.csv')
    n=pd.read_csv(BASE/'results/spectral_nulls.csv')
    synthetic=pd.read_csv(BASE/'results/spectral_synthetic.csv')
    scope=pd.read_csv(OLD/'results/all_run_status.csv')
    # Adversarial correction: criteria must coincide in each region/block and
    # all three regions. Preserve the previous marginal-gate result.
    saturation=r.set_index(['run','region']).saturation_fraction
    b['saturation_fraction']=[saturation.loc[(x,y)] for x,y in zip(b.run,b.region)]
    b['joint_supported']=b.transport_supported & b.band_ok & b.axial_ok & (b.saturation_fraction<=.20)
    for i,z in r.iterrows():
        q=b[(b.run==z.run)&(b.region==z.region)]
        if 'marginal_robust_supported' not in r:
            r['marginal_robust_supported']=r.robust_supported
        valid=q[q.joint_supported]
        r.loc[i,'joint_supported_blocks']=len(valid)
        r.loc[i,'joint_accepted_mean_px_s']=valid.signed_px_s.mean()
        r.loc[i,'robust_supported']=len(valid)>=4 and z.same_direction_blocks>=4
        r.loc[i,'joint_reasons']='' if len(valid)>=4 and z.same_direction_blocks>=4 else 'fewer_than_four_joint_blocks_or_mixed_sign'
    for i,z in p.iterrows():
        q=b[b.run==z.run]
        simultaneous=q.groupby('block').joint_supported.all()
        temporal_simultaneous=q.groupby('block').transport_supported.all()
        ids=simultaneous.index[simultaneous].tolist()
        valid=q[q.block.isin(ids)]
        if 'marginal_robust_supported' not in p:
            p['marginal_robust_supported']=p.robust_supported
        if 'marginal_temporal_supported' not in p:
            p['marginal_temporal_supported']=p.temporal_supported
        p.loc[i,'temporal_joint_blocks']=int(temporal_simultaneous.sum())
        p.loc[i,'temporal_supported']=temporal_simultaneous.sum()>=4
        p.loc[i,'joint_supported_blocks']=len(ids)
        p.loc[i,'joint_block_ids']=';'.join(map(str,ids))
        p.loc[i,'accepted_core_signed_px_s']=valid[valid.region=='core'].signed_px_s.mean()
        p.loc[i,'accepted_wall_signed_px_s']=valid[valid.region.isin(['upper_outer','lower_outer'])].signed_px_s.mean()
        p.loc[i,'robust_supported']=len(ids)>=4
        p.loc[i,'full_candidate_mean_supported']=len(ids)==5
        p.loc[i,'joint_reasons']='' if len(ids)>=4 else 'fewer_than_four_simultaneous_joint_blocks'
    b.to_csv(BASE/'results/spectral_blocks.csv',index=False)
    r.to_csv(BASE/'results/spectral_regions.csv',index=False)
    p.to_csv(BASE/'results/spectral_pairs.csv',index=False)
    archive=np.load(OLD/'intermediate/annular_profiles.npz')
    curves=np.load(BASE/'results/spectral_score_curves.npz')
    v=curves['velocities'];max_curve_error=0.;count=0
    assert len(p)==26 and set(p.run)==set(scope.run)
    assert np.isfinite(p[['core_signed_px_s','wall_signed_px_s']]).all().all()
    assert (p[['core_signed_px_s','wall_signed_px_s']]<0).all().all()
    for z in b.itertuples():
        a=archive[f'{z.run}__{z.region}'][z.first_frame:z.last_frame+1].astype(float)
        d=a-a.mean(axis=0)[None,:]
        d-=d.mean(axis=1)[:,None]
        power=abs(np.fft.rfft2(d*np.hanning(len(a))[:,None]*np.hanning(a.shape[1])[None,:]))**2
        power[0]=0;power[1]=0;power[-1]=0
        fx=np.fft.rfftfreq(a.shape[1]);ix=np.flatnonzero((fx>=1/256)&(fx<=1/32))
        totals=power[:,ix].sum(axis=0)
        normalized=power[:,ix]/totals
        weights=totals**.5*(totals>max(totals.max()*1e-10,1e-15))
        phase=(-v[:,None]*fx[ix][None,:])%1*len(a)
        floor=np.floor(phase);frac=phase-floor;indices=floor.astype(int)%len(a)
        line=np.zeros_like(phase)
        for k in (-1,0,1):
            for offset,w in [(0,1-frac),(1,frac)]:
                line+=w*normalized[(indices+k+offset)%len(a),np.arange(len(ix))[None,:]]
        recovered=np.sum(line*weights[None,:],axis=1)/weights.sum()
        saved=curves[f'{z.run}__{z.region}__block{z.block}'].astype(float)
        delta=float(np.max(abs(recovered-saved)));max_curve_error=max(max_curve_error,delta)
        assert delta<5e-8,(z.run,z.region,z.block,delta)
        at=int(np.argmax(recovered));assert v[at]==z.signed_px_frame
        assert abs(recovered[at]-z.score)<1e-12
        excluded=abs(v-z.signed_px_frame)>max(3,.15*abs(z.signed_px_frame))
        margin=recovered[at]-recovered[excluded].max()
        assert abs(margin-z.margin)<1e-12
        nn=n[(n.run==z.run)&(n.region==z.region)&(n.block==z.block)]
        assert len(nn)==19 and abs(nn.score.max()-z.shuffle_max)<1e-12
        support=(z.texture>=1 and abs(z.signed_px_frame)>=1 and not z.boundary
                 and z.score>=.10 and z.score>=z.shuffle_max+.03 and z.margin>=.03
                 and abs(z.signed_px_frame+z.reverse_signed_px_frame)<=.25)
        assert support==z.transport_supported
        assert np.isclose(z.signed_px_s,z.signed_px_frame*z.fps,rtol=1e-13)
        count+=len(recovered)
    for z in r.itertuples():
        q=b[(b.run==z.run)&(b.region==z.region)]
        assert len(q)==5 and sum(q.last_frame-q.first_frame+1)==501
        assert np.isclose(q.signed_px_s.mean(),z.signed_px_s,rtol=1e-13)
        sup=q.transport_supported.sum();bands=q.band_ok.sum();axial=q.axial_ok.sum()
        direction=max((q.signed_px_frame<0).sum(),(q.signed_px_frame>0).sum())
        assert sup==z.supported_blocks and bands==z.band_consistent_blocks and axial==z.axial_consistent_blocks
        assert bool(sup>=4 and direction>=4)==z.temporal_supported
        assert bool(sup>=4 and bands>=4 and axial>=4 and direction>=4 and z.saturation_fraction<=.20)==z.marginal_robust_supported
        joint=q.transport_supported&q.band_ok&q.axial_ok&(q.saturation_fraction<=.20)
        assert joint.sum()==z.joint_supported_blocks
        assert bool(joint.sum()>=4 and direction>=4)==z.robust_supported
    for z in p.itertuples():
        q=r[r.run==z.run].set_index('region')
        assert q.loc['core','signed_px_s']==z.core_signed_px_s
        assert np.isclose(q.loc[['upper_outer','lower_outer'],'signed_px_s'].mean(),z.wall_signed_px_s,rtol=1e-13)
        assert q.marginal_robust_supported.all()==z.marginal_robust_supported
        assert q.temporal_supported.all()==z.marginal_temporal_supported
        qq=b[b.run==z.run]
        joint=qq.groupby('block').joint_supported.all()
        temporal=qq.groupby('block').transport_supported.all()
        assert temporal.sum()==z.temporal_joint_blocks and (temporal.sum()>=4)==z.temporal_supported
        assert joint.sum()==z.joint_supported_blocks
        assert (joint.sum()>=4)==z.robust_supported
        assert (joint.sum()==5)==z.full_candidate_mean_supported
        sel=qq[qq.block.isin(joint.index[joint])]
        if len(sel):
            assert np.isclose(sel[sel.region=='core'].signed_px_s.mean(),z.accepted_core_signed_px_s)
            assert np.isclose(sel[sel.region!='core'].signed_px_s.mean(),z.accepted_wall_signed_px_s)
    checks=dict(status='passed',candidate_runs=len(p),regions=len(r),blocks=len(b),
        saved_null_trials=len(n),independently_recomputed_curve_values=count,
        maximum_numpy_vs_saved_scipy_float32_curve_difference=max_curve_error,
        region_and_pair_aggregations_checked=True,unit_and_scope_checked=True,
        limitation='Independent implementation confirms arithmetic and declared gates; no experimental phase identity or accuracy validation.')
    checks['jointly_supported_pairs']=int(p.robust_supported.sum())
    checks['all_five_block_means_supported']=int(p.full_candidate_mean_supported.sum())
    metadata=json.loads((BASE/'results/spectral_metadata.json').read_text())
    metadata['postprocessing_script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    metadata['marginal_pair_support']=int(p.marginal_robust_supported.sum())
    metadata['final_joint_pair_support']=int(p.robust_supported.sum())
    metadata['full_candidate_mean_support']=int(p.full_candidate_mean_supported.sum())
    metadata['marginal_temporal_pair_support']=int(p.marginal_temporal_supported.sum())
    metadata['temporal_pair_support']=int(p.temporal_supported.sum())
    metadata['final_joint_gate']='At least4 simultaneous blocks satisfying transport, both alternative bands, both axial halves and saturation<=20% in ALL THREE regions. All-five-block mean supported only if all5 simultaneous blocks pass.'
    metadata['accepted_aggregation']='Separate accepted coordinates use only simultaneously passing blocks; never present a failed block as accepted. Candidate means of all5 remain distinct. Initial marginal-gate fields preserved.'
    (BASE/'results/spectral_metadata.json').write_text(json.dumps(metadata,indent=2))
    (BASE/'results/spectral_verification.json').write_text(json.dumps(checks,indent=2))
    # Reuse layout only; numerical verification above imports no estimator.
    spec=importlib.util.spec_from_file_location('old_paper',OLD/'scripts/paper_figures.py')
    layout=importlib.util.module_from_spec(spec);spec.loader.exec_module(layout)
    colors=layout.COLORS;marks=layout.MARKERS;pt=layout.PT
    fig,(ax,key)=plt.subplots(1,2,figsize=(9.3,6.1),gridspec_kw={'width_ratios':[2.4,1]})
    annotations=[]
    for i,z in enumerate(p.itertuples()):
        ax.scatter(z.core_magnitude_px_s,z.wall_magnitude_px_s,marker=marks[z.pattern],s=42,
                   edgecolors=colors[z.pattern],facecolors=colors[z.pattern] if z.full_candidate_mean_supported else 'none',lw=1,zorder=3)
        annotations.append(ax.annotate(str(i+1),(z.core_magnitude_px_s,z.wall_magnitude_px_s),xytext=(7,7),
                                       textcoords='offset points',fontsize=7))
        key.text(.02,1-(i+.8)/27,f'{i+1:02d}  {z.run}',fontsize=8,color=colors[z.pattern],va='center')
    ax.set(xlim=(25000,280000),ylim=(25000,280000),
        xlabel='Propagação óptica no núcleo (px/s)',ylabel='Propagação óptica nas faixas externas (px/s)',
        title='Alternativa espectral: todas as 26 gravações')
    ax.set_aspect('equal');ax.grid(alpha=.18)
    formatter=FuncFormatter(lambda x,pos:f'{x:,.0f}'.replace(',','.'))
    ax.xaxis.set_major_formatter(formatter);ax.yaxis.set_major_formatter(formatter)
    key.axis('off');key.set_title('Identificação')
    handles=[Line2D([],[],linestyle='',marker=marks[x],color=colors[x],label=pt[x]) for x in colors]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.085),ncol=5,frameon=False,fontsize=8)
    fig.text(.5,.015,'Média de cinco velocidades dominantes locais. Preenchido: cinco blocos com suporte conjunto (2/26); vazio: candidato rejeitado.\n'
             f'{int(p.temporal_supported.sum())}/26 passaram pelos critérios temporais simultâneos. Os eixos não representam velocidades materiais de gás e líquido.',
             ha='center',va='bottom',fontsize=8)
    fig.tight_layout(rect=(0,.16,1,1));layout.place_labels(ax,annotations)
    for ext in ('png','pdf','svg'):fig.savefig(BASE/'figures'/f'spectral_all_candidates_PT_BR.{ext}',dpi=600,bbox_inches='tight')
    plt.close(fig)
    original=pd.read_csv(OLD/'results/annular_variant_pairs.csv')
    original=original[original.variant=='outer_bands'].set_index('run')
    comparison=p.set_index('run')[['core_signed_px_s','wall_signed_px_s']].copy()
    comparison['old_core_signed_px_s']=original.core_signed_px_s
    comparison['old_wall_signed_px_s']=original.nearwall_signed_px_s
    comparison['relative_core_difference']=abs(comparison.core_signed_px_s-comparison.old_core_signed_px_s)/abs(comparison.old_core_signed_px_s)
    comparison['relative_wall_difference']=abs(comparison.wall_signed_px_s-comparison.old_wall_signed_px_s)/abs(comparison.old_wall_signed_px_s)
    comparison.to_csv(BASE/'results/spectral_vs_ncc.csv')
    report=f'''# Alternativa espectral adversarial: propagação regional

**Resultado:** foram calculadas 26 coordenadas completas, mas somente **{int(p.temporal_supported.sum())}/26 pares** passaram pelos controles temporais e **{int(p.robust_supported.sum())}/26** permaneceram após os controles de frequência, posição axial e saturação. A alternativa não completou 26 pares sustentados nem identificou velocidades materiais das fases. Os dois pares finais foram m05q01_001 (Plug) e m07q01_001 (bolhas dispersas).

![Todos os candidatos espectrais](../figures/spectral_all_candidates_PT_BR.png)

## Hipótese nova

Em vez do pico da correlação média de toda a gravação, estimou-se uma velocidade dominante em cada um de cinco blocos temporais disjuntos. Para um campo I(x,t)=F(x-vt), a energia de Fourier se concentra em f_t=-v f_x, incluindo a periodicidade temporal decorrente da amostragem. Procurou-se essa direção no espectro de cada bloco. Essa definição permite variação de velocidade entre blocos, sem exigir a antiga concordância de velocidade constante nos intervalos de 1/2/4 quadros.

Usaram-se as mesmas regiões em todas as gravações: núcleo 35:54, faixa superior externa 18:29 e inferior externa 60:71. Os perfis foram derivados das imagens nativas anteriormente e foram reutilizados com hash registrado. Cada região foi descrita pela média de cinco máximos espectrais locais, com peso igual por bloco. A parede foi definida como a média das duas regiões externas: não se exigiu igualdade entre as duas, e sua diferença ficou salva. Essa média descreve duas regiões; não pressupõe um filme simétrico ou velocidade líquida média.

## Regras fixas e controles

Nenhum parâmetro, faixa ou método foi selecionado por gravação, classe ou proximidade ao gráfico desejado. Por isso não houve seleção de hiperparâmetros que precisasse de separação treino/teste. Todos os resultados continuam exploratórios; o conjunto não é uma validação externa.

- Busca de −80 a +80 px/quadro, passo de 0,25; faixa espacial principal 1/256–1/32 ciclo/pixel.
- Remoção da média temporal de cada bloco e da média espacial de cada quadro; janela Hann nos dois eixos; exclusão do zero temporal e seus vizinhos imediatos.
- Score: fração de energia em uma faixa de três bins temporais sobre a direção candidata, normalizada por frequência espacial e ponderada pela raiz da energia dessa frequência.
- Suporte por bloco: desvio dinâmico ≥1 nível de cinza; |v|≥1; máximo fora do limite da busca; score≥0,10 e ≥máximo de 19 embaralhamentos+0,03; margem≥0,03 contra direções fora de max(3 px/quadro,15% de |v|); reversão com erro≤0,25 px/quadro.
- Suporte temporal por região: pelo menos quatro de cinco blocos com suporte e quatro com o mesmo sinal.
- Critérios completos: também pelo menos quatro blocos concordando com **ambas** as faixas alternativas 1/512–1/64 e 1/128–1/16, e com **ambas** as metades axiais da imagem, dentro de max(3 px/quadro,25% de |v|); saturação≤20%.

**Correção adversarial:** a regra inicial contava os critérios separadamente e poderia aceitar conjuntos de blocos diferentes. A regra final exige que transporte, bandas, metades axiais e saturação passem **no mesmo bloco e simultaneamente nas três regiões** em pelo menos quatro de cinco blocos. Os resultados marginais originais permanecem em `marginal_temporal_supported` e `marginal_robust_supported`. A triagem temporal passou de 13 pares marginais para 8 simultâneos. As duas gravações finais passaram nos cinco blocos; portanto, essa correção não alterou a contagem de dois pares com critérios completos neste conjunto.

O CSV preserva separadamente a média candidata dos cinco blocos e a média apenas dos blocos conjuntamente aceitos (`accepted_core_signed_px_s`, `accepted_wall_signed_px_s`). A média dos cinco blocos só recebe suporte quando **todos os cinco** passam; esse é o preenchimento da figura. Médias parciais com menos de quatro blocos são apenas diagnósticos numéricos.

Esses limites são critérios empíricos de diagnóstico, não erro experimental conhecido nem nível de significância. Os 19 embaralhamentos são controles negativos; não se declarou p-valor. Faixas e metades axiais selecionam estruturas diferentes: discordância limita a representatividade de uma única velocidade, sem provar que todos os candidatos estejam errados.

## Resultados dos testes

Foram avaliadas 78 regiões, 390 blocos, **7.410 controles de ordem temporal embaralhada**, três faixas espectrais e duas posições axiais. Dos 390 blocos, **{int(b.transport_supported.sum())}** passaram pelos critérios temporais, **{int(b.band_ok.sum())}** concordaram entre faixas e **{int(b.axial_ok.sum())}** entre metades axiais. As **390 reversões** concordaram com a troca de sinal. Essa última propriedade confirma a convenção matemática, sem provar movimento material.

As translações sintéticas conhecidas −51, −31,8, −17,3, −6,7 e +23,75 px/quadro apresentaram erro máximo **0,05 px/quadro**; duas camadas com texturas independentes a −35 e −14 px/quadro foram recuperadas separadamente. O sinal estacionário, o constante e os quadros aleatórios foram rejeitados como transporte. Um controle com cinco velocidades sucessivas recuperou a média −21 px/quadro.

No controle com movimentos simultâneos a −12 e −38 px/quadro, o algoritmo escolheu −12 px/quadro; isso não equivale à média −25, nem identifica uma fase. Os sinais usaram multiplicadores iguais, mas não tinham energia exatamente igual. A pontuação depende da energia, das janelas e da largura espectral: mesmo uma direção bem definida pode esconder outra velocidade. Esse é um limite essencial do observável dominante.

## Interpretação adversarial

A mudança de estimador produziu uma nuvem visualmente próxima da diagonal para muitas gravações. Isso é compatível com influência compartilhada de interfaces, ondas ou padrões de brilho nas regiões e não demonstra igualdade das velocidades materiais. Em particular, não é legítimo chamar o núcleo de gás e as faixas externas de líquido somente por localização.

Os valores mudaram em relação ao pico da correlação média; ambas as estatísticas também atribuem pesos diferentes às estruturas e ao tempo. A tabela `spectral_vs_ncc.csv` mantém todas essas diferenças, sem escolher o resultado mais conveniente. Exemplo: m04q01 mudou de aproximadamente −34.492 para −88.200 px/s no núcleo. A alternativa passou pelos controles temporais nesse caso, mas falhou nos controles completos. Não há base para promover a nova estimativa como correção física da anterior.

**Conclusão:** o modelo temporalmente variável é uma opção legítima para estudar propagação óptica local. Os resultados executados não justificam substituir as 26 coordenadas por velocidades de gás e líquido aceitas. A cobertura só parece completa quando se mantêm candidatos que falharam nos diagnósticos.

## Verificação e reprodução

`spectral_finish.py` reimplementou o cálculo espectral usando NumPy, sem importar o estimador: **{count:,} valores** de curvas foram conferidos, com diferença máxima **{max_curve_error:.3g}**, compatível com o armazenamento float32 das curvas. Também conferiu critérios dos 390 blocos, máximos dos controles salvos, 78 agregações regionais, 26 pares, unidades e escopo. Não refez as imagens experimentais nem validou a identidade material.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B src2/v15.1_gas_liquid_velocity/all_scenarios_adversarial/scripts/spectral_transport.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B src2/v15.1_gas_liquid_velocity/all_scenarios_adversarial/scripts/spectral_finish.py
```

Dados completos: [pares](../results/spectral_pairs.csv), [regiões](../results/spectral_regions.csv), [blocos](../results/spectral_blocks.csv), [sintéticos](../results/spectral_synthetic.csv), [metadados](../results/spectral_metadata.json), [verificação](../results/spectral_verification.json).
'''
    (BASE/'findings/spectral_review.md').write_text(report)
    print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
