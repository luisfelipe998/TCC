"""Publish canonical bilingual summaries from the final numerical tables."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
BASE=Path(__file__).resolve().parents[1]
P={'dispersed bubbles':'Bolhas dispersas','plug':'Plug','slug':'Slug','churn':'Churn','wavy annular':'Anular ondulado'}

def num(x,d=1,pt=False):
    if not np.isfinite(x):return '—'
    s=f'{x:.{d}f}';return s.replace('.',',') if pt else s

def run():
    per=pd.read_csv(BASE/'results/object_relative_per_run.csv')
    pairs=pd.read_csv(BASE/'results/paired_object_consensus.csv')
    counts=pd.read_csv(BASE/'results/residual_counts.csv')
    fps=pd.read_csv(BASE.parent/'v15.1_gas_liquid_velocity/results/scenario_inventory.csv');fps=fps[fps.included]
    status=per.merge(pairs[['run','baseline_supported','threshold_stable','ratio','common_intervals']],on='run')
    status['identified_gas_material_velocity_px_s']=np.nan;status['identified_liquid_material_velocity_px_s']=np.nan
    status['identified_material_velocity_ratio']=np.nan
    status['conclusion']='object/interface motion only; material phase pair not identified'
    status.to_csv(BASE/'results/final_per_run_status.csv',index=False)
    metrics=dict(recordings=len(per),native_frames=int(fps.frames.sum()),small_reference=int(per.small_supported.sum()),
        small_stable=int(per.small_threshold_stable.sum()),meniscus_reference=int(per.meniscus_supported.sum()),
        meniscus_stable=int(per.meniscus_threshold_stable.sum()),common_object_pairs=int(pairs.threshold_stable.sum()),
        identified_material_phase_pairs=0,residual_stable_all_background_options=0)
    (BASE/'results/report_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    for pt in (True,False):
        lines=['# EXECUTIVE SUMMARY — v15.2'+(' (PT-BR)' if pt else ''),'']
        if pt:
            lines += [
                '**Conclusão: a investigação recuperou movimentos de objetos e interfaces, mas não identificou dois movimentos materiais independentes de líquido e gás. O objetivo de obter velocidades líquido/gás, mesmo relativas para todos os padrões, não foi alcançado com os dados e métodos testados. Nenhuma velocidade ausente foi preenchida.**','',
                f'Foram analisados {metrics["recordings"]} cenários e {metrics["native_frames"]:,} imagens nativas, preservando as exclusões solicitadas. As unidades principais são px/s; valores negativos indicam movimento para colunas menores. Os gráficos de velocidade usam magnitudes, pois os resultados selecionados se deslocam nesse sentido. O FPS é documentado; a escala axial em metros por pixel e a magnificação relativa entre gravações permanecem sem confirmação independente.','',
                '## Resultados principais','',
                '| Resultado | Cenários | Significado |','|---|---:|---|',
                f'| Objetos pequenos com movimento resolvido | {metrics["small_reference"]}/26 na referência; **{metrics["small_stable"]}/26** estáveis nos três limiares | Translação de componentes segmentados; não são traçadores líquidos validados. |',
                f'| Envelopes alongados, ajuste robusto das extremidades | {metrics["meniscus_reference"]}/26 na referência; **{metrics["meniscus_stable"]}/26** estáveis | Celeridade de estruturas candidatas a gás. |',
                f'| Comparação de objetos pequenos e envelopes nos mesmos intervalos | **{metrics["common_object_pairs"]}/26** estáveis | Razão entre movimentos de duas populações de objetos, não escorregamento de fases. |',
                '| Textura após subtração condicionada à interface | 2/26 com erosão de 3 px; 4/26 com erosão de 1 px | Nenhum cenário manteve cinco blocos comuns nas quatro opções de fundo/erosão. |',
                '| Pares materiais líquido/gás identificados | **0/26** | Não há gráfico completo validado com essas duas grandezas. |','',
                'Nesta análise, **objeto pequeno significa um componente pequeno da máscara**, não necessariamente uma bolha inteira. A inspeção dos exemplos nativos mostrou componentes que cobrem somente parte do contorno de uma bolha maior. Desabilitar o preenchimento de envelopes abertos não garante um contorno físico fechado. Portanto, até a comparação pequena/envelope pode compartilhar informação da interface; não representa dois traçadores materiais independentes.','',
                '## Gráficos e resultados auditáveis','',
                '![Cobertura completa dos 26 cenários](../figures/coverage_all_scenarios_PT_BR.png)','',
                'Todos os cenários aparecem nessa figura de cobertura. Traço indica ausência de suporte, não velocidade zero. Símbolos preenchidos indicam suporte e concordância nos limiares 0,60/0,70/0,80; os abertos têm suporte somente na referência. As barras representam dispersão entre trajetórias, não incerteza experimental. Não se misturou textura anular com velocidade de bolhas para completar pontos.','',
                '| Cenário | Padrão | Objetos pequenos (10³ px/s) | Envelopes alongados (10³ px/s) | Razão pequena/envelope nos intervalos comuns |','|---|---|---:|---:|---:|']
        else:
            lines += [
                '**Conclusion: the investigation recovered object and interface motion, but did not identify two independent liquid/gas material velocities. The requested liquid/gas velocity result, including relative values across all patterns, was not achieved with the available data and tested methods. Missing velocities were not filled.**','',
                f'The study covered {metrics["recordings"]} eligible recordings and {metrics["native_frames"]:,} native images, preserving the requested exclusions. Primary units are px/s. Negative velocities move toward decreasing image columns; speed figures show magnitudes of those selected motions. Native FPS is documented. Axial metres per pixel and relative magnification across recordings remain independently unverified.','',
                '## Main results','',
                '| Result | Recordings | Meaning |','|---|---:|---|',
                f'| Motion-resolved small objects | {metrics["small_reference"]}/26 at reference; **{metrics["small_stable"]}/26** stable across thresholds | Translation of segmented components, not validated liquid tracers. |',
                f'| Robust elongated-envelope endpoints | {metrics["meniscus_reference"]}/26 at reference; **{metrics["meniscus_stable"]}/26** stable | Celerity of gas-candidate structures. |',
                f'| Small objects and elongated envelopes on common intervals | **{metrics["common_object_pairs"]}/26** stable | Two object-population motions; not phase slip. |',
                '| Interface-conditioned residual texture | 2/26 at erosion 3 px; 4/26 at erosion 1 px | No recording retained five common blocks across all four background/erosion options. |',
                '| Identified material liquid/gas pairs | **0/26** | A validated complete chart of those quantities is unavailable. |','',
                'Here, **small object means a small mask component**, not necessarily an entire bubble. Native examples include components covering only part of a larger bubble outline. Disabling open-envelope filling does not guarantee a physically closed outline. Consequently, small/envelope comparisons may still share interface information rather than follow two independent material tracers.','',
                '## Figures and auditable results','',
                '![All-recording coverage; Portuguese publication figure](../figures/coverage_all_scenarios_PT_BR.png)','',
                'All recordings appear in the coverage figure. A dash means unsupported, not zero speed. Filled markers pass support and agreement at thresholds 0.60/0.70/0.80; open markers pass only the reference. Bars show between-track spread, not experimental uncertainty. Annular texture was not mixed with bubble translation to complete coordinates.','',
                '| Recording | Pattern | Small objects (10³ px/s) | Elongated envelopes (10³ px/s) | Small/envelope ratio on common intervals |','|---|---|---:|---:|---:|']
        for r in status.itertuples():
            sv=abs(r.small_px_s)/1000 if r.small_threshold_stable else np.nan
            mv=abs(r.meniscus_px_s)/1000 if r.meniscus_threshold_stable else np.nan
            lines.append(f'| {r.run} | {P[r.pattern] if pt else r.pattern} | {num(sv,2,pt)} | {num(mv,2,pt)} | {num(r.ratio,3,pt)} |')
        if pt:
            lines += ['',
                'A tabela principal mostra apenas resultados estáveis nos três limiares. As medianas de objetos/envelopes usam suas próprias trajetórias selecionadas. A razão usa outra agregação: médias das velocidades ajustadas que estão presentes nos mesmos intervalos temporais. Por isso, não deve ser calculada dividindo as duas medianas dessa tabela.','',
                '![Comparação de objetos em intervalos comuns](../figures/paired_objects_PT_BR.png)','',
                'Em **m05q05_001 (slug)**, a razão foi **0,9345**, com sensibilidade de **0,9345–0,9924** e 120/500 intervalos na referência. Em **m07q05_001 (plug)**, foi **0,9639**, sensibilidade **0,9444–0,9817**, em 133/500 intervalos. São envelopes de sensibilidade à segmentação, não intervalos de confiança. Os ajustes usam as trajetórias inteiras; apenas sua agregação foi restrita aos intervalos comuns. Proximidade de 1 não demonstra ausência de escorregamento líquido/gás: os pequenos candidatos também podem ser bolhas.','',
                '![Razões entre gravações](../figures/inter_recording_ratios_PT_BR.png)','',
                'A mesma definição de movimento de objetos pequenos permitiu dez comparações entre cinco gravações: três dentro do grupo plug e sete entre gravações de padrões diferentes. São comparações descritivas de populações selecionadas de objetos. Não são velocidades líquidas, efeitos causais do padrão nem médias de toda a fase gasosa. Interpretá-las como razões físicas requer magnificação axial comum; a gravação única de bolhas dispersas não permite estimar variabilidade entre repetições desse padrão.','',
                '## O que foi testado e qual opção funcionou melhor','',
                'O detector de contornos fechados recuperou objetos pequenos, mas rejeitou as trajetórias grandes de plug/slug. Sombras conectadas à parede comprometiam sua geometria. A alternativa de rastrear intervalos axiais do núcleo recuperou meniscos, mas os saltos de localização reprovavam o ajuste original. O ajuste robusto Theil–Sen das duas extremidades, com controles de deformação e concordância temporal, recuperou os três resultados estáveis apresentados. É a opção mais útil aqui para **translação de envelopes alongados**, dentro desse subconjunto; não foi validada como média material de gás.','',
                'O filtro adicional de movimento removeu 29 candidatos de referência quase imóveis ou com percurso ajustado insuficiente. A inspeção encontrou exemplos persistentes perto da coluna 370 em diferentes vídeos, compatíveis com estrutura óptica fixa. Foram preservados nas tabelas originais e excluídos da comparação de objetos em movimento. Mudanças de população, limiar e ponderação permanecem fontes de viés.','',
                'A subtração da contribuição prevista pela distância à interface removeu halos sintéticos especificados e preservou duas texturas independentes. Porém, uma iluminação móvel compartilhada passou pelos critérios nas duas máscaras mesmo sem transporte material. Em imagens reais, erros retangulares da máscara também foram transferidos ao resíduo. O novo método deu o mesmo número de pares apoiados que o fundo fixo (2/26 e 4/26), com conjuntos diferentes; não demonstrou identificação adicional de líquido. Nenhum dos nove anulares ondulados teve par residual apoiado na regra principal.','',
                'Também foram testadas extrapolações de velocidade para tamanho aparente zero: constante, linear, quadrática e quadrática com posição radial. A validação retirou grupos temporais e eliminou do treinamento trajetórias que compartilhassem intervalos com o grupo de teste. A associação tamanho/velocidade variou entre cenários e modelos. Mesmo um bom ajuste não identifica seu intercepto como líquido: tamanho tridimensional, deslizamento, refração e profundidade não estão determinados. Os interceptos ficaram fora dos resultados de velocidade de fase.','',
                'A conservação oferece o limite **α_b ≤ U_l,pistão/c_bolha ≤ 1**, condicionado a célula rígida, pistão sem gás, geometria correta e filme com 0 ≤ U_f ≤ c_bolha. A aplicação prespecificada a contornos grandes encontrou **0 trajetórias elegíveis**; nenhum limite líquido foi atribuído. Os meniscos robustos melhoraram a translação axial, mas mantiveram contato de máscara com as paredes e, em alguns casos, α_b=1 aparente. Isso não confirma espessura do filme nem autoriza reaproveitar a mesma máscara como geometria validada. A média de toda a célula ainda exigiria sua identificação completa.','',
                '## Validação, limitações e decisão','',
                'Foram usados controles sintéticos de translação, ruído de localização, crescimento, aceleração, imagem estacionária, ambiguidade, reversão temporal e iluminação compartilhada. O resíduo recuperou as velocidades sintéticas prescritas com erro máximo de 0,015875 px/quadro; isso é erro nesse controle, não acurácia experimental. Quarenta e cinco picos reais foram recalculados diretamente por Pearson, com diferença máxima de 8,05×10⁻⁸. A auditoria independente também verificou trajetórias, correspondências, decisões e redetecções; os relatórios completos preservam suas coberturas e limitações.','',
                'As gravações duram apenas 71–125 ms. Exposição ocupa 43,4–66,0% de um intervalo. Há refração, sobreposição fora do plano, saturação e seleção dos objetos mais fáceis. Concordância de algoritmos, máscaras e razões não substitui um traçador material. de Oliveira e Dias calculam velocidades/escorregamento mediante redução termodinâmica ou fechamento; os PDFs não validam essa atribuição a movimentos das imagens atuais.','',
                '**Recomendação:** conservar translação de bolhas/envelopes e comparações ópticas bem definidas como resultados desta base. Para líquido/gás materiais, manter os campos ausentes. Preencher todos os pontos exigiria hipóteses adicionais não validadas; não seria correto apresentar esse preenchimento como medição ou como razão de fase determinada.','',
                '## Arquivos','',
                '- [Tabela final dos 26 cenários](../results/final_per_run_status.csv), [razões em intervalos comuns](../results/paired_object_consensus.csv), [comparações entre gravações](../results/inter_recording_object_ratios.csv).',
                '- [Protocolo e histórico das opções](PROTOCOL.md), [revisão de resíduos](residual_review_PT_BR.md), [revisão física e literatura](physical_identifiability_review_PT_BR.md), [aplicação dos limites](conditional_piston_PT_BR.md).',
                '- [Objetos pequenos nos quadros reais](../figures/moving_object_examples_PT_BR.png), [auditoria de objetos quase imóveis](../figures/object_stationary_audit.png), [extrapolações completas](../results/zero_size_extrapolation.csv).',
                '- [Reprodução](../README.md), [verificação final](../results/final_verification.json). Figuras principais: PNG de 300 dpi e versões vetoriais PDF/SVG; texto em PT-BR.']
        else:
            lines += ['',
                'The main table retains threshold-stable results only. Object/envelope columns are medians over their separately selected trajectories. The ratio instead aggregates fitted velocities present on common time intervals; it is not the quotient of those two medians.','',
                '![Common-interval object comparison](../figures/paired_objects_PT_BR.png)','',
                '**m05q05_001 (slug):** ratio **0.9345**, segmentation sensitivity **0.9345–0.9924**, using 120/500 reference intervals. **m07q05_001 (plug):** **0.9639**, sensitivity **0.9444–0.9817**, using 133/500 intervals. These ranges are not confidence intervals. Fits use whole trajectories; only aggregation is restricted to common intervals. Near-unity ratios do not establish negligible phase slip: the small candidates can themselves be bubbles.','',
                '![Between-recording object-speed ratios](../figures/inter_recording_ratios_PT_BR.png)','',
                'A common small-object definition supports ten comparisons between five recordings: three within plug and seven between recordings of different patterns. These are descriptive comparisons of selected object populations, not liquid speeds, causal pattern effects, or full gas-phase means. Physical interpretation requires equal axial magnification. A single dispersed-bubble recording supplies no between-recording replication for that pattern.','',
                '## Tested options and selection','',
                'Closed-contour tracking recovered small objects but failed on large plug/slug contours. Wall-connected shadows compromised geometry. Axial core-envelope tracking recovered menisci, but localization jumps defeated the original trajectory gates. Robust Theil–Sen fits to both endpoints, with deformation and temporal-half agreement checks, recovered three threshold-stable recordings. This is the most useful tested option for **elongated-envelope translation** in that subset, not a validated gas material mean.','',
                'The additional motion screen removed 29 reference candidates with near-zero motion or insufficient fitted travel. Visual examples recur near image column 370 across recordings, consistent with fixed optical structure. Original tables retain them. Selection, segmentation, weighting and object-population changes remain potential biases.','',
                'The interface-distance predictor removed specified synthetic halos while preserving independent interior textures. Shared moving illumination nevertheless passed both masks without material transport. Real rectangular mask errors also propagated into residual texture. Residual subtraction supported the same pair counts as the fixed background (2/26 and 4/26), in different subsets, and did not identify an additional liquid observable. None of the nine wavy-annular recordings had a supported residual pair under the main rule.','',
                'Zero-size extrapolations compared constant, linear, quadratic and radial-adjusted quadratic models. Temporal-group validation purged training trajectories that shared native intervals with test trajectories. Size/speed associations depended on model and recording. Even predictive fit cannot identify the intercept as liquid motion without a justified slip model, size calibration and optical depth information. No intercept was promoted to a phase velocity.','',
                'Conservation gives **alpha_body <= U_liquid_piston/c_bubble <= 1** only for a rigid cell, a gas-free piston, valid geometry and film mean speed between zero and bubble speed. The prespecified closed-large-object application found **zero eligible trajectories**, so no liquid bounds were assigned. Robust menisci improve axial translation but retain substantial radial mask contact, sometimes apparent alpha_body=1. This does not establish film thickness or validate the geometry. A whole-cell liquid mean additionally requires identification of the complete cell.','',
                '## Validation and limitations','',
                'Controls covered known translation, localization noise, growth, acceleration, stationarity, ambiguous correspondence, time reversal and shared illumination. Maximum residual-method error on the prescribed synthetic velocities was 0.015875 px/frame; this is not experimental accuracy. Forty-five saved real peaks were checked by direct Pearson calculation, maximum difference 8.05e-8. Independent trajectory and detection checks document their exact scope separately.','',
                'Recordings last only 71–125 ms; exposures span 43.4–66.0% of a frame interval. Refraction, out-of-plane overlap, saturation and selection of easy objects remain. Agreement of algorithms or ratios does not identify material transport. The de Oliveira and Dias papers obtain velocities/slip through thermodynamic reduction or closure, and do not validate these image motions as the two material phase velocities.','',
                '**Recommendation:** retain well-defined bubble/envelope translation and optical comparisons as the results supported by this dataset. Leave material liquid/gas fields missing. Filling every point would require additional unvalidated assumptions and must not be presented as measured or identified phase slip.','',
                '## Files','',
                '- [Final 26-recording table](../results/final_per_run_status.csv), [common-interval ratios](../results/paired_object_consensus.csv), [between-recording comparisons](../results/inter_recording_object_ratios.csv).',
                '- [Protocol and option history](PROTOCOL.md), [residual review](residual_review_PT_BR.md), [physical/literature review](physical_identifiability_review_PT_BR.md), [conditional-bound application](conditional_piston_PT_BR.md).',
                '- [Native moving-object examples](../figures/moving_object_examples_PT_BR.png), [stationary-candidate audit](../figures/object_stationary_audit.png), [complete extrapolations](../results/zero_size_extrapolation.csv).',
                '- [Reproduction](../README.md), [final verification](../results/final_verification.json). Main figures are provided in Portuguese as 300-dpi PNG and vector PDF/SVG.']
        name='EXECUTIVE_SUMMARY_PT_BR.md' if pt else 'EXECUTIVE_SUMMARY.md'
        (BASE/'findings'/name).write_text('\n'.join(lines)+'\n')

if __name__=='__main__':run()
