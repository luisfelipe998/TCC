#!/usr/bin/env python3
"""Generate both executive summaries and paper captions from saved evidence."""
from pathlib import Path
import json
import math
import pandas as pd
from publish_latest_summary import write_legacy_report, publish_latest

BASE=Path(__file__).resolve().parents[1]
PT={'dispersed bubbles':'Bolhas dispersas','plug':'Plug','slug':'Slug','churn':'Churn','wavy annular':'Anular ondulado'}


def number(value,pt=False,digits=0):
    if pd.isna(value): return '—'
    text=f'{value:,.{digits}f}'
    return text.replace(',','X').replace('.',',').replace('X','.') if pt else text


def emit(path,lines):
    path.write_text(chr(10).join(lines)+chr(10))


def run():
    passage=pd.read_csv(BASE/'results/passage_per_run.csv')
    spatial=pd.read_csv(BASE/'results/annular_spatial_pairs.csv')
    regions=pd.read_csv(BASE/'results/annular_region_summary.csv')
    variants=pd.read_csv(BASE/'results/annular_variant_pairs.csv')
    stress=json.loads((BASE/'results/mass_balance_stress_checks.json').read_text())
    p=passage[passage.eligible_passage_comparison]
    s=spatial[spatial.robust_spatial_descriptor]
    outer=variants[(variants.variant=='outer_bands')&variants.recommended_supported]
    annular=regions[(regions.pattern=='wavy annular')&(regions.region=='core')]
    core=annular[annular.robust_ensemble_supported]
    counts=[]
    for pattern in PT:
        counts.append(dict(pattern=pattern,recordings=int((passage.pattern==pattern).sum()),
                           passage_pairs=int((p.pattern==pattern).sum()),regional_pairs=int((s.pattern==pattern).sum()),outer_regional_pairs=int((outer.pattern==pattern).sum()),
                           supported_core_regions=int(((regions.pattern==pattern)&(regions.region=='core')&regions.robust_ensemble_supported).sum())))
    pd.DataFrame(counts).to_csv(BASE/'results/support_by_pattern.csv',index=False)
    status=[]
    for row in passage.itertuples():
        sr=spatial[spatial.run==row.run].iloc[0]
        regs=regions[regions.run==row.run].set_index('region')
        inner_variant=variants[(variants.run==row.run)&(variants.variant=='inner_bands')].iloc[0]
        outer_variant=variants[(variants.run==row.run)&(variants.variant=='outer_bands')].iloc[0]
        assert bool(sr.robust_spatial_descriptor)==bool(inner_variant.recommended_supported)
        reasons=inner_variant.recommended_withheld_reasons if pd.notna(inner_variant.recommended_withheld_reasons) else ''
        status.append(dict(run=row.run,pattern=row.pattern,passage_comparison=row.eligible_passage_comparison,
                           passage_reasons=row.withheld_reasons,regional_comparison=bool(sr.robust_spatial_descriptor),
                           regional_reasons=reasons,outer_regional_comparison=bool(outer_variant.recommended_supported),
                           outer_regional_reasons=outer_variant.recommended_withheld_reasons,
                           core_optical_propagation=bool(regs.loc['core','robust_ensemble_supported']),
                           gas_liquid_material_pair=False,superficial_velocity_pair=False))
    pd.DataFrame(status).to_csv(BASE/'results/all_run_status.csv',index=False)
    for pt in (False,True):
        lang='PT_BR' if pt else 'EN'
        lines=[
            '# RESUMO EXECUTIVO — v15.1: movimento óptico e comparação gás–líquido' if pt else '# EXECUTIVE SUMMARY — v15.1: optical motion and gas–liquid comparison',
            '',
            ('**A unidade principal recomendada é px/s.** A taxa de aquisição foi confirmada, mas não foi encontrada calibração axial em m/pixel. A análise produziu comparações por passagem em **'+str(len(p))+'/26 gravações** e comparações entre regiões em **'+str(len(s))+'/26**. No anular, **'+str(len(core))+'/9 núcleos** apresentaram propagação óptica com suporte nos critérios finais. Nenhum desses resultados identifica, por si só, duas velocidades materiais ou superficiais independentes.' if pt else
             '**Use px/s for the principal results.** Native frame timing was confirmed, but no axial m/pixel calibration was found. The analysis supports passage comparisons in **'+str(len(p))+'/26 recordings**, regional comparisons in **'+str(len(s))+'/26**, and core optical propagation in **'+str(len(core))+'/9 annular recordings** under the final criteria. None independently identifies both material or superficial phase velocities.'),
            '',
            ('Foram mantidas as exclusões: m05q05_002, anular liso e classificações duplas. As 26 gravações contêm 501 imagens nativas cada, totalizando 13.026 imagens e 13.000 intervalos adjacentes. Os registros duram somente 71–125 ms.' if pt else
             'The original exclusions remain: m05q05_002, smooth annular and double classifications. The 26 recordings contain 501 native images each: 13,026 images and 13,000 adjacent intervals. Each record spans only 71–125 ms.'),
            '',
            '## 1. Comparação solicitada: um ponto por gravação com suporte' if pt else '## 1. Requested comparison: one point per supported recording',
            '',
            f'![{"Comparação por passagem em px/s" if pt else "Passage comparison in px/s"}](../figures/gas_liquid_passages_px_s_{lang}.png)',
            '',
            ('O eixo horizontal representa o módulo da média com sinal do movimento durante passagens com aparência de gás; o vertical, durante intervalos com aparência de líquido. Todos os deslocamentos aceitos das quatro comparações apontam para a esquerda. A linha tracejada indica módulos iguais e não delimita padrões. Os pontos não representam Jg e Jl.' if pt else
             'The horizontal axis is the magnitude of signed mean image motion during gas-like passages; the vertical axis uses liquid-looking gaps. All accepted displacements in these four comparisons are leftward. The dashed line marks equal magnitudes and is not a pattern boundary. The coordinates are not Jg and Jl.'),
            '',
            ('| Gravação | Padrão | Média no gás aparente, px/s | Média no líquido aparente, px/s | Razão das médias | Pares gás/líquido |' if pt else
             '| Recording | Pattern | Gas-appearance mean, px/s | Liquid-appearance mean, px/s | Ratio of means | Gas/liquid pairs |'),
            '|---|---|---:|---:|---:|---:|']
        for r in p.itertuples():
            lines.append(f'| {r.run} | {PT[r.pattern] if pt else r.pattern} | {number(r.gas_mean_signed_px_s,pt)} | {number(r.liquid_mean_signed_px_s,pt)} | {number(r.ratio_of_signed_means,pt,3)} | {r.gas_pairs}/{r.liquid_pairs} |')
        lines += ['',
            ('A referência foi preservada e reagregada diretamente dos pares: concordância dos limiares v14 0,60/0,70/0,80, verificação da aparência líquida, episódios de pelo menos cinco quadros e consistência no intervalo de dois quadros. Exigiram-se pelo menos 20 pares e três episódios por tipo de passagem. A classificação ocorreu em uma janela de 24×53 pixels; a correlação do movimento utilizou uma região de 536 pixels de comprimento que pode conter ambas as fases.' if pt else
             'The existing reference was preserved and reaggregated directly from pairs: agreement across v14 thresholds 0.60/0.70/0.80, a liquid-appearance guard, episodes of at least five frames and lag-2 consistency. Each passage type required at least 20 pairs and three episodes. Appearance was classified in a 24×53-pixel window; motion correlation used a 536-pixel template that may include both phases.'),
            '',
            ('A média com peso igual por episódio fornece razões '+', '.join(number(r.episode_equal_ratio,True,3) for r in p.itertuples())+', na mesma ordem. A última cruza a unidade. Essa sensibilidade impede concluir qual fase material foi mais rápida. Os testes anteriores de regiões internas e janelas adaptativas continuam sem um par material validado; pixels não marcados pelo detector não são líquido confirmado.' if pt else
             'Equal weighting of episodes gives ratios '+', '.join(number(r.episode_equal_ratio,False,3) for r in p.itertuples())+' in the same order. The last crosses unity. This sensitivity prevents a conclusion about which material phase is faster. Previous interior/adaptive-patch tests still establish no validated material pair; unmarked detector pixels are not verified liquid.'),
            '',
            '## 2. Alternativa espacial: núcleo e regiões junto às paredes' if pt else '## 2. Spatial alternative: core and near-wall regions',
            '',
            f'![{"Comparação regional em px/s" if pt else "Regional comparison in px/s"}](../figures/regional_motion_px_s_{lang}.png)',
            '',
            ('A nova correlação normalizada foi calculada para todas as imagens em cinco faixas radiais. A referência usou núcleo 35:54, faixa superior 24:35 e inferior 54:65; são fatias Python com limite final exclusivo. A coordenada junto à parede é a média das celeridades superior e inferior, somente quando ambas têm suporte e concordam. A estimativa de cada faixa é o deslocamento do pico da curva média de correlação, dividido pelo intervalo temporal: **não é a média aritmética das velocidades dos pares**.' if pt else
             'New normalized correlations used every native image in five radial bands. The reference uses core rows 35:54, upper 24:35 and lower 54:65, with Python-exclusive stop indices. The near-wall coordinate averages upper/lower propagation estimates only when both are supported and agree. Each band estimate is the peak displacement of the mean correlation curve divided by elapsed time: **it is not the arithmetic mean of pair velocities**.'),
            '',
            ('Foram testados intervalos de 1, 2 e 4 quadros, cinco blocos temporais e 19 embaralhamentos da ordem temporal em cada intervalo. A revisão adversarial identificou picos fracos em alguns intervalos maiores. O critério final também exigiu força e separação dos controles em todos os intervalos, consistência da reversão temporal e de duas posições axiais adicionais. As decisões do critério inicial foram preservadas para comparação. Esses critérios avaliam suporte e consistência, não exatidão física.' if pt else
             'Tests used lags of 1, 2 and 4 frames, five temporal blocks and 19 shuffled-order controls at each lag. Adversarial review found weak peaks at some longer lags. Final reporting additionally requires peak strength/control separation at every lag, time-reversal consistency and consistency at two additional axial template positions. Initial decisions remain saved for comparison. These are support/consistency criteria, not physical accuracy bounds.'),
            '',
            '## 3. Tratamento especial dos anulares' if pt else '## 3. Special handling of annular flow',
            '',
            f'![{"Alternativa espacial com faixas externas" if pt else "Spatial alternative with outer bands"}](../figures/regional_outer_motion_px_s_{lang}.png)',
            '',
            ('A alternativa fixa de faixas externas forneceu **'+str(len(outer))+'/26 pares**, incluindo todos os cinco padrões e uma gravação anular: m04q80, com núcleo −233.283 px/s e média das faixas externas −221.429 px/s. Ela usou a mesma definição em todas as gravações; não combinou faixas diferentes caso a caso. A referência com faixas internas forneceu '+str(len(s))+' pares e nenhum anular. Cobertura maior não estabelece maior exatidão física.' if pt else
             'The fixed outer-band alternative supplied **'+str(len(outer))+'/26 pairs**, covering all five patterns and one annular recording: m04q80, with core −233,283 px/s and outer-band average −221,429 px/s. It used the same definition for every recording, without mixing bands case by case. The inner-band reference supplied '+str(len(s))+' pairs and no annular pair. Higher coverage does not establish greater physical accuracy.'),
            '',
            f'![{"Resultados anulares por região" if pt else "Annular results by region"}](../figures/annular_region_support_{lang}.png)',
            '',
            ('A correlação média recuperou propagação coerente do brilho/estrutura no núcleo de '+str(len(core))+' gravações anulares, com módulos entre '+number(abs(core.signed_velocity_px_s).min(),True)+' e '+number(abs(core.signed_velocity_px_s).max(),True)+' px/s. As faixas junto às paredes foram mais sensíveis ao intervalo e à posição. A falta de textura exclusiva de uma fase continua impedindo chamar esses valores de velocidade do gás ou do líquido.' if pt else
             'Mean correlation recovered coherent core brightness/structure propagation in '+str(len(core))+' annular recordings, with magnitudes from '+number(abs(core.signed_velocity_px_s).min())+' to '+number(abs(core.signed_velocity_px_s).max())+' px/s. Near-wall bands were more sensitive to lag and location. The absence of phase-exclusive texture prevents assigning these values to gas or liquid material velocity.'),
            '',
            ('Também foram testadas faixas externas 18:29 e 60:71 e posições de gradientes superior/inferior como candidatos de interface. As alternativas foram aplicadas por uma regra fixa, e não escolhidas individualmente para completar o gráfico. A celeridade de uma onda de filme não equivale à velocidade média do filme: em um modelo h_t+[q(h)]_x=0, c=dq/dh, enquanto a média material é q/h. Não foi adotado um fator universal de conversão. m03q40 permaneceu sinalizado por saturação.' if pt else
             'Alternatives also tested outer bands 18:29 and 60:71 and upper/lower gradient positions as interface candidates. Each alternative follows a fixed rule rather than being selected separately to fill missing plot points. Film-wave celerity differs from film mean material speed: in h_t+[q(h)]_x=0, c=dq/dh while the material mean is q/h. No universal conversion factor was imposed. m03q40 remains flagged for saturation.'),
            '',
            ('**Limitação remanescente:** a filtragem passa-altas com escala gaussiana de 8 pixels mudou quatro dos oito resultados anulares de núcleo além da tolerância declarada; com 32 pixels, os oito concordaram. Em m03q10, a celeridade mudou de −31,81 para −10,19 px/quadro. A celeridade estimada depende da escala óptica enfatizada pelo pré-processamento; um pico persistente não identifica automaticamente a velocidade material. Nenhum par baseado nas posições de gradiente permaneceu após exigir evidência forte em todos os intervalos.' if pt else
             '**Remaining limitation:** high-pass filtering at Gaussian scale 8 pixels changed four of the eight annular core estimates beyond the declared tolerance; scale 32 agreed in all eight. In m03q10 the celerity changed from −31.81 to −10.19 px/frame. The estimated celerity depends on the optical scale emphasized by preprocessing; a persistent peak does not automatically identify material velocity. No gradient-location pair survived the requirement for strong evidence at every lag.'),
            '',
            '## 4. Cobertura e casos difíceis' if pt else '## 4. Coverage and difficult cases',
            '',
            ('| Padrão | Gravações | Pares por passagem | Faixas internas | Faixas externas | Núcleos com suporte |' if pt else
             '| Pattern | Recordings | Passage pairs | Inner-band pairs | Outer-band pairs | Supported core regions |'),
            '|---|---:|---:|---:|---:|---:|']
        for c in counts:
            lines.append(f'| {PT[c["pattern"]] if pt else c["pattern"]} | {c["recordings"]} | {c["passage_pairs"]} | {c["regional_pairs"]} | {c["outer_regional_pairs"]} | {c["supported_core_regions"]} |')
        lines += ['',
            ('Plug oferece o melhor suporte para a comparação por passagem. Em slug, m04q01 e m06q05 não fornecem intervalos líquidos suficientes pela regra de referência; relaxar o limiar muda a mistura amostrada. Em bolhas dispersas, o rastreamento de objetos isolados foi mantido como observável separado: 38 trajetórias candidatas, 1.038 intervalos de objetos e média −43.904 px/s no limiar 0,70. Vários objetos podem contribuir no mesmo intervalo entre quadros. Em churn, contornos sobrepostos e mistura visual prejudicam a identificação de fases.' if pt else
             'Plug offers the strongest passage-comparison support. Slug runs m04q01 and m06q05 lack sufficient liquid-looking intervals under the reference rule; relaxing the gate changes the sampled mixture. Dispersed-bubble tracking remains a separate observable: 38 candidate chains, 1,038 object intervals and mean −43,904 px/s at threshold 0.70. Several objects can contribute to the same frame interval. Churn has overlapping contours and mixed appearance that obstruct phase identification.'),
            '',
            ('As 26 gravações permanecem na [tabela de situação](../results/all_run_status.csv) e no [painel de cobertura](../figures/all_run_availability_PT_BR.png). O [gráfico diagnóstico com todas as candidatas](../figures/all_regional_candidates_PT_BR.png) apresenta resultados rejeitados com marcadores vazios e saturação com ×. Nesse diagnóstico, a coordenada de parede conserva a média dos dois picos até quando discordam; esses pontos não são estimativas aceitas. Não foram preenchidas velocidades ausentes por valores de outro padrão ou por zero.' if pt else
             'All 26 recordings remain in the [status table](../results/all_run_status.csv) and [coverage panel](../figures/all_run_availability_EN.png). The [all-candidate diagnostic](../figures/all_regional_candidates_EN.png) marks rejected results with hollow symbols and saturation with ×. In that diagnostic, the wall coordinate retains the average of two peaks even when they disagree; these points are not accepted estimates. Missing velocities were not filled from another pattern or with zero.'),
            '',
            '## 5. Conversão para m/s e inferência rejeitada' if pt else '## 5. Conversion to m/s and rejected inference',
            '',
            ('A auditoria nova inspecionou 78 arquivos de câmera, metadados de 13.026 JPEGs e 104 planilhas. FPS e proporção 1:1 dos pixels estão registrados; parâmetros de lente, referência e correção óptica não fornecem m/pixel. O campo ScaleFactor=1 pertence a canais de tensão desativados. A escala 13,68 µm/pixel de Dias pertence a uma campanha diferente e não foi transferida.' if pt else
             'A fresh audit inspected 78 camera files, metadata from 13,026 JPEGs and 104 workbooks. Native FPS and 1:1 pixel aspect are documented; lens/reference/correction metadata supply no m/pixel calibration. ScaleFactor=1 belongs to disabled voltage channels. Dias reports 13.68 µm/pixel for a different campaign, so that value was not transferred.'),
            '',
            'v_imagem = deslocamento_px × FPS / intervalo_em_quadros' if pt else 'v_image = displacement_px × FPS / frame_interval_count',
            '',
            ('A conversão geométrica opcional usa s=0,001/53=18,8679 µm/pixel. Ela pressupõe que o diâmetro físico de 1 mm corresponda aos 53 pixels radiais e que essa escala também valha axialmente, apesar da refração no tubo. O [gráfico condicional em m/s](../figures/gas_liquid_passages_conditional_m_s_PT_BR.png) multiplica os mesmos quatro pares por esse fator; não oferece validação adicional.' if pt else
             'The optional geometric conversion uses s=0.001/53=18.8679 µm/pixel. It assumes that the 1 mm physical diameter corresponds to the 53 radial cells and that this scale also applies axially despite tube refraction. The [conditional m/s chart](../figures/gas_liquid_passages_conditional_m_s_EN.png) multiplies the same four pairs by that factor; it supplies no additional validation.'),
            '',
            ('Para testar uma alternativa sem entalpia de entrada, rearranjou-se G=ρg α ug+ρl(1−α)ul usando vazão registrada, geometria nativa e movimento óptico substituído por ug. Foram calculados **702 casos**, todos rejeitados como velocidades materiais. No cenário de referência, a velocidade líquida algébrica atingiu **'+number(stress['baseline_candidate_liquid_max_m_s'],True,2)+' m/s**; a variação apenas do limiar alterou alguns resultados em até **'+number(stress['max_threshold_liquid_factor'],True,2)+' vezes**. Esses são resultados de um teste de hipóteses inadequadas, não medidas de escoamento. A área líquida quase nula amplifica o resultado; médias não sincronizadas e identidade material não verificada invalidam a interpretação.' if pt else
             'An inlet-enthalpy-free alternative rearranged G=ρg α ug+ρl(1−α)ul using recorded mass flow, native geometry and an optical proxy substituted for ug. **702 cases** were evaluated; all were rejected as material velocities. The baseline algebraic liquid value reached **'+number(stress['baseline_candidate_liquid_max_m_s'],False,2)+' m/s**, and changing the segmentation threshold alone varied some results by up to **'+number(stress['max_threshold_liquid_factor'],False,2)+' times**. These are outputs of an invalidated assumption test, not measured flow speeds. A nearly zero liquid-area denominator amplifies the result; unmatched averaging and unverified material identity invalidate the interpretation.'),
            '',
            ('O gráfico de referência usa velocidades superficiais: Jg=Qg/A e Jl=Ql/A. FPS não fornece Qg nem Ql. Dias calculou essas quantidades a partir de G e título de vapor; o artigo de de Oliveira descreve entrada sub-resfriada. As equações e as páginas originais foram verificadas na [revisão de calibração e literatura](calibration_review.md).' if pt else
             'The supplied example uses superficial velocities: Jg=Qg/A and Jl=Ql/A. FPS supplies neither Qg nor Ql. Dias computed those quantities from mass flux and vapor quality; de Oliveira describes a subcooled inlet. Equations and original pages were checked in the [calibration and literature review](calibration_review.md).'),
            '',
            '## 6. Recomendação e limites da validação' if pt else '## 6. Recommendation and validation limits',
            '',
            ('Recomenda-se publicar o gráfico por passagem em px/s com sua definição completa, acompanhado do mapa regional como análise distinta. Para anular, a correlação média com verificações ampliadas é a opção mais sustentada para estudar propagação óptica no núcleo; a velocidade média do líquido continua sem identificação. Não foi estabelecido um método vencedor em exatidão física nem uma razão universal por padrão.' if pt else
             'Publish the passage plot in px/s with its full definition, and keep the regional map as a separate analysis. For annular flow, mean correlation with the extended checks provides the strongest support for studying core optical propagation; mean liquid material velocity remains unidentified. No physical-accuracy winner or universal per-pattern ratio was established.'),
            '',
            ('A verificação independente conferiu classificações, critérios de aceitação, conversões e agregações dos 13.000 registros de pares salvos; não refez o rastreamento original de imagens. Também recalculou a ocupação e a geometria das 39.078 máscaras salvas por outro algoritmo e os 702 casos do teste de massa, sem importar os geradores correspondentes. Os testes sintéticos, de reversão, de posição e de embaralhamento avaliam consistência numérica e correspondência. Não validam fase material, escala óptica ou erro experimental. As figuras usam os mesmos dados em inglês e PT-BR e foram exportadas em PDF/SVG e PNG de 600 dpi. Os intervalos e blocos não são repetições experimentais independentes.' if pt else
             'Independent verification checked classifications, acceptance rules, unit conversions and aggregation of the 13,000 saved pair records; it did not independently rematch the original images. It also recomputed occupancy and geometry for all 39,078 saved masks with a separate algorithm and the 702 mass-stress cases, without importing their generating modules. Synthetic, reversal, location and shuffle checks assess numerical consistency and correspondence. They do not validate material identity, optical scale or experimental error. English and PT-BR figures use identical numerical data and are exported as PDF/SVG and 600-dpi PNG. Frame intervals and blocks are not independent experimental repetitions.'),
            '',
            ('Evidências: [revisão adversarial](adversarial_method_review.md), [verificação independente](../results/independent_verification.json), [dados dos pontos por passagem](../results/paper_passage_plot_data.csv), [dados dos pontos regionais](../results/paper_regional_plot_data.csv), [legendas para o artigo](FIGURE_CAPTIONS_PT_BR.md) e [reprodução](../README.md).' if pt else
             'Evidence: [adversarial review](adversarial_method_review.md), [independent verification](../results/independent_verification.json), [passage plot data](../results/paper_passage_plot_data.csv), [regional plot data](../results/paper_regional_plot_data.csv), [paper captions in PT-BR](FIGURE_CAPTIONS_PT_BR.md), and [reproduction](../README.md).')]
        write_legacy_report('EXECUTIVE_SUMMARY_PT_BR.md' if pt else 'EXECUTIVE_SUMMARY.md',lines)
    captions=[
        '# Legendas e arquivos para publicação — PT-BR','',
        'As figuras finais possuem versões EN e PT_BR, PNG de 600 dpi, PDF vetorial e SVG com texto. Os valores dos pontos são os mesmos nos dois idiomas. Os eixos lineares foram escolhidos para a faixa observada; a organização gás no eixo horizontal e líquido no vertical foi preservada somente para os descritores por passagem, sem copiar as grandezas superficiais do exemplo.','',
        '## Figura 1 — Comparação por passagem','',
        '[PNG](../figures/gas_liquid_passages_px_s_PT_BR.png) · [PDF](../figures/gas_liquid_passages_px_s_PT_BR.pdf) · [SVG](../figures/gas_liquid_passages_px_s_PT_BR.svg)','',
        f'Comparação dos módulos das médias com sinal do movimento da imagem durante passagens com aparência de gás e intervalos com aparência de líquido. Cada ponto corresponde a uma gravação que satisfez os critérios de classificação e correspondência ({len(p)} de 26). Cores e símbolos indicam o padrão fornecido; a linha tracejada representa igualdade dos módulos. Todos os deslocamentos aceitos desses pares ocorreram para a esquerda. A classificação local não garante pureza de fase na região maior usada para estimar movimento; os eixos não representam velocidades superficiais ou materiais.','',
        '## Figura 2 — Propagação regional','',
        '[PNG](../figures/regional_motion_px_s_PT_BR.png) · [PDF](../figures/regional_motion_px_s_PT_BR.pdf) · [SVG](../figures/regional_motion_px_s_PT_BR.svg)','',
        f'Celeridades ópticas no núcleo e nas faixas junto às paredes, obtidas do pico das curvas médias de correlação entre imagens. Foram apresentados {len(s)} pares que passaram pelos critérios ampliados de força de correlação, intervalos temporais, embaralhamento, blocos, reversão e posição axial. A estimativa junto às paredes é a média dos picos superior e inferior, exigindo suporte e concordância de ambos. As regiões espaciais não equivalem a máscaras materiais de gás e líquido.','',
        '## Figura 3 — Alternativa com faixas externas','',
        '[PNG](../figures/regional_outer_motion_px_s_PT_BR.png) · [PDF](../figures/regional_outer_motion_px_s_PT_BR.pdf) · [SVG](../figures/regional_outer_motion_px_s_PT_BR.svg)','',
        f'Comparação das celeridades ópticas no núcleo e nas faixas externas 18:29 e 60:71. A definição foi mantida nas 26 gravações e forneceu {len(outer)} pares com suporte, incluindo m04q80 do padrão anular ondulado. A coordenada vertical é a média dos picos das duas faixas após verificação de suporte e concordância. A mudança das faixas internas para externas altera o observável; estes resultados não devem ser tratados como velocidade material de gás e líquido nem combinados seletivamente com a Figura 2.','',
        '## Figura 4 — Anulares','',
        '[PNG](../figures/annular_region_support_PT_BR.png) · [PDF](../figures/annular_region_support_PT_BR.pdf) · [SVG](../figures/annular_region_support_PT_BR.svg)','',
        'Propagação óptica no núcleo e nas faixas internas superior e inferior das nove gravações anulares onduladas. Somente resultados que satisfizeram os critérios ampliados foram apresentados. Pontos ausentes indicam insuficiência de suporte, e não velocidade zero. O resultado mede propagação de estruturas de brilho ou interfaces; não foi aplicada uma conversão de velocidade de onda para velocidade material do filme. A filtragem passa-altas com escala gaussiana de 8 pixels mudou quatro das oito estimativas de núcleo além da tolerância; a escala de 32 pixels concordou nas oito. Essa sensibilidade ao pré-processamento permanece mesmo nos resultados apresentados.','',
        '## Figura suplementar — Todas as candidatas','',
        '[PNG](../figures/all_regional_candidates_PT_BR.png) · [PDF](../figures/all_regional_candidates_PT_BR.pdf) · [SVG](../figures/all_regional_candidates_PT_BR.svg)','',
        'Diagnóstico das estimativas regionais calculadas para as 26 gravações. Marcadores preenchidos indicam aprovação nos critérios; vazios indicam rejeição; × assinala saturação. Os números identificam cada gravação na chave lateral. Os valores rejeitados permanecem visíveis apenas para auditar o algoritmo. A coordenada junto à parede conserva a média dos picos superior e inferior mesmo quando discordam, tornando-a inadequada como velocidade representativa nesses casos. Não se inferiram fronteiras, tendências ou razões de velocidade material a partir desses pontos.','',
        '## Figura suplementar — Conversão geométrica condicional','',
        '[PNG](../figures/gas_liquid_passages_conditional_m_s_PT_BR.png) · [PDF](../figures/gas_liquid_passages_conditional_m_s_PT_BR.pdf) · [SVG](../figures/gas_liquid_passages_conditional_m_s_PT_BR.svg)','',
        'Reescala da Figura 1 por s=0,001/53 m/pixel, assumindo diâmetro interno de 1 mm e equivalência da escala radial aparente com a escala axial. A proporção 1:1 dos pixels não calibra os efeitos ópticos do tubo. A figura não estabelece velocidades físicas medidas nem corrige a incerteza de identificação de fases. Recomenda-se manter a Figura 1 em px/s no resultado principal.','',
        '[Dados por passagem](../results/paper_passage_plot_data.csv) · [Dados regionais aceitos](../results/paper_regional_plot_data.csv) · [Todas as candidatas](../results/paper_all_candidate_plot_data.csv) · [Situação das 26 gravações](../results/all_run_status.csv).']
    write_legacy_report('FIGURE_CAPTIONS_PT_BR.md',captions)
    published=publish_latest()
    report=dict(original_result_counts=dict(passage_pairs=len(p),regional_pairs=len(s),annular_core=len(core)),
                published_latest_reports=published,summaries=['EN','PT_BR'])
    metrics=BASE/'all_scenarios_adversarial/results/final_summary_metrics.json'
    if published and metrics.exists():report['latest_counts']=json.loads(metrics.read_text())['counts']
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    run()
