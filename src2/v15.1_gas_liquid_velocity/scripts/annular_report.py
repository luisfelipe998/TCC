#!/usr/bin/env python3
"""Summarize spatial/annular diagnostics and verify their saved arithmetic."""
from pathlib import Path
import csv
import hashlib
import json
import os
import numpy as np
import pandas as pd

BASE=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(BASE/'intermediate/annular_mplcache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D


def variants(regions):
    rows=[]
    for run,group in regions.groupby('run',sort=False):
        d=group.set_index('region');core=d.loc['core']
        for variant,u,l in [('inner_bands','upper_inner','lower_inner'),
                            ('outer_bands','upper_outer','lower_outer'),
                            ('gradient_locations','upper_edge','lower_edge')]:
            if u not in d.index or l not in d.index:continue
            upper,lower=d.loc[u],d.loc[l]
            wall=(upper.signed_velocity_px_frame+lower.signed_velocity_px_frame)/2
            tolerance=max(2.,.25*abs(wall))
            agreement=bool(np.isfinite(wall) and np.sign(upper.signed_velocity_px_frame)==np.sign(lower.signed_velocity_px_frame)
                and abs(upper.signed_velocity_px_frame-lower.signed_velocity_px_frame)<=tolerance)
            direction=bool(np.isfinite(wall) and abs(wall)>=1 and np.sign(wall)==np.sign(core.signed_velocity_px_frame))
            base_reasons=[name+':'+str(row.withheld_reasons) for name,row in [('core',core),(u,upper),(l,lower)] if not row.ensemble_supported]
            if not agreement:base_reasons.append('upper_lower_disagreement')
            if not direction:base_reasons.append('core_wall_direction_or_zero')
            strict_reasons=list(base_reasons)
            for name,row in [('core',core),(u,upper),(l,lower)]:
                for field,reason in [('time_reversal_consistent','reversal'),('axial_template_consistent','axial_template'),('all_lag_peaks_supported','lag_peak_evidence')]:
                    if not row[field]:strict_reasons.append(name+':'+reason)
            rows.append(dict(run=run,pattern=core['pattern'],variant=variant,fps=core.fps,
                core_region='core',upper_region=u,lower_region=l,
                baseline_supported=not base_reasons,recommended_supported=not strict_reasons,
                baseline_withheld_reasons=';'.join(base_reasons),recommended_withheld_reasons=';'.join(strict_reasons),
                core_signed_px_frame=core.signed_velocity_px_frame,nearwall_signed_px_frame=wall,
                core_signed_px_s=core.signed_velocity_px_s,nearwall_signed_px_s=wall*core.fps,
                upper_signed_px_s=upper.signed_velocity_px_s,lower_signed_px_s=lower.signed_velocity_px_s,
                upper_lower_agreement=agreement,core_wall_same_direction=direction,
                core_to_nearwall_celerity_ratio=core.signed_velocity_px_frame/wall if abs(wall)>=1 else np.nan,
                interpretation='Peak-derived spatial optical celerities; wall value is the equal average of two regional celerities, not liquid material velocity'))
    frame=pd.DataFrame(rows)
    frame.to_csv(BASE/'results/annular_variant_pairs.csv',index=False)
    return frame


def verify(regions,lag,blocks,nulls,checks,pairs,metadata):
    packed=np.load(BASE/'results/annular_pair_results.npz')
    checked=0;max_error=0.
    for row in lag.itertuples():
        key=f'{row.run}__{row.region}'
        curve=packed[key+f'__curve{row.lag}'].astype(float)
        k=int(np.nanargmax(curve));offset=0.
        if 0<k<480:
            curvature=curve[k-1]-2*curve[k]+curve[k+1]
            if curvature < -1e-12:offset=np.clip(.5*(curve[k-1]-curve[k+1])/curvature,-.5,.5)
        delta=abs(k-240+offset-row.ensemble_shift_px)
        # The packed diagnostic curves retain float32; CSV estimates are float64.
        assert delta<.002,(key,row.lag,delta)
        max_error=max(max_error,delta)
        a=packed[key+f'__lag{row.lag}'].astype(float)
        threshold=.25 if 'edge' in row.region else 1.
        accepted=(a[:,1]>=.6)&(a[:,2]>=.03)&(a[:,3]>=threshold)&(np.abs(np.rint(a[:,0]))<238)
        assert np.array_equal(accepted,a[:,4].astype(bool)),(key,row.lag,'pair gates')
        assert accepted.sum()==row.individual_accepted_pairs
        n=nulls[(nulls.run==row.run)&(nulls.region==row.region)&(nulls.lag==row.lag)]
        assert len(n)==19
        assert np.isclose(n.peak.max(),row.shuffle_max_peak,rtol=1e-12,atol=1e-12)
        expected=row.ensemble_peak>=.10 and row.ensemble_peak>=row.shuffle_max_peak+.03 and row.valid_texture_pairs>=250 and not row.boundary
        assert expected==row.peak_and_shuffle_supported
        checked+=len(a)+4
    for row in regions.itertuples():
        ls=lag[(lag.run==row.run)&(lag.region==row.region)].set_index('lag')
        bs=blocks[(blocks.run==row.run)&(blocks.region==row.region)]
        cs=checks[(checks.run==row.run)&(checks.region==row.region)].set_index('check')
        assert len(ls)==3 and len(bs)==5 and len(cs)==6
        speed=row.signed_velocity_px_frame;tol=max(2.,.2*abs(speed))
        assert np.isclose(row.signed_velocity_px_s,speed*row.fps,rtol=1e-13)
        assert row.block_agree_count==int((abs(bs.signed_px_frame-speed)<=tol).sum())
        lag_agree=all(abs(ls.loc[k,'ensemble_velocity_px_frame']-speed)<=tol and not ls.loc[k,'boundary'] for k in (2,4))
        assert row.lag_consistent==lag_agree
        assert row.all_lag_peaks_supported==bool(ls.peak_and_shuffle_supported.all())
        reverse=abs(cs.loc['time_reversed','signed_px_frame']+speed)<=tol and not cs.loc['time_reversed','boundary']
        axial=all(abs(cs.loc[k,'signed_px_frame']-speed)<=tol and not cs.loc[k,'boundary'] for k in ('source_x160','source_x320'))
        assert row.time_reversal_consistent==reverse and row.axial_template_consistent==axial
        assert row.robust_ensemble_supported==bool(row.ensemble_supported and reverse and axial and row.all_lag_peaks_supported)
        checked+=8
    original=pd.read_csv(BASE.parent/'v15_gas_liquid_velocity/results/scenario_inventory.csv')
    selected=original[original.included]
    assert set(selected.run)==set(regions.run)==set(pairs.run)
    assert len(regions)==26*5+9*2 and len(pairs)==26*2+9
    assert metadata['script_sha256']==hashlib.sha256((BASE/'scripts/annular_motion.py').read_bytes()).hexdigest()
    result=dict(status='passed',comparisons_checked=checked,region_rows=len(regions),lag_rows=len(lag),
        null_trials=len(nulls),additional_checks=len(checks),variant_rows=len(pairs),
        maximum_shift_difference_from_float32_curves_px=max_error,
        limitation='Saved arithmetic, gates and provenance only; no physical accuracy or material phase identity verified')
    (BASE/'results/annular_verification.json').write_text(json.dumps(result,indent=2)+chr(10))
    return result


def figures(regions):
    ann=regions[regions.pattern=='wavy annular']
    runs=list(dict.fromkeys(ann.run))
    fig,ax=plt.subplots(figsize=(9.2,6.3))
    for k,(name,label,color,marker) in enumerate([('core','Núcleo: linhas 35–53','#2769a7','o'),
            ('upper_inner','Faixa superior: linhas 24–34','#9e731b','^'),('lower_inner','Faixa inferior: linhas 54–64','#218563','s')]):
        for i,run in enumerate(runs):
            r=ann[(ann.run==run)&(ann.region==name)].iloc[0]
            ok=r.robust_ensemble_supported
            ax.scatter(r.signed_velocity_px_s,i+(k-1)*.17,marker=marker,s=40,
                color=color if ok else '#bbbbbb',facecolors=color if ok else 'none',
                label=label if i==0 else None)
    ax.set_yticks(range(len(runs)),runs);ax.invert_yaxis()
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos:f'{x:,.0f}'.replace(',','.')))
    ax.set_xlabel('Celeridade aparente da estrutura óptica (px/s)')
    ax.grid(axis='x',alpha=.2)
    legend=[Line2D([],[],marker=marker,color=color,linestyle='',label=label) for label,color,marker in
        [('Núcleo: linhas 35–53','#2769a7','o'),('Faixa superior: linhas 24–34','#9e731b','^'),('Faixa inferior: linhas 54–64','#218563','s')]]
    ax.legend(handles=legend,loc='upper center',bbox_to_anchor=(.5,-.15),ncol=1,frameon=False,fontsize=9)
    fig.text(.5,.025,'Marcadores vazios: estimativa retida apenas como diagnóstico, sem suporte nos critérios finais.',ha='center',fontsize=8)
    ax.set_title('Escoamento anular ondulado: sensibilidade à região da imagem')
    fig.tight_layout(rect=(0,.07,1,1))
    for ext in ('png','svg','pdf'):fig.savefig(BASE/'figures'/f'annular_region_celerities_pt_br.{ext}',dpi=300)
    plt.close(fig)
    packed=np.load(BASE/'results/annular_pair_results.npz')
    chosen=['m03q10_001','m04q80_001','m06q80_001']
    fig,axes=plt.subplots(3,3,figsize=(10.7,7.2),sharex=True)
    for i,run in enumerate(chosen):
        for j,name in enumerate(['core','upper_inner','lower_inner']):
            ax=axes[i,j]
            for lag,color in [(1,'#2769a7'),(2,'#b26321'),(4,'#218563')]:
                curve=packed[f'{run}__{name}__curve{lag}']
                ax.plot(np.arange(-240,241)/lag,curve,color=color,lw=1,label='1 quadro' if lag==1 else f'{lag} quadros')
            ax.set_xlim(-80,5);ax.set_ylim(-.15,1)
            ax.grid(alpha=.2)
            if j==0:ax.set_ylabel(run+chr(10)+'Correlação normalizada')
            if i==0:ax.set_title({'core':'Núcleo','upper_inner':'Faixa superior','lower_inner':'Faixa inferior'}[name])
            if i==2:ax.set_xlabel('Deslocamento / intervalo (px/quadro)')
    axes[0,2].legend(fontsize=8,loc='upper left')
    fig.suptitle('Picos de correlação agregada em três intervalos temporais',fontsize=12)
    fig.tight_layout(rect=(0,0,1,.95))
    for ext in ('png','svg','pdf'):fig.savefig(BASE/'figures'/f'annular_correlation_examples_pt_br.{ext}',dpi=300)
    plt.close(fig)


def run():
    path=BASE/'results'
    regions=pd.read_csv(path/'annular_region_summary.csv')
    lag=pd.read_csv(path/'annular_lag_summary.csv')
    blocks=pd.read_csv(path/'annular_block_summary.csv')
    nulls=pd.read_csv(path/'annular_null_trials.csv')
    checks=pd.read_csv(path/'annular_additional_checks.csv')
    metadata=json.loads((path/'annular_metadata.json').read_text())
    pairs=variants(regions)
    verification=verify(regions,lag,blocks,nulls,checks,pairs,metadata)
    figures(regions)
    ann=regions[regions.pattern=='wavy annular']
    core=ann[ann.region=='core']
    inner=pairs[pairs.variant=='inner_bands'];outer=pairs[pairs.variant=='outer_bands'];edges=pairs[pairs.variant=='gradient_locations']
    lines=['# Spatial propagation and special annular handling','',
        f'**Aggregating normalized image correlations recovers a supported core optical celerity in {int(core.robust_ensemble_supported.sum())}/9 annular recordings under the final checks. It does not establish gas or liquid material velocity.** Spatial core/near-wall pairs are available in {int(inner.recommended_supported.sum())}/26 recordings with the fixed inner-band definition, and {int(outer.recommended_supported.sum())}/26 with the fixed outer-band alternative. Neither descriptor is a superficial gas/liquid velocity map.','',
        '## What was calculated','',
        'All 501 native images were used in each of the 26 retained recordings. The source configurations were independently checked for matching native FPS and consecutive frame IDs. The nine wavy-annular recordings also received upper/lower brightness-gradient tracking. No masks, inlet assumptions, flow-rate fit or pixel-to-length conversion enters these computations.','',
        '| Profile | Half-open image-row interval | Observable |','|---|---|---|',
        '| Core | 35:54 | Mean grayscale profile along the image |',
        '| Inner upper/lower bands | 24:35 and 54:65 | Separate mean grayscale profiles |',
        '| Outer upper/lower bands | 18:29 and 60:71 | Separate mean grayscale profiles closer to lumen limits |',
        '| Upper/lower gradient locations, annular only | 22:42 and 45:67 | Locations of strongest signed radial intensity gradients after smoothing |','',
        'The five grayscale-region definitions are fixed across all 26 runs. Their spatial names do not certify which material fills the projected image region. In annular flow, a central optical ray can also traverse the near and far liquid films. A gradient location can follow a shadow or reflection instead of the material interface.','',
        'For each profile, its 501-frame temporal mean at each axial coordinate was subtracted. A source template x=240:776 (536 pixels) was normalized and compared against all equal-width target windows over ±240 pixels. Per-pair NCC uses peak≥0.60, competing-peak margin≥0.03 beyond ±5 pixels, adequate texture and an interior search peak.','',
        'The alternative averages the full normalized correlation curves over all texture-qualified pairs before locating a peak. Its output is a **peak-derived propagation celerity**, not an arithmetic average of instantaneous displacement estimates. This can resolve coherent broad structures rejected by individual competing-peak margins; it can also privilege the dominant optical mode.','',
        '## Reporting checks and adversarial tightening','',
        'Baseline gates require at least 250 texture-qualified adjacent pairs, ensemble peak≥0.10, a peak at least 0.03 above the largest of 19 shuffled-order peaks, nonzero interior-search motion, ≤20% saturated regional pixels, lag-2/4 celerities within max(2 px/frame,20%) of lag 1, and agreement in at least four of five temporal blocks. The edge candidates additionally exclude >20% weak gradients or edge-search-boundary picks. These are diagnostic gates, not confidence limits.','',
        'The adversarial review identified that matching lag speeds alone was insufficient: the longer-lag correlation peak itself could be weak. Final recommended support therefore also requires peak≥0.10 and peak≥the same-lag maximum shuffled peak+0.03 at each lag (1, 2 and 4 frames), correct reversed-time celerity within tolerance, and agreement under both alternative axial source starts x=160 and x=320. Those templates keep width 536 and use their full valid asymmetric shift ranges. Baseline decisions remain saved separately.','',
        'The wall descriptor is the equal arithmetic average of the upper and lower **regional celerities**, available only if both regions pass, both directions agree, and their difference is ≤max(2 px/frame,25% of their average magnitude). A spatial pair additionally needs supported core motion in the same direction. No run-specific mixture of inner/outer/gradient definitions is used.','',
        '## Scope and support','',
        '| Pattern | Runs | Inner-band pairs, baseline/final | Outer-band pairs, baseline/final |','|---|---:|---:|---:|']
    for pattern in ['dispersed bubbles','plug','slug','churn','wavy annular']:
        a=inner[inner.pattern==pattern];b=outer[outer.pattern==pattern]
        lines.append(f'| {pattern} | {len(a)} | {int(a.baseline_supported.sum())}/{int(a.recommended_supported.sum())} | {int(b.baseline_supported.sum())}/{int(b.recommended_supported.sum())} |')
    lines+=['','## Every annular recording','',
        '| Run | Core celerity px/s | Core final support | Upper/lower inner final support | Outer-band spatial pair | Gradient spatial pair |','|---|---:|---|---|---|---|']
    for run,group in ann.groupby('run',sort=False):
        d=group.set_index('region');c=d.loc['core'];u=d.loc['upper_inner'];l=d.loc['lower_inner']
        o=outer[outer.run==run].iloc[0];e=edges[edges.run==run].iloc[0]
        value=f'{c.signed_velocity_px_s:,.0f}' if c.robust_ensemble_supported else 'withheld'
        lines.append(f'| {run} | {value} | {bool(c.robust_ensemble_supported)} | {bool(u.robust_ensemble_supported)}/{bool(l.robust_ensemble_supported)} | {bool(o.recommended_supported)} | {bool(e.recommended_supported)} |')
    lines+=['', '![Annular region celerities](../figures/annular_region_celerities_pt_br.png)','',
        'The hollow markers retain nonreportable estimates as diagnostic evidence only. Missing supported coordinates must not be replaced with these estimates in a main velocity scatter. Negative values indicate leftward propagation.','',
        '![Correlation peaks at different temporal intervals](../figures/annular_correlation_examples_pt_br.png)','',
        'The gradient-location alternative has one baseline candidate (m04q40), but **zero** final supported pairs: the longer-lag evidence is too weak. The final outer-band alternative supplies one annular spatial pair, m04q80, with core celerity −233,283 px/s and the mean of two outer-band celerities −221,429 px/s. Its inner-band alternative fails. This dependence on region definition must remain visible in any comparison.','',
        'High-pass scale 8 changes the selected core celerity beyond the declared 20%/2 px-per-frame tolerance in four of the eight otherwise supported annular recordings: m03q10, m03q20, m04q40 and m07q80. For example, m03q10 changes from −31.81 to −10.19 px/frame. Scale 32 retains agreement in all eight. These results identify an optical-scale dependence; they do not determine which optical mode follows the gas material.','',
        '## Interpretation and limitations','',
        '- The aggregate method recovers strong, persistent core brightness propagation where the earlier pair-level sharpness rule supplied very little support. For example, m06q80 has zero accepted individual core pairs, yet its aggregate peak is stable across lag, temporal block, source position and reversal checks. This is evidence for a coherent optical structure, not validated gas material tracking.',
        '- The inner near-wall bands can contain multiple scales with different apparent propagation. Their maxima may move from one mode to another when the temporal gap changes. This is why core support alone does not create a second phase velocity.',
        '- The outer-band alternative and gradient-location alternative are reported consistently as separate analyses. Any annular pair appearing only under an alternative definition is definition-sensitive, not a repaired gas/liquid slip measurement.',
        '- High-pass sensitivity at Gaussian scales 8 and 32 and a raw-profile alternative are saved. They change the optical structure being emphasized. A different celerity is method/scale sensitivity, not an error against physical ground truth.',
        '- Core and wall image signals can respond to the same moving boundary. Different or similar regional celerities do not prove independent phase-material velocities. Bubble/film/interface wave propagation, projected illumination and refraction remain possible sources.',
        '- m03q40 remains withheld for core saturation. A surviving dark outer strip does not restore the lost core information.',
        '- Native recording durations are only 71–125 ms; overlapping pairs and blocks are not independent experimental repeats. The shuffle maxima are null diagnostics, not calibrated significance levels or uncertainty intervals.','',
        '## Reproduction and arithmetic checks','',
        'Run python3 -B src2/v15.1_gas_liquid_velocity/scripts/annular_motion.py, then scripts/annular_report.py at the same absolute base. The --pilot option retains the two original pilot cases separately.','',
        f'Saved checks reproduce {verification["comparisons_checked"]:,} pair/summary operations, {verification["region_rows"]} region rows, {verification["lag_rows"]} lag rows, {verification["null_trials"]} shuffled trials and {verification["additional_checks"]} supplementary checks. Diagnostic correlation arrays are packed as float32; full-precision reported values remain in CSVs. The maximum reconstructed shift difference from this packing is {verification["maximum_shift_difference_from_float32_curves_px"]:.6f} px. These checks do not establish phase truth or physical accuracy.','',
        'Tables: [all regional outputs](../results/annular_region_summary.csv), [consistent pair variants](../results/annular_variant_pairs.csv), [lag support](../results/annular_lag_summary.csv), [all shuffle trials](../results/annular_null_trials.csv), [source/reversal/preprocessing sensitivity](../results/annular_additional_checks.csv), [temporal blocks](../results/annular_block_summary.csv), [metadata and synthetic controls](../results/annular_metadata.json), and [verification](../results/annular_verification.json).']
    (BASE/'findings/annular_methods.md').write_text(chr(10).join(lines)+chr(10))
    print(json.dumps(verification,indent=2))
    print(pairs.groupby('variant')[['baseline_supported','recommended_supported']].sum().to_string())


if __name__=='__main__':run()
