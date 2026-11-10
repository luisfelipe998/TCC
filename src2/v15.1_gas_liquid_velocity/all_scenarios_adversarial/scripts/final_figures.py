#!/usr/bin/env python3
"""All-recording optical figures with unchanged numerical support labels."""
from pathlib import Path
import sys,os,json,hashlib
import numpy as np
import pandas as pd
BASE=Path(__file__).resolve().parents[1]
OLD=BASE.parent
os.environ.setdefault('MPLCONFIGDIR',str(BASE/'intermediate/mplcache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap
sys.path.insert(0,str(OLD/'scripts'))
from paper_figures import place_labels,COLORS,MARKERS,PT
plt.rcParams.update({'font.size':9,'axes.labelsize':10,'axes.titlesize':11,'svg.fonttype':'none','pdf.fonttype':42})


def save(fig,stem):
    for ext in ('png','pdf','svg'):
        fig.savefig(BASE/'figures'/f'{stem}.{ext}',dpi=600,bbox_inches='tight')
    plt.close(fig)


def all_plot(data,pt):
    fig,(ax,key)=plt.subplots(1,2,figsize=(9.6,6.8),gridspec_kw={'width_ratios':[2.45,1]})
    notes=[]
    for r in data.itertuples():
        x=abs(r.plot_core_signed_px_s)/1000;y=abs(r.plot_wall_signed_px_s)/1000
        color=COLORS[r.pattern]
        ax.scatter(x,y,s=50,marker=MARKERS[r.pattern],edgecolors=color,linewidth=1.2,
                   facecolors=color if r.stable_all_variants else 'white',zorder=3)
        if not r.internal_supported:ax.scatter(x,y,s=48,marker='x',c='black',zorder=4)
        point_label=str(r.figure_id)+('†' if r.severe_saturation else '')
        notes.append(ax.annotate(point_label,(x,y),xytext=(6,7),textcoords='offset points',fontsize=7,zorder=5))
    ax.set(xlim=(0,290),ylim=(0,290))
    ax.set_aspect('equal',adjustable='box');ax.grid(alpha=.16)
    ax.set_xlabel('Velocidade óptica no núcleo (10³ px/s)' if pt else 'Core optical propagation speed (10³ px/s)')
    ax.set_ylabel('Velocidade óptica na parede selecionada (10³ px/s)' if pt else 'Selected-wall optical propagation speed (10³ px/s)')
    ax.set_title('Todas as 26 gravações — comparação óptica' if pt else 'All 26 recordings — optical comparison')
    key.axis('off');key.set_title('Gravação' if pt else 'Recording')
    for i,r in enumerate(data.itertuples()):
        key.text(.0,1-(i+.8)/27,f'{r.figure_id:02d}  {r.run}'+(' †' if r.severe_saturation else ''),fontsize=8,color=COLORS[r.pattern],va='center')
    patterns=[Line2D([],[],linestyle='',marker=MARKERS[p],color=COLORS[p],label=PT[p] if pt else p) for p in COLORS]
    fig.legend(handles=patterns,ncol=5,loc='lower center',bbox_to_anchor=(.5,.135),frameon=False,fontsize=8)
    if pt:
        caption=('Preenchido: estável nas variantes testadas. Vazio: sensível ou com suporte insuficiente. ×: cobertura insuficiente. †: saturação.'
                 '\nMódulos das médias sobre blocos aceitos da metade de avaliação; uma parede escolhida na outra metade.'
                 '\nAs regiões não identificam velocidades materiais de gás e líquido. Estabilidade não representa exatidão experimental.')
    else:
        caption=('Filled: stable across tested variants. Hollow: sensitive or insufficiently supported. ×: insufficient coverage. †: saturation.'
                 '\nMagnitudes of accepted assessment-block means; one wall selected using the other temporal half.'
                 '\nRegions do not identify material gas/liquid velocities. Stability does not establish experimental accuracy.')
    fig.text(.5,.02,caption,ha='center',va='bottom',fontsize=7.5,linespacing=1.55)
    fig.tight_layout(rect=(0,.19,1,1));place_labels(ax,notes)
    save(fig,'all_26_selected_wall_'+('PT_BR' if pt else 'EN'))


def coverage_plot(table,pt):
    columns=['old_outer','block_two_walls','selected_wall','phase_radius3','phase_radius1','spectral_joint','sensitivity_stable']
    labels=(['Anterior: duas paredes','Novo: blocos locais','Novo: parede selecionada','Máscaras: erosão 3 px','Máscaras: erosão 1 px','Espectro espaço–tempo','Parede: todas sensibilidades'] if pt else
            ['Previous: two walls','New: local blocks','New: selected wall','Masks: 3 px erosion','Masks: 1 px erosion','Space–time spectrum','Wall: all sensitivities'])
    values=table[columns].to_numpy(int)
    fig,ax=plt.subplots(figsize=(9.3,8.5))
    ax.imshow(values,aspect='auto',cmap=ListedColormap(['#f0f0f0','#0072B2']),vmin=0,vmax=1)
    ax.set_xticks(range(len(columns)),labels,rotation=30,ha='left');ax.xaxis.tick_top()
    ax.set_yticks(range(26),[f'{r.run} · {PT[r.pattern] if pt else r.pattern}' for r in table.itertuples()])
    for i in range(26):
        for j in range(len(columns)):
            ax.text(j,i,('Sim' if pt else 'Yes') if values[i,j] else '—',ha='center',va='center',fontsize=7,color='white' if values[i,j] else '.6')
    ax.tick_params(length=0)
    fig.text(.5,.013,'Critérios e grandezas diferem entre métodos; nenhuma coluna valida velocidades materiais.' if pt else
             'Methods have different criteria and observables; no column validates material phase velocities.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.03,1,1))
    save(fig,'method_coverage_'+('PT_BR' if pt else 'EN'))


def main():
    split=pd.read_csv(BASE/'results/split_background_pairs.csv')
    prev=pd.read_csv(BASE/'results/selected_wall_per_run.csv')
    reference=split[(split.variant=='training_background')&(split.selection_direction=='forward')].set_index('run').sort_index()
    fields=['core_signed_px_s','wall_signed_px_s']
    tolerance=(reference[fields].abs()*.2).clip(lower=reference.fps*2,axis=0)
    variants={'reverse':split[(split.variant=='training_background')&(split.selection_direction=='reverse')].set_index('run'),
              'raw':prev[(prev.variant=='raw_blocks10')&(prev.selection_direction=='forward')].set_index('run'),
              'highpass32':split[(split.variant=='training_background_hp32')&(split.selection_direction=='forward')].set_index('run'),
              'highpass8':split[(split.variant=='training_background_hp8')&(split.selection_direction=='forward')].set_index('run')}
    stability=reference.supported.copy();sens=[]
    for name,q in variants.items():
        match=((q[fields]-reference[fields]).abs()<=tolerance).all(axis=1)
        stability &= q.supported&match
        for run in reference.index:
            sens.append(dict(run=run,variant=name,reference_supported=bool(reference.loc[run,'supported']),
                alternative_supported=bool(q.loc[run,'supported']),within_tolerance=bool(match[run]),
                alternative_core_px_s=q.loc[run,'core_signed_px_s'],alternative_wall_px_s=q.loc[run,'wall_signed_px_s'],
                delta_core_px_s=q.loc[run,'core_signed_px_s']-reference.loc[run,'core_signed_px_s'],
                delta_wall_px_s=q.loc[run,'wall_signed_px_s']-reference.loc[run,'wall_signed_px_s'],
                same_selected_wall=bool(q.loc[run,'selected_wall']==reference.loc[run,'selected_wall'])))
    pd.DataFrame(sens).to_csv(BASE/'results/final_sensitivity.csv',index=False)
    regions=pd.read_csv(OLD/'results/annular_region_summary.csv').set_index(['run','region'])
    chart=reference.reset_index();chart.insert(0,'figure_id',range(1,27))
    chart['internal_supported']=chart.supported
    chart['stable_all_variants']=chart.run.map(stability)
    chart['severe_saturation']=[max(regions.loc[(r.run,'core'),'saturation_fraction'],regions.loc[(r.run,r.selected_wall),'saturation_fraction'])>.2 for r in chart.itertuples()]
    chart['plot_core_signed_px_s']=chart.core_signed_px_s.fillna(chart.all_test_core_candidate_px_s)
    chart['plot_wall_signed_px_s']=chart.wall_signed_px_s.fillna(chart.all_test_wall_candidate_px_s)
    chart['coordinate_source']=np.where(chart.core_signed_px_s.notna()&chart.wall_signed_px_s.notna(),'jointly_accepted_assessment_blocks','all_assessment_blocks_unaccepted_candidate')
    chart['observable']='core / training-selected wall optical propagation, not material phase speed'
    assert len(chart)==26 and np.isfinite(chart[['plot_core_signed_px_s','plot_wall_signed_px_s']].to_numpy()).all()
    assert (chart[['plot_core_signed_px_s','plot_wall_signed_px_s']]<0).all().all()
    chart.to_csv(BASE/'results/final_all_26_plot_data.csv',index=False)
    old=pd.read_csv(OLD/'results/all_run_status.csv').set_index('run')
    blocks=pd.read_csv(BASE/'results/block_ncc_per_run.csv');blocks=blocks[blocks.variant=='blocks10'].set_index('run')
    phase=pd.read_csv(BASE/'results/phase_local_pairs.csv')
    spectral=pd.read_csv(BASE/'results/spectral_pairs.csv').set_index('run')
    ledger=chart[['run','pattern']].copy()
    ledger['old_outer']=ledger.run.map(old.outer_regional_comparison)
    ledger['block_two_walls']=ledger.run.map(blocks.supported)
    ledger['selected_wall']=ledger.run.map(reference.supported)
    for radius in (3,1):ledger[f'phase_radius{radius}']=ledger.run.map(phase[phase.erosion_radius==radius].set_index('run').pair_optical_supported)
    ledger['spectral_joint']=ledger.run.map(spectral.robust_supported)
    ledger['sensitivity_stable']=ledger.run.map(stability)
    ledger['material_phase_pair']=False
    ledger.to_csv(BASE/'results/final_method_coverage.csv',index=False)
    counts={k:int(ledger[k].sum()) for k in ledger.columns if k not in ('run','pattern')}
    no_support=chart[~chart.internal_supported].run.tolist()
    sensframe=pd.DataFrame(sens)
    senscounts={v:{'coordinate_agreement':int((g.reference_supported&g.within_tolerance).sum()),
                   'support_and_coordinate_agreement':int((g.reference_supported&g.alternative_supported&g.within_tolerance).sum())}
                for v,g in sensframe.groupby('variant')}
    for pt in (True,False):all_plot(chart,pt);coverage_plot(ledger,pt)
    info={'counts':counts,'sensitivity_counts':senscounts,'insufficient_support':no_support,
          'all_points':26,'units':'px/s','plot_display_scale':1000,
          'background_fitting':'Training half only, frozen before assessment',
          'stability_rule':'Primary support AND alternative support AND both coordinates within max(2*FPS,20% magnitude) for reversed halves, raw, sigma32, sigma8.',
          'limitations':'Sensitivity compares different temporal/content samples; failure is not proof of erroneous matching. Agreement is not physical accuracy or phase identification.',
          'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (BASE/'results/final_summary_metrics.json').write_text(json.dumps(info,indent=2))
    print(json.dumps(info,indent=2))


if __name__=='__main__':main()
