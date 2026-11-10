#!/usr/bin/env python3
"""Bilingual, publication-size figures with explicit optical observables."""
from pathlib import Path
import json
import os
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(BASE/'intermediate/paper_mplcache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MaxNLocator
from matplotlib.colors import ListedColormap
from matplotlib.transforms import Bbox

COLORS={'dispersed bubbles':'#E69F00','plug':'#0072B2','slug':'#009E73','churn':'#D55E00','wavy annular':'#CC79A7'}
MARKERS={'dispersed bubbles':'o','plug':'s','slug':'^','churn':'D','wavy annular':'v'}
PT={'dispersed bubbles':'Bolhas dispersas','plug':'Plug','slug':'Slug','churn':'Churn','wavy annular':'Anular ondulado'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':10,
                     'axes.titlesize':11,'legend.fontsize':8,'xtick.labelsize':8,'ytick.labelsize':8,
                     'svg.fonttype':'none','pdf.fonttype':42,'savefig.facecolor':'white'})


def save(fig,stem):
    for ext in ('png','svg','pdf'):
        fig.savefig(BASE/'figures'/f'{stem}.{ext}',dpi=600,bbox_inches='tight',
                    metadata={'Creator':'v15.1 reproducible image-motion analysis'} if ext=='pdf' else None)
    plt.close(fig)


def fmt(pt, decimals=0):
    def tick(value,pos=None):
        text=f'{value:,.{decimals}f}'
        return text.replace(',','X').replace('.',',').replace('X','.') if pt else text
    return FuncFormatter(tick)


def place_labels(ax,annotations):
    """Keep labels inside the axes and away from other labels/markers."""
    fig=ax.figure
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    point_xy=np.array([ax.transData.transform(a.xy) for a in annotations])
    point_boxes=[Bbox.from_extents(x-5,y-5,x+5,y+5) for x,y in point_xy]
    occupied=[]
    if ax.get_legend() is not None:
        occupied.append(ax.get_legend().get_window_extent(renderer))
    area=ax.get_window_extent(renderer)
    offsets=[(x,y) for x in (-72,-48,-24,-8,8,24,48,72) for y in (-48,-36,-24,-12,0,12,24,36,48)]
    order=sorted(range(len(annotations)),key=lambda i:-sum(np.linalg.norm(point_xy-point_xy[i],axis=1)<65))
    for i in order:
        annotation=annotations[i]
        candidates=[]
        for dx,dy in offsets:
            annotation.set_position((dx,dy))
            annotation.set_ha('left' if dx>0 else 'right')
            annotation.set_va('center')
            box=annotation.get_window_extent(renderer).expanded(1.04,1.14)
            overlaps=sum(box.overlaps(old) for old in occupied)
            covered=sum(box.overlaps(point) for point in point_boxes)
            outside=not(area.contains(box.x0,box.y0) and area.contains(box.x1,box.y1))
            cost=100000*outside+10000*overlaps+1000*covered+np.hypot(dx,dy)
            candidates.append((cost,dx,dy,box))
        _,dx,dy,box=min(candidates,key=lambda item:item[0])
        annotation.set_position((dx,dy));annotation.set_ha('left' if dx>0 else 'right')
        occupied.append(box)
        ax.annotate('',annotation.xy,xytext=(dx,dy),textcoords='offset points',
                    arrowprops={'arrowstyle':'-','color':'.5','lw':.45,'shrinkA':3,'shrinkB':5},zorder=2)


def scatter(rows, kind, pt, conditional=False):
    language='PT_BR' if pt else 'EN'
    is_passage=kind=='passage'
    is_outer=kind=='regional_outer'
    scale=.001/53 if conditional else 1
    fig,ax=plt.subplots(figsize=(7.1,5.7))
    if is_passage:
        xfield,yfield='gas_plot_magnitude_px_s','liquid_plot_magnitude_px_s'
    else:
        xfield,yfield='core_magnitude_px_s','nearwall_magnitude_px_s'
    x,y=rows[xfield].to_numpy()*scale,rows[yfield].to_numpy()*scale
    lo=min(x.min(),y.min())*.83
    hi=max(x.max(),y.max())*1.13
    ax.plot([lo,hi],[lo,hi],color='.7',linestyle='--',linewidth=.8,zorder=1)
    label_offsets={
        'm05q01_001':(-10,15), 'm05q05_001':(8,9),'m06q01_001':(8,-14),
        'm07q05_001':(-8,14),'m05q10_001':(-10,12),'m05q20_001':(-10,-17),
        'm05q40_001':(8,-14),'m07q20_001':(8,10),'m07q40_001':(-8,11),
        'm07q01_001':(8,8)}
    if not is_passage:
        label_offsets.update({'m05q10_001':(-20,22),'m05q20_001':(-15,-32),
                              'm07q05_001':(-18,28),'m06q01_001':(15,-21),
                              'm05q05_001':(8,-18),'m07q20_001':(12,17)})
    annotations=[]
    for (_,r),xx,yy in zip(rows.iterrows(),x,y):
        ax.scatter(xx,yy,s=48,marker=MARKERS[r['pattern']],color=COLORS[r['pattern']],
                   edgecolors='white',linewidth=.5,zorder=3)
        dx,dy=label_offsets.get(r.run,(7,8))
        annotations.append(ax.annotate(r.run,(xx,yy),xytext=(dx,dy),textcoords='offset points',
                           ha='left' if dx>0 else 'right',fontsize=8,zorder=4))
    unit='m/s' if conditional else 'px/s'
    if is_passage:
        title='Movimento da imagem condicionado à passagem' if pt else 'Image motion conditioned on passage appearance'
        xlabel=('Passagens com aparência de gás'+chr(10)+'Módulo do movimento médio' if pt else 'Gas-like passages'+chr(10)+'Magnitude of mean image motion')+f' ({unit})'
        ylabel=('Intervalos com aparência de líquido'+chr(10)+'Módulo do movimento médio' if pt else 'Liquid-looking gaps'+chr(10)+'Magnitude of mean image motion')+f' ({unit})'
        caption=(f'{len(rows)}/26 gravações com os dois descritores; velocidades materiais não identificadas.' if pt
                 else f'{len(rows)}/26 recordings with both descriptors; material velocities remain unidentified.')
        stem='gas_liquid_passages_'+('conditional_m_s' if conditional else 'px_s')+'_'+language
    else:
        title='Propagação de estruturas em regiões da imagem' if pt else 'Structure propagation in image regions'
        xlabel=('Módulo da celeridade óptica no núcleo' if pt else 'Magnitude of optical propagation speed in core')+f' ({unit})'
        ylabel=('Módulo da celeridade óptica junto às paredes' if pt else 'Magnitude of optical propagation speed near walls')+f' ({unit})'
        caption=(f'{len(rows)}/26 gravações; regiões espaciais não identificam gás e líquido separadamente.' if pt
                 else f'{len(rows)}/26 recordings; spatial regions do not separately identify gas and liquid.')
        stem='regional_motion_px_s_'+language
        if is_outer:
            title='Propagação no núcleo e nas faixas externas' if pt else 'Propagation in core and outer bands'
            caption=(f'{len(rows)}/26 gravações; faixas 18:29 e 60:71; descritores ópticos, sem identificação material.' if pt
                     else f'{len(rows)}/26 recordings; bands 18:29 and 60:71; optical descriptors without material identification.')
            stem='regional_outer_motion_px_s_'+language
    ax.set(xlim=(lo,hi),ylim=(lo,hi),xlabel=xlabel,ylabel=ylabel,title=title)
    ax.set_aspect('equal',adjustable='box')
    ax.grid(alpha=.18,linewidth=.6)
    ax.xaxis.set_major_formatter(fmt(pt,1 if conditional else 0))
    ax.yaxis.set_major_formatter(fmt(pt,1 if conditional else 0))
    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(6))
    handles=[Line2D([],[],linestyle='',marker=MARKERS[p],color=COLORS[p],label=PT[p] if pt else p)
             for p in COLORS if p in set(rows.pattern)]
    ax.legend(handles=handles,loc='lower right',framealpha=.95)
    if conditional:
        caption=('ESCALA CONDICIONAL: 1 mm / 53 px; equivalência radial–axial e óptica não calibradas.' if pt
                 else 'CONDITIONAL SCALE: 1 mm / 53 px; radial–axial equivalence and optics are uncalibrated.')+chr(10)+caption
    fig.text(.5,.01,caption,ha='center',va='bottom',fontsize=8)
    fig.tight_layout(rect=(0,.065 if conditional else .045,1,1))
    place_labels(ax,annotations)
    save(fig,stem)


def annular_plot(regions, pt):
    language='PT_BR' if pt else 'EN'
    a=regions[regions.pattern=='wavy annular']
    names=sorted(a.run.unique())
    fig,ax=plt.subplots(figsize=(7.1,5.3))
    types=[('core',-.2,'o','#0072B2','Núcleo' if pt else 'Core'),
           ('upper_inner',0,'^','#009E73','Faixa superior interna' if pt else 'Upper inner band'),
           ('lower_inner',.2,'v','#D55E00','Faixa inferior interna' if pt else 'Lower inner band')]
    for region,offset,mark,color,label in types:
        for i,name in enumerate(names):
            row=a[(a.run==name)&(a.region==region)].iloc[0]
            if row.robust_ensemble_supported:
                ax.scatter(abs(row.signed_velocity_px_s),i+offset,marker=mark,color=color,s=40)
    for i,name in enumerate(names):
        supported=a[(a.run==name)&a.region.isin(['core','upper_inner','lower_inner'])].robust_ensemble_supported.any()
        if not supported:
            ax.text(.03,i,'Sem suporte' if pt else 'Unsupported',transform=ax.get_yaxis_transform(),va='center',fontsize=8,color='.4')
    ax.set_yticks(range(len(names)),names)
    ax.invert_yaxis()
    ax.set_xlim(0,290000)
    ax.xaxis.set_major_formatter(fmt(pt))
    ax.grid(axis='x',alpha=.2)
    ax.set_xlabel('Módulo da celeridade óptica (px/s)' if pt else 'Magnitude of optical propagation speed (px/s)')
    ax.set_title('Anular ondulado: suporte por região' if pt else 'Wavy annular flow: support by image region',pad=37)
    handles=[Line2D([],[],marker=m,color=c,linestyle='',label=l) for _,_,m,c,l in types]
    ax.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,1.01),ncol=3,frameon=False)
    caption=('Somente estimativas que passaram pelos critérios; posição ausente não significa velocidade zero.' if pt else
             'Only estimates passing the declared criteria; an absent point does not mean zero velocity.')
    fig.text(.5,.015,caption,ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.045,1,1))
    save(fig,'annular_region_support_'+language)


def availability(passage, spatial, regions, variants, pt):
    names=passage.run.tolist()
    values=[]
    for name in names:
        p=passage[passage.run==name].iloc[0]
        s=spatial[spatial.run==name].iloc[0]
        core=regions[(regions.run==name)&(regions.region=='core')].iloc[0]
        outer=variants[(variants.run==name)&(variants.variant=='outer_bands')].iloc[0]
        values.append([bool(p.eligible_passage_comparison),bool(s.robust_spatial_descriptor),bool(outer.recommended_supported),bool(core.robust_ensemble_supported)])
    values=np.array(values,dtype=int)
    fig,ax=plt.subplots(figsize=(7.1,8.2))
    ax.imshow(values,aspect='auto',cmap=ListedColormap(['#eeeeee','#0072b2']),vmin=0,vmax=1)
    labels=['Passagem','Faixas internas','Faixas externas','Núcleo'] if pt else ['Passage pair','Inner bands','Outer bands','Core alone']
    ax.set_xticks(range(4),labels)
    ax.set_yticks(range(len(names)),[n+'  ·  '+(PT[p] if pt else p) for n,p in zip(names,passage.pattern)])
    for i in range(len(names)):
        for j in range(4):
            ax.text(j,i,('Sim' if pt else 'Yes') if values[i,j] else '—',ha='center',va='center',fontsize=8,color='white' if values[i,j] else '.5')
    ax.set_title('Disponibilidade dos descritores por gravação' if pt else 'Descriptor availability by recording')
    ax.tick_params(axis='both',length=0)
    fig.text(.5,.01,'Critérios de suporte não representam validação da velocidade material.' if pt else 'Support criteria do not validate material phase velocity.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.025,1,1))
    save(fig,'all_run_availability_'+('PT_BR' if pt else 'EN'))


def all_candidate_plot(spatial, pt):
    """Show real algorithm outputs and failed gates, without endorsing them."""
    fig,(ax,key)=plt.subplots(1,2,figsize=(9.2,6.1),gridspec_kw={'width_ratios':[2.4,1]})
    table=spatial.sort_values('run').reset_index(drop=True)
    assert len(table)==26
    assert np.isfinite(table[['core_signed_px_s','nearwall_signed_px_s']].to_numpy()).all()
    assert (table.core_signed_px_s<0).all() and (table.nearwall_signed_px_s<0).all()
    ax.set(xlim=(20000,285000),ylim=(20000,250000))
    annotations=[]
    for i,r in table.iterrows():
        x,y=abs(r.core_signed_px_s),abs(r.nearwall_signed_px_s)
        ok=bool(r.robust_spatial_descriptor)
        ax.scatter(x,y,s=42,marker=MARKERS[r['pattern']],edgecolors=COLORS[r['pattern']],
                   facecolors=COLORS[r['pattern']] if ok else 'none',linewidth=1,zorder=3)
        if r.run=='m03q40_001':
            ax.scatter(x,y,s=24,marker='x',color='black',linewidth=1,zorder=4)
        annotations.append(ax.annotate(str(i+1),(x,y),xytext=(7,7),textcoords='offset points',
                                      fontsize=7,ha='left',va='bottom'))
    ax.set_xlabel('Candidato de propagação óptica no núcleo (px/s)' if pt else 'Candidate core optical propagation speed (px/s)')
    ax.set_ylabel('Candidato de propagação óptica junto às paredes (px/s)' if pt else 'Candidate near-wall optical propagation speed (px/s)')
    ax.set_title('Diagnóstico: todas as 26 gravações' if pt else 'Diagnostic: all 26 recordings')
    ax.xaxis.set_major_formatter(fmt(pt));ax.yaxis.set_major_formatter(fmt(pt))
    ax.xaxis.set_major_locator(MaxNLocator(6));ax.yaxis.set_major_locator(MaxNLocator(6))
    ax.grid(alpha=.18)
    key.axis('off')
    key.set_title('Identificação' if pt else 'Recording key')
    for i,r in table.iterrows():
        key.text(.02,1-(i+.8)/27,f'{i+1:02d}  {r.run}',fontsize=8,va='center',color=COLORS[r['pattern']])
    handles=[Line2D([],[],linestyle='',marker=MARKERS[p],color=COLORS[p],label=PT[p] if pt else p) for p in COLORS]
    fig.legend(handles=handles,ncol=5,loc='lower center',bbox_to_anchor=(.5,.085),frameon=False,fontsize=8)
    caption=('Preenchido: passou pelos critérios. Vazio: rejeitado. ×: saturação. Módulos das estimativas com sinal.' if pt else
             'Filled: passed the criteria. Hollow: rejected. ×: saturation. Magnitudes of signed estimates.')
    caption+=chr(10)+('Parede = média dos picos superior/inferior, inclusive quando discordam; não são velocidades de fase.' if pt else
                      'Wall = average of upper/lower peaks, including disagreements; these are not phase velocities.')
    fig.text(.5,.015,caption,ha='center',va='bottom',fontsize=8)
    fig.tight_layout(rect=(0,.16,1,1))
    place_labels(ax,annotations)
    table.insert(0,'figure_id',np.arange(1,27))
    table.to_csv(BASE/'results/paper_all_candidate_plot_data.csv',index=False)
    save(fig,'all_regional_candidates_'+('PT_BR' if pt else 'EN'))


def run():
    passage=pd.read_csv(BASE/'results/passage_per_run.csv')
    spatial=pd.read_csv(BASE/'results/annular_spatial_pairs.csv')
    regions=pd.read_csv(BASE/'results/annular_region_summary.csv')
    variants=pd.read_csv(BASE/'results/annular_variant_pairs.csv')
    chosen=passage[passage.eligible_passage_comparison].copy()
    assert 'robust_spatial_descriptor' in spatial and 'robust_ensemble_supported' in regions
    regional=spatial[spatial.robust_spatial_descriptor].copy()
    outer=variants[(variants.variant=='outer_bands')&variants.recommended_supported].copy()
    assert (chosen.gas_mean_signed_px_s<0).all() and (chosen.liquid_mean_signed_px_s<0).all()
    assert (regional.core_signed_px_s<0).all() and (regional.nearwall_signed_px_s<0).all()
    regional['core_magnitude_px_s']=abs(regional.core_signed_px_s)
    regional['nearwall_magnitude_px_s']=abs(regional.nearwall_signed_px_s)
    assert (outer.core_signed_px_s<0).all() and (outer.nearwall_signed_px_s<0).all()
    outer['core_magnitude_px_s']=abs(outer.core_signed_px_s)
    outer['nearwall_magnitude_px_s']=abs(outer.nearwall_signed_px_s)
    chosen.to_csv(BASE/'results/paper_passage_plot_data.csv',index=False)
    regional.to_csv(BASE/'results/paper_regional_plot_data.csv',index=False)
    outer.to_csv(BASE/'results/paper_regional_outer_plot_data.csv',index=False)
    for pt in (False,True):
        scatter(chosen,'passage',pt)
        scatter(chosen,'passage',pt,conditional=True)
        scatter(regional,'regional',pt)
        scatter(outer,'regional_outer',pt)
        annular_plot(regions,pt)
        availability(passage,spatial,regions,variants,pt)
        all_candidate_plot(spatial,pt)
    meta=dict(passage_points=len(chosen),regional_points=len(regional),outer_regional_points=len(outer),all_candidate_points=len(spatial),languages=['EN','PT_BR'],
              formats=['PNG 600dpi','SVG','PDF'],phase_axis_units='px/s',regional_axis_units='px/s',
              conditional_scale_m_px=.001/53,log_axes=False,
              display='magnitude of signed mean for passages; magnitude of ensemble propagation peak for regional descriptors; all plotted directions leftward',
              no_pattern_boundaries_drawn=True,regional_descriptors_are_not_phase_velocities=True)
    (BASE/'results/paper_figure_metadata.json').write_text(json.dumps(meta,indent=2)+chr(10))
    print(json.dumps(meta,indent=2))


if __name__=='__main__':
    run()
