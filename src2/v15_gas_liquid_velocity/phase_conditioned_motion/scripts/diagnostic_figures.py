#!/usr/bin/env python3
"""Figures from saved phase-conditioned evidence; no estimator fitting."""
from pathlib import Path
import csv
import os
import numpy as np
from PIL import Image
BASE=Path(__file__).resolve().parents[1]
DATA=BASE.parent.parent/'data_processed'
os.environ.setdefault('MPLCONFIGDIR',str(BASE/'intermediate'/'figure_mplcache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


def read(name):return list(csv.DictReader((BASE/'results'/name).open()))


def run():
    out=BASE/'figures';out.mkdir(exist_ok=True)
    frames=read('passage_frames.csv');pairs=read('passage_pairs.csv')
    selected=['m05q01_001','m05q05_001','m06q05_001','m07q01_001','m05q20_001','m04q80_001']
    colors={'unresolved':'#aaaaaa','liquid_candidate':'#298c63','gas_candidate':'#326bc2'}
    cmap=ListedColormap([colors[k] for k in ['unresolved','liquid_candidate','gas_candidate']])
    fig,axes=plt.subplots(len(selected),2,figsize=(13,11))
    for name,axs in zip(selected,axes):
        fs=[r for r in frames if r['run']==name and r['config']=='consensus']
        ps=[r for r in pairs if r['run']==name and r['config']=='consensus' and r['motion_gate']=='adjacent']
        fps=float(fs[0]['fps']);t=np.arange(501)/fps*1000
        occ=np.array([[float(r[k]) for r in fs] for k in ['occupancy_060','occupancy_070','occupancy_080']])
        axs[0].fill_between(t,occ.min(axis=0),occ.max(axis=0),color='#a1b6d6',alpha=.7)
        axs[0].plot(t,occ[1],color='#24519a',lw=.7)
        axs[0].axhline(.6,color='0.5',ls='--',lw=.6);axs[0].axhline(.05,color='0.5',ls='--',lw=.6)
        state=np.array([['unresolved','liquid_candidate','gas_candidate'].index(r['state']) for r in fs])
        axs[0].imshow(state[None,:],aspect='auto',extent=[0,500/fps*1000,-.16,-.06],cmap=cmap,vmin=0,vmax=2)
        axs[0].set_ylim(-.18,1.05);axs[0].set_ylabel(name+'\nProbe mask fraction',fontsize=9)
        for phase in colors:
            data=[r for r in ps if r['phase']==phase and (r['accepted']=='True' or phase=='unresolved')]
            axs[1].scatter([int(r['frame_index'])/fps*1000 for r in data],
                [float(r['velocity_px_s'])/1000 for r in data],s=4,color=colors[phase],alpha=.7,
                label=phase.replace('_',' ') if name==selected[0] else None)
        axs[1].set_ylabel('Structure motion\n1000 px/s',fontsize=9)
        for ax in axs:ax.grid(alpha=.2)
    axes[0,0].set_title('Fixed probe: thresholds 0.60–0.80 and retained state')
    axes[0,1].set_title('Global-profile motion conditioned on probe state')
    axes[0,1].legend(fontsize=8,loc='lower left')
    axes[-1,0].set_xlabel('Time (ms)');axes[-1,1].set_xlabel('Time (ms)')
    fig.suptitle('Gas-candidate and liquid-looking passages: retained ambiguity\nProbe: x=496:520, y=18:71; global motion uses a wider 536-pixel template',fontsize=12)
    fig.tight_layout(rect=(0,0,1,.95));fig.savefig(out/'passage_time_diagnostics.png',dpi=180);plt.close(fig)
    ratios=read('passage_ratios.csv');local=read('phase_patch_ratios.csv')
    names=['m05q01_001','m05q05_001','m06q01_001','m07q05_001']
    fig,ax=plt.subplots(figsize=(10,5.5))
    methods=[('consensus','adjacent','o',-.20,'#326bc2','Passage-conditioned, adjacent'),
        ('consensus','lag2_confirmed','s',-.05,'#66a3d1','Passage-conditioned, lag-2 screen')]
    for config,gate,mark,dx,col,label in methods:
        for i,name in enumerate(names):
            r=next(r for r in ratios if r['run']==name and r['config']==config and r['motion_gate']==gate)
            ax.scatter(i+dx,float(r['descriptive_ratio']),marker=mark,color=col,label=label if i==0 else None)
    for i,name in enumerate(names):
        r=next(r for r in local if r['run']==name and r['variant']=='majority_070' and r['method']=='background_subtracted')
        if r['eligible_descriptive_ratio']=='True':
            ax.scatter(i+.15,float(r['descriptive_ratio']),marker='^',color='#d36a27',label='Local patches, ≥80% phase candidate' if i==0 else None)
    ax.axhline(1,color='0.5',ls='--',lw=.8)
    ax.set_xticks(range(4),names);ax.set_ylabel('Ratio of signed conditional mean image motions')
    ax.set_ylim(.55,1.07);ax.grid(axis='y',alpha=.2)
    ax.set_title('The estimated ratio changes with motion support\nThese are optical descriptors, not measured gas/liquid slip')
    ax.legend(fontsize=9,loc='lower right')
    fig.text(.5,.02,'No run meets the two-phase reporting gates with 95% consensus or 100% eroded interiors.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.05,1,1))
    for ext in ['png','svg']:fig.savefig(out/('conditional_ratio_comparison.'+ext),dpi=180)
    plt.close(fig)
    summary=read('phase_patch_summary.csv');variants=['majority_070','consensus_95','interior_consensus_100']
    fig,axes=plt.subplots(1,2,figsize=(11,4.6))
    for ax,phase in zip(axes,['gas_candidate','liquid_candidate']):
        counts=[sum(int(r['n_accepted_tracks']) for r in summary if r['method']=='background_subtracted' and r['variant']==v and r['phase']==phase) for v in variants]
        bars=ax.bar(range(3),counts,color=colors[phase])
        ax.bar_label(bars,padding=3);ax.set_xticks(range(3),['≥80% baseline','≥95% consensus','100% eroded\nconsensus interior'],fontsize=9)
        ax.set_title(phase.replace('_',' '));ax.set_ylabel('Accepted local tracks, all 26 recordings')
        ax.set_ylim(0,max(counts)*1.18);ax.grid(axis='y',alpha=.2)
    fig.suptitle('Removing ambiguous/boundary regions sharply reduces liquid-candidate support')
    fig.tight_layout(rect=(0,0,1,.92));fig.savefig(out/'phase_support_attrition.png',dpi=180);plt.close(fig)
    # Inspect every distinct surviving strict liquid-candidate example.
    tr=read('phase_patch_tracks.csv')
    candidates=[r for r in tr if r['variant']=='interior_consensus_100' and r['method']=='background_subtracted' and r['assigned_phase']=='liquid_candidate' and r['accepted']=='True']
    assert len(candidates)==4
    fig,axes=plt.subplots(len(candidates),2,figsize=(12,6))
    for r,axs in zip(candidates,axes):
        i=int(r['frame_index']);j=int(r['next_frame_index']);x=int(r['window_x']);y=int(r['window_y'])
        bx=int(r['matched_x']);by=int(r['matched_y']);w=int(r['width']);h=int(r['height'])
        files=sorted((DATA/r['run']).glob('ImgA*.jpg'))
        for ax,frame,xx,yy,label in [(axs[0],i,x,y,'source'),(axs[1],j,bx,by,'matched target')]:
            a=np.asarray(Image.open(files[frame]).convert('L'))
            lo=max(0,xx-20);hi=min(1016,xx+w+20);top=max(0,yy-8);bottom=min(80,yy+h+8)
            ax.imshow(a[top:bottom,lo:hi],cmap='gray',vmin=0,vmax=255,interpolation='nearest',aspect='auto',extent=[lo,hi,bottom,top])
            from matplotlib.patches import Rectangle
            ax.add_patch(Rectangle((xx,yy),w,h,fill=False,edgecolor='#39b987',lw=1))
            ax.set_title(f'{r["run"]}, frame {frame}, {label}; dx={float(r["dx_px"]):.2f} px',fontsize=9)
            ax.set_ylabel('y (px)');ax.set_xlabel('x (px)')
    fig.suptitle('Only four strict liquid-candidate matches remain, all in the dispersed run\nGreen boxes are mask-interior matches; feature material identity is unverified',fontsize=11)
    fig.tight_layout(rect=(0,0,1,.90));fig.savefig(out/'strict_liquid_candidate_examples.png',dpi=200);plt.close(fig)
    print('Saved phase-passage, ratio, support-attrition and strict-patch diagnostics.')


if __name__=='__main__':run()
