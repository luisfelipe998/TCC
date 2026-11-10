#!/usr/bin/env python3
"""Native-frame visual diagnostics for the interface-distance residual."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parents[1]/'intermediate/residual_mplcache'))
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import residual_motion as study

def run():
    metadata=[]
    for run in ('m05q10_001','m04q80_001','m07q01_001'):
        files=sorted((study.SRC2/'data_processed'/run).glob('ImgA*.jpg'))
        im=np.stack([np.array(Image.open(f).convert('L'),dtype=np.float32) for f in files])
        with np.load(study.MASKS/(run+'.npz'),allow_pickle=False) as z:
            mm=[np.unpackbits(z[k],axis=-1,count=1016,bitorder='big').astype(bool) for k in ('masks_060','masks_070','masks_080')]
        gas=np.logical_and.reduce(mm);liquid=~np.logical_or.reduce(mm)
        liquid[:,:18]=False;liquid[:,71:]=False
        gas=study.erode(gas,3)[:,18:71];liquid=study.erode(liquid,3)[:,18:71]
        data,diagnostics=study.crossfit_predictions(im[:,18:71],study.distance_bins(mm[1])[:,18:71])
        frame=200;native=im[frame,18:71]
        overlay=np.repeat(native[...,None]/255,3,axis=-1)
        for m,color in [(gas[frame],np.array([.05,.45,1.])),(liquid[frame],np.array([1.,.5,.05]))]:
            overlay[m]=.55*overlay[m]+.45*color
        fig,axes=plt.subplots(5,1,figsize=(15,6.3),sharex=True)
        panels=[(native,'Imagem nativa: lúmen nas linhas 18:71','gray',0,255),
            (overlay,'Máscaras interiores, erosão 3 px: azul = gás candidato; laranja = líquido candidato',None,None,None),
            (data['fixed_background'][frame],'Resíduo após fundo fixo treinado na outra metade','RdBu_r',-50,50),
            (native-data['interface_distance'][frame],'Intensidade prevista pela distância à segmentação e posição radial','gray',0,255),
            (data['interface_distance'][frame],'Resíduo após predição por distância à segmentação','RdBu_r',-50,50)]
        for ax,(values,title,cmap,vmin,vmax) in zip(axes,panels):
            ax.imshow(values,cmap=cmap,vmin=vmin,vmax=vmax,aspect='auto',interpolation='nearest')
            ax.set_yticks([]);ax.set_title(title,fontsize=10,loc='left')
        axes[-1].set_xlabel('Coordenada horizontal (pixels nativos)')
        fig.suptitle(f'{run}, quadro {frame}: inspeção do tratamento óptico',fontsize=13)
        fig.text(.5,.012,'Resíduos: escala comum −50 a +50 unidades de intensidade. Cores das máscaras não certificam a identidade física.',ha='center',fontsize=9)
        fig.tight_layout(rect=(0,.04,1,.96))
        for ext in ('png','pdf'):
            fig.savefig(study.BASE/'figures'/f'residual_example_{run}.{ext}',dpi=220,bbox_inches='tight')
        plt.close(fig)
        metadata.append(dict(run=run,frame=frame,full_native_frames_for_predictor=501,erosion_radius=3,
            gas_candidate_fraction=float(gas[frame].mean()),liquid_candidate_fraction=float(liquid[frame].mean()),
            residual_limits=[-50,50],crossfit=diagnostics))
        print(run,flush=True)
    study.dump('residual_visual_examples.json',metadata)

if __name__=='__main__':run()
