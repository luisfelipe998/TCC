"""Figura do pico WPS usando apenas os resultados já calculados na estatistica6."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/estatistica7-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT.parent/'estatistica6/results/a3_features.csv'
OUT=ROOT/'figures/selecionadas'
CONDITIONS=['prefix','center','suffix']


def main():
    data=pd.read_csv(SOURCE,float_precision='round_trip')
    rows=data[data.pattern.isin(['plug','slug'])&data.condition.isin(CONDITIONS)]
    assert len(rows)==18 and rows.run.nunique()==6
    assert rows.groupby('run').condition.nunique().eq(3).all()
    for condition,direction in [('prefix',-1),('suffix',1)]:
        frame=rows[rows.condition==condition]
        plug=frame[frame.pattern=='plug'].peak_px
        slug=frame[frame.pattern=='slug'].peak_px
        assert (plug.max()<slug.min()) if direction==-1 else (plug.min()>slug.max())
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':13,
        'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.dpi':180})
    fig,ax=plt.subplots(figsize=(10,6.1),layout='constrained')
    colors={'plug':'#8064a8','slug':'#31977b'}
    for pattern in ['plug','slug']:
        for k,run in enumerate(rows[rows.pattern==pattern].run.drop_duplicates()):
            frame=rows[rows.run==run].set_index('condition').loc[CONDITIONS]
            ax.plot(range(3),frame.peak_px,color=colors[pattern],marker=['o','s','^'][k],
                    ms=7,lw=1.8,ls='-' if pattern=='plug' else '--',label=f'{pattern} · {run}')
    ax.set(xticks=range(3),xticklabels=['Começo','Centro','Fim'],
           ylabel='Comprimento de onda do pico WPS (px)',ylim=(280,1320),xlim=(-.15,2.15),
           title='O pico de plug ficou abaixo de slug no começo e acima no fim')
    ax.grid(axis='y',alpha=.2)
    ax.legend(loc='upper center',bbox_to_anchor=(.5,-.09),ncol=3,fontsize=8.5,frameon=False)
    fig.get_layout_engine().set(rect=(0,.12,1,.88))
    fig.text(.02,.025,'Seis ensaios, três por padrão; trechos de 3.171 px e a mesma grade de 200 escalas entre 50 e 3.000 px.\nCada linha conecta regiões do mesmo ensaio. As janelas podem se sobrepor e não são novas aquisições.',fontsize=9,linespacing=1.5)
    OUT.mkdir(exist_ok=True)
    for ext in ['png','pdf']:
        fig.savefig(OUT/f'pico_wps_plug_slug.{ext}',bbox_inches='tight',facecolor='white')
    plt.close(fig)
    print('Figura selecionada criada a partir de 18 observações já arquivadas; nenhum resultado recalculado.')


if __name__=='__main__':main()
