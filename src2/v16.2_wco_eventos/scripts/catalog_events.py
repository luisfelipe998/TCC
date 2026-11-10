"""Anota propostas da imagem e gera pares sem acesso a wavelets."""
from pathlib import Path
import hashlib
import itertools
import json
import os
import time
os.environ.setdefault('MPLCONFIGDIR','/tmp/wco-eventos-mpl')
import numpy as np
import pandas as pd
from PIL import Image
from scipy.ndimage import median_filter, binary_closing
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from inspect_inputs import ROOT, BASE, DATA, LABELS, digest

CFG=json.loads((ROOT/'config.json').read_text())
# A lista registra a inspeção óptica; não é um classificador por rótulo global.
OBSERVABLE={'m05q01_001','m06q01_001','m07q05_001',
            'm04q01_001','m05q05_001','m06q05_001'}

def intervals(pano, threshold=60):
    f=(pano[18:72]<threshold).mean(axis=0)
    mask=median_filter(f,size=CFG['mediana_deteccao_px'],mode='nearest')>=CFG['fracao_escura_minima']
    size=CFG['fechamento_deteccao_px']
    mask=binary_closing(np.pad(mask,size,mode='edge'),structure=np.ones(size))[size:-size]
    diff=np.diff(np.r_[False,mask,False].astype(int))
    return [(int(a),int(b)) for a,b in zip(np.where(diff==1)[0],np.where(diff==-1)[0]) if b-a>=CFG['extensao_minima_px']],f

def frame_for(x, shifts):
    bounds=np.r_[0,np.cumsum(shifts)]
    i=int(np.searchsorted(bounds,x,side='right')-1)
    return min(i,len(shifts)-1), int(x-bounds[i])

def main():
    t0=time.perf_counter();out=ROOT/'results';figs=ROOT/'figures/catalogo';figs.mkdir(parents=True,exist_ok=True)
    inv=pd.read_csv(BASE/'results/input_audit.csv');inv=inv[inv.included].sort_values('run')
    reviewpath=out/'revisao_visual.csv'
    review=pd.read_csv(reviewpath).set_index('event_id').to_dict('index') if reviewpath.exists() else {}
    records=[];fixed=[];observability=[];frames_hashes={}
    for r in inv.itertuples(index=False):
        pano=np.asarray(Image.open(DATA/r.run/'__panorama_v11.png').convert('L'))
        cache=np.load(BASE/'intermediate'/f'{r.run}_input.npz');shifts=cache['shifts'];seams=cache['seams']
        raw=pano[18:72].mean(axis=0);proposals,f=intervals(pano)
        alternatives={q:intervals(pano,q)[0] for q in [55,65]}
        local=[];N=CFG['comprimento_segmento'];jitter=CFG['margem_sensibilidade_px']
        observable=r.run in OBSERVABLE
        for i,(a,b) in enumerate(proposals):
            eid=f'{r.run}_E{i+1:03d}';anchor=(a+b)//2;s=anchor-N//2;e=s+N
            reasons=[]
            if not observable:reasons.append('sem_extremidades_morfologicas_confirmadas')
            if a==0 or b==r.n:reasons.append('evento_truncado_no_panorama')
            if b-a>CFG['extensao_maxima_px']:reasons.append('extensao_maior_1536')
            if s<jitter or e>r.n-jitter:reasons.append('contexto_2048_mais_margem_indisponivel')
            if raw[max(0,s):min(r.n,e)].std()<1e-10:reasons.append('sinal_degenerado')
            if observable and eid in review and not review[eid]['aceito']:
                reasons.append(str(review[eid]['motivo']))
            reviewed=eid in review
            if observable and not reviewed:reasons.append('revisao_visual_pendente')
            anchor_frame,offset=frame_for(anchor,shifts)
            left_frame,_=frame_for(a,shifts);right_frame,_=frame_for(b-1,shifts)
            framefile=DATA/r.run/f'ImgA{anchor_frame:06d}.jpg'
            frames_hashes[str(framefile.relative_to(ROOT.parent))]=digest(framefile)
            frame_width=Image.open(framefile).width
            anchor_col=frame_width//2-int(shifts[anchor_frame])//2+offset
            # Um pixel no lado esquerdo da costura pode ser misturado com a faixa seguinte.
            mixed=bool(np.any((seams-anchor>0)&(seams-anchor<=2)))
            border_changes=[]
            for q,alt in alternatives.items():
                overlaps=[max(0,min(b,d)-max(a,c))/(max(b,d)-min(a,c)) for c,d in alt]
                if overlaps and max(overlaps)>=.5:
                    c,d=alt[int(np.argmax(overlaps))];change=max(abs(a-c),abs(b-d))
                else:change=np.nan
                border_changes.append(change)
            ctx_a=max(0,s);ctx_b=min(r.n,e)
            cf0,_=frame_for(ctx_a,shifts);cf1,_=frame_for(ctx_b-1,shifts)
            nonpos=np.r_[cache['raw_shifts']<=0,cache['raw_shifts'][-1]<=0]
            row=dict(event_id=eid,aquisicao=r.run,padrao=r.pattern,
                categoria='estrutura_alongada_delimitada' if observable else 'textura_ambigua',
                inicio=a,fim=b,extensao=b-a,ancora=anchor,linhas='18–71',
                segmento_inicio=s,segmento_fim=e,contexto_esquerdo=a-s,contexto_direito=e-b,
                distancia_costura_ancora=int(np.min(np.abs(seams-anchor))),
                costuras_evento=int(((seams>a)&(seams<b)).sum()),
                costuras_contexto=int(((seams>ctx_a)&(seams<ctx_b)).sum()),
                deslocamentos_nao_positivos_evento=int(nonpos[left_frame:right_frame+1].sum()),
                deslocamentos_nao_positivos_contexto=int(nonpos[cf0:cf1+1].sum()),
                quadro_ancora=framefile.name,coluna_quadro_ancora=anchor_col,
                ancora_em_zona_mistura=mixed,quadro_inicial=f'ImgA{left_frame:06d}.jpg',
                quadro_final=f'ImgA{right_frame:06d}.jpg',
                mudanca_fronteira_limiar55=border_changes[0],mudanca_fronteira_limiar65=border_changes[1],
                qualidade='revisada' if reviewed else 'ambigua' if not observable else 'pendente',
                incluida=not reasons,motivo=';'.join(reasons) if reasons else 'aceita',
                versao_regra=CFG['versao'],modo='limiar_de_imagem_e_revisao_visual_unica')
            local.append(row);records.append(row)
        starts=[0,(r.n-N)//2,r.n-N]
        for k,s in enumerate(starts):
            fixed.append(dict(segment_id=f'{r.run}_F{k}',aquisicao=r.run,padrao=r.pattern,
                              inicio=s,fim=s+N,ancora=s+N//2,posicao=['inicio','centro','fim'][k]))
        observability.append(dict(aquisicao=r.run,padrao=r.pattern,alongados_observaveis=observable,
            candidatos=len(local),eventos_incluidos=sum(x['incluida'] for x in local),
            motivo=('extremidades_visiveis_revisadas' if observable else
                    'objetos_pequenos_sem_segunda_aquisicao' if r.pattern=='dispersed bubbles' else
                    'textura_continua_fragmentada_sem_fronteiras_repetiveis'),
            atlas=f'figures/observabilidade/{r.run}.png'))
        if observable:
            chunks=int(np.ceil(r.n/2048));fig,axs=plt.subplots(chunks,1,figsize=(16,chunks*2.1),constrained_layout=True)
            for ax,k in zip(np.atleast_1d(axs),range(chunks)):
                x0=k*2048;x1=min((k+1)*2048,r.n)
                ax.imshow(pano[:,x0:x1],cmap='gray',vmin=0,vmax=255,aspect='auto',extent=[x0,x1,80,0])
                ax.set_ylim(90,-18)
                for rec in local:
                    a,b=rec['inicio'],rec['fim']
                    if b<=x0 or a>=x1:continue
                    col='#168b9a' if rec['incluida'] else '#c44832'
                    ax.axvline(a,color=col,lw=1);ax.axvline(b,color=col,lw=1)
                    ax.text(max(a,x0),-4,rec['event_id'].split('_')[-1]+f' ({b-a} px)',color=col,fontsize=8)
                ax.set_xlim(x0,x1)
                ax.set_xlabel('Posição no panorama (px)');ax.set_ylabel('Linha')
            fig.suptitle(f'{r.run} | {LABELS[r.pattern]} | Propostas de imagem; nenhuma WCO consultada')
            fig.savefig(figs/f'{r.run}.png',dpi=120);plt.close(fig)
    cat=pd.DataFrame(records)
    # Conta vizinhos e sobreposição somente depois de conhecer todos os eventos aceitos.
    cat['outros_eventos_no_contexto']=0;cat['segmentos_sobrepostos']=0
    for run,g in cat[cat.incluida].groupby('aquisicao'):
        for i,r in g.iterrows():
            cat.loc[i,'outros_eventos_no_contexto']=int(((g.inicio<r.segmento_fim)&(g.fim>r.segmento_inicio)).sum()-1)
            cat.loc[i,'segmentos_sobrepostos']=int(((g.segmento_inicio<r.segmento_fim)&(g.segmento_fim>r.segmento_inicio)).sum()-1)
    cat.to_csv(out/'catalogo_eventos.csv',index=False)
    pd.DataFrame(fixed).to_csv(out/'segmentos_fixos.csv',index=False)
    pd.DataFrame(observability).to_csv(out/'observabilidade.csv',index=False)
    (out/'quadros_catalogo_sha256.json').write_text(json.dumps(frames_hashes,indent=2)+'\n')
    if not reviewpath.exists():
        print(cat[cat.aquisicao.isin(OBSERVABLE)][['event_id','inicio','fim','extensao','motivo']].to_string(index=False));return
    eligible=cat[cat.incluida];allpairs=[]
    for rx,ry in itertools.combinations(sorted(eligible.aquisicao.unique()),2):
        gx=eligible[eligible.aquisicao==rx];gy=eligible[eligible.aquisicao==ry]
        for x in gx.itertuples(index=False):
            for y in gy.itertuples(index=False):
                ratio=max(x.extensao,y.extensao)/min(x.extensao,y.extensao)
                relation='equivalentes' if ratio<=1.25 else 'distintos' if ratio>=1.5 else 'intermediarios'
                scope='intra' if x.padrao==y.padrao else 'entre'
                group=f'{scope}_{relation}';pid=f'{x.event_id}__{y.event_id}'
                priority=hashlib.sha256(f'{CFG["semente"]}|{pid}|{group}'.encode()).hexdigest()
                allpairs.append(dict(pair_id=pid,event_x=x.event_id,event_y=y.event_id,
                    run_x=rx,run_y=ry,pattern_x=x.padrao,pattern_y=y.padrao,
                    grupo=group,razao_extensoes=ratio,prioridade=priority))
    pairs=pd.DataFrame(allpairs);pairs['selecionado']=False
    for _,g in pairs.groupby(['run_x','run_y','grupo']):
        if g.grupo.iloc[0].endswith('intermediarios'):continue
        ix=g.sort_values('prioridade').head(CFG['max_pares_por_aquisicoes_e_grupo']).index
        pairs.loc[ix,'selecionado']=True
    pairs.to_csv(out/'pares_candidatos.csv',index=False)
    pairs[pairs.selecionado].to_csv(out/'pares_eventos.csv',index=False)
    cat['parceiros_equivalentes_outra_aquisicao']=0
    for i,r in cat[cat.incluida].iterrows():
        eq=pairs[pairs.grupo.str.endswith('equivalentes')]
        cat.loc[i,'parceiros_equivalentes_outra_aquisicao']=int(((eq.event_x==r.event_id)|(eq.event_y==r.event_id)).sum())
    cat.to_csv(out/'catalogo_eventos.csv',index=False)
    meta=dict(candidatos=len(cat),incluidos=int(cat.incluida.sum()),aquisicoes_com_eventos=int(eligible.aquisicao.nunique()),
              pares_admissiveis=len(pairs),pares_selecionados=int(pairs.selecionado.sum()),
              grupos=pairs[pairs.selecionado].grupo.value_counts().to_dict(),segundos=time.perf_counter()-t0)
    (out/'catalogo_metadados.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(meta,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
