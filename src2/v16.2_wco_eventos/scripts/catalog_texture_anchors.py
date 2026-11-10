"""Catálogo óptico e pares da extensão; não importa wavelets."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib
import itertools
import json
import os
import time
os.environ.setdefault('MPLCONFIGDIR','/tmp/wco-eventos-mpl')
import numpy as np
import pandas as pd
from PIL import Image
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT.parent/'v16_spatial2_18_71'
DATA=ROOT.parent/'data_processed'
OUT=ROOT/'results/extensao_churn_anular'
FIG=ROOT/'figures/extensao_churn_anular'
CFG=json.loads((ROOT/'config_ancoras_textura.json').read_text())
LABELS={'churn':'Agitado (churn)','wavy annular':'Anular ondulado (wavy annular)'}

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def detect(raw,sigma=16):
    sm=gaussian_filter1d(np.asarray(raw,float),sigma,mode='reflect')
    amp=float(np.diff(np.percentile(sm,[5,95]))[0])
    threshold=max(CFG['prominencia_minima_cinza'],CFG['fracao_amplitude_robusta_prominencia']*amp)
    rows=[];margin=CFG['comprimento_segmento']//2+CFG['margem_sensibilidade_px']
    for kind,sign in [('minimo',-1),('maximo',1)]:
        peaks,props=find_peaks(sign*sm,prominence=threshold,width=CFG['largura_meia_prominencia_px'],
            distance=CFG['distancia_candidatos_px'],wlen=CFG['janela_prominencia_px'])
        for j,c in enumerate(peaks):
            eligible=margin<=c<=len(raw)-margin
            stratum=int((c-margin)*CFG['estratos_por_aquisicao']/(len(raw)-2*margin+1)) if eligible else -1
            rows.append(dict(tipo=kind,ancora=int(c),largura=float(props['widths'][j]),
                inicio_meia_prominencia=float(props['left_ips'][j]),fim_meia_prominencia=float(props['right_ips'][j]),
                prominencia=float(props['prominences'][j]),limiar_prominencia=threshold,amplitude_robusta=amp,
                intensidade_suavizada=float(sm[c]),estrato=stratum,elegivel_contexto=eligible))
    columns=['tipo','ancora','largura','inicio_meia_prominencia','fim_meia_prominencia','prominencia',
             'limiar_prominencia','amplitude_robusta','intensidade_suavizada','estrato','elegivel_contexto']
    return pd.DataFrame(rows,columns=columns),sm

def select(candidates):
    selected=set()
    for kind in ['minimo','maximo']:
        positions=[]
        for block in range(CFG['estratos_por_aquisicao']):
            g=candidates[(candidates.tipo==kind)&(candidates.estrato==block)].sort_values(['prominencia','ancora'],ascending=[False,True])
            for idx,r in g.iterrows():
                if all(abs(r.ancora-c)>=CFG['distancia_ancoras_mesmo_tipo_px'] for c in positions):
                    selected.add(idx);positions.append(int(r.ancora));break
    return selected

def frame_at(x,shifts):
    bounds=np.r_[0,np.cumsum(shifts)]
    i=min(int(np.searchsorted(bounds,int(x),side='right')-1),len(shifts)-1)
    return i,int(x-bounds[i])

def counterpart(c,kind,other):
    g=other[(other.tipo==kind)&other.elegivel_contexto].copy()
    if len(g):
        g['distancia']=abs(g.ancora-c);r=g.sort_values(['distancia','ancora']).iloc[0]
        if r.distancia<=CFG['tolerancia_correspondencia_detector_px']:return int(r.ancora)
    return np.nan

def render_run(run,pattern,pano,raw,sm,rows):
    chosen=[r for r in rows if r['selecionada']]
    chosen.sort(key=lambda r:(r['tipo']!='minimo',r['ancora']))
    fig=plt.figure(figsize=(17,3+len(chosen)*1.55),layout='constrained')
    grid=fig.add_gridspec(2+len(chosen),3,height_ratios=[1,1]+[1.25]*len(chosen))
    ax=fig.add_subplot(grid[0,:]);ax.imshow(pano,cmap='gray',vmin=0,vmax=255,aspect='auto',extent=[0,len(raw),80,0])
    ax.set_title('Panorama completo; marcas = âncoras selecionadas por imagem');ax.set_xlabel('Coluna do panorama (px)')
    signalax=fig.add_subplot(grid[1,:]);signalax.plot(raw,color='#b3b3b3',lw=.5,label='Média 18–71 sem filtro adicional')
    signalax.plot(sm,color='#303030',lw=.7,label='Perfil para detecção: σ=16 px')
    for r in chosen:
        col='#137f8b' if r['tipo']=='minimo' else '#b75a2a';c=r['ancora']
        ax.axvline(c,color=col,lw=1);signalax.scatter(c,sm[c],color=col,marker='v' if r['tipo']=='minimo' else '^',s=35)
        signalax.text(c,sm[c],r['anchor_id'].split('_')[-1],fontsize=7)
    signalax.set_ylabel('Intensidade média');signalax.set_xlabel('Coluna do panorama (px)');signalax.legend(fontsize=8)
    for k,r in enumerate(chosen):
        c=r['ancora'];color='#137f8b' if r['tipo']=='minimo' else '#b75a2a'
        a=fig.add_subplot(grid[k+2,0]);a.imshow(pano[:,c-256:c+256],cmap='gray',vmin=0,vmax=255,aspect='auto',extent=[c-256,c+256,80,0])
        a.axvline(c,color=color,lw=1.3);a.set_title(f'{r["anchor_id"]} | {r["tipo"]}',fontsize=9);a.set_xlabel('Detalhe do panorama (px)',fontsize=8)
        frame=np.asarray(Image.open(DATA/run/r['quadro_ancora']).convert('L'));col=r['coluna_quadro_ancora']
        lo=max(0,col-256);hi=min(frame.shape[1],col+256)
        a=fig.add_subplot(grid[k+2,1]);a.imshow(frame[:,lo:hi],cmap='gray',vmin=0,vmax=255,aspect='auto',extent=[lo,hi,80,0])
        a.axvline(col,color=color,lw=1.3);a.set_title(f'Quadro de origem: {r["quadro_ancora"]}',fontsize=9);a.set_xlabel('Coluna no quadro (px)',fontsize=8)
        a=fig.add_subplot(grid[k+2,2]);xx=np.arange(c-256,c+256)
        a.plot(xx,raw[c-256:c+256],color='#aaa',lw=.5);a.plot(xx,sm[c-256:c+256],color=color,lw=1)
        a.axvline(c,color=color);a.axvspan(r['inicio_meia_prominencia'],r['fim_meia_prominencia'],color=color,alpha=.12)
        a.set_title(f'Largura óptica={r["largura"]:.1f} px; proeminência={r["prominencia"]:.1f}'.replace('.',','),fontsize=9)
        a.set_xlabel('Coluna do panorama (px)',fontsize=8)
    fig.suptitle(f'{LABELS[pattern]} | {run}\nÂncoras ópticas; detalhes de 512 px; o recorte wavelet terá 2.048 px',fontsize=13)
    fig.savefig(FIG/'ancoras'/f'{run}.png',dpi=125);plt.close(fig)

def main(render=True):
    t0=time.perf_counter();OUT.mkdir(parents=True,exist_ok=True);(FIG/'ancoras').mkdir(parents=True,exist_ok=True)
    inv=pd.read_csv(BASE/'results/input_audit.csv');inv=inv[inv.included&inv.pattern.isin(CFG['padroes'])].sort_values('run')
    assert len(inv)==19
    records=[];sources={};allraw={}
    for r in inv.itertuples():
        imagefile=DATA/r.run/'__panorama_v11.png';cachefile=BASE/'intermediate'/f'{r.run}_input.npz'
        assert sha(imagefile)==r.panorama_sha256
        pano=np.asarray(Image.open(imagefile).convert('L'));raw=pano[18:72].mean(axis=0);allraw[r.run]=raw
        cache=np.load(cachefile);assert np.array_equal(raw,cache['raw'])
        shifts=cache['shifts'];seams=cache['seams'];nonpos=np.r_[cache['raw_shifts']<=0,cache['raw_shifts'][-1]<=0]
        for p in [imagefile,cachefile]:sources[str(p.relative_to(ROOT.parent))]=sha(p)
        cand,sm=detect(raw,CFG['sigma_deteccao_px']);selected=select(cand)
        alternatives={f'sigma{s}':detect(raw,s)[0] for s in CFG['sigmas_sensibilidade_px']}
        alternatives['faixa30_49']=detect(pano[30:50].mean(axis=0),CFG['sigma_deteccao_px'])[0]
        local=[]
        for j,x in cand.iterrows():
            d=x.to_dict();c=int(x.ancora);s=c-1024;e=c+1024;frame,offset=frame_at(c,shifts)
            file=DATA/r.run/f'ImgA{frame:06d}.jpg';width=Image.open(file).width;col=width//2-int(shifts[frame])//2+offset
            mixed=bool(np.any((seams-c>0)&(seams-c<=2)))
            cf0,_=frame_at(max(0,s),shifts);cf1,_=frame_at(min(r.n,e)-1,shifts)
            d.update(anchor_id=f'{r.run}_{"MIN" if x.tipo=="minimo" else "MAX"}{c:05d}',aquisicao=r.run,padrao=r.pattern,
                selecionada=j in selected,motivo='selecionada_por_imagem' if j in selected else 'sem_contexto' if not x.elegivel_contexto else 'limite_por_estrato_ou_distancia',
                segmento_inicio=s,segmento_fim=e,linhas='18–71',quadro_ancora=file.name,coluna_quadro_ancora=col,
                ancora_em_zona_mistura=mixed,distancia_costura=int(np.min(abs(seams-c))),
                costuras_contexto=int(((seams>max(0,s))&(seams<min(r.n,e))).sum()),
                deslocamentos_nao_positivos_contexto=int(nonpos[cf0:cf1+1].sum()),versao=CFG['versao'])
            if j in selected:
                sources[str(file.relative_to(ROOT.parent))]=sha(file)
                if not mixed:assert np.array_equal(pano[:,c],np.asarray(Image.open(file).convert('L'))[:,col])
            for name,alternative in alternatives.items():
                pos=counterpart(c,x.tipo,alternative);d['ancora_'+name]=pos;d['deslocamento_'+name]=pos-c
            local.append(d);records.append(d)
        if render:render_run(r.run,r.pattern,pano,raw,sm,local)
        print(f'Catálogo {r.run}: {len(cand)} candidatos, {len(selected)} selecionados',flush=True)
    cat=pd.DataFrame(records);cat['vizinhos_selecionados_no_contexto']=0;cat['segmentos_sobrepostos']=0
    for run,g in cat[cat.selecionada].groupby('aquisicao'):
        for i,r in g.iterrows():
            cat.loc[i,'vizinhos_selecionados_no_contexto']=int(((g.ancora>=r.segmento_inicio)&(g.ancora<r.segmento_fim)).sum()-1)
            cat.loc[i,'segmentos_sobrepostos']=int(((g.segmento_inicio<r.segmento_fim)&(g.segmento_fim>r.segmento_inicio)).sum()-1)
    cat.to_csv(OUT/'catalogo_ancoras.csv',index=False)
    chosen=cat[cat.selecionada];candidate_pairs=[]
    for rx,ry in itertools.combinations(sorted(chosen.aquisicao.unique()),2):
        for x,y in itertools.product(chosen[chosen.aquisicao==rx].itertuples(),chosen[chosen.aquisicao==ry].itertuples()):
            ratio=max(x.largura,y.largura)/min(x.largura,y.largura)
            if x.tipo==y.tipo=='maximo':continue
            group=('minimos_compativeis' if ratio<=1.5 else 'minimos_distintos' if ratio>=2 else 'minimos_intermediarios') if x.tipo==y.tipo else 'polaridades_opostas'
            pid=x.anchor_id+'__'+y.anchor_id
            priority=hashlib.sha256(f'{CFG["semente"]}|{group}|{pid}'.encode()).hexdigest()
            candidate_pairs.append(dict(pair_id=pid,anchor_x=x.anchor_id,anchor_y=y.anchor_id,run_x=rx,run_y=ry,
                pattern_x=x.padrao,pattern_y=y.padrao,grupo=group,tipo_x=x.tipo,tipo_y=y.tipo,
                escopo='intra' if x.padrao==y.padrao else 'entre',razao_larguras=ratio,prioridade=priority))
    pairs=pd.DataFrame(candidate_pairs);pairs['selecionado']=False
    for (_,_,group),g in pairs.groupby(['run_x','run_y','grupo']):
        if group=='minimos_intermediarios':continue
        if group=='polaridades_opostas':
            for _,h in g.groupby('tipo_x'):pairs.loc[h.sort_values('prioridade').head(1).index,'selecionado']=True
        else:pairs.loc[g.sort_values('prioridade').head(CFG['max_pares_por_aquisicoes_e_grupo']).index,'selecionado']=True
    pairs.to_csv(OUT/'pares_candidatos.csv',index=False);pairs[pairs.selecionado].to_csv(OUT/'pares_ancoras.csv',index=False)
    coverage=[]
    for r in inv.itertuples():
        g=chosen[chosen.aquisicao==r.run];coverage.append(dict(aquisicao=r.run,padrao=r.pattern,
            candidatos=len(cat[cat.aquisicao==r.run]),minimos=int((g.tipo=='minimo').sum()),maximos=int((g.tipo=='maximo').sum()),
            sigma12_disponiveis=int(g.ancora_sigma12.notna().sum()),sigma20_disponiveis=int(g.ancora_sigma20.notna().sum()),
            nao_positivos_aquisicao=int(r.nonpositive_shifts)))
    pd.DataFrame(coverage).to_csv(OUT/'cobertura_ancoras.csv',index=False)
    (OUT/'fontes_sha256.json').write_text(json.dumps(sources,indent=2)+'\n')
    meta=dict(data_utc=datetime.now(timezone.utc).isoformat(),candidatos=len(cat),ancoras=len(chosen),
        pares_candidatos=len(pairs),pares_selecionados=int(pairs.selecionado.sum()),
        grupos=pairs[pairs.selecionado].groupby(['escopo','grupo']).size().to_dict(),segundos=time.perf_counter()-t0)
    meta['grupos']={' / '.join(k):int(v) for k,v in meta['grupos'].items()}
    (OUT/'catalogo_metadados.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(meta,indent=2,ensure_ascii=False))

if __name__=='__main__':main()
