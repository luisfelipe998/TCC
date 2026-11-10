"""Pearson dos perfis fixos e resumos descritivos dos 325 pares da v16."""
import ast
from itertools import combinations
import json
import sys
sys.dont_write_bytecode = True
from common import *
from scipy.stats import pearsonr


def residuals(power, wavelengths):
    design = np.column_stack((np.ones(len(wavelengths)), np.log(wavelengths)))
    y = np.log(power+1e-12)
    coefficients = np.linalg.lstsq(design, y.T, rcond=None)[0].T
    z = y-coefficients@design.T
    t = np.log(wavelengths)-np.log(wavelengths).mean()
    slopes = (y-y.mean(axis=1,keepdims=True))@t/np.dot(t,t)
    independent = y-y.mean(axis=1,keepdims=True)-slopes[:,None]*t
    np.testing.assert_allclose(z,independent,atol=1e-12,rtol=1e-12)
    return z,coefficients


def pearson_matrix(z):
    centered = z-z.mean(axis=1,keepdims=True)
    norms = np.linalg.norm(centered,axis=1,keepdims=True)
    if np.any(norms<=1e-12*np.maximum(1.,np.linalg.norm(z,axis=1,keepdims=True))):
        raise ValueError("Descritor constante ou degenerado.")
    unit = centered/norms
    r = np.clip(unit@unit.T,-1,1)
    np.fill_diagonal(r,1.)
    return r,unit


def describe(values):
    a = np.asarray(values,float)
    assert len(a)>0 and np.isfinite(a).all()
    return dict(n=len(a),media=float(a.mean()),mediana=float(np.median(a)),
                minimo=float(a.min()),maximo=float(a.max()),
                q25=float(np.quantile(a,.25)),q75=float(np.quantile(a,.75)))


def calculate():
    acquisitions = records()
    pair_rows = pairs()
    runs = [r["ensaio"] for r in acquisitions]
    n = len(runs)
    assert n==26 and len(pair_rows)==325
    powers,means=[],[]
    for run in runs:
        with np.load(INTER/"individuais"/f"{run}.npz") as saved:
            np.testing.assert_array_equal(saved["wavelengths"],WAVES)
            powers.append(saved["wps_zero"]);means.append(saved["wps_mean_all"])
    powers,means=np.asarray(powers),np.asarray(means)
    z,coefficients=residuals(powers,WAVES)
    r,unit=pearson_matrix(z)
    np.testing.assert_allclose(r,np.corrcoef(z),atol=1e-13,rtol=1e-13)
    errors=[abs(r[i,j]-pearsonr(z[i],z[j]).statistic) for i,j in combinations(range(n),2)]
    assert max(errors)<1e-12 and np.linalg.eigvalsh(r).min()>-1e-10
    # Verify the exact descriptor convention of v16.3 on the v16 wavelength grid.
    ref=ROOT.parent/"v16.3_ancoras/scripts/statistics_anchors.py"
    nodes=[node for node in ast.parse(ref.read_text()).body if isinstance(node,ast.FunctionDef)
           and node.name in ("residuals","pearson_matrix")]
    assert len(nodes)==2
    scope={"np":np,"EPS":1e-12}
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(ref),"exec"),scope)
    reference_z,_=scope["residuals"](powers,WAVES)
    reference_r,_=scope["pearson_matrix"](reference_z)
    np.testing.assert_allclose(z,reference_z,atol=1e-13)
    np.testing.assert_allclose(r,reference_r,atol=1e-13)
    score=np.full((n,n),np.nan)
    records_by_id={v["ensaio"]:v for v in acquisitions}
    rows=[]
    for pair in pair_rows:
        a,b=pair["run_x"],pair["run_y"]
        i,j=runs.index(a),runs.index(b)
        assert (pair["tipo"]=="intra") == (records_by_id[a]["padrao"]==records_by_id[b]["padrao"])
        score[i,j]=score[j,i]=pair["escore_piso05"]
        rows.append(dict(numero_x=i+1,numero_y=j+1,ensaio_x=a,ensaio_y=b,
                         padrao_x=pair["padrao_x"],padrao_y=pair["padrao_y"],tipo=pair["tipo"],
                         comprimento_comum=pair["comprimento"],pearson_r=float(r[i,j]),
                         escore_piso05=pair["escore_piso05"],escore_piso01=pair["escore_piso01"],
                         escore_piso10=pair["escore_piso10"],wco_maxima=pair["wco_maxima"],
                         wco_media_integral=pair["wco_media_integral"],
                         wco_ponderada_zero=pair["wco_ponderada_zero"],arquivo=pair["arquivo"]))
    assert sum(v["tipo"]=="intra" for v in rows)==87
    write_csv(TABLES/"pares_pearson_wco.csv",rows)
    groups={}
    for row in rows:
        key=tuple(sorted((row["padrao_x"],row["padrao_y"]),key=PATTERNS.index))
        groups.setdefault(key,[]).append(row)
    group_rows=[]
    for (a,b),subset in sorted(groups.items(),key=lambda v:(PATTERNS.index(v[0][0]),PATTERNS.index(v[0][1]))):
        group_rows.append(dict(padrao_x=a,padrao_y=b,tipo="intra" if a==b else "inter",pares=len(subset),
            pearson_media=float(np.mean([v["pearson_r"] for v in subset])),
            pearson_mediana=float(np.median([v["pearson_r"] for v in subset])),
            escore_media=float(np.mean([v["escore_piso05"] for v in subset])),
            escore_mediana=float(np.median([v["escore_piso05"] for v in subset])),
            wco_maxima_media=float(np.mean([v["wco_maxima"] for v in subset]))))
    write_csv(TABLES/"resumo_por_padrao.csv",group_rows)
    aggregate=[]
    for kind in ("intra","inter"):
        subset=[v for v in rows if v["tipo"]==kind]
        same_groups=[g for (a,b),g in groups.items() if (a==b)==(kind=="intra")]
        for metric in ("pearson_r","escore_piso05","escore_piso01","escore_piso10","wco_maxima","wco_ponderada_zero"):
            aggregate.append(dict(tipo=kind,metrica=metric,**describe([v[metric] for v in subset]),
                media_balanceada=float(np.mean([np.mean([v[metric] for v in g]) for g in same_groups])),
                grupos=len(same_groups)))
    write_csv(TABLES/"resumo_intra_inter.csv",aggregate)
    lookup={(v["tipo"],v["metrica"]):v for v in aggregate}
    original=next(v for v in read_csv(OLD/"results/preprocessing_confirm_rankings.csv")
                  if v["variant"]=="baseline" and v["metric"]=="weighted_scale_peak_power05")
    for kind,old_key in (("intra","same"),("inter","different")):
        assert abs(lookup[kind,"escore_piso05"]["media_balanceada"]-float(original[old_key]))<1e-12
    averages=[]
    for i,acquisition in enumerate(acquisitions):
        others=np.arange(n)!=i
        same=np.array([v["padrao"]==acquisition["padrao"] for v in acquisitions]) & others
        different=others & ~same
        averages.append(dict(numero=i+1,ensaio=runs[i],padrao=acquisition["padrao"],
            pearson_media_todos=float(r[i,others].mean()),
            pearson_media_intra=float(r[i,same].mean()) if same.any() else "",
            pearson_media_inter=float(r[i,different].mean()),
            escore_media_todos=float(score[i,others].mean()),
            escore_media_intra=float(score[i,same].mean()) if same.any() else "",
            escore_media_inter=float(score[i,different].mean()),pares_intra=int(same.sum()),pares_inter=int(different.sum())))
    write_csv(TABLES/"medias_por_ensaio.csv",averages)
    for name,matrix in (("matriz_pearson_todos",r),("matriz_escore_wco",score)):
        write_csv(TABLES/f"{name}.csv",[dict(ensaio=run,**{other:float(matrix[i,j]) if np.isfinite(matrix[i,j]) else ""
                  for j,other in enumerate(runs)}) for i,run in enumerate(runs)])
    write_csv(TABLES/"perfis_espectrais.csv",[dict(ensaio=run,padrao=acquisitions[i]["padrao"],
        lambda_px=float(w),wps_zero=float(powers[i,k]),potencia_media_todos_lags=float(means[i,k]),
        residuo=float(z[i,k]),descritor_unitario=float(unit[i,k]))
        for i,run in enumerate(runs) for k,w in enumerate(WAVES)])
    np.savez_compressed(INTER/"descritores.npz",runs=np.asarray(runs),wavelengths=WAVES,power=powers,
                        power_mean_all=means,residuals=z,unit=unit,coefficients=coefficients,r=r,score=score)
    contrasts={metric:dict(por_par=lookup["intra",metric]["media"]-lookup["inter",metric]["media"],
                  balanceado=lookup["intra",metric]["media_balanceada"]-lookup["inter",metric]["media_balanceada"])
               for metric in ("pearson_r","escore_piso05")}
    result=dict(ensaios=26,pares=325,intra=87,inter=238,grupos=group_rows,resumos=aggregate,
                contrastes=contrasts,lambda_px=[50,3000],escalas=200,
                pearson_geral=describe([v["pearson_r"] for v in rows]))
    write_json(TABLES/"resultados.json",result)
    write_json(INTER/"auditoria/estatistica.json",dict(status="PASS",pares_scipy=325,
        maior_erro_scipy=max(errors),menor_autovalor=float(np.linalg.eigvalsh(r).min()),
        referencia_funcao_sha256=sha(ref),identidade_descritor_v163=True,
        medias_balanceadas_conferidas=True,descritores_sha256=sha(INTER/"descritores.npz")))
    print(json.dumps(dict(status="PASS",contrastes=contrasts,resumos=[v for v in aggregate
          if v["metrica"] in ("pearson_r","escore_piso05")]),ensure_ascii=False,indent=2),flush=True)
    return acquisitions,rows,averages,result


def render():
    import figuras
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.ticker import NullLocator, FuncFormatter
    acquisitions=records()
    pair_rows=read_csv(TABLES/"pares_pearson_wco.csv")
    with np.load(INTER/"descritores.npz") as saved:
        arrays={key:saved[key] for key in saved.files}
    cmap=LinearSegmentedColormap.from_list("pearson",[
        (0,"#5e3c99"),(.25,"#b2abd2"),(.5,"#fffdf0"),(.625,"#d6eef1"),
        (.75,"#99d4de"),(.875,"#3b93bd"),(1,"#164b72")],N=256)
    layouts=[]
    for key,name,title,palette,vmin,vmax,ticks,label in [
        ("r","matriz_pearson_todos","Correlação de Pearson · panoramas completos\n26 ensaios · λ de 50 a 3.000 px",cmap,-1,1,[-1,-.5,0,.5,1],"Pearson r"),
        ("score","matriz_escore_wco","Escore de coerência wavelet (WCO)\nMáximos por escala · ponderação pela potência · piso de 5%",plt.get_cmap("magma").copy(),0,1,[0,.25,.5,.75,1],"Escore ponderado"),
    ]:
        palette.set_bad("#dedede")
        fig=plt.figure(figsize=(6.3,7.4))
        ax=fig.add_axes([.13,.24,.81,.70])
        picture=ax.imshow(np.ma.masked_invalid(arrays[key]),vmin=vmin,vmax=vmax,cmap=palette,interpolation="nearest")
        ax.set_xticks(range(26),range(1,27),rotation=90,fontsize=8.5)
        ax.set_yticks(range(26),range(1,27),fontsize=8.5)
        ax.tick_params(length=0,pad=3)
        ax.set_xlabel("Ensaio",labelpad=4);ax.set_ylabel("Ensaio")
        for boundary in (9.5,18.5,21.5,24.5):
            ax.axhline(boundary,color="white",lw=1.4)
            ax.axvline(boundary,color="white",lw=1.4)
        bar=fig.add_axes([.22,.15,.62,.02])
        cb=fig.colorbar(picture,cax=bar,orientation="horizontal",ticks=ticks,format=FuncFormatter(figuras.decimal))
        cb.set_label(label,fontsize=9);cb.ax.tick_params(labelsize=9)
        fig.text(.53,.055,"1–10: churn · 11–19: wavy annular",ha="center",fontsize=9)
        fig.text(.53,.028,"20–22: plug · 23–25: slug · 26: dispersed bubble",ha="center",fontsize=9)
        fig.suptitle(title,y=.995,fontsize=10)
        layouts.append(figuras.save_figure(fig,STATS/f"{name}.png"))

    fig,axes=plt.subplots(2,1,figsize=(6.3,7.1))
    fig.suptitle("Correlações dentro e entre padrões",y=.99,fontsize=11)
    for ax,metric,title,ylabel,limits in zip(axes,("pearson_r","escore_piso05"),
        ("Pearson dos perfis espectrais","Escore de coerência wavelet (WCO)"),
        ("Pearson r","Escore ponderado"),((-1,1),(0,1))):
        for x,kind,color in ((1,"intra","#087f8c"),(2,"inter","#b66b27")):
            values=np.array([float(row[metric]) for row in pair_rows if row["tipo"]==kind])
            rng=np.random.default_rng(2718+x)
            ax.scatter(x+rng.uniform(-.17,.17,len(values)),values,s=13,alpha=.45,color=color,edgecolor="none")
            q1,q3=np.quantile(values,[.25,.75])
            ax.plot([x-.23,x+.23],[values.mean()]*2,color="#182d37",lw=2)
            ax.plot([x,x],[q1,q3],color="#182d37",lw=2)
        ax.set_xticks([1,2],["Mesmo padrão\n87 pares","Padrões diferentes\n238 pares"],fontsize=9)
        ax.set(xlim=(.45,2.55),ylim=limits,ylabel=ylabel)
        ax.set_title(title,loc="left",fontsize=10)
        ax.yaxis.set_major_formatter(FuncFormatter(figuras.decimal))
        ax.grid(axis="y",color="#e6ecef",lw=.6)
    fig.subplots_adjust(left=.14,right=.97,top=.92,bottom=.11,hspace=.43)
    fig.text(.53,.025,"Traço: média · barra vertical: intervalo interquartil",ha="center",fontsize=9)
    layouts.append(figuras.save_figure(fig,STATS/"correlacoes_dentro_entre.png"))

    fig,axes=plt.subplots(5,2,figsize=(6.3,8.5))
    fig.suptitle("Perfis espectrais e descritores de Pearson",y=.99,fontsize=11)
    fig.text(.32,.945,"Espectro de potência\nwavelet (WPS)",ha="center",va="center",fontsize=9)
    fig.text(.78,.945,"Resíduos unitários\npara Pearson",ha="center",va="center",fontsize=9)
    colors=["#d34a33","#2377b4","#328652","#9466ac","#b0851c"]
    for i,(pattern,color) in enumerate(zip(PATTERNS,colors)):
        indices=[j for j,row in enumerate(acquisitions) if row["padrao"]==pattern]
        for col,key,ylabel in ((0,"power","WPS / máximo"),(1,"unit","Resíduo")):
            ax=axes[i,col]
            data=arrays[key][indices].copy()
            if key=="power":data/=data.max(axis=1,keepdims=True)
            for profile in data:ax.plot(WAVES,profile,color=color,alpha=.35,lw=.8)
            ax.plot(WAVES,np.median(data,axis=0),color=color,lw=1.6)
            ax.set_title(pattern+f" · n = {len(indices)}" if col==0 else "",loc="left",fontsize=9)
            ax.set(xscale="log",xlim=(50,3000),ylabel=ylabel)
            ax.set_xticks([50,300,3000],["50","300","3000"])
            ax.xaxis.set_minor_locator(NullLocator())
            ax.yaxis.set_major_formatter(FuncFormatter(figuras.decimal))
            if i==4:ax.set_xlabel("λ (px)")
            ax.grid(color="#edf0f0",lw=.5)
    fig.subplots_adjust(left=.14,right=.97,top=.885,bottom=.065,wspace=.47,hspace=.85)
    layouts.append(figuras.save_figure(fig,STATS/"perfis_comparados.png"))
    write_json(INTER/"apresentacao/estatistica.json",layouts)
    print(f"Estatística: {len(layouts)} figuras geradas.",flush=True)


if __name__ == "__main__":
    calculate()
    render()
