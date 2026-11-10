"""Pearson of fixed zero-lag WPS descriptors, following estatistica5.

Primary: frozen 2048 px anchor contexts for all 25 acquisitions.
Sensitivity control: full panoramas, with the identical 200 wavelengths.
All 300 unordered pairs are retained; no significance test treats scales or
shared-acquisition pairs as independent experimental replicates.
"""
import sys
sys.dont_write_bytecode = True
import argparse
import ast
import json
from itertools import combinations

import numpy as np
from scipy.signal import correlate
from scipy.stats import pearsonr

from common import INTER, ROOT, SCRIPTS, image_path, load_image, sha, write_csv, write_json
from prepare_comparisons import load_catalog, load_individual
from wavelet_cache import get_grid
from audit_method_schuck import independent_prepare, max_error
import wavelet_core as core

OUT = INTER/"estatistica"
DEST = ROOT/"04_estatistica"
DATA = DEST/"dados"
PATTERNS = ["churn", "wavy annular", "plug", "slug"]
COLORS = {"churn":"#d34a33", "wavy annular":"#2377b4", "plug":"#328652", "slug":"#9466ac"}
REFERENCE = ROOT.parent/"estatistica5/todos_contra_todos/run_all_pairs.py"
EPS = 1e-12
FACTOR = 4*np.pi/(6+np.sqrt(38))


def residuals(power, waves):
    power = np.asarray(power, dtype=float)
    if power.ndim != 2 or not np.isfinite(power).all() or np.any(power<0):
        raise ValueError("Perfis invalidos.")
    design = np.column_stack((np.ones(len(waves)), np.log(waves)))
    y = np.log(power+EPS)
    coefficients = np.linalg.lstsq(design, y.T, rcond=None)[0].T
    z = y-coefficients@design.T
    assert np.max(abs(z@design)) < 1e-9
    # Independent scalar regression checks intercept, slope and log axis.
    t = np.log(waves)
    tc = t-t.mean()
    slopes = (y-y.mean(axis=1,keepdims=True))@tc/np.dot(tc,tc)
    independent = y-y.mean(axis=1,keepdims=True)-slopes[:,None]*tc
    np.testing.assert_allclose(z, independent, atol=1e-12, rtol=1e-12)
    return z, coefficients


def pearson_matrix(z):
    z = np.asarray(z, dtype=float)
    if z.ndim != 2 or not np.isfinite(z).all():
        raise ValueError("Descritores invalidos.")
    centered = z-z.mean(axis=1, keepdims=True)
    norm = np.linalg.norm(centered, axis=1, keepdims=True)
    if np.any(norm <= 1e-12*np.maximum(1., np.linalg.norm(z, axis=1, keepdims=True))):
        raise ValueError("Pearson indefinido para descritor constante ou degenerado.")
    u = centered/norm
    r = u@u.T
    assert np.isfinite(r).all() and r.min()>=-1-1e-12 and r.max()<=1+1e-12
    r = np.clip(r,-1,1)  # Floating-point roundoff only, as in estatistica5.
    np.fill_diagonal(r, 1.)
    return r, u


def temporal_zero(rxx, waves):
    lags = np.arange(-(len(rxx)//2),len(rxx)//2+1)
    assert len(lags) == len(rxx)
    out = []
    for scale in waves/FACTOR:
        u = lags/scale
        kernel = np.pi**(-.25)/np.sqrt(scale)*np.exp(-u*u/2)*np.cos(6*u)
        out.append(abs(np.sum(rxx*kernel)))
    return np.asarray(out)


def formulas(waves):
    known = np.array([[-1.,0,1], [1.,0,-1], [1.,-2,1], [6.,13,20]])
    r,_ = pearson_matrix(known)
    np.testing.assert_allclose(r[0], [1,-1,0,1], atol=1e-12)
    try:
        pearson_matrix(np.array([[1,2,3], [1,1,1]]))
    except ValueError:
        pass
    else:
        raise AssertionError("Vetor constante aceito.")
    # Execute only the two pure mathematical functions from the requested
    # reference, without its legacy paths, 26-run dataset, or processing main.
    tree = ast.parse(REFERENCE.read_text())
    funcs = [n for n in tree.body if isinstance(n,ast.FunctionDef)
             and n.name in ("residuals", "pearson_matrix")]
    assert len(funcs) == 2
    namespace = {"np":np, "WAVELENGTHS":waves}
    exec(compile(ast.Module(body=funcs,type_ignores=[]),str(REFERENCE),"exec"),namespace)
    rng = np.random.default_rng(20261025)
    power = np.exp(rng.normal(size=(7,len(waves))))
    z,_ = residuals(power,waves)
    zr,_ = namespace["residuals"](power,waves)
    rr,_ = namespace["pearson_matrix"](zr)
    own,_ = pearson_matrix(z)
    np.testing.assert_allclose(z,zr,atol=1e-13)
    np.testing.assert_allclose(own,rr,atol=1e-13)
    zscale,_ = residuals(power,waves/FACTOR)
    np.testing.assert_allclose(z,zscale,atol=1e-12)
    assert abs(core.FOURIER_FACTOR-FACTOR)<1e-15
    return dict(status="PASS", casos_menos_um_zero_um=True, constante_recusado=True,
                invariancia_afim_positiva=True, igualdade_com_estatistica5=True,
                residuos_equivalentes_escala_lambda=True, referencia_sha256=sha(REFERENCE))


def describe(values):
    a = np.asarray(values, dtype=float)
    assert a.ndim==1 and len(a)>0 and np.isfinite(a).all()
    return dict(n=len(a), minimo=float(a.min()), q25=float(np.quantile(a,.25)),
                mediana=float(np.median(a)), q75=float(np.quantile(a,.75)),
                maximo=float(a.max()), negativos=int((a<0).sum()))


def calculate():
    for folder in (OUT,DEST,DATA):
        folder.mkdir(exist_ok=True)
    catalog = load_catalog()
    waves = get_grid(2048)["wavelengths"]
    controls = formulas(waves)
    records, primary, complete, verification = [], [], [], []
    source_hashes = {}
    for pattern in PATTERNS:
        for run in sorted(r for r,v in catalog["aquisicoes"].items() if v["padrao"]==pattern):
            record, auto = load_individual(catalog,run,"comparavel")
            pano, raw = load_image(run)
            # Independent inclusive ROI extraction, retaining the original pixels.
            np.testing.assert_array_equal(raw,pano[18:72].sum(axis=0,dtype=np.uint64)/54)
            start, stop = record["segmento_inicio"], record["segmento_fim"]
            np.testing.assert_array_equal(raw[start:stop],auto["raw"])
            x = independent_prepare(raw[start:stop])
            error_prepare = max_error(x,auto["x"])
            rxx = np.correlate(x,x,"full")/len(x)
            p = abs(auto["wxx"][:,len(x)-1])
            error_primary = max_error(temporal_zero(rxx,waves),p)
            full_x = core.prepare(raw)
            error_full_prepare = max_error(full_x,independent_prepare(raw))
            full_r = correlate(full_x,full_x,mode="full",method="fft")/len(full_x)
            n = len(full_x)
            direct_errors = [abs(full_r[n-1+m]-np.dot(full_x[m:],full_x[:n-m])/n)
                             for m in (0,1,17,127,n//4,n//2,n-1)]
            assert max(direct_errors)<1e-12
            pf = temporal_zero(full_r,waves)
            # All scales, different numerical implementation (frequency multiplier).
            transformed,_ = core.transform_correlation(full_r,waves)
            error_full = max_error(pf,abs(transformed[:,n-1]))
            del transformed
            native = catalog["aquisicoes"][run]["perfis"]["nativo"]["contexto"]
            records.append(dict(numero=len(records)+1, ensaio=run, padrao=pattern,
                ancora_px=record["ancora"], janela_nativa_px=native["comprimento"],
                janela_comparavel_px=2048, inicio_px=start, fim_exclusivo_px=stop,
                panorama_px=n, cache_individual=auto["key"], panorama_sha256=sha(image_path(run))))
            primary.append(p); complete.append(pf)
            source_hashes[str(image_path(run))] = sha(image_path(run))
            verification.append(dict(ensaio=run, erro_preparo=error_prepare,
                erro_wps_zero=error_primary, erro_preparo_completo=error_full_prepare,
                erro_correlacao_completa=max(direct_errors), erro_wps_zero_completo=error_full))
            print(f"Pearson: entrada e WPS verificadas {len(records)}/25: {run}",flush=True)
    assert len(records)==25
    p, pf = np.asarray(primary), np.asarray(complete)
    z, coef = residuals(p,waves); zf, coeff = residuals(pf,waves)
    r,u = pearson_matrix(z); rf,uf = pearson_matrix(zf)
    validation = []
    for condition, descriptors, matrix in (("ancora_2048",z,r), ("panorama_completo",zf,rf)):
        np.testing.assert_allclose(matrix,np.corrcoef(descriptors),atol=1e-13)
        eig = np.linalg.eigvalsh(matrix)
        assert eig.min()>=-1e-10
        errors = [abs(matrix[i,j]-pearsonr(descriptors[i],descriptors[j]).statistic)
                  for i,j in combinations(range(25),2)]
        assert max(errors)<1e-12
        validation.append(dict(condicao=condition, pares_scipy=300,
            maior_erro_scipy=max(errors), menor_autovalor=float(eig.min()),
            celulas_finitas=int(np.isfinite(matrix).sum()), simetrica=bool(np.allclose(matrix,matrix.T))))
    # Match WCO by unordered acquisition identity, never by array order or label.
    wco_rows = []
    for batch in ("lote_intra", "lote_inter"):
        wco_rows += json.loads((INTER/batch/"metricas_wavelets.json").read_text())
    wco_lookup = {frozenset((v["run_x"],v["run_y"])):v for v in wco_rows}
    assert len(wco_lookup)==len(wco_rows)==300
    pairs=[]
    for i,j in combinations(range(25),2):
        a,b = records[i],records[j]
        wave = wco_lookup[frozenset((a["ensaio"],b["ensaio"]))]
        assert wave["perfil"]=="comparavel" and wave["comprimento"]==2048
        assert wave["escalas_validas"]==200 and wave["lags_resumo"]==65
        assert {wave["pattern_x"],wave["pattern_y"]}=={a["padrao"],b["padrao"]}
        pairs.append(dict(numero_x=i+1,numero_y=j+1,ensaio_x=a["ensaio"],ensaio_y=b["ensaio"],
            padrao_x=a["padrao"],padrao_y=b["padrao"],
            tipo="intra" if a["padrao"]==b["padrao"] else "inter",
            pearson_r=float(r[i,j]),pearson_r_completo=float(rf[i,j]),
            diferenca_absoluta_r=float(abs(r[i,j]-rf[i,j])),
            wco_media_central=wave["wco_media"], mapa_sha256=wave["mapa_sha256"]))
    assert sum(v["tipo"]=="intra" for v in pairs)==87
    aggregate=[]
    for kind in ("intra","inter"):
        subset=[v for v in pairs if v["tipo"]==kind]
        for metric in ("pearson_r","pearson_r_completo","wco_media_central"):
            aggregate.append(dict(tipo=kind,metrica=metric,**describe([v[metric] for v in subset])))
    groups=[]
    for i,px in enumerate(PATTERNS):
        for py in PATTERNS[i:]:
            subset=[v for v in pairs if (v["padrao_x"],v["padrao_y"])==(px,py)]
            groups.append(dict(padrao_x=px,padrao_y=py,pares=len(subset),
                **{f"r_{k}":v for k,v in describe([v["pearson_r"] for v in subset]).items() if k!="n"},
                wco_mediana=float(np.median([v["wco_media_central"] for v in subset]))))
    intra = [v for v in pairs if v["tipo"]=="intra"]
    inter = [v for v in pairs if v["tipo"]=="inter"]
    limits = dict(inter_r_acima_mediana_intra=int(sum(v["pearson_r"]>=np.median([x["pearson_r"] for x in intra]) for v in inter)),
        inter_wco_acima_mediana_intra=int(sum(v["wco_media_central"]>=np.median([x["wco_media_central"] for x in intra]) for v in inter)),
        delta_r=describe([v["diferenca_absoluta_r"] for v in pairs]),
        maiores_r_inter=sorted(inter,key=lambda v:v["pearson_r"],reverse=True)[:5],
        menores_wco_inter=sorted(inter,key=lambda v:v["wco_media_central"])[:5],
        maiores_wco_inter=sorted(inter,key=lambda v:v["wco_media_central"],reverse=True)[:5])
    np.savez_compressed(OUT/"descritores.npz",wavelengths=waves,p=p,p_full=pf,z=z,z_full=zf,
        coefficients=coef,coefficients_full=coeff,unit=u,unit_full=uf,r=r,r_full=rf,
        runs=np.asarray([v["ensaio"] for v in records]))
    write_csv(DATA/"indice_ensaios.csv",records)
    write_json(OUT/"indice_ensaios.json",records)
    write_csv(DATA/"pares_pearson_wco.csv",pairs)
    write_json(OUT/"pares.json",pairs)
    write_csv(DATA/"resumo_intra_inter.csv",aggregate)
    write_csv(DATA/"resumo_por_padrao.csv",groups)
    for name,matrix in (("matriz_pearson_todos",r),("matriz_pearson_completo",rf)):
        write_csv(DATA/(name+".csv"),[dict(ensaio=a["ensaio"],**{b["ensaio"]:float(matrix[i,j])
            for j,b in enumerate(records)}) for i,a in enumerate(records)])
    profiles = [dict(ensaio=a["ensaio"],padrao=a["padrao"],lambda_px=float(w),
        wps_zero=float(p[i,k]),residuo=float(z[i,k]),descritor_unitario=float(u[i,k]),
        wps_zero_completo=float(pf[i,k]),residuo_completo=float(zf[i,k]))
        for i,a in enumerate(records) for k,w in enumerate(waves)]
    write_csv(DATA/"perfis_espectrais.csv",profiles)
    write_csv(OUT/"verificacao_entradas.csv",verification)
    result=dict(ensaios=25,pares=300,intra=87,inter=213,janela=2048,
                lambda_px=[50,1024],escalas=200,resumos=aggregate,grupos=groups,**limits)
    write_json(DATA/"resultados.json",result)
    checks=dict(status="PASS",controles=controls,matrizes=validation,
        erro_maximo_perfil_primario=max(v["erro_wps_zero"] for v in verification),
        erro_maximo_perfil_completo=max(v["erro_wps_zero_completo"] for v in verification),
        preparacao_sha256=catalog["congelamento_sha256"],script_sha256=sha(__file__),
        descritores_sha256=sha(OUT/"descritores.npz"),fontes_sha256=source_hashes,
        numpy=np.__version__)
    write_json(OUT/"verificacao.json",checks)
    return dict(ensaios=25,pares=300,status=checks["status"],resumos=aggregate)


def save(fig,name):
    import plots
    plots.INTER=OUT
    plots.save_figure(fig,DEST/name)


def matrix_figure(matrix,records,name,title):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    cmap=LinearSegmentedColormap.from_list("purple_cream_blue",[
        (0,"#5e3c99"),(.25,"#b2abd2"),(.5,"#fffdf0"),(.625,"#d6eef1"),
        (.75,"#99d4de"),(.875,"#3b93bd"),(1,"#164b72")],N=256)
    fig=plt.figure(figsize=(6.3,6.7))
    ax=fig.add_axes([.16,.16,.77,.73])
    im=ax.imshow(matrix,vmin=-1,vmax=1,cmap=cmap,interpolation="nearest")
    ax.set_xticks(range(25),range(1,26),fontsize=8.5,rotation=90)
    ax.set_yticks(range(25),range(1,26),fontsize=8.5)
    ax.tick_params(length=0,pad=3)
    start=0
    for pattern in PATTERNS:
        count=sum(v["padrao"]==pattern for v in records)
        center=start+(count-1)/2
        ax.text(center,-1.1,pattern,ha="center",va="bottom",fontsize=9)
        ax.text(-2.1,center,pattern,ha="center",va="center",rotation=90,fontsize=9)
        if start:
            ax.axhline(start-.5,color="white",lw=1.6)
            ax.axvline(start-.5,color="white",lw=1.6)
        start+=count
    bar=fig.add_axes([.24,.075,.60,.022])
    cb=fig.colorbar(im,cax=bar,orientation="horizontal",ticks=[-1,-.5,0,.5,1])
    cb.set_label("Pearson r",fontsize=9)
    cb.ax.tick_params(labelsize=9)
    fig.suptitle(title,y=.98,fontsize=10)
    save(fig,name)


def figures():
    import plots
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullLocator
    checked=json.loads((OUT/"verificacao.json").read_text())
    assert checked["status"]=="PASS" and checked["script_sha256"]==sha(__file__)
    assert checked["descritores_sha256"]==sha(OUT/"descritores.npz")
    with np.load(OUT/"descritores.npz") as a:
        arrays={k:a[k] for k in a.files}
    records=json.loads((OUT/"indice_ensaios.json").read_text())
    pairs=json.loads((OUT/"pares.json").read_text())
    matrix_figure(arrays["r"],records,"matriz_pearson_todos",
                  "Matriz de correlação de Pearson\nÂncoras · 2.048 px · λ de 50 a 1.024 px")
    matrix_figure(arrays["r_full"],records,"matriz_pearson_completo",
                  "Matriz de correlação de Pearson\nPanoramas completos · λ de 50 a 1.024 px")
    # Same descriptive comparison as estatistica5, extended with WCO.
    fig,axs=plt.subplots(2,1,figsize=(6.3,7.2))
    fig.suptitle("Pearson e WCO · intra e inter padrão",fontsize=11)
    for ax,metric,label,limits in zip(axs,("pearson_r","wco_media_central"),
        ("Pearson r","WCO média central"),((-1,1),(0,1))):
        for x,kind,color in ((1,"intra","#087f8c"),(2,"inter","#b66b27")):
            values=np.array([v[metric] for v in pairs if v["tipo"]==kind])
            rng=np.random.default_rng(2718+x)
            ax.scatter(x+rng.uniform(-.17,.17,len(values)),values,s=12,alpha=.48,
                       color=color,edgecolor="none")
            q1,median,q3=np.quantile(values,[.25,.5,.75])
            ax.plot([x-.23,x+.23],[median,median],color="#182d37",lw=2)
            ax.plot([x,x],[q1,q3],color="#182d37",lw=2.5)
        ax.set_xticks([1,2],["Mesmo padrão (87)","Padrões diferentes (213)"])
        ax.set_xlim(.45,2.55);ax.set_ylim(limits)
        ax.set_ylabel(label)
        ax.grid(axis="y",color="#e6ecef",lw=.6)
    fig.subplots_adjust(left=.13,right=.97,top=.92,bottom=.07,hspace=.29)
    save(fig,"correlacoes_dentro_entre")
    fig,axs=plt.subplots(4,2,figsize=(6.3,8.6))
    fig.suptitle("Perfis WPS e descritores · 2.048 px",fontsize=11)
    waves=arrays["wavelengths"]
    for i,pattern in enumerate(PATTERNS):
        indices=[j for j,v in enumerate(records) if v["padrao"]==pattern]
        for col,key,ylabel in ((0,"p","WPS / máximo"),(1,"unit","Resíduo unitário")):
            ax=axs[i,col]
            data=arrays[key][indices].copy()
            if key=="p":
                data/=data.max(axis=1,keepdims=True)
            for v in data:
                ax.plot(waves,v,color=COLORS[pattern],alpha=.35,lw=.8)
            ax.plot(waves,np.median(data,axis=0),color=COLORS[pattern],lw=1.6)
            ax.set_title(pattern,loc="left",fontsize=9)
            ax.set_xscale("log");ax.set_xlim(50,1024)
            ax.set_xticks([64,256,1024],["64","256","1024"])
            ax.xaxis.set_minor_locator(NullLocator())
            ax.set_xlabel("λ (px)");ax.set_ylabel(ylabel)
            ax.grid(color="#edf0f0",lw=.5)
    fig.subplots_adjust(left=.13,right=.97,top=.93,bottom=.06,wspace=.43,hspace=.77)
    save(fig,"perfis_comparados")
    fig,ax=plt.subplots(figsize=(6.3,5.9))
    for kind,color,title in (("intra","#087f8c","Mesmo padrão"),("inter","#b66b27","Padrões diferentes")):
        group=[v for v in pairs if v["tipo"]==kind]
        ax.scatter([v["pearson_r"] for v in group],[v["pearson_r_completo"] for v in group],
                   s=18,alpha=.65,color=color,label=title,edgecolor="none")
    ax.plot([-1,1],[-1,1],ls="--",color="#697780",lw=.8)
    ax.set(xlim=(-1,1),ylim=(-1,1),xlabel="Pearson r · âncora de 2.048 px",
           ylabel="Pearson r · panorama completo")
    ax.legend(frameon=False,loc="lower right",fontsize=9)
    ax.grid(color="#edf0f0",lw=.5)
    fig.suptitle("Sensibilidade ao contexto · Pearson",fontsize=11)
    fig.subplots_adjust(left=.14,right=.97,bottom=.11,top=.91)
    save(fig,"sensibilidade_comprimento")
    files=sorted(DEST.glob("*.png"))
    write_json(OUT/"figuras.json",[dict(arquivo=p.name,sha256=sha(p)) for p in files])
    return dict(figuras=len(files))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage",choices=["calculate","figures"])
    args=parser.parse_args()
    print(json.dumps(globals()[args.stage](),ensure_ascii=False,indent=2),flush=True)
