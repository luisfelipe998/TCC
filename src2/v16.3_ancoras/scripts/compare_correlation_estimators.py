"""Isolated pilot figures and known-answer controls for N versus N-|lag|.

No production anchors, images, configuration or reusable spectra are changed.
Display uses shared WPS/WCS limits across estimators and WCO 0..1 with a
distinct over-range color; saved numerical maps are never clipped.
"""
import sys
sys.dont_write_bytecode = True
import json
from datetime import datetime, timezone

import numpy as np

from common import INTER, ROOT, SCRIPTS, sha, write_csv, write_json
from matplotlib.ticker import FuncFormatter, MaxNLocator
from pair_identity import pair_key
from prepare_comparisons import load_catalog, load_individual
from wavelet_core import pair as production_pair
from audit_method_schuck import (FACTOR, direct_correlation, temporal_cwt,
                                independent_pair, independent_prepare,
                                explicit_cwt, max_error, protected_files)

OUT = INTER/"comparacao_estimadores"
METHODS = (("biased","01_enviesado","N"),("unbiased","02_nao_enviesado","N−|lag|"))


def figures_and_maps(catalog):
    import plots
    selection = json.loads((INTER/"selecao.json").read_text())
    anchors = {r["anchor_id"]:r for r in selection["ancoras"]}
    results = []
    original_heatmap, original_save = plots.heatmap, plots.save_figure
    original_inter, original_wavelets, original_cmap = plots.INTER, plots.WAVELETS, plots.COHERENCE_CMAP
    try:
        for comparison in selection["pares"]:
            _,a = load_individual(catalog,comparison["run_x"],"nativo")
            _,b = load_individual(catalog,comparison["run_y"],"nativo")
            x,y = independent_prepare(a["raw"]),independent_prepare(b["raw"])
            n = len(x)
            waves = np.geomspace(50,n/2,200)
            lags = np.arange(1-n,n)
            support = n-1-abs(lags)[None,:] >= np.sqrt(2)*(waves/FACTOR)[:,None]
            all_maps = {}
            for estimator,_,_ in METHODS:
                maps = production_pair(x,y,waves,estimator=estimator)
                reference = independent_pair(x,y,waves,unbiased=estimator=="unbiased")
                for key in ("wxx","wyy","wxy","wco"):
                    max_error(reference[key],maps[key])
                maps.update(x=x,y=y,raw_x=np.array(a["raw"]),raw_y=np.array(b["raw"]),
                            wavelengths=waves,support=support)
                for key,name in (("wxx","power_mean_x"),("wyy","power_mean_y")):
                    maps[name] = np.asarray([np.abs(row)[mask].mean() for row,mask in zip(maps[key],support)])
                all_maps[estimator] = maps
            common_power = max(float(np.abs(m[k]).max()) for m in all_maps.values() for k in ("wxx","wyy"))
            common_cross = max(float(np.abs(m["wxy"]).max()) for m in all_maps.values())
            common_mean = 1.05*max(float(m[k].max()) for m in all_maps.values() for k in ("power_mean_x","power_mean_y"))
            for estimator,directory,divisor in METHODS:
                folder = OUT/directory
                (folder/"dados").mkdir(parents=True,exist_ok=True)
                maps = all_maps[estimator]
                center = maps["wco"][:,n-33:n+32]
                valid = maps["denominator_valid"]
                summary = dict(padrao=comparison["padrao"],estimador=estimator,divisor=divisor,comprimento=n,
                    run_x=comparison["run_x"],run_y=comparison["run_y"],
                    wco_media=float(center[valid].mean()),wco_mediana=float(np.median(center[valid])),
                    wco_maxima=float(np.nanmax(maps["wco"])),
                    wco_maxima_centro=float(np.nanmax(center)),
                    wco_maxima_suporte=float(np.nanmax(np.where(support,maps["wco"],np.nan))),
                    fracao_acima_1_centro=float(np.mean(center[valid]>1+1e-10)),
                    fracao_acima_1_suporte=float(np.mean(maps["wco"][support]>1+1e-10)),
                    escalas_com_auto_nao_positiva_x=int(np.sum(maps["wxx"][:,n-1].real<=0)),
                    escalas_com_auto_nao_positiva_y=int(np.sum(maps["wyy"][:,n-1].real<=0)),
                    limite_wps_comum=common_power,limite_wcs_comum=common_cross,limite_potencia_comum=common_mean)
                metadata = dict(comparacao=comparison,metricas=summary,estudo="sensibilidade do estimador; mesmos quatro pares",
                    selecao_sha256=sha(INTER/"selecao.json"),nucleo_sha256=sha(SCRIPTS/"wavelet_core.py"),
                    escala_visual_wco="0 a 1; azul indica >1; arrays completos sem clipping",
                    data_utc=datetime.now(timezone.utc).isoformat())
                np.savez_compressed(folder/"dados"/f"{pair_key(comparison)}.npz",**maps,
                                    meta_json=np.array(json.dumps(metadata,ensure_ascii=False)))
                plots.INTER = folder/"dados"
                plots.WAVELETS = folder
                plots.COHERENCE_CMAP = original_cmap.copy()
                plots.COHERENCE_CMAP.set_over("#00a7c5")

                def comparison_heatmap(fig,slot,maps,key,title,upper,mean=None,power_upper=None):
                    limit = 1. if key=="wco" else common_cross if key=="wxy" else common_power
                    over = key=="wco" and bool(np.nanmax(maps[key])>1+1e-10)
                    first_axis = len(fig.axes)
                    original_heatmap(fig,slot,maps,key,"WCO · azul: >1" if over else title,
                                     limit,mean,common_mean if mean is not None else None)
                    if key=="wco" and over:
                        ax,bar = fig.axes[first_axis:first_axis+2]
                        picture = ax.collections[0]
                        bar.clear()
                        cbar = fig.colorbar(picture,cax=bar,extend="max",format=FuncFormatter(plots.decimal))
                        cbar.ax.tick_params(labelsize=8.5,length=2,pad=2)
                        cbar.locator = MaxNLocator(4)
                        cbar.update_ticks()

                def comparison_save(fig,stem):
                    n_label = f"{n:,}".replace(",",".")
                    fig._suptitle.set_text(f"{comparison['padrao']} · {divisor} · {comparison['run_x']} × {comparison['run_y']} · {n_label} px")
                    original_save(fig,stem)

                plots.heatmap,plots.save_figure = comparison_heatmap,comparison_save
                plots.pair_figure(comparison,anchors,maps,summary,label=estimator)
                results.append(summary)
                print(f"Figura isolada: {comparison['padrao']} / {estimator}; max centro={summary['wco_maxima_centro']:.3f}",flush=True)
    finally:
        plots.heatmap,plots.save_figure = original_heatmap,original_save
        plots.INTER,plots.WAVELETS,plots.COHERENCE_CMAP = original_inter,original_wavelets,original_cmap
    for estimator,directory,_ in METHODS:
        write_csv(OUT/directory/"dados/metricas.csv",[r for r in results if r["estimador"]==estimator])
    write_csv(OUT/"comparacao_piloto.csv",results)
    return results


def known_answer_controls():
    n = 2048
    waves = np.array([128.])
    x,y = np.zeros(n),np.zeros(n)
    x[600],y[1100] = 1,1
    impulses = []
    for estimator,_,_ in METHODS:
        result = independent_pair(x,y,waves,unbiased=estimator=="unbiased")
        measured = int(np.argmax(result["wco"][0])-(n-1))
        peak = float(np.nanmax(result["wco"]))
        expected = 1. if estimator=="biased" else (n/(n-500))**2
        assert abs(peak-expected)<1e-12 and measured==-500
        impulses.append(dict(estimador=estimator,atraso_real=500,lag_estimado=measured,
                             wco_pico=peak,wco_pico_prevista_algebricamente=expected))
    write_json(OUT/"controle_impulsos.json",impulses)

    rng = np.random.default_rng(20261026)
    axis = np.arange(n)
    base = np.exp(-.5*((axis-600)/50)**2)*np.cos(2*np.pi*(axis-600)/128)
    lags = np.arange(1-n,n)
    supported = n-1-abs(lags)>=np.sqrt(2)*(128/FACTOR)
    trials = []
    for delay in (0,500):
        target = np.exp(-.5*((axis-(600+delay))/50)**2)*np.cos(2*np.pi*(axis-(600+delay))/128)
        for snr in (-10,0,10):
            noise_std = np.sqrt(np.mean(base**2)/10**(snr/10))
            for trial in range(48):
                a,b = [independent_prepare(128+z+noise_std*rng.normal(size=n)) for z in (base,target)]
                for estimator,_,_ in METHODS:
                    rr = [direct_correlation(u,v,unbiased=estimator=="unbiased") for u,v in ((a,a),(b,b),(a,b))]
                    aa,bb = [explicit_cwt(r,waves,[n-1])[0,0] for r in rr[:2]]
                    cross = temporal_cwt(rr[2],waves)[0]
                    coherence = abs(cross)**2/(abs(aa)*abs(bb))
                    peak_at = int(np.argmax(np.where(supported,coherence,-np.inf)))
                    lag = int(lags[peak_at])
                    trials.append(dict(estimador=estimator,atraso_real=delay,snr_db=snr,repeticao=trial,
                                       lag_estimado=lag,erro_absoluto_px=abs(lag+delay),
                                       acerto_16px=abs(lag+delay)<=16,
                                       wco_no_lag_real=float(coherence[n-1-delay]),
                                       wco_maxima_suporte=float(coherence[supported].max())))
            print(f"Controle de transiente: atraso={delay} px, SNR={snr} dB, 48 realizacoes por estimador",flush=True)
    grouped = []
    for estimator,_,_ in METHODS:
        for delay in (0,500):
            for snr in (-10,0,10):
                selected = [r for r in trials if r["estimador"]==estimator and r["atraso_real"]==delay and r["snr_db"]==snr]
                grouped.append(dict(estimador=estimator,atraso_real=delay,snr_db=snr,
                                    erro_mediano_px=float(np.median([r["erro_absoluto_px"] for r in selected])),
                                    taxa_acerto_16px=float(np.mean([r["acerto_16px"] for r in selected])),
                                    fracao_real_acima_1=float(np.mean([r["wco_no_lag_real"]>1+1e-10 for r in selected])),
                                    fracao_maximo_acima_1=float(np.mean([r["wco_maxima_suporte"]>1+1e-10 for r in selected]))))
    write_csv(OUT/"controle_transientes_repeticoes.csv",trials)
    write_csv(OUT/"controle_transientes_resumo.csv",grouped)
    return impulses,grouped


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    before = protected_files()
    catalog = load_catalog()
    pilot = figures_and_maps(catalog)
    impulses,transients = known_answer_controls()
    assert before==protected_files(), "Um artefato de producao foi modificado."
    write_json(OUT/"verificacao.json",dict(status="PASS",figuras_png=8,pares_reais_piloto=4,
        pares_inter_calculados=0,estimadores=[r[0] for r in METHODS],producao_preservada=True,
        criterio="Limite 0..1, identificacao do atraso conhecido, estabilidade do denominador; nao maximizar a WCO.",
        controle_impulsos=impulses,controle_transientes=transients,
        limitacao_controles="Simulacoes com transiente gaussiano modulado e ruido branco; nao validam a classificacao de escoamentos.",
        criterio_visual="Mesmas escalas WPS/WCS entre estimadores em cada par; WCO 0..1 e azul para >1, sem clipping numerico.",
        metricas_piloto=pilot,script_sha256=sha(SCRIPTS/"compare_correlation_estimators.py")))
    print("PASS: oito figuras isoladas; quatro pares originais preservados; nenhum par inter calculado.",flush=True)


if __name__=="__main__":
    main()
