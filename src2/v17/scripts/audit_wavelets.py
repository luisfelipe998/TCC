"""Independent numerical audit of active v16; all output stays within v17.

Run from any directory with Python 3 and numpy/scipy/Pillow/PyWavelets.
The production module is imported only for comparison; no production main runs.
Exploratory checks are descriptive, use acquisition-level labels, and fit no ML.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parents[1]
REPO = HERE.parents[1]
SOURCE = REPO / "src2/v16_spatial2_18_71"
OUT = HERE / "resultados/wavelets"
OUT.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(HERE / "logs/mpl_wavelets"))
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import scipy
from scipy.fft import fft, ifft, fftfreq, next_fast_len
from scipy.signal import correlate
from scipy.stats import spearmanr
from PIL import Image

WAVES = np.geomspace(50.0, 3000.0, 200)
FACTOR = 4*np.pi/(6+np.sqrt(38))
SEED = 27012026


def read_json(path):
    return json.loads(path.read_text())


def save_json(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+"\n")


def read_csv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def save_csv(name, rows):
    with (OUT / name).open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def independent_prepare(raw):
    # Scalar closed-form OLS rather than scipy.signal.detrend.
    raw = np.asarray(raw, dtype=float)
    t = np.arange(len(raw), dtype=float)
    t -= t.mean()
    x = 255-raw
    x -= x.mean()
    x -= t*np.dot(x,t)/np.dot(t,t)
    x -= x.mean()
    return x/np.sqrt(np.dot(x,x)/len(x))


def zero_quadrature(raw, wavelengths):
    x = independent_prepare(raw)
    r = correlate(x, x, mode="full", method="fft")/len(x)
    m = np.arange(1-len(x),len(x))
    values = []
    for a in wavelengths/FACTOR:
        u = m/a
        values.append(abs(np.dot(r,np.pi**(-.25)/np.sqrt(a)*np.exp(-u*u/2)*np.cos(6*u))))
    return np.array(values)


def independent_residuals(powers, wavelengths):
    t = np.log(wavelengths)
    t -= t.mean()
    y = np.log(powers+1e-12)
    y -= y.mean(axis=1,keepdims=True)
    return y-(y@t/np.dot(t,t))[:,None]*t


def aggregates(matrix, labels):
    groups = {}
    within, between = [], []
    for i in range(len(labels)):
        for j in range(i+1,len(labels)):
            key = tuple(sorted((labels[i],labels[j])))
            groups.setdefault(key,[]).append(float(matrix[i,j]))
            (within if labels[i]==labels[j] else between).append(float(matrix[i,j]))
    a = [np.mean(v) for (x,y),v in groups.items() if x==y]
    b = [np.mean(v) for (x,y),v in groups.items() if x!=y]
    return dict(intra=float(np.mean(within)), inter=float(np.mean(between)),
                delta=float(np.mean(within)-np.mean(between)),
                intra_balanceado=float(np.mean(a)),inter_balanceado=float(np.mean(b)),
                delta_balanceado=float(np.mean(a)-np.mean(b)))


def direct_spectral_pair(x, y, wavelengths, keep_maps=False):
    """Alternative computation without constructing the correlation sequence.

    W_Rxy(a,m) = IFFT(H_a FFT(x) conj(FFT(y)))[m]/N, with zero padding.
    H_a>=0 gives the Cauchy-Schwarz bound for the biased estimator.
    """
    n = len(x)
    scales = wavelengths/FACTOR
    nfft = next_fast_len(2*n-1+2*int(np.ceil(8*scales.max())))
    omega = 2*np.pi*fftfreq(nfft)
    fx,fy = fft(x,nfft),fft(y,nfft)
    xy = fx*fy.conj()
    lags = np.arange(1-n,n)
    index = lags % nfft
    zx,zy,peak,zero,mean,peaklag = [],[],[],[],[],[]
    co_maps = []
    for a in scales:
        h = np.sqrt(2*np.pi*a)*np.pi**(-.25)*np.exp(-.5*(a*omega-6)**2)
        px = np.sum(h*abs(fx)**2)/nfft/n
        py = np.sum(h*abs(fy)**2)/nfft/n
        wc = ifft(h*xy)[index]/n
        co = abs(wc)**2/(px*py)
        zx.append(px);zy.append(py)
        peak.append(co.max());zero.append(co[n-1]);mean.append(co.mean())
        peaklag.append(lags[np.argmax(co)])
        if keep_maps: co_maps.append(co)
    zx,zy,peak,zero,mean = map(np.asarray,(zx,zy,peak,zero,mean))
    weight = np.sqrt(zx*zy)
    active = (zx>=.05*zx.max()) & (zy>=.05*zy.max())
    weighted = weight*active
    result = dict(zx=zx,zy=zy,peak=peak,zero=zero,mean=mean,
                  peaklag=np.asarray(peaklag), score=float(np.dot(weighted,peak)/weighted.sum()),
                  weighted_zero=float(np.dot(weight,zero)/weight.sum()),
                  full_mean=float(mean.mean()),max=float(peak.max()),lags=lags)
    if keep_maps: result["wco"] = np.asarray(co_maps)
    return result


def equation_checks(core):
    rng = np.random.default_rng(SEED)
    x,y = rng.normal(size=(2,257))
    m = np.arange(1-len(x),len(x))
    correlation_errors = {}
    for mode in ["biased","unbiased"]:
        expected = np.correlate(x,y,"full")/(len(x) if mode=="biased" else len(x)-abs(m))
        correlation_errors[mode] = float(np.max(abs(expected-core.correlation(x,y,mode))))
    waves = np.array([50.,125.,512.,3000.])
    r = np.correlate(x,y,"full")/len(x)
    actual = core.morlet_cwt(r,waves/FACTOR)
    quadrature_errors = []
    for k,a in enumerate(waves/FACTOR):
        for b in [0,1,77,256,410,512]:
            u = (np.arange(len(r))-b)/a
            expected = np.dot(r,np.pi**(-.25)/np.sqrt(a)*np.exp(-u*u/2-6j*u))
            quadrature_errors.append(abs(actual[k,b]-expected))
    forward = core.pair(x,y,waves)
    reverse = core.pair(y,x,waves)
    independent = direct_spectral_pair(x,y,waves,keep_maps=True)
    auto = core.pair(x,x,waves)
    direct_error = float(np.max(abs(independent["wco"]-forward["wco"])))
    swap_error = float(np.max(abs(forward["wco"]-reverse["wco"][:,::-1])))
    identity_error = float(np.max(abs(auto["wco"][:,256]-1)))
    sx = np.zeros(2048);sx[500:756] = rng.normal(size=256)
    sy = np.r_[np.zeros(100),sx[:-100]]
    delayed = direct_spectral_pair(sx,sy,np.array([50.,100.,300.]))
    # Padding must remove wraparound, but does not create new observations.
    padded = core.morlet_cwt(r,waves/FACTOR,padding_sigma=16)
    padding_error = float(np.max(abs(actual-padded)))
    result = dict(correlation_max_abs_error=correlation_errors,
        cwt_complex_quadrature_max_abs_error=float(max(quadrature_errors)),
        independent_spectral_identity_max_abs_error=direct_error,
        swap_reflection_max_abs_error=swap_error, self_zero_max_abs_error=identity_error,
        known_delay_peak_lags=delayed["peaklag"].tolist(),
        known_delay_peak_errors=(delayed["peak"]-1).tolist(),
        double_padding_max_abs_error=padding_error,
        fourier_factor=float(FACTOR),morlet_dc_ratio=float(np.exp(-18)))
    assert max(correlation_errors.values())<1e-12
    assert result["cwt_complex_quadrature_max_abs_error"]<1e-12
    assert max(direct_error,swap_error,identity_error,padding_error)<1e-12
    assert delayed["peaklag"].tolist()==[-100,-100,-100]
    save_json("equacoes.json",result)
    return result


def data_and_profiles(core, records):
    inventory = read_csv(SOURCE/"intermediario/auditoria/entradas.csv")
    inputs = [];raw_by_run = {};profile_rows = [];powers = []
    for row in inventory:
        run = row["ensaio"]
        image_path = REPO/"src2/data_processed"/run/"__panorama_v11.png"
        with Image.open(image_path) as image:
            pixels = np.asarray(image.convert("L"))
        raw = pixels[18:72].astype(float).sum(axis=0)/54
        with np.load(SOURCE/"intermediario/entradas"/f"{run}.npz") as saved:
            stored = saved["raw"]
        raw_by_run[run]=raw
        digest = hashlib.sha256(raw.astype("<f8").tobytes()).hexdigest()
        inputs.append(dict(ensaio=run,n=len(raw),raw_error=float(np.max(abs(raw-stored))),
            raw_hash_match=digest==row["raw_sha256"],
            image_hash_match=hashlib.file_digest(image_path.open("rb"),"sha256").hexdigest()==row["panorama_sha256"]))
    assert all(r["raw_error"]==0 and r["raw_hash_match"] and r["image_hash_match"] for r in inputs)
    for row in records:
        run=row["ensaio"];raw=raw_by_run[run]
        x = independent_prepare(raw)
        p = zero_quadrature(raw,WAVES)
        with np.load(SOURCE/"intermediario/individuais"/f"{run}.npz") as saved:
            baseline=saved["wps_zero"]
            prepare_error=float(np.max(abs(x-saved["prepared"])))
        profile_rows.append(dict(ensaio=run,padrao=row["padrao"],n=len(raw),
            prepare_error=prepare_error,power_abs_error=float(np.max(abs(p-baseline))),
            power_relative_error=float(np.max(abs((p-baseline)/baseline))),
            cycles_at_lambda3000=len(raw)/3000,cycles_at_lambda1024=len(raw)/1024,
            mean_x=float(x.mean()),std_x=float(x.std())))
        powers.append(p)
    powers=np.array(powers)
    z=independent_residuals(powers,WAVES)
    r=np.corrcoef(z)
    with np.load(SOURCE/"intermediario/descritores.npz") as saved:
        r_error = float(np.max(abs(r-saved["r"])))
        z_error = float(np.max(abs(z-saved["residuals"])))
    assert r_error<1e-10 and max(v["power_abs_error"] for v in profile_rows)<1e-10
    save_csv("entradas.csv",inputs)
    save_csv("perfis_reproduzidos.csv",profile_rows)
    save_json("descritores_reproduzidos.json",dict(pearson_max_abs_error=r_error,
        residuals_max_abs_error=z_error,summary=aggregates(r,[v["padrao"] for v in records])))
    return raw_by_run


def actual_pair_checks(core,raw_by_run,records):
    pairs = read_csv(SOURCE/"04_estatistica/dados/pares_pearson_wco.csv")
    pairs.sort(key=lambda r:(int(r["comprimento_comum"]),r["ensaio_x"],r["ensaio_y"]))
    chosen = [pairs[0],pairs[len(pairs)//2],pairs[-1]]
    rows=[]
    for row in chosen:
        a,b=row["ensaio_x"],row["ensaio_y"]
        n=int(row["comprimento_comum"])
        x,y=[independent_prepare(raw_by_run[v][:n]) for v in (a,b)]
        result=direct_spectral_pair(x,y,WAVES,keep_maps=True)
        # All-lag independent result checked against saved exact summaries;
        # stored maps retain only the figure sampling positions.
        with np.load(SOURCE/"intermediario/pares"/f"{a}__{b}.npz") as saved:
            indices=saved["lags"].astype(int)+n-1
            map_error=float(np.max(abs(result["wco"][:,indices]-saved["wco"])))
            zero_errors=[float(np.max(abs(result[k]-saved[field]))) for k,field in
                [("zx","wps_zero_x"),("zy","wps_zero_y"),("peak","scale_peak_wco")]]
        errors=[abs(result["score"]-float(row["escore_piso05"])),
                abs(result["max"]-float(row["wco_maxima"])),
                abs(result["weighted_zero"]-float(row["wco_ponderada_zero"])),
                abs(result["full_mean"]-float(row["wco_media_integral"]))]
        # Scale-specific boundary screening is a sensitivity, not a validated COI.
        overlap_scores=[]
        w=np.sqrt(result["zx"]*result["zy"])
        active=(result["zx"]>=.05*result["zx"].max())&(result["zy"]>=.05*result["zy"].max())
        for ratio in [.25,.5]:
            peak=result["wco"][:,abs(result["lags"])<=ratio*n].max(axis=1)
            overlap_scores.append(float(np.sum(w[active]*peak[active])/w[active].sum()))
        rows.append(dict(ensaio_x=a,ensaio_y=b,n=n,map_abs_error=map_error,
            zero_or_peak_abs_error=max(zero_errors),summary_abs_error=max(errors),
            score_full=result["score"],score_abs_lag_le_N4=overlap_scores[0],
            score_abs_lag_le_N2=overlap_scores[1]))
        assert max(map_error,*zero_errors,*errors)<1e-10
        print(f"Pair reproduced {a} {b} n={n}",flush=True)
    save_csv("pares_reproduzidos.csv",rows)
    by_run={r["ensaio"]:r for r in records}
    pair_geometry=[]
    for row in pairs:
        a,b=row["ensaio_x"],row["ensaio_y"]
        n=int(row["comprimento_comum"])
        with np.load(SOURCE/"intermediario/pares"/f"{a}__{b}.npz") as saved:
            zx,zy,peak=[saved[k] for k in ("wps_zero_x","wps_zero_y","scale_peak_wco")]
        active=(zx>=.05*zx.max())&(zy>=.05*zy.max())
        w=np.sqrt(zx*zy)*active
        score=np.dot(w,peak)/w.sum()
        assert abs(score-float(row["escore_piso05"]))<1e-12
        pair_geometry.append(dict(ensaio_x=a,ensaio_y=b,tipo=row["tipo"],n=n,
            padrao_x=by_run[a]["padrao"],padrao_y=by_run[b]["padrao"],
            score=float(score),max=float(row["wco_maxima"]),pearson=float(row["pearson_r"]),
            weighted_fraction_lt3_cycles=float(w[WAVES>n/3].sum()/w.sum()),
            weighted_fraction_lt5_cycles=float(w[WAVES>n/5].sum()/w.sum()),
            scale_at_max=float(WAVES[np.argmax(peak)]),
            contribution_lt3_cycles=float(np.dot(w[WAVES>n/3],peak[WAVES>n/3])/np.dot(w,peak))))
    save_csv("suporte_pares.csv",pair_geometry)
    association={}
    for key in ["score","max","pearson"]:
        association[key]=float(spearmanr([r["n"] for r in pair_geometry],[r[key] for r in pair_geometry]).statistic)
    save_json("associacao_comprimento.json",dict(spearman_descritivo=association,
        note="Correlações entre pares dependentes; nenhum valor-p é utilizado. Comprimento e padrão são confundidos.",
        fraction_weight_lt3_cycles_median=float(np.median([r["weighted_fraction_lt3_cycles"] for r in pair_geometry])),
        fraction_weight_lt5_cycles_median=float(np.median([r["weighted_fraction_lt5_cycles"] for r in pair_geometry])),
        pairs_peak_at_max_wavelength=int(sum(r["scale_at_max"]==WAVES[-1] for r in pair_geometry)),
        pairs_lt3cycles_weight_over_half=sum(r["weighted_fraction_lt3_cycles"]>.5 for r in pair_geometry)))


def crop_sensitivity(raw_by_run,records):
    lengths={r["ensaio"]:len(raw_by_run[r["ensaio"]]) for r in records}
    n=min(lengths.values())
    labels=[r["padrao"] for r in records]
    configs={"full_50_3000":(None,WAVES),"full_50_1024":(None,np.geomspace(50,1024,200))}
    for position in ["prefix","center","suffix"]:
        for end in [1024,3000]:configs[f"{position}_{n}_50_{end}"]=(position,np.geomspace(50,end,200))
    rows=[];within_run=[];store={};slope_rows=[]
    for name,(position,waves) in configs.items():
        powers=[]
        for record in records:
            raw=raw_by_run[record["ensaio"]]
            start=0 if position in [None,"prefix"] else (len(raw)-n)//2 if position=="center" else len(raw)-n
            fragment=raw if position is None else raw[start:start+n]
            powers.append(zero_quadrature(fragment,waves))
        powers=np.asarray(powers)
        z=independent_residuals(powers,waves)
        t=np.log(waves);t-=t.mean()
        logp=np.log(powers+1e-12)
        slopes=(logp-logp.mean(axis=1,keepdims=True))@t/np.dot(t,t)
        for record,slope in zip(records,slopes):
            slope_rows.append(dict(config=name,ensaio=record["ensaio"],padrao=record["padrao"],slope=float(slope)))
        matrix=np.corrcoef(z)
        store[name]=z
        rows.append(dict(config=name,window_n="variable" if position is None else n,
            wavelength_min=float(waves[0]),wavelength_max=float(waves[-1]),**aggregates(matrix,labels)))
        np.savez_compressed(OUT/f"{name}.npz",runs=np.array(list(lengths)),labels=np.array(labels),
                            wavelengths=waves,power=powers,residuals=z,pearson=matrix)
        print("Crop sensitivity",name,rows[-1]["delta"],flush=True)
    for end in [1024,3000]:
        for i,record in enumerate(records):
            for a,b in [("prefix","center"),("prefix","suffix"),("center","suffix")]:
                r=float(np.corrcoef(store[f"{a}_{n}_50_{end}"][i],store[f"{b}_{n}_50_{end}"][i])[0,1])
                within_run.append(dict(ensaio=record["ensaio"],padrao=record["padrao"],n=n,
                    lambda_max=end,window_a=a,window_b=b,pearson=r,
                    identical_full_window=lengths[record["ensaio"]]==n))
    save_csv("sensibilidade_janelas.csv",rows)
    save_csv("repetibilidade_interna_janelas.csv",within_run)
    save_csv("inclinacoes_espectrais.csv",slope_rows)
    save_json("janelas_metodo.json",dict(common_n=n,placements=["prefix","center","suffix"],
        grid_points=200,bands=[[50,1024],[50,3000]],
        warning="Janelas da mesma aquisição são dependentes e frequentemente sobrepostas; não são réplicas independentes. O menor ensaio tem três janelas idênticas. Controle exploratório, sem seleção da melhor configuração."))


def white_noise_controls():
    rng=np.random.default_rng(SEED+1)
    fixed=np.array([50.,125.,250.,500.,750.,1024.,1500.,3000.])
    ratios=np.array([.02,.04,.08,.16,.32,.64,.95])
    rows=[]
    for n in [3171,8192,23079]:
        for replicate in range(48):
            x,y=[independent_prepare(v) for v in rng.normal(size=(2,n))]
            for band,waves in [("fixed_pixels",fixed),("fixed_lambda_over_n",ratios*n)]:
                result=direct_spectral_pair(x,y,waves)
                for k,w in enumerate(waves):
                    rows.append(dict(n=n,replicate=replicate,band=band,lambda_px=float(w),
                        lambda_over_n=float(w/n),peak=float(result["peak"][k]),
                        zero=float(result["zero"][k]),mean=float(result["mean"][k])))
        print("White-noise controls n",n,flush=True)
    save_csv("nulos_brancos.csv",rows)
    summaries=[]
    for n in [3171,8192,23079]:
        for band in ["fixed_pixels","fixed_lambda_over_n"]:
            selected=[r for r in rows if r["n"]==n and r["band"]==band]
            for w in sorted({r["lambda_px"] for r in selected}):
                group=[r for r in selected if r["lambda_px"]==w]
                peak=np.array([r["peak"] for r in group])
                summaries.append(dict(n=n,band=band,lambda_px=w,lambda_over_n=w/n,
                    replicates=len(group),peak_mean=float(peak.mean()),
                    peak_median=float(np.median(peak)),peak_q05=float(np.quantile(peak,.05)),
                    peak_q95=float(np.quantile(peak,.95)),peak_gt09_fraction=float(np.mean(peak>.9)),
                    zero_mean=float(np.mean([r["zero"] for r in group]))))
    save_csv("resumo_nulos_brancos.csv",summaries)
    save_json("nulos_metodo.json",dict(seed=SEED+1,replicates_per_n=48,
        independence="Duas sequências gaussianas independentes por replicação; mesmo par usado nas duas grades.",
        estimator="Divisor N, Morlet omega0=6, preparo com remoção de reta e variância unitária.",
        scope="Teste de mecanismo, não modelo nulo calibrado dos panoramas. Não fornece valores-p para pares reais.",
        fixed_grid=fixed.tolist(),relative_grid=ratios.tolist()))


def main():
    source_path=SOURCE/"scripts/wavelet_core.py"
    spec=importlib.util.spec_from_file_location("v17_audit_wavelet_core",source_path)
    core=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(core)
    records=read_json(SOURCE/"intermediario/aquisicoes.json")
    eq=equation_checks(core)
    raw_by_run=data_and_profiles(core,records)
    actual_pair_checks(core,raw_by_run,records)
    crop_sensitivity(raw_by_run,records)
    white_noise_controls()
    save_json("execucao.json",dict(status="PASS",python=platform.python_version(),
        numpy=np.__version__,scipy=scipy.__version__,seed=SEED,
        original_core_sha256=hashlib.file_digest(source_path.open("rb"),"sha256").hexdigest(),
        audit_script_sha256=hashlib.file_digest(Path(__file__).open("rb"),"sha256").hexdigest(),
        inputs_checked=31,profiles_checked=26,pair_summaries_checked=325,
        full_pair_maps_independently_recomputed=3,synthetic_replicates=144,
        all_outputs=str(OUT.relative_to(REPO))))
    print("PASS wavelet audit",flush=True)


if __name__=="__main__":
    main()
