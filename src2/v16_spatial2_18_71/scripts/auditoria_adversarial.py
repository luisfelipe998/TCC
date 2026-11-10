"""Auditoria independente das entradas, equações e sensibilidade ao ensaio m03q40.

Não substitui panoramas, rótulos, coeficientes ou o conjunto de 26 aquisições.
Reconstrói os 31 panoramas em memória a partir dos JPEGs, reestima deslocamentos
e verifica a extração. Controles estatísticos são explicitamente exploratórios.
"""
import argparse
from io import BytesIO
from itertools import combinations
import json
import sys
sys.dont_write_bytecode = True
from common import *
from scipy.fft import rfft, irfft, next_fast_len
from scipy.signal import correlate
import openpyxl
import wavelet_core as core

TARGET = "m03q40_001"
OUT = INTER / "auditoria" / "adversarial"


def direct_zero(raw, wavelengths=WAVES, linear=True):
    """CWT no lag zero por quadratura temporal, sem chamar a CWT de produção."""
    x = core.prepare(raw, linear=linear)
    r = correlate(x, x, mode="full", method="fft") / len(x)
    lag = np.arange(1-len(x), len(x))
    result = []
    for a in wavelengths / core.FOURIER_FACTOR:
        u = lag/a
        kernel = np.pi**(-.25)/np.sqrt(a)*np.exp(-u*u/2)*np.cos(6*u)
        result.append(abs(np.dot(r, kernel)))
    return np.array(result)


def descriptor_matrix(power, wavelengths):
    """Regressão escalar independente, seguida de Pearson por produto escalar."""
    t = np.log(wavelengths)
    t -= t.mean()
    y = np.log(power+1e-12)
    y -= y.mean(axis=1, keepdims=True)
    residual = y-(y@t / np.dot(t,t))[:,None]*t
    assert (np.linalg.norm(residual,axis=1)>1e-10).all()
    return np.corrcoef(residual)


def aggregates(matrix, acquisitions, excluded=()):
    groups, within, between = {}, [], []
    for i,j in combinations(range(len(acquisitions)),2):
        a,b = acquisitions[i],acquisitions[j]
        if a["ensaio"] in excluded or b["ensaio"] in excluded:
            continue
        key = tuple(sorted((a["padrao"],b["padrao"])))
        value = float(matrix[i,j])
        groups.setdefault(key,[]).append(value)
        (within if key[0]==key[1] else between).append(value)
    wi = [float(np.mean(v)) for (a,b),v in groups.items() if a==b]
    be = [float(np.mean(v)) for (a,b),v in groups.items() if a!=b]
    return dict(ensaios=len(acquisitions)-len(excluded),pares_intra=len(within),
        pares_inter=len(between),intra=float(np.mean(within)),inter=float(np.mean(between)),
        contraste=float(np.mean(within)-np.mean(between)),
        intra_balanceado=float(np.mean(wi)),inter_balanceado=float(np.mean(be)),
        contraste_balanceado=float(np.mean(wi)-np.mean(be)),
        wavy_annular_intra=float(np.mean(groups['wavy annular','wavy annular'])))


def estimate_shifts(profiles):
    """Pearson por sobreposição: FFT dos produtos e somas cumulativas por lag.

    Implementação independente do laço de np.corrcoef do pipeline v11.
    O sinal é o lag a[n+s] b[n], e não o deslocamento físico para a direita.
    """
    a,b = profiles[:-1],profiles[1:]
    width = profiles.shape[1]
    limit = min(200,width//4)
    lags = np.arange(-limit,limit+1)
    n = width-abs(lags)
    nfft = next_fast_len(2*width-1)
    f = rfft(profiles,n=nfft,axis=1)
    products = irfft(f[:-1]*f[1:].conj(),n=nfft,axis=1)[:,lags % nfft]
    starts_a,stops_a = np.maximum(lags,0),width+np.minimum(lags,0)
    starts_b,stops_b = np.maximum(-lags,0),width-np.maximum(lags,0)
    def moments(x,starts,stops):
        sums = np.c_[np.zeros(len(x)),np.cumsum(x,axis=1)]
        squares = np.c_[np.zeros(len(x)),np.cumsum(x*x,axis=1)]
        return sums[:,stops]-sums[:,starts],squares[:,stops]-squares[:,starts]
    sa,qa = moments(a,starts_a,stops_a)
    sb,qb = moments(b,starts_b,stops_b)
    den = np.sqrt(np.maximum(qa-sa*sa/n,0)*np.maximum(qb-sb*sb/n,0))
    scores = np.divide(products-sa*sb/n,den,out=np.full_like(products,-np.inf),where=den>0)
    best = np.argmax(scores,axis=1)
    return lags[best],scores[np.arange(len(best)),best]


def literal_shift(a,b):
    width = len(a)
    candidates = np.arange(-min(200,width//4),min(200,width//4)+1)
    values = []
    for lag in candidates:
        aa = a[max(lag,0):width+min(lag,0)]
        bb = b[max(-lag,0):width-max(lag,0)]
        values.append(np.corrcoef(aa,bb)[0,1])
    i = int(np.nanargmax(values))
    return int(candidates[i]),float(values[i])


def audit_frames():
    audit = read_csv(INTER/"auditoria/entradas.csv")
    manual = (DATA/"CLASSIFICACAO_PADRAOES_MANUAL.md").read_text()
    rows,shift_rows,frame_manifest = [],[],[]
    for row in audit:
        run = row["ensaio"]
        folder = DATA/run
        files = sorted(folder.glob("ImgA*.jpg"))
        assert len(files)==501,(run,len(files))
        assert [p.stem for p in files]==[f"ImgA{i:06d}" for i in range(501)]
        frames,hashes = [],[]
        for path in files:
            data = path.read_bytes()
            hashes.append(hashlib.sha256(data).hexdigest())
            with Image.open(BytesIO(data)) as image:
                frames.append(np.asarray(image.convert("L")))
        frames = np.stack(frames)
        assert frames.shape[1]==80
        frame_manifest.append(dict(ensaio=run,arquivos=[p.name for p in files],sha256=hashes))
        profiles = frames[:,20:60].mean(axis=1)
        shifts,scores = estimate_shifts(profiles)
        with np.load(OLD/"intermediate"/f"{run}_input.npz") as archived:
            np.testing.assert_array_equal(shifts,archived['raw_shifts'])
            np.testing.assert_allclose(scores,archived['scores'],atol=1e-10,rtol=1e-10)
        samples = set((0,250,499))
        if run==TARGET:
            samples.update(np.flatnonzero(shifts<=0).tolist())
        literal_error = 0.
        for i in sorted(samples):
            expected,score = literal_shift(profiles[i],profiles[i+1])
            assert shifts[i]==expected,(run,i,shifts[i],expected)
            literal_error = max(literal_error,abs(scores[i]-score))
        widths = np.maximum(np.r_[shifts,shifts[-1]],1)
        center = frames.shape[2]//2
        strips = [frame[:,center-w//2:center-w//2+w] for frame,w in zip(frames,widths)]
        canvas = np.concatenate(strips,axis=1).copy()
        seams = np.cumsum(widths)[:-1]
        for i,seam in enumerate(seams,1):
            feather = min(2,min(widths[i-1],widths[i])//2)
            if feather:
                w = np.linspace(0,1,feather)[None,:]
                mixed = strips[i-1][:,-feather:]*(1-w)+strips[i][:,:feather]*w
                canvas[:,seam-feather:seam] = mixed.astype(np.uint8)
        panorama = load_image(run)
        np.testing.assert_array_equal(canvas,panorama)
        roi = panorama[18:72]
        raw = roi.sum(axis=0,dtype=np.uint64)/54
        np.testing.assert_array_equal(raw,load_raw(run))
        wb = openpyxl.load_workbook(folder/'__dados_cropped_v11.xlsx',read_only=True,data_only=True)
        sheet = wb['Sheet']
        values = sheet.iter_rows(values_only=True)
        assert next(values)==('position_px','Voltage_0','Voltage_1','Voltage_2')
        values = list(values)
        wb.close()
        np.testing.assert_array_equal([v[0] for v in values],np.arange(len(raw)))
        xlsx_error = float(np.max(abs(np.array([v[3] for v in values])-panorama[20:60].mean(axis=0))))
        assert xlsx_error<1e-10
        label = [p.stem.lstrip('_') for p in folder.glob('*.txt')]
        assert len(label)==1 and pattern_name(label[0])==row['padrao']
        condition = run.split('_')[0]
        manual_rows=[line.split('|')[2].strip() for line in manual.splitlines()
                     if line.startswith(f'| {condition} |')]
        assert len(manual_rows)==1,(run,manual_rows)
        label_matches=pattern_name(manual_rows[0])==row['padrao']
        if row['incluido']=='True':
            assert label_matches,(run,label,manual_rows)
        bad = shifts<=0
        prepared = core.prepare(raw)
        starts = np.r_[0,seams]
        bad_pixels = np.zeros(len(raw),bool)
        for i in np.flatnonzero(np.r_[bad,bad[-1]]):
            bad_pixels[starts[i]:starts[i]+widths[i]]=True
        rows.append(dict(ensaio=run,padrao=row['padrao'],incluido=row['incluido']=='True',
            rotulo_manual_condicao=manual_rows[0],rotulo_condicao_concordante=label_matches,
            quadros=501,intervalos=500,deslocamentos_nao_positivos=int(bad.sum()),
            fracao_substituida=float(bad.mean()),deslocamentos_negativos=int((shifts<0).sum()),
            deslocamentos_nulos=int((shifts==0).sum()),primeiro_intervalo_substituido=int(np.flatnonzero(bad)[0]) if bad.any() else '',
            lag_mediano_px=float(np.median(shifts)),comprimento_px=len(raw),
            intensidade_media=float(roi.mean()),fracao_pixels_255=float((roi==255).mean()),
            fracao_quadros_pixels_255=float((frames[:,18:72]==255).mean()),
            fracao_quadros_pixels_ge250=float((frames[:,18:72]>=250).mean()),
            desvio_sinal=float(raw.std()),minimo_normalizado=float(prepared.min()),maximo_normalizado=float(prepared.max()),
            colunas_substituidas=int(bad_pixels.sum()),fracao_colunas_substituidas=float(bad_pixels.mean()),
            erro_reconstrucao=0,erro_extracao=0,erro_planilha=xlsx_error,
            verificacoes_literais=len(samples),erro_pearson_deslocamento=literal_error,
            panorama_sha256=sha(folder/'__panorama_v11.png')))
        for i,(lag,score) in enumerate(zip(shifts,scores)):
            shift_rows.append(dict(ensaio=run,quadro_i=i,quadro_j=i+1,lag_estimado_px=int(lag),
                largura_usada_px=int(widths[i]),substituido=bool(lag<=0),pearson_perfil=float(score),
                coluna_inicio=int(starts[i]),coluna_fim_exclusiva=int(starts[i]+widths[i])))
        print(f'Entrada {run}: reconstrução exata; {int(bad.sum())}/500 substituições.',flush=True)
    write_csv(TABLES/'auditoria_entradas.csv',rows)
    write_csv(OUT/'deslocamentos_reestimados.csv',shift_rows)
    write_json(OUT/'quadros_sha256.json',frame_manifest)
    return rows


def audit_metadata():
    """Confere os parâmetros nominais da aquisição em duas fontes locais."""
    rows=[]
    keys=['SensorGain','GammaCorrection','Brightness','Contrast','PixelGain','XDRContrast']
    for source in read_csv(INTER/'auditoria/entradas.csv'):
        run=source['ensaio'];folder=DATA/run
        settings={};inside=False
        for line in (folder/'_config.xsv').read_text(errors='strict').splitlines():
            if line.startswith('['):inside=line=='[Acquisition]'
            elif inside and '=' in line:
                key,value=line.split('=',1);settings[key]=value
        movie={}
        for line in (folder/'_config.mii').read_text(errors='strict').splitlines():
            if ':' in line:
                key,value=line.split(':',1);movie[key.strip()]=value.strip()
        rate=float(settings['Rate']);exposure=int(settings['ExposureNs'])
        assert settings['Name']==movie['Origin 1']==run
        assert rate==float(movie['Film speed 1'])
        assert abs(exposure*1e-9-float(movie['Shutter time 1']))<1e-12
        with Image.open(folder/'ImgA000000.jpg') as image:
            width,height=image.size
        assert width==int(settings['ROIW'])==int(movie['Width of image 1'])
        assert height==int(settings['ROIH'])==int(movie['Height of image 1'])
        assert int(settings['Frames'])==int(movie['Number of images 1'])==501
        rows.append(dict(ensaio=run,incluido=source['incluido']=='True',
            taxa_quadros_s=rate,exposicao_ns=exposure,profundidade_bits=int(settings['PixelDepth']),
            largura_quadro_px=width,altura_quadro_px=height,
            **{key:int(settings[key]) for key in keys},
            config_xsv_sha256=sha(folder/'_config.xsv'),config_mii_sha256=sha(folder/'_config.mii')))
    write_csv(TABLES/'metadados_aquisicao.csv',rows)
    chosen=[r for r in rows if r['incluido']]
    lookup={r['ensaio']:r for r in rows}
    compare=['taxa_quadros_s','exposicao_ns','profundidade_bits',*keys]
    return dict(status='PASS',arquivos_conferidos=62,
        taxas_incluidas=sorted({r['taxa_quadros_s'] for r in chosen}),
        exposicoes_ns_incluidas=sorted({r['exposicao_ns'] for r in chosen}),
        brilho_codigos_incluidos=sorted({r['Brightness'] for r in chosen}),
        m03q40_igual_m03q20_nos_campos_verificados=all(lookup[TARGET][k]==lookup['m03q20_001'][k] for k in compare),
        campos_comparados=compare,
        observacao='Valores nominais dos metadados; códigos de ajustes não foram convertidos em funções fotométricas.')


def audit_equations():
    rng = np.random.default_rng(27102026)
    x,y = rng.normal(size=(2,513))
    lags = np.arange(1-len(x),len(x))
    errors = []
    for estimator in ['biased','unbiased']:
        direct = np.correlate(x,y,'full')/(len(x) if estimator=='biased' else len(x)-abs(lags))
        errors.append(float(abs(core.correlation(x,y,estimator)-direct).max()))
    wavelengths = np.array([50.,125.,512.,3000.])
    scales = wavelengths/core.FOURIER_FACTOR
    r = np.correlate(x,y,'full')/len(x)
    actual = core.morlet_cwt(r,scales)
    selected = [0,1,123,len(x)-1,len(r)-2,len(r)-1]
    cwt_error = 0.
    for k,a in enumerate(scales):
        for b in selected:
            u = (np.arange(len(r))-b)/a
            expected = np.dot(r,np.pi**(-.25)/np.sqrt(a)*np.exp(-u*u/2)*np.exp(-6j*u))
            cwt_error = max(cwt_error,float(abs(actual[k,b]-expected)))
    assert max(errors)<1e-12 and cwt_error<1e-12
    maps = core.pair(x,y,wavelengths)
    reverse = core.pair(y,x,wavelengths)
    swap_error = float(abs(maps['wco']-reverse['wco'][:,::-1]).max())
    assert swap_error<1e-12
    identity = core.pair(x,x,wavelengths)
    self_error = float(abs(identity['wco'][:,len(x)-1]-1).max())
    assert self_error<1e-12 and identity['wco'].max()<=1+1e-12
    sx = np.zeros(2048); sx[500:756]=rng.normal(size=256)
    sy = np.r_[np.zeros(100),sx[:-100]]
    delay_maps = core.pair(sx,sy,np.array([50.,100.,300.]))
    peak_lags = delay_maps['lags'][np.argmax(delay_maps['wco'],axis=1)]
    np.testing.assert_array_equal(peak_lags,[-100,-100,-100])
    assert float(abs(delay_maps['wco'].max(axis=1)-1).max())<1e-12
    # An affine change of grayscale must leave the prepared signal unchanged.
    affine_error = float(abs(core.prepare(x)-core.prepare(17+2.3*x)).max())
    assert affine_error<1e-12
    rejected = 0
    for bad in [np.ones(128),np.r_[np.nan,np.zeros(127)],np.zeros((8,8))]:
        try: core.prepare(bad)
        except ValueError: rejected+=1
    assert rejected==3
    return dict(status='PASS',erro_correlacao_maximo=max(errors),erro_cwt_temporal=cwt_error,
        erro_troca_par=swap_error,erro_autocoerencia=self_error,lag_atraso_100=peak_lags.tolist(),
        erro_invariancia_afim=affine_error,entradas_invalidas_rejeitadas=rejected,
        pdf_sha256=sha(ROOT.parents[1]/'assets/schuck-bardo-wavelet-coherence.pdf'),
        paginas_pdf_conferidas=[3,5,6],paginas_impressas=[411,413,414],
        definicoes='32.4, 32.6, 32.8, 32.14–32.17',
        adaptacao='divisor N; o artigo usa N−|m| em 32.21–32.23',
        ambiguidade='32.16 média complexa versus 32.20 média de módulos; K=1 não altera a WCO')


def audit_statistics(inputs):
    acquisitions = records()
    runs = [a['ensaio'] for a in acquisitions]
    j = runs.index(TARGET)
    ann = [i for i,a in enumerate(acquisitions) if a['padrao']=='wavy annular']
    others = [i for i in ann if i!=j]
    with np.load(INTER/'descritores.npz') as saved:
        p = saved['power'].copy(); r = saved['r'].copy(); score = saved['score'].copy()
    reconstructed = np.array([direct_zero(load_raw(run)) for run in runs])
    error = float(abs(p-reconstructed).max())
    np.testing.assert_allclose(p,reconstructed,atol=1e-11,rtol=1e-11)
    np.testing.assert_allclose(r,descriptor_matrix(reconstructed,WAVES),atol=1e-12,rtol=1e-12)
    sensitivity = []
    for metric,matrix in [('pearson',r),('escore_wco',score)]:
        for excluded in [(),(TARGET,)]:
            sensitivity.append(dict(metrica=metric,conjunto='26_originais' if not excluded else '25_sem_m03q40',
                **aggregates(matrix,acquisitions,excluded)))
    write_csv(TABLES/'sensibilidade_exclusao_m03q40.csv',sensitivity)
    peer_r = [r[i,k] for i,k in combinations(others,2)]
    controls = []
    def append_control(name,matrix,power,wavelengths,indices=ann):
        jj = indices.index(j)
        peer = [k for k in range(len(indices)) if k!=jj]
        controls.append(dict(controle=name,ensaios=len(indices),
            pearson_m03q40_outros=float(matrix[jj,peer].mean()),
            pearson_entre_outros=float(np.mean([matrix[a,b] for a,b in combinations(peer,2)])),
            pico_m03q40_px=float(wavelengths[np.argmax(power[jj])]),
            pico_outros_min_px=float(min(wavelengths[np.argmax(power[i])] for i in peer)),
            pico_outros_max_px=float(max(wavelengths[np.argmax(power[i])] for i in peer))))
    append_control('panoramas_completos_18_71',r[np.ix_(ann,ann)],p[ann],WAVES)
    with np.load(OLD/'intermediate'/f'{TARGET}_input.npz') as old:
        first_bad = int(np.flatnonzero(old['raw_shifts']<=0)[0])
        prefix = int(old['shifts'][:first_bad].sum())
    assert all(acquisitions[i]['comprimento']>=prefix for i in ann)
    experiments = [('prefixo_antes_da_primeira_substituicao',WAVES),
        ('faixa_lambda_50_1024',np.geomspace(50,1024,200)),
        ('sem_remocao_tendencia_linear',WAVES),('media_linhas_20_59',WAVES),
        ('equalizacao_global_histograma',WAVES)]
    for name,waves in experiments:
        powers = []
        for i in ann:
            run = runs[i]
            raw = load_raw(run)
            if name.startswith('prefixo'):raw=raw[:prefix]
            elif name=='media_linhas_20_59':raw=load_image(run)[20:60].mean(axis=0)
            elif name=='equalizacao_global_histograma':
                roi=load_image(run)[18:72]
                hist=np.bincount(roi.ravel(),minlength=256)
                lut=255*(np.cumsum(hist)-.5*hist)/hist.sum()
                raw=lut[roi].mean(axis=0)
            powers.append(direct_zero(raw,waves,linear=name!='sem_remocao_tendencia_linear'))
        powers=np.array(powers)
        append_control(name,descriptor_matrix(powers,waves),powers,waves)
    write_csv(TABLES/'controles_m03q40.csv',controls)
    input_target = next(a for a in inputs if a['ensaio']==TARGET)
    lengths = np.array([a['comprimento'] for a in acquisitions])
    windows = []
    target_image = load_image(TARGET)[18:72]
    raw = load_raw(TARGET)
    for a,b in [(0,prefix),(prefix,len(raw))]:
        windows.append(dict(inicio=a,fim_exclusivo=b,comprimento=b-a,
            media=float(raw[a:b].mean()),desvio=float(raw[a:b].std()),
            fracao_pixels_255=float((target_image[:,a:b]==255).mean())))
    historical = read_csv(OLD/'displacement/results/run_summary.csv')
    h = next(a for a in historical if a['run']==TARGET)
    return dict(status='PASS',erro_perfis_por_quadratura=error,ensaio_13=TARGET,
        pearson_intra_m03q40=float(r[j,others].mean()),pearson_intra_demais_wavy=float(np.mean(peer_r)),
        entrada_m03q40=input_target,comparacoes_com_sem=sensitivity,controles=controls,
        prefixo_sem_substituicoes_px=prefix,primeiro_intervalo_substituido=first_bad,
        janelas_m03q40=windows,
        ciclos_na_escala_3000_minimo=float(lengths.min()/3000),
        ciclos_na_escala_3000_maximo=float(lengths.max()/3000),
        suporte_deslocamento_historico=dict(intervalos=int(h['intervals']),suportados=int(h['supported_count']),
            substituicoes_suportadas=int(h['supported_old_clamps']),
            fonte='_OLD/displacement/results/run_summary.csv',sha256=sha(OLD/'displacement/results/run_summary.csv'),
            reexecutado_nesta_auditoria=False))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--reuse-input-audit',action='store_true',help='Usa auditoria de quadros já concluída nesta revisão.')
    args=parser.parse_args()
    OUT.mkdir(exist_ok=True)
    equations = audit_equations()
    write_json(OUT/'equacoes.json',equations)
    if args.reuse_input_audit:
        inputs=read_csv(TABLES/'auditoria_entradas.csv')
        text_fields={'ensaio','padrao','rotulo_manual_condicao','panorama_sha256'}
        for row in inputs:
            for key,value in row.items():
                if key not in text_fields and value!='':row[key]=json.loads(value.lower())
            assert sha(DATA/row['ensaio']/'__panorama_v11.png')==row['panorama_sha256']
        manifest=json.loads((OUT/'quadros_sha256.json').read_text())
        assert len(manifest)==31
        for record in manifest:
            for name,digest in zip(record['arquivos'],record['sha256']):
                assert sha(DATA/record['ensaio']/name)==digest
    else:
        inputs=audit_frames()
    metadata=audit_metadata()
    statistics=audit_statistics(inputs)
    result=dict(status='PASS_COM_RESSALVAS_DE_ENTRADA_E_METODO',quadros=31*501,
        panoramas_reconstruidos=31,equacoes=equations,metadados=metadata,estatistica=statistics,
        scripts_sha256={n:sha(ROOT/'scripts'/n) for n in ['wavelet_core.py','estatistica.py','auditoria_adversarial.py']},
        decisao='Manter o conjunto original de 26 como referência; sinalizar m03q40 e apresentar sensibilidade com 25. Não excluir em função de Pearson.')
    write_json(OUT/'resultado.json',result)
    print(json.dumps(dict(status=result['status'],estatistica=statistics),ensure_ascii=False,indent=2),flush=True)


if __name__=='__main__':
    main()
