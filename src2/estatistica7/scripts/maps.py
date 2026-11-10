"""Calcula WCO dos pares completos em comprimentos e posições fixos."""
from datetime import datetime,timezone
import time
from core import *


def main():
    assert json.loads((ROOT/'validation/engine.json').read_text())['status']=='PASS'
    runs,labels,inv,raw=load_inputs()
    coordinates=[];validation=[];started=time.perf_counter()
    csv_save('inventory.csv',inv)
    for condition in CFG['conditions']:
        rr,ll,x,coords=condition_inputs(condition,runs,labels,raw)
        n=x.shape[1];wave=Wavelet(n);coordinates.extend(coords)
        p=[];self_area=[];self_width=[];censor=[];auto_rows=[]
        for run,label,signal in zip(rr,ll,x):
            pp,co=wave.auto(signal);aa,ww,cc=self_descriptors(co,n)
            p.append(pp);self_area.append(aa);self_width.append(ww);censor.append(cc)
            auto_rows.append(dict(run=run,pattern=label,condition=condition,n=n,
                                 mean_self_wco=float((H/H.sum())@aa),censored_width_scales=int(cc.sum()),
                                 self_zero_max_error=float(np.max(abs(co[:,n-1]-1)))))
        p=np.array(p);pairs=np.array(list(combinations(range(len(rr)),2)))
        rows=[];mean=[];peak=[];zero=[];peak_lag=[];examples=[];example_pair=[];displays=np.linspace(0,2*n-2,301,dtype=int)
        # Exemplos definidos por identificador, antes de observar os mapas.
        examples_wanted=set()
        for g in PATTERNS:
            members=np.flatnonzero(ll==g)
            if len(members)>=2:examples_wanted.add((int(members[0]),int(members[1])))
        for a,b in [('churn','wavy annular'),('plug','slug')]:
            if a in ll and b in ll:examples_wanted.add(tuple(sorted([int(np.flatnonzero(ll==a)[0]),int(np.flatnonzero(ll==b)[0])])))
        for k,(i,j) in enumerate(pairs):
            co,valid=wave.cross(x[i],x[j],p[i],p[j])
            scores,profiles=map_summary(co,p[i],p[j],n)
            rows.append(dict(condition=condition,n=n,run_x=rr[i],run_y=rr[j],pattern_x=ll[i],pattern_y=ll[j],same_pattern=ll[i]==ll[j],**scores))
            mean.append(profiles['mean']);peak.append(profiles['peak']);zero.append(profiles['zero']);peak_lag.append(profiles['peak_lag'])
            if (i,j) in examples_wanted:
                examples.append(co[:,displays].astype(np.float32));example_pair.append([i,j])
            validation.append(dict(condition=condition,run_x=rr[i],run_y=rr[j],valid_scales=int(valid.sum()),
                                   max_wco=float(np.nanmax(co)),min_wco=float(np.nanmin(co)),
                                   shared_peak_le_scale_peak=bool(scores['common_lag_peak']<=scores['scale_peak_weighted']+1e-12)))
            if (k+1)%100==0:print(f'{condition}: {k+1}/{len(pairs)} pares; total {time.perf_counter()-started:.1f}s',flush=True)
        csv_save(f'scores_{condition}.csv',rows);csv_save(f'auto_{condition}.csv',auto_rows)
        np.savez_compressed(ROOT/f'results/profiles_{condition}.npz',runs=rr,labels=ll,pairs=pairs,n=n,wavelengths_px=LAM,
                            zero_wps=p,mean=np.asarray(mean),peak=np.asarray(peak),zero=np.asarray(zero),peak_lag=np.asarray(peak_lag),
                            self_area=np.asarray(self_area),self_width=np.asarray(self_width),self_width_censored=np.asarray(censor))
        np.savez_compressed(ROOT/f'results/examples_{condition}.npz',runs=rr,labels=ll,pairs=np.asarray(example_pair),
                            display_lags=displays-n+1,wavelengths_px=LAM,display_wco=np.asarray(examples),
                            note='Somente apresentação: 301 lags. Todos os resumos foram calculados nos 2N-1 lags integrais.')
        print(f'{condition}: concluído, {len(rr)} ensaios e {len(pairs)} pares.',flush=True)
    csv_save('window_coordinates.csv',coordinates)
    frame=pd.DataFrame(validation)
    assert frame.shared_peak_le_scale_peak.all()
    frame.to_csv(ROOT/'validation/maps.csv',index=False,float_format='%.17g')
    manifest=json.loads((ROOT/'results/manifest_started.json').read_text())
    assert all(sha(SRC/k)==v for k,v in manifest['sources'].items())
    json_write(ROOT/'results/maps_completed.json',dict(completed_utc=datetime.now(timezone.utc).isoformat(),
               elapsed_seconds=time.perf_counter()-started,n_conditions=len(CFG['conditions']),n_pair_maps=len(frame),
               total_auto_maps=len(coordinates),all_sources_unchanged=True,all_lags_used=True))


if __name__=='__main__':main()
