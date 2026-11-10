"""Publish the independently checked object-motion findings in Portuguese."""
from pathlib import Path
import gzip,hashlib,json
import numpy as np
import pandas as pd
from PIL import Image,ImageDraw
import object_motion as obj
from object_comparisons import eligible

BASE=Path(__file__).resolve().parents[1]

def run():
    a=json.loads((BASE/'results/object_independent_verification.json').read_text())
    b=json.loads((BASE/'results/comparison_independent_verification.json').read_text())
    assert not a['cache_without_objects_function_hash']
    # The final protocol was read after the added analyses and purge were audited.
    # Refresh that document snapshot only; all numerical inputs and scripts must
    # still have the hashes from the full numerical verification.
    for payload,name in ((a,'object_independent_verification.json'),(b,'comparison_independent_verification.json')):
        for relative,previous in list(payload['manifest'].items()):
            current=hashlib.sha256((BASE/relative).read_bytes()).hexdigest()
            if current!=previous:
                assert relative=='findings/PROTOCOL.md',('Numerical audit input changed',relative)
                payload['document_manifest_review']=dict(path=relative,previous_sha256=previous,reviewed_sha256=current,
                    reason='Final added-method protocol reviewed; numerical inputs and calculation scripts unchanged.')
                payload['manifest'][relative]=current
        (BASE/'results'/name).write_text(json.dumps(payload,indent=2)+'\n')
    mask=np.zeros((80,1016),bool);mask[30:40,300:310]=True
    tracks=obj.track([obj.objects(mask) for _ in range(21)])
    assert len(tracks)==1 and tracks[0]['accepted'] and abs(tracks[0]['signed_px_frame'])<1e-12
    assert not eligible(tracks[0])
    stationary=dict(synthetic=True,frames=21,known_velocity_px_frame=0.,original_tracker_accepted=tracks[0]['accepted'],
        recovered_velocity_px_frame=tracks[0]['signed_px_frame'],motion_resolved_filter_accepts=eligible(tracks[0]),
        interpretation='A stationary segmented square is trackable; it is not evidence of material transport.')
    (BASE/'results/stationary_object_control.json').write_text(json.dumps(stationary,indent=2)+'\n')
    examples=[('m04q01_001',4),('m05q01_001',30),('m05q05_001',91),('m07q01_001',0)]
    overview=Image.new('RGB',(1016,420),'white');draw=ImageDraw.Draw(overview)
    for k,(run,tid) in enumerate(examples):
        with gzip.open(BASE/'intermediate'/f'tracks_{run}_0.7.json.gz','rt') as f:sequence=json.load(f)
        candidate=next(t for t in sequence if t['track_id']==tid);node=candidate['nodes'][0]
        assert candidate['accepted']
        im=Image.open(BASE.parent/'data_processed'/run/f"ImgA{node['frame']:06d}.jpg").convert('RGB');box=ImageDraw.Draw(im)
        box.rectangle((node['xmin']-2,node['ymin']-2,node['xmax']+2,node['ymax']+2),outline='red',width=2)
        overview.paste(im,(0,k*105+24))
        draw.text((3,k*105+5),f"{run} | trajetória {tid} | quadro {node['frame']} | {candidate['signed_px_frame']:.4f} px/quadro",fill='black')
    overview.save(BASE/'figures/object_stationary_audit.png')
    table=pd.DataFrame(a['support'])
    names={('object','small'):'Objetos pequenos — protocolo original',('object','isolated_small'):'Pequenos afastados de candidatos grandes',
        ('object','large'):'Objetos grandes completamente fechados',('meniscus','large'):'Meniscos — regressão OLS original'}
    rows=['| Observável | Registros com suporte a 0,70 | Registros estáveis nos três limiares |','|---|---:|---:|']
    for r in table.to_dict('records'):rows.append(f"| {names[r['method'],r['group']]} | {r['baseline_supported']} | {r['threshold_robust']} |")
    counts=a['counts'];comp=b['counts'];paired=pd.read_csv(BASE/'results/paired_object_consensus.csv');paired=paired[paired.threshold_stable]
    pairs=['| Registro | Intervalos comuns | Pequenos/envoltória gasosa | Sensibilidade aos limiares |','|---|---:|---:|---:|']
    for r in paired.to_dict('records'):
        pairs.append(f"| {r['run']} | {int(r['common_intervals'])} | {r['ratio']:.6f} | {r['ratio_sensitivity_min']:.6f}–{r['ratio_sensitivity_max']:.6f} |")
    error=max(a['max_errors'][k] for k in ('OLS_velocity','OLS_left','OLS_right','OLS_volume'))
    report=f'''# Revisão independente dos movimentos de objetos e meniscos

## Resultado

Os cálculos de deslocamento, conversão pelo FPS, suporte e razões foram reproduzidos. **A auditoria não identificou um par de velocidades materiais líquido/gás.** Os resultados sustentam movimentos de candidatos ópticos e de envoltórias gasosas, com diferentes critérios de seleção. A concordância entre esses observáveis não valida microbolhas como traçadores do líquido.

Foram auditados os 26 registros selecionados, com 501 quadros cada e três limiares. As exclusões permaneceram as mesmas. Não foi escrito nem modificado dado em versões anteriores da investigação.

## Verificação numérica independente

O verificador reconstruiu os ajustes por mínimos quadrados com uma matriz de regressão e solução por decomposição numérica, em vez de chamar o ajuste do autor. Foram reavaliadas **{counts['object_tracks']:,} trajetórias de objetos** e **{counts['meniscus_tracks']:,} trajetórias de meniscos**, incluindo rejeitadas, em todos os registros e limiares. Esses totais não representam bolhas físicas únicas: incluem segmentações alternativas e fragmentos de trajetórias.

Foram reconstruídas {counts['object_links']:,} ligações de objetos e {counts['meniscus_links']:,} ligações de meniscos. Cada ligação esperada apareceu exatamente uma vez; cada nó armazenado pertenceu a uma única trajetória dentro da mesma combinação de método, registro e limiar. Não houve duplicação de ligações ou nós nessas combinações. Objetos fisicamente iguais em limiares diferentes não foram contados como observações independentes.

O detector foi reexecutado nos quadros 0, 125, 250, 375 e 500 de cada registro: **130 imagens nativas, {counts['native_frames_redetected_by_threshold']} combinações imagem/limiar e {counts['native_objects_redetected']} candidatos**, comparados com o cache. Os componentes, centroides, dimensões e pesos geométricos foram recalculados por uma implementação independente da extração de objetos. Para meniscos, foram reprocessadas todas as {counts['meniscus_masks_reprocessed']:,} combinações quadro/limiar a partir das máscaras de entrada. O detector de sombras em si é o mesmo algoritmo de v14; esta repetição verifica implementação e proveniência, não a verdade física da segmentação.

Os dois campos nativos de FPS foram novamente lidos em cada gravação. Os 13.026 arquivos de imagem tiveram nomes, sequência e conteúdo verificados contra os hashes do cache. A maior diferença entre velocidades OLS recalculadas foi {error:.3g} px/quadro; a maior diferença na conversão para px/s foi {a['max_errors']['native_fps_conversion']:.3g} px/s. Esses resíduos são aritméticos, não incertezas experimentais.

Houve {counts.get('roundoff_boundary_rejection_differences',0)} diferenças de enumeração do motivo de rejeição em comparações exatamente na fronteira do corte de deformação, explicadas por arredondamento numérico inferior a 10⁻⁹. Nenhuma alterou a aceitação final das trajetórias ou o suporte dos registros. O verificador registra essa equivalência numérica explicitamente; não oculta discrepâncias de decisão fora dessa tolerância.

O protocolo inicial foi confirmado:

{chr(10).join(rows)}

O cache original não vinculava os candidatos à versão da função de extração. O autor corrigiu isso incluindo o hash das funções de objetos e pesos geométricos e verificando novamente o conteúdo das imagens antes de reutilizar o cache. Após regeneração, a auditoria confirmou o hash em **todos os 26 caches**. Os manifestos representam os arquivos efetivamente auditados.

## Achado adversarial: objetos quase estacionários

O rastreador original permite acompanhar uma estrutura estacionária. Um controle executado com um quadrado segmentado imóvel durante 21 quadros foi aceito, com velocidade calculada nula. Isso é correto como rastreamento geométrico, mas não sustenta transporte material. Nos resultados reais a 0,70, 29 das 243 trajetórias pequenas aceitas tinham magnitude de velocidade ≤2 px/quadro: 8 em m04q01_001, 12 em m05q01_001 e 9 em m05q05_001.

Uma inspeção de exemplos encontrou pequenas regiões quase fixas próximas da mesma coluna, em torno de x=370–372, em gravações diferentes. Elas são compatíveis com uma contribuição óptica fixa; a inspeção não certifica a origem do artefato. Os exemplos estão na [figura de auditoria](../figures/object_stationary_audit.png). A comparação entre medianas com e sem esse grupo está na [tabela de sensibilidade](../results/object_stationary_sensitivity.csv).

A alternativa posterior exige magnitude superior a 2 px/quadro e deslocamento ajustado total de pelo menos 20 px. Sua seleção, todas as contagens, medianas e cálculo de deslocamento durante a exposição foram reproduzidos. Ela retém **{b['motion_resolved_small_baseline_tracks']} trajetórias pequenas no limiar de referência**, distribuídas em nove registros; possuir alguma trajetória não significa obter suporte no registro. O filtro elimina deslocamentos pequenos nesta seleção, mas também pode excluir movimentos físicos lentos; não certifica a fase de origem. Os limites foram declarados como análise exploratória, preservando a referência original.

## Meniscos robustos e razões em janelas comuns

As {comp['independent_robust_tracks']} regressões robustas foram recalculadas usando explicitamente a mediana de todas as inclinações entre pares de pontos de cada extremidade. A velocidade é a média das inclinações robustas das duas extremidades. Os resíduos, ajustes nas metades, critérios de aceitação e conversão para px/s coincidiram com os resultados publicados. O método usa seus próprios critérios e não exige aprovação pelo ajuste OLS anterior; nos dados atuais, todas as trajetórias robustas aceitas também satisfazem os diagnósticos anteriores de correspondência em atraso dois e reversão.

O ponto médio dos meniscos não é o centroide da área projetada. A auditoria verificou que sua velocidade OLS é a média das velocidades das duas extremidades. O centroide ponderado por diâmetros ao quadrado é uma sensibilidade geométrica condicionada à reconstrução circular; não demonstra volume tridimensional ou diâmetro físico corretos.

Foram reconstruídos **{b['common_intervals_checked']} intervalos pareados**, somados sobre os três limiares e todos os registros. As mesmas identidades de intervalos contribuíram para as duas coordenadas de cada par. Dois registros mantiveram suporte e estabilidade aos limiares:

{chr(10).join(pairs)}

Essas razões são entre **pequenos objetos e envoltórias gasosas**, não entre líquido e gás materiais. As coordenadas em cada intervalo são médias das velocidades ajustadas às trajetórias presentes. Cada ajuste utiliza a trajetória inteira, incluindo nós fora do conjunto de intervalos comuns. Portanto, há pareamento da janela de agregação, mas não uma estimação inteiramente restrita aos quadros dessa janela; não são velocidades instantâneas independentes de quadro a quadro.

Cinco registros sustentaram o mesmo observável de pequenos objetos em todos os limiares. Foram verificadas as **{b['inter_recording_pairs']} razões entre esses registros** e seus envelopes de sensibilidade. A razão física entre gravações ainda exige equivalência de magnificação axial. A variação entre registros de padrões diferentes não estabelece um efeito causal do padrão, e o único registro de bolhas dispersas não fornece replicação independente dentro desse padrão.

## Extrapolação para tamanho zero

Os {comp['extrapolation_models_recalculated']} ajustes, interceptos e erros de predição foram recalculados. O intercepto em tamanho zero permanece uma extrapolação óptica, sem validação como velocidade do líquido. O tamanho aparente, a posição radial, o desfoque e a seleção dos objetos podem influenciar o ajuste.

Os grupos de validação foram definidos pelo bloco temporal do ponto médio de cada trajetória. Na implementação final, qualquer trajetória de treinamento com intervalos nativos em comum com o teste é removida antes do ajuste. A auditoria reconstruiu essa remoção por interseção explícita dos conjuntos de intervalos e confirmou **{b['extrapolation_midpoint_folds_with_temporal_overlap']} sobreposições remanescentes em {b['extrapolation_midpoint_folds']} grupos**. Também conferiu o número de trajetórias que receberam predição quando o treinamento conservou posto suficiente. Trajetórias fragmentadas e condições comuns ainda impedem tratar os grupos como experimentos independentes. A [auditoria das sobreposições](../results/extrapolation_holdout_overlap_audit.csv) registra os casos e a quantidade removida.

## Limites e reprodução

O suporte depende de persistência, resolução, forma e possibilidade de correspondência. Não inclui necessariamente eventos rápidos, desfocados, coalescentes, fragmentados ou cortados. O lúmen foi definido como 18:71; o fechamento radial estrito exclui sombras que tocam essa fronteira, enquanto a alternativa de meniscos preserva esse contato como diagnóstico. Não é válido substituir um método pelo outro para preencher um valor ausente e manter o mesmo nome físico.

Esta revisão verificou a aritmética, a implementação e a proveniência dos métodos de objetos. Não demonstra ausência de mudança de fase, movimento sem escorregamento das microbolhas, validade de reconstrução circular, escala axial em m/pixel ou identidade material do líquido. A investigação de resíduos ópticos possui uma revisão separada.

Scripts: [auditoria dos objetos](../scripts/independent_object_verify.py), [auditoria das comparações](../scripts/independent_comparison_verify.py), [publicação deste relatório](../scripts/write_object_audit_report.py). Resultados completos: [objetos](../results/object_independent_verification.json), [comparações](../results/comparison_independent_verification.json), [controle estacionário](../results/stationary_object_control.json). Execute os dois verificadores após os cálculos e, em seguida, o publicador do relatório.
'''
    (BASE/'findings/object_independent_review_PT_BR.md').write_text(report)

if __name__=='__main__':run()
