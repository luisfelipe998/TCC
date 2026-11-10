"""Generate the local technical README from the completed comparison tables."""
from pathlib import Path
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def value(x, digits=3):
    return '—' if pd.isna(x) else f'{x:.{digits}f}'


def main():
    summary = pd.read_csv(ROOT / 'results/run_summary.csv')
    pairs = pd.read_csv(ROOT / 'results/displacements.csv')
    eligible = pairs[pairs.included]
    supported = eligible[eligible.supported]
    improved = int((supported.dynamic_correlation_change > 1e-6).sum())
    worsened = int((supported.dynamic_correlation_change < -1e-6).sum())
    tied = int((supported.dynamic_correlation_change.abs() <= 1e-6).sum())
    target = summary.set_index('run').loc['m03q40_001']
    control = summary.set_index('run').loc['m07q01_001']
    numerical = json.loads((ROOT / 'results/numerical_validation.json').read_text())
    text = f'''# Comparação do deslocamento — antes e agora

**Faz sentido melhorar a estimativa e explicitar onde ela falha. Os resultados não sustentam substituir automaticamente todos os panoramas.** Esta extensão trabalha com os quadros nativos e escreve somente em `displacement/`.

Abra a [galeria com todos os registros](comparison.html), a [visão geral](figures/overview.png), o [caso m03q40_001](figures/m03q40_001_comparison.png) e o [controle dispersed bubbles](figures/m07q01_001_comparison.png).

## Resultado observado

- Foram processados **31 registros × 500 intervalos = 15.500 pares de quadros**, conservando os rótulos da tabela fornecida pelo usuário e os indicadores de inclusão anteriores. Cinco registros continuam fora do escopo principal.
- Nos **26 registros elegíveis**, **{len(supported):,}/{len(eligible):,} intervalos ({len(supported)/len(eligible):.1%})** passaram os critérios espaciais e temporais. Os demais permanecem sem estimativa aceita (`new_dx = NaN`).
- Nos mesmos intervalos aceitos, a correlação 2-D média das estruturas após remoção do fundo passou de **{supported.old_dynamic_correlation.mean():.6f} para {supported.new_dynamic_correlation.mean():.6f}**. Houve **{improved} aumentos ({improved/len(supported):.1%}), {worsened} reduções e {tied} empates** com tolerância de 1e−6. O ganho médio foi pequeno. A média atribui peso igual por intervalo aceito, não por aquisição.
- **m03q40_001:** {int(target.supported_count)}/500 intervalos passaram; as **60 substituições antigas por 1 px continuam sem suporte suficiente**. Nenhuma foi preenchida por interpolação ou pela velocidade de trechos vizinhos. No subconjunto aceito, a correlação passou de {target.old_check_correlation:.6f} para {target.new_check_correlation:.6f}.
- **m07q01_001:** {int(control.supported_count)}/500 intervalos passaram; a correlação passou de {control.old_check_correlation:.6f} para {control.new_check_correlation:.6f}. Serve como exemplo de estruturas com correspondência mais clara.
- Foram gerados **{int((summary.segment_intervals > 0).sum())} pares de panoramas**, somente para segmentos contínuos com pelo menos oito intervalos verificados. Esses segmentos não representam automaticamente o registro completo.

Melhor alinhamento de imagem não determina o erro físico do deslocamento. O critério de suporte seleciona eventos, as linhas vizinhas e os quadros são correlacionados e a avaliação 2-D não constitui um conjunto experimental independente. Os campos relativos à correlação bruta, sem remoção do fundo, também foram preservados no CSV.

## O que mudou

1. **Fundo estacionário:** subtrai-se a média temporal em cada posição. Assim, reflexos e estruturas imóveis têm menor influência sobre a busca. A média usa todos os 501 quadros: é um procedimento offline.
2. **Faixa 18–71:** as 27 linhas pares (18, 20, …, 70) estimam o deslocamento; as 27 linhas ímpares (19, 21, …, 71) verificam a concordância. Essa separação não torna as observações estatisticamente independentes.
3. **Correlação normalizada com suporte fixo:** um modelo central de 616 colunas é comparado com janelas de igual comprimento. A busca permanece em ±200 px. Uma parábola de três pontos refina o pico para precisão subpixel.
4. **Sinal preservado:** positivo significa movimento para a direita na imagem; negativo, para a esquerda. O `s` antigo do v11 tem convenção oposta; os gráficos mostram `old_dx = −max(s, 1)` para comparar o mesmo sentido. Valores não positivos não são transformados em 1 px na nova estimativa.
5. **Suporte explícito:** correlação ≥0,60, diferença para pico concorrente fora de ±5 px ≥0,03, textura ≥1 nível de cinza e pico fora da borda da busca. Exige-se concordância com as outras linhas e com a busca em ordem inversa, dentro de `max(1 px, 15% × |dx|)`. A busca inversa usa outro modelo central e verifica consistência espacial; não é um retorno exato do mesmo recorte.
6. **Consistência temporal:** um deslocamento medido diretamente entre quadros separados por dois intervalos deve concordar com a soma dos dois adjacentes, dentro de `max(2 px, 15% × |soma|)`. Cada intervalo aceito precisa participar de pelo menos uma comparação válida; nenhuma comparação válida associada pode contradizê-lo. O CSV distingue suporte apenas no par (`base_supported`) de suporte também temporal (`supported`).

Os limites são heurísticos, herdados em grande parte da investigação anterior de quadros nativos. Foram fixados antes da comparação experimental desta pasta. Não são intervalos de confiança, e não foram ajustados para aumentar WCO, selecionar um padrão ou recuperar obrigatoriamente as 60 falhas.

O núcleo matemático está em [motion_core.py](scripts/motion_core.py). O estudo anterior consultado está em [frame_motion_evaluation](../../v15_gas_liquid_velocity/frame_motion_evaluation/README.md); nenhum script desse estudo foi executado com efeitos de escrita.

## Como ler as figuras lado a lado

- **Primeira linha:** deslocamentos antigos e novos ao longo de todos os 500 intervalos. Pontos cinza à direita são candidatos rejeitados, não movimentos aprovados. Pontos vermelhos à esquerda identificam substituições antigas por 1 px.
- **Segunda linha:** um exemplo de alinhamento de dois quadros após remoção do fundo. O intervalo aceito mais próximo do quadro 250 é escolhido, com desempate pelo primeiro; se não houver intervalo aceito, o quadro 250 é mostrado como candidato sem suporte. Não se escolhe o maior ganho visual. `r` mede correlação 2-D nas linhas não usadas no ajuste.
- **Terceira linha:** recorte exato do panorama v11 e nova reconstrução para os mesmos quadros do maior segmento contínuo aprovado. Os dois painéis usam a mesma escala horizontal e a mesma escala de cinza; comprimentos diferentes não são esticados para parecer iguais. Sem oito intervalos contínuos aprovados, não há novo panorama.
- **Quarta linha:** média das linhas 18–71 de cada PNG desse segmento. Não é uma nova análise wavelet, nem uma medida de fração de vazio.

Na nova construção, a posição de cada quadro resulta da soma dos deslocamentos subpixel. As fronteiras de suas células ficam a meio caminho entre posições vizinhas. A imagem é amostrada a cada 1 px por interpolação linear, termina nos centros do primeiro/último quadro observado e é arredondada para PNG de 8 bits. Não há estimativa inventada para um quadro futuro, mistura de emendas ou ponte sobre intervalos rejeitados. O mapa de origem de cada coluna é salvo em `intermediate/`.

A comparação de panoramas inclui mudanças tanto no deslocamento quanto na construção: o v11 concatena tiras inteiras, modifica emendas e reutiliza o último passo; a reconstrução nova usa posições cumulativas e células de ponto médio. Portanto, a diferença visual não isola um único fator. A comparação do alinhamento em pares nativos avalia o deslocamento sem depender da reconstrução do panorama.

## Resultado por registro

As medianas anterior e atual abaixo usam **exatamente os mesmos intervalos aceitos**. Um valor ausente indica falta de suporte, não velocidade zero. A amplitude da mudança e a cobertura precisam ser lidas juntas; dois intervalos aprovados não validam toda uma aquisição.

| Registro | Padrão | Escopo principal | Pares aceitos | Mediana anterior, px/quadro | Mediana atual, px/quadro | Variação média de r | Segmento, intervalos |
|---|---|---|---:|---:|---:|---:|---:|
'''
    for r in summary.itertuples():
        text += (f'| {r.run} | {r.pattern} | {"sim" if r.included else "não"} | {r.supported_count}/500 | '
                 f'{value(r.old_median_on_supported)} | {value(r.new_supported_median)} | '
                 f'{value(r.mean_check_correlation_change, 6)} | {r.segment_intervals} |\n')
    text += f'''
## Limites e recomendação

Este método estima movimento de estruturas ópticas sob uma aproximação de translação horizontal. Deformação, movimento vertical, velocidades diferentes entre interfaces, saturação e falta de textura continuam limitando a correspondência. Em m03q40_001, {target.saturation_fraction:.1%} dos pixels nativos da faixa 18–71 têm intensidade ≥250; remover o fundo não recupera informação saturada. A ausência de suporte no escoamento smooth annular não demonstra ausência de escoamento.

**Usar estes resultados como diagnóstico e como candidatos para segmentos verificáveis. Não substituir a entrada wavelet completa por concatenação desses segmentos:** isso criaria novas emendas e omitiria intervalos sem uma malha espacial validada. Os sinais dos segmentos têm comprimentos diferentes e não devem ser usados para concluir melhora de WCO sem o controle de comprimento discutido na revisão.

O caso m03q01_001, fora do escopo principal por ser plug-slug, recuperou correspondência para suas 116 substituições antigas. Isso não altera sua exclusão original nem resolve as 60 falhas do registro m03q40_001.

## Validação e reprodução

Os controles com translações conhecidas em ambos os sentidos tiveram erro máximo de **{numerical['maximum_translation_error_px']:.6f} px**; o controle com fundo estacionário mais intenso teve erro máximo de **{numerical['static_background_sequence']['max_error_px']:.6f} px**. Perfis constantes, periodicidade ambígua e uma sequência estacionária foram rejeitados. A reconstrução passou por um teste geométrico independente de intensidades lineares. Esses números não são uma estimativa de erro nos filmes reais.

O [verificador independente](scripts/verify_comparison.py) não importa o estimador. Confere janelas explícitas de Pearson em quadros reais, decisões e contagens, 15.500 intervalos, 15.469 verificações temporais, recortes antigos exatos, reconstrução dos novos PNGs, sinais e preservação por hash de todos os artefatos anteriores. Veja [verification.json](results/verification.json), [numerical_validation.json](results/numerical_validation.json) e [manifesto](results/artifact_manifest.json).

Execute a partir da raiz do projeto, com as mesmas dependências NumPy, SciPy, pandas, Pillow e Matplotlib da análise principal:

```sh
export MPLCONFIGDIR=/tmp/tcc-displacement-mpl
python3 -B src2/v16_spatial2_18_71/displacement/scripts/validate_displacement.py
python3 -B src2/v16_spatial2_18_71/displacement/scripts/run_comparison.py
python3 -B src2/v16_spatial2_18_71/displacement/scripts/write_report.py
python3 -B src2/v16_spatial2_18_71/displacement/scripts/verify_comparison.py
```

O inventário congelado `results/parent_snapshot.json` é a referência anterior à criação desta extensão. Não o substitua para mascarar alterações na análise original.

Saídas numéricas: [deslocamentos e motivos de rejeição](results/displacements.csv), [consistência temporal](results/temporal_checks.csv), [resumo por registro](results/run_summary.csv) e [parâmetros/procedência](results/method.json). Os PNGs dos segmentos estão em `panoramas/`, e seus sinais, em `signals/`.
'''
    (ROOT / 'README.md').write_text(text)
    result = dict(eligible_intervals=len(eligible), supported_intervals=len(supported),
                  supported_fraction=len(supported) / len(eligible),
                  improved=improved, worsened=worsened, tied=tied,
                  old_dynamic_correlation=float(supported.old_dynamic_correlation.mean()),
                  new_dynamic_correlation=float(supported.new_dynamic_correlation.mean()),
                  unresolved_m03q40_clamps=60, recommendation='Use verified segments as diagnostics; retain gaps and do not replace full wavelet inputs.')
    (ROOT / 'results/conclusion.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__': main()
