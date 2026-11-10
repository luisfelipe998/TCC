# Comparação do deslocamento — antes e agora

**Faz sentido melhorar a estimativa e explicitar onde ela falha. Os resultados não sustentam substituir automaticamente todos os panoramas.** Esta extensão trabalha com os quadros nativos e escreve somente em `displacement/`.

Abra a [galeria com todos os registros](comparison.html), a [visão geral](figures/overview.png), o [caso m03q40_001](figures/m03q40_001_comparison.png) e o [controle dispersed bubbles](figures/m07q01_001_comparison.png).

## Resultado observado

- Foram processados **31 registros × 500 intervalos = 15.500 pares de quadros**, conservando os rótulos da tabela fornecida pelo usuário e os indicadores de inclusão anteriores. Cinco registros continuam fora do escopo principal.
- Nos **26 registros elegíveis**, **4,888/13,000 intervalos (37.6%)** passaram os critérios espaciais e temporais. Os demais permanecem sem estimativa aceita (`new_dx = NaN`).
- Nos mesmos intervalos aceitos, a correlação 2-D média das estruturas após remoção do fundo passou de **0.874534 para 0.878614**. Houve **3690 aumentos (75.5%), 1191 reduções e 7 empates** com tolerância de 1e−6. O ganho médio foi pequeno. A média atribui peso igual por intervalo aceito, não por aquisição.
- **m03q40_001:** 171/500 intervalos passaram; as **60 substituições antigas por 1 px continuam sem suporte suficiente**. Nenhuma foi preenchida por interpolação ou pela velocidade de trechos vizinhos. No subconjunto aceito, a correlação passou de 0.528059 para 0.536088.
- **m07q01_001:** 481/500 intervalos passaram; a correlação passou de 0.938421 para 0.947095. Serve como exemplo de estruturas com correspondência mais clara.
- Foram gerados **22 pares de panoramas**, somente para segmentos contínuos com pelo menos oito intervalos verificados. Esses segmentos não representam automaticamente o registro completo.

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
| m03q01_001 | plug-slug | não | 341/500 | -11.000 | -3.224 | 0.003630 | 177 |
| m03q05_001 | churn | sim | 94/500 | -21.000 | -20.259 | 0.002727 | 9 |
| m03q10_001 | wavy annular | sim | 7/500 | -23.000 | -11.166 | 0.006075 | 0 |
| m03q20_001 | wavy annular | sim | 23/500 | -42.000 | -42.953 | 0.002716 | 0 |
| m03q40_001 | wavy annular | sim | 171/500 | -49.000 | -52.295 | 0.008029 | 30 |
| m03q80_001 | smooth annular | não | 0/500 | — | — | — | 0 |
| m04q01_001 | slug | sim | 325/500 | -6.000 | -5.160 | 0.011799 | 105 |
| m04q05_001 | slug-churn | não | 379/500 | -10.000 | -8.656 | 0.008524 | 65 |
| m04q10_001 | churn | sim | 201/500 | -17.000 | -16.467 | 0.006110 | 23 |
| m04q20_001 | churn | sim | 183/500 | -18.000 | -18.386 | 0.005269 | 23 |
| m04q40_001 | wavy annular | sim | 0/500 | — | — | — | 0 |
| m04q80_001 | wavy annular | sim | 0/500 | — | — | — | 0 |
| m05q01_001 | plug | sim | 302/500 | -12.000 | -11.965 | 0.000457 | 93 |
| m05q05_001 | slug | sim | 395/500 | -13.000 | -12.831 | 0.001968 | 191 |
| m05q05_002 | dispersed bubbles | não | 500/500 | -9.000 | -9.121 | 0.000505 | 500 |
| m05q10_001 | churn | sim | 312/500 | -14.000 | -14.171 | 0.003837 | 53 |
| m05q20_001 | churn | sim | 302/500 | -15.000 | -14.872 | 0.002343 | 41 |
| m05q40_001 | churn | sim | 144/500 | -20.000 | -19.947 | 0.005291 | 16 |
| m05q80_001 | wavy annular | sim | 0/500 | — | — | — | 0 |
| m06q01_001 | plug | sim | 375/500 | -16.000 | -15.712 | 0.000790 | 96 |
| m06q05_001 | slug | sim | 308/500 | -17.000 | -16.834 | 0.005638 | 31 |
| m06q10_001 | churn | sim | 258/500 | -14.000 | -13.537 | 0.001682 | 41 |
| m06q20_001 | churn | sim | 153/500 | -20.000 | -19.774 | 0.004655 | 20 |
| m06q40_001 | wavy annular | sim | 2/500 | -18.000 | -15.064 | 0.003639 | 0 |
| m06q80_001 | wavy annular | sim | 0/500 | — | — | — | 0 |
| m07q01_001 | dispersed bubbles | sim | 481/500 | -6.000 | -6.309 | 0.008674 | 332 |
| m07q05_001 | plug | sim | 366/500 | -15.000 | -14.996 | 0.000638 | 85 |
| m07q10_001 | slug-churn | não | 301/500 | -19.000 | -18.580 | 0.006167 | 42 |
| m07q20_001 | churn | sim | 302/500 | -16.000 | -15.748 | 0.002223 | 29 |
| m07q40_001 | churn | sim | 184/500 | -20.000 | -20.098 | 0.001819 | 11 |
| m07q80_001 | wavy annular | sim | 0/500 | — | — | — | 0 |

## Limites e recomendação

Este método estima movimento de estruturas ópticas sob uma aproximação de translação horizontal. Deformação, movimento vertical, velocidades diferentes entre interfaces, saturação e falta de textura continuam limitando a correspondência. Em m03q40_001, 17.7% dos pixels nativos da faixa 18–71 têm intensidade ≥250; remover o fundo não recupera informação saturada. A ausência de suporte no escoamento smooth annular não demonstra ausência de escoamento.

**Usar estes resultados como diagnóstico e como candidatos para segmentos verificáveis. Não substituir a entrada wavelet completa por concatenação desses segmentos:** isso criaria novas emendas e omitiria intervalos sem uma malha espacial validada. Os sinais dos segmentos têm comprimentos diferentes e não devem ser usados para concluir melhora de WCO sem o controle de comprimento discutido na revisão.

O caso m03q01_001, fora do escopo principal por ser plug-slug, recuperou correspondência para suas 116 substituições antigas. Isso não altera sua exclusão original nem resolve as 60 falhas do registro m03q40_001.

## Validação e reprodução

Os controles com translações conhecidas em ambos os sentidos tiveram erro máximo de **0.006513 px**; o controle com fundo estacionário mais intenso teve erro máximo de **0.034105 px**. Perfis constantes, periodicidade ambígua e uma sequência estacionária foram rejeitados. A reconstrução passou por um teste geométrico independente de intensidades lineares. Esses números não são uma estimativa de erro nos filmes reais.

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
