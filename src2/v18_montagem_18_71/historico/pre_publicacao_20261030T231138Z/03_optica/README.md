# Análise óptica com faixa única de 18 a 71

Esta execução lê os 26 panoramas ópticos de `../01_panoramas/manifesto.json`,
formados por concatenação das tiras originais e alinhados com as linhas de 18 a 71.
Todas as operações de intensidade começam pelo recorte `18:72`: normalização,
referência líquida, suavização, segmentação e triagem de saturação. A máscara e
as duas medidas usam as mesmas 54 linhas, com diâmetro geométrico de 54 pixels.

## Reprodução

Na raiz do repositório, com os panoramas já existentes:

```sh
python3 -B src2/v14.1_new_voidfraction/scripts/run_all.py
python3 -B src2/v18_montagem_18_71/04_analise_simplificada/scripts/atualizar_optica.py
python3 -B TCC/figuras/v18_optica/gerar_figuras_opticas.py
python3 -B src2/v18_montagem_18_71/scripts/finalize_optical.py
```

O primeiro comando executa `scripts/process_optical.py` desta versão. Há uma
execução canônica compartilhada com a entrada v14.1, não duas replicações
independentes. Nenhum cálculo wavelet ou coeficiente de Pearson é reexecutado.
A junção usa os pares já selecionados em `04_analise_simplificada/dados`.
O fechamento verifica os máximos nessa matriz e os hashes das figuras.
A opção `--verification CAMINHO` associa um relatório independente com status PASS.

## Procedimento atual

O modelo compartilhado está em `../v14.1_new_voidfraction/scripts/panorama_void_alternatives.py`
(a partir do diretório da versão); o loader verifica seu SHA256 e isola as funções.
Após recortar as 54 linhas, cada coluna é normalizada pelo percentil 75 de toda
essa faixa. A seleção de referência usa o desvio padrão nas mesmas 54 linhas:
colunas no percentil 5 ou abaixo, nos ensaios `m05q01_001`, `m06q01_001` e
`m07q05_001`; medianas por linha e entre ensaios fornecem a referência de 54 valores.
A suavização radial gaussiana de desvio padrão de 1 pixel recebe somente esse
recorte, sem ler intensidades externas.

As buscas das sombras usam as metades derivadas da faixa: linhas 18–44 e 45–71.
Permanecem o limiar de 0,70, o fechamento de 3×3 pixels, a duração mínima de
20 colunas para interfaces abertas e o mínimo de 4 pixels por componente
isolado. A saturação é a proporção de pixels 255 nas 54 linhas; valores acima
de 20% excluem o ensaio dos resumos. Nenhum total de elegíveis é imposto.
A área circular condicional e as transições usam diâmetro e altura de 54 pixels.

## Resultados e verificações

- `intermediate/features/panorama_liquid_reference.json`: referência e colunas selecionadas.
- `intermediate/phase_masks/`: 26 máscaras, com zeros fora das 54 linhas analisadas.
- `results/void_fraction/` e `results/flip_rate/`: medidas, elegibilidade e resumos.
- `results/summary/optical_group_summary.csv`: medianas e médias por padrão.
- `results/validation/change_from_previous_optical.csv`: comparação antes/depois nas mesmas entradas v18.
- `results/validation/single_roi_checks.json`: invariância a intensidades externas, inclusão da linha 71 e normalização da máscara cheia.
- `results/validation/independent_arithmetic_checks.json`: intervalos e transições conferidos por cálculos alternativos.
- `results/validation/hull_intersections.json`: interseções verificadas por polígonos e viabilidade linear.
- `results/validation/input_hashes.json` e `artifact_manifest.json`: proveniência de entradas e produtos.
- `resumo_tcc.json`: síntese, concluída depois da junção e das figuras.

O estado anterior está preservado em `../historico/before_roi_unification_20261030`.
A comparação atual usa os mesmos panoramas antes e depois; alteraram-se as
janelas de processamento e triagem. As medidas permanecem condicionadas à
segmentação e à hipótese de seções de gás circulares sem sobreposição em
profundidade, sem calibração física independente.
