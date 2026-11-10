# Análise óptica após a remontagem nas linhas 18 a 71

Esta etapa lê os 26 panoramas ópticos do manifesto `../01_panoramas/manifesto.json`, reconstruídos pela concatenação das tiras na ordem dos quadros. O perfil de alinhamento desses panoramas usa `18:72`, isto é, 54 linhas de 18 a 71 inclusive. A análise óptica usa a mesma faixa para as máscaras, com altura e diâmetro geométrico de 54 pixels.

## Reprodução

Na raiz do repositório, depois da conclusão da reconstrução:

```sh
python3 -B src2/v18_montagem_18_71/scripts/process_optical.py
python3 -B TCC/figuras/v18_optica/gerar_figuras_opticas.py
python3 -B src2/v18_montagem_18_71/scripts/finalize_optical.py
```

O ambiente utilizado contém Python 3.13.2, NumPy, SciPy, Pillow e Matplotlib. As versões são registradas nos metadados. A execução sobrescreve somente os resultados ópticos desta versão e suas figuras novas.

A última etapa requer a matriz de Pearson e a seleção intrapadrão da nova análise espectral em `../02_wavelet/04_estatistica/dados/`. Confere os máximos na matriz, associa as medidas ópticas pelo identificador de aquisição e fecha o manifesto. A opção `--verification CAMINHO` inclui no resumo o hash de um relatório independente com status PASS.

## Método preservado e dados recalculados

`../scripts/optical_model.py` carrega somente cinco funções da implementação histórica, verificando o SHA-256 do arquivo de origem e compilando as definições em um namespace isolado. Os argumentos padrão recebem `LUMEN_ROWS=(18,72)` e `DIAMETER_PX=54.0`. Nenhum módulo histórico é importado ou modificado.

O limiar de 0,70, a suavização radial, o fechamento morfológico 3×3, as regras dos componentes, as janelas iniciais de busca e o comprimento mínimo de 20 pixels permanecem iguais aos da análise óptica de 54 linhas anterior. A referência líquida é **recalculada** dos novos panoramas de `m05q01_001`, `m06q01_001` e `m07q05_001`, com o mesmo critério de seleção de colunas. A referência anterior é usada apenas para registrar a comparação.

A triagem de saturação permanece na faixa auxiliar `20:60` e exclui dos resumos registros com mais de 20% dos pixels no valor 255. A saturação e a elegibilidade são recalculadas. Essa faixa de triagem não é a faixa das máscaras, e nenhum total de elegíveis é imposto.

## Resultados e verificações

- `intermediate/phase_masks/`: 26 máscaras com as dimensões dos novos panoramas.
- `intermediate/features/panorama_liquid_reference.json`: referência nova e seleção de colunas dos três panoramas de plug.
- `results/void_fraction/` e `results/flip_rate/`: medidas por aquisição, resumos por padrão e coordenadas do gráfico.
- `results/summary/optical_group_summary.csv`: medianas das tabelas e médias aritméticas das estrelas, com pesos iguais por aquisição.
- `results/validation/change_from_v14_2.csv`: comparação com a análise óptica anterior, incluindo largura, saturação, elegibilidade, fração de vazio e densidade de transições.
- `results/validation/independent_arithmetic_checks.json`: conferência dos intervalos contíguos por coluna e das contagens de transições, sem reutilizar o cálculo geométrico principal.
- `results/validation/hull_intersections.json`: envoltórias e interseções verificadas por geometria de polígonos e viabilidade linear.
- `results/validation/input_hashes.json`: hashes das entradas antes e depois da execução, incluindo os panoramas verificados contra o manifesto da reconstrução.
- `results/validation/artifact_manifest.json`: hashes dos artefatos e dos scripts de processamento.
- `resumo_tcc.json`: síntese para atualização do manuscrito. Os exemplos da Discussão são preenchidos depois da nova comparação espectral.

Não se comparam máscaras antigas e novas pixel a pixel, pois a remontagem pode alterar a largura e a correspondência das colunas. A comparação histórica é feita pelas medidas por aquisição. Os resultados permanecem descritivos e a reconstrução volumétrica, condicional à máscara e à hipótese circular.

As figuras em `TCC/figuras/v18_optica/` mantêm as funções de desenho da apresentação anterior. O gerador primeiro reproduz os PNGs antigos pixel a pixel, depois utiliza exclusivamente os panoramas, máscaras e tabelas novos. Não executa segmentação.
O rótulo vertical foi corrigido para indicar transições por 1.000 pares de pixels, conforme o denominador da medida, preservando cores e disposição gráfica.
