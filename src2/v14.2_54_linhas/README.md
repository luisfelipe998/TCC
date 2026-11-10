# Reprocessamento óptico com 54 linhas

Esta versão aplica a correção para as linhas **18 a 71, inclusive**: recorte Python `18:72`, altura de 54 pixels e diâmetro geométrico de 54 pixels. A versão histórica `v14.1_new_voidfraction` permanece intacta. Seus panoramas corrigidos são reutilizados sem reconstrução ou alteração de pixels.

## Reprodução

Na raiz do repositório:

```sh
/Users/i521862/.pyenv/versions/3.13.2/bin/python3 -B src2/v14.2_54_linhas/scripts/process_optical.py
```

Também é possível usar outro Python com NumPy, SciPy e Pillow instalados. As versões efetivamente usadas são registradas em `results/validation/metadata.json`. A reprodução sobrescreve somente os arquivos gerados nesta versão.

## Escopo e fidelidade ao algoritmo

`scripts/optical_model.py` lê o código histórico e verifica seu SHA-256. Compila exclusivamente as definições de `normalize`, `make_reference`, `close_shadow_crop`, `silhouette` e `geometry` em um namespace novo. Não importa módulos históricos, não modifica seus atributos e não usa seus caches. A geometria de 54 linhas é definida antes da compilação das funções, incluindo seus argumentos padrão. A única função substituída é a leitura da imagem, que abre os mesmos panoramas corrigidos sem efeitos colaterais.

O limiar permanece 0,70. Permanecem os mesmos três panoramas de referência, suas colunas selecionadas, a suavização, o fechamento morfológico 3×3, as regras dos componentes e o comprimento mínimo de 20 pixels. As janelas de busca inicial dos contornos permanecem `19:38` e `42:71`. A extensão radial e o recorte passam a poder incluir a linha 71. A referência reconstruída é comparada exatamente à referência salva da v14.1.

São processadas as mesmas 26 aquisições. A triagem de saturação permanece nas linhas `20:60`, com valor 255 e limite de 20%. A mesma aquisição saturada é excluída dos resumos, que contêm 25 aquisições. Não há novas varreduras de limiar, ajuste de parâmetros, inferência estatística ou processamento wavelet.

## Resultados e conferência

- `intermediate/phase_masks/`: máscaras binárias de todas as 26 aquisições, com as dimensões originais e primeiro plano restrito às 54 linhas.
- `results/void_fraction/`: fração de vazio condicional por ensaio, resumo por padrão e coordenadas do gráfico.
- `results/flip_rate/`: contagens de transições, denominadores `54*(W-1)` e resumos por padrão.
- `results/summary/optical_group_summary.csv`: medianas usadas na tabela e médias aritméticas usadas nas estrelas do gráfico, com pesos iguais por aquisição.
- `results/summary/discussion_example_runs.csv`: valores ópticos dos quatro ensaios dos exemplos da discussão. Os coeficientes de Pearson não são recalculados aqui.
- `results/validation/change_from_v14_1.csv`: comparação numérica e contagem de pixels alterados nas máscaras em relação à versão de 53 linhas.
- `results/validation/independent_arithmetic_checks.json`: reconstrução independente dos intervalos contíguos por coluna, releitura das máscaras salvas e conferência das contagens e da elegibilidade.
- `results/validation/hull_intersections.json`: interseções dos fechos convexos observados, verificadas por inclusão/cruzamento geométrico e por viabilidade de meias-planícies.
- `results/validation/input_hashes.json`: hashes dos panoramas, máscaras históricas, rótulos e demais entradas, iguais antes e depois da execução.
- `results/validation/artifact_manifest.json`: hashes dos scripts, deste README e dos artefatos produzidos, exceto o próprio manifesto.

A geometria é conferida por uma implementação independente que identifica cada intervalo radial contíguo e soma os quadrados de suas larguras divididos por `54²`. A densidade de transições é conferida linha a linha. A mudança de recorte pode alterar a conectividade e as máscaras, portanto os resultados não são obtidos por simples reescala dos valores históricos. A interpretação física permanece condicional à segmentação óptica e à reconstrução circular.
