# Análise óptica atual: faixa única de 18 a 71

A análise óptica vigente usa exclusivamente as 54 linhas de 18 a 71, inclusive
(`18:72` em Python), desde a normalização e referência líquida até a suavização,
segmentação, triagem de saturação e medidas. O diâmetro adotado é 54 pixels.
O modelo compartilhado está em [scripts/panorama_void_alternatives.py](scripts/panorama_void_alternatives.py).

A execução canônica usa os panoramas ópticos de
[01_panoramas da v18](../v18_montagem_18_71/01_panoramas/manifesto.json), montados
com perfis nas mesmas 54 linhas. Seus resultados estão somente em
[03_optica da v18](../v18_montagem_18_71/03_optica/README.md).
O [manifesto atual](execucao_atual.json) identifica entradas, código e resultados
por SHA256. São entradas e cálculo compartilhados, não duas replicações independentes.

```sh
python3 -B src2/v14.1_new_voidfraction/scripts/run_all.py
```

Esse comando reproduz referência, máscaras, medidas dos 26 ensaios, triagem,
resumos e envoltórias, com verificações de invariância fora da ROI e aritmética
independente. Os consumidores do manuscrito e a junção com os pares espectrais
são etapas posteriores; nenhum cálculo wavelet ou classificador é reexecutado.

## Histórico preservado

A execução anterior desta pasta usava `18:71`, portanto 53 linhas, além de
janelas auxiliares distintas e alinhamento em 20:60. Seus scripts, resultados,
figuras, alternativas nativas, classificadores e relatórios foram preservados
sem alteração no [snapshot anterior à unificação](intermediate/history/before_roi_unification_20261030/README.md),
com [hashes](intermediate/history/before_roi_unification_20261030/snapshot_manifest.json).
Esses materiais documentam a execução anterior e não são resultados atuais.
Os panoramas antigos permanecem também em `intermediate/reconstructed_inputs`
para preservar referências históricas. O snapshot v14 mais antigo e a literatura
permaneceram intactos.

As estimativas permanecem condicionadas à máscara e à hipótese de seções de gás
circulares sem sobreposição em profundidade, sem calibração física independente.
