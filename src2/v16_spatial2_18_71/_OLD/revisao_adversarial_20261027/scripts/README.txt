ANÁLISE COMPLETA V16 · LINHAS 18–71

Executar a partir da raiz do projeto, nesta ordem:

python3 src2/v16_spatial2_18_71/scripts/analise.py --workers 2
python3 src2/v16_spatial2_18_71/scripts/estatistica.py
python3 src2/v16_spatial2_18_71/scripts/relatorio.py
python3 src2/v16_spatial2_18_71/scripts/verificar.py

analise.py conferiu 31 entradas, recalculou os 26 perfis individuais e os 325 pares,
e gerou 26 panoramas, 87 figuras intrapadrão e 238 entre padrões. Pares já concluídos
podem ser reutilizados se as assinaturas de código e dados coincidirem. --force
força o recálculo dos pares. O núcleo wavelet_core.py é idêntico ao original.

estatistica.py calculou Pearson, resumos e médias por ensaio, gerando quatro figuras.
relatorio.py escreveu a metodologia, o índice, as galerias e o resumo executivo.
verificar.py conferiu arquivos, fontes, cobertura, equações nos mapas, formatação,
matrizes, links e o conteúdo de _OLD. Os resultados da verificação ficam em
intermediario/auditoria/verificacao_final.json.

Os cálculos usaram todos os lags. Os mapas guardados foram reduzidos somente no
eixo dos lags para apresentação. As curvas médias e máximos salvos vieram das
matrizes integrais e não devem ser recalculados apenas das matrizes reduzidas.

01_panoramas: 26 figuras dos panoramas completos e dos sinais.
02_wavelets: 87 pares intrapadrão no formato aprovado.
03_wavelets_inter: 238 pares entre padrões no mesmo formato.
04_estatistica: quatro figuras, metodologia, índice, médias e tabelas numéricas.
05_conclusao: Resumo Executivo.md, resultados e interpretação.
intermediario: entradas, perfis, mapas, galerias e auditorias necessárias à reprodução.
_OLD: 583 arquivos originais e o piloto aprovado, preservados.

Os controles históricos de pré-processamento, extração e fases aleatórias foram
mantidos em _OLD e citados como resultados anteriores. Os cálculos principais da
versão atual foram reexecutados e confrontados com esses resultados de referência.

Dependências: numpy, scipy, matplotlib e Pillow; PyWavelets é importado pelo núcleo
original para manter os controles da implementação histórica disponíveis.
