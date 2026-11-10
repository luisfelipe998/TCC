estatistica8 — avaliação preliminar além de Pearson, usando v16.3_ancoras

Leia RELATORIO.txt para a decisão, os números e seus limites. As figuras ficam
em figuras/, as tabelas/NPZ em dados/ e os cálculos reproduzíveis em scripts/.
Nenhum arquivo da fonte v16.3_ancoras ou do texto do TCC é alterado.

Reprodução, a partir da raiz TCC:

  export PYTHONDONTWRITEBYTECODE=1
  export MPLCONFIGDIR=/tmp/estatistica8-mpl
  export OPENBLAS_NUM_THREADS=1
  export OMP_NUM_THREADS=1
  export VECLIB_MAXIMUM_THREADS=1
  python3 src2/estatistica8/scripts/extract.py
  python3 src2/estatistica8/scripts/analyze.py
  python3 src2/estatistica8/scripts/followup.py
  python3 src2/estatistica8/scripts/validate.py
  python3 src2/estatistica8/scripts/report.py

As etapas reescrevem apenas seus próprios resultados em estatistica8. A extração
reutiliza os 300 mapas congelados e calcula as duas sensibilidades sem salvar
novos mapas volumosos. analyze.py enumera exatamente as alocações binárias e
salva as distribuições nulas, com correção Holm conjunta. followup.py executa
controles de ruído, permutações por m/q e sensibilidades das combinações.

O protocolo contém os critérios definidos antes dos novos resultados e um adendo
posterior que substitui 9999 sorteios por enumeração exata das 92378 alocações de
churn/annular, preservando todos os candidatos, estatísticas e famílias Holm.

Dependências já disponíveis: NumPy, pandas, SciPy, Pillow e Matplotlib.
Versões registradas em dados/validacao_final.json. Não requer scikit-learn.

src2 está ignorado pelo .gitignore existente: os artefatos ficam disponíveis no
disco, mas não são automaticamente incluídos no controle de versão.
