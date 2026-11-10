"""Gera a apresentação espectral vigente com os módulos locais da etapa 04.

Alias de 04_analise_simplificada/scripts/analisar.py. Requer as etapas 02 e 03.
Não publica arquivos no manuscrito nem reproduz análises históricas.
"""
from pathlib import Path
import sys

sys.dont_write_bytecode = True
from runtime import PACKAGE_ROOT

sys.path.insert(0, str(PACKAGE_ROOT / "04_analise_simplificada/scripts"))
from analisar import main


if __name__ == "__main__":
    main()
