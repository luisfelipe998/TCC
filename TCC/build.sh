#!/usr/bin/env bash
# Build the TCC PDF with the abntex2 template (bibtex backend).
# Usage: ./build.sh   (run from anywhere; paths are absolute-safe via script dir)
set -euo pipefail

# TeX Live binaries from MacTeX-no-gui live here on Apple Silicon:
export PATH="/Library/TeX/texbin:$PATH"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "==> Using latexmk: $(command -v latexmk)"
echo "==> pdflatex:      $(command -v pdflatex)"
echo "==> bibtex:        $(command -v bibtex)"

# -pdf       : produce PDF via pdflatex
# -bibtex    : use bibtex (abntex2cite[alf] is bibtex-based, NOT biber)
# -outdir    : write ALL generated files (aux, log, bbl, pdf, ...) into dist/
# -halt-on-error / -interaction : fail loudly, no interactive prompts
# latexmk runs pdflatex/bibtex as many passes as needed to resolve refs/citations.
latexmk -pdf -bibtex -outdir=dist -halt-on-error -interaction=nonstopmode -file-line-error main.tex

echo "==> Done. Output: $SCRIPT_DIR/dist/main.pdf"
