# Panoramas dos 26 ensaios

Os 501 quadros de cada ensaio foram alinhados pela média das linhas 18 a 71. As variantes `wavelet` e `optica` usaram as mesmas tiras centrais; somente a variante wavelet substituiu a coluna da emenda nas condições descritas na metodologia. `deslocamentos.csv` conserva os deslocamentos e larguras; `manifesto.json` identifica as 52 imagens. A montagem lê diretamente `data_processed`.

Para reproduzir toda a análise, use `python3 -B executar.py --workers 2` a partir da raiz da v18. Consulte `../METODOLOGIA.md` e `../README.md`.
