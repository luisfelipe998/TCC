# Máscaras e medidas ópticas

O modelo local `../scripts/optical_model.py` aplica o recorte 18:72 antes da normalização, referência, suavização, segmentação e triagem. Usa referência dos três plugs, percentis 75 e 5, limiar 0,70, filtro gaussiano radial de um pixel e fechamento de 3×3. O diâmetro e a altura usados nas medidas são 54 pixels.

As 26 máscaras ficam em `intermediate/phase_masks`; os resumos incluem 25 ensaios, conforme a triagem de saturação. `results/summary/optical_group_summary.csv` reúne médias e medianas. As figuras estão em `figuras/`, as verificações em `results/validation/` e a comparação conjunta em `results/summary/discussion_example_pairs.csv`.

Para reproduzir toda a análise, use `python3 -B executar.py --workers 2` a partir da raiz da v18. Consulte `../METODOLOGIA.md` e `../README.md`.
