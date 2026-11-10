# Sinais, perfis e métricas wavelet

Esta etapa contém os sinais completos, os 26 perfis de potência e as métricas dos 325 pares. As comparações usaram prefixos de comprimento comum, cortados antes do preparo. O núcleo está em `../scripts/wavelet_core.py`: Morlet própria, frequência central 6, divisor N, 200 comprimentos de onda entre 50 e 3000 pixels e todos os atrasos.

A estatística vigente e os cinco mapas selecionados estão em `../04_analise_simplificada`. Esta pasta não contém a antiga análise por resíduos logarítmicos. `manifesto.json` registra fontes e produtos; `intermediario/auditoria/wavelet.json` registra a verificação numérica.

Para reproduzir toda a análise, use `python3 -B executar.py --workers 2` a partir da raiz da v18. Consulte `../METODOLOGIA.md` e `../README.md`.
