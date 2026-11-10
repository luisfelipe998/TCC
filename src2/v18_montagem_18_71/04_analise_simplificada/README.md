# Pearson direto da potência no atraso zero

Esta é a única comparação espectral vigente. O coeficiente de Pearson foi calculado diretamente entre os 200 valores dos perfis completos, sem logaritmo nem retirada de tendência espectral. O preparo dos sinais de intensidade continua incluindo retirada de tendência linear.

`scripts/analisar.py` lê as etapas 02 e 03, calcula os 325 coeficientes, seleciona o maior Pearson de cada padrão com repetição e o maior pico de coerência interpadrão, recalcula os cinco mapas selecionados e produz oito figuras locais. A seleção e as métricas usam os campos completos antes da redução de colunas para exibição. Os exemplos ópticos dos mesmos ensaios ficam associados nas tabelas de `dados/`.

`atualizar_optica.py` permite refazer somente a junção óptica depois de uma alteração dessa etapa; a execução completa não precisa dessa chamada separada. Modificações científicas precisam passar novamente pelas verificações do comando principal.

Para reproduzir toda a análise, use `python3 -B executar.py --workers 2` a partir da raiz da v18. Consulte `../METODOLOGIA.md` e `../README.md`.
