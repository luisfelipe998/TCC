# Método e reprodução da estatistica6

Execução exploratória do [plano original](PLANO.md), autorizada pela solicitação posterior do usuário, com triagem prévia e comparação direta plug/slug. O plano permanece preservado como registro. Configuração congelada antes do cálculo: [triagem_v1.json](configs/triagem_v1.json); seu SHA-256 está nos manifestos. Limites inclusivos do filtro usam tolerância numérica de 10⁻¹².

## Amostra e comparabilidade

Unidade amostral: ensaio. Conjuntos: 26 no inventário e matrizes; 19 churn/wavy annular; 6 plug/slug; 25 nos quatro padrões com réplicas. Dispersed bubbles é mantido nas comparações entre ensaios, mas excluído da avaliação de compactação. Nenhum registro foi excluído pelo resultado. O controle de 7.001 px contém exatamente os mesmos 19 churn/wavy annular, nunca as quatro classes. Os códigos m/q são usados como códigos; documentação histórica de unidades não substitui calibração das aquisições.

Perfil principal: magnitude da alternative WPS no lag zero, com autocorrelação enviesada, Morlet ω₀=6, 200 comprimentos de onda entre 50–3.000 px. Cortar antes de inverter intensidade, retirar tendência linear e dividir pelo desvio-padrão populacional. Preservar a banda e as equações da v16. Sem novo COI, alinhamento, reconstrução, suavização, alteração de faixa óptica ou quadratura da magnitude.

## Representações e métricas

- **Pearson:** retirar a reta de log(P+10⁻¹²) em log(λ), centrar o resíduo e normalizar sua norma. Distância euclidiana d=√[2(1−r)], não 1−r para a geometria. A1 relata r e mediana(r intra)−mediana(r inter).
- **Distribuição WPS:** t=log(λ/1 px); h são pesos trapezoidais em t; p=P·h/Σ(P·h). Hellinger=||√pᵢ−√pⱼ||/√2. Wasserstein=Σ|Fᵢ−Fⱼ|Δt. Não se interpreta p como energia física ou probabilidade de padrão.
- **Características:** log da escala do máximo, inclinação log-linear e largura Q75(t)−Q25(t), por inversa da CDF. Centro, entropia e participações [50,200), [200,800), [800,3000] são complementares. Empates, máximos locais e bordas foram registrados. Trio: distância euclidiana após mediana/IQR dos 26 prefixos, preservados nas demais posições; não há interpretação preditiva desse ajuste global.
- **Mean power:** média de |CWT(Rxx)| sobre todos os 2N−1 lags. Prefixos recalculados para os 26; centro/fim recalculados somente nos seis plug/slug aprovados preliminarmente. Comparação por Pearson residual, Hellinger e Wasserstein.
- **WCS no lag zero:** todos os 325 pares com N=3.171. Similaridade Σh|CWT(Rxy)(0,a)| / Σh√(PᵢPⱼ), nas escalas com denominadores válidos; dissimilaridade 1−similaridade. A CWT foi calculada por soma temporal complexa e conferida contra a FFT original. Este resumo não cobre a informação de todo o mapa WCS.
- **WCO:** apenas variant=baseline, 200 escalas, pesos nativos 18–71 do arquivo preprocessing_confirm_pairs.csv. Máximo global, média ponderada no lag zero e média ponderada de máximos por escala com piso 5%. Dissimilaridade 1−escore, sem clipping; não se presume métrica euclidiana. Comprimento mínimo por par, diferentemente dos demais descritores. Sem busca entre variantes/pisos, novas simulações ou extrapolação para WCO de comprimento fixo.

## Triagem e aprofundamento

Silhueta de cada ensaio s=(b−a)/max(a,b), com a=distância média aos outros ensaios do mesmo padrão e b=menor distância média a outra classe. Resumo principal: média das médias por classe, para não favorecer as maiores. Vizinho mais próximo exclui diagonal; empates recebem a proporção dos vizinhos empatados do mesmo padrão. Concordância balanceada é a média por classe. Referência aleatória: média de (n_g−1)/(n−1), igual a 47,222% nos 19 e 40% nos seis. Não são testes, previsões fora da amostra ou classificadores.

Filtro inicial: silhueta macro ≥0,10; média de silhueta de cada classe >0; concordância balanceada de vizinhança pelo menos 0,10 acima do acaso. Filtro geométrico seguinte: silhueta macro e vantagem de vizinhança positivas em centro/fim; mínimo leave-one-run-out da silhueta positivo no prefixo; silhueta positiva no controle 7.001 quando aplicável. Limites são convenções de prioridade, não margens científicas validadas. Todos os candidatos, falhas e condições permanecem nos CSVs.

**Passar o filtro geométrico não significa manter a direção do efeito, compactar cada classe em toda janela ou recuperar os rótulos por agrupamento.** Essas limitações são verificadas nos efeitos de Cliff, silhuetas por classe e agrupamentos. O pico plug/slug exemplifica a inversão que uma distância por si só não detecta. Não se alterou o filtro após observar essa limitação.

Agrupamento: ligação média e número de grupos fixo (k=2 nos dois contrastes que avançaram). Rótulos são usados somente para comparar a partição por ARI; 1 indica coincidência e 0 a expectativa ajustada ao acaso. Nenhuma escolha de k, linkage ou melhor janela. Os quatro padrões não passaram a triagem, portanto seu agrupamento não foi aprofundado. As janelas preservam a ordem dos ensaios; não são replicações experimentais.

## Blocos executados e limites

| Bloco | Execução |
| --- | --- |
| A1 | Contrastes de r, pares por classe, remoção do singleton, peso igual entre classes, retirada de cada ensaio e 5.000 perturbações de pesos por contraste binário. |
| A2 | Começo/centro/fim; 76 janelas distintas e 75 pares de janelas do mesmo ensaio; sobreposição registrada. Controle de comprimento nos mesmos 19. |
| A3 | Características, medianas e delta de Cliff por condição; sem testes inferenciais. |
| A4 | Curvas/diferenças descritivas completas; sem teste de máximo por escala. |
| A5 | Distâncias Hellinger/Wasserstein, 325 pares cada; aprofundamento conforme triagem. |
| A6 | Cobertura de m/q, Spearman, efeitos dentro de q; sem regressões que confundam padrão e operação. |
| A7 | Raios ao centroide dos vetores residuais unitários; ajuste √[n/(n−1)], sem PERMDISP. |
| A8 | Não ampliado: não resolve a separação e a triagem WCO foi desfavorável. |

Perturbações A1: pesos exponenciais positivos normalizados para média um dentro de cada classe; peso do par é o produto; somente pares distintos. Mediana ponderada pela inversa da CDF em 0,5, média dos adjacentes em empate (tolerância 10⁻¹²). Quantis 2,5–97,5% são sensibilidade à ponderação, não intervalos de confiança. Semente 2026102406, fluxos SeedSequence independentes. Os 10.000 resultados estão salvos.

Não se fizeram permutações inferenciais: permutabilidade não estabelecida. Existem 200 arranjos dentro de q em churn/wavy annular e 9 em plug/slug; dentro de (m,q), apenas um em ambos. Sem p-valores novos, família inferencial ou correção Holm a aplicar. Não houve bootstrap espacial, ANOSIM ou classificadores. Comparações de pares não foram tratadas como amostras independentes.

## Reprodução

Da raiz do projeto, com Python e dependências já instaladas:

```sh
python3 src2/estatistica6/scripts/screen.py
python3 src2/estatistica6/scripts/advance.py
python3 src2/estatistica6/scripts/validate_results.py
python3 src2/estatistica6/scripts/report.py
python3 src2/estatistica6/scripts/finalize.py
```

Os programas escrevem apenas em estatistica6 e no cache temporário do Matplotlib. Versões efetivas: ver manifest_started.json; PyWavelets 1.8.0 difere do manifesto histórico (1.9.0), mas os cálculos novos usam a Morlet própria, sem motor legacy. As fontes foram conferidas antes/depois por SHA-256. Executar novamente substitui somente os produtos desta análise. Relatórios e gráficos são derivados dos CSVs/NPZs salvos.

## Arquivos para auditoria

- [Triagem completa](results/triagem.csv), [decisões preliminares](results/decisoes_preliminares.csv), [decisões geométricas finais](results/decisoes_finais.csv).
- [Métricas por condição](results/metricas_finais.csv), [características](results/a3_features.csv), [efeitos](results/a3_effects.csv), [agrupamentos](results/cluster_summary.csv).
- [Influência dos candidatos](results/candidatos_leave_one_out.csv), [sensibilidade por pesos](results/a1_weight_sensitivity.csv), [janelas](results/a2_within_run_windows.csv).
- [Covariáveis](results/a6_associations.csv), [cobertura](results/a6_coverage.csv), [efeitos dentro de q](results/a6_within_q_effects.csv).
- [Validação de fontes e fórmulas](validation/numerical.json), [verificação independente final](validation/final_numeric.json), [manifesto](results/manifest.json).
