# Perfis de WPS alternativa: análise e leitura guiada

**Comece pelo [guia visual em HTML](RESUMO_EXECUTIVO.html)**. Ele explica como uma imagem virou uma curva, como duas curvas foram comparadas e o que significam os percentuais, com cinco figuras comentadas. As imagens estão incorporadas ao HTML, que pode ser lido sem conexão.

A análise investigou se aquisições do mesmo padrão apresentaram perfis mais parecidos. A tendência apareceu, com exceções e sobreposição. **Os 76,5% apresentados são concordância de ordenação entre pares, não acurácia de classificação.**

Foram verificadas 31 entradas e mantidas 26 aquisições da base v16_spatial2_18_71. A média das linhas 18–71 forneceu os sinais. Seus perfis formaram 87 pares do mesmo padrão e 238 entre padrões, comparados por distância de Hellinger.

## Leitura recomendada

1. [Leitura guiada em HTML](RESUMO_EXECUTIVO.html) ou [a mesma explicação em Markdown](EXECUTIVE_SUMMARY.md): conceitos, exemplos reais, distribuição dos pares, conta de A e limites.
2. [Conclusão da comparação entre v16 e v16.1](conclusão.md): o que mudou entre as análises.
3. [Relatório técnico](REPORT.md): equações, todas as variantes e resultados detalhados.
4. [Histórico da execução](PASSO_A_PASSO.md), [protocolo fixado](PROTOCOL.md) e [decisões](DECISIONS.md): procedência e reprodução. O passo a passo registra a execução; a introdução didática está no guia.
5. [Verificação final](results/final_verification.json), [manifesto](results/artifact_manifest.json) e [fontes das novas figuras](results/reading_guide.json).

As cinco figuras centrais têm uma leitura comentada. O final do guia contém um glossário e uma tabela que explica os demais gráficos da pasta. Todos os gráficos científicos continuam disponíveis em PNG e PDF.

## Atualização somente da apresentação

Para gerar as novas figuras e os documentos a partir das tabelas existentes, sem executar novamente o processamento dos sinais:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 src2/v16.1_perfis_wps/scripts/reading_guide.py
PYTHONDONTWRITEBYTECODE=1 python3 src2/v16.1_perfis_wps/scripts/write_reports.py
PYTHONDONTWRITEBYTECODE=1 python3 src2/v16.1_perfis_wps/scripts/export_html.py
PYTHONDONTWRITEBYTECODE=1 python3 src2/v16.1_perfis_wps/scripts/verify_final.py
```

O guia também utiliza a figura `figures/conclusao/01_imagem_sinal_perfil.png`, já incluída na entrega. Se precisar regenerá-la, execute `scripts/render_conclusion.py` depois de `scripts/verify_conclusion.py`.

## Reprodução

Executar a partir da raiz do projeto, com as dependências de [requirements.txt](requirements.txt). Os scripts resolvem a base pela localização do próprio arquivo e funcionam após a renomeação. A base e data_processed devem continuar disponíveis como pastas irmãs; esta entrega não duplica os 1,7 GB da análise de origem.

```sh
PYTHONDONTWRITEBYTECODE=1 python3 src2/v16.1_perfis_wps/scripts/preflight.py
PYTHONDONTWRITEBYTECODE=1 python3 src2/v16.1_perfis_wps/scripts/validate_operations.py
PYTHONDONTWRITEBYTECODE=1 python3 src2/v16.1_perfis_wps/scripts/run_profiles.py
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/v16-perfis-mpl python3 src2/v16.1_perfis_wps/scripts/render_results.py
PYTHONDONTWRITEBYTECODE=1 python3 src2/v16.1_perfis_wps/scripts/verify_conclusion.py
PYTHONDONTWRITEBYTECODE=1 python3 src2/v16.1_perfis_wps/scripts/render_conclusion.py
PYTHONDONTWRITEBYTECODE=1 python3 src2/v16.1_perfis_wps/scripts/reading_guide.py
PYTHONDONTWRITEBYTECODE=1 python3 src2/v16.1_perfis_wps/scripts/write_reports.py
PYTHONDONTWRITEBYTECODE=1 python3 src2/v16.1_perfis_wps/scripts/export_html.py
PYTHONDONTWRITEBYTECODE=1 python3 src2/v16.1_perfis_wps/scripts/verify_final.py
```

Esses comandos escrevem somente na nova análise e na pasta temporária de configuração do Matplotlib. Não executam lotes de reprodução do v5 ou v16. O protocolo e a configuração estão congelados; mudar o protocolo exige registrar uma nova decisão e uma nova versão da configuração, sem alegar que a mudança foi anterior à avaliação.

O lote numérico desta execução levou aproximadamente 106 s; a duração não inclui auditoria, gráficos e redação. PyWavelets apresenta metadados de distribuição 1.9.0 e módulo 1.8.0 neste ambiente. Os hashes e ambas as versões foram registrados. A Morlet principal usa NumPy/SciPy diretamente.

## Artefatos principais

| Arquivo | Conteúdo |
| --- | --- |
| results/input_verification.csv | Verificação das 31 extrações e diagnósticos de deslocamento |
| results/profiles_raw.csv | Perfis brutos, incluindo sensibilidades e valores sem suporte |
| results/profiles_normalized.csv | Massas, densidades, áreas e pesos por domínio |
| results/pair_prefix_profiles.csv | 350 combinações únicas de aquisição/comprimento para comparação justa por par |
| results/lag_support_counts.csv | Contagens integrais de lags por escala e aquisição |
| results/all_pairs.csv | 325 pares em cada uma das 21 representações/referências |
| results/ranking_summary.csv | Concordância A por par e balanceada |
| results/ranking_by_pattern.csv | Desfechos por padrão e controles que o envolvem |
| results/pair_group_summary.csv | Quatro grupos internos e dez combinações cruzadas |
| results/acquisition_deletion.csv | 26 retiradas de aquisição por método, sem intervalo de confiança |
| results/ranking_effects.csv | Diferenças absolutas e relativas de A |
| results/within_acquisition_sensitivity.csv | Mudanças de perfil da mesma aquisição |
| results/selected_cases.csv | Casos escolhidos após os resultados, para interpretação |
| results/illustrated_maps.csv | Origem, equação, denominadores e máximos integrais dos mapas |
| results/operations_validation.json | Controles sintéticos e testes das novas operações |
| results/source_snapshot.json | Hashes dos 583 arquivos da base antes da execução |
| results/execution.log | Registro de progresso e tempo do lote |
| results/reading_guide.json | Contagens explicadas, faixas do histograma e procedência das figuras didáticas |
| figures/ | 17 gráficos técnicos em PNG e PDF, em PT-BR |
| figures/conclusao/ | Três figuras da comparação entre v16 e v16.1 |
| figures/leitura/ | Quatro figuras didáticas adicionais em PNG e PDF |

As colunas de identificadores e os nomes dos arquivos preservam termos do código para reprodução. Textos e gráficos foram redigidos em PT-BR. Os rótulos de origem foram mapeados para bolhas dispersas, bolhas alongadas (plug), golfadas (slug), agitado (churn) e anular ondulado.

A pasta final é v16.1_perfis_wps. Os caminhos antigos v16.1 e v16_perfis_wps foram removidos por solicitação do usuário. src2 continua ignorado pelo Git; nenhum arquivo foi adicionado à força e nenhum commit foi feito.
