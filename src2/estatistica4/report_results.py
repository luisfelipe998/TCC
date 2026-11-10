"""Relatório cujos números são lidos das tabelas efetivamente calculadas."""
import json
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE/'results'


def table(frame, digits=4):
    def cell(value):
        if pd.isna(value):
            return 'não estimável'
        if isinstance(value,(float,np.floating)):
            return f'{value:.{digits}f}'
        return str(value)
    lines = ['| '+' | '.join(map(str,frame.columns))+' |',
             '| '+' | '.join(['---']*len(frame.columns))+' |']
    lines += ['| '+' | '.join(cell(v) for v in row)+' |'
              for row in frame.itertuples(index=False,name=None)]
    return chr(10).join(lines)


def generate_report():
    s = json.loads((RESULTS/'summary.json').read_text())
    v = json.loads((RESULTS/'source_validation.json').read_text())
    robustness = pd.read_csv(RESULTS/'robustness.csv').set_index('scenario')
    ecdf = pd.read_csv(RESULTS/'distance_summary.csv')
    old = pd.read_csv(RESULTS/'previous_analysis_rechecked.csv')
    dispersion = pd.read_csv(RESULTS/'dispersion_summary.csv').set_index('pattern')
    all_summary = pd.read_csv(RESULTS/'all26_descriptive_summary.csv')
    correlations = pd.read_csv(RESULTS/'feature_covariate_correlations.csv')
    influence = pd.read_csv(RESULTS/'leave_one_run_out.csv')
    windows = pd.read_csv(RESULTS/'within_run_window_distances.csv')
    q_counts = pd.read_csv(RESULTS/'pattern_by_q.csv')
    tests = pd.DataFrame([
        dict(Teste='ANOSIM',Estatística=s['ANOSIM_R'],p_exato_q=s['ANOSIM_p_exact_q'],p_Holm=s['ANOSIM_p_holm']),
        dict(Teste='PERMANOVA',Estatística=s['PERMANOVA_F'],p_exato_q=s['PERMANOVA_p_exact_q'],p_Holm=s['PERMANOVA_p_holm'])])
    names = {
        'principal_7001_inicio':'18–71, 7.001 px, começo (principal)',
        '7001_centro':'18–71, 7.001 px, centro',
        '7001_fim':'18–71, 7.001 px, fim',
        'linha40_7001_inicio':'Linha 40, 7.001 px, começo',
        'referencia20_59_7001_inicio':'20–59, mesmo motor, 7.001 px, começo',
        'comprimento_completo':'18–71, comprimento completo de cada ensaio',
        '3171_inicio':'18–71, 3.171 px, começo',
        '3171_centro':'18–71, 3.171 px, centro',
        '3171_fim':'18–71, 3.171 px, fim',
        'Spearman_residuos':'Spearman dos resíduos, 7.001 px, começo',
        'Pearson_log_sem_remover_inclinacao':'Pearson dos log-perfis, sem retirar inclinação'}
    sensitivities = pd.DataFrame({'Condição':[names[k] for k in robustness.index],
                                 'ANOSIM R':robustness.ANOSIM_R.to_numpy()})
    ref_r = robustness.loc['referencia20_59_7001_inicio','ANOSIM_R']
    weakest = influence.loc[influence.ANOSIM_R.idxmin()]
    strongest = influence.loc[influence.ANOSIM_R.idxmax()]
    ratio = dispersion.loc['churn','mean_adjusted_radius']/dispersion.loc['wavy annular','mean_adjusted_radius']
    peak_rho = correlations[(correlations.feature == 'peak_wavelength_px') & (correlations.covariate == 'n')].spearman_rho.iloc[0]
    ps = windows[(windows.window_a == 'principal_7001_inicio') & (windows.window_b == '7001_fim')]
    overlap_ch = ps[ps.pattern == 'churn'].overlap_fraction
    overlap_an = ps[ps.pattern == 'wavy annular'].overlap_fraction
    source = HERE.parent/'v16_spatial2_18_71/results'
    pairs = pd.read_csv(source/'pair_methods.csv')
    pairs = pairs[pairs.method == 'v16_biased']
    null = pd.read_csv(source/'null_diagnostic.csv')
    null_table = null[['run_x','run_y','observed_max','conditional_monte_carlo_p']]
    peaks = pd.read_csv(source/'per_run_methods.csv')
    peaks = peaks[peaks.included & (peaks.method == 'v16_biased')]
    peaks = peaks.groupby('pattern').peak_zero_lag_px.agg(n='size',mínimo='min',mediana='median',máximo='max').reset_index()
    report = f'''# Análise dos resultados v16_spatial2_18_71

## Resposta

**Sim, faz sentido reaplicar Pearson para comparar a forma dos perfis wavelet por ensaio. Nos novos sinais, o descritor principal não sustentou separação robusta entre churn e wavy annular.** ANOSIM R = {s['ANOSIM_R']:.6f}, p exato condicionado a q = {s['ANOSIM_p_exact_q']:.3f}. Após Holm para os dois testes definidos no plano, p = {s['ANOSIM_p_holm']:.3f} para ANOSIM e {s['PERMANOVA_p_holm']:.3f} para PERMANOVA.

Isso não prova equivalência entre padrões ou ausência de informação útil nos sinais. Este conjunto pequeno e observacional, com este descritor, não forneceu evidência convincente de separação estável. O efeito observado foi menor que o histórico e variou com o trecho escolhido. Não se selecionou o trecho ou a transformação com resultado mais favorável.

Os [métodos e fórmulas](METODOLOGIA.md) foram definidos antes do cálculo novo. Todos os números vêm dos CSVs e JSONs desta pasta ou das tabelas de origem identificadas. Os controles sintéticos validam o código e não integram os resultados experimentais.

## Verificação dos dados e das fórmulas

- Os 31 sinais salvos e seus CSVs coincidiram com a média das 54 linhas 18–71 dos PNGs; hashes dos panoramas e CSVs foram conferidos.
- Os 31 perfis completos de alternative WPS no lag zero e os 78 perfis de janelas existentes foram recalculados independentemente. A maior diferença relativa L2 foi {v['full_profile_max_relative_l2']:.3e} nos perfis completos e {v['window_profile_max_relative_l2']:.3e} nas janelas.
- A Eq. 32.17 foi conferida nos coeficientes complexos dos 87 mapas intrapadrão salvos, com diferença máxima {v['saved_maps_equation_max_abs_error']:.1f}. A cobertura e os resumos foram cruzados com as tabelas.
- Um trecho experimental novo de 7.001 px foi conferido nas 200 escalas contra a FFT da v16. Diferença absoluta máxima: {s['experimental_direct_profile_vs_fft_max_abs_error']:.3e}.
- Os controles verificaram Pearson contra SciPy; ANOSIM por uma identidade independente baseada em AUC, incluindo empates; PERMANOVA contra ANOVA univariada e somas de quadrados de distâncias; enumeração exata; Holm; ortogonalidade dos resíduos; autocorrelações por somas literais; Morlet por convolução temporal; e uma senoide de comprimento de onda conhecido.

Não se refez a reconstrução dos panoramas nem a transformada completa de cada par. Os mapas foram conferidos a partir dos coeficientes armazenados; o recálculo independente foi dos perfis em zero. Os arquivos de origem consultados permaneceram inalterados, com hashes no [manifesto](results/manifest.json). [Controles de fórmulas](results/formula_validation.json) e [validação das fontes](results/source_validation.json) tiveram status PASS.

Na v16 enviesada, a mediana do máximo de alternative WCO entre os 87 pares foi {pairs.max_wco.median():.4f}; o maior máximo foi {pairs.max_wco.max():.4f}. Um máximo procurado em todo o plano de lag e escala não demonstra similaridade uniforme. Autocoerência igual a um em lag zero é uma identidade algébrica, não uma validação experimental. A versão não enviesada permanece uma sensibilidade metodológica e não foi usada como fonte do descritor principal.

Os picos dos perfis completos conferidos foram:

{table(peaks,1)}

São comprimentos de onda em pixels do panorama reconstruído, não medidas automaticamente equivalentes a comprimento de bolha, frequência temporal ou dimensão em milímetros.

## Relação com a estatística anterior

A análise anterior usou Pearson entre resíduos de log alternative WPS no lag zero, seguido de ANOSIM. Seus descritores arquivados reproduziram:

{table(old[['variant','n_runs','ANOSIM_R','p_exact_q','n_q_arrangements']])}

Essa conferência reexecutou a estatística sobre os descritores existentes; não refez a extração v8 das planilhas. Foram preservados os mesmos 19 ensaios, o comprimento comum de 7.001 px e a ideia do descritor. Para os novos dados, usaram-se a média 18–71 e a Morlet da v16, selecionando o trecho antes de retirar a tendência e normalizar.

As 200 escalas são coordenadas correlacionadas de um vetor, e os 171 pares compartilham ensaios. Nenhum deles foi tratado como observação independente. Não se aplicou p paramétrico de Pearson às escalas nem testes de duas amostras aos pares. A unidade de permutação foi o rótulo de ensaio.

O whitening anterior é interpretado estritamente como retirada de uma reta de log(P) contra log(escala). Não demonstra independência entre escalas nem eliminação de um fundo físico 1/f. Uma costura nominal abaixo de 50 px não prova ausência de contaminação dentro da banda; uma máscara que não cobre escala alguma não oferece um controle experimental adicional.

## Resultado principal

Foram usados 10 churn e 9 wavy annular, com 81 pares intrapadrão e 90 entre padrões. O perfil foi a magnitude da alternative WPS no lag zero, transformada em log e com tendência log-linear removida. A dissimilaridade foi D = 1 − Pearson, nos primeiros 7.001 px.

{table(tests)}

O ANOSIM R pequeno indica pouca vantagem na ordenação das distâncias para pares do mesmo padrão. Pela identidade R = 2 AUC − 1, a proporção empírica de uma distância inter exceder uma intra, com meio crédito para empates, foi {(s['ANOSIM_R']+1)/2:.2%}. **Não é acurácia de um classificador.**

A PERMANOVA usou os descritores centrados e normalizados para norma um, cuja distância euclidiana é sqrt(2D). Seu R² marginal foi {s['PERMANOVA_R2_marginal']:.4f}, ou {s['PERMANOVA_R2_marginal']:.2%} da soma de quadrados nessa geometria. Não é R² ajustado por q, porcentagem física explicada pelo regime ou capacidade preditiva. O p bruto de {s['PERMANOVA_p_exact_q']:.3f} não foi interpretado isoladamente: após Holm, ficou em {s['PERMANOVA_p_holm']:.3f}.

As distribuições de distância apresentam forte sobreposição:

{table(ecdf[['category','n_pairs','median','mean']])}

Pearson alto também ocorre entre padrões diferentes; similaridade absoluta alta não torna o descritor discriminativo. A divisão q-vizinho/q-distante conserva o limiar histórico de 20 unidades do código q como descrição, sem controlar todos os efeitos de q.

![Matriz de dissimilaridade e ECDF](figures/pearson_ecdf.png)

## Alcance das permutações

As contagens por código q foram:

{table(q_counts)}

Os estratos q = 5 e q = 80 não trocam rótulos, pois contêm só uma classe. Os demais geram 4 × 5 × 10 = 200 arranjos distintos, incluindo o observado. O p exato é a fração de estatísticas pelo menos tão grandes quanto a observada.

Permutar somente dentro de m gerou {s['m_arrangements']} arranjos e p ANOSIM = {s['ANOSIM_p_exact_m_sensitivity']:.6f}. Essa sensibilidade não controla m e q simultaneamente. Dentro das combinações (m,q) houve apenas {s['joint_m_q_arrangements']} arranjo, sem troca possível. O esquema não permite isolar o padrão dos dois códigos operacionais ao mesmo tempo.

A resolução de 0,005 é conhecida, mas a principal limitação é a exigência de permutabilidade dos ensaios dentro dos estratos, não garantida por randomização experimental. Os valores-p são condicionais e exploratórios. Não existem réplicas suficientes em condições idênticas para demonstrar repetibilidade experimental independente.

## Sensibilidade ao observável e ao trecho

Todas as condições usam os mesmos 19 ensaios. Os R são diagnósticos, sem selecionar uma nova análise principal e sem acrescentar testes por condição.

{table(sensitivities)}

Com motor v16, ensaios e trecho iguais, a referência 20–59 teve R = {ref_r:.6f}, contra {s['ANOSIM_R']:.6f} em 18–71, diferença de {s['ANOSIM_R']-ref_r:.6f}. **Ampliar a faixa diminuiu a separação observada deste descritor neste trecho.** Não se testou a significância dessa diferença; não se conclui que 18–71 seja globalmente pior. A comparação com a v8 antiga, isoladamente, não separa mudanças de faixa, motor e ordem de pré-processamento.

Nas janelas de 3.171 px, R foi positivo no começo e no centro, mas negativo no fim. Isso enfraquece uma afirmação de separação espacialmente estável. O resultado maior na linha 40 não autoriza substituir retrospectivamente a média escolhida.

As janelas não são aquisições novas. Para começo versus fim de 7.001 px, a sobreposição em churn variou de {overlap_ch.min():.1%} a {overlap_ch.max():.1%}; em wavy annular foi {overlap_an.max():.1%}. Portanto, maior concordância entre janelas churn não implica maior estacionariedade física. O controle de 3.171 px mantém comprimentos de onda até 3.000 px, e o registro finito limita a interpretação das maiores escalas.

![Sensibilidade do ANOSIM](figures/robustness.png)

Ao retirar um ensaio por vez, R variou de {s['leave_one_out_R_min']:.6f} a {s['leave_one_out_R_max']:.6f}. O menor valor ocorreu sem {weakest.removed_run}; o maior, sem {strongest.removed_run}. Nenhum ensaio foi excluído por esse diagnóstico. A faixa mede influência, não é intervalo de confiança.

## Análises complementares realizadas

### Localização e dispersão

A PERMANOVA acrescenta uma leitura da diferença entre centroides. Foi acompanhada de dispersão descritiva porque dispersões distintas também podem contribuir para diferenças entre distribuições.

{table(dispersion.reset_index())}

A distância média ajustada ao centroide foi {ratio:.2f} vezes maior em churn. Isso exige cautela antes de descrever a diferença como puro deslocamento entre grupos compactos. Não se executou teste formal de PERMDISP; a razão não demonstra significância de dispersões diferentes. A PCA reteve {s['PCA_two_components_fraction']:.2%} da variância nos dois eixos mostrados, com sobreposição entre classes. Os demais eixos permaneceram nos testes. Os dois pontos mais afastados na projeção receberam rótulos; todos os identificadores estão no CSV.

![PCA e dispersão](figures/pca_dispersion.png)

### Consistência espacial e condições operacionais

Foram exportadas as distâncias entre janelas de cada ensaio e suas frações de sobreposição. Esse diagnóstico avalia dependência da região do panorama sem converter janelas em réplicas.

![Perfis e consistência entre janelas](figures/profiles_windows.png)

Nos 19 registros, Pearson entre comprimento original e deslocamento médio da reconstrução foi {s['shift_length_pearson_r']:.6f}; não são evidências independentes neste conjunto. Spearman entre o pico do perfil de 7.001 px e o comprimento original foi {peak_rho:.6f}. Truncar os registros iguala o número de pixels, mas não elimina associações com condições operacionais e reconstrução. A associação não pode ser atribuída exclusivamente a um artefato nem exclusivamente ao escoamento.

As associações descritivas de pico, inclinação e magnitude residual com m, q, comprimento e deslocamento estão em [feature_covariate_correlations.csv](results/feature_covariate_correlations.csv), sem busca de significância por variável.

### Cinco classes: descrição com comprimento igual

Os 26 ensaios foram comparados nos primeiros 3.171 px. Os picos abaixo se referem a essa janela, não ao perfil completo apresentado antes.

{table(all_summary,3)}

Plug e slug têm três ensaios cada; dispersed bubbles tem apenas um, impossibilitando estimar sua distância intraclasse. Curvas e intervalos de pico apresentam sobreposições. Não se aplicou teste confirmatório entre cinco classes nem se estimou acurácia de classificação com os mesmos dados.

![Perfis dos 26 ensaios](figures/all26_profiles.png)

### Diagnóstico nulo de alternative WCO já disponível

A v16 já contém máximos de WCO comparados a 199 pares de fases aleatórias por exemplo. Os valores-p foram conferidos nas replicações salvas:

{table(null_table)}

Nenhum dos cinco exemplos foi incomum ao nível de 5% sob esse modelo. Um pico alto de WCO não demonstra sincronização ou acoplamento entre aquisições. O diagnóstico é condicional ao modelo estacionário, não prova ausência de relação física e não representa um teste de todos os 87 pares.

## Conclusão

Pearson foi adequado como medida descritiva de forma, com teste por rótulos de ensaio sem pseudorreplicação. A versão 18–71 apresentou efeito pequeno e sensível à janela. Os resultados não autorizam afirmar separação robusta, repetibilidade experimental, classificação validada ou causalidade do padrão. Também não demonstram equivalência ou ausência de informação discriminativa nos sinais.

PERMANOVA, dispersão, PCA, influência por ensaio e variação espacial acrescentaram informação relevante: sobreposição entre classes, dispersões amostrais distintas e dependência do recorte. A prioridade para evidência mais forte é obter aquisições independentes repetidas em condições operacionais comparáveis. Se houver desenvolvimento posterior de classificadores, seleção de descritores e pré-processamento precisam ocorrer dentro da validação, mantendo todos os pixels e janelas de um ensaio no mesmo grupo de treino ou teste.

Reprodução computacional e fórmulas verificadas não eliminam limitações do delineamento.
'''
    (HERE/'RELATORIO.md').write_text(report)
    readme = f'''# estatistica4 — Pearson nos dados 18–71

**Resultado:** separação fraca entre churn e wavy annular. ANOSIM R = {s['ANOSIM_R']:.4f}; p exato por q = {s['ANOSIM_p_exact_q']:.3f}; p ajustado por Holm = {s['ANOSIM_p_holm']:.3f}. Não houve evidência convincente de separação robusta pelo descritor principal.

- [Relatório completo](RELATORIO.md)
- [Métodos, fórmulas e plano anterior ao cálculo](METODOLOGIA.md)
- [Resumo numérico](results/summary.json)
- [Verificação das fórmulas](results/formula_validation.json) e [das fontes](results/source_validation.json)
- [Matriz de Pearson](results/pearson_r_matrix.csv), [distâncias por par](results/pair_distances.csv) e [200 permutações exatas](results/exact_q_permutations.csv)
- [Sensibilidades](results/robustness.csv), [influência](results/leave_one_run_out.csv), [dispersão](results/dispersion_summary.csv) e [26 ensaios](results/all26_descriptive_summary.csv)
- [Manifesto, versões e hashes](results/manifest.json)

## Reproduzir

Na raiz do repositório:

    PYTHONDONTWRITEBYTECODE=1 python3 src2/estatistica4/run_analysis.py

O comando verifica as fórmulas antes dos dados, confere as fontes e gera descritores, tabelas, figuras e relatório. Escreve em estatistica4 e usa /tmp/estatistica4-matplotlib para o cache gráfico. As análises anteriores não são alteradas. As dependências estão em requirements.txt e as versões utilizadas, no manifesto.

Somente controles de fórmulas:

    PYTHONDONTWRITEBYTECODE=1 python3 src2/estatistica4/validate_analysis.py

O segundo comando atualiza apenas o JSON de controles. Para renovar o manifesto após mudanças, executar novamente o comando completo. As permutações experimentais são enumeradas integralmente; somente os controles sintéticos usam semente fixa. Campos vazios em comparações não aplicáveis ou dispersão não estimável não representam zero.
'''
    (HERE/'README.md').write_text(readme)


if __name__ == '__main__':
    generate_report()
