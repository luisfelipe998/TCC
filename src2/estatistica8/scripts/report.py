"""Relatório preliminar e figuras exportáveis, sem alterar o texto do TCC."""
from core8 import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from datetime import datetime, timezone

COLORS = dict(zip(PATTERNS, ['#cb4b35', '#277eb3', '#279467', '#8c5ca6']))
NAMES = {
    'wps_pearson': 'Pearson · WPS residual', 'wps_width': 'Largura WPS',
    'wps_hellinger': 'Hellinger · WPS', 'wps_wasserstein': 'Wasserstein · WPS',
    'wps_log_peak': 'Pico WPS', 'wps_slope': 'Inclinação WPS', 'wps_trio': 'Trio WPS',
    'mean_power_pearson': 'Pearson · potência média', 'mean_power_hellinger': 'Hellinger · potência média',
    'mean_power_wasserstein': 'Wasserstein · potência média', 'mean_power_log_peak': 'Pico · potência média',
    'mean_power_width': 'Largura · potência média', 'mean_power_slope': 'Inclinação · potência média',
    'mean_power_trio': 'Trio · potência média',
    'wco_global_max': 'WCO · máximo global', 'wco_zero_weighted': 'WCO · lag zero ponderado',
    'wco_scale_peak_weighted': 'WCO · pico por escala', 'wco_scale_peak_power05': 'WCO · pico / piso 5%',
    'wco_common_lag_peak': 'WCO · pico em lag comum', 'wco_map_mean': 'WCO · média integral',
    'wco_map_mean_weighted': 'WCO · média ponderada', 'wco_lag_q90_weighted': 'WCO · percentil 90',
    'wco_area05_weighted': 'WCO · área ≥ 0,5', 'wco_mean_short': 'WCO · média 50–200 px',
    'wco_mean_middle': 'WCO · média 200–800 px', 'wco_mean_long': 'WCO · média 800–1024 px',
    'wco_common_peak_short': 'WCO · pico comum 50–200 px',
    'wco_common_peak_middle': 'WCO · pico comum 200–800 px',
    'wco_common_peak_long': 'WCO · pico comum 800–1024 px',
    'wco_central_mean': 'WCO · média central ±32 px',
    'self_area': 'Autocoerência · área', 'self_width': 'Autocoerência · meia largura',
    'self_bands': 'Autocoerência · três bandas',
    'kernel_area_full': 'WCO normalizada · área integral',
    'kernel_area_short': 'WCO normalizada · área 50–200 px',
    'kernel_zero': 'WCO normalizada · lag zero'
}


def read(name): return pd.read_csv(ROOT/'dados'/name)
def num(x, digits=3): return f'{x:.{digits}f}'.replace('.', ',')
def pct(x): return num(100*x, 1)+'%'


def save_figure(fig, name):
    (ROOT/'figuras').mkdir(exist_ok=True)
    for ext in ['png', 'pdf']:
        fig.savefig(ROOT/f'figuras/{name}.{ext}', dpi=180, bbox_inches='tight', facecolor='white')
    plt.close(fig)


def figures():
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'figure.dpi': 120, 'savefig.dpi': 180, 'axes.titleweight': 'bold'})
    features, metrics_frame, gain = read('caracteristicas.csv'), read('metricas.csv'), read('ganho_sensibilidades.csv')
    f = features[(features.family == 'wps') & (features.condition == 'ancora')]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3), layout='constrained')
    for k, pattern in enumerate(PATTERNS):
        rows = f[f.pattern == pattern].sort_values('run')
        x = k+np.linspace(-.15, .15, len(rows))
        axes[0].scatter(x, rows.width, color=COLORS[pattern], s=38, zorder=3, edgecolor='white', linewidth=.5)
        axes[0].hlines(rows.width.median(), k-.24, k+.24, color=COLORS[pattern], lw=2.5)
    axes[0].set(xticks=range(4), xticklabels=['churn\nn=10', 'wavy annular\nn=9', 'plug\nn=3', 'slug\nn=3'],
                ylabel='Largura WPS = ln(λ75 / λ25)', title='Maior evidência em churn × wavy annular')
    axes[0].grid(axis='y', alpha=.2)
    for pattern in PATTERNS[:2]:
        rows = features[(features.family == 'wps') & (features.pattern == pattern)]
        for run, group in rows.groupby('run'):
            values = group.set_index('condition').loc[['desloc_m64', 'ancora', 'desloc_p64'], 'width']
            axes[1].plot([-64, 0, 64], values, color=COLORS[pattern], alpha=.4, marker='o', ms=3, lw=1)
        med = rows.groupby('condition').width.median().reindex(['desloc_m64', 'ancora', 'desloc_p64'])
        axes[1].plot([-64, 0, 64], med, color=COLORS[pattern], lw=3, marker='o', label=pattern)
    axes[1].set(xticks=[-64, 0, 64], xlabel='Deslocamento do recorte em relação à âncora (px)',
                ylabel='Largura WPS', title='Mesmos 19 ensaios nas três posições')
    axes[1].legend(frameon=False); axes[1].grid(alpha=.2)
    fig.suptitle('Largura espectral: sinal promissor, ainda exploratório', fontsize=14)
    save_figure(fig, '01_largura_wps')

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3), layout='constrained')
    for name, combination, label, color in [
            ('wps_pearson', False, 'Pearson', '#777777'),
            ('wps_width', False, 'Largura WPS', '#277eb3'),
            ('wps_width', True, 'Pearson + largura', '#279467')]:
        subset = gain[(gain.scope == 'churn_annular') & (gain.metric == name) & (gain.combination == combination)]
        pos = subset[subset.control == 'leave_run_out'].set_index('condition').loc[['desloc_m64', 'ancora', 'desloc_p64']]
        axes[0].plot([-64, 0, 64], 100*pos.balanced_neighbor, color=color, marker='o', lw=2, label=label)
        for xx, yy in zip([-64, 0, 64], 100*pos.balanced_neighbor):
            axes[0].annotate(f'{yy:.1f}', (xx, yy), xytext=(0, 6), textcoords='offset points', ha='center', fontsize=9, color=color)
        cv = subset[subset.condition == 'ancora'].set_index('control').loc[['leave_run_out', 'leave_m_out', 'leave_q_out']]
        axes[1].plot(range(3), 100*cv.balanced_neighbor, color=color, marker='o', lw=2, label=label)
    axes[0].axhline(100*17/36, color='#999999', ls='--', lw=1, label='Referência aleatória 47,2%')
    axes[0].set(xticks=[-64, 0, 64], xlabel='Deslocamento do recorte (px)', ylabel='Concordância balanceada de vizinhos (%)',
                title='Sensibilidade espacial', ylim=(40, 104))
    axes[1].set(xticks=range(3), xticklabels=['Excluir\nensaio', 'Excluir\nmesmo m', 'Excluir\nmesmo q'],
                ylabel='Concordância balanceada de vizinhos (%)', title='Controles na janela original', ylim=(40, 104))
    for ax in axes: ax.grid(alpha=.2)
    axes[0].legend(frameon=True, facecolor='white', framealpha=.95, fontsize=8, loc='lower left')
    fig.suptitle('Churn × wavy annular: ganho observado, sem validação externa', fontsize=14)
    save_figure(fig, '02_vizinhos_e_sensibilidade')

    main = metrics_frame[metrics_frame.condition == 'ancora']
    names = list(NAMES)
    scopes = ['quatro_padroes', 'churn_annular', 'plug_slug']
    values = main.pivot(index='metric', columns='scope', values='min_class_silhouette').loc[names, scopes]
    fig, ax = plt.subplots(figsize=(9, 11), layout='constrained')
    im = ax.imshow(values, cmap='RdBu', norm=TwoSlopeNorm(vcenter=0, vmin=-.6, vmax=.45), aspect='auto')
    ax.set(yticks=range(len(names)), yticklabels=[NAMES[v] for v in names], xticks=range(3),
           xticklabels=['Quatro padrões', 'Churn ×\nwavy annular', 'Plug × slug'])
    ax.tick_params(axis='y', labelsize=8)
    for i in range(len(names)):
        for j in range(3):
            v = values.iloc[i, j]
            ax.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=8,
                    color='white' if v < -.38 or v > .32 else '#222222')
    ax.set_title('Menor silhueta média de classe: todas as 36 representações\nPositivo exige que ambas/todas as classes tenham separação média favorável', fontsize=11)
    cb = fig.colorbar(im, ax=ax, shrink=.55, pad=.03); cb.set_label('Menor silhueta de classe')
    save_figure(fig, '03_triagem_completa')

    obs = read('resumos_wco.csv').query('condition == "ancora"')
    noise = read('ruidos_independentes.csv')
    groups = [obs.loc[obs.pattern_x == obs.pattern_y, 'global_max'],
              obs.loc[obs.pattern_x != obs.pattern_y, 'global_max'],
              noise.loc[noise.model == 'branco', 'global_max'], noise.loc[noise.model == 'ar1_095', 'global_max']]
    fig, ax = plt.subplots(figsize=(8.5, 4), layout='constrained')
    ax.boxplot(groups, tick_labels=['Reais intra\n87 pares', 'Reais inter\n213 pares', 'Branco independente\n99 pares', 'AR(1) independente\n99 pares'],
               showfliers=True, widths=.5)
    ax.set(ylabel='Máximo global WCO', ylim=(0, 1.04), title='Máximos altos também ocorrem entre sinais independentes')
    ax.grid(axis='y', alpha=.2)
    save_figure(fig, '04_controle_maximos_wco')


def report():
    met, tests, gains = read('metricas.csv'), read('permutacoes.csv'), read('ganho_pearson.csv')
    effects, cov, conditional = read('efeitos_individuais.csv'), read('largura_covariaveis_por_escopo.csv'), read('permutacoes_condicionais_exploratorias.csv')
    get = lambda metric: met[(met.condition == 'ancora') & (met.scope == 'churn_annular') & (met.metric == metric)].iloc[0]
    width, pearson, trio = get('wps_width'), get('wps_pearson'), get('wps_trio')
    tw = tests[(tests.scope == 'churn_annular') & (tests.metric == 'wps_width')].iloc[0]
    tt = tests[(tests.scope == 'churn_annular') & (tests.metric == 'wps_trio')].iloc[0]
    gw = gains[(gains.scope == 'churn_annular') & (gains.metric == 'wps_width')].iloc[0]
    ew = effects[(effects.condition == 'ancora') & (effects.scope == 'churn_annular') & (effects.family == 'wps') & (effects.feature == 'width')].iloc[0]
    validation = json.loads((ROOT/'dados/validacao_final.json').read_text())
    text = f'''AVALIAÇÃO PRELIMINAR — estatistica8
Data de geração: {datetime.now(timezone.utc).isoformat()}
Fonte: src2/v16.3_ancoras. Comparação metodológica: src2/estatistica7.

RECOMENDAÇÃO

Não descartar toda a extensão além de Pearson. A largura espectral da WPS
apresentou separação prática interessante entre churn e wavy annular e merece
ser mantida como candidata, especialmente junto com Pearson. Entretanto, nenhuma
hipótese passou a correção global de múltiplos testes a 5%: o benefício permanece
exploratório, sem demonstração confirmatória de generalização.

Para separar plug de slug ou os quatro padrões simultaneamente, esta extensão
não entregou um separador consistente. Eu descartaria, para esse objetivo, os
resumos adicionais de WCO/autocoerência e a análise do pico plug/slug. O controle
de máximos WCO é um diagnóstico metodológico e não foi contado como sucesso de
identificação, conforme o critério solicitado.

O QUE FOI TESTADO

Foram preservados os 25 ensaios e seus rótulos: churn=10, wavy annular=9, plug=3,
slug=3. Cada ensaio teve uma janela fixa de 2048 px centrada na âncora existente.
Usaram-se os 300 pares únicos, sem contar diagonal ou pares espelhados como novas
observações. Cada aquisição, e não cada par/escala/janela, foi a unidade de análise.

As 200 escalas disponíveis cobrem 50–1024 px. Não foram importados os comprimentos
3171/5071/7001 px nem o teto de 3000 px da base antiga. A banda longa foi adaptada
para 800–1024 px; os dados não sustentam uma avaliação das escalas ausentes.

A triagem incluiu 36 representações admissíveis: 15 resumos cruzados WCO do estudo
anterior; WCO média central da fonte; 3 descritores de autocoerência; 3 kernels
WCO admissíveis; 7 representações WPS e 7 da potência média. Estas últimas foram
Pearson residual, Hellinger, Wasserstein, pico, inclinação, largura e trio
(pico/inclinação/largura). As análises WPS são as selecionadas em estatistica7,
com implementação original em estatistica6. Dois kernels de área (bandas média
e longa) foram recusados por autovalores negativos nas três condições, sem
reparar a matriz para obter artificialmente uma distância euclidiana.

Todos os candidatos foram avaliados nos seis contrastes binários e nas quatro
classes. Mediram-se silhueta macro e por classe, vizinhança balanceada, agrupamento
por ligação média com k fixo e ARI, influência de cada ensaio e ganho sobre Pearson.
Não se escolheram k, pesos ou bandas para maximizar os resultados.

RESULTADO MAIS ÚTIL: LARGURA WPS, CHURN × WAVY ANNULAR

A largura é L = ln(λ75/λ25), obtida dos quantis da WPS no lag zero normalizada
com pesos trapezoidais em ln(λ). Os quantis usam o primeiro ponto da grade que
atinge 25% e 75%, sem interpolação. Isso descreve distribuição espectral relativa;
não é energia física nem tamanho de bolha em milímetros.

Mediana de L: churn={num(ew.median_a)}; wavy annular={num(ew.median_b)}.
Faixas: churn {num(ew.min_a)}–{num(ew.max_a)}; annular {num(ew.min_b)}–{num(ew.max_b)}.
Delta de Cliff = {num(ew.cliff_delta)}; AUC de ordenação descritiva = {num(ew.rank_auc)}.
A AUC resume as 90 comparações dos MESMOS 19 ensaios, não 90 réplicas independentes.

                                      Pearson       Largura WPS
Silhueta macro                        {num(pearson.macro_silhouette):>7}        {num(width.macro_silhouette):>7}
Menor silhueta média entre classes     {num(pearson.min_class_silhouette):>7}        {num(width.min_class_silhouette):>7}
Vizinhança balanceada                  {pct(pearson.balanced_neighbor):>7}        {pct(width.balanced_neighbor):>7}

A combinação 50:50 de distância de Pearson e distância da largura atingiu
{pct(gw.combination_neighbor)} de concordância balanceada (18 de 19 ensaios corretos;
9/10 churn e 9/9 annular). As distâncias foram normalizadas pela mediana dos pares
do conjunto de treinamento em cada exclusão. A escala/IQR dos vetores também foi
estimada somente nos ensaios de treinamento para os testes de vizinhança.

São resultados de 1 vizinho com exclusão do próprio ensaio para candidatos fixos.
A escolha do candidato promissor ocorreu após a exploração desta amostra: 95%
não é uma estimativa imparcial de acurácia futura do processo de seleção, nem
validação externa. Não há conjunto de teste independente neste estudo preliminar.

O agrupamento sem supervisão não recuperou dois grupos limpos: ARI da largura foi
{num(width.ari)} e do trio {num(trio.ari)}. Proximidade útil para classificação por
vizinhos e recuperação de grupos por ligação média são resultados diferentes.

ESTABILIDADE E CONDIÇÕES OPERACIONAIS

Recalcularam-se TODOS os 300 pares após deslocar cada janela em -64 e +64 px,
mantendo N=2048 e a mesma grade. Os mesmos ensaios apareceram nas três condições;
isso não aumentou o tamanho amostral. A largura manteve silhueta macro de
0,417 / 0,435 / 0,357 e médias positivas nas duas classes. Sua vizinhança foi
78,9% / 81,1% / 66,7%, indicando alguma sensibilidade local.

Pearson + largura manteve 94,4% / 95,0% / 95,0%. Excluir todos os ensaios com o
mesmo código m do ensaio avaliado produziu 95,0% na janela original; excluir o
mesmo q também produziu 95,0%. Esses controles evitam a reutilização daquele
código operacional no treinamento, mas não eliminam toda a confusão física.

Ao retirar individualmente cada um dos 19 ensaios, os mínimos da silhueta macro
e da menor média de classe foram 0,406 e 0,370. Não foi necessário um único
ensaio específico para manter a diferença. No controle com panoramas completos,
a largura também permaneceu favorável (silhueta macro 0,330), com vizinhança menor
(64,6%). Esse controle muda conteúdo e comprimento juntos, não isola comprimento.

A largura ainda se associou a q (Spearman -0,591) e ao comprimento original do
panorama (-0,729) em churn/annular. Somente q=40 teve >=2 ensaios de cada classe
(2 churn e 3 annular); nesse subconjunto a vizinhança da largura foi 33,3%, embora
as silhuetas médias fossem ligeiramente positivas. Não se pode atribuir o efeito
exclusivamente ao padrão de escoamento.

SIGNIFICÂNCIA ESTATÍSTICA: RESULTADO AINDA NÃO CONFIRMADO

Os testes permutaram rótulos de AQUISIÇÕES, preservando contagens e todas as
dependências de pares. Enumeraram-se as 92378 alocações de churn/annular e todas
as alocações dos demais contrastes binários (20 em plug/slug, 286 ou 220 nos
outros); usaram-se 9999 permutações de semente fixa nas quatro classes.
O refinamento para enumeração exata está documentado no adendo do protocolo.

Para largura/churn-annular, p exato da silhueta={num(tw.p_silhouette, 6)};
após Holm sobre 504 testes (36 candidatos × 7 escopos × 2 desfechos),
p={num(tw.p_silhouette_holm, 4)}. O trio WPS teve o menor p corrigido da família,
{num(tt.p_silhouette_holm, 4)}, também acima de 0,05. A correção é conservadora para
descritores correlacionados, mas não foi reduzida após a escolha do resultado.

O ganho da combinação sobre Pearson teve p bruto={num(gw.p_gain_combination, 4)}
sob a referência de permutação de rótulos, porém p Holm={num(gw.p_gain_combination_holm, 3)}
na família própria de 504 comparações de ganho. Esse teste usa a hipótese de
permutabilidade dos rótulos; não é um teste de equivalência nem demonstra que um
preditor não possa superar outro em dados futuros. Ausência de significância
não significa efeito inexistente. Aqui há efeito descritivo forte e pouca
evidência confirmatória depois de considerar toda a busca realizada.

O aprofundamento por blocos m/q foi explicitamente posterior à triagem. Para
largura, as permutações exatas dentro de q deram p bruto=0,025 da silhueta, mas
p corrigido=1 na família exploratória conjunta de 216 testes. Dentro de m, o
p corrigido foi 0,125. Esses controles não converteram o achado em confirmação.

O QUE DESCARTARIA OU NÃO PRIORIZARIA

1. WCO cruzada, seus kernels e autocoerência como extensão para identificar os
   quatro padrões ou separar churn/annular e plug/slug: nenhum passou a triagem.
   Máximos, médias, percentis, áreas e bandas não deram separação estável útil
   nesses escopos. Dois kernels nem sequer tiveram geometria admissível.

2. Houve sinais WCO locais em outros contrastes: WCO central atingiu 83,3% de
   vizinhança para churn/plug, e lag zero ponderado atingiu 83,3% para annular/plug.
   Eles não foram escondidos. A retirada de um ensaio tornou a silhueta macro
   negativa; excluir o mesmo q reduziu a vizinhança para 45,0% e 44,4%. Não são
   bons candidatos de continuidade diante da fragilidade e da confusão operacional.

3. Pico WPS para plug/slug: silhueta macro -0,137, vizinhança 33,3% na âncora.
   Não reproduziu uma separação útil. Com apenas 3+3 aquisições, testes exatos de
   separação invariantes à troca dos nomes das classes têm resolução mínima de
   2/20=0,10; nem uma separação perfeita atingiria 5% nesse teste de separação.

4. Hellinger e Wasserstein: tiveram sinais em contrastes secundários, sobretudo
   churn/slug, mas não foram superiores à largura no contraste prioritário.
   Hellinger WPS em churn/annular teve uma classe com silhueta negativa na âncora.
   No panorama completo ele foi mais favorável, o que confirma a dependência do
   contexto, sem substituir silenciosamente a condição principal.

5. Nos quatro padrões, todas as 36 silhuetas macro da janela original foram
   negativas. A melhor combinação de vizinhos foi Pearson+largura, 42,2%, contra
   26,7% de Pearson sozinho. Esse ganho parcial não resolve a identificação das
   quatro classes e não passou a correção dos testes de ganho.

O trio e a largura da potência média passaram a triagem churn/annular, mas não
justificam ampliar a análise principal diante da simplicidade, magnitude e
interpretabilidade da largura WPS. Os resultados completos, inclusive os
contrastes secundários favoráveis, ficaram nas tabelas; não foram descartados.

CONTROLE DE MÁXIMOS E TESTE DE FASES

Nos 99 pares independentes de ruído branco e 99 AR(1)=0,95, os máximos WCO tiveram
medianas 0,864 e 0,846. Cerca de 37,4% e 29,3% atingiram pelo menos 0,9. São
simulações de independência conhecida, não réplicas experimentais nem um limiar
universal de significância. Esse resultado enfraquece o uso de um máximo alto
como identificador, mas não é considerado um resultado útil de classificação.

O nulo de fases não foi recalculado. O protocolo condicionou esse custo a algum
resumo WCO passar a triagem prioritária, o que não ocorreu. O teste de fases
investigaria o modelo espectral/temporal, sem recuperar uma separação que os
resumos observados não apresentaram. A conclusão negativa para WCO é restrita
aos resumos e dados testados, não uma prova de impossibilidade universal.

VERIFICAÇÃO E REPRODUÇÃO

Auditoria: {validation['status']}. Todos os {validation['sources_unchanged']} arquivos-fonte usados
mantiveram seus SHA-256. A WPS e autocoerência dos 25 ensaios foram confrontadas
com a fonte; seis mapas cruzados completos foram conferidos contra o motor novo
(erro máximo de WCO 1,78e-15). Os 300 resumos centrais coincidiram com os originais.
As silhuetas/vizinhos/ARI foram confrontados com funções puras de estatistica7;
as 1008 caudas de permutação e ambas as correções Holm foram rechecadas.
O resultado destacado da combinação foi reconstruído por cálculo independente,
com normalização exclusiva no treinamento. Nenhum arquivo da fonte/TCC foi editado.

Arquivos principais:
  protocolo.json — critérios, hipóteses e adendo de precisão.
  dados/metricas.csv — todos os métodos, escopos e posições.
  dados/decisoes.csv — triagem e sensibilidade por candidato.
  dados/caracteristicas.csv — valores por ensaio, sem pares artificiais.
  dados/ganho_pearson.csv — comparação e combinação com Pearson.
  dados/ganho_sensibilidades.csv — posições e exclusão por m/q.
  dados/permutacoes.csv — valores-p brutos e corrigidos de todos os candidatos.
  dados/permutacoes_condicionais_exploratorias.csv — controles posteriores por m/q.
  dados/validacao_final.json — auditoria, versões e preservação das fontes.
  figuras/01_largura_wps.png — efeito e estabilidade local da largura.
  figuras/02_vizinhos_e_sensibilidade.png — ganho e controles de vizinhança.
  figuras/03_triagem_completa.png — todas as 36 representações, sem seleção visual.
  figuras/04_controle_maximos_wco.png — diagnóstico de máximos independentes.

Os comandos completos estão em README.txt. Os gráficos também têm versão PDF.
'''
    (ROOT/'RELATORIO.txt').write_text(text, encoding='utf-8')
    save_csv('nomes_metricas.csv', [dict(metric=k, descricao=v) for k, v in NAMES.items()])
    json_write(ROOT/'dados/resumo.json', dict(
        recommendation='Manter largura WPS/combinação Pearson para churn versus wavy annular como candidato exploratório; descartar extensão WCO e pico plug/slug como separadores nesta base.',
        significant_after_global_holm=False, n_runs=25, candidates_admissible=36,
        churn_annular=dict(n_runs=19, silhouette_width=float(width.macro_silhouette),
            neighbor_pearson=float(pearson.balanced_neighbor), neighbor_width=float(width.balanced_neighbor),
            neighbor_combined=float(gw.combination_neighbor), width_cliff=float(ew.cliff_delta),
            width_p_exact=float(tw.p_silhouette), width_p_holm=float(tw.p_silhouette_holm),
            gain_p_exact=float(gw.p_gain_combination), gain_p_holm=float(gw.p_gain_combination_holm)),
        validation=validation['status']))


if __name__ == '__main__':
    figures()
    report()
    manifest = {str(p.relative_to(ROOT)): sha(p) for p in sorted(ROOT.rglob('*'))
                if p.is_file() and p.name != 'artefatos_sha256.json'
                and not any(part in ('__pycache__', 'cache_matplotlib') for part in p.parts)}
    json_write(ROOT/'dados/artefatos_sha256.json', manifest)
    print('Relatório e quatro figuras PNG/PDF gerados.', flush=True)
