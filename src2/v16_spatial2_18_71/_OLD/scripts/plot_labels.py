"""Nomes dos padrões em inglês; demais rótulos em PT-BR.

As chaves preservam os identificadores das tabelas. Na apresentação,
"dispersed bubbles" usa o nome "dispersed bubble" solicitado pelo usuário.
PATTERN_PT mantém o nome da interface existente dos geradores.
"""

PATTERN_PT = {
    'churn': 'churn',
    'dispersed bubbles': 'dispersed bubble',
    'dispersed bubble': 'dispersed bubble',
    'plug': 'plug',
    'slug': 'slug',
    'wavy annular': 'wavy annular',
    'smooth annular': 'smooth annular',
    'plug-slug': 'plug–slug',
    'slug-churn': 'slug–churn',
    'cross-pattern control': 'Controle entre padrões',
}

METHOD_PT = {
    'v5': 'v5 aplicado a 18–71',
    'v5_recomputed': 'v5 aplicado a 18–71',
    'v16_biased': 'v16 enviesado, 18–71',
    'schuck_unbiased': 'Schuck não enviesado, 18–71',
}

VARIANT_PT = {
    'baseline': 'Referência: linhas 18–71',
    'demean_only': 'Remoção apenas da média',
    'center_row': 'Linha 40 (controle fixo)',
    'central_band20': 'Faixa 30–49 (controle fixo)',
    'core_row_median': 'Mediana das linhas 18–71',
    'core_trimmed_mean20': 'Média aparada das linhas 18–71',
    'spatial_gaussian2': 'Filtro gaussiano espacial (σ = 2 px)',
    'spatial_median5': 'Mediana espacial (5 px)',
    'rank_uniform': 'Transformação uniforme por postos',
    'winsor_1pct': 'Winsorização de 1%',
    'v5': 'v5 aplicado a 18–71',
}

WINDOW_PT = {'prefix': 'Início', 'center': 'Centro', 'suffix': 'Fim'}


def pattern_label(name):
    return PATTERN_PT.get(name, name)


def method_label(name):
    return METHOD_PT.get(name, name)


def variant_label(name):
    return VARIANT_PT.get(name, name)
