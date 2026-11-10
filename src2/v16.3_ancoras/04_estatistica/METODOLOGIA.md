# Método e leitura dos resultados

A análise principal utilizou os **25 ensaios** do catálogo congelado de `v16.3_ancoras`: 10 de churn, 9 de wavy annular, 3 de plug e 3 de slug. Foram calculados **300 pares únicos**, sendo 87 intra padrão e 213 inter padrão. A matriz contém 625 células; a diagonal unitária e a metade espelhada não foram contadas nos resumos.

## Entradas e compatibilidade dos recortes

Cada ensaio foi representado pelo mesmo recorte de **2.048 px**, centrado na sua âncora individual, em todas as comparações. Esse perfil já havia sido preparado e auditado para os pares de `02_wavelets`. Os contextos nativos de 4.096 px de churn/wavy annular permaneceram nas imagens de âncoras. Não houve mudança de centro, reamostragem, preenchimento artificial do sinal medido nem escolha de janela em função de Pearson ou WCO.

Os pixels das linhas **18–71, inclusive**, foram conferidos novamente. A média por coluna foi invertida por `255 − intensidade`, teve a tendência linear removida e foi centralizada e dividida pelo desvio padrão populacional. Essa normalização numérica do sinal é distinta da intensidade média mostrada nas imagens de âncoras.

Foram mantidos 200 comprimentos de onda igualmente espaçados em log, entre **50 e 1.024 px**, Morlet com ω₀ = 6, amostragem de 1 px e uma realização (**K = 1**). As unidades são espaciais do panorama montado; não representam milímetros nem frequência temporal calibrada.

## Wavelets e descritor de Pearson

A autocorrelação e a correlação cruzada bilaterais utilizaram o divisor **N**, estimador enviesado adotado na análise anterior:

$$\widehat R_{xy}[m]=\frac{1}{N}\sum_{n\;\mathrm{válidos}}x[n+m]\,\overline{y[n]}.$$

A alternativa de Schuck–Bodmann foi mantida, com a adaptação explícita desse estimador. O livro usa N − |m|; a auditoria anterior documentou por que essa opção não foi adotada. A WPS foi o **módulo** da CWT da autocorrelação, sem elevar esse módulo novamente ao quadrado. A WCS foi armazenada complexa e exibida em módulo:

$$WPS_x(m,a)=|\mathrm{CWT}\{\widehat R_{xx}\}(m,a)|,$$

$$WCS_{xy}(m,a)=\mathrm{CWT}\{\widehat R_{xy}\}(m,a),\qquad
WCO_{xy}(m,a)=\frac{|WCS_{xy}(m,a)|^2}{WPS_x(0,a)\,WPS_y(0,a)}.$$

Foi usada a conversão λ = 4πa/(6 + √38). Denominadores abaixo do piso relativo de 10⁻¹⁰ seriam recusados; as 200 escalas foram válidas nos 300 pares. A WCO não foi truncada artificialmente em 1. As linhas curvas pontilhadas dos mapas delimitam o suporte geométrico adotado, não uma região de significância estatística.

O Pearson reproduziu o descritor de [`estatistica5`](../../estatistica5/RESUMO_EXECUTIVO.md): a **forma da WPS no lag zero após retirada da tendência log-linear**. Para cada ensaio i:

$$P_{ik}=WPS_i(0,a_k),\qquad
z_{ik}=\log(P_{ik}+10^{-12})-\widehat c_i-\widehat d_i\log(\lambda_k).$$

Os coeficientes c e d foram ajustados por mínimos quadrados sobre as mesmas 200 escalas. O Pearson assinado foi calculado entre os resíduos:

$$r_{ij}=\frac{\sum_k(z_{ik}-\bar z_i)(z_{jk}-\bar z_j)}
{\sqrt{\sum_k(z_{ik}-\bar z_i)^2}\sqrt{\sum_k(z_{jk}-\bar z_j)^2}}.$$

Cada descritor foi fixo por ensaio, independentemente de seu parceiro. **Não foi usado o Pearson da WCO, da curva potência média ou diretamente das intensidades.** Na matriz, roxo indica r negativo, creme indica proximidade de zero e azul indica r positivo; os limites cromáticos foram fixados em −1 e 1, como na referência.

A curva **potência média**, ao lado de cada WPS, continua sendo a média do módulo por escala nos lags com suporte geométrico. Já a **WCO média central**, apresentada somente nas tabelas e no gráfico estatístico, resume as 200 escalas e os 65 lags de −32 a +32 px, todos com suporte válido. Cada célula teve o mesmo peso; portanto, o peso por intervalo em log λ é uniforme. Esse resumo facilita comparar os pares, mas pode esconder faixas locais de alta coerência.

## Controle de sensibilidade e limites

Repetiu-se somente o Pearson com os **panoramas completos**, mantendo a faixa de 50–1.024 px e o mesmo método. Cada panorama completo também gerou um único descritor fixo. Esse controle alterou comprimento e conteúdo espacial simultaneamente; não isolou o efeito de apenas um deles. As matrizes principal e de controle foram conservadas, sem seleção pelo maior contraste.

Em relação à referência, mudaram a amostra e a janela: aqui há 25 ensaios dos quatro padrões do catálogo atual, recortes centrados de 2.048 px e λ até 1.024 px. `estatistica5` usou 26 ensaios, incluindo dispersed bubbles, os primeiros 3.171 px e λ até 3.000 px. Os números das duas análises não devem ser comparados como se só uma variável tivesse mudado.

Os pontos de `correlacoes_dentro_entre.png` são pares únicos; o traço horizontal indica a mediana e o vertical, os quartis 25%–75%. Os pares compartilham aquisições, e as escalas são correlacionadas. Por isso, não foram produzidos valores-p, intervalos de confiança ou alegações de acurácia tomando esses pontos como réplicas independentes. Os perfis de `perfis_comparados.png` mostram cada ensaio em linha fina e a mediana por padrão em linha grossa; a divisão da WPS pelo máximo é apenas visual.

Nos nomes de arquivos inter padrão, **favoravel** indica predominância de baixa coerência junto ao alinhamento; **desfavoravel** destaca faixas de coerência elevada que dificultam distinguir os padrões. Os mapas completos foram inspecionados; os rótulos são descritivos e não obedecem a um limiar estatístico. Picos restritos ou deslocados podem existir nos pares favoráveis. Em `02_wavelets`, favorável mantém o sentido anterior: semelhança dentro do mesmo padrão.

## Verificação e arquivos

Os perfis foram confrontados com somas temporais independentes da Morlet; o controle completo foi verificado também por CWT via FFT. A regressão foi conferida por uma segunda formulação. Os 300 coeficientes de cada condição foram comparados com SciPy e NumPy. Também foram verificados r = −1, 0 e 1, recusa de descritores constantes, simetria, diagonal unitária e semidefinição positiva das matrizes. As funções matemáticas de `estatistica5` foram confrontadas diretamente, sem importar sua base antiga.

- [Matriz principal](matriz_pearson_todos.png) e [matriz do controle](matriz_pearson_completo.png).
- [Índice dos ensaios](INDICE_ENSAIOS.md), [300 pares e seus coeficientes](dados/pares_pearson_wco.csv), [perfis espectrais](dados/perfis_espectrais.csv) e [resumos por padrão](dados/resumo_por_padrao.csv).
- [Verificação estatística](../intermediario/estatistica/verificacao.json), [verificação dos pares inter](../intermediario/lote_inter/verificacao.json) e [auditoria anterior das âncoras e pares intra](../intermediario/auditoria_lote_completo/relatorio.txt).
- Scripts: [comparações inter](../scripts/run_full_inter.py), [estatísticas](../scripts/statistics_anchors.py) e [registro da revisão visual](../scripts/record_inter_review.py).

As precomputações individuais permaneceram em `intermediario/precalculo_wavelets`; os novos mapas cruzados ficaram em `intermediario/lote_inter/mapas`. Os descritores, matrizes e coeficientes de ajuste reutilizáveis ficaram em `intermediario/estatistica/descritores.npz`. As pastas finais contêm PNGs e tabelas/documentação, sem PDFs. As imagens foram dimensionadas para inserção com 16 cm de largura em uma página A4 retrato.
