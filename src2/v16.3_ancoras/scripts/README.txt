v16.3 — LOTE INTRA PADRÃO COMPLETO

Escopo atual: 25 aquisições, uma âncora individual por aquisição e todos os
87 pares intra padrão sem repetição: plug 3 ensaios/3 pares, slug 3/3,
churn 10/45 e wavy annular 9/36. Dispersed bubbles, transições e aquisições
terminadas em _002 permanecem fora do escopo dos quatro padrões definidos.
As comparações inter padrão ainda não foram executadas.

Expansão a partir do piloto validado, na raiz TCC:
  python3 src2/v16.3_ancoras/scripts/run_full_intra.py preparar
  python3 src2/v16.3_ancoras/scripts/run_full_intra.py ancoras
  [inspecionar as 25 figuras; registrar lote_intra/inspecao_ancoras.json]
  python3 src2/v16.3_ancoras/scripts/run_full_intra.py calcular
  python3 src2/v16.3_ancoras/scripts/run_full_intra.py galerias
  [inspecionar os 87 WCO; registrar lote_intra/rotulos_visuais.json]
  python3 src2/v16.3_ancoras/scripts/run_full_intra.py figuras
  python3 src2/v16.3_ancoras/scripts/run_full_intra.py verificar

A preparação do lote já foi executada; não repeti-la sobre os mapas congelados.
O comando preparar recusa reiniciar um lote já existente. As etapas seguintes
podem ser chamadas separadamente, respeitando os hashes e as inspeções exigidas.
As evidências visuais são registros do agente; não são aprovações automáticas.
As galerias servem à conferência e ficam apenas em intermediario/lote_intra.

Reprodução do piloto anterior (o reset abaixo exige solicitação de limpeza):
  python3 src2/v16.3_ancoras/scripts/reset_outputs.py --executar
  python3 src2/v16.3_ancoras/scripts/run_pilot.py --etapa ancoras
  [inspecionar as oito figuras e registrar inspecao_visual_piloto.json]
  python3 src2/v16.3_ancoras/scripts/run_pilot.py --etapa preparar
  python3 src2/v16.3_ancoras/scripts/run_pilot.py --etapa wavelets
  [inspecionar as quatro prévias e registrar rotulos_visuais.json]
  python3 src2/v16.3_ancoras/scripts/run_pilot.py --etapa figuras
  python3 src2/v16.3_ancoras/scripts/run_pilot.py --etapa verificar
  python3 src2/v16.3_ancoras/scripts/verify_rebuild.py

reset_outputs.py só deve ser executado quando a limpeza for solicitada.
Apaga 01_ancoras, 02_wavelets e intermediario, incluindo caches, históricos,
auditorias e comparações antigas. Preserva scripts e as fontes originais.
O registro da limpeza é criado em intermediario/reinicio.json.
Nenhum arquivo calculado da execução anterior é copiado para a execução nova.

A inspeção visual é registrada como evidência humana/do agente, sem fingir uma
aprovação automática. Os oito IDs e hashes dos novos PNGs devem constar em
inspecao_visual_piloto.json antes das wavelets. rotulos_visuais.json contém
padrão, ensaios, favoravel/desfavoravel, justificativa, data e hash do mapa novo.
Os rótulos são descritivos e não mudam âncoras, janelas, parceiros ou cálculos.

PASTAS
01_ancoras: 25 PNGs de panorama completo, contexto nativo e média por coluna.
02_wavelets: 87 PNGs com contextos A/B de 2048 px, sinais, WPS, potência média,
WCS e WCO. Um arquivo por par; favoravel/desfavoravel é só orientação visual.
03_conclusoes: permanece vazia, conforme solicitado.
intermediario: dados, prévias, caches, inspeções e verificações. O catálogo e as
figuras anteriores do piloto estão em lote_intra/piloto_anterior. Os quatro
mapas nativos do piloto continuam em intermediario/mapas; não misturá-los com
os 87 mapas comparáveis novos em intermediario/lote_intra/mapas.

As imagens mantêm o formato compacto em retrato, PNG a 300 dpi, para inserção
com 16 cm de largura, até 24 cm de altura e fonte mínima de 8 pt. Não são páginas
A4 rasterizadas e não há PDFs. Padrões aparecem em inglês; demais textos são
curtos e em português. WPS/WCS usam viridis; WCO usa magma. Favoravel/desfavoravel
aparece somente no nome do arquivo. O eixo vertical dos mapas é comprimento
de onda lambda (px); a média da WCO fica apenas nos dados, fora das figuras.

SELEÇÃO ÓPTICA
Plug/slug: pixels <60 nas linhas 18–71 inclusivas, fração por coluna >=25%,
mediana de 5 colunas e fechamento de 13. Intervalos de 96–1536 px, fora das
bordas. Aprovação morfológica herdada da revisão individual da v16.2; ausência
de revisão não é aprovação. A detecção e as coordenadas são recalculadas do
panorama original e conferidas contra as anotações usadas como entrada.

Exige correspondência única nos limiares 55 e 65 (IoU >=0,5). Ordena pela menor
variação máxima das fronteiras; em empate, maior contraste da extremidade
menos definida; depois, menor coordenada. O contraste usa a mediana externa
menos a interna em faixas de 32 colunas do perfil médio. A âncora é o centro
inteiro do intervalo, com fim exclusivo. O ranking não consulta wavelets.

Churn/wavy annular: vales na média suavizada com sigma=16; proeminência mínima
max(3; 0,25*(P95-P5)), largura 16–512 px, separação 128 px e janela de
proeminência 1025 px. O mesmo vale precisa corresponder unicamente em sigma=12
e sigma=20, até 16 px. Ordena por maior proeminência em sigma=16, menor variação
posicional e menor coordenada. A suavização serve ao detector; o preparo da
wavelet parte da média original recortada.

As janelas nativas preferidas são 2048 px para plug/slug e 4096 px para
churn/wavy annular. Exigem margem de aquisição de 16 px e suporte óptico que
inclui a união das fronteiras/bases das três variantes. Além desse suporte,
usam margens de 256 px para eventos e 128 px para vales. Se a janela preferida
inviabilizar os candidatos, é reduzida em passos de 256; mínimo 1024 px.

PRÉ-COMPUTAÇÃO REUTILIZÁVEL
prepare_comparisons.py é limitado ao piloto por padrão; o lote completo chama
explicitamente escopo=completo. Agora há 25 ensaios, 50 perfis nomeados,
44 contextos numéricos distintos e duas grades Morlet. Os 12 contextos numéricos
do piloto foram reutilizados sem reescrita; foram acrescentados outros 32.
Os perfis nativo e comparavel compartilham arquivos quando seus recortes são
idênticos. A preparação individual nunca calcula pares inter.

nativo: 2048 px em plug/slug e 4096 px em churn/wavy annular.
comparavel: 2048 px para todos, com as mesmas 200 escalas entre 50 e 1024 px,
resolução de 1 px e lags de -2047 a +2047 px. A janela comum é explícita na
configuração; a expansão não pode alterá-la silenciosamente.

O centro escolhido no contexto nativo é preservado entre perfis. Os diagnósticos
são recalculados para cada janela: costuras, vizinhos e repetibilidade do ranking
não são herdados da janela maior. O catálogo distingue ranking nativo e posição
no ranking do contexto atual. Repetibilidade da seleção completa não é o mesmo
que encontrar o vale correspondente nas três suavizações.

Cada contexto é recortado diretamente da média das linhas 18–71 antes do preparo.
São guardados: média recortada, sinal preparado, autocorrelação, CWT individual
complexa (WPS em módulo), potência média, comprimentos de onda, escalas, lags,
suporte geométrico e núcleos Morlet. WCS e WCO dependem do parceiro e ficam para
a comparação solicitada. Foram calculados os 87 pares intra comparáveis.

intermediario/preparacao_comparacoes contém seleção, contextos por perfil,
congelamento anterior aos espectros, especificações esperadas e índice.
intermediario/precalculo_wavelets contém os arrays imutáveis reutilizáveis.
O plano registra os 300 pares possíveis entre os 25 ensaios (87 intra e
213 inter). A verificação da compatibilidade dos 213 pares inter usa somente
metadados: não executa suas correlações cruzadas, WCS ou WCO.

Cada NPZ do lote guarda WCS complexa, WCO, máscara de denominador válido e
metadados que apontam para as duas entradas individuais. WPS e sinais não são
duplicados em cada par. run_full_intra.read_pair reconstitui o conjunto para
plotagem a partir desses arquivos, sem recalcular a comparação. Os hashes do
catálogo, código, mapas e PNGs vinculam cada resultado aos respectivos dados.

API compartilhada pelas futuras comparações intra e inter:
  from prepare_comparisons import load_catalog, comparison_spec
  from wavelet_cache import pair_from_autos
  catalogo = load_catalog()
  comparacao, contextos, individuais = comparison_spec(
      catalogo, "m05q01_001", "m03q05_001", profile="comparavel")

Até aí não há cálculo do par. pair_from_autos(*individuais) executará a WCS/WCO
quando a análise correspondente for solicitada. Use o perfil comparavel tanto
na intra quanto na inter para comparar resultados no mesmo domínio. O perfil
nativo permite examinar as escalas mais longas e deve ser interpretado à parte.

Preparação individual completa usada nesta expansão autorizada:
  python3 src2/v16.3_ancoras/scripts/prepare_comparisons.py --escopo completo
Acrescenta os dados individuais restantes e reaproveita os arquivos atuais se
fonte, centro, recorte, método e
implementação coincidirem. O nome do parceiro e o escopo do lote não fazem parte
da identidade do cache individual. Repreparar substitui o congelamento: não
fazê-lo sobre pares existentes sem regenerar seus vínculos.

PROTEÇÕES CONTRA DADOS DESATUALIZADOS
A identidade do cache inclui ensaio, fonte, centro, início/fim, linhas, janela,
grade, preparo e hashes do núcleo, do leitor de entrada e da implementação.
O carregador verifica arrays, especificações e vínculos contra o congelamento;
um arquivo íntegro de outro ensaio não é aceito. Mudanças nos arquivos após a
primeira leitura invalidam a reutilização em memória. Metadados retornados não
alteram o objeto interno compartilhado. As fontes e o código também são conferidos.

O método principal exige explicitamente estimador biased, K=1 e omega0=6. Outros
valores são rejeitados, não ignorados. Comprimentos fracionários, ímpares ou fora
do passo de 256 são rejeitados antes de recortar; faixas de linhas inválidas não
são truncadas silenciosamente. Grades de 2048 e 4096 não podem ser misturadas.

FÓRMULAS E LIMITES
wavelet_core.py permanece byte a byte igual ao núcleo da v16.2. O preparo depois
do recorte usa 255-I, remoção da tendência linear/média e desvio padrão unitário.
Correlação com divisor N; Morlet omega0=6. São as definições alternativas de
Schuck–Bodmann com K=1, adaptadas para correlação enviesada:
WPS = |CWT(Rxx)|; WCS = CWT(Rxy);
WCO = |WCS|² / (WPSx(lag0)*WPSy(lag0)).

Schuck usa divisor N-|lag| nas Eqs. 32.21–32.23. A adaptação com N foi mantida
porque conserva a estrutura de energia não negativa e o limite unitário da WCO
com este núcleo. "Não enviesado" para a correlação não assegura uma WCO melhor.
A curva potência média é a média da WPS nos lags com suporte geométrico, sem novo
quadrado; não é a média sobre realizações do capítulo nem potência física.

A WCO descreve o contexto inteiro e o eixo horizontal é lag. Não localiza
exclusivamente a estrutura central. Coerência alta, especialmente nas escalas
longas, não é significância estatística nem validação de classificação. O suporte
geométrico e o piso numérico do denominador também não são testes estatísticos.
Os panoramas permanecem montagens em pixels; conferir a origem das colunas não
certifica uma reconstrução física em mm ou Hz.

VERIFICAÇÃO
run_full_intra.py verificar confere todos os 44 espectros individuais, os
50 perfis e os 87 pares contra regressão, correlação direta e convolução com
a Morlet temporal. Confere também fontes, recortes, máscaras, limite unitário
da WCO, ausência de pares duplicados/inter, hashes e preservação do cache.
As 112 imagens finais são verificadas quanto a dimensão inserida e fonte.
Resultados: intermediario/lote_intra/verificacao/verificacao.json.

Os verificadores abaixo são específicos do piloto e não devem ser usados para
validar o escopo completo. Seus relatórios anteriores foram preservados.
verify_pilot.py verifica o ranking, as fontes, os recortes, as fórmulas dos pares
e a apresentação. verify_rebuild.py confere todos os doze espectros individuais
e os quatro pares por regressão, correlação direta e convolução com a Morlet
temporal, sem usar os núcleos pré-calculados. Verifica também diagnósticos por
perfil, entradas inválidas, vínculo trocado, corrupção após leitura, reutilização
sem reescrita, tamanho das figuras e criação dos arquivos depois da limpeza.
Os resultados ficam em intermediario/verificacao_reconstrucao.
validate_reusable_analysis.py agora chama essa verificação, sem exigir históricos.

Os scripts dos estudos anteriores foram preservados como código de referência;
seus resultados e históricos foram removidos. analyze_windows.py,
audit_anchor_selection.py e study_anchor_consensus.py dependem de estados
históricos e não fazem parte da reprodução atual. As variantes de estimadores
também não foram executadas novamente nesta reconstrução do método escolhido.
