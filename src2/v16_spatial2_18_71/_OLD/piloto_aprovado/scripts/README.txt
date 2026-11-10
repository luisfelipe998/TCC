PILOTO DA REORGANIZAÇÃO

Executar a partir da raiz do projeto:
python3 src2/v16_spatial2_18_71/scripts/gerar_piloto.py

Para atualizar somente a apresentação com os dados existentes:
python3 src2/v16_spatial2_18_71/scripts/gerar_piloto.py --somente-figuras

Os títulos por extenso e os rótulos dos eixos são compartilhados pela função
pair_figure, que aceita pares intra ou entre padrões. O piloto continua restrito
ao único par abaixo até a revisão pelo usuário.

O script processa somente m05q01_001 × m06q01_001 (plug).
O processamento dos demais ensaios depende da revisão deste piloto pelo usuário.

wavelet_core.py foi preservado sem alterações em relação a _OLD/scripts.
extraction_config.py apenas atualizou o caminho da configuração para _OLD/results.
A entrada foi novamente extraída dos PNGs
originais em src2/data_processed e conferida contra os sinais arquivados.

01_panoramas: dois panoramas completos e seus sinais.
02_wavelets: figura do par, com prefixos comuns, sinais, WPS, WCS e WCO.
04_estatistica/dados: métricas e perfis espectrais do piloto.
05_conclusao: descrição do piloto e seus limites.
intermediario: dados numéricos, verificação do arquivo antigo e auditoria do piloto.
_OLD: todo o conteúdo anterior, preservado integralmente.

Metodologia: média das linhas 18–71; 255−intensidade, remoção de tendência linear,
centralização e divisão pelo desvio padrão populacional; no par, corte no prefixo
comum antes do preparo; correlação enviesada / N; Morlet omega0=6; 200 comprimentos
de onda de 50 a 3000 px; plano integral de lags, sem COI nem máscara geométrica.
Potência média: média de abs(CWT(Rxx)) em todos os lags, sem elevar ao quadrado.
Pearson: resíduos de log(WPS no lag zero + 1e-12), após ajuste contra 1 e log(lambda),
usando os panoramas completos de cada aquisição e a grade nativa de 50–3000 px.
