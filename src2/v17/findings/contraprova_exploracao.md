# Contraprova independente da exploração espectral v17

Revisor: auditor independente de estatistica8. Revisão realizada depois de encerrar `auditoria_estatistica8.md`, sem alterar o script do coordenador. Fontes: `PROTOCOLO.txt`, `scripts/exploracao_espectral.py` e resultados em `resultados/exploracao`.

**O resultado descritivo forte da largura Welch foi reproduzido por um cálculo independente e resistiu às contraprovas de normalização, interpolação de quantis e agregação. A separação perfeita dos valores observados depende da resolução/suavização do estimador: usar um único segmento de 2048 px introduziu pequena sobreposição em 50–512 px. Não foi encontrado erro que explique artificialmente a separação original.**

Isso não é validação externa: os mesmos 19 ensaios e os mesmos panoramas sustentam o resultado e suas sensibilidades. As contraprovas foram definidas depois de conhecer o resultado positivo e antes de executar os novos cálculos; não fornecem p-valores confirmatórios, desempenho preditivo ou novas réplicas.

## Fórmulas examinadas

**Jacobiano da densidade.** Se λ = 1/f, então `|df/dlnλ| = f`. Portanto `PSD(f)/λ = f PSD(f)` é a densidade correspondente em log λ. A operação de `exploracao_espectral.py:190–191` está correta. O sentido invertido entre f e λ não gera massa negativa porque se usa o módulo do Jacobiano. A normalização de `mass()` é uma quadratura trapezoidal nessa variável, coerente com a definição declarada.

**Reponderação da WPS alternativa.** Para uma Morlet estreita em frequência, a CWT da autocorrelação no atraso zero se comporta aproximadamente como `S(f₀)/√a`: o fator √a da transformada combina com uma banda de largura proporcional a 1/a. Como λ é proporcional a a, dividir o perfil por √λ aproxima uma ponderação proporcional a `f₀ S(f₀)`. A interpretação deve continuar aproximada; Morlet e Welch têm janelas, larguras de banda e efeitos de suporte distintos. Não é uma identidade universal entre os dois estimadores e não converte intensidade óptica em energia física das fases.

**Mistura de três janelas.** `exploracao_espectral.py:208–215` primeiro normaliza a distribuição de cada janela na banda e depois calcula a média aritmética das três massas. Dividir por `quadrature = h/Σh` reconstrói um perfil proporcional à densidade dessa mistura. A renormalização posterior recupera exatamente as massas desejadas. Não há dupla aplicação indevida dos pesos. O resultado é largura da mistura de distribuições; não é média das três larguras e não é o espectro de um panorama concatenado. Essa distinção precisa acompanhar o resultado.

**Silhueta e contraste.** `separation():108–127` exclui a diagonal no termo intraclasse, usa a menor distância média às demais classes e balanceia as médias por classe. O contraste de distâncias balanceia classes no termo intra e pares de classes no termo inter. Essas são definições coerentes, desde que nomeadas; não equivalem a dar o mesmo peso a cada par individual.

**Janelas dependentes.** Prefixo, centro e sufixo são definidos sem consultar o conteúdo/rótulo, mas não são réplicas. A mistura produz uma distribuição por aquisição e mantém n = 19 no contraste churn/anular. Registros mais curtos têm sobreposição entre janelas. O protocolo e os resultados tratam isso como sensibilidade dependente, de forma adequada.

## Reconstrução independente do Welch

O script [contraprova_exploracao.py](../scripts/contraprova_exploracao.py) releu os pixels, retirou tendência por ajuste polinomial e calculou manualmente os periodogramas: Hann periódica, segmentos de 1024 px, passo 512, FFT 8192 e normalização `|FFT|² / Σw²`, duplicando somente as frequências positivas internas. A média de três periodogramas reproduziu o Welch original sem chamar `scipy.signal.welch` nem importar funções do coordenador.

Foram comparados **150 perfis**: 25 aquisições × 3 janelas × 2 bandas. O erro absoluto máximo foi **2,44 × 10⁻¹⁵**. Em cada cálculo, a integral discreta da PSD foi conferida contra a média da potência ponderada pela janela, pela identidade de Parseval. As larguras da mistura coincidiram com os CSV do coordenador com erro inferior a 10⁻¹².

O resultado original foi confirmado:

| Faixa λ | Silhueta macro | Menor média de classe | Faixa churn | Faixa anular |
| --- | ---: | ---: | --- | --- |
| 50–512 px | 0,590546 | 0,533726 | 0,980278–1,281469 | 0,529590–0,908251 |
| 50–1024 px | 0,595699 | 0,550448 | 1,080573–1,474044 | 0,693644–1,064524 |

Em ambas, todos os 10 valores churn superaram os nove anulares: 90/90 comparações ordenadas, delta de Cliff = 1. São comparações dos mesmos 19 ensaios; não são 90 observações independentes e não expressam acurácia de um classificador.

## Tentativas de refutação

Todas as variantes abaixo foram executadas após conhecer o resultado principal. A configuração principal permaneceu a definida pelo coordenador; não se escolheu a variante de maior silhueta.

| Contraprova em 50–512 px | Silhueta macro | Menor média de classe | Probabilidade empírica churn > anular |
| --- | ---: | ---: | ---: |
| Original, mistura e CDF interpolada | 0,590546 | 0,533726 | 1,0000 |
| Quantil pela integral do interpolante linear da densidade, 200 pontos | 0,590583 | 0,533785 | 1,0000 |
| Mesmo quantil, grade de 2000 pontos | 0,590537 | 0,533706 | 1,0000 |
| Média dos perfis antes da normalização da banda | 0,610250 | 0,560724 | 1,0000 |
| Hann 1024 com remoção de média por segmento | 0,590899 | 0,534174 | 1,0000 |
| Hann 512, sobreposição 50% | 0,629869 | 0,610008 | 1,0000 |
| Hann 2048, um segmento por janela | 0,425204 | 0,289909 | 0,9778 |

**Normalização e quantis não explicaram o efeito.** O pooling antes da normalização modifica o peso das janelas conforme sua massa na banda, enquanto a mistura original lhes dá peso igual. Ambos conservaram a ordenação completa nesta amostra. A mudança de quadratura e quantil alterou pouco as silhuetas.

**A escolha do estimador importa para a separação perfeita.** Com segmentos de 2048 px, as faixas em 50–512 passaram a churn 0,951095–1,319530 e anular 0,453177–0,990095, com 88/90 comparações ordenadas e delta de Cliff 0,9556. A separação continuou positiva nas duas classes, mas a ausência completa de sobreposição deixou de valer. Essa variante tem menor suavização por médias e maior resolução; não é automaticamente mais correta. Na banda 50–1024, essa mesma variante manteve a ordenação completa e silhueta 0,552422.

**Resolução do Welch.** O segmento principal tem 1024 px. Uma escala de 512 px corresponde a apenas duas oscilações por segmento, e 1024 px a uma. A FFT 8192 fornece interpolação da grade, não oito vezes mais informação independente. Por isso a forma no extremo de longos comprimentos de onda deve ser interpretada com cautela, especialmente sem incerteza espectral ou replicação experimental.

## Objeções concretas encontradas

### CEXP-01 — Baixa relevância para os dados atuais — Pearson da WPS e da WPS reponderada é redundante por construção

`log(P/√λ) = log P − 0,5 log λ`. Como o descritor residual remove justamente uma reta em log λ, a transformação não deve mudar o Pearson de cada janela, exceto pelos efeitos do ε = 10⁻¹². As matrizes realmente diferiram no máximo aproximadamente **4,16 × 10⁻¹¹**.

Assim, a concordância dessas duas linhas na triagem não é corroboração independente. A observação não exige remover as variantes de Hellinger, Wasserstein ou largura: elas medem objetos diferentes e respondem à reponderação. Para a mistura de distribuições, normalizar cada janela antes de misturar altera pesos relativos; a invariância algébrica simples não se aplica necessariamente entre as misturas dos dois métodos.

Evidência: `exploracao_espectral.py:76–81,190–191`; [redundância numérica](../resultados/contraprova_exploracao/redundancia_pearson_reponderacao.csv).

### CEXP-02 — Baixa relevância para os resultados atuais — Falta guarda para perfil sem forma residual

Na versão inspecionada, `residual():76–81` normaliza pela norma sem testar se o resíduo é numericamente nulo. Um perfil constante gerou norma de aproximadamente **1,77 × 10⁻²⁷** e retornou um vetor unitário de arredondamento. Para esse caso, uma correlação de forma residual não possui informação identificável.

**Impacto observado:** nenhum nos perfis reais individuais examinados; a menor norma foi **1,3522**. Tampouco afeta a largura Welch. Recomenda-se sinalizar norma residual quase nula antes de normalizar, em vez de transformar erro de arredondamento em direção unitária. A contraprova preservou o script do coordenador; qualquer correção posterior deve registrar nova versão/hash.

Evidência: [caso degenerado](../resultados/contraprova_exploracao/caso_degenerado_residual.json).

**Resolvido durante a revisão.** O coordenador acrescentou uma rejeição explícita quando a norma residual é menor ou igual a `10⁻¹² max(1, ||log P||)`. A contraprova independente da versão corrigida confirmou `ValueError` para perfil constante e lei de potência exata; os 600 perfis individuais reais arquivados foram aceitos e permaneceram finitos. A alteração ocorreu somente em v17. Registro e hash da versão conferida: [verificação após correção](../resultados/contraprova_exploracao/caso_degenerado_apos_correcao.json).

## Conclusão de revisão

O candidato Welch é mais forte descritivamente que a largura WPS original nas janelas avaliadas, mas não constitui prova geral de separabilidade dos padrões. O efeito sobreviveu a várias tentativas razoáveis de refutação computacional; sua interpretação física continua dependente dos mesmos panoramas, das condições operacionais e da definição do espectro. A ausência de sobreposição é uma propriedade dessas observações sob certos estimadores, não uma fronteira validada para novas aquisições.

Para continuar sem aprendizado de máquina, a largura de uma distribuição espectral agregada por aquisição é um caminho concreto e verificável. Devem permanecer explícitos: banda, estimador, normalização por janela, agregação e variabilidade espacial. Novas aquisições e controle das condições experimentais são necessários para distinguir padrão, operação e resposta óptica.

Artefatos da contraprova: [verificação](../resultados/contraprova_exploracao/verificacao.json), [todas as sensibilidades](../resultados/contraprova_exploracao/contraprovas_welch.csv), [larguras reconstruídas](../resultados/contraprova_exploracao/reproducao_largura_welch.csv), [protocolo posterior ao resultado](../resultados/contraprova_exploracao/protocolo_contraprova.json) e [log](../logs/contraprova_exploracao.log).
