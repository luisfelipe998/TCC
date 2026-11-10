


## 1. As imagems por exemplo de m03q05 são naturalmente mais claras que a m04q05. Como normalizar?

![m03q05](image-1.png)
![m04q05](image.png)

## 2. Tamanho da filmagem e a dominio do tempo ao inves do dominio do espaço. F_min possivel não é suficente, o que da pra fazer?

| Cenário | Padrão | Tempo total (s) | FPS câmera | fs nativa (Hz) | Nº imagens | Deslocamento médio (px/frame) | f_min wavelet (Hz) |
|---------|--------|----------------|------------|----------------|------------|-------------------------------|-------------------|
| m03q01_001 | plug-slug | 0.1002 | 5000 | 51377.2 | 501 | 10.28 | 9.98 |
| m03q05_001 | churn | 0.1002 | 5000 | 95668.7 | 501 | 19.13 | 9.98 |
| m03q10_001 | wavy annular | 0.1002 | 5000 | 141956.1 | 501 | 28.39 | 9.98 |
| m03q20_001 | wavy annular | 0.1253 | 4000 | 184263.5 | 501 | 46.07 | 7.98 |
| m03q40_001 | wavy annular | 0.1253 | 4000 | 165612.8 | 501 | 41.40 | 7.98 |
| m03q80_001 | smooth annular | 0.1253 | 4000 | 4000.0 | 501 | 1.00 | 7.98 |
| m04q01_001 | slug | 0.0716 | 7000 | 70852.3 | 501 | 10.12 | 13.97 |
| m04q05_001 | slug-churn | 0.0716 | 7000 | 76287.4 | 501 | 10.90 | 13.97 |
| m04q10_001 | churn | 0.1002 | 5000 | 88113.8 | 501 | 17.62 | 9.98 |
| m04q20_001 | churn | 0.1002 | 5000 | 93253.5 | 501 | 18.65 | 9.98 |
| m04q40_001 | wavy annular | 0.1002 | 5000 | 149301.4 | 501 | 29.86 | 9.98 |
| m04q80_001 | wavy annular | 0.1002 | 5000 | 225908.2 | 501 | 45.18 | 9.98 |
| m05q01_001 | plug | 0.0716 | 7000 | 82994.0 | 501 | 11.86 | 13.97 |
| m05q05_001 | slug | 0.0716 | 7000 | 88345.3 | 501 | 12.62 | 13.97 |
| m05q05_002 | slug | 0.0716 | 7000 | 63014.0 | 501 | 9.00 | 13.97 |
| m05q10_001 | churn | 0.0716 | 7000 | 97818.4 | 501 | 13.97 | 13.97 |
| m05q20_001 | churn | 0.0716 | 7000 | 103267.5 | 501 | 14.75 | 13.97 |
| m05q40_001 | churn | 0.0716 | 7000 | 138840.3 | 501 | 19.83 | 13.97 |
| m05q80_001 | wavy annular | 0.0716 | 7000 | 244846.3 | 501 | 34.98 | 13.97 |
| m06q01_001 | plug | 0.0716 | 7000 | 108954.1 | 501 | 15.56 | 13.97 |
| m06q05_001 | slug | 0.0716 | 7000 | 108800.4 | 501 | 15.54 | 13.97 |
| m06q10_001 | churn | 0.0716 | 7000 | 105964.1 | 501 | 15.14 | 13.97 |
| m06q20_001 | churn | 0.0716 | 7000 | 140335.3 | 501 | 20.05 | 13.97 |
| m06q40_001 | wavy annular | 0.0716 | 7000 | 206395.2 | 501 | 29.49 | 13.97 |
| m06q80_001 | wavy annular | 0.0716 | 7000 | 260788.4 | 501 | 37.26 | 13.97 |
| m07q01_001 | dispersed bubbles | 0.0716 | 7000 | 44305.4 | 501 | 6.33 | 13.97 |
| m07q05_001 | plug | 0.0716 | 7000 | 102806.4 | 501 | 14.69 | 13.97 |
| m07q10_001 | slug-churn | 0.0716 | 7000 | 127159.7 | 501 | 18.17 | 13.97 |
| m07q20_001 | churn | 0.0716 | 7000 | 115339.3 | 501 | 16.48 | 13.97 |
| m07q40_001 | churn | 0.0716 | 7000 | 144596.8 | 501 | 20.66 | 13.97 |
| m07q80_001 | wavy annular | 0.0716 | 7000 | 232746.5 | 501 | 33.25 | 13.97 |

> **f_min** é a menor frequência que a análise wavelet consegue resolver, dada por `f_min = 1 / T_total`, onde `T_total = N_imagens / FPS` é a duração total da filmagem. A wavelet precisa de pelo menos um ciclo completo dentro da janela do sinal para detectar aquela frequência. Com 100 ms de gravação, `f_min ≈ 10 Hz` — ou seja, padrões mais lentos como slug (~2,8 Hz) ou plug (<1 Hz) exigiriam filmagens de pelo menos 0,36 s e 1 s, respectivamente.
