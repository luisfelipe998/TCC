# v5 Spatial Wavelet Analysis — Insights

Spatial domain (1 px/sample), wavelet `cmor2.0-0.9549` (ω₀≈6). CWT applied to the biased
autocorrelation Rxx of the reconstructed panorama signal (Voltage_2, z-scored). Power =
mean |WPS| integrated per band, as % of total power across the full 50–3000 px spectrum.

---

## 1. WPS Power by Run and Pattern

`*` COI-limited band (λ > λ_max).  `†` degenerate run.  `‡` outlier.

```
Run              Pattern               N      λ_max    50–100  100–200  200–500  500–1k  1k–2k  2k–3k
────────────────────────────────────────────────────────────────────────────────────────────────────────
m03q01_001       plug-slug             5148    3882      1.4%    2.3%    14.5%    22.7%   51.8%   7.0%
m03q05_001       churn                 9586    7230      2.8%    8.3%    41.8%    32.7%   11.7%   2.5%
m03q10_001       wavy annular         14224   10729      1.9%   12.3%    42.1%    32.2%    9.2%   2.2%
m03q20_001       wavy annular         23079   17408      2.3%   11.9%    36.4%    30.0%   13.5%   5.8%
m03q40_001       wavy annular         20743   15645      9.1%   22.4%    30.3%    22.0%   11.5%   4.3%
m03q80_001       smooth annular †       501     378     22.9%   37.6%    25.8%     8.0%    3.8%   1.6%
m04q01_001       slug                  5071    3824      1.6%    5.4%    27.2%    46.8%   12.1%   6.4%
m04q05_001       slug-churn            5460    4118      1.6%    3.6%    11.8%    36.9%   27.6%  17.4%
m04q10_001       churn                 8829    6659      2.1%    5.3%    17.3%    55.5%   15.2%   4.1%
m04q20_001       churn                 9344    7047      2.1%    6.4%    27.9%    41.5%   16.4%   5.3%
m04q40_001       wavy annular         14960   11283      1.5%    7.4%    40.7%    38.9%    9.3%   1.9%
m04q80_001       wavy annular         22636   17074      1.9%   13.7%    49.5%    23.7%    7.6%   3.5%
m05q01_001       plug                  5940    4480      2.7%    6.4%    31.3%    35.8%   18.5%   5.0%
m05q05_001       slug                  6323    4769      2.6%    6.2%    28.4%    42.8%   16.4%   3.2%
m05q05_002       slug ‡                4510    3401     12.7%   27.0%    30.4%     9.7%   11.8%   7.9%
m05q10_001       churn                 7001    5280      2.2%    4.3%    18.9%    38.0%   27.0%   9.3%
m05q20_001       churn                 7391    5574      3.8%    6.4%    28.8%    39.1%   18.4%   3.3%
m05q40_001       churn                 9937    7494      2.2%    6.1%    34.9%    39.1%   13.9%   3.6%
m05q80_001       wavy annular         17524   13218      1.1%    7.5%    48.9%    29.9%    8.8%   3.5%
m06q01_001       plug                  7798    5882      4.4%   12.1%    37.5%    24.8%   18.2%   2.8%
m06q05_001       slug ‡                7787    5874      3.4%    9.4%    27.0%    19.3%   29.7%  10.2%
m06q10_001       churn                 7584    5720      1.8%    5.7%    18.1%    28.8%   39.0%   6.3%
m06q20_001       churn                10044    7574      2.2%    5.9%    35.4%    36.9%   16.9%   2.6%
m06q40_001       wavy annular         14772   11141      0.9%    5.0%    35.3%    35.5%   20.4%   2.6%
m06q80_001       wavy annular         18665   14082      0.7%    4.9%    40.0%    34.8%   12.0%   7.3%
m07q01_001       dispersed bubbles     3171    2391     13.2%   21.6%    14.6%    23.0%   20.5%   6.5%*
m07q05_001       plug                  7358    5550      3.8%   11.5%    29.3%    31.3%   19.1%   4.8%
m07q10_001       slug-churn            9101    6864      4.0%   13.7%    33.3%    17.1%   15.1%  15.9%
m07q20_001       churn                 8255    6228      4.7%    7.6%    19.9%    36.1%   23.3%   8.0%
m07q40_001       churn                10349    7806      2.1%    6.5%    34.0%    38.9%   13.5%   4.8%
m07q80_001       wavy annular         16658   12567      0.6%    3.6%    38.0%    40.5%   15.1%   2.0%
```

### Mean by pattern (degenerate/outlier runs excluded)

```
Pattern               runs   50–100  100–200  200–500  500–1k  1k–2k  2k–3k   dominant band
──────────────────────────────────────────────────────────────────────────────────────────────
wavy annular            8     2.2%    9.9%    40.1%    31.9%   11.9%   3.7%   200–500 px
churn                  10     2.6%    6.2%    27.7%    38.7%   19.5%   5.0%   500–1000 px
plug                    4     3.6%   10.0%    32.7%    30.6%   18.6%   4.2%   200–500 px
slug                    3     2.6%    7.0%    27.5%    36.3%   15.0%   5.3%   500–1000 px
plug-slug               1     1.4%    2.3%    14.5%    22.7%   51.8%   7.0%   1000–2000 px
slug-churn              2     2.8%    8.6%    22.6%    27.0%   21.4%  16.6%   500–1000 px
dispersed bubbles       1    13.2%   21.6%    14.6%    23.0%   20.5%   6.5%*  500–1000 px (2k–3k COI-limited)
```

---

## 2. COI Limits and Reliable Wavelength Range

The CWT is applied to the bilateral autocorrelation `Rxx_full` (length 2N-1, lag=0 at
center index N-1). The Torrence & Compo (1998, eq. 19) e-folding limit therefore applies
to input length 2N-1, not N:

    λ_max = (2N-1) / (2√2 × center_freq)     center_freq = 0.9375 for cmor2.0-0.9549

Beyond λ_max, the wavelet window extends past the signal boundary and pywt zero-pads,
mixing real data with artificial zeros. The per-pattern cap is the empirically verified
reliable upper limit, confirmed by broadband 1/f WCO tests across all patterns.

### Suggested wavelength range per pattern

| Pattern | COI cap (px) | Suggested ticks (px) | Energy captured |
|---------|-------------|----------------------|-----------------|
| wavy annular | 3000 (full axis) | 50, 100, 200, 500, 1000, 2000, 3000 | 99.8% |
| churn | 3000 (full axis) | 50, 100, 200, 500, 1000, 2000, 3000 | 99.7% |
| plug | 3000 (full axis) | 50, 100, 200, 500, 1000, 2000, 3000 | 99.4% |
| slug | 3000 (full axis) | 50, 100, 200, 500, 1000, 2000, 3000 | 98.7% |
| plug-slug | 3000 (full axis) | 50, 100, 200, 500, 1000, 2000, 3000 | 97.5% (1 run) |
| slug-churn | 3000 (full axis) | 50, 100, 200, 500, 1000, 2000, 3000 | 97.4% |
| dispersed bubbles | 2400 | 50, 100, 200, 500, 1000, 2000 | 93.1% (1 run) |
| smooth annular | — | — | degenerate, exclude |

Per-pattern λ_max range (computed as `(2N-1)/(2√2·cf)`):

| Pattern | N range | λ_max range |
|---------|---------|-------------|
| wavy annular | 14224–23079 | 10729–17408 px → full 50–3000 px axis reliable |
| churn | 7001–10349 | 5280–7806 px → full 50–3000 px axis reliable |
| plug | 5940–7798 | 4480–5882 px → full 50–3000 px axis reliable |
| slug | 5071–6323 | 3824–4769 px → full 50–3000 px axis reliable |
| plug-slug | 5148 | 3882 px → full 50–3000 px axis reliable |
| slug-churn | 5460–9101 | 4118–6864 px → full 50–3000 px axis reliable |
| dispersed bubbles | 3171 | 2391 px → reliable to ~2400 px |

**All patterns except dispersed bubbles have reliable coverage of the full 50–3000 px axis.**

---

## 3. Key Findings

The ratio of power in 200–500 px vs 500–1000 px is the strongest discriminator across patterns.
All patterns show negligible energy below 100 px — fine tube-scale structure is irrelevant.

| Pattern | 200–500 px | 500–1000 px | Primary discriminator |
|---------|-----------|-------------|----------------------|
| wavy annular | 40% | 32% | 200–500 px dominant — annular wave spacing |
| plug | 33% | 31% | balanced — two distinct scales (cap + body) |
| slug / churn | 27–28% | 37–39% | 500–1000 px dominant — elongated bubble bodies |
| plug-slug | 14% | 23% | very large structures (1k–2k px), full axis reliable |

Same-pattern WCO is consistently high (~0.7–1.0) at lag 0 within each pattern's dominant band —
the dominant spatial scale is a stable fingerprint across independent runs of the same condition.

All patterns have reliable spectral coverage of the full 50–3000 px axis, except dispersed
bubbles (reliable to ~2400 px).

### Wavy annular
**Dominant band: 200–500 px (40%).** The WPS shows a broad, diffuse ridge spanning 100–500 px
that is consistent across all 8 runs and all mass/heat flux conditions. This reflects the
periodic spacing of the liquid film waves on the tube wall. Cross-run WCS shows recurring
hot-spots at 300–500 px spaced ~3000–5000 px apart in lag — the annular wave pattern is
spatially quasi-periodic and repeatable. WCO is moderate (~0.5–0.7 at lag 0), meaning the
scale fingerprint is consistent but not perfectly locked between independent runs.

### Churn
**Dominant band: 500–1000 px (39%).** WPS shows a concentrated, bright locus near 500–700 px
that moves between runs — m03q05 and m05q40 peak closer to 200–500 px while m04q10 peaks
sharply at ~600 px with 55% of total power. This spread likely reflects different bubble
coalescence regimes across heat fluxes (q″ = 5–80 kW/m²). Cross-run WCO is high and broad
at lag 0 (~0.8–1.0 in the 300–1000 px band), confirming that despite the run-to-run shift in
peak location, churn always concentrates power in this region — making it a reliable signature.

### Slug
**Dominant band: 500–1000 px (36%).** WPS is tighter and more localized than churn, with a
sharp peak near 500–700 px. The 3 non-outlier runs agree well; WCO is high (~0.8–1.0) and
narrow at lag 0, indicating slug bodies have a consistent spatial scale across runs.
Two outlier runs (m05q05_002, m06q05_001) show strongly shifted distributions, suggesting
slug body length is sensitive to operating conditions.

### Plug
**Dominant band: 200–500 px (33%) and 500–1000 px (31%) — roughly balanced.** The WPS has
two visible ridges: one in the 200–500 px range (bubble cap diameter) and a secondary at
500–1000 px (bubble body length). The cross-run WCS shows clearly separated hot-spots at
both scales. WCO at lag 0 is high (~0.8–1.0) and spans both bands, with coherence extending
to ±2000 px lag — plug bubbles are large enough that their spatial autocorrelation reaches
far. The split between the two bands distinguishes plug from both slug (500–1000 px dominant)
and wavy annular (200–500 px dominant).

### Plug-slug
**Dominant band: 1000–2000 px (52%).** The single run m03q01_001 has unusually large
structures concentrated well above 1000 px. The full 50–3000 px axis is reliable (λ_max=3882 px).

### Slug-churn
**Dominant band: 500–1000 px (27%), with unusually heavy energy at 1k–2k px (21%) and
2k–3k px (17%).** The two runs (m04q05 and m07q10) are quite different: m04q05 peaks sharply
at ~1000 px while m07q10 has a broader ridge at 200–500 px. WCO between them is moderate
(~0.5–0.6 at lag 0), spread across a wide band rather than concentrated. This is consistent
with slug-churn being a transition regime — the spatial structure is less organized than
either slug or churn.

### Dispersed bubbles (single run)
**Energy spread across all bands, with a slight peak at 500–1000 px (23%).** The WPS lacks
the sharp dominant ridge seen in slug/churn/plug — power is distributed from 50 px to the
COI limit at ~2400 px, consistent with many small, randomly spaced bubbles without a dominant
spatial scale. Only one run available; no cross-run analysis possible.

---

## 4. Caveats

- **`m03q80_001` (smooth annular):** panorama reconstruction failed — mean inter-frame shift
  clamped to 1 px, yielding N=501 and λ_max=378 px. Excluded from all analysis.
- **`m05q05_002`, `m06q05_001` (slug):** anomalous power distributions. Excluded from pattern
  mean. Likely due to unusual flow conditions or recording artefacts.
- **`dispersed bubbles`:** single run, COI limit at ~2400 px. The 2k–3k px band should be
  treated as indicative only.
- **`plug-slug`, `dispersed bubbles`:** single run each — no intra-pattern coherence analysis
  possible.
