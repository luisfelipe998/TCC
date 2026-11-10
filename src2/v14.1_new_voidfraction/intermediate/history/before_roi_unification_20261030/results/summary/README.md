# Did the fixes help?

**The column-order bug is removed. Improved physical void accuracy is not established.**

For the 25 eligible acquisitions, maximum absolute α change is 0.000976. Flip-rate changes range from -7.26% to 6.76%. Among 12 fixed classifier comparisons, balanced accuracy improves in 1, worsens in 4, and stays equal in 7.

Threshold-dependent hull overlap remains a limitation; the primary threshold was not retuned. No lower/higher flip-rate is interpreted as intrinsically more accurate.

- [Full discussion](../../findings/fix_effectiveness.md)
- [Machine-readable results](fix_effectiveness.json)
- [Every acquisition, including frozen-reference control](fix_comparison_per_run.csv)
- [Group medians before and after](fix_comparison_by_pattern.csv)
- [All classifier models, validations and variants](classifier_comparison.csv)
- [All threshold hull intersections in both versions](threshold_hull_intersections.csv)
- [Exact threshold coordinates in both versions](threshold_plot_data.csv)
- [English summary](../../findings/EXECUTIVE_SUMMARY.md) · [Resumo executivo](../../findings/EXECUTIVE_SUMMARY_PT_BR.md)
