"""Navigation for the isolated v14.1 correction and full v14 comparison."""
import json
from v14_paths import ROOT, artifact_path, write_markdown
from run_all import STAGES


def write_navigation():
    write_markdown(ROOT/'README.md', '''# v14.1 — void fraction and flip-rate

[Executive summary](EXECUTIVE_SUMMARY.md) · [Resumo executivo](EXECUTIVE_SUMMARY_PT_BR.md) · [Effectiveness of the correction](fix_effectiveness.md)

This version repeats the v14 analysis with the same sample, formulas, thresholds, geometry and classifiers. It corrects the panorama seam overwrite that reordered/duplicated columns, rebuilds the liquid reference under the same rule, and exposes threshold dependence in the main scatter analysis. The primary threshold remains 0.70 and the full lumen remains rows 18:71, D=53.

The original v14 folder, original panoramas and native frames are preserved. A frozen [v14 baseline](intermediate/history/v14_baseline/README.md) and its [hash manifest](intermediate/history/v14_baseline_manifest.json) support the comparison. Historical snapshot documents preserve original bytes and may contain links intended for their original locations.

| Start here | Contents |
|---|---|
| [findings](findings/README.md) | Summaries, method reviews and adversarial checks |
| [results](results/README.md) | Current measurements, classifier predictions and verification |
| [results/summary](results/summary/README.md) | Whether the fixes helped, including unfavorable changes |
| [figures](figures/README.md) | Regenerated figures and sensitivity comparisons |
| [scripts](scripts/README.md) | Ordered reproduction stages |
| [intermediate](intermediate/README.md) | Corrected inputs, masks, features, sensitivity and logs |
| [literature](literature/README.md) | Preserved literature evidence |

![Primary analysis](void_fraction_vs_flip_rate.png)

Void estimates remain conditional optical measurements, without independent physical calibration. Flip-rate is spatial transition density, not a temporal frequency. Observed convex hulls are not classification boundaries.

## Reproduce

From the repository root, with NumPy, SciPy, pandas, Pillow, Matplotlib and PyWavelets installed:

```sh
python3 -B src2/v14.1_new_voidfraction/scripts/run_all.py
```

The runner includes all original analysis stages and native-frame alternatives, regenerates comparisons and reports, stops on failure, and records [logs and status](pipeline_run.json). All inputs are local. Reconstruction requires the original native JPEGs, configurations, markers and pipeline_v11.py. No source data or shared pipeline is overwritten.
''')
    documents={
        'findings/README.md': ('Findings', [
            ('EXECUTIVE_SUMMARY.md','English executive summary'),('EXECUTIVE_SUMMARY_PT_BR.md','Resumo executivo'),
            ('fix_effectiveness.md','Before/after results and whether the fixes helped'),
            ('seam_adversarial_review.md','Adversarial source-order and method-preservation audit'),
            ('post_regeneration_adversarial_review.md','Independent audit of completed results and both summaries'),
            ('deep_method_review.md','Independent segmentation and arithmetic'),('lumen_normalization_review.md','Preserved D53 geometry'),
            ('panorama_void_alternatives.md','Void method and assumptions'),('phase_mask_flip_review.md','Spatial mask flip definitions'),
            ('spatial_flip_review.md','Optical spatial sensitivity'),('flip_rate_review.md','Native temporal controls'),
            ('pattern_identification.md','Morphology classifiers'),('silhouette_pattern_review.md','Silhouette classifiers'),
            ('adversarial_review.md','Phase anchors, reference and crop sensitivity'),('source_audit.md','Literature audit'),
            ('prior_findings_audit.md','Selection and independent optical arithmetic'),('scope_verification.md','Files, scope and links'),
            ('GRAFICOS_PARA_DISCUSSAO_PT_BR.md','Leitura dos gráficos')]),
        'results/README.md': ('Numerical results', [
            ('results/summary/README.md','Correction effectiveness and full before/after comparison'),
            ('panorama_void_per_run.csv','Conditional void per acquisition'),('panorama_void_by_pattern.csv','Void group medians'),
            ('phase_mask_flip_per_run.csv','Mask flip counts and denominators'),('phase_mask_flip_by_pattern.csv','Mask flip group medians'),
            ('flip_rate_per_run.csv','Optical spatial and native temporal measures'),('void_fraction_flip_plot_data.csv','Exact chart coordinates'),
            ('pattern_identification_predictions.csv','Morphology predictions'),('silhouette_pattern_predictions.csv','Silhouette predictions'),
            ('final_verification.json','Arithmetic and fold checks'),('seam_adversarial_checks.json','Adversarial correction checks'),
            ('post_regeneration_adversarial_checks.json','Independent post-summary checks'),('layout_verification.json','Final layout and provenance checks')]),
        'figures/README.md': ('Figures', [
            ('void_fraction_vs_flip_rate.png','Primary scatter at threshold 0.70'),('threshold_sensitivity.png','All five existing thresholds, identical axes'),
            ('fix_comparison.png','Changes per eligible acquisition'),('void_fraction_vs_center_flips.png','Fixed row y=40 comparison'),
            ('void_crop_sensitivity.png','Remove known interior rows at D53'),('panorama_phase_masks.png','Gas mask overlays'),
            ('panorama_liquid_reference_candidates.png','Inferred liquid-reference intervals'),('pattern_panorama_contact_sheet.png','All selected panoramas'),
            ('pattern_confusion_matrices.png','Morphology classification'),('silhouette_pattern_confusion_matrices.png','Silhouette classification'),
            ('adversarial_phase_anchor_patches.png','Internal visual controls'),('adversarial_slug_gap_mask.png','Liquid gap control'),
            ('flip_rate_spatial_vs_native.png','Spatial versus native optical measures'),('flip_rate_native_kymographs.png','Native time records')]),
        'intermediate/README.md': ('Intermediate evidence', [
            ('intermediate/reconstructed_inputs','Corrected panoramas, preserving source-strip order'),('phase_masks','Current binary gas masks'),
            ('intermediate/features','Liquid reference, morphology and silhouette features'),('intermediate/sensitivity','All original perturbations'),
            ('intermediate/exploratory','Original supplemental native-image alternatives'),('intermediate/history/v14_baseline','Unmodified v14 snapshot'),
            ('intermediate/history/v14_baseline_manifest.json','Frozen file hashes'),('pipeline_run.json','Latest pipeline status'),
            ('pipeline_logs','Per-stage logs')]),
    }
    for relative,(title,items) in documents.items():
        write_markdown(ROOT/relative,'# '+title+'\n\n[Back to v14.1](../README.md)\n\n'+'\n'.join(f'- [{description}]({target})' for target,description in items)+'\n')
    stages='\n'.join(f'{i}. [{name}.py]({name}.py) — {description}' for i,(name,description) in enumerate(STAGES,1))
    write_markdown(ROOT/'scripts/README.md','# Reproduction stages\n\n[Back to v14.1](../README.md)\n\n'+stages+'\n\nThe ordered runner is [run_all.py](run_all.py). Numerical estimator and classifier functions are checked against the frozen baseline by [seam_adversarial.py](seam_adversarial.py).\n')
    effect=json.loads(artifact_path('fix_effectiveness.json').read_text())
    write_markdown(ROOT/'results/summary/README.md',f'''# Did the fixes help?

**The column-order bug is removed. Improved physical void accuracy is not established.**

For the {effect['n_eligible']} eligible acquisitions, maximum absolute α change is {effect['max_absolute_alpha_change']:.6f}. Flip-rate changes range from {effect['min_relative_flip_change_percent']:.2f}% to {effect['max_relative_flip_change_percent']:.2f}%. Among {effect['fixed_classifier_comparisons']} fixed classifier comparisons, balanced accuracy improves in {effect['fixed_classifiers_improved']}, worsens in {effect['fixed_classifiers_worsened']}, and stays equal in {effect['fixed_classifiers_unchanged']}.

Threshold-dependent hull overlap remains a limitation; the primary threshold was not retuned. No lower/higher flip-rate is interpreted as intrinsically more accurate.

- [Full discussion](fix_effectiveness.md)
- [Machine-readable results](fix_effectiveness.json)
- [Every acquisition, including frozen-reference control](fix_comparison_per_run.csv)
- [Group medians before and after](fix_comparison_by_pattern.csv)
- [All classifier models, validations and variants](classifier_comparison.csv)
- [All threshold hull intersections in both versions](threshold_hull_intersections.csv)
- [Exact threshold coordinates in both versions](threshold_plot_data.csv)
- [English summary](EXECUTIVE_SUMMARY.md) · [Resumo executivo](EXECUTIVE_SUMMARY_PT_BR.md)
''')
    write_markdown(artifact_path('GRAFICOS_PARA_DISCUSSAO_PT_BR.md'),'''# Gráficos para discussão — v14.1

O [gráfico principal](void_fraction_vs_flip_rate.png) mantém a definição e a apresentação da v14, agora com as entradas corrigidas. As estrelas são médias por aquisição; as envoltórias descrevem apenas a amostra.

O [painel dos cinco limiares](threshold_sensitivity.png) revela a dependência da separação em relação à segmentação. O [gráfico de mudanças](fix_comparison.png) mostra os efeitos por aquisição. Leia os [resultados comparativos](fix_effectiveness.md) antes de afirmar melhora.

A [sensibilidade ao recorte](void_crop_sensitivity.png) mantém D=53 e remove somente linhas internas. Nenhuma destas figuras constitui calibração física de volume, validação externa ou intervalo de confiança.
''')
