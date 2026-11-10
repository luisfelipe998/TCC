# Confirmed-lumen chart and summary handoff

Completed source changes in `summary_figures.py`, `summary_text.py` and `write_summary.py`. Syntax parsing passed for all three. No pipeline or result-writing script was run by this agent.

- Shared geometry now supplies D=53, full-lumen rows and valid internal crop variants. Primary input rows must carry `lumen_start_row`, `lumen_stop_row_exclusive`, `lumen_diameter_px` and `boundary_source`, with exact agreement with the shared configuration.
- Both scatter axes now extend through alpha=1, with x limits −0.025 to 1.03. Original filenames, the approved Portuguese Y-axis label and strict 25-eligible-point selection remain intact. No ceiling annotations or saturated points were added.
- Plot provenance removes missing-width fields and uses `full_lumen_maximum=1`; all 26 retained acquisitions remain traceable. Chart metadata schema is now 3.
- Crop chart accepts only `INTERNAL_CROP_VARIANTS`, with heights 50/51/52/53 and fixed D=53. It describes removal of known-interior rows, not alternative lumen boundaries. Its dashed diagnostic curve is `(h/53)^2`, reaching 1 for the complete lumen.
- Both summary generators retain seven named methods, formulas/pipelines, embedded figures, the shared numerical table and nearby limitations. They explain the user-supplied complete lumen, the mistaken former D60 normalization and the exact correction `alpha53=alpha60*(60/53)^2`. Full mask implies computed alpha=1, not proof of pure gas. No claims of missing lumen or valid outside-lumen extension remain.
- Historical morphology changes are computed by comparing the two frozen D60 snapshots (`before_deep_method_review` and `before_confirmed_lumen_D53`), so the reported morphology delta does not mix in the D53 normalization. The summary writer also asserts the current alpha rescaling and unchanged mask files against the latter snapshot.
- Actual hull overlap is still calculated from saved polygon vertices; no pattern overlap pair is hardcoded.

Root-owned follow-ups before regeneration:

1. `v14_navigation.py`: update former clipping/missing-width prose and list the shared geometry module.
2. `findings/GRAFICOS_PARA_DISCUSSAO_PT_BR.md`: the crop recommendation should say **removal of known-interior rows with D=53**, and the proposed per-column pipeline should explicitly use **D=53**. Its existing relevant recommendations and status as suggestions remain suitable.
3. The last direct collaboration message was unreadable ciphertext; if it contained another requirement, put it in the readable task file.
4. After the main pipeline, visually review regenerated scatter and crop figures plus both summaries; this handoff does not claim that regenerated outputs have already been checked.

## Latest summary addendum completed

Applied the new requirements in confirmed_lumen_summary_addendum.md. Both summaries now include the same concise three-item scientific-limitations passage: light/shadow identification and possible liquid/reflection mislabeling; circular-area reconstruction with side-view/depth analogy and refraction; and explicitly superseded D60 crop history. Earlier repeated explanations were shortened while retaining nearby limitations and all seven method steps.

Historical values are loaded by write_summary.py directly from the frozen D60 adversarial JSON, not hardcoded or recomputed with the current D53 detector. Verified source values for m06q80_001: 0.6825585171295058 to 0.8575148822811562, rounded in text to 0.6826–0.8575; the old 53/60 squared maximum is 0.7802777777777777, rounded to 0.7803 with an approximation symbol. PT-BR formatting uses decimal commas. The writer asserts that these historical crop intervals are exactly 20:70, 18:71, 16:73 and 15:75 and that their analytic maxima use D60.

Both passages link the frozen historical experiment and the current normalization review. They explicitly say that the wider old crops exceeded the confirmed lumen and that the historical range is not current uncertainty. Current internal-row-removal sensitivity remains next to its current chart. Syntax checks passed after the addendum; no pipeline stage was run.
