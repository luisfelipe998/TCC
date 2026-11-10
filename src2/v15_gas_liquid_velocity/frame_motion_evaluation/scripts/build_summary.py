#!/usr/bin/env python3
"""Build reviewable frame-motion findings from saved numerical evidence."""
from pathlib import Path
import csv
import json

BASE=Path(__file__).resolve().parents[1]


def read(name):
    return list(csv.DictReader((BASE/'results'/name).open()))


def pct(value):
    return '—' if value=='' else f'{100*float(value):.1f}'


def run():
    rows=read('temporal_per_run.csv')
    local=read('local_patch_summary.csv')
    sensitivity={r['run']:r for r in json.loads((BASE/'results'/'independent_existing_sensitivity.json').read_text())}
    meta=json.loads((BASE/'results'/'temporal_metadata.json').read_text())
    stable=[r for r in rows if r['stable_under_lag_diagnostics']=='True']
    totals=meta['totals']
    assert len(rows)==26 and len(stable)==7 and totals['base_supported']==17
    bg=[r for r in local if r['method']=='background_subtracted']
    annular=[r for r in bg if r['pattern']=='wavy annular']
    assert sum(int(r['accepted_patch_tracks']) for r in bg)==3880
    assert sum(int(r['accepted_patch_tracks']) for r in annular)==434
    assert sum(int(r['coherent_pairs']) for r in annular)==10
    lines=[
        '# EXECUTIVE SUMMARY — Motion between native frames',
        '',
        '**Yes: the 501 native frames support estimates of moving image structures in px/frame and px/s. They do not establish two independent liquid and gas material velocities.** This evaluation uses the existing images and native camera rates, without inlet-state assumptions, mass-flow matching, or a physical pixel scale.',
        '',
        'The original scope is preserved: **26 recordings**, excluding m05q05_002, smooth annular and double classifications. All retained recordings contain consecutive frame identifiers 0–500. The two camera configuration formats agree: 18 recordings at 7,000 fps, six at 5,000 fps and two at 4,000 fps. Each recording therefore spans **0.07143–0.125 s**, with 500 adjacent comparisons. Those overlapping comparisons are not 500 independent experimental repeats.',
        '',
        '## Results of the feasibility tests',
        '',
        '| Test | Actual result | Interpretation |',
        '|---|---|---|',
        f'| Whole-profile adjacent-frame correlation | {totals["adjacent_pass"]:,}/13,000 accepted pairs; {totals["base_supported"]}/26 recordings pass coverage and saturation gates | Many image structures have measurable displacement; this is initial numerical support. |',
        f'| Compare frame gaps of 1, 2 and 4 | {totals["temporal_pairs"]:,} comparisons; {len(stable)}/26 recordings pass the declared temporal diagnostics | Stronger evidence for persistent structure motion; two of these seven remain materially sensitive to region/selection. |',
        '| Local 2-D patches with reciprocal matching | 3,880/7,800 background-subtracted tracks accepted | Reveals local motion and spatial variability; individual matches do not imply one common velocity across the image. |',
        '| Annular local patches | 434/2,700 accepted tracks, but only 10/450 sampled pairs with distributed agreement; no annular recording meets the reporting gate | Isolated annular events can be matched; a representative annular velocity remains unsupported by these tests. |',
        '| Visual multi-frame review | Six recordings, covering all five retained patterns | Clear bubble/plug/slug examples translate left; churn and annular examples show stronger deformation or correspondence ambiguity. |',
        '',
        '**The 17-recording initial-support result must not be described as 17 validated velocities.** Longer frame gaps reduce correspondence support, especially in churn. The local multi-window gate asks a stricter, different question; its low reporting count does not mean the other recordings contain no trackable features.',
        '',
        '![Motion and comparison coverage](../figures/temporal_motion_overview.png)',
        '',
        'Plotted values are signed accepted-pair medians. Horizontal bars are the accepted-pair interquartile range, not measurement confidence intervals. Positive x is rightward, so these negative values describe motion toward smaller image coordinates. Unsupported annular estimates are withheld from the velocity panel. The right panel shows joint comparison coverage, including its denominator; agreement based on very few accepted events must not appear to be strong evidence.',
        '',
        '## Recordings passing the temporal diagnostic',
        '',
        'These are structure-motion candidates. The table retains all seven temporal passes, including the two with large selection sensitivity. Numeric rounding is for readability and does not assert instrument precision.',
        '',
        '| Recording | Pattern | Accepted adjacent pairs | Median px/frame | Median px/s | Accepted-variant median span (%) |',
        '|---|---|---:|---:|---:|---:|',
    ]
    for r in stable:
        s=sensitivity[r['run']]
        lines.append(f'| {r["run"]} | {r["pattern"]} | {r["adjacent_pass_pairs"]}/500 | {float(r["median_signed_px_frame"]):.3f} | {float(r["median_signed_px_s"]):,.0f} | {s["supported_range_pct_baseline"]:.2f} |')
    lines += [
        '',
        'The final column is the maximum minus minimum median among previously tested ROI/lag variants with at least 50% accepted pairs, divided by the magnitude of the baseline median. The variants accept different frame subsets. This is sensitivity to analysis/selection, not velocity error or an uncertainty bound.',
        '',
        '**Recommended first reporting candidates:** m07q01_001, m05q01_001, m05q05_001, m06q01_001 and m07q05_001. Their supported-variant median spans are approximately 2.1–4.5%. This recommendation is based on the observed diagnostics, not a new accuracy-validation gate. m05q01 has eight of nine ROI/lag variants above the 50% coverage gate; the remaining variant has 49%, so it is not uniformly supported under every setting.',
        '',
        '**Keep m04q01_001 and m06q05_001 explicitly qualified:** their supported-variant spans are 36.37% and 19.41%. Passing pairwise temporal checks does not make a selected per-recording median invariant to which interfaces and frame pairs are accepted.',
        '',
        '## Concrete evidence of frame displacement',
        '',
        '![Native plug sequence](../figures/visual_sequence_m05q01_001.png)',
        '',
        'The plug sequence shows recognizable interfaces at frames 200, 202, 204, 206, 208 and 210. The displayed 10-frame interval lasts 1.429 ms at 7,000 fps. Across the complete recording, the accepted adjacent-pair median is −11.886 px/frame, giving −83,199 px/s. The short contact sheet illustrates correspondence; the reported median uses all accepted adjacent comparisons across the 501-frame recording.',
        '',
        'For any matched feature, image velocity is displacement in pixels multiplied by capture rate and divided by the number of frame intervals. A 2-frame displacement must be divided by two. Capture rate is taken from each native configuration, not from the resampled 4,000-Hz workbook signals or camera playback rate.',
        '',
        '## What can be claimed',
        '',
        '- Report **observed axial structure/interface motion in px/s** by scenario, with accepted coverage, temporal checks, region sensitivity and spreads.',
        '- Describe identifiable objects as bubble translation or interface/wave propagation when their correspondence is supported. A global brightness-profile shift is more generally an aggregate structure-motion estimate.',
        '- Preserve rejected pairs and unsupported scenarios. No representative annular speed is reported from the present tested methods.',
        '- Keep cross-scenario comparisons in image coordinates. Physical speed rankings require consistent magnification; equal pixel dimensions do not prove an equal length scale.',
        '- Treat the observed motion as a short-record summary. The available sequence durations do not establish long-term operating variability.',
        '',
        '## What remains unidentified',
        '',
        '**Separate liquid and gas material velocities remain unidentified.** The visible boundaries belong to both phase masks. Tracking a liquid plug boundary measures an interface, not a liquid parcel. The visual review did not establish liquid tracers in its six examples. A successful patch match can follow a bubble, reflection, refractive texture or wave, and does not supply phase mass/volume weighting. These quantities also differ from superficial gas/liquid velocities.',
        '',
        'Changing to px/s removes the need for a verified mm/pixel scale for the reported observable. It does not establish phase identity, physical speed, or measurement accuracy. Finite exposure, changing shape, saturation and spatially varying motion remain relevant. The apparent displacement during the exposure is stored as a diagnostic estimate, not a measured blur length.',
        '',
        '## Verification and deliverables',
        '',
        'Independent reaggregation checked **307,619 numeric values** across the baseline, temporal and local-patch outputs. Scope, rates, all pixel-rate conversions, summary statistics and declared decision rules reproduce. Synthetic translated textures verify sign and basic numerical behavior. These checks establish computational consistency; there is no independent experimental velocity reference.',
        '',
        'Details: [temporal methods and all 26 results](temporal_review.md), [local 2-D checks](local_patch_review.md), [native visual review](visual_review.md), [independent review](independent_review.md), and [reproduction](../README.md). All scripts, CSVs and figures are stored within this evaluation. Earlier v15 results and the source data are preserved.',
    ]
    (BASE/'findings'/'EXECUTIVE_SUMMARY.md').write_text('\n'.join(lines)+'\n')
    detail=[
        '# Temporal consistency of native-frame motion',
        '',
        'The preferred v15 normalized-correlation implementation was reapplied directly to all 501 native images per retained recording. Its original 13,000 adjacent comparisons reproduce exactly, including all 7,528 acceptance decisions and per-run medians. Images remain at 80×1016 pixels. Profiles average rows 20:60; each axial pixel has its 501-frame temporal mean subtracted. No high-pass filter, panorama, forced-positive displacement or physical length conversion is used.',
        '',
        '## Fixed correspondence and diagnostic rules',
        '',
        'A fixed 536-pixel central template is compared against equal-length candidate windows over ±240 pixels. An accepted comparison requires normalized correlation ≥0.60, margin ≥0.03 over the strongest candidate outside ±5 pixels, profile standard deviation ≥1 gray level, and a peak more than 2 pixels inside the search limit. Native central-band saturation at grayscale ≥250 must occupy no more than 20% of the recording for a run-level reporting candidate.',
        '',
        'The same profile method is applied at lags 1, 2 and 4 to every possible pair: 500+499+497=1,496 comparisons per recording and 38,896 total. For a longer gap, its displacement divided by the lag is compared with the mean of the intervening adjacent displacements. Both the longer comparison and every adjacent comparison on that path must pass to enter the joint diagnostic.',
        '',
        'A joint comparison agrees if the difference is ≤max(1 px/frame, 15% of the absolute adjacent-path mean). A run passes the temporal diagnostic if its adjacent coverage is ≥50%, saturation ≤20%, each of lags 2 and 4 has joint coverage ≥25%, and at least 80% of each lag’s joint comparisons agree. These thresholds were fixed before this computation. They are transparent heuristic diagnostics, not accuracy levels.',
        '',
        'Long-gap and adjacent templates cover changing spatial content. Their disagreement can represent genuine spatial heterogeneity, deformation or different structures as well as correspondence failure. Fixed-window reversed comparisons are also saved; because their central templates are spatially distinct, these are reversal-sensitivity diagnostics, not reciprocal tracking of an identical patch. The separate local 2-D experiment performs matched-location reciprocal checks.',
        '',
        '## Results for every retained recording',
        '',
        '| Recording | Pattern | Adjacent accepted (%) | Lag-2 joint coverage (%) | Lag-4 joint coverage (%) | Lag-2 agreement (%) | Lag-4 agreement (%) | Temporal diagnostic |',
        '|---|---|---:|---:|---:|---:|---:|---|',
    ]
    for r in rows:
        fields=[pct(r[k]) for k in ('adjacent_pass_fraction','lag2_joint_fraction','lag4_joint_fraction','lag2_within_tolerance_fraction','lag4_within_tolerance_fraction')]
        note='pass' if r['stable_under_lag_diagnostics']=='True' else 'not passed'
        if r['high_saturation']=='True': note+='; saturated'
        detail.append('| '+r['run']+' | '+r['pattern']+' | '+' | '.join(fields)+' | '+note+' |')
    detail += [
        '',
        'A dash means no jointly accepted comparisons. Agreement percentages use only the jointly accepted subset; coverage uses all possible comparisons for that lag. The 15% tolerance and gate percentages must not be described as confidence or measured error.',
        '',
        '![Full-sequence temporal examples](../figures/temporal_examples.png)',
        '',
        'The space–time panels use all 501 mean profiles after background subtraction. Their slopes visualize structure motion; no independent slope estimator is fitted here. Gray time-series points/lines retain rejected estimates; blue points indicate accepted adjacent comparisons.',
        '',
        'Five deterministic synthetic motion sequences at three lags verify the expected shift sign and px/frame calculation, with maximum error below 0.1 px/frame on these particular translated textures. Constant images are rejected. The synthetic check does not model multiphase deformation or establish real-data error.',
        '',
        'Full results: [all pairs](../results/temporal_pairs.csv), [per-run summaries](../results/temporal_per_run.csv), [input hashes, thresholds and checks](../results/temporal_metadata.json), [independent verification](../results/independent_verification.json).',
    ]
    (BASE/'findings'/'temporal_review.md').write_text('\n'.join(detail)+'\n')
    readme=[
        '# Native-frame motion feasibility evaluation',
        '',
        'Start with the [EXECUTIVE SUMMARY](findings/EXECUTIVE_SUMMARY.md). This extends v15 using the existing 501-frame recordings to evaluate structure displacement in px/frame and px/s. It does not estimate thermodynamic quality or independently measured liquid/gas material velocities.',
        '',
        'The 26-run scope and original data are unchanged. Earlier v15 numerical artifacts are preserved. New computations, figures and findings reside in this directory.',
        '',
        '## Reproduce',
        '',
        'From the repository root, using the existing Python environment:',
        '',
        '    python3 -B src2/v15_gas_liquid_velocity/frame_motion_evaluation/scripts/temporal_motion.py',
        '    python3 -B src2/v15_gas_liquid_velocity/frame_motion_evaluation/scripts/local_patch_motion.py',
        '    python3 -B src2/v15_gas_liquid_velocity/frame_motion_evaluation/scripts/local_patch_report.py',
        '    python3 -B src2/v15_gas_liquid_velocity/frame_motion_evaluation/scripts/visual_sequences.py',
        '    python3 -B src2/v15_gas_liquid_velocity/frame_motion_evaluation/scripts/independent_verify.py',
        '    python3 -B src2/v15_gas_liquid_velocity/frame_motion_evaluation/scripts/build_summary.py',
        '',
        'The temporal script imports the existing v15 NCC functions; their hash is recorded. The independent verifier does not import either analysis implementation. The summary script reads completed numerical outputs and the independent ROI/lag audit. No new package installation or network access is needed.',
        '',
        '## Outputs',
        '',
        '- [Executive summary](findings/EXECUTIVE_SUMMARY.md): feasibility, actual results and recommendation.',
        '- [Temporal evaluation](findings/temporal_review.md): all 26 scenarios and fixed correspondence rules.',
        '- [Local 2-D evaluation](findings/local_patch_review.md): independent spatial correspondence test.',
        '- [Native visual sequences](findings/visual_review.md): six contact sheets and six full-duration space–time figures.',
        '- [Independent verification](findings/independent_review.md): arithmetic/provenance checks and interpretation limits.',
        '',
        'Signed positive motion is rightward in the stored image. Medians and interquartile ranges describe accepted image-motion events; they are not measurement uncertainty or phase-weighted velocities. Rejected events remain in the CSVs. The reports separate numerical support, temporal consistency, region sensitivity and physical identification.',
    ]
    (BASE/'README.md').write_text('\n'.join(readme)+'\n')
    print('Generated executive summary, temporal report and README from verified saved results.')


if __name__=='__main__':
    run()
