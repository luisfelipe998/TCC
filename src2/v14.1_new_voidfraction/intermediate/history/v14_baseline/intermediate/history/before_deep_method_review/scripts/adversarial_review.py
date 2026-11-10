"""Independent, panorama-only phase anchors and crop/reference perturbations.

Visual anchors validate phase labeling in selected interiors, not true gas volume.
No gate or reconstruction parameter is fitted to these checks.
"""
from v14_paths import ROOT as HERE, DATA, CACHE, artifact_path, write_markdown
from pathlib import Path
import hashlib
import json
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter1d, binary_closing, label
import panorama_void_alternatives as estimator
from v14_selection import selected_runs, scope_metadata

ANCHORS = [
    ("m07q01_001", 180, 250, 24, 65, "liquid, held-out run"),
    ("m07q01_001", 580, 650, 24, 65, "liquid, held-out run"),
    ("m06q01_001", 455, 515, 24, 65, "liquid plug gap"),
    ("m06q01_001", 150, 300, 30, 55, "gas body"),
]
SELECTED = ["m07q01_001", "m06q01_001",
            "m06q05_001", "m06q20_001", "m06q80_001", "m03q20_001"]


def main():
    assert set(SELECTED) <= set(selected_runs())
    assert all(a[0] in selected_runs() for a in ANCHORS)
    reference, _ = estimator.make_reference()
    checks = []
    for run, x0, x1, y0, y1, visual in ANCHORS:
        im = estimator.read(run)
        mask, accepted, _ = estimator.silhouette(im, reference)
        values = []
        for gate in (.60, .65, .70, .75, .80):
            alternative, _, _ = estimator.silhouette(im, reference, threshold=gate)
            values.append(float(alternative[y0:y1, x0:x1].mean()))
        checks.append(dict(run=run, x_interval=[x0, x1], y_interval=[y0, y1],
                           visual_anchor=visual, gas_mask_fraction=float(mask[y0:y1, x0:x1].mean()),
                           gas_fraction_at_gates_060_065_070_075_080=values,
                           accepted_center_fraction=float(accepted[x0:x1].mean())))
    refs = {"automatic_three_plug_reference": reference}
    for run, x0, x1 in [("m06q01_001", 455, 515), ("m07q01_001", 180, 250)]:
        refs[f"manual_liquid_{run}_{x0}_{x1}"] = gaussian_filter1d(
            np.median(estimator.normalize(estimator.read(run))[:, x0:x1], axis=1), 1)
    manual, crops = [], []
    for run in SELECTED:
        im = estimator.read(run)
        rvalues = {}
        for name, ref in refs.items():
            mask, _, _ = estimator.silhouette(im, ref)
            rvalues[name] = float(estimator.geometry(mask).mean())
        manual.append(dict(run=run, conditional_alpha_by_reference=rvalues))
        cvalues = {}
        for bounds in [(15, 75), (16, 73), (18, 71), (20, 70)]:
            mask, _, _ = estimator.silhouette(im, reference, rows=bounds)
            alpha = estimator.geometry(mask)
            cvalues[str(bounds)] = dict(conditional_alpha=float(alpha.mean()),
                                        maximum_column_alpha=float(alpha.max()),
                                        analytic_crop_ceiling=(bounds[1]-bounds[0])**2 / 60**2)
        crops.append(dict(run=run, crop_variants=cvalues))
    im = estimator.read("m06q05_001")
    mask, accepted, ratio = estimator.silhouette(im, reference)
    cols = np.arange(im.shape[1])
    top = 19 + np.argmin(ratio[19:38], axis=0)
    bottom = 42 + np.argmin(ratio[42:71], axis=0)
    components, _ = label(binary_closing(ratio < .70, structure=np.ones((3, 3))))
    same = (components[top, cols] == components[bottom, cols]) & (components[top, cols] > 0)
    region = np.zeros(im.shape[1], bool); region[300:560] = True
    disconnected = np.flatnonzero(accepted & ~same & region).tolist()
    x0, x1, scale = 280, 600, 3
    panel = Image.new("RGB", (2*(x1-x0)*scale, 80*scale), "white")
    panel.paste(Image.fromarray(im[:, x0:x1].astype("uint8")).convert("RGB")
                .resize(((x1-x0)*scale, 80*scale)), (0, 0))
    rgb = np.stack([im, im, im], axis=-1)
    rgb[mask] = .5*rgb[mask] + .5*np.array([0, 180, 255])
    panel.paste(Image.fromarray(rgb[:, x0:x1].astype("uint8")).resize(((x1-x0)*scale, 80*scale)),
                ((x1-x0)*scale, 0))
    panel.save(artifact_path("adversarial_slug_gap_mask.png"))
    output = dict(estimator_sha256=hashlib.sha256(Path(estimator.__file__).read_bytes()).hexdigest(),
                  basis="Retained src2 panoramas only; existing visual anchors selected independently of threshold tuning",
                  dataset_scope=scope_metadata(),
                  anchors=checks, manual_liquid_reference_sensitivity=manual, crop_sensitivity=crops,
                  disconnected_paired_acceptances_in_slug_gap_x300_560=disconnected,
                  limits="Selected phase interiors are not an unbiased segmentation test set or volumetric calibration; perturbation intervals are not confidence intervals")
    (artifact_path("adversarial_checks.json")).write_text(json.dumps(output, indent=2)+"\n")
    write_current_review(output)
    print(json.dumps(dict(anchor_gas_mask_fractions=[x["gas_mask_fraction"] for x in checks],
                          disconnected_slug_gap_pairs=disconnected, estimator_sha256=output["estimator_sha256"])))


def write_current_review(output):
    """Regenerate current scope evidence, keeping the earlier visual review archived."""
    records = output['anchors']
    lines = ['# Adversarial phase, reference and crop review', '',
             f"Current scope: {len(selected_runs())} acquisitions in the five approved single patterns. Only retained acquisitions contribute to the checks below. The [selection manifest](dataset_selection.json) documents exclusions.", '',
             f"Current segmentation source SHA-256: {output['estimator_sha256']}.", '',
             '## Selected visual interior checks', '',
             'These existing gas/liquid interpretations come from the previous visual review; excluded-run patches are removed. Their current mask fractions are recomputed below.', '',
             '| Run | x interval | y interval | Visual interpretation | Gas-map fraction |',
             '|---|---|---|---|---:|']
    for r in records:
        lines.append(f"| {r['run']} | {r['x_interval']} | {r['y_interval']} | {r['visual_anchor']} | {r['gas_mask_fraction']:.6f} |")
    lines += ['', '**Limit:** selected clear patches are not an unbiased boundary-validation set. The retained gas-positive patch is from a plug reference panorama; independent positive-control coverage is reduced by the requested exclusion.', '',
              '![Retained visual interior patches](adversarial_phase_anchor_patches.png)', '',
              '## Liquid-reference alternatives', '',
              'Automatic reference construction still uses the same three retained plug panoramas. Manual alternatives now use only retained liquid patches.', '',
              '| Run | Conditional alpha range across retained references |', '|---|---:|']
    for r in output['manual_liquid_reference_sensitivity']:
        values = list(r['conditional_alpha_by_reference'].values())
        lines.append(f"| {r['run']} | {min(values):.6f}–{max(values):.6f} |")
    lines += ['', '## Crop sensitivity', '',
              '| Run | Rows 20:70 | Rows 18:71, default | Rows 16:73 | Rows 15:75 |', '|---|---:|---:|---:|---:|']
    for r in output['crop_sensitivity']:
        values = [r['crop_variants'][str(bounds)]['conditional_alpha'] for bounds in [(20, 70), (18, 71), (16, 73), (15, 75)]]
        lines.append('| ' + r['run'] + ' | ' + ' | '.join(f'{v:.6f}' for v in values) + ' |')
    lines += ['', '**Limit:** the default 53-row crop imposes the conditional maximum (53/60)^2 = 0.7803. Wider crops can include wall/reflection shadows. These ranges describe sensitivity, not physical confidence intervals.', '',
              '## Corrected false-bridge case', '',
              f"The liquid slug in m06q05_001, x300:560, has {len(output['disconnected_paired_acceptances_in_slug_gap_x300_560'])} accepted paired fills between different closed-outline components. The original false bridge and its correction remain documented in the [historical review](intermediate/history/pre_scope_exclusions/findings/adversarial_review.md).", '',
              '![Current corrected liquid-gap overlay](adversarial_slug_gap_mask.png)', '',
              'The connection/interval rule can still miss fragmented true gas outlines or retain a long false bridge. Both void and flip-rate come from the same inferred map; their agreement is not independent physical validation.', '',
              'Use the [current machine checks](adversarial_checks.json), [final arithmetic verification](final_verification.json) and [post-regeneration adversarial review](post_regeneration_adversarial_review.md). The previous wider-scope evidence is historical, not a current result.']
    write_markdown(artifact_path('adversarial_review.md'), '\n'.join(lines) + '\n')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(12, 6))
    for ax, r in zip(axes.flat, records):
        x0, x1 = r['x_interval']; y0, y1 = r['y_interval']
        patch = estimator.read(r['run'])[y0:y1, x0:x1]
        ax.imshow(patch, cmap='gray', vmin=0, vmax=255, aspect='auto')
        ax.set_title(f"{r['run']} — {r['visual_anchor']}\nx{x0}:{x1}, y{y0}:{y1}; inferred gas {r['gas_mask_fraction']:.3f}", fontsize=10)
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle('Retained visual interior checks: selected examples, not volumetric calibration')
    fig.tight_layout()
    fig.savefig(artifact_path('adversarial_phase_anchor_patches.png'), dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
