# Native phase-mask transfer from v14

**Completed:** the current v14 detector was applied to all 13,026 native frames in the 26 retained recordings, at shadow thresholds 0.60, 0.70 and 0.80. These 39,078 masks are gas-silhouette candidates. Unmarked pixels are not independently confirmed liquid.

The exact three-panorama liquid-like radial reference was recomputed and frozen. Native images, v14 scripts and previous outputs were not changed. The current owner-confirmed lumen is rows 18:71, or 53 pixel cells; no geometric void fraction or physical pixel scale is used in this transfer.

## Per-recording projected fractions and sensitivity

Fractions below average inferred gas pixels over all 501 native frames and all 53 lumen rows. They are projected optical fractions, not vapor quality or volumetric void. Disagreement is the fraction of lumen pixels whose labels differ among the three thresholds.

| Recording | Pattern | Threshold 0.60 (%) | Threshold 0.70 (%) | Threshold 0.80 (%) | Disagreement (%) |
|---|---|---:|---:|---:|---:|
| m03q05_001 | churn | 97.607 | 98.290 | 98.736 | 1.132 |
| m03q10_001 | wavy annular | 99.593 | 99.731 | 99.821 | 0.227 |
| m03q20_001 | wavy annular | 99.940 | 99.967 | 99.981 | 0.041 |
| m03q40_001 | wavy annular | 98.052 | 98.855 | 99.412 | 1.360 |
| m04q01_001 | slug | 80.950 | 83.315 | 87.948 | 7.699 |
| m04q10_001 | churn | 93.223 | 94.579 | 95.667 | 2.518 |
| m04q20_001 | churn | 93.953 | 95.292 | 96.335 | 2.407 |
| m04q40_001 | wavy annular | 98.842 | 99.184 | 99.392 | 0.551 |
| m04q80_001 | wavy annular | 99.929 | 99.961 | 99.979 | 0.050 |
| m05q01_001 | plug | 64.730 | 66.654 | 74.606 | 10.302 |
| m05q05_001 | slug | 58.584 | 61.321 | 69.705 | 11.621 |
| m05q10_001 | churn | 82.942 | 86.783 | 89.895 | 7.180 |
| m05q20_001 | churn | 85.505 | 89.459 | 92.312 | 6.963 |
| m05q40_001 | churn | 90.872 | 94.916 | 97.012 | 6.168 |
| m05q80_001 | wavy annular | 94.878 | 98.500 | 99.471 | 4.596 |
| m06q01_001 | plug | 64.039 | 70.152 | 75.319 | 11.932 |
| m06q05_001 | slug | 69.140 | 73.975 | 79.269 | 10.911 |
| m06q10_001 | churn | 77.307 | 81.628 | 85.920 | 9.150 |
| m06q20_001 | churn | 82.750 | 88.441 | 92.840 | 10.260 |
| m06q40_001 | wavy annular | 87.755 | 93.972 | 97.300 | 9.571 |
| m06q80_001 | wavy annular | 90.426 | 96.560 | 98.769 | 8.348 |
| m07q01_001 | dispersed bubbles | 0.455 | 15.184 | 23.679 | 24.369 |
| m07q05_001 | plug | 54.764 | 60.663 | 66.539 | 12.653 |
| m07q20_001 | churn | 73.688 | 78.791 | 83.486 | 10.519 |
| m07q40_001 | churn | 80.204 | 86.460 | 91.075 | 11.078 |
| m07q80_001 | wavy annular | 85.731 | 93.114 | 97.146 | 11.444 |

## Concrete transfer difficulties

**Dispersed-bubble m07q01_001 is a major sensitivity case:** the native projected gas fraction changes from 0.455% at threshold 0.60 to 15.184% at 0.70 and 23.679% at 0.80. The existing current-v14 panorama projected fraction is 0.821%. Native and panorama sampling differ, so this discrepancy does not identify which is correct; it demonstrates that panorama evidence cannot automatically validate native-frame labels.

In contrast, m05q01_001 and m03q05_001 have close native and saved-panorama projected fractions. The all-run comparison is saved with explicit temporal/spatial definitions. Do not compare projected fractions to circular-volume alpha, which is a different quantity.

The mask operation is not guaranteed monotonic in threshold: thresholding changes component connectivity, crop-contact rejection and hole filling. Observed nonmonotonic pixels reach approximately 1.061% of lumen pixels in a recording. Conservative consensus must intersect/union all three masks rather than assume one contains the others.

Gas-dominated churn and annular masks can leave very little detected liquid at a fixed probe. That is a phase-support limitation for a two-class conditional average. It does not establish that no liquid film or entrained liquid is physically present.

## Saved schema and verification

Each intermediate/masks/RUN.npz contains masks_060, masks_070 and masks_080 as uint8 arrays with shape (501,80,127), packed along the axial axis with big-endian bit order. It also stores mask_shape=(501,80,1016), original_width=1016, frame_ids, fps, thresholds, lumen_rows, bitorder and unpacked_dtype=bool.

Unpack with numpy.unpackbits(npz[key], axis=-1, count=int(npz["original_width"]), bitorder="big").astype(bool). Value 1 is inferred gas silhouette; value 0 is unmarked. Pixels outside rows 18:71 are always unmarked and must not become liquid labels.

Uniform images under a neutral reference and exact reference-shaped images produce no inferred gas at every threshold. Packed representations were compared with direct detector outputs, then reopened from disk and checked byte-for-byte. Fixed frames 0, 200 and 500 were also unpacked and compared exactly. Source, reference and packed-file hashes are retained. These tests establish faithful application and storage, not phase-label accuracy.

Reproduce with python3 -B scripts/native_phase_masks.py --workers 3 from the repository path for this script.

Evidence: [inventory](../results/native_phase_inventory.csv), [per-frame fractions](../results/native_phase_frame_fractions.csv), [reference](../results/native_phase_reference.json), [metadata and checks](../results/native_phase_metadata.json), [native/panorama comparison](../results/native_panorama_comparison.csv).
