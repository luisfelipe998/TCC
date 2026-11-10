# Independent centroid review

**PASS:** all five saved centroids match independent arithmetic. Rendered SVG marker coordinates and the PNG were reviewed after regeneration.

The intended centroid is the arithmetic mean of the displayed acquisition coordinates, with equal weight per acquisition: c_p = (sum(alpha_r)/n_p, sum(1000 F_r)/n_p). It is neither a median nor the centroid of the shaded convex polygon. Panorama length does not provide additional weight.

Expected coordinates from the saved per-run void and flip tables, calculated with math.fsum; selected by supplied pattern and saved core saturation <=0.20, independently of plotting code:

| Pattern | n | Alpha mean | Flips / 1,000 pixels mean |
|---|---:|---:|---:|
| dispersed bubbles | 1 | 0.001768204847710411 | 2.464139039342896 |
| plug | 3 | 0.5680338071082838 | 9.167131841634912 |
| slug | 3 | 0.6108514892109951 | 11.04687107352387 |
| churn | 10 | 0.8247855970067979 | 9.766872898630911 |
| wavy annular | 8 | 0.9715669024716174 | 2.237754367976236 |

All 25 eligible acquisitions are included once. The saturated m03q40_001 is excluded. The singleton dispersed-bubble centroid must exactly coincide with m07q01_001; it must not be offset to make both markers separately visible.

Input hashes:
- results/void_fraction/panorama_void_per_run.csv: 4b5bc79b37df34d2c2e1c1ea501b2d086c5bedd753b018b20a32ea593ec7899c
- results/flip_rate/phase_mask_flip_per_run.csv: 9e7d88ebbb47925a65f03fa22c2bb9f232a427d127717109b138b3f24192db8a

## Final verification

- Metadata means: maximum absolute difference from independent math.fsum means = 1.11e-16.
- All five named SVG centroid groups contain one plotted marker at the expected transformed mean; maximum SVG-coordinate difference = 4.53e-07 points (SVG serialization tolerance).
- Original PNG opened: five stars are clearly visible, each uses its pattern color; the legend identifies Centroide (média). No clipping or misleading coordinate offset was observed. The singleton bubble star coincides with its original point.
- The included acquisitions match the 25 displayed points exactly. No acquisition excluded by the saturation criterion or by the user contributes.
- Both summaries distinguish centroid means from the median values in their tables and specify equal weight per acquisition.
- No production files were edited and no pipeline was run by this reviewer.
- Reviewed primary PNG SHA-256: 7b6969b99377bd3188fcab20ca68db424b9c1fd0708fa2f20620d59d43b3339d.
