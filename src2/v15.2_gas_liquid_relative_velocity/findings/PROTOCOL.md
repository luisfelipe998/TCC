# v15.2 protocol: gas/liquid relative velocity

Protocol recorded before the new numerical results. The study uses the 26 eligible recordings from the v15.1 inventory, all 501 native frames per recording. Mixed classifications, smooth annular and m05q05_002 remain excluded. Native camera FPS is used. Original images and earlier analyses are read-only inputs.

## Questions and estimands

1. Can persistent, isolated closed-outline objects provide gas-bubble translation estimates that distinguish translation from changing length? Can smaller objects outside large gas candidates provide independent natural-tracer evidence for liquid motion?
2. Does removing an image contribution predicted from the instantaneous interface geometry leave reproducible, phase-interior optical transport? Removing a fixed temporal background is a separate comparator.
3. What additional assumptions would be required to turn those measurements into liquid/gas material velocities, or their ratio? Test non-identifiability and physical tracer limitations explicitly.

The primary units are px/s, signed positive toward increasing image columns. Within-recording ratios cancel a common spatial scale. Cross-recording px/s comparisons concern image motion and require common magnification to become physical speed comparisons. A ratio does not repair phase identification.

## Object tracking

Use the v14 shadow detector with open-envelope acceptance disabled and thresholds 0.60, 0.70 and 0.80. A closed outline remains a candidate: it can include shadows or unresolved spots. Keep fully visible components inside rows 20:69 and columns 10:1006. Small candidates: area 15–1000 px, width 3–60 px, height 3–35 px, aspect ratio at most 4. Large candidates: area >1000 px, width 60–900 px, height at least 20 px. No per-pattern tuning.

Match reciprocal nearest centroids with displacement at most 60 px/interval, radial displacement at most 5 px, area ratio 0.67–1.5, and a 3 px advantage over a competing eligible match. Keep trajectories with at least five intervals, 95th-percentile changes in interval displacement no greater than 3 px, linear centroid-fit residual no greater than 2 px, and consistent front/rear translation: their fitted velocities must differ by at most max(2 px/frame, 25% of centroid speed). Also examine the stricter 10% deformation criterion. These gates select persistent, sufficiently rigid objects; they do not characterize all gas.

For small objects, measure the distance to large candidates in the same frame; the candidate's own outline is not counted as the host interface. An object at least 5 px from large candidates and tube boundaries is called an isolated small object, not a liquid material tracer.

Summaries weight each trajectory equally (median fitted centroid velocity). Also retain interval-weighted values to expose weighting sensitivity. Recording-level support requires at least three accepted tracks, at least 50 distinct intervals and observations in at least three of five time blocks. Threshold robustness requires support at all three thresholds and differences within max(1.5 px/frame, 20% of reference speed), with a consistent sign. These thresholds are declared screening choices, not confidence guarantees.

Synthetic controls cover known translation, growth, competing objects and missing texture. Time reversal and a longer matching lag check direction and correspondence. Inspect actual overlays and representative accepted/rejected tracks. Do not adjust gates to fill missing runs.

## Interface-conditioned residuals

Compare raw intensities, a cross-half temporal background, and an interface-geometry predictor estimated on the other temporal half. Use source AND displaced destination phase-candidate masks, erosion sensitivity, several lags, block consistency and shuffled controls. The predictor must not contain an imposed transport velocity. Its own script records exact gates before execution. Synthetic interface-only and shared-illumination controls distinguish numerical performance from physical identifiability.

## Additional analyses specified during implementation

Alongside the projected centroid, compute an axisymmetry-dependent centroid using the sum of squared contiguous radial gas spans as axial weights. This is a geometric sensitivity, not a confirmed reconstruction of three-dimensional volume. The median area fraction over the central half of an object's length is recorded for possible conditional conservation bounds. These additions were specified during implementation, before inspecting accepted trajectories.

The physical review supplies a separate, assumption-dependent branch for isolated elongated bubbles: if a traveling cell has an all-liquid piston, negligible phase change, constant densities, and a film mean velocity between zero and bubble speed, conservation gives alpha_body <= U_liquid_piston / c_bubble <= 1. This is not the mean liquid velocity of the whole cell or video. The empirical geometry and adjacent mask-clear sections will be screened before assigning any conditional bounds. They do not verify the film hypothesis or the absence of gas outside the viewing plane. The exact implementation records additional eligibility checks; it is an exploratory branch, not a validated phase measurement.

The first object results showed that wall-connected shadows defeated the fully closed large-object detector, including visually elongated bubbles. A separate exploratory meniscus method was therefore specified: contiguous columns with more than 50% candidate-gas occupancy in rows 35:54, width 60–900 px, both ends inside columns 10:1006. Track the midpoint and the two ends with the same correspondence/deformation gates. Retain and report radial boundary contact instead of rejecting it. This changes the observable to axial gas-envelope/interface translation. It does not establish a correct gas diameter or validate the liquid-film model. Run all three segmentation thresholds on all 26 recordings, including failures. Do not merge its output with the small-object velocity as if they were the same observable.

That endpoint experiment exposed localization jumps despite similar mean endpoint slopes. A further exploratory estimator fits Theil–Sen slopes to each end of the SAME matched chains. Its gates are: >=20 intervals, speed magnitude >2 px/frame, net fitted displacement >=50 px, endpoint median absolute residual <=3 px and 90th percentile <=12 px, and agreement of both endpoint speeds and separate temporal-half speeds within max(1.5 px/frame,20%). The result is the mean of the two endpoint slopes. Known-translation, localization-jitter, growth, stationary and accelerated synthetic sequences are checked before real processing. No claim is made that synthetic jitter reproduces the experimental errors. This is a disclosed methodological development on the available recordings, not independent external validation.

The audit also found stationary small candidates at nearly fixed image columns. Retain original track tables, then add a motion-resolved sensitivity requiring |v|>2 px/frame and fitted travel >=20 px. Compare it against the original median; do not rename the excluded stationary candidates as slow liquid. Ratios between small-object and large-envelope velocities require supported, threshold-stable estimates of both, and remain ratios of these different optical objects.

Two further exploratory comparisons are documented separately. Common-interval ratios average the fitted velocities of all eligible objects present at each shared interval, then average over those intervals. At least 50 common intervals, three time blocks and three trajectories from each population are required, with threshold support and speed agreement. Fits themselves use whole trajectories, so the restriction concerns aggregation rather than exclusive fitting on shared intervals. Between-recording comparisons use only the same small-object definition in five threshold-stable recordings.

The zero-size probe compares a constant with linear height, squared height, and squared height plus squared radial position models. It uses at least 12 motion-resolved small trajectories and three temporal groups. Held-out groups are assigned by trajectory midpoint; training trajectories sharing any interval with a test trajectory are purged. All models and their zero-size intercepts are retained, and none is promoted to liquid speed. An early diagnostic used unpurged groups; the final output and independent verifier use the purged version. Separate tracks can still be statistically dependent, and this internal predictive check does not validate a material tracer or the extrapolation to zero size.

## Reporting

Report every eligible recording, including failures, with the same definitions. No zeros for unavailable velocities; no pooling of incompatible observables. Temporal blocks and trajectories are not independent experimental replications. Within-pattern spread and cross-pattern comparisons are descriptive of these recordings, not causal effects of pattern. A single dispersed recording does not establish within-pattern variation between recordings.

The success criterion for liquid/gas velocities requires evidence that both observables follow the respective material phases. If no method establishes that, state failure explicitly while preserving useful gas-object or phase-candidate optical results. Numerical/synthetic validation is not experimental velocity ground truth. Publish English and PT-BR executive summaries, PT-BR figures, numerical tables, scripts, limitations and provenance in this directory.
