# Independent check of the original findings

This script reads the saved panoramas directly and implements Otsu independently 
using cumulative histogram arithmetic. All 31 reported spatial flip counts and 
denominators are reproduced; the maximum absolute rate difference is zero.

There are 31 acquisitions: 2 dispersed bubbles, 3 plug, 3 slug, 10 churn, 
9 wavy annular, 2 slug-churn, 1 plug-slug and 1 smooth annular. 
Keeping the five pure groups and removing >20% clipped core images leaves 
26 acquisitions, including both dispersed-bubble runs. Excluding m07q01 
then leaves 25; the original findings table sums to 25 although it says 26.

The three clipped records are m03q01, m03q40 and m03q80. m03q01 is also 
a transition label and m03q80 is smooth annular, so exclusions overlap. 
Do not subtract these categories twice.

The original file conflates two definitions. v11.1 counts Schmitt-trigger 
switches of a one-dimensional dark occupancy. v11.3 counts individual binary 
pixel changes across neighboring panorama columns, averaged over 40 rows. 
The v11.3 values near 0.02 are fractions of adjacent comparisons, 
not automatically switches per 1000 px or temporal frequencies.

Complementing the optical binary mask preserves its flip count exactly. 
Thus swapping dark and bright phase names cannot correct or validate flip-rate. 
Both complement and spatial-reversal invariance pass for all 31 images.

The one-view/scale argument also needs correction: an axisymmetric physical 
gas diameter can yield alpha=(d_g/D)^2 with absolute scale canceling. 
Unvalidated optical shadows and unknown cross-sectional shapes are the actual 
obstacles here; the 60-pixel analysis window is not an independent tube calibration.

The data inventory records SHA-256 for every input panorama to make the source 
images used in this check identifiable. The pipeline does not change them.

Reproduce: `python3 src2/v14_new_voidfraction/scripts/inventory_and_legacy_audit.py`.
