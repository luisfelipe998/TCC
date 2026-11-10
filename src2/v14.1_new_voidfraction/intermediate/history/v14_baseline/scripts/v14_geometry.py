"""Image geometry explicitly confirmed by the user on 2026-10-23.

Pixel widths count occupied rows (pixel cells), not distances between the
first and last pixel centers. The half-open slice 18:71 contains 53 cells.
These image boundaries do not calibrate refraction or gas cross-section shape.
"""

LUMEN_ROWS = (18, 71)
LUMEN_SLICE = slice(*LUMEN_ROWS)
LUMEN_HEIGHT_PX = LUMEN_ROWS[1] - LUMEN_ROWS[0]
DIAMETER_PX = float(LUMEN_HEIGHT_PX)
FIXED_PROBE_ROW = 40
INTERNAL_CROP_VARIANTS = ((20, 70), (19, 70), (19, 71), LUMEN_ROWS)
BOUNDARY_PROVENANCE = (
    "User-confirmed full tube interior on 2026-10-23: Python slice 18:71, "
    "rows 18 through 70, 53 pixel cells. Optical profile minima are not used "
    "as lumen boundaries. Refraction mapping remains uncalibrated."
)
