"""Current single optical ROI, corrected on 2026-10-30.

54 occupied pixel cells: rows 18 through 71 inclusive, Python slice18:72.
Historical D53 definitions are preserved in the pre-unification snapshot.
"""
LUMEN_ROWS = (18, 72)
LUMEN_SLICE = slice(*LUMEN_ROWS)
LUMEN_HEIGHT_PX = LUMEN_ROWS[1] - LUMEN_ROWS[0]
DIAMETER_PX = float(LUMEN_HEIGHT_PX)
FIXED_PROBE_ROW = 40
INTERNAL_CROP_VARIANTS = (LUMEN_ROWS,)
BOUNDARY_PROVENANCE = "User correction2026-10-30: sole optical ROI18:72, rows18..71 inclusive, H=D=54 pixels."
