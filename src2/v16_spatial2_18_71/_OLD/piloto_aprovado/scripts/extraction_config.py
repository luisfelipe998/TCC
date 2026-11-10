"""Faixas de extração declaradas para a análise 18–71 (índices a partir de zero)."""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
CONFIG=json.loads((ROOT/'_OLD/results/analysis_config.json').read_text())
ROW_START,ROW_STOP=CONFIG['numpy_row_slice']
ROW_END=ROW_STOP-1
ROW_COUNT=ROW_STOP-ROW_START
ROW_SLICE=slice(ROW_START,ROW_STOP)
REFERENCE_ROW_SLICE=slice(20,60)
FIXED_ROW=CONFIG['fixed_row_control']
FIXED_BAND_START,FIXED_BAND_END=CONFIG['fixed_band_control_inclusive']
FIXED_BAND_SLICE=slice(FIXED_BAND_START,FIXED_BAND_END+1)
UPPER_START,UPPER_END=CONFIG['upper_diagnostic_rows_inclusive']
LOWER_START,LOWER_END=CONFIG['lower_diagnostic_rows_inclusive']
UPPER_SLICE=slice(UPPER_START,UPPER_END+1)
LOWER_SLICE=slice(LOWER_START,LOWER_END+1)
TRIM_FRACTION=CONFIG['trimmed_mean_fraction_requested']
TRIM_COUNT=int(TRIM_FRACTION*ROW_COUNT)
TRIM_SLICE=slice(TRIM_COUNT,ROW_COUNT-TRIM_COUNT)

assert (ROW_START,ROW_STOP,ROW_COUNT)==(18,72,54)
assert CONFIG['extraction_rows_inclusive']==[ROW_START,ROW_END]
assert CONFIG['row_count']==ROW_COUNT
assert (FIXED_ROW,FIXED_BAND_START,FIXED_BAND_END)==(40,30,49)
assert UPPER_END+1==LOWER_START and UPPER_START==ROW_START and LOWER_END==ROW_END
assert TRIM_COUNT==CONFIG['trimmed_rows_each_end']==10
assert ROW_COUNT-2*TRIM_COUNT==CONFIG['retained_trimmed_rows']==34
