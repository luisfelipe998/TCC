"""Names distinguish patterns, profiles and contexts without altering pilot filenames."""
from common import slug


def pair_key(record):
    pattern_x = record.get("pattern_x", record["padrao"])
    pattern_y = record.get("pattern_y", pattern_x)
    pattern = slug(pattern_x) if pattern_x == pattern_y else f"{slug(pattern_x)}_vs_{slug(pattern_y)}"
    key = f"{pattern}_{record['run_x']}_{record['run_y']}"
    if record.get("perfil"):
        key += f"_{record['perfil']}_n{record['comprimento']}"
    return key
