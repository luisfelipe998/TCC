"""Single dataset selection policy for every current v14 analysis stage."""
from collections import Counter
from pathlib import Path

from v14_paths import DATA

PATTERN_ORDER = ('dispersed bubbles', 'plug', 'slug', 'churn', 'wavy annular')
PATTERNS = PATTERN_ORDER
EXCLUDED_RUNS = frozenset({'m05q05_002'})
SATURATION_THRESHOLD = 0.20


def folder_pattern(folder):
    folder = Path(folder)
    markers = list(folder.glob('__*.txt'))
    if len(markers) != 1:
        raise ValueError(f'Expected one authoritative pattern marker in {folder}')
    return markers[0].stem[2:].replace('_', ' ')


def all_folders():
    return [folder for folder in sorted(DATA.glob('m??q??_???')) if folder.is_dir()]


def exclusion_reason(run, pattern):
    run = Path(run).name
    pattern = pattern.replace('_', ' ').strip()
    if run in EXCLUDED_RUNS:
        return 'Acquisition excluded by user: m05q05_002'
    if pattern == 'smooth annular':
        return 'Smooth annular excluded; only wavy annular retained'
    if '-' in pattern:
        return 'Compound pattern label excluded by user'
    if pattern not in PATTERN_ORDER:
        return 'Pattern outside the five approved single-pattern classes'
    return None


def is_in_scope(run, pattern):
    return exclusion_reason(run, pattern) is None


def is_eligible(run, pattern, saturation):
    return is_in_scope(run, pattern) and float(saturation) <= SATURATION_THRESHOLD


def selected_folders():
    return [folder for folder in all_folders()
            if is_in_scope(folder.name, folder_pattern(folder))]


def iter_runs():
    """Yield selected folder Paths; compatibility name for analysis loops."""
    return iter(selected_folders())


def selected_runs():
    return tuple(folder.name for folder in selected_folders())


def scope_metadata():
    selected, excluded = [], []
    for folder in all_folders():
        pattern = folder_pattern(folder)
        reason = exclusion_reason(folder.name, pattern)
        item = {'run': folder.name, 'supplied_pattern': pattern}
        if reason is None:
            selected.append(item)
        else:
            excluded.append({**item, 'reason': reason})
    return {
        'policy': 'Five single-pattern classes; wavy annular only; m05q05_002 excluded',
        'allowed_patterns': list(PATTERN_ORDER),
        'explicitly_excluded_acquisitions': sorted(EXCLUDED_RUNS),
        'source_acquisitions': len(selected) + len(excluded),
        'selected_acquisitions': len(selected),
        'selected_runs': [item['run'] for item in selected],
        'selected_pattern_counts': dict(Counter(item['supplied_pattern'] for item in selected)),
        'excluded_acquisitions': excluded,
        'statistical_rule': 'Selected acquisitions with core saturation fraction <= 0.20',
        'scope_does_not_relabel_inputs': True,
    }
