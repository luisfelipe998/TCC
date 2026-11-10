#!/usr/bin/env python3
"""Publish the saved phase-conditioned findings as the main v15 summary.

This only rebases report links and appends provenance navigation. It does not
convert units, recompute motion, or substitute the historical model report.
"""
from pathlib import Path
import os
import re
from urllib.parse import unquote, urlsplit, urlunsplit

BASE = Path(__file__).resolve().parents[1]
SOURCE = BASE / 'phase_conditioned_motion/findings/EXECUTIVE_SUMMARY.md'
DESTINATION = BASE / 'findings/EXECUTIVE_SUMMARY.md'


def rebase_links(markdown, source=SOURCE, destination=DESTINATION):
    """Preserve link meaning when moving a Markdown report between folders."""
    def replace(match):
        target = match.group(2)
        angle = target.startswith('<') and target.endswith('>')
        value = target[1:-1] if angle else target
        parts = urlsplit(value)
        if parts.scheme or parts.netloc or not parts.path:
            return match.group(0)
        path = (source.parent / unquote(parts.path)).resolve()
        relative = Path(os.path.relpath(path, destination.parent)).as_posix()
        rebased = urlunsplit(('', '', relative, parts.query, parts.fragment))
        if angle:
            rebased = '<' + rebased + '>'
        return match.group(1) + rebased + match.group(3)

    return re.sub(r'(\[[^\]\n]*\]\()([^\)\n]+)(\))', replace, markdown)


def current_summary_text():
    """Return the current report with its source evidence and units preserved."""
    markdown = SOURCE.read_text()
    if not markdown.startswith('# EXECUTIVE SUMMARY') or 'px/s' not in markdown:
        raise ValueError('Expected the saved pixel-based phase-conditioned summary')
    published = rebase_links(markdown).rstrip()
    appendix = [
        '', '', '## Earlier investigations and provenance', '',
        'This main report is synchronized from the [phase-conditioned executive summary](../phase_conditioned_motion/findings/EXECUTIVE_SUMMARY.md). Its velocity comparisons use px/s and dimensionless ratios. The current plots and numerical values are generated from the saved phase-conditioned results.', '',
        'The [historical conditional model summary](conditional_model_summary.md) preserves the earlier thermodynamic investigation, including its m/s predictions and inlet/property/closure assumptions. Those predictions are a different calculation; they are not physical-unit conversions of the image-motion results above.', '',
        'Supporting earlier work: [native-frame motion evaluation](../frame_motion_evaluation/findings/EXECUTIVE_SUMMARY.md), [conditional model tables](velocity_results.md), [model methods](method_comparison.md), [v12/v12.1 audit](legacy_audit.md), [recorded-flow audit](massflow_review.md), [literature review](literature_review.md), [independent model review](independent_review.md), and [reproduction](../README.md).',
    ]
    return published + '\n'.join(appendix) + '\n'


def publish():
    DESTINATION.write_text(current_summary_text())
    print(f'Published current px/s findings to {DESTINATION.relative_to(BASE)}')


if __name__ == '__main__':
    publish()
