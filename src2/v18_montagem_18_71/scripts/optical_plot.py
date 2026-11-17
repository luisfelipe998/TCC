"""Rendering helpers for the v18 optical scatter plot."""
from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
import numpy as np

from runtime import OUTPUT_ROOT

PATTERNS = ("dispersed bubbles", "plug", "slug", "churn", "wavy annular")


def artifact_path(name):
    path = OUTPUT_ROOT / "03_optica/figuras" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


COLORS = {
    "dispersed bubbles": "#e69f00", "plug": "#0072b2",
    "slug": "#009e73", "churn": "#d55e00", "wavy annular": "#cc79a7",
}


PATTERN_LABELS = {
    "dispersed bubbles": r"$\mathit{dispersed\ bubbles}$", "plug": r"$\mathit{plug}$",
    "slug": r"$\mathit{slug}$", "churn": r"$\mathit{churn}$",
    "wavy annular": r"$\mathit{wavy\ annular}$",
}


MARKERS = {
    "dispersed bubbles": "o", "plug": "o", "slug": "o", "churn": "o",
    "wavy annular": "o",
}


def hull_indices(points):
    """Return the observed convex polygon; reject fewer than three/collinear points."""
    unique = {}
    for index, (x, y) in enumerate(points):
        unique.setdefault((x, y), index)
    ordered = [(x, y, index) for (x, y), index in sorted(unique.items())]
    if len(ordered) < 3:
        return []

    def cross(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    lower, upper = [], []
    for point in ordered:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    for point in reversed(ordered):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    hull = lower[:-1] + upper[:-1]
    return [point[2] for point in hull] if len(hull) >= 3 else []


def save_figure(fig, basename):
    fig.savefig(artifact_path(basename + ".png"), dpi=180, facecolor="white")
    fig.savefig(artifact_path(basename + ".svg"), facecolor="white")
    plt.close(fig)


def portuguese_ticks(ax):
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_formatter(FuncFormatter(lambda value, position: f"{value:g}".replace(".", ",")))


def alternation_chart(rows, metric, basename, center=False, show_centroids=False):
    """Show eligible observations only; retain all input rows in the provenance CSV."""
    displayed = [r for r in rows if r["plotted"]]
    assert displayed and all(r["eligible_pure_unsaturated"] and not r["saturated"] for r in displayed)
    fig, ax = plt.subplots(figsize=(13.7, 7.3))
    fig.subplots_adjust(left=0.11, right=0.77, top=0.88, bottom=0.15)
    ymax = max(float(r[metric]) for r in displayed) * 1.10
    ax.set_xlim(-0.025, 1.03)
    ax.set_ylim(-0.025 * ymax, ymax)
    ax.set_xticks(np.arange(0, 1.01, 0.1))
    ax.grid(color="#d5d5d5", linestyle="--", linewidth=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    hull_checks = {}
    centroid_checks = {}
    for pattern in PATTERNS:
        selected = [r for r in displayed if r["supplied_pattern"] == pattern]
        points = [(float(r["conditional_circle_alpha_spatial_mean"]), float(r[metric]))
                  for r in selected]
        hull = hull_indices(points)
        if hull:
            polygon = np.asarray([points[i] for i in hull])
            ax.fill(polygon[:, 0], polygon[:, 1], facecolor=COLORS[pattern],
                    edgecolor=COLORS[pattern], linewidth=1.4, alpha=0.16, zorder=1)
        hull_checks[pattern] = {
            "eligible_runs": [r["run"] for r in selected],
            "polygon_vertex_runs": [selected[i]["run"] for i in hull],
            "shaded": bool(hull),
        }
        for row in selected:
            ax.scatter(float(row["conditional_circle_alpha_spatial_mean"]),
                       float(row[metric]), marker=MARKERS[pattern], s=83,
                       facecolors=COLORS[pattern], edgecolors="#222222",
                       linewidths=0.75, zorder=4)
        if show_centroids:
            # Each displayed acquisition has equal weight, regardless of its
            # panorama length. Compute the center of the points themselves.
            x_mean, y_mean = np.mean(np.asarray(points), axis=0)
            centroid_checks[pattern] = {
                "eligible_runs": [r["run"] for r in selected],
                "n": len(selected),
                "alpha_mean": float(x_mean),
                "flip_rate_per1000px_mean": float(y_mean),
                "weighting": "equal weight per eligible acquisition",
            }
            # White outline separates the star from nearby observations.
            # With n=1 the star correctly coincides with its sole observation.
            ax.scatter(x_mean, y_mean, marker="*", s=365,
                       facecolors=COLORS[pattern], edgecolors="white",
                       linewidths=3.1, zorder=6)
            ax.scatter(x_mean, y_mean, marker="*", s=270,
                       facecolors=COLORS[pattern], edgecolors="#222222",
                       linewidths=0.8, zorder=7,
                       gid="centroid_" + pattern.replace(" ", "_"))
    handles = []
    for pattern in PATTERNS:
        n = sum(r["supplied_pattern"] == pattern for r in displayed)
        handles.append(Line2D([], [], color="none", marker=MARKERS[pattern],
                              markerfacecolor=COLORS[pattern], markeredgecolor="#222222",
                              markersize=8, label=f"{PATTERN_LABELS[pattern]} (n={n})"))
    if show_centroids:
        handles.append(Line2D([], [], color="none", marker="*",
                              markerfacecolor="#bbbbbb", markeredgecolor="#222222",
                              markersize=13, label="Centroide (média)"))
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.01),
              title="Padrão informado", fontsize=10.2,
              title_fontsize=10.2, frameon=False, labelspacing=0.85)
    ax.set_xlabel("Fração de vazio estimada, α [−]", fontsize=13)
    ax.set_ylabel("Flip-rate na linha y = 40 (transições por 1.000 pixels)" if center else
                  "Flip-rate espacial\n(transições por 1.000 pares de pixels)", fontsize=12)
    ax.tick_params(labelsize=11)
    portuguese_ticks(ax)
    fig.suptitle("Fração de vazio × flip-rate na linha y = 40" if center else
                 "Fração de vazio × flip-rate espacial", fontsize=18, y=0.97)
    save_figure(fig, basename)
    return hull_checks, centroid_checks
