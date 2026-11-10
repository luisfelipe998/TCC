#!/usr/bin/env python3
"""Native-frame visual evidence; no inferred phase velocities or resampling.

Fixed examples span the retained patterns and include an annular saturation
failure. Frame indices were fixed before viewing these newly generated sheets.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re

import numpy as np
from PIL import Image, ImageDraw, ImageFont

HERE = Path(__file__).resolve().parents[1]
DATA = HERE.parents[1] / "data_processed"
os.environ.setdefault("MPLCONFIGDIR", str(HERE / "intermediate" / "visual_mplcache"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RUNS = (
    "m07q01_001", "m05q01_001", "m05q05_001", "m05q20_001",
    "m04q80_001", "m03q40_001",
)
INDICES = (200, 202, 204, 206, 208, 210)


def main():
    for name in ("figures", "results"):
        (HERE / name).mkdir(parents=True, exist_ok=True)
    metadata = []
    fontpath = Path("/System/Library/Fonts/Supplemental/Arial.ttf")
    font = ImageFont.truetype(str(fontpath), 17) if fontpath.exists() else ImageFont.load_default()
    heading = ImageFont.truetype(str(fontpath), 20) if fontpath.exists() else ImageFont.load_default()
    for run in RUNS:
        folder = DATA / run
        marker, = folder.glob("__*.txt")
        pattern = marker.stem[2:].replace("_", " ")
        config = (folder / "_config.mii").read_text()
        fps = float(re.search(r"Film speed 1\s*:([0-9.]+)", config).group(1))
        files = sorted(folder.glob("ImgA*.jpg"))
        ids = [int(re.search(r"(\d+)\.jpg$", f.name).group(1)) for f in files]
        assert ids == list(range(501)), run
        stack = np.stack([np.asarray(Image.open(f).convert("L")) for f in files])
        assert stack.shape == (501, 80, 1016), (run, stack.shape)
        canvas = Image.new("RGB", (1176, 650), "white")
        draw = ImageDraw.Draw(canvas)
        draw.text((16, 12), f"{run} | {pattern} | {fps:g} frames/s", fill="black", font=heading)
        draw.text((16, 39), "Native grayscale pixels; unchanged brightness; time relative to frame 200", fill="black", font=font)
        for row, index in enumerate(INDICES):
            y = 75 + row * 90
            draw.text((12, y + 15), f"Frame {index}", fill="black", font=font)
            draw.text((12, y + 39), f"{(index-200)/fps*1000:.3f} ms", fill="black", font=font)
            canvas.paste(Image.fromarray(stack[index]).convert("RGB"), (150, y))
        canvas.save(HERE / "figures" / f"visual_sequence_{run}.png")
        profile = stack[:, 20:60].mean(axis=1)
        fig, ax = plt.subplots(figsize=(11, 5.5), layout="constrained")
        im = ax.imshow(profile, cmap="gray", vmin=0, vmax=255, aspect="auto", origin="upper",
                       extent=(-0.5, 1015.5, 500.5/fps*1000, -0.5/fps*1000),
                       interpolation="nearest")
        ax.set_title(f"{run}: {pattern} — raw central-band x–time image")
        ax.set_xlabel("Axial image coordinate (px)")
        ax.set_ylabel("Time after frame 0 (ms)")
        fig.colorbar(im, ax=ax, label="Mean stored grayscale, rows 20:60")
        fig.savefig(HERE / "figures" / f"visual_xt_{run}.png", dpi=150)
        plt.close(fig)
        metadata.append(dict(
            run=run, pattern=pattern, fps=fps, frames=len(files),
            first_last_duration_s=500/fps, contact_frame_indices=list(INDICES),
            contact_relative_times_s=[(i-200)/fps for i in INDICES],
            mean_core_fraction_ge250=float((stack[:, 20:60] >= 250).mean()),
            raw_dimensions=[80, 1016], xt_rows=[20, 60],
            xt_interpolation="nearest", brightness_range=[0, 255],
            sequence_sha256={files[i].name: hashlib.sha256(files[i].read_bytes()).hexdigest()
                             for i in INDICES},
            description="Selected visual examples; no segmentation, motion fitting, or phase-velocity inference.",
        ))
        print(run, fps, metadata[-1]["mean_core_fraction_ge250"], flush=True)
    (HERE / "results" / "visual_sequences_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
