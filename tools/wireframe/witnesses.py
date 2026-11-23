"""Render the source, human labels, masks and typed contours at two witness crops."""
import argparse
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from score import read_csv

# Four locations printed in issue 3 and four floor intersections from our labels.
CROPS = [
    (
        "00079423",
        (175, 56),
        [(170.2, 51.4), (178.8, 51.6), (170.1, 60.5), (178.7, 60.7)],
    ),
    (
        "00078318",
        (224, 139),
        [(225.4, 134.1), (224.1, 141.3), (204.9, 137.2), (244.6, 138.3)],
    ),
]
STAGES = [
    ("Resized photograph", "source"),
    ("Human line annotations", "labels"),
    ("FDoG gap 6 mask", "fdog6.mask"),
    ("FDoG gap 6 traces", "fdog6.traces"),
    ("M-LSD mask", "mlsd_mask.mask"),
    ("M-LSD traces", "mlsd_mask.traces"),
    ("M-LSD with local contrast", "mlsd_clahe.traces"),
    ("Traced human lines", "annotation_mask.traces"),
]


def render(folder):
    detections = read_csv(folder / "junctions.csv")
    canvas = Image.new("RGB", (198 * len(STAGES), 447), "white")
    labels = ImageDraw.Draw(canvas)
    for row, (image_id, (cx, cy), points) in enumerate(CROPS):
        # PIL boxes are zero-based and right/bottom exclusive.
        box = (cx - 24, cy - 22, cx + 22, cy + 20)

        def screen(p):
            return tuple((np.asarray(p) - 1 - box[:2]) * 4 + 1.5)

        top = row * 222
        labels.text(
            (2, top),
            image_id
            + "    Coral rings: reference X.    Green: detected X.    Yellow: detected T.",
            fill="black",
        )
        for col, (title, stage) in enumerate(STAGES):
            if stage == "source" or stage.endswith(".mask"):
                suffix = "small.png" if stage == "source" else stage + ".png"
                with Image.open(folder / f"{image_id}.{suffix}") as raw:
                    tile = (
                        raw.convert("RGB")
                        .crop(box)
                        .resize((184, 168), Image.Resampling.NEAREST)
                    )
            else:
                tile = Image.new("RGB", (184, 168), (24, 26, 30))
            draw = ImageDraw.Draw(tile)
            if stage == "labels" or stage.endswith(".traces"):
                suffix = "annot.csv" if stage == "labels" else stage + ".csv"
                lines = np.loadtxt(
                    folder / f"{image_id}.{suffix}", delimiter=",", skiprows=1
                )
                for a, b in lines.reshape(-1, 2, 2):
                    draw.line([screen(a), screen(b)], fill="white", width=2)
            for point in points:
                x, y = screen(point)
                draw.ellipse((x - 6, y - 6, x + 6, y + 6), outline="#ff8267")
            for junction in detections:
                if (
                    stage == junction["method"] + ".traces"
                    and junction["id"] == image_id
                    and junction["type"] in ("T", "X")
                ):
                    x, y = screen([float(junction["x"]), float(junction["y"])])
                    color = "lime" if junction["type"] == "X" else "yellow"
                    draw.ellipse((x - 3, y - 3, x + 3, y + 3), fill=color)
            left = 198 * col
            labels.text((left + 2, top + 20), title, fill="black")
            canvas.paste(tile, (left, top + 38))
    canvas.save(folder / "witness-inspection.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    render(parser.parse_args().directory)
