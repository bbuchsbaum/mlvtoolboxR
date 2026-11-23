"""Predict straight segments from four overlapping crops of each resized RGB image."""
import argparse
import csv
import hashlib
from pathlib import Path

import numpy as np
from PIL import Image

from prepare import load_model, run_model


def write_tiles(path, pixels, model):
    h, w = pixels.shape[:2]
    cw, ch = (2 * w + 2) // 3, (2 * h + 2) // 3
    with path.open("w", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["left", "top", "width", "height", "x1", "y1", "x2", "y2"])
        for top in (0, h - ch):
            for left in (0, w - cw):
                crop = pixels[top : top + ch, left : left + cw]
                for line in run_model(model, crop):
                    writer.writerow([left, top, cw, ch, *(line + 1)])


def predict(folder):
    model = load_model(folder / "mlsd.tflite")
    with (folder / "panel.csv").open() as roster:
        for row in csv.DictReader(roster):
            image_id = row["id"]
            with Image.open(folder / f"{image_id}.small.png") as image:
                pixels = np.array(image.convert("RGB"))
            if hashlib.sha256(pixels.tobytes()).hexdigest() != row["pixel_sha256"]:
                raise ValueError(f"Resized RGB checksum mismatch for {image_id}")
            write_tiles(folder / f"{image_id}.tiles.csv", pixels, model)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    predict(parser.parse_args().directory)
