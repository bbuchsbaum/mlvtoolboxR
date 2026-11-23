"""Prepare a deterministic subset of the official Wireframe training split."""
import argparse
import csv
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from PIL import Image, ImageOps

from prepare import DATA, download, resize_photo


def prepare(folder):
    folder.mkdir(parents=True, exist_ok=True)
    url = DATA.removesuffix("/test") + "/train"
    meta = download(folder / "metadata.jsonl", url + "/metadata.jsonl")
    if (
        hashlib.sha256(meta.read_bytes()).hexdigest()
        != "2d08abcd824cfbe43b0de9b72ab395f599d1c447610f963f32ca9fee1d5a8c46"
    ):
        raise ValueError("Training annotations have changed")
    records = sorted(
        map(json.loads, meta.read_text().splitlines()),
        key=lambda r: hashlib.sha256(
            ("junction-arms-v1/" + r["image_file"]).encode()
        ).digest(),
    )[:320]

    def fetch(record):
        name = record["image_file"]
        i = Path(name).stem
        path = download(folder / name, url + "/images/" + name)
        with Image.open(path) as image:
            rgb = ImageOps.exif_transpose(image).convert("RGB")
        if list(rgb.size) != record["image_size"]:
            raise ValueError(f"Image dimension mismatch for {name}")
        small = resize_photo(rgb)
        small.save(folder / (i + ".small.png"))
        lines = (np.asarray(record["lines"]) + 0.5) * (
            np.array(small.size) / rgb.size
        ) + 0.5
        np.savetxt(
            folder / (i + ".annot.csv"),
            lines.reshape(-1, 4),
            delimiter=",",
            header="x1,y1,x2,y2",
            comments="",
        )
        return [i, "train", *small.size, hashlib.sha256(small.tobytes()).hexdigest()]

    with ThreadPoolExecutor(max_workers=6) as pool:
        rows = list(pool.map(fetch, records))
    with (folder / "panel.csv").open("w") as output:
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(["id", "panel", "width", "height", "pixel_sha256"])
        writer.writerows(rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    prepare(parser.parse_args().folder)
