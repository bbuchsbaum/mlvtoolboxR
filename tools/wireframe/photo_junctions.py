"""Write experimental photographic T/X detections for a single local image."""
import argparse
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from PIL import Image, ImageOps

from evaluate_support import extend, rank_crops
from fragment_crossings import crossings
from predict_tiles import write_tiles
from prepare import load_model, resize_photo, write_full_views
from score import ROOT, read_csv, write_csv


def detect(image, weights, method):
    with Image.open(image) as source:
        rgb = ImageOps.exif_transpose(source).convert("RGB")
    small = resize_photo(rgb)
    pixels = np.asarray(small)
    model = load_model(weights)
    with TemporaryDirectory() as temp:
        folder = Path(temp)
        write_full_views(folder, "photo", pixels, model)
        write_tiles(folder / "photo.tiles.csv", pixels, model)
        write_csv(
            folder / "panel.csv",
            [dict(id="photo", width=small.width, height=small.height)],
            ["id", "width", "height"],
        )
        subprocess.run(
            [
                "Rscript",
                "--vanilla",
                "-e",
                'pkgload::load_all(".", quiet=TRUE); source("tools/wireframe/trace_tiles.R")',
                str(folder),
                "--include-full",
            ],
            check=True,
            cwd=ROOT,
        )
        views = read_csv(folder / "all-views.csv")
    primary = [r for r in views if r["method"] == "mlsd_mask"]
    enhanced = [r for r in views if r["method"] == "mlsd_clahe"]
    crops = [r for r in views if r["method"].startswith("crop_")]
    proposals = rank_crops(enhanced, crops)
    if method == "supported":
        fragments = [dict(x=x, y=y, type="X") for _, x, y in crossings(pixels)]
        proposals = fragments + proposals
    detections = extend(primary, proposals, method)
    scale = np.array(rgb.size) / small.size
    for event in detections:
        x, y = float(event["x"]), float(event["y"])
        original_x, original_y = (np.array([x, y]) - 0.5) * scale + 0.5
        yield dict(
            type=event["type"], x=x, y=y, original_x=original_x, original_y=original_y
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "image", type=Path, help="EXIF orientation is applied before resizing"
    )
    parser.add_argument(
        "output", type=Path, help="CSV of one-based pixel centers at both sizes"
    )
    parser.add_argument(
        "--weights",
        type=Path,
        required=True,
        help="Pinned M-LSD large FP32 TFLite file",
    )
    parser.add_argument(
        "--method",
        choices=("supported", "crops"),
        default="supported",
        help="Include four-arm X proposals (supported) or use only traced contours (crops)",
    )
    args = parser.parse_args()
    if args.output.resolve() in (args.image.resolve(), args.weights.resolve()):
        parser.error("Choose an output file distinct from the image and weights")
    rows = list(detect(args.image, args.weights, args.method))
    write_csv(args.output, rows, ["type", "x", "y", "original_x", "original_y"])
