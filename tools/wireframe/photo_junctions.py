"""Write experimental photographic T/X detections for a single local image."""
import argparse
import gzip
import json
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


def detect(image, weights, method, hawp_source=None, hawp_weights=None):
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
    if method in ("supported", "topology", "context"):
        fragments = [dict(x=x, y=y, type="X") for _, x, y in crossings(pixels)]
        proposals = fragments + proposals
    detections = extend(primary, proposals, method)
    fitted_method = method in ("topology", "context")
    if fitted_method:
        from predict_hawp import load_predictor, crop_lines
        from topology import model_features, probabilities, resolve

        filename = (
            "context-scorer.json.gz"
            if method == "context"
            else "junction-scorer.json.gz"
        )
        with gzip.open(Path(__file__).with_name(filename), "rt") as f:
            fitted = json.load(f)
        predictor = load_predictor(hawp_source, hawp_weights)
        detail = (
            crop_lines(pixels, predictor) if fitted["variant"] == "multiscale" else None
        )
        if method == "context":
            from context_features import descriptors as context_descriptors

            lines, maps = predictor(pixels, with_context=True)
            xy, geometry = model_features(pixels, lines, fitted["variant"], detail)
            features = context_descriptors(geometry, maps, xy, small.size)
        else:
            xy, features = model_features(
                pixels, predictor(pixels), fitted["variant"], detail
            )
        detections = resolve(
            xy,
            probabilities(features, fitted["models"]),
            detections,
            fitted["cutoffs"],
        )
    scale = np.array(rgb.size) / small.size
    for event in detections:
        x, y = float(event["x"]), float(event["y"])
        original_x, original_y = (np.array([x, y]) - 0.5) * scale + 0.5
        row = dict(
            type=event["type"], x=x, y=y, original_x=original_x, original_y=original_y
        )
        if fitted_method:
            row["score"] = event["confidence"]
        yield row


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
        choices=("supported", "crops", "topology", "context"),
        default="supported",
        help="Use the capped baseline, contour-only crops, or a fitted topology or context scorer",
    )
    parser.add_argument(
        "--hawp-source", type=Path, help="Pinned HAWP source with CPU patch"
    )
    parser.add_argument(
        "--hawp-weights", type=Path, help="Supervised HAWPv2 checkpoint"
    )
    args = parser.parse_args()
    if args.method in ("topology", "context") and (
        not args.hawp_source or not args.hawp_weights
    ):
        parser.error(
            f"--method {args.method} requires --hawp-source and --hawp-weights"
        )
    if args.output.resolve() in {
        p.resolve() for p in [args.image, args.weights, args.hawp_weights] if p
    }:
        parser.error("Choose an output file distinct from the image and weights")
    rows = list(
        detect(
            args.image, args.weights, args.method, args.hawp_source, args.hawp_weights
        )
    )
    fields = ["type", "x", "y", "original_x", "original_y"]
    if args.method in ("topology", "context"):
        fields.append("score")
    write_csv(args.output, rows, fields)
