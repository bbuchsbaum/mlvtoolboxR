"""Cache the frozen context scorer's photo-only inputs for the contour experiment."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
from PIL import Image
from photo_junctions import contour_prior, scored_candidates, traced_views
from predict_hawp import load_predictor
from prepare import load_model
from score import positions, read_csv


def cache(pixels, mlsd, predict, fitted):
    segments, views = traced_views(pixels, mlsd)
    prior = contour_prior(pixels, views, "context")
    xy, prob = scored_candidates(pixels, predict, fitted, "context")
    return dict(
        xy=xy,
        prob=prob,
        mlsd=segments,
        prior_xy=positions(prior),
        prior_type=np.array([r["type"] for r in prior], dtype="U1"),
    )


def run(args):
    fitted = json.loads(
        gzip.decompress(Path(__file__).with_name("context-scorer.json.gz").read_bytes())
    )
    mlsd = load_model(args.assets / "mlsd.tflite")
    predict = load_predictor(
        args.assets / "hawp-src", args.assets / "hawpv2-edb9b23f.pth"
    )
    args.output.mkdir(parents=True, exist_ok=True)
    rows = read_csv(args.panel)
    for row in rows[args.shard :: args.shards]:
        dest = args.output / (row["id"] + ".npz")
        with Image.open(args.images / (row["id"] + ".small.png")) as image:
            rgb = np.asarray(image.convert("RGB"))
        digest = hashlib.sha256(rgb.tobytes()).hexdigest()
        if row.get("pixel_sha256") and digest != row["pixel_sha256"]:
            raise ValueError(f"Pixel hash mismatch at {row['id']}")
        if dest.exists():
            with np.load(dest) as prior:
                if digest != str(prior["pixel_sha256"]):
                    raise ValueError(
                        f"Cached inputs belong to another photograph at {dest}"
                    )
            continue
        inputs = cache(rgb, mlsd, predict, fitted)
        tmp = dest.with_suffix(".partial.npz")
        np.savez_compressed(tmp, pixel_sha256=digest, **inputs)
        tmp.replace(dest)
        print(row["id"], flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("assets", "panel", "images", "output"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--shards", type=int, default=1)
    run(parser.parse_args())
