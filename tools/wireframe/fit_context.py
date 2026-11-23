"""Refit the context scorer on the frozen 320-photo training roster."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from PIL import Image

from context_features import descriptors
from fit_topology import cache, fit_trees
from predict_hawp import load_predictor
from score import ROOT, read_csv


def fit(folder, source, weights, output):
    roster = read_csv(folder / "panel.csv")
    expected = read_csv(ROOT / "docs/evidence/wireframe/topology-train-panel.csv")
    if roster != expected:
        raise ValueError(
            "Training requires the frozen Wireframe roster in its original order"
        )
    predict = load_predictor(source, weights)
    records = []
    for row, (geometry, labels, hard) in zip(
        roster, cache(folder, "multiscale", training=True), strict=True
    ):
        i = row["id"]
        with Image.open(folder / (i + ".small.png")) as image:
            rgb = np.asarray(image.convert("RGB"))
        if hashlib.sha256(rgb.tobytes()).hexdigest() != row["pixel_sha256"]:
            raise ValueError(f"Training pixels differ for {i}")
        _, maps = predict(rgb, with_context=True)
        with np.load(folder / (i + ".multiscale.npz")) as record:
            f = descriptors(geometry, maps, record["xy"], rgb.shape[1::-1])
        records.append((f, labels, hard))
    f, y, hard = (np.concatenate(v) for v in zip(*records))
    models = fit_trees(f, y, hard, iterations=150, leaves=20, learning_rate=0.08)
    output.write_text(json.dumps(models, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument("hawp_source", type=Path)
    parser.add_argument("hawp_weights", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    fit(args.folder, args.hawp_source, args.hawp_weights, args.output)
