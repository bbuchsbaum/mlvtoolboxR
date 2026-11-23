"""Run the supervised HAWPv2 checkpoint on the benchmark's resized RGB photos."""
import argparse
import csv
import hashlib
import sys
from pathlib import Path

import numpy as np
import torch
from PIL import Image

WEIGHTS_SHA256 = "c9548d67dd1a60a4c3c0271eaf249adc234409d8c5512eaaa18a66ba64d2f670"


def predict(folder, source, weights):
    if hashlib.sha256(weights.read_bytes()).hexdigest() != WEIGHTS_SHA256:
        raise ValueError("HAWPv2 weight checksum mismatch")
    sys.path.insert(0, str(source.resolve()))
    from hawp.fsl.config import cfg
    from hawp.fsl.dataset.build import build_transform
    from hawp.fsl.model.build import build_model

    cfg.merge_from_file(str(source / "configs/hawpv2.yaml"))
    cfg.MODEL.DEVICE = "cpu"
    torch.set_num_threads(4)
    model = build_model(cfg)
    model.use_residual = 2  # Official HAWPv2 predict and benchmark CLI default.
    model.load_state_dict(
        torch.load(weights, map_location="cpu", weights_only=True)["model"]
    )
    model.eval()
    transform = build_transform(cfg)
    with (folder / "panel.csv").open() as roster:
        photos = list(csv.DictReader(roster))
    with torch.inference_mode():
        for photo in photos:
            image_id = photo["id"]
            with Image.open(folder / f"{image_id}.small.png") as image:
                rgb = np.array(image.convert("RGB"))
            if hashlib.sha256(rgb.tobytes()).hexdigest() != photo["pixel_sha256"]:
                raise ValueError(f"Resized RGB checksum mismatch for {image_id}")
            height, width = rgb.shape[:2]
            meta = dict(filename=image_id, width=width, height=height)
            output, _ = model(transform(rgb).unsqueeze(0), [meta])
            # Upstream scales back to the input photo and plots on zero-based centers.
            coords = output["lines_pred"].numpy() + 1
            scores = output["lines_score"].numpy()
            np.savetxt(
                folder / f"{image_id}.hawp.csv",
                np.column_stack([coords, scores]),
                delimiter=",",
                header="x1,y1,x2,y2,score",
                comments="",
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument(
        "source", type=Path, help="Pinned HAWP checkout with the CPU patch"
    )
    parser.add_argument("weights", type=Path)
    args = parser.parse_args()
    predict(args.directory, args.source, args.weights)
