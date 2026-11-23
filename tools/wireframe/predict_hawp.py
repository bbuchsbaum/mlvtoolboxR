"""Run the supervised HAWPv2 checkpoint on the benchmark's resized RGB photos."""
import argparse
import csv
import hashlib
import sys
from pathlib import Path

import numpy as np
import torch
from PIL import Image
from torch.nn import functional as F

WEIGHTS_SHA256 = "c9548d67dd1a60a4c3c0271eaf249adc234409d8c5512eaaa18a66ba64d2f670"


def load_predictor(source, weights):
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

    @torch.inference_mode()
    def lines(rgb, *, with_context=False):
        height, width = rgb.shape[:2]
        args = (
            transform(rgb).unsqueeze(0),
            [dict(filename="input", width=width, height=height)],
        )
        if with_context:
            # Capture only the full-photo pass. Crop predictions use the ordinary path.
            captured = []
            with model.backbone.register_forward_hook(
                lambda m, args, result: captured.append(result)
            ):
                output, _ = model(*args)
            heads, features = captured[0]
            head = heads[0]
            junction = head[:, 6:7] - head[:, 5:6]
            thin = model.fc3(features)
            maps = torch.cat(
                [
                    model.fc1(features),
                    thin,
                    head,
                    F.max_pool2d(junction, 5, 1, 2),
                    F.avg_pool2d(thin, 11, 1, 5),
                ],
                dim=1,
            )
        else:
            output, _ = model(*args)
        segments = np.column_stack(
            [output["lines_pred"].numpy() + 1, output["lines_score"].numpy()]
        )
        return (segments, maps) if with_context else segments

    return lines


def sample_context(maps, xy, size):
    """Bilinearly sample activations at one-based pixel centers, as during fitting."""
    grid = 2 * (xy - 0.5) / np.asarray(size) - 1
    grid = torch.as_tensor(grid, dtype=torch.float32).reshape(1, 1, -1, 2)
    sampled = F.grid_sample(maps, grid, align_corners=False, padding_mode="border")
    return sampled[0, :, 0, :].T.numpy().astype(np.float16)


def crop_lines(rgb, predict):
    from prepare import enhance_contrast

    # Give faint short edges both a local contrast adjustment and a larger
    # effective scale. All crops originate from the same benchmark pixels.
    image = enhance_contrast(rgb)
    h, w = rgb.shape[:2]
    cw, ch = (2 * w + 2) // 3, (2 * h + 2) // 3
    all_lines = []
    for top in sorted({0, h - ch}):
        for left in sorted({0, w - cw}):
            rows = predict(image[top : top + ch, left : left + cw])
            keep = rows[:, 4] > 0.06
            for axis, cut, size, extent in [(0, left, cw, w), (1, top, ch, h)]:
                ends = rows[:, [axis, axis + 2]]
                if cut:
                    keep &= ends.min(axis=1) > 2
                if cut + size < extent:
                    keep &= ends.max(axis=1) < size - 1
            keep_rows = rows[keep].copy()
            keep_rows[:, :4] += [left, top, left, top]
            all_lines.append(keep_rows)
    return np.vstack(all_lines)


def predict(folder, source, weights, tiles=False):
    lines = load_predictor(source, weights)
    with (folder / "panel.csv").open() as roster:
        for photo in csv.DictReader(roster):
            image_id = photo["id"]
            with Image.open(folder / f"{image_id}.small.png") as image:
                rgb = np.array(image.convert("RGB"))
            if hashlib.sha256(rgb.tobytes()).hexdigest() != photo["pixel_sha256"]:
                raise ValueError(f"Resized RGB checksum mismatch for {image_id}")
            result = crop_lines(rgb, lines) if tiles else lines(rgb)
            suffix = "hawp-tiles" if tiles else "hawp"
            np.savetxt(
                folder / f"{image_id}.{suffix}.csv",
                result,
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
    parser.add_argument("--tiles", action="store_true")
    args = parser.parse_args()
    predict(args.directory, args.source, args.weights, args.tiles)
