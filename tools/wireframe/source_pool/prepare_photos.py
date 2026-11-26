"""Prepare the two screened sources' photo inputs for the development photos.

Reads only the development roster's members of the Toronto archive. For each
photo it checks the staged working-image hash, writes a 512-pixel LANCZOS
resize for FDoG, and runs OpenCV's structured edge detector on the original
photo with the pinned BSDS model. Outputs are private pixels and edge maps.

Run with the dedicated environment that has opencv-contrib-python-headless:

    .agent-work/venv-sed/bin/python -I tools/wireframe/source_pool/prepare_photos.py \
        .agent-work/toronto/CP_images.zip docs/evidence/wireframe/curve-rank/development.csv \
        .agent-work/source-pool/model/model.yml.gz .agent-work/source-pool
"""

import argparse
import csv
import hashlib
import json
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile

import cv2
import numpy as np
import PIL
from PIL import Image, ImageOps

MODEL_COMMIT = "5e3a56880fb115e757855f8e01e744c154791144"
MODEL_PATH = "testdata/cv/ximgproc/model.yml.gz"
MODEL_SHA256 = "fe3586bca5ec124b73fa769c2b54e195e4526d6ceeddf161fad96bf68a65bf17"
THRESHOLD = 0.1


def resize_long(rgb, long_side):
    """The staging convention of toronto_scenes.resize_photo, at any long side."""
    size = tuple(max(1, round(d * long_side / max(rgb.size))) for d in rgb.size)
    return rgb.resize(size, Image.Resampling.LANCZOS)


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("roster", type=Path)
    parser.add_argument("model", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if sha256(args.model.read_bytes()) != MODEL_SHA256:
        raise ValueError("structured edge model checksum differs from the pinned file")
    with args.roster.open() as f:
        roster = list(csv.DictReader(f))
    if len(roster) != 60 or {r["panel"] for r in roster} != {"development"}:
        raise ValueError("expected the 60-photo development roster")
    (args.output / "fdog512").mkdir(parents=True, exist_ok=True)
    (args.output / "sed").mkdir(parents=True, exist_ok=True)
    detector = cv2.ximgproc.createStructuredEdgeDetection(str(args.model))
    rows = []
    with ZipFile(args.archive) as archive:
        for r in roster:
            name = r["id"]
            with Image.open(BytesIO(archive.read("CP_images/" + name))) as image:
                rgb = ImageOps.exif_transpose(image).convert("RGB")
            if rgb.size != (800, 600):
                raise ValueError(f"{name}: expected 800 x 600, found {rgb.size}")
            small = resize_long(rgb, 256)
            if small.size != (256, 192) or sha256(small.tobytes()) != r["pixel_sha256"]:
                raise ValueError(
                    f"{name}: working image does not reproduce the staged hash"
                )
            large = resize_long(rgb, 512)
            if large.size != (512, 384):
                raise ValueError(f"{name}: unexpected 512 resize {large.size}")
            large.save(args.output / "fdog512" / f"{name}.512.png")
            # The detector expects RGB channel order, float32 in [0, 1]. Pillow
            # already gives RGB, so no BGR conversion is applied.
            src = np.asarray(rgb, dtype=np.float32) / np.float32(255)
            edges = detector.detectEdges(src)
            orientation = detector.computeOrientation(edges)
            nms = detector.edgesNms(edges, orientation)  # module defaults r=2, s=0, m=1
            nms = np.ascontiguousarray(nms, dtype="<f4")
            mask = nms >= THRESHOLD
            Image.fromarray(mask.astype(np.uint8) * 255).save(
                args.output / "sed" / f"{name}.mask.png"
            )
            (args.output / "sed" / f"{name}.strength.f32").write_bytes(nms.tobytes())
            rows.append(
                dict(
                    id=name,
                    category=r["category"],
                    working_sha256_ok=1,
                    sed_raw_max=float(edges.max()),
                    sed_nms_nonzero=int((nms > 0).sum()),
                    sed_mask_pixels=int(mask.sum()),
                    sed_strength_sha256=sha256(nms.tobytes()),
                    fdog512_sha256=sha256(large.tobytes()),
                )
            )
    with (args.output / "prepare-manifest.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    meta = dict(
        model_repository="https://github.com/opencv/opencv_extra",
        model_commit=MODEL_COMMIT,
        model_path=MODEL_PATH,
        model_sha256=MODEL_SHA256,
        opencv=cv2.__version__,
        numpy=np.__version__,
        pillow=PIL.__version__,
        nms="edgesNms defaults: r=2, s=0, m=1",
        threshold=THRESHOLD,
        input="original 800 x 600 RGB after EXIF transpose, float32 / 255, RGB order",
        photos=len(rows),
    )
    (args.output / "prepare-meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
