"""Prepare the official York Urban test photos and paired line endpoints."""
import argparse
import hashlib
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from PIL import Image
from scipy.io import loadmat

from prepare import resize_photo
from score import ROOT, read_csv, write_csv


ARCHIVE_SHA256 = "09c922b9a87b6c9385fa8a8b20b9900fb39d43928baa963bd4790ab5a8d0afc5"


def prepare(archive, folder):
    evidence = ROOT / "docs/evidence/wireframe"
    rosters = (
        "panel",
        "replication-panel",
        "confirmation-panel",
        "prospective-panel",
        "remainder-panel",
        "topology-train-panel",
    )
    seen = {
        row["pixel_sha256"]
        for name in rosters
        for row in read_csv(evidence / (name + ".csv"))
    }
    with archive.open("rb") as stream:
        if hashlib.file_digest(stream, "sha256").hexdigest() != ARCHIVE_SHA256:
            raise ValueError("York Urban archive checksum differs")
    folder.mkdir(parents=True, exist_ok=True)
    rows = []
    with ZipFile(archive) as source:

        def mat(name, field):
            return loadmat(BytesIO(source.read("YorkUrbanDB/" + name)))[field]

        names = mat("Manhattan_Image_DB_Names.mat", "Manhattan_Image_DB_Names").ravel()
        indices = mat("ECCV_TrainingAndTestImageNumbers.mat", "testSetIndex").ravel()
        if len(indices) != 51 or len(set(indices)) != 51:
            raise ValueError("Unexpected York Urban test split")
        for idx in indices:
            image_id = str(names[idx - 1][0]).rstrip("/\\")
            prefix = image_id + "/" + image_id
            with Image.open(
                BytesIO(source.read("YorkUrbanDB/" + prefix + ".jpg"))
            ) as image:
                rgb = image.convert("RGB")
                small = resize_photo(rgb)
            digest = hashlib.sha256(small.tobytes()).hexdigest()
            if digest in seen:
                raise ValueError(f"Duplicate photo {image_id}")
            seen.add(digest)
            small.save(folder / (image_id + ".small.png"))
            lines = mat(prefix + "LinesAndVP.mat", "lines")
            # Each consecutive pair of rows forms one segment, at one-based pixel centers.
            resized = (lines - 0.5) * np.array(small.size) / rgb.size + 0.5
            np.savetxt(
                folder / (image_id + ".annot.csv"),
                resized.reshape(-1, 4),
                delimiter=",",
                header="x1,y1,x2,y2",
                comments="",
            )
            rows.append(
                dict(
                    id=image_id,
                    panel="york51",
                    width=small.width,
                    height=small.height,
                    pixel_sha256=digest,
                )
            )
    write_csv(folder / "panel.csv", rows, list(rows[0]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("folder", type=Path)
    args = parser.parse_args()
    prepare(args.archive, args.folder)
