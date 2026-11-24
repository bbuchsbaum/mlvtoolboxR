"""Prepare the frozen Toronto reserve only after the curve-rank gate passes.

The command has deliberately separate phases. ``photos`` consumes the pinned
photo archive and model weights but never reads artist vectors. ``references``
only adds the selected artist geometry after a prediction has been frozen.
Neither phase has a command-line roster argument: the 120-item reserve roster
is the fixed evidence file below.
"""

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
from io import BytesIO
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile

import numpy as np
from PIL import Image, ImageOps

from prepare import load_model, write_full_views
from score import ROOT, write_csv
from toronto_scenes import (
    ARCHIVE_SHA256,
    assign_splits,
    file_sha256,
    geometry,
    load_records,
    resize_photo,
)


ROSTER = ROOT / "docs/evidence/wireframe/curve-rank/reserve-test.csv"
PHOTO_INDEX = ROOT / "docs/evidence/wireframe/toronto/photos.csv"
SOURCES = ("fdog", "mlsd")


def read_rows(path):
    with Path(path).open(newline="") as source:
        return list(csv.DictReader(source))


def frozen_roster():
    rows = read_rows(ROSTER)
    required = {"id", "category", "panel", "fold", "pixel_sha256"}
    if not rows or set(rows[0]) != required:
        raise ValueError("reserve roster schema differs from the frozen evidence")
    if len(rows) != 120 or len({row["id"] for row in rows}) != len(rows):
        raise ValueError("expected 120 distinct frozen reserve IDs")
    if any(row["panel"] != "reserve-test" for row in rows):
        raise ValueError("reserve roster contains a non-reserve-test row")
    return rows


def photo_index():
    rows = read_rows(PHOTO_INDEX)
    required = {
        "id",
        "category",
        "panel",
        "rank",
        "split_sha256",
        "original_width",
        "original_height",
        "width",
        "height",
        "pixel_sha256",
    }
    if (
        not rows
        or set(rows[0]) != required
        or len({row["id"] for row in rows}) != len(rows)
    ):
        raise ValueError("Toronto photo index schema or IDs differ from its freeze")
    return {row["id"]: row for row in rows}


def checked_rows():
    roster = frozen_roster()
    index = photo_index()
    rows = []
    for row in roster:
        photo = index.get(row["id"])
        if photo is None:
            raise ValueError(
                f"reserve ID is absent from Toronto photo index: {row['id']}"
            )
        for field in ("category", "pixel_sha256"):
            if row[field] != photo[field]:
                raise ValueError(
                    f"reserve roster disagrees with photo index for {row['id']}: {field}"
                )
        if photo["panel"] != "reserve" or not 26 <= int(photo["rank"]) <= 45:
            raise ValueError(
                f"reserve roster is not rank 26--45 in its category: {row['id']}"
            )
        # Staged photo and geometry contracts retain the source ``reserve``
        # panel; ``reserve-test`` is an evaluation-roster-only alias.
        rows.append({**row, **photo})
    return rows


def write_mlsd_full_photo(path, image_id, pixels, model):
    """Copy the exact full-photo file emitted by the established M-LSD helper."""
    with TemporaryDirectory() as temporary:
        temporary = Path(temporary)
        # ``traced_views`` uses this same helper and unmodified RGB array.
        write_full_views(temporary, image_id, pixels, model)
        shutil.copyfile(temporary / f"{image_id}.segments.csv", path)


def trace_photo_drawings(stage, base, rows):
    plan = base / "reserve-contour-plan.csv"
    write_csv(
        plan,
        [dict(id=row["id"], method=source) for row in rows for source in SOURCES],
        ["id", "method"],
    )
    try:
        subprocess.run(
            [
                "Rscript",
                "--vanilla",
                "tools/wireframe/contour_drawings.R",
                str(stage),
                str(base),
                str(plan),
            ],
            cwd=ROOT,
            check=True,
        )
    finally:
        plan.unlink(missing_ok=True)


def prepare_photos(archive, weights, output):
    archive, weights, output = Path(archive), Path(weights), Path(output)
    rows = checked_rows()
    if file_sha256(archive) != ARCHIVE_SHA256:
        raise ValueError("Toronto archive checksum differs from the frozen archive")
    stage, base = output / "staged", output / "base"
    stage.mkdir(parents=True, exist_ok=True)
    base.mkdir(parents=True, exist_ok=True)
    model = load_model(weights)
    staged = []
    with ZipFile(archive) as source:
        names = set(source.namelist())
        for row in rows:
            member = "CP_images/" + row["id"]
            if member not in names:
                raise ValueError(f"archive does not contain {member}")
            with Image.open(BytesIO(source.read(member))) as image:
                rgb = ImageOps.exif_transpose(image).convert("RGB")
            expected_original = (
                int(row["original_width"]),
                int(row["original_height"]),
            )
            if rgb.size != expected_original:
                raise ValueError(f"archive dimensions differ for {row['id']}")
            small = resize_photo(rgb)
            expected_working = (int(row["width"]), int(row["height"]))
            if small.size != expected_working:
                raise ValueError(f"working dimensions differ for {row['id']}")
            pixels = np.asarray(small)
            if hashlib.sha256(pixels.tobytes()).hexdigest() != row["pixel_sha256"]:
                raise ValueError(f"fixed RGB checksum differs for {row['id']}")
            small.save(stage / f"{row['id']}.small.png")
            write_mlsd_full_photo(
                base / f"{row['id']}.mlsd.segments.csv", row["id"], pixels, model
            )
            staged.append(
                {
                    key: row[key]
                    for key in (
                        "id",
                        "panel",
                        "category",
                        "width",
                        "height",
                        "pixel_sha256",
                    )
                }
            )
    write_csv(
        stage / "panel.csv",
        staged,
        ["id", "panel", "category", "width", "height", "pixel_sha256"],
    )
    trace_photo_drawings(stage, base, rows)


def prepare_references(archive, weights, output):
    """Add artist geometry only after frozen photo-only predictions exist.

    ``archive`` and ``weights`` remain positional for a symmetric, auditable CLI;
    reference export intentionally does not open either one.
    """
    del archive, weights
    output = Path(output)
    rows = checked_rows()
    stage, base = output / "staged", output / "base"
    staged = {row["id"]: row for row in read_rows(stage / "panel.csv")}
    if set(staged) != {row["id"] for row in rows}:
        raise ValueError(
            "staged photo IDs do not exactly match the frozen reserve roster"
        )
    for row in rows:
        actual = staged[row["id"]]
        for field in ("panel", "category", "width", "height", "pixel_sha256"):
            if actual.get(field) != row[field]:
                raise ValueError(
                    f"staged photo metadata differs for {row['id']}: {field}"
                )
    predictions = output / "predicted"
    selection = read_rows(predictions / "selection.csv")
    if len(selection) != len(rows) or {row["id"] for row in selection} != set(staged):
        raise ValueError(
            "Reference scoring requires a prediction for each frozen photo"
        )
    for item in selection:
        path = predictions / (item["id"] + ".curve_rank.segments.csv")
        if file_sha256(path) != item["segments_sha256"]:
            raise ValueError(f"Frozen prediction differs for {item['id']}")
    records = {record["id"]: record for record in assign_splits(load_records())}
    if not {row["id"] for row in rows} <= set(records):
        raise ValueError("artist vectors do not cover the frozen reserve")
    base.mkdir(parents=True, exist_ok=True)
    for row in rows:
        record = records[row["id"]]
        doc = geometry(
            record,
            (int(row["original_width"]), int(row["original_height"])),
            (int(row["width"]), int(row["height"])),
            "CP_images/" + row["id"],
        )
        doc["working"]["pixel_sha256"] = row["pixel_sha256"]
        (stage / f"{row['id']}.json").write_text(
            json.dumps(doc, separators=(",", ":"), allow_nan=False) + "\n"
        )
    plan = base / "reserve-artist-plan.csv"
    write_csv(
        plan, [dict(id=row["id"], method="artist") for row in rows], ["id", "method"]
    )
    try:
        subprocess.run(
            [
                "Rscript",
                "--vanilla",
                "tools/wireframe/contour_drawings.R",
                str(stage),
                str(base),
                str(plan),
            ],
            cwd=ROOT,
            check=True,
        )
    finally:
        plan.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("photos", "references"))
    parser.add_argument("archive", type=Path)
    parser.add_argument("weights", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    if args.mode == "photos":
        prepare_photos(args.archive, args.weights, args.output_dir)
    else:
        prepare_references(args.archive, args.weights, args.output_dir)


if __name__ == "__main__":
    main()
