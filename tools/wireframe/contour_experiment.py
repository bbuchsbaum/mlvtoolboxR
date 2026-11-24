"""Build four photo-only contour treatments and trace their full-scene drawings."""

import argparse
import gzip
import hashlib
import json
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from PIL import Image

from contour_paths import resolve_all_modes
from contour_render import result_segments
from score import ROOT, read_csv, write_csv


def photo_outputs(photo, source, inputs, output, cutoffs, york=False):
    name = photo["id"]
    with Image.open(source / (name + ".small.png")) as image:
        rgb = np.asarray(image.convert("RGB"))
    digest = hashlib.sha256(rgb.tobytes()).hexdigest()
    with np.load(inputs / (name + ".npz")) as record:
        if digest != photo["pixel_sha256"] or digest != str(record["pixel_sha256"]):
            raise ValueError(f"Photo hash mismatch for {name}")
        lines = record["mlsd"]
        prior = [
            dict(x=x, y=y, type=t)
            for (x, y), t in zip(record["prior_xy"], record["prior_type"], strict=True)
        ]
        variants = resolve_all_modes(
            record["xy"], record["prob"], prior, cutoffs, rgb, lines
        )
    candidates = variants["audit"]["both"]["candidates"]
    ts = [r for r in candidates if r["type"] == "T"]
    xs = [r for r in candidates if r["type"] == "X"]
    audit = dict(
        T_candidates=len(ts),
        T_with_stem=sum("track_id" in r["stem_path"] for r in ts),
        T_vetoed=sum(r["stem_path"].get("veto", False) for r in ts),
        X_candidates=len(xs),
        X_with_two_tracks=sum(len(r.get("track_signature", ())) == 2 for r in xs),
    )
    (output / (name + ".audit.json")).write_text(json.dumps(audit) + "\n")
    keep_keys = ("x", "y", "type", "confidence")
    detections = {}
    masks = {}
    for mode in ("baseline", "stem", "tracks", "both"):
        method = "context" if mode == "baseline" else mode
        events = variants[mode]
        rejected = [
            r
            for r in variants["audit"][mode]["candidates"]
            if r.get("selection") == "rejected_contour_gap"
        ]
        detections[method] = [{key: r[key] for key in keep_keys} for r in events]
        comparator = variants["tracks" if method == "both" else "baseline"]
        counterfactual = {(r["x"], r["y"]) for r in comparator if r["type"] == "T"}
        cuts = [r for r in rejected if (r["x"], r["y"]) in counterfactual]
        masks[method] = result_segments(lines, events, cuts)
    if not york:
        masks["mlsd"] = lines
    for method, mask in masks.items():
        np.savetxt(
            output / f"{name}.{method}.segments.csv",
            mask,
            delimiter=",",
            header="x1,y1,x2,y2",
            comments="",
        )
    (output / (name + ".json")).write_text(
        json.dumps(detections, separators=(",", ":")) + "\n"
    )
    return [dict(id=name, method=m) for m in masks]


def evaluate(args):
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    cutoffs = json.loads(
        gzip.decompress(Path(__file__).with_name("context-scorer.json.gz").read_bytes())
    )["cutoffs"]
    rows = [
        r
        for r in read_csv(args.source / "panel.csv")
        if args.panel == "all" or args.panel == r.get("panel", r.get("split"))
    ]
    rows = rows[args.shard :: args.shards]
    records = []
    for r in rows:
        records.extend(
            photo_outputs(r, args.source, args.inputs, output, cutoffs, args.york)
        )
    if not args.york:
        records.extend(dict(id=r["id"], method="artist") for r in rows)
    if not args.york and not args.skip_fdog:
        records.extend(dict(id=r["id"], method="fdog") for r in rows)
    # Many configurations retain exactly the same segments. Trace identical
    # pixels once, retaining distinct scores for all four experimental arms.
    masks = {}
    jobs = []
    copies = []
    for row in records:
        prefix = output / (row["id"] + "." + row["method"])
        if row["method"] in ("artist", "fdog"):
            jobs.append(row)
            continue
        key = (
            row["id"],
            hashlib.sha256(Path(str(prefix) + ".segments.csv").read_bytes()).digest(),
        )
        if key in masks:
            copies.append((masks[key], prefix))
        else:
            jobs.append(row)
            masks[key] = prefix
    if jobs:
        import shutil

        with TemporaryDirectory() as temp:
            plan = Path(temp) / "plan.csv"
            write_csv(plan, jobs, ["id", "method"])
            subprocess.run(
                [
                    "Rscript",
                    "--vanilla",
                    "tools/wireframe/contour_drawings.R",
                    str(args.source.resolve()),
                    str(output.resolve()),
                    str(plan),
                ],
                cwd=ROOT,
                check=True,
            )
        for src, dst in copies:
            for suffix in (".drawing.png", ".sites.csv", ".contours.json"):
                shutil.copyfile(str(src) + suffix, str(dst) + suffix)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "inputs", "output"):
        p.add_argument(name, type=Path)
    p.add_argument("--panel", default="all")
    p.add_argument("--york", action="store_true")
    p.add_argument(
        "--skip-fdog",
        action="store_true",
        help="FDoG drawings have already been traced separately",
    )
    p.add_argument("--shard", type=int, default=0)
    p.add_argument("--shards", type=int, default=1)
    evaluate(p.parse_args())
