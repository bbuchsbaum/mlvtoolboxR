"""Select cutoffs on earlier photos, then evaluate frozen photographic predictions."""
import argparse
import json
import itertools
import gzip
import hashlib
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist

from evaluate_support import assess
from score import ROOT, match, positions, read_csv, score, supplement, totals, write_csv
from topology import anchor, probabilities, resolve, suppress


def check_rosters(work, include_remainder=False):
    used_ids = set()
    pixel_owner = {}
    excluded = set()
    test_ids = set()
    inputs = [
        (work / "topology-train", 320),
        (work, 47),
        (work / "replication", 24),
        (work / "confirmation", 48),
        (work / "prospective", 64),
    ]
    if include_remainder:
        inputs.append((work / "remainder", 279))
    for folder, count in inputs:
        rows = read_csv(folder / "panel.csv")
        ids = {r["id"] for r in rows}
        pixels = {r["pixel_sha256"] for r in rows}
        if (
            len(rows) != count
            or len(ids) != count
            or len(pixels) != count
            or ids & used_ids
        ):
            raise ValueError(
                f"Unexpected size, duplicate photo or overlapping pixels at {folder}"
            )
        overlap = {
            r["id"]: pixel_owner[r["pixel_sha256"]]
            for r in rows
            if r["pixel_sha256"] in pixel_owner
        }
        # The pinned test release repeats one confirmation photo under a new ID
        # and a different annotation. Score only genuinely new photographs.
        expected = {"00036808": "00037266"} if folder == work / "remainder" else {}
        if overlap != expected:
            raise ValueError(
                f"Unexpected cross-panel pixel overlap at {folder}: {overlap}"
            )
        excluded.update(overlap)
        used_ids.update(ids)
        pixel_owner.update({r["pixel_sha256"]: r["id"] for r in rows})
        if folder != work / "topology-train":
            test_ids.update(ids)
    if include_remainder:
        metadata = (work / "metadata.jsonl").read_text().splitlines()
        official_ids = {Path(json.loads(row)["image_file"]).stem for row in metadata}
        if test_ids != official_ids:
            raise ValueError(
                "Test rosters do not cover the official test split exactly"
            )
    manifest = work / "issue-3.json"
    if (
        hashlib.sha256(manifest.read_bytes()).hexdigest()
        != "a65df9d883d62e611a6affadac0e0929202db7e56193a8879e9797a751301f14"
    ):
        raise ValueError("Original issue manifest differs from the published freeze")
    return excluded


def freeze(output, variant, cutoffs, model):
    payload = json.dumps(
        dict(variant=variant, cutoffs=cutoffs, models=model), separators=(",", ":")
    )
    output.write_bytes(gzip.compress(payload.encode(), mtime=0))


def groups(folder, panel):
    reviews = ROOT / "docs/evidence/wireframe"
    if panel in ("prospective", "remainder"):
        photos = supplement(folder, (panel,))
        gaps = read_csv(reviews / f"{panel}-gap-review.csv")
    else:
        sites = read_csv(folder / f"support-{panel}-sites.csv")
        ids = dict.fromkeys(
            r["id"] for r in read_csv(folder / f"support-{panel}-scores.csv")
        )
        photos = [
            (
                i,
                panel,
                [
                    (r["kind"], [float(r["x"]), float(r["y"])])
                    for r in sites
                    if r["id"] == i
                ],
            )
            for i in ids
        ]
        gaps = [
            dict(r, include="1")
            for r in read_csv(folder / f"support-{panel}-near-miss-scores.csv")
        ]
    return photos, gaps


def prepared(folder, photos, models, variant):
    for i, _, _ in photos:
        with np.load(folder / (i + "." + variant + ".npz")) as data:
            yield i, (data["xy"], probabilities(data["f"], models))


def select(work, model, variant, anchored):
    truth = {}
    predictions = {}
    limits = np.zeros(2, dtype=int)
    for panel, folder in [
        ("development", work),
        ("replication", work / "replication"),
        ("confirmation", work / "confirmation"),
    ]:
        photos, gaps = groups(folder, panel)
        if truth.keys() & {p[0] for p in photos}:
            raise ValueError("Duplicate calibration image")
        previous = read_csv(folder / f"support-{panel}-predictions.csv")
        for i, (xy, p) in prepared(folder, photos, model, variant):
            if anchored:
                base = [
                    r
                    for r in previous
                    if r["id"] == i
                    and r["method"] == "supported"
                    and r["type"] in ("T", "X")
                ]
                predictions[i] = anchor(xy, p, base)
            else:
                predictions[i] = (xy, p, None)
        for i, _, refs in photos:
            truth[i] = (
                [
                    np.array([p for k, p in refs if k == kind]).reshape(-1, 2)
                    for kind in ("T", "X")
                ],
                positions([r for r in gaps if r["id"] == i]),
            )
        scores = read_csv(folder / f"support-{panel}-scores.csv")
        limits += [
            sum(int(r["detected"]) for r in scores if r["method"] == "supported"),
            sum(r["supported"] != "-" for r in gaps),
        ]
    nref = np.sum([[len(r) for r in refs] for refs, _ in truth.values()], axis=0)
    grid = []
    settings = (
        itertools.product(
            (0.03, 0.065, 0.10, 0.16),
            (0.015, 0.045, 0.075),
            (0.24, 0.36, 0.46),
            (0.055, 0.12, 0.25),
        )
        if anchored
        else itertools.product(
            (0.055, 0.08, 0.12, 0.17, 0.23, 0.31, 0.42),
            (0.035, 0.055, 0.08, 0.12, 0.17, 0.23, 0.31),
        )
    )

    def measure(setting):
        n = flags = 0
        hits = np.zeros(2, dtype=int)
        for i, (xy, p, mask) in predictions.items():
            refs, gaps = truth[i]
            cut = (
                setting
                if mask is None
                else np.where(mask[:, None], [setting[:2]], [setting[2:]])
            )
            events = suppress(xy, p, cut)
            n += len(events)
            if len(gaps):
                flags += (cdist(gaps, positions(events)) <= 4).any(axis=1).sum()
            for j, kind in enumerate(("T", "X")):
                hits[j] += len(
                    match(refs[j], positions([r for r in events if r["type"] == kind]))
                )
        recall = hits / nref
        grid.append(
            dict(
                cutoffs=json.dumps(setting),
                T=int(hits[0]),
                X=int(hits[1]),
                T_reference=int(nref[0]),
                X_reference=int(nref[1]),
                n=n,
                flags=int(flags),
                weak=min(recall),
                mean=np.mean(recall),
                eligible=int(n <= limits[0] and flags <= limits[1]),
            )
        )

    def strongest():
        return max(
            (r for r in grid if r["eligible"]),
            key=lambda r: (r["weak"], r["mean"], -r["n"]),
        )

    for setting in settings:
        measure(setting)
    center = json.loads(strongest()["cutoffs"])
    for factors in itertools.product((0.85, 1.0, 1.15), repeat=len(center)):
        if any(k != 1 for k in factors):
            measure(tuple(round(v * k, 6) for v, k in zip(center, factors)))
    write_csv(work / f"topology-{variant}-{anchored}-grid.csv", grid, list(grid[0]))
    best = strongest()
    print("Selected", best, "previous density/flags", limits, flush=True)
    return json.loads(best["cutoffs"])


def evaluate(folder, panel, model, cutoffs, variant, excluded=()):
    photos, gaps = groups(folder, panel)
    if excluded:
        photos = [photo for photo in photos if photo[0] not in excluded]
        gaps = [r for r in gaps if r["id"] not in excluded]
    if panel in ("prospective", "remainder"):
        assess(folder, panel, photos, gaps)
    events = [
        r
        for r in read_csv(folder / f"support-{panel}-predictions.csv")
        if r["method"] == "supported"
    ]
    base_by_id = {
        i: [r for r in events if r["id"] == i and r["type"] in ("T", "X")]
        for i, _, _ in photos
    }
    for i, (xy, p) in prepared(folder, photos, model, variant):
        for event in resolve(xy, p, base_by_id[i], cutoffs):
            events.append(
                dict(
                    id=i,
                    method="topology",
                    x=event["x"],
                    y=event["y"],
                    type=event["type"],
                    confidence=event["confidence"],
                )
            )
    methods = ("supported", "topology")
    scores, gap_scores = score(
        folder,
        photos,
        gaps,
        prefix=f"topology-{panel}-",
        methods=methods,
        detections=events,
    )
    write_csv(
        folder / f"topology-{panel}-predictions.csv",
        events,
        ["id", "method", "x", "y", "type", "confidence"],
    )
    rows = [totals(panel, m, scores, gap_scores) for m in methods]
    write_csv(folder / f"topology-{panel}-summary.csv", rows, list(rows[0]))
    print(panel, rows, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("model", type=Path)
    parser.add_argument("cutoffs", type=Path)
    parser.add_argument("--select", action="store_true")
    parser.add_argument(
        "--panel",
        choices=(
            "development",
            "replication",
            "confirmation",
            "issue",
            "prospective",
            "remainder",
        ),
    )
    parser.add_argument(
        "--variant", choices=("arms", "photoarms", "multiscale"), default="multiscale"
    )
    parser.add_argument("--direct", action="store_true")
    parser.add_argument(
        "--freeze",
        type=Path,
        help="Write the numerical scorer for the single-photo runner",
    )
    args = parser.parse_args()
    model = json.loads(args.model.read_text())
    excluded = set()
    if args.freeze or args.panel in ("prospective", "remainder"):
        excluded = check_rosters(args.work, include_remainder=args.panel == "remainder")
    if args.select:
        args.cutoffs.write_text(
            json.dumps(select(args.work, model, args.variant, not args.direct)) + "\n"
        )
    if args.freeze:
        freeze(args.freeze, args.variant, json.loads(args.cutoffs.read_text()), model)
    if args.panel:
        folder = (
            args.work
            if args.panel in ("issue", "development")
            else args.work / args.panel
        )
        evaluate(
            folder,
            args.panel,
            model,
            json.loads(args.cutoffs.read_text()),
            args.variant,
            excluded,
        )
