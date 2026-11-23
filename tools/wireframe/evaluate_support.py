"""Evaluate M-LSD additions supported by overlapping crops or local line fragments."""
import argparse
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist

from score import (
    ROOT,
    positions,
    original,
    read_csv,
    score,
    supplement,
    totals,
    write_csv,
)


def rank_crops(enhanced, crops):
    votes = crops + [r for r in enhanced if r["type"] in ("T", "X")]
    xy = positions(crops)
    distances = cdist(xy, positions(votes))
    ranked = []
    for i, row in enumerate(crops):
        nearest = {}
        for j in np.flatnonzero(distances[i] <= 4):
            other = votes[j]
            view = other["method"]
            if view != row["method"] and other["type"] == row["type"]:
                nearest[view] = min(nearest.get(view, np.inf), distances[i, j])
        if nearest:
            ranked.append((-len(nearest), np.mean(list(nearest.values())), i))
    return [crops[i] for _, _, i in sorted(ranked)]


def extend(primary, proposals, method):
    """Keep primary sites; add mutually separated new sites in proposal rank order."""
    primary = [r for r in primary if r["type"] in ("T", "X")]
    xy = positions(proposals)
    blocked = (cdist(xy, positions(primary)) <= 4).any(axis=1)
    kept = []
    for i in np.flatnonzero(~blocked):
        if len(kept) >= len(primary) // 5:
            break
        if all(np.linalg.norm(xy[i] - xy[j]) > 4 for j in kept):
            kept.append(i)
    return [dict(r, method=method) for r in primary + [proposals[i] for i in kept]]


def assess(folder, split, photos, reviews):
    full = read_csv(folder / "junctions.csv")
    tiles = read_csv(folder / "crop-junctions.csv")
    fragments = read_csv(folder / "fragment-candidates.csv")
    events = []
    for i, _, _ in photos:
        base = [r for r in full if r["id"] == i and r["method"] == "mlsd_mask"]
        assisted = [r for r in full if r["id"] == i and r["method"] == "mlsd_clahe"]
        crops = [r for r in tiles if r["id"] == i]
        events.extend(base)
        repeated = rank_crops(assisted, crops)
        events.extend(extend(base, repeated, "crop_consensus"))
        xs = [
            dict(id=i, x=r["x"], y=r["y"], type="X") for r in fragments if r["id"] == i
        ]
        events.extend(extend(base, xs, "fragment_x"))
        events.extend(extend(base, xs + repeated, "supported"))
    methods = ("mlsd_mask", "crop_consensus", "fragment_x", "supported")
    write_csv(
        folder / f"support-{split}-predictions.csv",
        events,
        ["id", "method", "x", "y", "type"],
    )
    scores, negatives = score(
        folder,
        photos,
        reviews,
        prefix=f"support-{split}-",
        methods=methods,
        detections=events,
    )
    return [totals(split, method, scores, negatives) for method in methods]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("replication", type=Path)
    parser.add_argument("--confirmation", type=Path, help="48 further unused photos")
    args = parser.parse_args()
    issue, issue_negatives = original(args.directory, args.manifest)
    issue_ids = {i for i, _, _ in issue}
    dev = [p for p in supplement(args.directory, ("sample",)) if p[0] not in issue_ids]
    dev_ids = {p[0] for p in dev}
    reviews = ROOT / "docs/evidence/wireframe"
    dev_negatives = [
        r for r in read_csv(reviews / "near-miss-review.csv") if r["id"] in dev_ids
    ]
    rows = assess(args.directory, "development", dev, dev_negatives)
    rows += assess(args.directory, "issue", issue, issue_negatives)
    fresh = supplement(args.replication, ("replication",))
    rows += assess(
        args.replication,
        "replication",
        fresh,
        read_csv(reviews / "replication-gap-review.csv"),
    )
    if args.confirmation:
        heldout = supplement(args.confirmation, ("confirmation",))
        ids = {p[0] for p in heldout}
        if len(ids) != 48 or ids & (issue_ids | dev_ids | {p[0] for p in fresh}):
            raise ValueError(
                "The confirmation panel must contain 48 distinct unused photos"
            )
        rows += assess(
            args.confirmation,
            "confirmation",
            heldout,
            read_csv(reviews / "confirmation-gap-review.csv"),
        )
    write_csv(args.directory / "support-summary.csv", rows, list(rows[0]))
