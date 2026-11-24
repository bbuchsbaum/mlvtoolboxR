"""Count which member paths survive selection at matched M-LSD development sites."""

import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np

from contour_incidence import incident_correspondence
from score import ROOT, match, positions, read_csv, write_csv


def geometry_key(contour):
    return tuple(np.asarray(contour["segments"], dtype=float).ravel())


def audit(staged, base, predictions):
    counts = Counter()
    for row in read_csv(ROOT / "docs/evidence/wireframe/curve-rank/development.csv"):
        name = row["id"]
        artist = json.loads((staged / f"{name}.json").read_text())
        source = json.loads((base / f"{name}.mlsd.contours.json").read_text())
        candidate = json.loads(
            (predictions / f"{name}.curve_rank.contours.json").read_text()
        )
        kept = {geometry_key(c) for c in candidate["contours"]}
        original = {c["id"]: geometry_key(c) for c in source["contours"]}
        strict = incident_correspondence(artist, source["sites"], source["contours"])
        for kind in ("T", "X"):
            sites = [s for s in source["sites"] if s["type"] == kind]
            truth = [s["xy"] for s in artist["junctions"] if s["type"] == kind]
            for criterion, pairs in [
                ("position", match(truth, positions(sites))),
                ("incidence", strict[kind]["pairs"]),
            ]:
                for _, j in pairs:
                    members = {original[m["contour_id"]] for m in sites[j]["members"]}
                    count = len(members & kept)
                    retention = (
                        "all" if count == len(members) else "some" if count else "none"
                    )
                    counts[kind, criterion, retention] += 1
    return [
        dict(
            type=kind,
            criterion=criterion,
            **{k: counts[kind, criterion, k] for k in ("all", "some", "none")},
        )
        for kind in ("T", "X")
        for criterion in ("position", "incidence")
    ]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ("staged", "base", "predictions", "output"):
        parser.add_argument(argument, type=Path)
    args = parser.parse_args()
    records = audit(args.staged, args.base, args.predictions)
    write_csv(args.output, records, list(records[0]))
