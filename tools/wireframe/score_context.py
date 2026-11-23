"""Compare both frozen scorers on the official York test set, including equal-count rankings."""
import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist

from check_localization import cardinality
from score import ROOT, match, positions, read_csv, supplement, write_csv


def compare(photo_id, panel, refs, gaps, previous, candidate):
    row = dict(panel=panel, id=photo_id, reviewed=len(gaps))
    near = defaultdict(list)
    for kind, truth in zip(("T", "X"), refs, strict=True):
        ranked = [
            sorted(
                (r for r in events if r["type"] == kind),
                key=lambda r: -float(r["score"]),
            )
            for events in (previous, candidate)
        ]
        budget = min(map(len, ranked))
        row[kind + "_reference"] = len(truth)
        for name, events in zip(("baseline", "context"), ranked, strict=True):
            for suffix, subset in (("", events), ("_equal", events[:budget])):
                xy = positions(subset)
                hits = cardinality(truth, xy)
                if hits != len(match(truth, xy)):
                    raise ValueError(f"Matching disagreement for {photo_id}")
                row[f"{kind}_{name}{suffix}"] = hits
                row[f"{kind}_{name}{suffix}_detected"] = len(xy)
                near[name + suffix].append((cdist(gaps, xy) <= 4).any(axis=1))
    for key, masks in near.items():
        row[key + "_flags"] = int(np.logical_or.reduce(masks).sum())
    return row


def york(folder, events_file, output):
    evidence = ROOT / "docs/evidence/wireframe"
    if read_csv(folder / "panel.csv") != read_csv(evidence / "york-panel.csv"):
        raise ValueError("York test roster differs from the freeze")
    with gzip.open(events_file, "rt", newline="") as stream:
        detections = defaultdict(list)
        for row in csv.DictReader(stream):
            detections[row["id"], row["method"]].append(row)
    gaps = defaultdict(list)
    for row in read_csv(evidence / "york-gap-review.csv"):
        if row["include"] == "1":
            gaps[row["id"]].append(row)
    rows = []
    for i, panel, sites in supplement(folder, ("york51",)):
        truth = [
            np.asarray([p for k, p in sites if k == kind]).reshape(-1, 2)
            for kind in ("T", "X")
        ]
        rows.append(
            compare(
                i,
                panel,
                truth,
                positions(gaps[i]),
                detections[i, "topology"],
                detections[i, "context"],
            )
        )
    write_csv(output, rows, list(rows[0]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument("events", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    york(args.folder, args.events, args.output)
