"""Score blinded card ratings and exhaustive region marks against fixed predictions."""
import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist

from score import match, positions, read_csv


METHODS = ("context", "scores", "working", "native")
CLASSES = {"T", "X", "gap", "other", "unclear"}


def card_scores(key, ratings, events):
    if set(ratings) != {r["token"] for r in key} or set(ratings.values()) - CLASSES:
        raise ValueError("A valid decision is required on every card")
    out = {}
    for m in METHODS:
        for kind in "TX":
            den = sum(r[m] == kind for r in events)
            num = sum(
                int(r["population"]) / int(r["sample"])
                for r in key
                if r[m] == kind and ratings[r["token"]] == kind
            )
            out[m + "_" + kind] = num / den if den else None
    return out


def region_scores(key, regions, events):
    if set(regions) != {r["token"] for r in key}:
        raise ValueError("Every unmarked region needs a review")
    by_image = defaultdict(list)
    for r in events:
        by_image[r["id"]].append(r)
    out = {
        m: {
            k: 0
            for k in (
                "T",
                "X",
                "T_reference",
                "X_reference",
                "gaps",
                "flags",
                "X_flags",
            )
        }
        for m in METHODS
    }
    for r in key:
        review = regions[r["token"]]
        if review["complete"] is not True:
            raise ValueError(f"Incomplete region {r['token']}")
        points = defaultdict(list)
        for mark in review["marks"]:
            x, y, label = mark["x"], mark["y"], mark["label"]
            if not (
                12 <= x < 228
                and 12 <= y < 188
                and label in ("T", "X", "gap", "unclear")
            ):
                raise ValueError(f"Invalid mark in {r['token']}")
            # Canvas clicks are zero-based indices within the original-pixel crop.
            xy = np.array([int(r["left"]) + x + 0.5, int(r["top"]) + y + 0.5])
            scale = np.array(
                [
                    int(r["small_width"]) / int(r["width"]),
                    int(r["small_height"]) / int(r["height"]),
                ]
            )
            points[label].append(xy * scale + 0.5)
        for m in METHODS:
            sample = by_image[r["id"]]
            for kind in "TX":
                ref = np.asarray(points[kind]).reshape(-1, 2)
                pred = positions([e for e in sample if e[m] == kind])
                out[m][kind] += len(match(ref, pred))
                out[m][kind + "_reference"] += len(ref)
            gaps = np.asarray(points["gap"]).reshape(-1, 2)
            out[m]["gaps"] += len(gaps)
            for field, sel in (
                ("flags", sample),
                ("X_flags", [e for e in sample if e[m] == "X"]),
            ):
                out[m][field] += int((cdist(gaps, positions(sel)) <= 4).any(1).sum())
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for field in ("keys", "events", "ratings"):
        parser.add_argument(field, type=Path)
    args = parser.parse_args()
    data = json.loads(args.ratings.read_text())
    if data["version"] != 1:
        raise ValueError("Unknown rating format")
    events = read_csv(args.events)
    print(
        json.dumps(
            dict(
                card_precision=card_scores(
                    read_csv(args.keys / "cards.csv"), data["cards"], events
                ),
                regions=region_scores(
                    read_csv(args.keys / "regions.csv"), data["regions"], events
                ),
            ),
            indent=2,
        )
    )
